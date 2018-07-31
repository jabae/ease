function Aem_proj = extract_em_segments(obj, mscan, mslice)
%% load the projected data
%{
%}

%% inputs: 
%{
    mscan: scan ID 
    mslice: slice ID 
%}

%% outputs: 
%{
    Aem_proj: project of each EM segment on the selected slice 
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

%% code 
 % choose loading data from matfile directly or from the matlab struct variable
if obj.em_load_flag        
    if ~exist_in_workspace('em_data_mem', 'base')
        obj.load_em();
    end
end

% if mscan/mslice is not specified, then use the current scan/slice id. 
if ~exist('mscan', 'var')
    mscan = obj.scan_id;
end
if ~exist('mslice', 'var')
    mslice = obj.slice_id;
end

% EM info
blur_size = obj.em_zblur;   % blurring in z direction
shifts_em_ii = obj.em_shifts.ii;
shifts_em_jj = obj.em_shifts.jj;

% create a mat file for storing results
matfile_proj = fullfile(obj.output_folder, ...
    sprintf('Aem_proj_blur%d.mat', blur_size));
if ~exist(matfile_proj, 'file')
    flag_processed = false(obj.num_scans, obj.num_slices);
    save(matfile_proj, 'blur_size', 'flag_processed', '-v7.3');
else
    temp = load(matfile_proj, 'flag_processed');
    flag_processed = temp.flag_processed;
end

% check whether the data has been projected already
var_name = sprintf('scan%d_slice%d', mscan, mslice);
zplane = obj.video_zvals_updated(mscan, mslice);    % current planes

if flag_processed(mscan, mslice)
    fprintf('plane: scan=%d, slice=%d, z=%d...\n', mscan, mslice, zplane);
    temp = load(matfile_proj, var_name);  %#ok<NASGU>
    Aem_proj = eval(sprintf('temp.%s', var_name));
    return;
end

% start projecting plane by plane
fprintf('\nprojecting all EM masks to the selected plane.\n'); 

% stack information: cropped area and size
d1_2p = obj.d1*2;
d2_2p = obj.d2*2;
FOV_2p = obj.FOV_stack;

% pre-allocate a space for saving results
dims_Aem = [d1_2p*d2_2p, obj.K_em];
Aem_proj = zeros(dims_Aem);

% crop the EM masks
dy = shifts_em_ii(mscan, mslice);
dx = shifts_em_jj(mscan, mslice);

% overlapping nearby frames
z0 = max(zplane-2*blur_size, 1);
z1 = min(zplane+2*blur_size, obj.dims_stack(3));
fprintf('running......\n');
for z=z0:z1
    fprintf('|');
end
fprintf('\n\n');
for z=z0:z1
    if isempty(obj.em_ranges{z})
        % not in the EM volume
        continue;
    else
        % load the EM representation on the slected plane
        if obj.em_load_flag
            tmp_slice = evalin('base', sprintf('em_data_mem.slice_%d', z));
        else
            tmp_slice = eval(sprintf('obj.em_data.slice_%d', z));
        end
        
        % get pixel locations and represent its cooridiantes
        % in the cropped area.
        [ii, jj] = ind2sub(obj.dims_stack(1:2), tmp_slice(:,1));
        ii = ii - dy - FOV_2p(1) + 1;
        jj = jj - dx - FOV_2p(3) + 1;
        ind_cropped = sub2ind([d1_2p, d2_2p], ii, jj);
        
        % uncompress the data, find all nonzero pixels
        tmpA = cell(obj.PACK_SIZE, 1);
        pack0 = (tmp_slice(:,2)-1) * 64;
        for nn=1:obj.PACK_SIZE
            ind_nz = logical(bitget(tmp_slice(:,3), nn));
            tmpA{nn} = sub2ind(dims_Aem, ind_cropped(ind_nz, 1), pack0(ind_nz)+nn);
        end
        tmpA = cell2mat(tmpA);
        
        % project data to the selected scanning plane with some
        % blurring weights
        tmp_weight = exp(-(z-zplane)^2 / (2*blur_size^2) );
        Aem_proj(tmpA) = Aem_proj(tmpA) + tmp_weight;
    end
    fprintf('.');
end
fprintf('\n');

% save the results
Aem_proj = sparse(Aem_proj);
eval(sprintf('%s = Aem_proj;', var_name));
flag_processed(mscan, mslice) = true;  %#ok<NASGU>
save(matfile_proj, var_name, 'flag_processed','-append');
fprintf('Done! The results were saved into the mat file \n %s\nand named as %s\n', ...
    matfile_proj, var_name);
end