function Aem_all = project_em_to_stack(obj)
%% load the projected data
if obj.em_load_flag         % choose loading data from matfile directly or from the matlab struct variable
    if ~exist_in_workspace('em_data_mem', 'base')
        obj.load_em();
    end
end

% EM info
blur_size = obj.em_zblur;   % blurring in z direction
shifts_em_ii = obj.em_shifts.ii;
shifts_em_jj = obj.em_shifts.jj;

% create a mat file for storing results
matfile_proj = fullfile(obj.output_folder, ...
    sprintf('em_stack_proj.mat', blur_size));
if ~exist(matfile_proj, 'file')
    flag_processed = false; 
    save(matfile_proj, 'flag_processed', '-v7.3');
else
    temp = load(matfile_proj, 'flag_processed');
    flag_processed = temp.flag_processed;
end

% check whether the data has been projected already
if flag_processed
    % the data has been processed 
    fprintf('all EM masks were projected to 2p stack already...\n'); 
    temp = load(matfile_proj, 'Aem_all'); 
    Aem_all = temp.Aem_all; 
else
    % load data
    fprintf('\nprojecting all EM masks to the selected plane.\n');
    
    % stack information: cropped area and size
    d1_2p = obj.d1*2;
    d2_2p = obj.d2*2;
    d3_2p = obj.dims_stack(3); 
    FOV_2p = obj.FOV_stack;
    
    % pre-allocate a space for saving results
    dims_Aem = [d1_2p*d2_2p, obj.K_em];
    Aem_all = cell(d3_2p, 1); 
    
    % crop the EM masks
    dy = 0; %shifts_em_ii(mscan, mslice);
    dx = 0; %shifts_em_jj(mscan, mslice);
    
    % overlapping nearby frames
    z0 = 1;
    z1 = d3_2p;
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
            
            % uncompress the data find all nonzero pixels
            tmpA = cell(obj.PACK_SIZE, 1);
            pack0 = (tmp_slice(:,2)-1) * 64;
            for nn=1:obj.PACK_SIZE
                ind_nz = logical(bitget(tmp_slice(:,3), nn));
                tmpA{nn} = sub2ind(dims_Aem, ind_cropped(ind_nz, 1), pack0(ind_nz)+nn);
            end
            tmpA = cell2mat(tmpA);
            
            % project data to the selected scanning plane with some
            % blurring weights
            [r, c] = ind2sub(dims_Aem, tmpA);
            Aem_all{z} = sparse(r, c, 1, dims_Aem(1), dims_Aem(2)); 
        end
        fprintf('.');
    end
    fprintf('\n');

end
flag_processed = true; 
%% blurr data 
% Aem_blur = imfilter(Aem_proj, fspecial('gaussian', [1, 5], blur_size)); 
% Aem_blur = reshape(Aem_blur, d1_2p, d2_2p, obj.K_em, d3_2p); 
% start projecting plane by plane 
% save the results
% Aem_proj = sparse(Aem_proj);
save(matfile_proj, 'Aem_all', 'flag_processed', 'dims_Aem', 'd3_2p', '-append');
fprintf('Done! The results were saved into the mat file \n %s\n.\n', ...
    matfile_proj);
end