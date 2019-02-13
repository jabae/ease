%% project all EM meshes to the scanning planes of calcium imaging data 

%% code
% choose loading data from matfile directly or from the matlab struct variable
if ease.em_load_flag
    if ~exist_in_workspace('em_data_mem', 'base')
        ease.load_em();
    end
end

% stack information: cropped area and size
FOV_2p = ease.FOV_stack;
d1_2p = diff(FOV_2p(1:2))+1; 
d2_2p = diff(FOV_2p(3:4))+1;

% EM info
blur_size = ease.em_zblur;   % blurring in z direction
shifts_em_ii = ease.em_shifts.ii;
shifts_em_jj = ease.em_shifts.jj;

% create a mat file for storing results
matfile_proj = fullfile(ease.output_folder, ...
    sprintf('Aem_proj_blur%d.mat', blur_size));
if ~exist(matfile_proj, 'file')
    flag_processed = false(ease.num_scans, ease.num_slices);
    save(matfile_proj, 'blur_size', 'flag_processed', '-v7.3');
else
    load(matfile_proj, 'flag_processed'); 
    fprintf('The EM meshes have been created. If you want to project them again,\ndelete %s first\n', matfile_proj); 
end

%% load all nonzero voxels 
ssub =  ease.dims_stack(1) / ease.dims_video(1); 
K_em = ease.K_em; 
ease_connect_database; 
temp = rel_voxels & sprintf('segmentation=%d', ease.em_segmentation); 
indices = cell2mat(temp.fetchn('indices')); 
n_voxels = ease.em_data.EM_info(:, 2); 
ids = zeros(length(indices), 1); 
temp = cumsum(n_voxels) + 1; 
ids(1) = 1; 
ids(temp(1:(end-1))) = 1; 
ids = cumsum(ids);      % segment_ids for each voxel 
[y_voxel, x_voxel, z_voxel] = ind2sub(ease.dims_stack, indices); 

%% project EM meshes
norm_z = 2 *blur_size^2; 

for mscan = 1:ease.num_scans
    for mslice=1:ease.num_slices
        fprintf('\nprojecting all EM masks to (scan-%d, slice-%d).\n', mscan, mslice);
        if flag_processed(mscan, mslice)
            disp([mscan, mslice]); 
            continue;
        end
        var_name = sprintf('scan%d_slice%d', mscan, mslice); 
        % shift information 
        dx = shifts_em_jj(mscan, mslice); 
        dy = shifts_em_ii(mscan, mslice); 
        z0 = ease.video_zvals_updated(mscan, mslice); 
        
        %
        ii = y_voxel - dy - FOV_2p(1) + 1;
        jj = x_voxel - dx - FOV_2p(3) + 1;
        ind_cropped = sub2ind([d1_2p, d2_2p], ii, jj);
        
        % neuron by neuron 
        temp = abs(z_voxel-z0); 
        temp(temp > 2*blur_size) = inf; 
        
        % create a sparse matrix
        tmp_slice = sparse(ind_cropped, ids, exp(-temp.^2 /norm_z), d1_2p*d2_2p, K_em);
        eval(sprintf('scan%d_slice%d = sparse(tmp_slice);', mscan, mslice));
        flag_processed(mscan, mslice) = true;
        if ssub > 1
            d1 = d1_2p / ssub; 
            d2 = d2_2p / ssub; 
            ind_cropped = sub2ind([d1_2p, d2_2p]/ssub, ceil(ii/ssub), ceil(jj/ssub));
            tmp_slice_ds = sparse(ind_cropped, ids, exp(-temp.^2 /norm_z), d1*d2, K_em);
            var_name_ds = sprintf('%s_ds%d', var_name, ssub);
            eval(sprintf('%s=tmp_slice_ds;', var_name_ds));
            save(matfile_proj, var_name, var_name_ds, 'flag_processed','-append');
        else
            save(matfile_proj, var_name, 'flag_processed','-append');
        end
        disp([mscan, mslice]);      
    end
end
