function construct_Aem(obj)

%% construct Aem \in d*K for each scan

%% input

%% output

%% create a folder for storing Aem
dir_em = fullfile(obj.output_folder, sprintf('segmentation_%d_zblur_%d', ...
    obj.em_segmentation, obj.em_zblur));
if ~exist(dir_em, 'dir')
    mkdir(dir_em);
else
    fprintf('The folder for storing EM footprints has been created. \n delete the folder if you want to do start it over.\n');
    temp = input('delete the old one (y/n): ', 's');
    if strcmpi(temp, 'y')
        rmdir(dir_em, 's');
        mkdir(dir_em);
    else
        return;
    end
end

%% order neurons and divide neurons into multiple groups.
min_pixel_number = 0; 
if isempty(obj.blur_version)
    obj.set_projections_options(); 
end
if strcmpi(obj.data_name, 'pinky40')
    rel0 = obj.rel_footprints &...
        (ta3.Mesh & (ta3.VoxelizedMesh & 'n_vertices>2500')) &...
        sprintf('segmentation=%d', obj.em_segmentation) & ...
        sprintf('version=%d', obj.blur_version);
    sprintf('version=%d', obj.blur_version);
    obj.em_shifts.ii = obj.stack_shifts.ii - 3;
    obj.em_shifts.jj = obj.stack_shifts.jj - 3;
else
    rel0 = obj.rel_footprints & ...
        sprintf('segmentation=%d', obj.em_segmentation) & ...
        sprintf('version=%d', obj.blur_version);
end
n_voxels = rel0.fetchn('n_voxels');  % fetch number of voxels
n_voxels(n_voxels<min_pixel_number) = [];
n_voxels = sort(n_voxels, 'descend');
Kem = length(n_voxels);

% find the break points to split data
sum_voxels = cumsum(n_voxels);
total_voxels = sum_voxels(end);
batch_size = 10^7;
break_points = batch_size:batch_size:total_voxels;
nbatch = length(break_points)+1;
indices = ones(nbatch, 2);
for m=2:nbatch
    idx = find(sum_voxels>break_points(m-1), 1, 'first');
    indices(m-1, 2) = idx;
    indices(m, 1) = idx + 1;
end
indices(end, 2) = Kem;
bp_voxels = n_voxels(indices(:,2));

%% fetch and data and construct Aem for each scan
[~, dims_new] = obj.convert_idx([]);
dvoxel = prod(dims_new);

% create indices for selecting given planes
scan_select = false(prod(dims_new), obj.num_scans);
for m=1:obj.num_scans
    ind = false(dims_new);
    ind(:, :, (1:obj.num_slices)+(m-1)*obj.num_slices) = true;
    scan_select(:, m) = ind(:);
end
%%
for mbatch=1:nbatch
    fprintf('batch %d/%d\n', mbatch, nbatch);
    
    % select the data
    if mbatch==1
        rel = rel0 & sprintf('n_voxels>=%d', bp_voxels(1));
    else
        rel = rel0 & ...
            sprintf('n_voxels<%d', bp_voxels(mbatch-1)) & ...
            sprintf('n_voxels>=%d', bp_voxels(mbatch));
    end
    
    % fetch data
    [segment_id, idx_values] = fetchn(rel, 'segment_id', 'idx_value');
    num_cell = length(segment_id);
    
    % save segment_ids
    file_name = fullfile(dir_em, sprintf('batch_%d_ids.mat', mbatch));
    save(file_name, 'segment_id', '-v7.3');
    
    %% construct Aem
    ii = cell(num_cell,1);
    jj = cell(num_cell, 1);
    vv = cell(num_cell, 2);
    for mcell=1:num_cell
        [ii{mcell}, ~] = obj.convert_idx(idx_values{mcell}(:, 1));
        jj{mcell} = ones(size(ii{mcell})) * mcell;
        vv{mcell} = idx_values{mcell}(:, 2);
    end
    ii = cell2mat(ii);
    jj = cell2mat(jj);
    vv = cell2mat(vv);
    
    Aem_all = sparse(ii, jj, vv, dvoxel,num_cell);
    if mbatch==1
        em_volume = full(sum(Aem_all,2))>0;
    else
        em_volume = or(em_volume, full(sum(Aem_all,2))>0);
    end
    %% split Aem into the given number of scans
    for mscan=1:obj.num_scans
        idx = scan_select(:, mscan);
        Aem = Aem_all(idx, :);
        file_name = fullfile(dir_em, sprintf('batch_%d_scan_%d_Aem.mat', mbatch, mscan));
        save(file_name, 'Aem', '-v7.3');
        fprintf('--scan %d/%d\n', mscan, obj.num_scans);
    end
    % resolutions
end
fprintf('The footprints of all EM segments were projected & aligned onto the scanning planes of calcium imaging data.\n');

%% save the spatial range
em_volume = reshape(em_volume, [], obj.num_scans);

% run an morphological erosion to update em_volume 

file_name = fullfile(dir_em, 'em_volume.mat');
save(file_name, 'em_volume', '-v7.3');

ease.em_volume = em_volume; 