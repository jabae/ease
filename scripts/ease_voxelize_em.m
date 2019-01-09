% setup path and connect to a database 
addpath(fullfile(EASE_dir, 'packages'));
addpath(genpath(fullfile(EASE_dir, 'packages', 'polygon2voxel')));

if ~exist('dj_connected', 'var') || ~dj_connected
    ease_connect_database;
end
%% get the transformation matrix between EM space and 2P space 
if ~exist(fullfile(ease.data_folder, ease.matfile_transformation), 'file')
    ease_get_transformation; 
else
    load(fullfile(ease.data_folder, ease.matfile_transformation)); 
end

% create a struct variable storing parameters
options.voxel_em = [256, 256, 240]/1000; % voxelize EM data with this resolution 
options.dims_2p = ease.dims_stack;
options.range_2p = ease.range_2p;
options.A = A_convert;
options.offset = offset;
options.user_parallel = true; 
if strcmpi(data_name, 'pinky40')  % convert the EM unit to um
    options.scale_factor = 0.001*3.58/4;      
else
    options.scale_factor = 0.001;      
end

%% voxelize all EM meshes 
parpopulate(rel_voxels);
%% collect EM information 
[segment_id, indices] = fetchn(rel_voxels, 'segment_id', 'indices');
EM_info = [segment_id, cellfun(@length, indices)];

% determine EM ranges 
[y, x, z] = ind2sub(options.dims_2p, unique(cell2mat(indices))); 
em_ranges = cell(options.dims_2p(3), 1);
for zz=1:options.dims_2p(3)
    ind = (z==zz); 
    if ~isempty(ind)
        tmp_y = y(ind);
        tmp_x = x(ind);
        k = boundary(tmp_x,tmp_y,0);
        em_ranges{zz} = [tmp_y(k), tmp_x(k)];       
    end
    if mod(zz, 10)==0
        disp(zz); 
    end
end


EM_data.EM_info = EM_info; 
EM_data.options = options; 
EM_data.em_ranges = em_ranges;

save(fullfile(ease.data_folder, ease.matfile_em), 'EM_info', 'options', 'em_ranges',...
    '-v7.3'); 