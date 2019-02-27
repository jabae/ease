%% create alias for database tables. 
if strcmpi(data_name, 'pinky40')
    rel_mesh = ta3.MeshFragment;
    rel_voxels = ta3.VoxelizedMesh;
elseif strcmpi(data_name, 'pinky100')
    rel_mesh = ta3p100.Mesh;
    rel_voxels = ta3p100.VoxelizedMesh;
else
    error('the selected data has not saved in the database.');
end

%% load stack data
ease.load_stack();      % stack data
if ~exist(fullfile(ease.data_folder, ease.matfile_em), 'file')
    ease_voxelize_em;
end

%% load EM data 
ease.load_em();         % EM data
ease_select_em; 
EM_info = ease.em_data.EM_info;

%% load video data 
ease.load_video();      % video data

%% get EM boundary
ease.get_em_boundaries();

%% choose FOV
show_fov = false;
ease.choose_FOV(show_fov);

%% crop video
ease_crop_video;

%% align video data and stack data
ease.align_video_stack();
if strcmpi(data_name, 'pinky40')
    ease.em_shifts.ii = ease.stack_shifts.ii - 3;
    ease.em_shifts.jj = ease.stack_shifts.jj+2;
else
    ease.em_shifts.ii = ease.stack_shifts.ii;
    ease.em_shifts.jj = ease.stack_shifts.jj;
end
Y_raw = cell(ease.num_scans, ease.num_blocks);
Y_denoised = cell(ease.num_scans, ease.num_blocks);

%% save the current configurations
ease.write_config(yaml_path);

%% GUI
if exist('with_GUI', 'var')  && with_GUI
    ease.startGUI();
end

%% create a flag indicating EASE is ready to be used.
flag_ease_running = true;
