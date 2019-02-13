%% setup packages and some configurations
close all; clear; clc;

with_GUI = false;

% folders related to the project
dir_scripts = fileparts(which('start_EASE.m'));
dir_project = fileparts(dir_scripts);
dir_data = fullfile(dir_project, 'data');   % place for storing data
dir_results = fullfile(dir_project, 'results'); % place for storing results
dir_fig = fullfile(dir_project, 'Figures'); % place for storing figures
dir_video = fullfile(dir_project, 'Videos'); % place fore storing videos

% packages to be used
fi.usepkg({'yaml', ... % YAML matlab
    'utils', ... % export matlab
    'idl', ... % ImagingDataLoader
    'oasis', ... % OASIS_matlab
    'ease'}... % EASE };
    );

%% choose the data to use
datasets = {'pinky40', 'pinky100'};
fprintf('\n**********choose the data to use**********\n');
for m=1:length(datasets)
    fprintf('%d: %s\n', m, datasets{m});
end
fprintf('********************************************\n');

data_id = input('data ID: ');
while true
    if any(data_id==[1, 2])
        data_name = datasets{data_id};
        fprintf('you selected data %s\n', data_name);
        break;
    else
        data_id = input('please type a valid data ID: ');
    end
end


%% connect to the database
ease_connect_database;
if strcmpi(data_name, 'pinky40')
    rel_mesh = ta3.MeshFragment;
    rel_voxels = ta3.VoxelizedMesh;
elseif strcmpi(data_name, 'pinky100')
    rel_mesh = ta3p100.Mesh;
    rel_voxels = ta3p100.VoxelizedMesh;
else
    error('the selected data has not saved in the database.');
end

%% create a class object to manage options and pipelines
yaml_path = fullfile(dir_scripts, sprintf('%s_config.yaml', data_name));
if ~exist(yaml_path, 'file')
    % create from default options 
    ease = EM2P();
    ease.write_config(yaml_path);
else
    ease = EM2P(yaml_path);
end

% update data information
ease.output_folder = fullfile(dir_results, data_name);
ease.data_folder = fullfile(dir_data, data_name);
ease.fig_folder = fullfile(dir_fig, data_name);
ease.video_folder = fullfile(dir_video, data_name);
ease.matfile_stack = 'stack_2p.mat';
ease.matfile_video = 'functional_data.mat';
ease.denoised_folder = 'cropped_denoised_video';
ease.raw_folder = 'cropped_raw_video';
ease.registration_csv = 'registration.csv';
ease.matfile_transformation = 'coor_convert.mat';
ease.matfile_em = 'em.mat';
ease.matfile_stack = 'stack_2p.mat';

if ~exist(ease.output_folder, 'dir')
    mkdir(ease.output_folder);
end
if ~exist(ease.fig_folder, 'dir')
    mkdir(ease.fig_folder);
end
if ~exist(ease.video_folder, 'dir')
    mkdir(ease.video_folder);
end

%% load data
ease_get_ready; 
