function init_project(obj, dir_project)
%% initialize a projct folder for running EASE
%{
	This function will create a project folder to store all
	data/code/results for running EASE analysis. The following folder will
	be created.
    - /scripts/
    - /data/
    - /results/
    - /Videos/
    - /Figures/
    - /schemas
%}

%% inputs
%{
	obj: class object
	dir_project: string; a folder for stroing everything. if it's missing
	or empty, EASE will ask you to select one.
%}

%% outputs
%{
%}

%% Author
%{
	Pengcheng Zhou
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License
%}


%% create a project folder
fprintf('create a project folder for storing the data and analysis output.\n'); 
if ~exist('dir_project', 'var') || isempty(dir_project)
    obj.dir_project = uigetdir('', 'choose or create a folder');
end

%% create subfolders
subfolder_list = {'scripts', 'data', 'results', 'Figures', 'Videos', 'schemas'};
for m=1:length(subfolder_list)
    tmp_folder = fullfile(obj.dir_project, subfolder_list{m});
    if ~exist(tmp_folder, 'dir')
        mkdir(tmp_folder);
    end
end

%% create a mat file to store important information
project_info = fullfile(obj.dir_project, 'metainfo.yaml');
if ~exist(project_info, 'file')
    metainfo.datasets_list = {'pinky100'};
    metainfo.databases_list = {'127.0.0.1:3306'};
    metainfo.example_data = struct('datajoint_name', 'ta3p100', ...
        'database_name', 'microns_ta3p100', 'rel_mesh', ta3p100.Mesh);
    yaml.WriteYaml(project_info, metainfo);
end

%% add this project to EASE directory
tmp_file = fullfile(fi.locate('ease', true), '.projects.mat');
if exist(tmp_file, 'file')
    load(tmp_file, 'projects');
else
    projects = {};
end
for m=1:length(projects)
    if strcmp(projects{m}, obj.dir_project)
        return;
    end
end
projects{end+1} = obj.dir_project;
if exist(tmp_file, 'file')
    save(tmp_file, 'projects', '-append');
else
    save(tmp_file, 'projects');
end
