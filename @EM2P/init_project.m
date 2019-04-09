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
if ~exist('dir_project', 'var') || isempty(dir_project) 
    obj.dir_project = uigetdir();
end

%% create subfolders
subfolder_list = {'scripts', 'data', 'results', 'Figures', 'Videos'};
for m=1:length(subfolder_list)
    tmp_folder = fullfile(obj.dir_project, subfolder_list{m});
    if ~exist(tmp_folder, 'dir')
        mkdir(tmp_folder);
    end
end

%% create a mat file to store important information
project_info = fullfile(obj.dir_project, 'metainfo.yaml');
if ~exist(project_info, 'file')
    metainfo.datasets_list = {'example_data'};
    metainfo.databases_list = {'example_host'};
    metainfo.example_data = struct('datajoint_name', 'ta3'); 
    yaml.WriteYaml(project_info, metainfo);
end

%% add this project to EASE directory
tmp_file = fullfile(fi.locate('ease', true), '.projects.mat');
load(tmp_file, 'projects');
for m=1:length(projects)
    if strcmp(projects{m}, obj.dir_project)
        return;
    end
end
projects{end+1} = obj.dir_project;
save(tmp_file, 'projects', '-append');

