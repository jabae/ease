function data_name = select_data(obj, datasets)
%% select data for running EASE
%{%}

%% inputs
%{
	obj: class object
	datasets: cell array; each element is a string identifying data
%}

%% outputs
%{
    data_name: string; the selected data name
%}

%% Author
%{
	Pengcheng Zhou
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License
%}

%% select data

if ~exist('datasets', 'var')
    datasets = obj.datasets_list;
else
    datasets = union(obj.datasets_list, datasets);
end

if isempty(datasets) || (length(datasets)== 1 && strcmpi(datasets, 'example_data'))
    edit(fullfile(obj.dir_project, 'metainfo.yaml'));
    warning('no datasets available for this project.');
end
fprintf('\n**********choose the data to use**********\n');
fprintf('0: add new dataset\n');
for m=1:length(datasets)
    fprintf('%d: %s\n', m, datasets{m});
end
fprintf('********************************************\n');

data_id = input('data ID: ');
while true
    if any(data_id==(1:length(datasets)))
        data_name = datasets{data_id};
        fprintf('you selected data %s\n', data_name);
        break;
    elseif data_id==0
        obj.add_dataset();
        
        % make the choice
        temp = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
        obj.datasets_list = temp.datasets_list;
        datasets = obj.datasets_list;
        fprintf('\n**********choose the data to use**********\n');
        fprintf('0: add new dataset\n');
        for m=1:length(datasets)
            fprintf('%d: %s\n', m, datasets{m});
        end
        fprintf('********************************************\n');
        
    else
        data_id = input('please type a valid data ID: ');
    end
end

obj.data_name = data_name;

%% configuration file
obj.yaml_path = fullfile(obj.dir_project, sprintf('%s_config.yaml', data_name));
if exist(obj.yaml_path, 'file')
    obj.read_config();
else
    % assign folders corresponding this this data
    obj.data_folder = fullfile(obj.dir_project, 'data', data_name);
    obj.fig_folder = fullfile(obj.dir_project, 'Figures', data_name);
    obj.video_folder = fullfile(obj.dir_project, 'Videos', data_name);
end

folder_list = {'scripts', 'data', 'results', 'Figures', 'Videos'};
name_list = {'script_folder', 'data_folder', 'output_folder', 'fig_folder', 'video_folder'};
for m=1:length(name_list)
    % folder for scripts
    if isempty(obj.(name_list{m}))
        obj.(name_list{m}) = fullfile(obj.dir_project, folder_list{m}, data_name);
    end
    if     ~exist( obj.(name_list{m}) , 'dir')
        mkdir( obj.(name_list{m}) );
    end
end

%% datajoint
metainfo = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
temp = metainfo.(obj.data_name);
obj.dj_name = temp.datajoint_name;
obj.db_name = temp.database_name;

% create a schema
if ~exist(fullfile(obj.dir_project, 'schemas', ['+',obj.dj_name]), 'dir')
    obj.create_schema();
end

if isfield(temp, 'em_scale_factor')
    obj.em_scale_factor = temp.em_scale_factor;
else
    obj.em_scale_factor = 0.001;
end

if isfield(temp, 'rel_mesh')  % high resolution meshes
    eval(sprintf('obj.rel_mesh=%s;', temp.rel_mesh));
else
    eval(sprintf('obj.rel_mesh=%s.Mesh;', temp.datajoint_name));
end
try  % voxelized meshes
    eval(sprintf('obj.rel_voxels=%s.VoxelizedMesh;', obj.dj_name));
catch
    fprintf('The meshes has not been voxelized yet.\n');
end
try  % footprints
    eval(sprintf('obj.rel_footprints=%s.FootprintsEM;', obj.dj_name));
catch
    fprintf('The EM footprints were not generated yet');
end

%% choose segmentation data
try
    %eval(sprintf('rel_seg = %s.Segmentation;', obj.dj_name));
    %temp = fetchn(rel_seg, 'segmentation');
    eval(sprintf('rel_seg = %s.Segment;', obj.dj_name));
    temp = unique(fetchn(rel_seg, 'version'));
catch
    segmentation = fetchn(obj.rel_mesh, 'version');
    temp = unique(segmentation);
end
if length(temp)>1
    fprintf('there are %d versions of mesh: ', length(temp));
    for m=1:length(temp)
        fprintf('%d; ', temp(m));
    end
    fprintf('\n');
    obj.em_segmentation = input('which one to use: ');
else
    obj.em_segmentation = temp(1);
end
obj.matfile_em = sprintf('em_%d.mat', obj.em_segmentation);

%% import data
% obj.import_data();

%% update transformation info
obj.transformation = [];
obj.get_transformation();