function add_dataset(obj)
%% select data for running EASE
%{%}

%% inputs
%{
	obj: class object
%}

%% outputs
%{
    dataname: string; the selected data name
%}

%% Author
%{
	Pengcheng Zhou
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License
%}

%% existing datasets
datasets = obj.datasets_list;

% name of the new dataset
dataname = input('name of the datasets: ', 's');
if any(strcmpi(datasets, dataname))
    fprintf('The name ''%s'' is already used by another dataset. Please use a different name\n', dataname);
    return;
else
    datasets{end+1} = dataname;
end

% dabase information
database_name = input('database name (e.g., microns_ta3): ', 's'); 
datajoint_name = input('database alias (e.g. ta3): ', 's');
mesh_tablename = input('table of the EM meshes (e.g., ta3.Mesh): ', 's');

% save the added dataset
obj.datasets_list = datasets;
temp = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
temp.datasets_list = datasets;
temp.(dataname) = struct('datajoint_name', datajoint_name,...
    'database_name', database_name, ...
    'rel_mesh', mesh_tablename);
yaml.WriteYaml(fullfile(obj.dir_project, 'metainfo.yaml'), temp);

%％ create subfolders for the new file 
tmp_datafolder = fullfile(obj.dir_project, 'data', dataname); 
if ~exist(tmp_datafolder, 'dir')
    mkdir(tmp_datafolder); 
end
tmp_yamlpath = fullfile(obj.dir_project, sprintf('%s_config.yaml', dataname));
temp = yaml.ReadYaml(fullfile(fi.locate('ease', true), 'config.yaml')); 
temp.yaml_path = tmp_yamlpath; 

obj.yaml_path = fullfile(obj.dir_project, sprintf('%s_config.yaml', dataname));
if exist(obj.yaml_path, 'file')
    obj.read_config();
end

% assign folders corresponding this this data
obj.script_folder = fullfile(obj.dir_project, 'scripts', dataname);
obj.output_folder = fullfile(obj.dir_project, 'results', dataname);
obj.data_folder = fullfile(obj.dir_project, 'data', dataname);
obj.fig_folder = fullfile(obj.dir_project, 'Figures', dataname);
obj.video_folder = fullfile(obj.dir_project, 'Videos', dataname);

obj.write_config();

%% things to do 
fprintf('%s: dataset info added. Next we need the actual data for this dataset.\n', dataname); 
fprintf('\t1: registration.csv. A csv file for aligning EM space and the CI space.\n'); 
fprintf('\t2: calcium imaging videos.\n'); 
fprintf('\t3: 2p structural imaging.\n'); 
fprintf('You can find the detailed description of these datasets on xxx.\n');
fprintf('Once you get everything ready, run \n\t>>ease.import_data();\nto import them.\n'); 