function add_dataset(obj)
%% select data for running EASE
%{%}

%% inputs
%{
	obj: class object
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

% create subfolders for the new file 
tmp_datafolder = fullfile(obj.dir_project, 'data', dataname); 
mkdir(tmp_datafolder); 
tmp_yamlpath = fullfile(obj.dir_project, sprintf('%s_config.yaml', dataname)); 
temp = yaml.ReadYaml(fullfile(fi.locate('ease', true), 'config.yaml')); 
temp.yaml_path = tmp_yamlpath; 
yaml.WriteYaml(tmp_yamlpath, temp); 

%% things to do 
fprintf('%s: dataset added. Here are things you need to do: \n', dataname); 
fprintf('\t1. add data files to folder: %s\n',  tmp_datafolder); 
fprintf('\t2. modify data options: %s\n', fullfile(obj.dir_project, sprintf('%s_config.yaml', dataname)));
