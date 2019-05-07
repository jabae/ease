function add_database(obj, databases_list)
%% add a database link
%{%}

%% inputs
%{
	obj: class object
%}

%% outputs
%{
    databases_list: a cell of strings or one string; the selected data name
%}

%% Author
%{
	Pengcheng Zhou
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License
%}

%% existing datasets
if ~exist('databases_list', 'var') || isempty(databases_list)
    databases_list = input('database link (e.g., 127.0.0.1:3306): ', 's'); 
end 

if ischar(databases_list) 
    databases_list = {databases_list}; 
end 

% save the added dataset
temp = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
obj.databases_list = union(temp.databases_list, databases_list); 
temp.databases_list = obj.databases_list; 
yaml.WriteYaml(fullfile(obj.dir_project, 'metainfo.yaml'), temp);
