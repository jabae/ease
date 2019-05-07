function del_dataset(obj)
%% select data for running EASE
%{%}

%% inputs
%{
	obj: class object
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

%% choose dataset 
metainfo = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
datasets = metainfo.datasets_list; 

if isempty(datasets) || (length(datasets)== 1 && strcmpi(datasets, 'example_data'))
    edit(fullfile(obj.dir_project, 'metainfo.yaml'));
    error('no datasets available for this project.');
end
fprintf('\n**********choose the data to use**********\n');
for m=1:length(datasets)
    fprintf('%d: %s\n', m, datasets{m});
end
fprintf('0: cancel\n'); 
fprintf('********************************************\n');

data_id = input('data ID: ');
while true
    if any(data_id==(1:length(datasets)))
        data_name = datasets{data_id};
        fprintf('you are going to delete data: %s \n', data_name);
        temp = input('are you sure? (y/n): ', 's'); 
        if strcmpi(temp, 'y')
            break;
        end
    elseif data_id==0
        fprintf('no dataset will be deleted.\n'); 
        return; 
    else
        data_id = input('please type a valid data ID: ');
    end
end

%% load data info 
metainfo.datasets_list(data_id) = [];
if isfield(metainfo, data_name)
    metainfo = rmfield(metainfo, data_name);
end
yaml.WriteYaml(fullfile(obj.dir_project, 'metainfo.yaml'), metainfo);
obj.datasets_list = metainfo.datasets_list; 

%% delete the corresponding subfolders for the new file 
tmp_yamlpath = fullfile(obj.dir_project, sprintf('%s_config.yaml', data_name)); 
if exist(tmp_yamlpath, 'file')
    delete(tmp_yamlpath); 
end 

temp = {'data', 'scripts', 'results', 'Figures', 'Videos'};
for m=1:length(temp)
    tmp_folder = fullfile(obj.dir_project, temp{m}, data_name);
    if exist(tmp_folder, 'dir')
        rmdir(tmp_folder, 's');
    end
end