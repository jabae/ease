function flag = read_config(obj, path_file)
%% what does this function do
%{
    read configurations for running EASE from a YAML file
%}

%% inputs:
%{
    path_file: the path of the yaml file. if it's empty, they use the
    default yaml file obj.yaml_path
%}

%% outputs:
%{
    flag: boolean variable, success or no
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code

% check the input arguments
if ~exist('path_file', 'var') || isempty(path_file) || ~exist(path_file, 'file')
    path_file = obj.yaml_path;
end

% load yaml file
configs = yaml.ReadYaml(path_file);

% pass configs to the class object
try
    temp = fieldnames(configs);
    for m=1:numel(temp)
        eval(sprintf('obj.%s=configs.%s;', temp{m}, temp{m}));
    end
    
    flag = true;
    fprintf('\nThe configuration of the EASE environment has been updated from \n%s\n\n', ...
        path_file); 
    obj.yaml_path = path_file; 
    
    % update few variables to correct the format 
    nams = {'dims_stack', ...
        'dims_video', ...
        'range_2p', ...
        'video_shifts.ii', ...
        'video_shifts.jj', ...
        'video_zvals', ...
        'video_zvals_updated', ...
        'FOV', ...
        'FOV_stack', ...
        'range_2p', ...
        'em_shifts.ii', ...
        'em_shifts.jj', ...
        'stack_shifts.ii', ...
        'stack_shifts.jj' };
    for m=1:length(nams)
        try
            eval(sprintf('obj.%s=cell2mat(obj.%s); ', nams{m}, nams{m}));
        end
    end
catch
    flag = false;
    fprintf('\bLoading configurations has been failed when passing variable %s\n\n', temp{m});
end