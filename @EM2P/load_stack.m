function load_stack(obj)
%% what does this function do
%{
    load 2p stack data into the working space
%}

%% inputs:
%{
%}

%% outputs:
%{
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code
if exist_in_workspace('stack_2p', 'base')
    fprintf('The 2P stack has been loaded already.\n');
    return;
end

stack_mat_file1 = fullfile(obj.data_folder, obj.matfile_stack);
stack_mat_file2 = fullfile(obj.data_folder, 'dl_stack.mat');

if exist(stack_mat_file1, 'file')
    stack_data = matfile(stack_mat_file1, 'Writable', false);
    fprintf('\nloading the 2p stack data...\n');
    assignin('base', 'stack_2p', stack_data.stack_2p);
    fprintf('Done!\n\n');
    return;
elseif exist(stack_mat_file2, 'file')
    load(stack_mat_file2, 'dl_stack');
    temp = dl_stack.load_tzrc(); 
    assignin('base', 'stack_2p', temp(:, :, end:-1:1));
else
    obj.import_data();
    obj.load_stack();
end