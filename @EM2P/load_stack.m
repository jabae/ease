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

stack_data = matfile(fullfile(obj.data_folder, obj.matfile_stack), 'Writable', false);

if exist_in_workspace('stack_2p', 'base')
    % exist, skip it 
    fprintf('The 2P stack has been loaded already.\n');
else
    % not exist, load it. 
    fprintf('\nloading the 2p stack data...\n');
    assignin('base', 'stack_2p', stack_data.stack_2p);
    fprintf('Done!\n\n');
end
end