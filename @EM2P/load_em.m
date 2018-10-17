function load_em(obj)
%% what does this function do 
%{
    load EM data and saved the results to obj 
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
%% get the matfile information
em_mat_file = fullfile(obj.data_folder, obj.matfile_em); 
if isempty(obj.em_data)
    % check the existance of the EM matfile
    if ~exist(em_mat_file, 'file')
        error('EM meshes have not been voxelized yet. You need stop everything and run ease_voxelize_em.m'); 
        edit('ease_voxelize_em.m'); 
    end
    
    % map the data to memory
    obj.em_data = matfile(em_mat_file, 'Writable', false);
    obj.em_variables = who(obj.em_data, 'slice_*');
    obj.em_info = obj.em_data.EM_info;
    obj.em_ranges = obj.em_data.em_ranges;
    obj.K_em = size(obj.em_info, 1); 
else
    fprintf('\nEM data has been mapped to memory already\n\n');
end


%% load data into memory for speed concern
if obj.em_load_flag
    if exist_in_workspace('em_data_mem', 'base')
        % loaded already
        fprintf('\nEM data has been loaded to memory already\n\n');
    else
        fprintf('loading the EM data into memory...\n');
        evalin('base', sprintf('em_data_mem=load(''%s'',''-regexp'',''^(?!cell_).*'');', ...
            em_mat_file));
        fprintf('Done!\n\n');
    end
else
    fprintf(['\nTips: The EM data are accessed from the hard drive. \n', ...
        'This requires less memory, but costs more time for \n', ...
        'fetching data. So if your memory is large enough, \n', ...
        'you can set ease.em_load_flag as true to preload the data\n\n']);
end
end
