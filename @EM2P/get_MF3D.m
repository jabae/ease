function neuron = get_MF3D(obj, T, create_new)
%% what does this function do
%{
    create a class object for storing the precessing results of the selected scan and block
%}

%% inputs:
%{
    T: number of frames
    create_new: boolean
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
% check the inputs
mscan = obj.scan_id;
mblock = obj.block_id;

if ~exist('T', 'var') || isempty(T)
    T = obj.video_T;
end
if ~exist('create_new', 'var')
    create_new = false;
end
if mblock==0
    var_name = sprintf('neuron_scan%d_all_blocks', mscan);
else
    var_name = sprintf('neuron_scan%d_block%d', mscan, mblock);
end

%% create a matfile for saving the results or load data from it directly
FOV_ = obj.FOV;
matfile_mf3d = fullfile(obj.output_folder, ...
    sprintf('neurons_%d_%d_%d_%d.mat',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));
folder_mf3d = fullfile(obj.output_folder, ...
    sprintf('neurons_%d_%d_%d_%d',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));
% a marker indicating where the data has been processed already
if ~exist(matfile_mf3d, 'file')
    flag_created = false(obj.num_scans, obj.num_blocks);
    save(matfile_mf3d, 'FOV_', 'flag_created','-v7.3');
else
    temp = load(matfile_mf3d, 'flag_created');
    flag_created = temp.flag_created;
end

%% load data directly
if ~create_new
    % existed already
    tmp_file = fullfile(folder_mf3d, sprintf('%s.mat', var_name)); 
    if exist(tmp_file, 'file')
        load(tmp_file, var_name); 
        neuron = eval(var_name);
        return;
    else
        fprintf('the desired MF3D wrapper was not created yet. create a new one.\n');
    end
end

%% load data directly, or create one and save it.
%% create MF3D objects
neuron = MF3D('d1', obj.d1, 'd2', obj.d2, 'd3', obj.d3,...
    'se', strel(ones(3,3,3)), 'search_method', 'dilate');
neuron.Fs = obj.video_Fs;
if mblock>0
    [dl_Yr, dl_Yd] = obj.create_dataloader(mscan, mblock, T);
    % dataloader for single blocks
    neuron.dataloader_denoised = dl_Yd;
    neuron.dataloader_raw = dl_Yr;
    if obj.use_denoise
        neuron.dataloader = neuron.dataloader_denoised;
    else
        neuron.dataloader = neuron.dataloader_raw;
    end
    %% save the results
    flag_created(mscan, mblock) = true; %#ok<NASGU>
    eval(sprintf('%s = neuron;', var_name));
    save(matfile_mf3d, var_name, 'flag_created', '-append');
else
    eval(sprintf('%s = neuron;', var_name));
    save(matfile_mf3d, var_name, '-append'); 
end


end