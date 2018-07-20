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
    T = ease.video_T;
end
if ~exist('create_new', 'var')
    create_new = false; 
end
var_name = sprintf('neuron_scan%d_block%d', mscan, mblock);

%% create a matfile for saving the results or load data from it directly
FOV_ = obj.FOV;
matfile_mf3d = fullfile(obj.output_folder, ...
    sprintf('neurons_%d_%d_%d_%d.mat',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));

% a marker indicating where the data has been processed already
if ~exist(matfile_mf3d, 'file')
    flag_created = false(obj.num_scans, obj.num_blocks);
    save(matfile_mf3d, 'FOV_', 'flag_created','-v7.3');
else
    temp = load(matfile_mf3d, 'flag_created');
    flag_created = temp.flag_created;
end

%% load data directly, or create one and save it.
if flag_created(mscan, mblock) && ~create_new
    % existed already
    temp = load(matfile_mf3d, var_name);  %#ok<NASGU>
    neuron = eval(sprintf('temp.%s', var_name));
    return;
else
    %% create MF3D objects
    neuron = MF3D('d1', obj.d1, 'd2', obj.d2, 'd3', obj.d3,...
        'se', strel(ones(3,3,3)), 'search_method', 'dilate');
    neuron.Fs = obj.video_Fs; 
    [dl_Yr, dl_Yd] = obj.create_dataloader(mscan, mblock, T);
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
end
end