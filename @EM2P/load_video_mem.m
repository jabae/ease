function neuron = load_video_mem(obj, mscan, mblock, create_new)
%% what does this function do
%{
    load calcium imaging videos into the memory for faster access
%}

%% inputs:
%{
    mscan: scan ID
    mblock: block ID
    create_new: boolean, create a new neuron object instead of loading it
    from the saved results
%}

%% outputs:
%{
    neuron: a class object of @MF3D. It stores information of one imaging
    scan

    in the base wokspace:
    Y_in_use: 4D array, the video data of the selected scan
    Y_denoised{mscan, mblock} / Y_denoised{mscan, mblock} = Y_in_use;
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code
% check the inputs
if ~exist('mscan', 'var') || isempty(mscan)
    mscan = obj.scan_id;
else
    obj.scan_id = mscan;
end
if ~exist('mblock', 'var') || isempty(mscan)
    mblock = obj.block_id;
else
    obj.block_id = mblock;
end
if ~exist('create_new', 'var') || isempty(create_new)
    create_new = false;
end
if obj.use_denoise
    var_name = 'Y_denoised';
else
    var_name = 'Y_raw';
end

% check the existance in the base workspace
if ~exist_in_workspace(var_name, 'base')
    assignin('base', var_name, cell(obj.num_scans, obj.num_blocks));
end

% load data
tmpY = evalin('base', sprintf('%s{%d,%d}', var_name, mscan, mblock));
neuron = obj.get_MF3D(obj.video_T, create_new);
if obj.use_denoise
    dl = neuron.dataloader_denoised;
else
    dl = neuron.dataloader_raw;
end
if isempty(tmpY)
    assignin('base', 'Y_in_use', dl.load_tzrc());
    evalin('base', sprintf('%s{%d,%d}=Y_in_use;', var_name, mscan, mblock));
else
    evalin('base', sprintf('Y_in_use=%s{%d,%d};', var_name, mscan, mblock));
end
evalin('base', sprintf('cell_id=1;'));
end