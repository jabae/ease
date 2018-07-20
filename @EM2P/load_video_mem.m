function tmpY = load_video_mem(obj, mscan, mblock)
%% what does this function do 
%{
    load calcium imaging videos into the memory for faster access 
%}

%% inputs: 
%{
    mscan: scan ID 
    mblock: block ID 
%}

%% outputs: 
%{
    tmpY: data corresponding to the selected scan ID and block ID 
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
if isempty(tmpY)
    neuron = obj.get_MF3D(obj.video_T);
    if obj.use_denoise
        dl = neuron.dataloader_denoised;
    else
        dl = neuron.dataloader_raw;
    end
    assignin('base', 'tmpY', dl.load_tzrc());
    evalin('base', sprintf('%s{%d,%d}=tmpY;', var_name, mscan, mblock));
else
    evalin('base', sprintf('tmpY=%s{%d,%d};', var_name, mscan, mblock)); 
end
end