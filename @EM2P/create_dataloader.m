function [dl_Yr, dl_Yd] = create_dataloader(obj, mscan, mblock, T)
if ~exist('mscan', 'var') || isempty(mscan)
    mscan = obj.scan_id; 
end 
if ~exist('mblock', 'var') || isempty(mblock)
    mblock = obj.block_id; 
end
if ~exist('T', 'var') || isempty(T)
    T = obj.video_T;
end

% data loader for denoised data
vars_denoised = {fullfile(obj.data_folder, obj.denoised_folder, ...
    sprintf('scan%d_block%d_gpca_p',mscan, mblock)),...
    '_data.npy'};
fname_denoised = @(vars, z) [vars{1}, int2str(z), vars{2}];
dl_Yd = IDL('vars', vars_denoised, 'type', 'npy', 'fname', fname_denoised, ...
    'dims', [obj.d1, obj.d2, obj.d3], 'num_frames', T);

% data loader for raw data
vars_raw = {fullfile(obj.data_folder, obj.raw_folder, sprintf('scan%d_block%d_complete.mat',...
    mscan, mblock))};
fname_raw = @(vars, z) vars{1};
dl_Yr = IDL('vars', vars_raw, 'type', 'mat4d', 'fname', fname_raw, 'dims', ...
    [obj.d1, obj.d2, obj.d3], 'num_frames', T, 'nfiles', false, 'shifts', ...
    [0, 0, 0, 0]);
end