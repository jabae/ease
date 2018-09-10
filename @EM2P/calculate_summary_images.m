function summary_images = calculate_summary_images(obj, mscan, mblock)
%% what does this function do
%{
    create a class object for storing the precessing results of the selected scan and block
%}

%% inputs:
%{
    mscan: scan ID
    mblock: block ID
%}

%% outputs:
%{
    summary_images: struct variable with all summary images for the
    corresponding data.
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code
% check the inputs
if ~exist('mscan', 'var')
    mscan = obj.scan_id;
end
if ~exist('mblock', 'var')
    mblock = obj.block_id;
end
% create a mat file for storing results
FOV_ = obj.FOV;
matfile_summary = fullfile(obj.output_folder, ...
    sprintf('summary_images_%d_%d_%d_%d.mat',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));
if ~exist(matfile_summary, 'file')
    flag_processed_raw = false(obj.num_scans, obj.num_blocks);
    flag_processed_denoised = flag_processed_raw; %#ok<NASGU>
    flag_processed = false(obj.num_scans, obj.num_blocks);
    save(matfile_summary, 'FOV_', 'flag_processed_raw',...
        'flag_processed_denoised', '-v7.3');
else
    if obj.use_denoise
        temp = load(matfile_summary, 'flag_processed_denoised');
        flag_processed = temp.flag_processed_denoised;
    else
        temp = load(matfile_summary, 'flag_processed_raw');
        flag_processed = temp.flag_processed_raw;
    end
end

%% when we process all blocks, just randomly select one block to get its results
if mblock==0
    mblock = find(flag_processed(mscan, :), 1, 'first');
    if isempty(mblock)
        mblock = 1;
    end
end

% check whether the dada has been processed or not
if obj.use_denoise
    var_name = sprintf('scan%d_block%d_denoised', mscan, mblock);
    data_name = 'Y_denoised';
else
    var_name = sprintf('scan%d_block%d_raw', mscan, mblock);
    data_name = 'Y_raw';
end
fprintf('video data: scan=%d, block=%d\n', mscan, mblock);

if flag_processed(mscan, mblock)
    temp = load(matfile_summary, var_name);  %#ok<NASGU>
    summary_images = eval(sprintf('temp.%s', var_name));
    fprintf('    The summary images for the selected data were loaded.\n\n');
    return;
end

% data loader
fprintf('\n    Loading the data of the selected scan.\n');
if exist_in_workspace(data_name, 'base')
    Y = evalin('base', sprintf('%s{%d,%d}', data_name,mscan, mblock));
else
    Y = [];
end

if isempty(Y)
    neuron = obj.get_MF3D(mscan, mblock);
    if obj.use_denoise
        dl = neuron.dataloader_denoised;
    else
        dl = neuron.dataloader_raw;
    end
    Y = dl.load_tzrc();
end
fprintf('    Done\n');
%% computing summary statistics
summary_images = struct('std', zeros(obj.d1, obj.d2, obj.num_slices), ...
    'cn', zeros(obj.d1, obj.d2, obj.num_slices), ...
    'max', zeros(obj.d1, obj.d2, obj.num_slices), ...
    'mean', zeros(obj.d1, obj.d2, obj.num_slices), ...
    'sn', zeros(obj.d1, obj.d2, obj.num_slices), ...
    'pnr', zeros(obj.d1, obj.d2, obj.num_slices));
for mslice=1:obj.num_slices
    fprintf('\nSlice %d\n', mslice);
    tmpY = squeeze(Y(:, :, mslice,:));
    
    fprintf('\n    Calculating summary images: ');
    % standard deviation
    fprintf('\n\t standard deviation');
    summary_images.std(:, :, mslice) = std(double(tmpY), 0, 3);
    
    % correlation image
    fprintf('\n\t correlation image');
    summary_images.cn(:, :, mslice) = correlation_image(tmpY);
    
    % maximum projection
    fprintf('\n\t maximum projection');
    summary_images.max(:, :, mslice)  = double(max(tmpY, [], 3));
    
    % mean projection
    fprintf('\n\t mean projection');
    summary_images.mean(:, :, mslice)  = mean(tmpY, 3);
    
    % noise level
    temp = reshape(tmpY, obj.d1 * obj.d2, []);
    fprintf('\n\t noise levels');
    summary_images.sn(:, :, mslice)  = reshape(GetSn(temp), obj.d1, obj.d2);
    fprintf('\n\t peak-to-noise ratio\n');
    summary_images.pnr(:, :, mslice)  = (summary_images.max(:, :, mslice) -...
        summary_images.mean(:, :, mslice)) ./ summary_images.sn(:, :, mslice);
end

%% save the results
flag_processed(mscan, mblock) = true;
eval(sprintf('%s=summary_images;', var_name));
if obj.use_denoise
    flag_processed_denoised = flag_processed;  %#ok<NASGU>
    save(matfile_summary, var_name, 'flag_processed_denoised', '-append');
else
    flag_processed_raw = flag_processed;  %#ok<NASGU>
    save(matfile_summary, var_name, 'flag_processed_raw', '-append');
end

end