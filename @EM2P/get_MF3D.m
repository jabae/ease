function neuron = get_MF3D(obj, create_new)
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
        
        fprintf('load an existing MF3D wrapper for (scan %d, block %d).\n', mscan, mblock);
        
        load(tmp_file, var_name);
        neuron = eval(var_name);
        
        if isempty(neuron.P.sn)
            summary_images = obj.calculate_summary_images();
            neuron.P.sn = neuron.reshape(summary_images.sn, 1);
            
        end
        obj.summary_images = obj.calculate_summary_images();

        return;
    else
        fprintf('create a MF3D wrapper for (scan %d, block %d).\n', mscan, mblock);
    end
end

%% load data directly, or create one and save it.
%% create MF3D objects
neuron = MF3D('d1', obj.d1, 'd2', obj.d2, 'd3', obj.d3,...
    'se', strel(ones(3,3,3)), 'search_method', 'dilate');
neuron.Fs = obj.video_Fs;

% EM volume 
img = neuron.reshape(obj.em_volume(:, mscan), 3); 
img = imerode(img, strel('disk', 1)); 
neuron.spatial_range = neuron.reshape(img, 1);
neuron.options.normalize_data  = obj.normalize_data; 
neuron.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100, ...    % maximum decay time (unit: frame);
    'remove_large_residuals', true); % remove large residuals

obj.summary_images = obj.calculate_summary_images();
neuron.P.sn = neuron.reshape(obj.summary_images.sn, 1);

fprintf('done\n');
end