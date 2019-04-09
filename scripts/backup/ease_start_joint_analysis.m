%% select scan IDs for joint analysis 
ease.block_id = 0; 
options_joint.block_id = ease.block_id; 
options_joint.scan_ids = input('scan IDs for joint analysis (e.g., [1, 2, 3, 4]): '); 
options_joint.nscan = length(options_joint.scan_ids); 

%% deconvolution options 
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100, ...    % maximum decay time (unit: frame);
    'remove_large_residuals', true); % remove large residuals

%% load the raw data for all these scans 
neurons_all = cell(options_joint.nscan, 1);  % wrapper for each scan
Y_all = cell(options_joint.nscan, 1);   % video data 
sn_all = cell(options_joint.nscan, 1);  % noise level 
indz_all = cell(options_joint.nscan, 1); 

for m=1:options_joint.nscan
    ease.scan_id = options_joint.scan_ids(m);
    neuron = ease.get_MF3D();  % create a wrapper for one scan 
     
    Y = neuron.reshape(ease.get_Y(options_joint.scan_ids(m)), 1);  % load data
    
    if isempty(neuron.P.sn)  % compute noise level 
        T = min(size(Y, 2), 5000);
        Y_sn = GetSn(bsxfun(@minus, double(Y(:, (end-T+1):end)), mean(Y(:, (end-T+1):end), 2)));
        neuron.P.sn = Y_sn;
    end 
    
    
    neuron.options.deconv_flag = true; 
    neuron.options.deconv_options = deconv_options; 
    neurons_all{m} = neuron; 
    Y_all{m} = Y;  % normalize data 
    
    ind = false(ease.d1, ease.d2, ease.num_scans*ease.num_slices); 
    ind(:, :, (1:ease.num_slices) + (ease.scan_id-1)*ease.num_slices) = true; 
    indz_all{m} = ind; 
end 

%% initialize neurons in different scans independently 
segment_ids = fetchn(ta3p100.AllenSomaClass &...
    'cell_class=''excitatory''', 'segment_id'); 
[Aem, segment_ids] = ease.get_em_footprints(segment_ids, options_joint.scan_ids); 
options_init = ease.options_init; 
options_init.clear_results = true; 
for m=1:options_joint.nscan 
    tmpAem = Aem{m}; 
    neuron = neurons_all{m}; 
    neuron.initialize_multiple(tmpAem, segment_ids, ...
        Y_all{m}(:, 201:end), options_init); 
    neurons_all{m} = neuron; 
end 

%% 
em_id = uint64(648518346349490624); 
Aem = ease.get_em_footprints(em_id, options_joint.scan_ids);
for m=1:options_joint.nscan
    neuron = neurons_all{m};
    ai_em = Aem{m}; 
    [ind_in, ind_out] = neuron.construct_in_out(ai_em);
    temp = struct('ind_in', ind_in, 'ind_out', ind_out); 
    temp.Yin = Y_all{m}(ind_in, 201:end); 
    temp.Yout = Y_all{m}(ind_out, 201:end); 

    [ai, ci, si, ci_raw] = neuron.initialize_one(ai_em, temp);
end
