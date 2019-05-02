function at_ease(obj, configs)
%% run EASE to automatically initialize neurons & update model variables
%{
	your computer is running EASE, so you are at ease. it's coffee time.
%}

%% inputs
%{
	obj: MF3D class object
    configs: struct variable; configurations for running EASE
%}

%% outputs
%{
%}

%% Author
%{
	Pengcheng Zhou
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	XXX License
%}

%% code
fprintf('\n******************************************************************\n');
fprintf('yay, it''s coffee time!\nThe computer is working very hard for you.\n');
fprintf('******************************************************************\n');

%% load data
% video
fprintf('preparing the video data...\n');
evalin('base', 'Y=ease.load_Y();');
Y = evalin('base', 'Y');
Y = obj.preprocess(Y);
obj.options.pre_process_data = false;
fprintf('done.\n');

% EM
fprintf('loading the EM info...\n');

scan_id = evalin('base', 'ease.scan_id');
try
    current_scan_id_for_em = evalin('base', 'current_scan_id_for_em');
catch
    current_scan_id_for_em = 0;
end
if current_scan_id_for_em ~=scan_id
    evalin('base', '[Aem, segment_ids] = ease.load_Aem();');
end
Aem = evalin('base', 'Aem');
segment_ids = evalin('base', 'segment_ids');
fprintf('done.\n');

%% configurations
n_iters = length(configs);

for miter=1:n_iters
    tmp_config = configs{miter};
    
    % initialize K_new neurons from the top K components
    options_init = tmp_config.options_init;
    if miter==2
        options_init.clear_results = false;
    elseif options_init.clear_results && ~isempty(obj.A)
        temp = input('Do you want to clear all existing results? (y/n)  ', 's');
        if ~strcmpi(temp, 'y')
            options_init.clear_results = false;
        end
    end
    
    % initialize neurons
    obj.ease_initialization(Y, Aem, segment_ids, options_init, ...
        tmp_config.black_list, tmp_config.white_list);
    obj.options.nb = tmp_config.nb;
    
    % run HALS to update
    fprintf('updating model variables\n'); 
    obj.hals(Y, [], [], false);
    fprintf('done\n'); 
    % evaluate the matching performance and correct the bad matches
    if ~isnan(tmp_config.min_rank)
        obj.evaluate_matching_confidence(Y, Aem, segment_ids);
        
        n_correction = obj.rematch(Aem, segment_ids, tmp_config.min_rank);
        if n_correction
            % remove duplicate matches
            obj.merge_repeats();
            % run hals to update model variables again
            obj.hals(Y, [], [], false);
        end
    end
end
obj.evaluate_matching_confidence(Y, Aem, segment_ids);
obj.options.pre_process_data = true;
