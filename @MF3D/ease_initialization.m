function ease_initialization(obj, Y, Aem, segment_ids, options, black_list, white_list)
%% initialize A & C in the CNMF model using spatial masks given by EM segments.
%{
    This function is used for initializing neurons given a pool of spatial
    masks. These masks for segmented from EM data.
%}

%% inputs:
%{
    Aem: d_em * K_em matrix, spatial footprints given by EM segments
    Y: d * T, data
    segment_ids: d*1, unique IDs for each EM segment
    options: struct variable for containing all options
    black_list: a list of EM indices that should be ignored during the
    initialization step
    white_list: a list of candidate EM components
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

%% check options for the initialization
options_default = struct('init_method', 'tf', ...%method for initializing single neurons
    'order_statistics', 'l3norm', ... % method for ordering neurons
    'min_similarity', 0.6,...   % threshold for rejecting bad initialization
    'clear_results', false,...  % remove previous result
    'save_fig', false, ...        % save figures for visualizing the initialization step
    'show_fig', true, ...        % show figures for visualizing the initializaito step
    'K_candidate', 3000,...     % number of neurons to be considered
    'K_new', 50, ...             % number of neurons to be added
    'min_pnr',3, ...             % the minimum peak-to-noise ratio for a good trace
    'preprocess_Y', true);       % preprocess Y or not

options = fix_missing_options(options_default, options); % fix missing values of the input options
init_method = options.init_method;
order_statistics = options.order_statistics;
min_similarity = options.min_similarity;
clear_results = options.clear_results;
save_fig = options.save_fig;
show_fig = options.show_fig;
K_candidate = options.K_candidate;
K_new = options.K_new;
min_pnr = options.min_pnr;
preprocess_Y = options.preprocess_Y;

%% pre-process Y
Y = obj.reshape(Y,1);   % reshape the data to a matrix

% preprocess results 
if obj.options.pre_process_data
    Y = obj.preprocess(Y); 
end

% clear existing results
if clear_results
    % delete all neurons
    if isempty(obj.A)
        obj.delete(1:length(obj.labels));
    end
    % delete the background
    obj.b = [];
    obj.f = [];
end

% roughly re-estimate the background
if isempty(obj.b)
    tmp_f = median(Y(obj.spatial_range,:), 1);
    tmp_f = tmp_f - mean(tmp_f);
    obj.b = (Y*tmp_f')/(tmp_f*tmp_f');
    obj.f = tmp_f;
end

% subtract background
if ~isempty(obj.b)  % keep
    Y = Y - obj.b*obj.f;
end

% subtract neural signal
if ~isempty(obj.A)
    Y = Y - obj.A*obj.C;
end

% subtract baseline
Y = bsxfun(@minus, Y, mean(Y, 2));

%% pre-allocate spaces for storing the results
[d, T] = size(Y);
if clear_results || isempty(obj.A)
    K_pre = 0;
    obj.A = zeros(d, K_new);
    obj.A_mask = zeros(d, K_new);
    obj.C = zeros(K_new, T);
    obj.C_raw = zeros(K_new, T);
    obj.S = zeros(K_new, T);
    obj.labels = zeros(K_new,1);
    obj.ids = zeros(K_new, 1);
    obj.match_status.status = zeros(1, K_new);
    obj.match_status.em_ids = cell(1, K_new);
    obj.match_status.confidence = zeros(1, K_new);
    obj.match_status.scores = zeros(1, K_new);
else
    K_pre = size(obj.A, 2);     % number of existing neurons
    obj.A = [obj.A, zeros(d, K_new)];
    obj.A_mask = [obj.A_mask, zeros(d, K_new)];
    obj.C = [obj.C; zeros(K_new, T)];
    obj.C_raw = [obj.C_raw; zeros(K_new, T)];
    obj.S = [obj.S; zeros(K_new, T)];
    obj.ids((end+1):(end+K_new)) = 0;
    obj.labels((end+1):(end+K_new)) = 0;
    obj.match_status.status((end+1):(end+K_new)) = 0;
    obj.match_status.confidence((end+1):(end+K_new)) = 0;
    obj.match_status.scores((end+1):(end+K_new)) = 0;
    obj.match_status.em_ids = [obj.match_status.em_ids, cell(1, K_new)];
end
k_new = K_pre;      % indices for the current cell ID
k_tried = 0;        % number of neurons tried for the initialization

%% pre-process Aem
if iscell(Aem)
    Aem =  cell2mat(Aem);
    segment_ids = cell2mat(segment_ids);
end

% only keep the top K_candidate components
K_total = size(Aem, 2);
if K_candidate < K_total
    temp = sum(Aem, 1);
    [~, idx] = sort(temp, 'descend');
    idx(1:K_candidate) = [];
    Aem(:, idx) = [];
    segment_ids(idx) = [];
end

% only select neurons from the whitelist
if exist('white_list', 'var') && ~isempty(white_list)
    idx = zeros(size(white_list));
    for m=1:length(white_list)
        temp = find(segment_ids==white_list(m), 1, 'first');
        if ~isempty(temp)
            idx(m) = temp;
        end
    end
    idx(idx==0) = [];
    Aem = Aem(:, idx);
    segment_ids = segment_ids(idx, 1);
    with_quality_control = false;
else
    with_quality_control = true;
end

K_total = length(segment_ids);
flag_use = true(K_total, 1);

% ignore the neurons in the black list
if exist('black_list', 'var') && ~isempty(black_list)
    for m=1:length(black_list)
        idx = find(segment_ids==black_list(m), 1, 'first');
        if ~isempty(idx)
            flag_use(idx) = false;
        end
    end
end

% ignore the neurons that have been found
if ~clear_results
    ind = (obj.match_status.status==1);
    tmp_list = cell2mat(obj.match_status.em_ids(ind));
    for m=1:length(tmp_list)
        idx = find(segment_ids==tmp_list(m), 1, 'first');
        if ~isempty(idx)
            flag_use(idx) = false;
        end
    end
end

% normalize Aem 
switch lower(order_statistics)
    case 'l3norm'
        norm_Aem = (sum(Aem.^3, 1)) ./ (sum(Aem.^3, 1).^(1.0/3));
    case 'ss3mad'
        norm_Aem = (sum(Aem.^3, 1)) ./ sqrt(sum(Aem.^2, 1));
    otherwise
        norm_Aem = (sum(Aem.^3, 1)) ./ sqrt(sum(Aem.^2, 1));
end
% remove segments that are empty
flag_use(isnan(norm_Aem)) = false;

% ignore the above EM segments
Aem = Aem(:, flag_use);
segment_ids = segment_ids(flag_use, 1);
norm_Aem = norm_Aem(1, flag_use);
Aem_l2norm = sqrt(sum(Aem.^2, 1)); 

%% compute the summay statistics for each neuron
batch_size = 1000;
tsub = 2;
K_total = size(Aem, 2);
nbatch = ceil(K_total/batch_size);
cc = zeros(1, size(Aem, 2));
Tmax = min(T, 10000);

% create a memmap file to store C
nbins = ceil(Tmax/tsub);
C = zeros(nbins, K_total);
C_fileID = fopen('temp_C.bin', 'w');
fwrite(C_fileID, C, 'double');
fclose(C_fileID);
clear C;
C_mmap = memmapfile('temp_C.bin', ...
    'Format', {'double' [nbins, K_total] 'C'}, ...
    'Writable', true);

for m=1:nbatch
    % select Aem
    if m==nbatch
        idx = ((m-1)*batch_size+1) : K_total;
    else
        idx = (m-1)*batch_size + (1:batch_size);
    end
    A_ = Aem(:, idx).^2;
    norm_ = reshape(norm_Aem(idx), [], 1);
    
    % compute temporal traces and their summary statistics
    C_ = bsxfun(@times, A_'*Y(:, 1:tsub:Tmax), 1./norm_);
    C_mmap.Data.C(:, idx) = C_';
    switch lower(order_statistics)
        case 'l3norm'
            % l3norm of each C
            cc(idx) = sum(C_.^3, 2);
        case 'ss3mad'
            % sum of squares of values above 3 mad
            C_med = median(C_, 2);
            C_mad = median(abs(bsxfun(@minus, C_, C_med)), 2);
            cc(idx) = sum((C_.^2).*bsxfun(@gt, C_, C_med+C_mad*3), 2);
        otherwise
                       cc(idx) = sum(C_.^3, 2);
    end
    clear C_;
end


ind_ignore = false(1, K_total);
ind_used = false(1, K_total);

%% run the initialization
% create a figure for exporting the initialization figures
if save_fig
    tmp_folder = evalin('base', 'ease.output_folder');
    output_folder = fullfile(tmp_folder, ['initialization_', get_date()]);
    mkdir(output_folder);
end

% canvas for ploting results
if show_fig
    figure('papersize', [10.08, 6.83]);
    init_fig;
    colormap jet;
    Nrows = 5;
end

deconv_options = obj.options.deconv_options;

%% greedy initialization
while (k_new < K_new+K_pre) && (k_tried < K_new*5)
    % find the best candidate components
    cc(ind_used) = 0;
    cc(ind_ignore) = 0;
    [cc_max, ind_max] = max(cc);
    if cc_max<=0
        obj.delete((k_new+1):size(obj.A,2));
        break;
    end
    k_tried = k_tried +1;
    
    % extract the corresponding EM segment ai and ci
    ai = obj.reshape(Aem(:, ind_max), 3);
    ci = reshape(ai, 1, []) * Y; %C_(ind_max, :);
    
    % visualize this EM segment
    if show_fig
        clf;
        % show ai and ci
        ai_max = max(ai(:));
        for m=1:3
            subplot(Nrows, 3, m);
            imagesc(ai(:, :, m), [0,ai_max/2]);
            axis equal off tight;
        end
        subplot(Nrows, 3, 4:6); cla;
        plot(ci, 'linewidth', 1);
        axis tight;
        set(gca, 'fontsize', 16);
    end
    
    %% refine (ai, ci)
    if strcmpi(init_method, 'semi_nmf')
        ci_new = deconvolveCa(ci, deconv_options, 'sn', mad(ci));
        ai = reshape(ai, [], 1);
        ind_mask = reshape(ai_mask>0, [], 1);
        for miter=1:2
            %% refine ai
            % threshold ci, replace it with a denoiser in the future
            if std(ci_new) ==0
                break;
            end
            ci_new = reshape(ci_new, 1, []) - mean(ci_new);
            ai_proj = (Y * ci_new')/(ci_new * ci_new'); %.*(ai(:)>0);
            ai_new = max(0, ai_proj.*ind_mask);
            
            %% refine ci
            ci_new = ((ai_new.*ai)' * Y) / ((ai_new.*ai)'*ai_new);
            ci_new_raw = ci_new;
            ci_new = deconvolveCa(ci_new, deconv_options, 'sn', mad(ci_new));
        end
    else % default method
        [ind_in, ind_out] = obj.construct_in_out(ai);
        Yin = Y(ind_in(:), :);
        Yout = Y(ind_out(:), :);
        
        % initialize one neuron
        ci0 = ci;
        ai_new = double(ind_in);
        [ai_new(ind_in), ci_new_raw] = initialize_ac_tf(Yin, Yout, ai(ind_in(:)), ci0, 15);
        
        % deconvolve/denoise the temporal trace
        [~, sn] = estimate_baseline_noise(ci_new_raw);
        ci_new = deconvolveCa(ci_new_raw, deconv_options, 'sn', sn);
        ci_new = reshape(ci_new, 1,[]);
        pnr = max(ci_new) / std(ci_new_raw - ci_new);
        
        %
        tmp_ci = ci_new - mean(ci_new);
        ai_proj = (Y * tmp_ci')/(tmp_ci*tmp_ci'); %.*(ai(:)>0);
        ai_proj = obj.reshape(ai_proj, 3);
    end
    
    % deal with silent activity cases
    if std(ci_new(2:end)) == 0
        continue;
    end
    
    %% show initialized components
    if show_fig
        % show new ai and ci
        ai_new_max = max(ai_new(:));
        ai_new = obj.reshape(ai_new, 3);
        ai_proj = obj.reshape(ai_proj, 3);
        for m=1:3
            subplot(Nrows, 3, m+6); cla;
            imagesc(ai_proj(:, :, m), [0, ai_new_max/2]);
            axis equal off tight;
            
            subplot(Nrows, 3, m+9); cla;
            imagesc(ai_new(:, :, m), [0, ai_new_max/2]);
            axis equal off tight;
        end
        subplot(Nrows, 3, 4:6);
        cla;
        plot(ci_new_raw); hold on;
        plot(ci_new, 'r', 'linewidth', 1);
        axis tight;
    end
    
    %% quality control
    if with_quality_control
        % check the temporal trace
        if pnr < min_pnr
            continue;
        end
        % check the spatial footprint
        if sum(ai_new(:).*ai(:))/norm(ai_new(ai>0), 2)/norm(ai(:),2) < min_similarity
            ind_ignore(ind_max) = true;
            cc(ind_max) = 0;
            fprintf('hmmm, try next\n');
            if save_fig
                saveas(gcf, fullfile(output_folder, sprintf('delete_%d.pdf', k_tried)));
            end
            
            % find the next candidate neuron
            tmp_ai = ai_proj; 
            tmp_ai(tmp_ai<0) = 0; 
            temp = (reshape(tmp_ai,1, [])*Aem) ./ Aem_l2norm / norm(tmp_ai(:), 2); 
            [tmp_max, tmp_ind] = max(temp);
            if tmp_max>=min_similarity
                cc(tmp_ind) = inf;  % find the next candidate neuron
            end
            continue;
        end
    end
    
    %% add the neurons
    k_new  = k_new + 1;
    obj.C_raw(k_new,:) = ci_new_raw(:);
    obj.A(:, k_new) = ai_new(:);
    obj.C(k_new, :) = ci_new(:);
    obj.A_mask(:, k_new) = ai(:);
    obj.ids(k_new) = segment_ids(ind_max);
    obj.match_status.status(k_new) = 1;
    obj.match_status.confidence(k_new) = 0;
    obj.match_status.scores(k_new) = 0;
    obj.match_status.em_ids{k_new} = segment_ids(ind_max);
    ind_used(ind_max) = true;
    fprintf('\t %d\n', k_new);
    
    %% subtract this components
    ci_new = reshape(ci_new, 1, []);
    ci_new = ci_new - mean(ci_new);
    Y = Y - reshape(ai_new, [], 1) * ci_new;
    if show_fig
        ai_res = ai_proj - ai_new;
        ai_res = obj.reshape(ai_res, 3);
        for m=1:3
            subplot(Nrows, 3, 12+m); cla;
            imagesc(ai_res(:, :, m), [0, ai_new_max/3]);
            axis equal off tight;
        end
    end
    drawnow();
    
    if save_fig
        saveas(gcf, fullfile(output_folder, sprintf('%d.pdf', k_new)));
    end
    
    %% update cc
    temp = (Aem.^2)' * reshape(ai_new, [], 1);
    idx = (temp>0);
    C_ = C_mmap.Data.C(:, idx);
    C_ = C_' - temp(idx)./norm_Aem(idx)'*ci_new(1:tsub:Tmax);
    C_mmap.Data.C(:, idx) = C_';
    switch lower(order_statistics)
        case 'l3norm'
            % l3norm of each C
            cc(idx) = sum(C_.^3, 2);
        case 'ss3mad'
            % sum of squares of values above 3 mad
            C_med = median(C_, 2);
            C_mad = median(abs(bsxfun(@minus, C_, C_med)), 2);
            cc(idx) = sum((C_.^2).*bsxfun(@gt, C_, C_med+C_mad*3), 2);
        otherwise
            cc(idx) = sum(C_.^3,2); 
    end
    clear C_;
end

delete('temp_C.bin'); 
obj.delete(find(obj.match_status.status==0)); 
