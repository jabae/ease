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
    black_list: a list of EM indices that should be ignored 
    white_list: a list of candidate EM components
%}

%% outputs:
%{
    all results are saved to the obj class object. 
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
    'use_wrss', true, ...    % use wrss or not
    'min_similarity', 0.6,...   % threshold for rejecting bad initialization
    'clear_results', false,...  % remove previous result
    'save_fig', false, ...        % save figures for visualizing the initialization step
    'show_fig', true, ...        % show figures for visualizing the initializaito step
    'K_candidate', 3000,...     % number of neurons to be considered
    'K_new', 50, ...             % number of neurons to be added
    'min_pnr',3);                % the minimum peak-to-noise ratio for a good trace

options = fix_missing_options(options_default, options); % fix missing values of the input options
init_method = options.init_method;
order_statistics = options.order_statistics;
use_wrss = options.use_wrss;
min_similarity = options.min_similarity;
clear_results = options.clear_results;
save_fig = options.save_fig;
show_fig = options.show_fig;
K_candidate = options.K_candidate;
K_new = options.K_new;
min_pnr = options.min_pnr;

%% pre-process Y
Y = obj.reshape(Y,1);   % reshape the data to a matrix
% preprocess data
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
if use_wrss
    temp1 = sum(Aem.^3,1);
else
    temp1 = sum(Aem.^2,1);
end
switch lower(order_statistics)
    case 'l3norm'
        temp2 = sum(Aem.^3,1).^(1.0/3);
    case 'ss3mad'
        temp2 = sqrt(sum(Aem.^2));
    case 'l3norm_xcorr'
        temp2 = sum(Aem.^3,1).^(1.0/3);
    otherwise
        error('the selected order statistics - %s is not supported yet', order_statistics);
end
norm_Aem = temp1./temp2;

% remove segments that are empty
flag_use(isnan(norm_Aem)) = false;
Aem = Aem(:, flag_use);
segment_ids = segment_ids(flag_use, 1);
norm_Aem = norm_Aem(1, flag_use);
Aem_l2norm = sqrt(sum(Aem.^2, 1));
K_new = min(K_new, size(Aem, 2)); 
if isempty(Aem)
    fprintf('no candidate EM components. delete the pre-specified space.\n');
    obj.delete(obj.match_status.confidence==0);
    return;
end

%% pre-allocate spaces for storing the results
[d, T] = size(Y);
if clear_results || isempty(obj.A)
    K_pre = 0;
    obj.A = zeros(d, K_new);
    obj.A_em = zeros(d, K_new);
    obj.C = zeros(K_new, T);
    obj.C_raw = zeros(K_new, T);
    obj.S = zeros(K_new, T);
    obj.labels = zeros(K_new,1);
    obj.ids = uint64(zeros(K_new, 1));
    obj.match_status.status = zeros(1, K_new);
    obj.match_status.em_ids = cell(1, K_new);
    obj.match_status.confidence = zeros(1, K_new);
    obj.match_status.scores = zeros(1, K_new);
else
    K_pre = size(obj.A, 2);     % number of existing neurons
    obj.A = [obj.A, zeros(d, K_new)];
    obj.A_em = [obj.A_em, zeros(d, K_new)];
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
    if use_wrss
        A_ = Aem(:, idx).^2;
    else
        A_ = Aem(:, idx);
    end
    norm_ = reshape(norm_Aem(idx), [], 1);
    
    % compute temporal traces and their summary statistics
    C_ = bsxfun(@times, A_'*Y(:, 1:tsub:Tmax), 1./norm_);
    switch lower(order_statistics)
        case 'l3norm'
            % l3norm of each C
            cc(idx) = sum(C_.^3, 2);
        case 'ss3mad'
            % sum of squares of values above 3 mad
            C_med = median(C_, 2);
            C_mad = median(abs(bsxfun(@minus, C_, C_med)), 2);
            cc(idx) = sum((C_.^2).*bsxfun(@gt, C_, C_med+C_mad*3), 2);
        case 'l3norm_xcorr'
            cc(idx) = sum(C_.^3, 2) .* (sum(C_(:,1:(end-1)).*C_(:, 2:end), 2)./sum(C_(:,2:end).^2,2)); 
    end
    % save the result
    C_mmap.Data.C(:, idx) = C_';
    clear C_;
end

ind_ignore = false(1, K_total);
ind_used = false(1, K_total);

%% run greedy initialization
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

while (k_new < K_new+K_pre) && (k_tried < K_new*10)
    % find the best candidate components
    cc(ind_used) = -inf;
    cc(ind_ignore) = -inf;
    [cc_max, ind_max] = max(cc);
    if cc_max==(-inf)
        obj.delete((k_new+1):size(obj.A,2));
        break;
    end
    k_tried = k_tried +1;
    
    % extract the corresponding EM segment ai and ci
    ai = obj.reshape(Aem(:, ind_max), 3);
    %     ai = imfilter(ai, fspecial('gaussian', 3, 1));
    ci = reshape(ai, 1, []) * Y; %C_(ind_max, :);
    temp = imfilter(obj.reshape(double(ai>max(ai(:))*0.1), 3), ones(10));
    if max(temp(:))<50
        is_soma = false; 
    else
        is_soma = true; 
    end
    
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
        [ind_in, ~] = obj.construct_in_out(ai);
        Yin = Y(ind_in(:), :);
        
        % initialize one neuron
        ci0 = ci;
        ai_new = double(ind_in);
        [ai_new(ind_in), ci_new_raw] = initialize_ac_seminmf(Yin, ai(ind_in(:)), ci0, 20);
        
        % deconvolve/denoise the temporal trace
        [~, sn] = estimate_baseline_noise(ci_new_raw);
        ci_new = deconvolveCa(ci_new_raw, deconv_options, 'sn', sn);
        ci_new = reshape(ci_new, 1,[]);
        pnr = max(ci_new) / std(ci_new_raw - ci_new);
        
        %
        tmp_ci = ci_new - mean(ci_new);
        ai_proj = (Y * tmp_ci')/(tmp_ci*tmp_ci'); %.*(ai(:)>0);
        ai_proj = obj.reshape(ai_proj, 3);
    else % default method 
        [ind_in, ind_out] = obj.construct_in_out(ai);
        Yin = Y(ind_in(:), :);
        Yout = Y(ind_out(:), :);
        
        % initialize one neuron
        ci0 = ci;
        ai_new = double(ind_in);
        [ai_new(ind_in), ci_new_raw] = initialize_ac_tf(Yin, Yout, ai(ind_in(:)), ci0, 20);
        
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
    
    if with_quality_control && (~is_soma)
        % delete traces with low PNR
        if pnr < min_pnr
            continue;
        end
        
        % check the spatial footprint
        tmp_ai = ai_proj;
        tmp_ai(tmp_ai<0) = 0;
        tmp_img = image_local_corr(tmp_ai, ai, 7);
        tmp_corr = corr(tmp_img(:).*(ai(:)>0), ai(:));
        %         tmp_corr = corr(tmp_img(:).*tmp_ai(:), tmp_img(:).*ai(:));  % use local correlation to weight pixels
        if tmp_corr <min_similarity
            ind_ignore(ind_max) = true;
            cc(ind_max) = 0;
            fprintf('hmmm, try next\n');
            if save_fig
                saveas(gcf, fullfile(output_folder, sprintf('delete_%d.pdf', k_tried)));
            end
            
            % find the next candidate neuron
            temp = (reshape(tmp_ai,1, [])*Aem) ./ Aem_l2norm; % / norm(tmp_ai(:), 2);
            temp(ind_ignore)= -inf;
            % check the first 3
            for m=1:3
                [~, tmp_ind] = max(temp);
                ai = obj.reshape(Aem(:, tmp_ind), 3);
                tmp_img = image_local_corr(tmp_ai, ai, 7);
                %             tmp_corr = corr(tmp_img(:).*tmp_ai(:), tmp_img(:).*ai(:));
                tmp_corr = corr(tmp_img(:).*(ai(:)>0), ai(:));
                
                if tmp_corr >= min_similarity
                    cc(tmp_ind) = inf;  % find the next candidate neuron
                    break;
                else
                    temp(tmp_ind) = -inf;
                end
            end
            
            continue;
        end
    end
    
    %% add the neurons
    k_new  = k_new + 1;
    obj.C_raw(k_new,:) = ci_new_raw(:);
    obj.A(:, k_new) = ai_new(:);
    obj.C(k_new, :) = ci_new(:);
    obj.A_em(:, k_new) = ai(:);
    obj.ids(k_new) = segment_ids(ind_max);
    obj.match_status.status(k_new) = 1;
    obj.match_status.confidence(k_new) = 0;
    obj.match_status.scores(k_new) = 0;
    obj.match_status.em_ids{k_new} = segment_ids(ind_max);
    ind_used(ind_max) = true;
    fprintf('\t %d\n', k_new);
    
    %% subtract this components from Y
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
    if use_wrss
        temp = (Aem.^2)' * reshape(ai_new, [], 1);
    else
        temp = Aem' * reshape(ai_new, [], 1);
    end
    idx = (temp>0);
    
    % there are some neurons' temporal traces to be updated.
    if any(idx)  
        C_ = C_mmap.Data.C(:, idx);
        C_ = C_' - temp(idx)./norm_Aem(idx)'*ci_new(1:tsub:Tmax);
        switch lower(order_statistics)
            case 'l3norm'
                % l3norm of each C
                cc(idx) = sum(C_.^3, 2);
            case 'ss3mad'
                % sum of squares of values above 3 mad
                C_med = median(C_, 2);
                C_mad = median(abs(bsxfun(@minus, C_, C_med)), 2);
                cc(idx) = sum((C_.^2).*bsxfun(@gt, C_, C_med+C_mad*3), 2);
            case 'l3norm_xcorr'
                cc(idx) = sum(C_(:,1:(end-1)).*C_(:, 2:end), 2)./sum(C_(:,2:end).^2,2);
        end
        C_mmap.Data.C(:, idx) = C_';
        
        clear C_;
    end
end

%% finish 
fprintf('initialized %d neurons after trying %d EM components.\n', k_new-K_pre, k_tried); 
delete('temp_C.bin');
ind_del = (obj.match_status.status==0);
ind_del(1:K_pre) = false; 
if any(ind_del) 
    fprintf('early termination. delete the pre-assigned spaces.\n');
    obj.delete(ind_del);
end
if show_fig
    close(gcf);
end
