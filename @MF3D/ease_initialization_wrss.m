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
    'order_statistics', 'ss3mad', ... % method for ordering neurons
    'min_similarity', 0.6,...   % threshold for rejecting bad initialization
    'clear_results', false,...  % remove previous result
    'save_fig', false, ...        % save figures for visualizing the initialization step
    'show_fig', true, ...        % show figures for visualizing the initializaito step
    'K_candidate', 30000,...     % number of neurons to be considered
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

if preprocess_Y
    % select frames to be analyzed
    if isempty(obj.frame_range)
        Y = double(Y);
    else
        t0 = obj.frame_range(1);
        t1 = obj.frame_range(2);
        Y = double(Y(:, t0:t1));
    end
    
    % normalize data
    if obj.options.normalize_data
        sn = obj.reshape(obj.P.sn, 1);
        Y = bsxfun(@times, Y, 1./sn);
    end
    
    % remove all pixels outside of the EM volume
    if ~isempty(obj.spatial_range)
        Y(~obj.spatial_range, :) = 0; % remove pixels outside of the EM volume
    end
    
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
    tmp_f = median(Y, 1);
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

% normalize Aem by the l2 norm of columns
norm_Aem = sqrt(sum(Aem.^2, 1));

% remove segments that are empty
flag_use(norm_Aem==0) = false;

% ignore the above EM segments
Aem = Aem(:, flag_use);
segment_ids = segment_ids(flag_use, 1);
norm_Aem = norm_Aem(1, flag_use);

%% compute the summay statistics for each neuron
batch_size = 300;
K_total = size(Aem, 2);
nbatch = ceil(K_total/batch_size);
wrss_1 = zeros(1, K_total);
wrss_2 = zeros(1, K_total);

Q1 = (Aem.*Aem);
Q1_norm = sum(Aem.^3, 1);
Q2 = (Aem);
Q2_norm = sum(Aem, 1);

%% create a shift operation 
% [ii, jj, vv] = find(Q1); 
% d1 = obj.options.d1; 
% d2 = obj.options.d2; 
% d3 = obj.options.d3; 
% [idx_r, idx_c, idx_z] = ind2sub([d1,d2,d3], ii); 
% idx_r = idx_r - 3; 
% idx_c = idx_c - 3; 
% ii = sub2ind([d1,d2, d3], idx_r, idx_c, idx_z); 
% Q2 = sparse(ii, jj, vv, d1*d2*d3, K_total); 
% Q2_norm = Q1_norm; 

%%
if exist('temp.mat', 'file')
    load temp.mat;
else
    for m=1:nbatch
        sprintf('batch %d / %d\n', m, nbatch);
        % select Aem
        if m==nbatch
            idx = ((m-1)*batch_size+1) : K_total;
        else
            idx = (m-1)*batch_size + (1:batch_size);
        end
        
        % case: ai = pi
        temp = Q1(:, idx)'*Y;
        wrss_1(idx) = sum(temp.^2, 2)'./Q1_norm(idx);
        
        % case: ai = (pi>0)
        temp = Q2(:, idx)'*Y;
        wrss_2(idx) = sum(temp.^2, 2)'./Q2_norm(idx);
    end
    save temp.mat wrss_1 wrss_2; 
end
cc = wrss_1 - wrss_2;
cc(wrss_1==0) = -inf;

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
while (k_new < K_new+K_pre)
    % find the best candidate components
    cc(ind_used) = 0;
    cc(ind_ignore) = 0;
    [cc_max, ind_max] = max(cc);
    if cc_max==0
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
        if sum(ai_new(:).*ai(:))/norm(ai_new(ai>0), 2) < min_similarity
            ind_ignore(ind_max) = true;
            cc(ind_max) = 0;
            fprintf('hmmm, try next');
            if save_fig
                saveas(gcf, fullfile(output_folder, sprintf('delete_%d.pdf', k_tried)));
            end
            
            % find the next candidate neuron
            temp = reshape(ai_new(:)./norm(ai_new(:), 2), 1, []) * Aem./norm_Aem;
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
    
    %% update cc
    ci_new = reshape(ci_new, 1, []);
    ci_new = ci_new - mean(ci_new);
    
    ai_proj = Y*ci_new'/(ci_new*ci_new');
    aQ1 = reshape(ai_new, 1, []) * Q1;
    aQ1_proj = ai_proj' * Q1;
    idx = (aQ1>0);
    wrss_1(idx) = wrss_1(idx) + (ci_new*ci_new') * aQ1(idx).*...
        (aQ1(idx)-2*aQ1_proj(idx))./Q1_norm(idx);
    
    aQ2 = reshape(ai_new, 1, []) * Q2;
    aQ2_proj = ai_proj' * Q2;
    wrss_2(idx) = wrss_2(idx) + (ci_new*ci_new') * aQ2(idx).*...
        (aQ2(idx)-2*aQ2_proj(idx))./Q2_norm(idx);
    
    cc(idx) = wrss_1(idx) - wrss_2(idx);
    
    %% subtract this components
    Y = Y - reshape(ai_new, [], 1) * ci_new;
    if show_fig
        ai_res = reshape(ai_proj, size(ai_new)) - ai_new;
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
    
    %% update model variables
    if mod(k_new, 50)==0
        obj.options.nb = max(k_new/50, 3);
        b_old = obj.b;
        f_old = obj.f;
        
        % add back the old background
        for m=1:size(b_old, 2)
            ci_new = reshape(-f_old(m,:), 1, []);
            ci_new = ci_new - mean(ci_new);
            ai_new = b_old(:, m);
            
            ai_proj = Y*ci_new'/(ci_new*ci_new');
            aQ1 = reshape(ai_new, 1, []) * Q1;
            aQ1_proj = ai_proj' * Q1;
            wrss_1 = wrss_1 + (ci_new*ci_new') * aQ1.*...
                (aQ1-2*aQ1_proj)./Q1_norm;
            
            aQ2 = reshape(ai_new, 1, []) * Q2;
            aQ2_proj = ai_proj' * Q2;
            wrss_2 = wrss_2 + (ci_new*ci_new') * aQ2.*...
                (aQ2-2*aQ2_proj)./Q2_norm;
            Y = Y + b_old(:, m) *f_old(m,:);
        end
        
        % subtract the new background
        [obj.b, obj.f] = fit_svd_model(Y+b_old*f_old, obj.options.nb);
        % add back the old background
        for m=1:size(obj.b, 2)
            ci_new = reshape(obj.f(m,:), 1, []);
            ci_new = ci_new - mean(ci_new);
            ai_new = obj.b(:, m);
            
            ai_proj = Y*ci_new'/(ci_new*ci_new');
            aQ1 = reshape(ai_new, 1, []) * Q1;
            aQ1_proj = ai_proj' * Q1;
            wrss_1 = wrss_1 + (ci_new*ci_new') * aQ1.*...
                (aQ1-2*aQ1_proj)./Q1_norm;
            
            aQ2 = reshape(ai_new, 1, []) * Q2;
            aQ2_proj = ai_proj' * Q2;
            wrss_2 = wrss_2 + (ci_new*ci_new') * aQ2.*...
                (aQ2-2*aQ2_proj)./Q2_norm;
            Y = Y - obj.b(:, m) *obj.f(m,:);
        end
        
        
    end
end
