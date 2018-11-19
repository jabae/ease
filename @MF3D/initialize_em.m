function ind_voxels_em = initialize_em(obj, Aem, Y, options, black_list, white_list)
%% initialize A & C in the CNMF model using spatial masks given by EM segments.
%{
    This function is used for initializing neurons given a pool of spatial
    masks. These masks for segmented from EM data.
%}

%% inputs:
%{
    Aem: d_em * K_em matrix, spatial footprints given by EM segments
    Y: d * T, data
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

%% choose candidate EM segments

% make sure that Aem is a matrix
if iscell(Aem)
    Aem =  obj.convert_matrix(Aem);
end

% options for the initialization
if ~exist('options', 'var') || isempty(options)
    init_method = 'tf';     %method for initializing single neurons
    order_statistics = 'thresh'; % method for ordering neurons
    min_similarity = 0.6;   % threshold for rejecting bad initialization
    clear_results = false;  % remove previous result
    save_fig = true;        % save figures for visualizing the initialization step
    show_fig = true;        % show figures for visualizing the initializaito step
    K_candidate = 3000;     % number of neurons to be considered
    K_new = 50;             % number of neurons to be added
    min_pnr = 5;                % the minimum peak-to-noise ratio for a good trace
else
    init_method = options.init_method;
    order_statistics = options.order_statistics;
    min_similarity = options.min_similarity;
    clear_results = options.clear_results;
    save_fig = options.save_fig;
    show_fig = options.show_fig;
    K_candidate = options.K_candidate;
    K_new = options.K_new;
    min_pnr = options.min_pnr;
end

% order neurons and determine the list of candidate neurons
ind_voxels_em = sparse(sum(Aem, 2)>0); % find voxels within EM volumes
Aem_sum = sum(Aem, 1);          % l1 norm of each ai_em
if ~clear_results
    %ignore existing matches
    ids_matched = cell2mat(obj.match_status.em_ids(obj.match_status.status==1));
    Aem_sum(ids_matched) = 0;
end
if exist('black_list', 'var') && ~isempty(black_list)
    Aem_sum(black_list) = 0;
end
if exist('white_list', 'var') && ~isempty(white_list)
    %     [v, idx] = sort(Aem_sum(white_list), 'descend');
    %     em_ids = white_list(idx);
    em_ids = white_list(Aem_sum(white_list)>0);
    K_new = length(em_ids);
    with_quality_control = false;
else
    K_em = length(Aem_sum>0);    % number of EM components.
    [~, temp] = sort(Aem_sum, 'descend');
    tmp_ind = temp(1:min(K_candidate, K_em));
    em_ids = tmp_ind;           % ids of the candidate neurons
    with_quality_control = true;
    white_list = [];
end

% normalize the spatial masks of the selected neurons
A_selected = full(Aem(:, em_ids));       % spatial masks of these neurons
A_norm = sqrt(sum(A_selected.^2,1));     % normalize spatial components
A_ = bsxfun(@times, A_selected, 1./A_norm);

if save_fig
    tmp_folder = evalin('base', 'ease.output_folder');
    output_folder = fullfile(tmp_folder, ['initialization_', get_date()]);
    mkdir(output_folder);
end

%% load data
if ~exist('Y', 'var') || isempty(Y)
    if isempty(obj.frame_range)
        Y = evalin('base', 'Y_cnmf');
    else
        temp = obj.frame_range;
        Y = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
    end
    %     dl = obj.dataloader_denoised;
    %     Y = dl.load_tzrc();
end

if clear_results
    Y = obj.reshape(Y,1);
    obj.A = [];
    obj.C = [];
    obj.labels = [];
    
    tmp_f = median(Y, 1);
    tmp_f = tmp_f - mean(tmp_f);
    tmp_b = (Y*tmp_f')/(tmp_f*tmp_f');
    Y = obj.reshape(Y,1) - tmp_b*tmp_f;
else
    Y = obj.reshape(Y, 1) - obj.A*obj.C - obj.b*obj.f;
end

%% get the spatial range
spatial_range = obj.spatial_range;
if ~isempty(spatial_range)
    Y(~spatial_range, :) = 0; % remove pixels outside of the EM volume
end
%% center data and remove the background
Y = bsxfun(@minus, Y, mean(Y, 2));

%% create matrices for storing results
[d, T] = size(Y);

if clear_results
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
    obj.match_status.confidence = ones(1, K_new)*5;
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
    obj.match_status.confidence((end+1):(end+K_new)) = 5;
    obj.match_status.em_ids = [obj.match_status.em_ids, cell(1, K_new)];
end
k_new = K_pre;
k_tried = 0;
ind_ignore = false(1, size(A_, 2));
ind_used = ind_ignore;

%% estimate temporal traces for each EM
tsub = 4;
C_ = A_' * Y(:, 1:tsub:end);
switch lower(order_statistics)
    case 'skewness'
        % skewness of each C
        cc = skewness(C_, 0, 2);
    case 'ss3mad'
        % sum of squares of values above 3 mad
        C_med = median(C_, 2);
        C_mad = median(abs(bsxfun(@minus, C_, C_med)), 2);
        cc = sum((C_.^2).*bsxfun(@gt, C_, C_med+C_mad*3), 2);
    otherwise
        cc = skewness(C_, 0, 2);
end

%% canvas for ploting results
if show_fig
    figure('papersize', [10.08, 6.83]);
    init_fig;
    colormap jet;
    Nrows = 5;
end

%% start initialization
deconv_options = obj.options.deconv_options;

while (k_new < K_new+K_pre) && (k_tried<size(A_,2))
    %% find the best components
    if isempty(white_list)
        cc(ind_ignore) = 0;
        cc(ind_used) = 0;
        [~, ind_max] = max(cc);
        k_tried = k_tried +1;
        ind_ignore(ind_max) = true;
    else
        k_tried = k_tried + 1;
        ind_max = k_tried;
    end
    %% extract the corresponding EM segment ai and ci
    ai = obj.reshape(A_(:, ind_max), 3);
    ai_mask = imdilate(ai>0, strel('square', 4)); %imdilate(repmat(sum(ai, 3)>0, [1, 1, 3]), strel('square', 3));
    ci = reshape(ai, 1, []) * Y; %C_(ind_max, :);
    
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
    if strcmpi(init_method, 'tf')
        % choose pixels within the EM masks and pixels surrounding the EM
        % masks
        ind_in = obj.reshape(ai_mask, 3);
        ind_out = xor(ind_in, imdilate(ai_mask, strel('square', 8)));
        Yin = Y(ind_in(:), :);
        Yout = Y(ind_out(:), :);
        ci0 = ci;
        ai_new = double(ind_in);
        [ai_new(ind_in), ci_new] = initialize_ac_tf(Yin, Yout, ai(ind_in(:)), ci0, 15);
        [~, sn] = estimate_baseline_noise(ci_new);
        
        ci_new_raw = ci_new;
        ci_new = deconvolveCa(ci_new, deconv_options, 'sn', sn);
        ci_new = reshape(ci_new, 1,[]);
        pnr = max(ci_new) / std(ci_new_raw - ci_new);
        
        tmp_ci = ci_new - mean(ci_new);
        ai_proj = (Y * tmp_ci')/(tmp_ci*tmp_ci'); %.*(ai(:)>0);
        ai_proj = obj.reshape(ai_proj, 3);
    else%if strcmpi(init_method, 'semi-nmf')
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
    end
    if std(ci_new(2:end)) == 0
        ind_ignore(ind_max) = true;
        if isempty(white_list)
            continue;
        else
            ai_new = ai;
            ci_new = zeros(size(ci_new));
        end
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
            temp = reshape(ai_new(:)./norm(ai_new(:), 2), 1, []) * A_;
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
    obj.ids(k_new) = em_ids(ind_max);
    obj.match_status.status(k_new) = 1;
    obj.match_status.confidence(k_new) = 5;
    obj.match_status.em_ids{k_new} = em_ids(ind_max);
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
    
    %% update C and cc
    C_ = C_ - A_' * reshape(ai_new, [], 1) * reshape(ci_new(1:tsub:end), 1, []);
    if isempty(white_list)
        switch lower(order_statistics)
            case 'skewness'
                % skewness of each C
                cc = skewness(C_, 0, 2);
            case 'ss3mad'
                % sum of squares of values above 3 mad
                C_med = median(C_, 2);
                C_mad = median(abs(bsxfun(@minus, C_, C_med)), 2);
                cc = sum((C_.^2).*bsxfun(@gt, C_, C_med+C_mad*3), 2);
            otherwise
                cc = skewness(C_, 0, 2);
        end
    end
end
