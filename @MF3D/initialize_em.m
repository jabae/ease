function ind_voxels_em = initialize_em(obj, Aem, Y, Kmax, options)
%% initialize A & C in the CNMF model using spatial masks given by EM segments.
%{
    This function is used for initializing neurons given a pool of spatial
    masks. These masks for segmented from EM data. 
%}

%% inputs: 
%{
    Aem: d_em * K_em matrix, spatial footprints given by EM segments
    Y: d * T, data
    Kmax: integer, the maximum number of neurons to be initialized.
    options: struct variable for containing all options 
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
if iscell(Aem)
    % convert cell array to a matrix
    Aem =  obj.convert_matrix(Aem);
end

if ~exist('options', 'var') || isempty(options)
    init_method = 'tf';
    order_statistics = 'thresh';
    min_similarity = 0.6;
    clear_results = false;
    save_fig = true;
    show_fig = true; 
    K_candidate = 3000;
else
    init_method = options.init_method;
    order_statistics = options.order_statistics;
    min_similarity = options.min_similarity;
    clear_results = options.clear_results;
    save_fig = options.save_fig;
    show_fig = options.show_fig; 
    K_candidate = options.K_candidate;
end

if ~exist('Kmax', 'var')    % number of neurons to be initialized 
    Kmax = 50;
end

ind_voxels_em = sparse(sum(Aem, 2)>0);     % find voxels within EM volumes
Aem_sum = sum(Aem, 1);          % l1 norm of each ai_em
if ~clear_results
    %ignore existing matches
    ids_matched = cell2mat(obj.match_status.em_ids(obj.match_status.status==1));
    Aem_sum(ids_matched) = 0;
end
K_em = length(Aem_sum>0);    % number of EM components.

[~, temp] = sort(Aem_sum, 'descend');    % find components within the scanning planes
tmp_ind = temp(1:min(K_candidate, K_em));
em_ids = tmp_ind;           % ids of the candidate neurons

% normalize the spatial masks of the selected neurons 
A_selected = full(Aem(:, tmp_ind));       % spatial masks of these neurons
A_norm = sqrt(sum(A_selected.^2,1));     % normalize spatial components
A_ = bsxfun(@times, A_selected, 1./A_norm);

if save_fig
    tmp_folder = evalin('base', 'ease.output_folder');
    output_folder = fullfile(tmp_folder, ['initialization_', get_date()]);
    mkdir(output_folder);
end

%% load data
if ~exist('Y', 'var') || isempty(Y)
    dl = obj.dataloader_denoised;
    Y = dl.load_tzrc();
end

if clear_results
    Y = obj.reshape(Y,1);
    obj.update_background(Y);
    Y = obj.reshape(Y,1) - obj.b*obj.f;
else
    Y = obj.reshape(Y, 1) - obj.A*obj.C - obj.b*obj.f;
end

%% center data and remove the background
Y = bsxfun(@minus, Y, mean(Y, 2));

%% create matrices for storing results
[d, T] = size(Y);

if clear_results
    K_pre = 0;
    obj.A = zeros(d, Kmax);
    obj.A_mask = zeros(d, Kmax);
    obj.C = zeros(Kmax, T);
    obj.C_raw = zeros(Kmax, T);
    obj.S = zeros(Kmax, T);
    obj.ids((end+1):(end+Kmax)) = 0;
    obj.match_status.status = zeros(1, Kmax);
    obj.match_status.em_ids = cell(1, Kmax);
else
    K_pre = size(obj.A, 2);     % number of existing neurons
    obj.A = [obj.A, zeros(d, Kmax)];
    obj.A_mask = [obj.A_mask, zeros(d, Kmax)];
    obj.C = [obj.C; zeros(Kmax, T)];
    obj.C_raw = [obj.C_raw; zeros(Kmax, T)];
    obj.S = [obj.S; zeros(Kmax, T)];
    obj.ids((end+1):(end+Kmax)) = 0;
    obj.match_status.status((end+1):(end+Kmax)) = 0;
    obj.match_status.em_ids = [obj.match_status.em_ids, cell(1, Kmax)];
end
k_new = K_pre;
k_tried = 0;
ind_ignore = false(1, size(A_, 2));
ind_used = ind_ignore;

%% estimate temporal traces for each EM
C_ = A_' * Y;
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

while (k_new < Kmax+K_pre) && (k_tried<size(A_,2))
    %% find the best components
    cc(ind_ignore) = 0; 
    cc(ind_used) = 0;
    [~, ind_max] = max(cc);
    k_tried = k_tried +1;
    ind_ignore(ind_max) = true; 
    
    %% extract the corresponding EM segment ai and ci
    ai = obj.reshape(A_(:, ind_max), 3);
    ai_mask = imdilate(ai>0, strel('square', 4)); %imdilate(repmat(sum(ai, 3)>0, [1, 1, 3]), strel('square', 3));
    ci = C_(ind_max, :); 
    
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
        [ai_new(ind_in), ci_new] = initialize_ac_tf(Yin, Yout, ai(ind_in(:)), ci0, 10);
        [~, sn] = estimate_baseline_noise(ci_new);
        ci_new_raw = ci_new; 
        ci_new = deconvolveCa(ci_new, deconv_options, 'sn', sn);
        ci_new = reshape(ci_new, 1,[]);
        
        tmp_ci = ci_new - mean(ci_new);
        ai_proj = (Y * tmp_ci')/(tmp_ci*tmp_ci'); %.*(ai(:)>0);
        ai_proj = obj.reshape(ai_proj, 3);
    elseif strcmpi(init_method, 'iteration')
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
%             ai_new = ai_new / norm(ai_new, 2);
            
            %% refine ci
            ci_new = ((ai_new.*ai)' * Y) / ((ai_new.*ai)'*ai_new);
            ci_new_raw = ci_new; 
            ci_new = deconvolveCa(ci_new, deconv_options, 'sn', mad(ci_new));
        end
        if std(ci_new) ==0
            ind_ignore(ind_max) = true;
            continue;
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
    if sum(ai_new(:).*ai(:))/norm(ai_new(:), 2) < min_similarity
        ind_ignore(ind_max) = true;
        cc(ind_max) = 0;
        display('hmmm, try next');
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
    
    %% add the neurons
    k_new  = k_new + 1;
    obj.C_raw(k_new,:) = ci_new_raw(:);
    obj.A(:, k_new) = ai_new(:);
    obj.C(k_new, :) = ci_new(:);
    obj.A_mask(:, k_new) = ai(:);
    obj.ids(k_new) = em_ids(ind_max);
    obj.match_status.status(k_new) = 1;
    obj.match_status.em_ids{k_new} = em_ids(ind_max);
    ind_used(ind_max) = true;
    fprintf('\t %d\n', k_new);
    
    %% subtract this components and compute correlation image
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
    if save_fig
        saveas(gcf, fullfile(output_folder, sprintf('%d.pdf', k_new)));
    end
    
    %% update C and cc
    C_ = C_ - A_' * reshape(ai_new, [], 1) * reshape(ci_new, 1, []);
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
    
    pause(0.1);
    if mod(k_new, 20)==0
        % update background
        tmpY = Y + obj.b*obj.f;
        [u, s, v] = svdsecon(tmpY, 1);
        obj.b = u;
        obj.f = s*v';
        Y = tmpY - u*s*v';
        ind_ignore = false(size(ind_ignore));
    end
    
end
