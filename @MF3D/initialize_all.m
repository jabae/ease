function ind_voxels_em = initialize_all(obj, Aem, Y, Kmax, options)
%% initialize A & C in the CNMF model using spatial masks given by EM segments.
% inputs:
%   Aem: d_em * K_em matrix, spatial footprints given by EM segments
%   Y: d * T, data
%   Kmax: integer, the maximum number of neurons to be initialized.
%   save_fig: boolean, save results or not.
%   clear_results: boolean, clear the earlier results

% outputs:
%   ind_voxels_em: spatial mask for identifying EM volume

% author:
%   Pengcheng Zhou, Columbia University, 2018

%% choose candidate EM segments
% convert cell array to a matrix
if iscell(Aem)
    Aem =  obj.convert_matrix(Aem);
end

if ~exist('options', 'var') || isempty(options)
    init_method = 'regression';
    clear_results = false;
    save_fig = true;
else
    init_method = options.init_method;
    clear_results = options.clear_results;
    save_fig = options.save_fig;
end

ind_voxels_em = sparse(sum(Aem, 2)>0);     % find voxels within EM volumes
temp = sum(Aem>0, 1);
K_em = length(temp);    % number of EM components.
if ~clear_results
    %ignore existing matches
    ids_matched = cell2mat(obj.match_status.em_ids(obj.match_status.status==1));
    temp(ids_matched) = 0;
end
if ~exist('Kmax', 'var')
    Kmax = 300;
end

[~, temp] = sort(temp, 'descend');    % find components within the scanning planes
tmp_ind = temp(1:min(Kmax, K_em));
em_ids = tmp_ind;           % ids of the selected neurons

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
Yres = bsxfun(@minus, Y, mean(Y, 2));
Yres(~ind_voxels_em, :) = 0;
clear Y;

%% create matrices for storing results
[d, T] = size(Yres);

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
k_tried = 0;
ind_ignore = false(1, size(A_, 2));
ind_used = ind_ignore;
C_ = A_' * Yres;

%%
figure('papersize', [10.08, 6.83]);
init_fig;
colormap jet;
Nrows = 5;
deconv_options = obj.options.deconv_options;
for k_new = (K_pre+1):(K_pre+Kmax) 
    ind_max = k_new + K_pre; 
    k_tried = k_tried +1;
    
    % extract the corresponding EM segment ai and ci
    ai = obj.reshape(A_(:, ind_max), 3);
    ai_mask = imdilate(ai>0, strel('square', 3));
    ci = C_(ind_max, :); %*sn(ind_max);
    ai_max = max(ai(:));
    
    % show ai and ci
    for m=1:3
        subplot(Nrows, 3, m);
        imagesc(ai(:, :, m), [0,ai_max/2]);
        axis equal off tight;
    end
    subplot(Nrows, 3, 4:6); cla;
    plot(ci, 'linewidth', 1);
    axis tight;
    set(gca, 'fontsize', 16);
    
    %% initialize (ai, ci)
    if strcmpi(init_method, 'tf')
        ind_in = obj.reshape(ai_mask, 3);
        ind_out = xor(ind_in, imdilate(ai_mask, strel('square', 8)));
        ind_out(~ind_voxels_em) = false;
        Yin = Yres(ind_in(:), :);
        Yout = Yres(ind_out(:), :);
        ci0 = ci;
        ai_new = double(ind_in);
        [ai_new(ind_in), ci_new] = initialize_ac_tf(Yin, Yout, ai(ind_in(:)), ci0, 10);
        [~, sn] = estimate_baseline_noise(ci_new);
        ci_new_raw = ci_new; 
        ci_new = deconvolveCa(ci_new, deconv_options, 'sn', sn);
        ci_new = reshape(ci_new, 1,[]);
        tmp_ci = ci_new - mean(ci_new);
        ai_proj = (Yres * tmp_ci')/(tmp_ci*tmp_ci'); %.*(ai(:)>0);
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
            ai_proj = (Yres * ci_new')/(ci_new * ci_new'); %.*(ai(:)>0);
            ai_new = max(0, ai_proj.*ind_mask);
%             ai_new = ai_new / norm(ai_new, 2);
            
            %% refine ci
            ci_new = ((ai_new.*ai)' * Yres) / ((ai_new.*ai)'*ai_new);
            ci_new_raw = ci_new; 
            ci_new = deconvolveCa(ci_new, deconv_options, 'sn', mad(ci_new));
        end
        if std(ci_new) ==0
            ind_ignore(ind_max) = true;
            continue;
        end
    end

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
    subplot(Nrows, 3, 4:6); hold on;
    plot(ci_new, 'r', 'linewidth', 1);
    
    %% add the neurons
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
    ai_res = obj.reshape(ai_proj, 3) - obj.reshape(ai_new, 3); 
    for m=1:3
        subplot(Nrows, 3, 12+m); cla;
        imagesc(ai_res(:, :, m), [0, ai_new_max/3]);
        axis equal off tight;
    end
    if save_fig
        saveas(gcf, fullfile(output_folder, sprintf('%d.pdf', k_new)));
    end
  
end
