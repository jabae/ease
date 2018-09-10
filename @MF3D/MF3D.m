classdef MF3D < handle
    %% properties
    properties
        % spatial
        A;          % spatial components of neurons
        A_prev;     % previous estimation of A
        A_mask;     % spatial support to neuron A.
        A_corr;     % correlation betwen the raw video and C
        % temporal
        C;          % temporal components of neurons
        C_prev;     % previous estimation of C
        C_raw;      % raw traces of temporal components
        S;          % spike counts
        kernel;     % calcium dynamics. this one is less used these days.
        % background
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        W;          % a sparse weight matrix matrix
        b0;         % constant baselines for each pixel
        b0_new;  % correct the changes in b0 when we update A & C
        % optiosn
        options;    % options for model fitting
        P;          % some estimated parameters or parameters relating to data
        % data info
        Fs = nan;    % frame rate
        file = '';
        frame_range;  % frame range of the data
        spatial_range; % spatial area to be processed
        %quality control
        ids;        % unique identifier for each neuron
        match_status = struct('status', [], 'em_ids', [], 'confidence', ...
            []);  % the status of neuron matching. It's a struct variable with two fields:
        %  status: an array indicating the status for each neuron
        % -1: no match because the neuron is outsize of EM volume
        % 0: all EM segments are potential matches
        % 1: bingo! we found a match and we are pretty sure the match
        % is right
        % n: we got n potential matches, but we don't know which one is
        % the correct. The good thing is that we shrinked the range.
        
        %  em_ids: a cell aray for stroing em IDs relating to each status
        % status = -1, 0, it's []
        % status = 1, it's a scalar storing the matching ID
        % status = n, it's an array of all EM IDs
        
        labels = [];   % label the neuron as soma(1) or dendrite(2)
        tags;       % tags bad neurons with multiple criterion using a 16 bits number
        % ai indicates the bit value of the ith digits from
        % right to the left
        % a1 = 1 means that neuron has too few nonzero pixels
        % a2 = 1 means that neuron has silent calcium transients
        % a3 = 1 indicates that the residual after being deconvolved has
        % much smaller std than the raw data, which is uaually the case
        % when temporal traces are random noise
        % a4 = 1 indicates that the CNMF-E fails in deconvolving temporal
        % traces. this usually happens when neurons are false positives.
        
        % a4 = 1 indicates that the neuron has small PNR
        %others
        Cn;         % correlation image
        PNR;        % peak-to-noise ratio image.
        Coor;       % neuron contours
        neurons_per_patch;
        Df;         % background for each component to normalize the filtered raw data
        C_df;       % temporal components of neurons and background normalized by Df
        S_df;       % spike counts of neurons normalized by Df
        batches = cell(0);  % results for each small batch data
        file_id = [];    % file id for each batch.
        dataloader = [];    % data loader
        dataloader_denoised = [];  %data loader for the denoised data
        dataloader_raw = [];  % data loader for the raw data
        
        scores = [];  % matching scores
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = MF3D(varargin)
            obj.options = CNMFSetParms();
            obj.P =struct('mat_file', [], 'mat_data', [], 'indicator', '', 'k_options', 0, ...
                'k_snapshot', 0, 'k_del', 0, 'k_merge', 0, 'k_trim', 0, 'sn', [], ...
                'kernel_pars',[], 'k_ids', 0);
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:});
            end
            obj.kernel = create_kernel('exp2');
        end
        %% update parameters
        function updateParams(obj, varargin)
            obj.options = CNMFSetParms(obj.options, varargin{:});
        end
        %% reshape spatial data
        function Y = reshape(obj, Y, dim)
            % reshape the imaging data into diffrent dimensions
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            d3 = obj.options.d3;
            
            if dim==1
                Y=reshape(Y, d1*d2*d3, []);  %each frame is a vector
            else
                Y = reshape(full(Y), d1, d2, d3, []);    %each frame is an image
            end
        end
        
        %% delete neurons
        function delete(obj, ind)
            % write the deletion into the log file
            if ~exist('ind', 'var') || isempty(ind)
                return;
            end
            if islogical(ind)
                n_del = sum(ind(:));
            else
                n_del = length(ind);
            end
            if size(obj.scores,1) == size(obj.A, 2)
                obj.scores(ind, :) = [];
            end
            obj.A(:, ind) = [];
            obj.A_mask(:, ind) = [];
            if ~isempty(obj.A_corr)
                try  obj.A_corr(:, ind) = []; catch;  end
            end
            obj.C(ind, :) = [];
            obj.labels(ind) = [];
            if ~isempty(obj.S)
                try obj.S(ind, :) = []; catch; end
            end
            if ~isempty(obj.C_raw)
                try obj.C_raw(ind, :) = []; catch;  end
            end
            if isfield(obj.P, 'kernel_pars')&&(  ~isempty(obj.P.kernel_pars))
                try obj.P.kernel_pars(ind, :) = []; catch; end
            end
            try  obj.ids(ind) = [];   catch;   end
            try obj.tags(ind) =[]; catch; end
            if ~isempty(obj.match_status.status)
                obj.match_status.status(ind) = [];
                obj.match_status.em_ids(ind) = [];
                obj.match_status.confidence(ind) = [];
            end
            
            % save the log
            fprintf('%d neurons were deleted. \n', n_del);
        end
        
        
        %% view neurons
        showNeuron(obj, ind, orientation);
        
        %% determine spatial support
        Amask = determine_spatial_support(obj);
        
        %% update spatial components
        update_spatial(obj, Y);
        
        %% update temporal components
        [C_offset] = update_temporal(obj, Y, allow_deletion, wight_em);
        
        %% update background
        update_background(obj, Yr);
        
        %% reconstruct background
        function Yb = reconstruct_background(obj)
            if strcmpi(obj.options.background_model, 'svd')
                Yb =  bsxfun(@plus, obj.b*obj.f, obj.b0);
            end
        end
        
        %% play movie of demixing
        showDemixing(obj, Y, min_max, col_map, avi_nm, t_pause, ind_neuron, rot_info);
        showDenoised(obj, min_max, col_map, avi_nm, t_pause);
        
        %% play movie
        playMovie(obj, Y, min_max, col_map, avi_nm, t_pause);
        
        %% loading results
        load_super_pixels(obj, sp);
        
        %% determine whether a neuron is out of EM volumes
        function find_out_of_range(obj, em_scan_mask, thr)
            A_ = bsxfun(@times, obj.A, 1./sum(obj.A, 1));
            if ~exist('thr', 'var')
                thr = 0.01;
            end
            temp = reshape(em_scan_mask, 1, []) * A_;
            obj.match_status.status(temp<thr) = -1;
        end
        
        %% label no matches
        function label_no_matches(obj, thr)
            max_score = max(obj.scores, [],2);
            if ~exist('thr', 'var')
                thr = 0.1;
            end
            obj.match_status.status(max_score<thr) = -1;
        end
        %% computate the correlation coefficients between the video and neurons' activity
        function A_corr = calculate_corr(obj, Yr, type)
            if ~exist('type', 'var') || isempty(type)
                type = 'residual';
            end
            if ~exist('Yr', 'var') || isempty(Yr)
                if isempty(obj.frame_range)
                    Yr = evalin('base', 'Y_cnmf');
                else
                    temp = obj.frame_range;
                    Yr = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
                end
            end
            if strcmpi(type, 'residual')
                Yres = obj.reshape(Yr, 1) - obj.A*obj.C - obj.b*obj.f;
                Yres = bsxfun(@minus, Yres, mean(Yres, 2));
                var_Yres = sum(Yres.^2, 2);
                A_corr = zeros(size(obj.A));
                for m=1:size(obj.A, 2)
                    ci = obj.C(m, :);
                    ai = obj.A(:, m);
                    A_corr(:, m) = (Yres*ci' + ai*(ci*ci')) ...
                        ./sqrt(var_Yres+ai.^2*sum(ci.^2))/norm(ci, 2);
                end
            else
                Yr = obj.reshape(Yr, 1);
                Yr_mean = mean(Yr, 2);
                Yr_std = std(Yr, 0, 2);
                C_ = obj.C;
                C_mean = mean(C_, 2);
                C_std = std(C_, 0, 2);
                T = size(C_, 2);
                A_corr = (Yr * C_'/T-Yr_mean*C_mean')./(Yr_std*C_std');
            end
        end
        
        %% convert EM masks to a matrix
        function Aem = convert_matrix(obj, Aem)
            % get the dimension of Aem
            num_slices = length(Aem);
            d_em = size(Aem{1}, 1);
            
            % convert Aem to a matrix
            d1_ = obj.options.d1;
            d2_ = obj.options.d2;
            ssub = sqrt(d_em/(d1_*d2_));
            if ssub > 1
                for m=1:num_slices
                    temp = reshape(full(Aem{m}), d1_*ssub, d2_*ssub, []);
                    temp = imresize(temp, [d1_, d2_], 'box');
                    Aem{m} = reshape(temp, d1_ * d2_, []);
                end
            end
            Aem = sparse(cell2mat(Aem));
        end
        
        %% copy objects
        function obj_new = copy(obj)
            % Instantiate new object of the same class.
            obj_new = feval(class(obj));
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                obj_new.(p{i}) = obj.(p{i});
            end
        end
        
        
        %% calculate the matching score between all neuron shapes and EM components
        function scores = calculate_matching_scores(obj, Aem, method)
            if iscell(Aem)
                Aem = obj.convert_matrix(Aem);
            end
            em_mask = (sum(Aem, 2)<=0);   % constrain to the area within the EM volume
            
            % compute matching score
            A_ = obj.A;
            A_(em_mask, :) = 0;
            A_(A_<=0) = 0;
            Asum_em = sum(Aem, 1);      %   1*K_em
            AAsum_em = sum(Aem.^2, 1);  % 	1*K_em
            Asum_2p = sum(A_, 1);    %   1*K_2p
            AAsum_2p = sum(A_.^2, 1);%   1*K_2p
            AAsum_em2p = A_' * Aem;  % K_2p * K_em
            d = size(obj.A, 1);
            
            if strcmpi(method, 'corr')
                scores = (AAsum_em2p - Asum_2p' * Asum_em/d) ./...
                    (sqrt(AAsum_2p'-(Asum_2p.^2)'/d) * sqrt(AAsum_em-(Asum_em.^2)/d) );
            else % strcmpi(method, 'sim')
                scores = AAsum_em2p./(sqrt(AAsum_2p)' * sqrt(AAsum_em));
            end
            scores(:, Asum_em==0) = 0;
            scores = sparse(scores);
        end
        
        %% initialization given EM masks
        ind_voxels_em = initialize_em(obj, Aem, Y, options);
        
        %% deconvolve all temporal components
        C_ = deconvTemporal(obj, use_parallel, method_noise)
        
        %% run HALS
        function hals(obj, Y)
            if ~exist('Y', 'var') || isempty(Y)
                if isempty(obj.frame_range)
                    Y = evalin('base', 'Y_cnmf');
                else
                    temp = obj.frame_range;
                    Y = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
                end
            end
            %% get the spatial range
            ind = obj.spatial_range;
            if ~isempty(ind)
                Y(~ind, :) = 0; % remove pixels outside of the EM volume
            end
            
            %% run HALS
            obj.update_temporal(Y, false, true);
            obj.update_spatial(Y);
            obj.update_background(Y);
        end
        
        %% function remove false positives
        function ids = remove_false_positives(obj)
            ids = find(sum(obj.S, 2)==0);
            if isempty(ids)
                fprintf('all components are good \n');
            else
                obj.delete(ids);
            end
        end
        
        %% order ROIs
        function [srt] = orderROIs(obj, srt)
            % srt: sorting order
            nA = sqrt(sum(obj.A.^2));
            nr = length(nA);
            if nargin<2
                srt='srt';
            end
            K = size(obj.C, 1);
            
            if ischar(srt)
                if strcmp(srt, 'mean')
                    if obj.options.deconv_flag
                        temp = mean(obj.C,2)'.*sum(obj.A);
                    else
                        temp = mean(obj.C_raw.*(obj.C_raw>0),2)'.*sum(obj.A);
                    end
                    [~, srt] = sort(temp, 'descend');
                elseif strcmp(srt, 'sparsity_spatial')
                    temp = sqrt(sum(obj.A.^2, 1))./sum(abs(obj.A), 1);
                    [~, srt] = sort(temp);
                elseif strcmp(srt, 'sparsity_temporal')
                    temp = sqrt(sum(obj.C_raw.^2, 2))./sum(abs(obj.C_raw), 2);
                    [~, srt] = sort(temp, 'descend');
                elseif  strcmpi(srt, 'pnr')
                    pnrs = max(obj.C, [], 2)./std(obj.C_raw-obj.C, 0, 2);
                    [~, srt] = sort(pnrs, 'descend');
                elseif strcmpi(srt, 'l3l2')
                    l3l2 = sum(obj.C_raw.^3, 2) ./ (sum(obj.C_raw.^2,2).^(1.5));
                    [~, srt] = sort(l3l2, 'descend');
                elseif strcmpi(srt, 'temporal_cluster')
                    obj.orderROIs('pnr'); 
                    dd = pdist(obj.C_raw, 'cosine');
                    tree = linkage(dd, 'average');
                    srt = optimalleaforder(tree, dd);
                elseif strcmpi(srt, 'spatial_cluster')
                    obj.orderROIs('pnr');
                    A_ = bsxfun(@times, obj.A, 1./sqrt(sum(obj.A.^2, 1))); 
                    temp = 1-A_' * A_; 
                    dd = temp(tril(true(size(temp)), -1)); 
                    dd = reshape(dd, 1, []); 
                    tree = linkage(dd, 'average');
                    srt = optimalleaforder(tree, dd);                
                else %if strcmpi(srt, 'snr')
                    snrs = var(obj.C, 0, 2)./var(obj.C_raw-obj.C, 0, 2);
                    [~, srt] = sort(snrs, 'descend');
                end
            end
            obj.A = obj.A(:, srt);
            obj.C = obj.C(srt, :);
            
            try
                obj.A_mask = obj.A_mask(:, srt);
                obj.C_raw = obj.C_raw(srt,:);
                obj.S = obj.S(srt, :);
                obj.ids = obj.ids(srt);
                obj.match_status.status = obj.match_status.status(srt);
                obj.match_status.em_ids = obj.match_status.em_ids(srt);
                if ~isempty(obj.scores)
                    obj.scores = obj.scores(srt, :);
                end
                obj.A_corr = obj.A_corr(:, srt);
                
            end
        end
        
        %% normalize data
        function [Y, Y_sn] = normalize_data(obj, Y, update_var)
            Y = double(obj.reshape(Y,1));
            Y_sn = GetSn(bsxfun(@minus, Y, mean(Y, 2)));
            Y = bsxfun(@times, Y, 1./Y_sn);
            
            if ~exist('update_var', 'var') || isempty(update_var)
                update_var = false;
            end
            if update_var
                if ~isempty(obj.b)
                    obj.b = bsxfun(@times, obj.b, 1./Y_sn);
                end
                if ~isempty(obj.b0)
                    obj.b0 = obj.b0 ./Y_sn;
                end
                if ~isempty(obj.A)
                    obj.A = bsxfun(@times, obj.A, 1./Y_sn);
                end
            end
        end
        
        %% compute the residual
        function Yres = compute_residual(obj, Y)
            if ~exist('Y', 'var') || isempty(Y)
                if isempty(obj.frame_range)
                    Y = evalin('base', 'Y_cnmf');
                else
                    temp = obj.frame_range;
                    Y = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
                end
            end
            Y = obj.reshape(Y, 1);
            Yres = bsxfun(@minus, obj.reshape(Y,1) - obj.b*obj.f-obj.A*obj.C, obj.b0);
        end
        
        function rss = compute_rss(obj, Y)
            if ~exist('Y', 'var') || isempty(Y)
                if isempty(obj.frame_range)
                    Y = evalin('base', 'Y_cnmf');
                else
                    temp = obj.frame_range;
                    Y = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
                end
            end
            temp = bsxfun(@minus, obj.reshape(Y,1) - obj.b*obj.f-obj.A*obj.C, obj.b0);
            rss = sum(temp(:).^2);
        end
        %% decorrelate neuron traces
        C_ = decorrTemporal(obj, wd)
        
        %% show image
        function showImage(obj, ai, orientation, vlim)
            if iscell(ai)  % images one each plane are elements of the cell array
                [d1, d2, ~] = size(ai{1});
                d3 = length(ai);
            else   %
                if numel(ai) == 1
                    ai = obj.reshape(obj.A(:, ai), 3);
                elseif ndims(ai) ~=3
                    ai = obj.reshape(ai, 3);
                end
                [d1, d2, d3] = size(ai);
                img_max = max(ai(:))*0.8;
                if ~exist('vlim', 'var') || isempty(vlim)
                    vlim = [0, img_max];
                end
            end
            if ~exist('orientation', 'var') || isempty(orientation)
                orientation = 'vertical';
            end
            if strcmpi(orientation, 'horizental')
                h = d1+2;
                w = d3*(d2+2);
                pos_ax = [1-(d2+1)/w, 1/h, d2/w, 1-2/h];
                dpos = [-(d2+2)/w, 0, 0, 0];
            else
                h = d3*(d1+2);
                w = d2+2;
                pos_ax = [1/w, 1/h, 1-2/w, d1/h];
                dpos = [0, (d1+2)/h, 0, 0];
            end
            figure('papersize', [w, h] / max(w, h)*7);
            set(gcf, 'color', 'w', ...
                'position', [200, 200, 100*get(gcf, 'papersize')], ...
                'paperposition', [0, 0, get(gcf, 'papersize')]);
            
            for m=1:d3
                axes('position', pos_ax);
                if iscell(ai)
                    image(ai{m});
                else
                    imagesc(ai(:, :, d3-m+1), vlim);
                end
                axis equal tight; %equal off tight;
                box on; 
                set(gca, 'xtick', [], 'ytick', []); 
                pos_ax = pos_ax + dpos;
            end
        end
        
        %% show video
        function showVideo(obj, Y)
            Y = obj.reshape(Y, 3);
            [d1, d2, d3, T] = size(Y);
            ai_max = max(Y(:));
            figure('papersize', [d2, d1*d3]*5/d2);
            init_fig;
            kt = 3;
            for t=1:kt:T
                for m=1:d3
                    subplot(d3, 1, m);
                    imagesc(squeeze(Y(:, :, m, t)), [0, ai_max*0.9]);
                    axis equal off tight;
                end
                pause(0.1);
            end
        end
        
        %% get the rank of the matched neuron
        function ind_rank = get_match_rank(obj)
            K = size(obj.A, 2);
            ind_rank = nan(K, 1);
            tmp_status = obj.match_status;
            
            for cell_id = 1:K
                tmp_scores = obj.scores(cell_id, :);
                [~, em_sort_id] = sort(tmp_scores, 'descend');
                if tmp_status.status(cell_id) == 1
                    ind_rank(cell_id) = find(em_sort_id==tmp_status.em_ids{cell_id});
                end
            end
        end
        
        %% export results
        function results = export_results(obj, em_ids, output_path)
            ind = find(obj.match_status.status==1);
            results.A = obj.reshape(obj.A(:, ind), 3);
            results.A_mask = obj.reshape(obj.A_mask(:, ind), 3);
            results.C_raw = obj.C_raw(ind, :);
            results.C = obj.C(ind, :);
            results.S = obj.S(ind, :);
            tmp_ids = cell2mat(obj.match_status.em_ids(ind));
            results.EM_IDs = em_ids(tmp_ids, 1);
            results.confidence = obj.match_status.confidence;
            results.deconv_options = obj.options.deconv_options;
            if exist('output_path', 'var')
                save(output_path, 'results');
            end
        end
        % end of methods section
        
        %% fuse neurons' spatial footprints for visualizing them together
        function img = overlap_neurons(obj, ind, var_name)
            if ~exist('var_name', 'var') || isempty(var_name)
                var_name = 'A';
            end
            K = size(obj.A, 2);
            if ~exist('ind', 'var') || isempty(ind)
                ind = 1:K;  %#ok<NASGU>
            end
            
            % extract neurons to be visualized.
            A_ = obj.reshape(eval(sprintf('obj.%s(:, ind)', var_name)), 1);
            
            % fuse them
            d3 = obj.options.d3;
            img = cell(obj.options.d3, 1);
            a1 = obj.reshape(A_(:, 1), 3);
            a2 = obj.reshape(A_(:, 2), 3);
            for m=1:d3
                % plane by plane
                img1 = a1(:, :, m);
                img2 = a2(:, :, m);
                img{m} = imfuse(img1, img2); 
            end
        end
        
        %% compress results
        function compress(obj, k_score)
            obj.A = sparse(obj.A);
            obj.A_mask = sparse(obj.A_mask);
            obj.S = sparse(obj.S);
            [K_2p, K_em] = size(obj.scores);
            if ~exist('k_score', 'var') || isempty(k_score)
                % keep the top 100 matches
                k_score = min(100, K_em);
            end
            if size(obj.A, 2) == K_2p
                x = obj.scores;
                z = quantile(x, 1-k_score/K_em, 2);
                x(bsxfun(@lt, x, z)) = 0;
                obj.scores = sparse(full(x));
            end
        end
        
        %% reset the class object and clear results
        function obj = reset(obj)
            neuron = MF3D();
            neuron.options = obj.options;
            neuron.dataloader = obj.dataloader;
            neuron.dataloader_denoised = obj.dataloader_denoised;
            neuron.dataloader_raw = obj.dataloader_raw;
            neuron.Fs = obj.Fs;
            neuron.frame_range = obj.frame_range;
            neuron.spatial_range = obj.spatial_range;
            obj = neuron.copy();
        end
        
        %%        
        function [img, col0, AA] = overlapA(obj, ind, ratio, bg_white)
            %merge all neurons' spatial components into one singal image
            if nargin<2 || isempty(ind)
                AA = obj.A;
            else
                AA = obj.A(:, ind);
            end
            if nargin<3
                ratio = 0.3;
            end
            if ~exist('bg_white', 'var') || isempty(bg_white)
                bg_white = false; 
            end
            d3 = obj.options.d3; 
            v_max = 1;
            AA = bsxfun(@times, AA, v_max./max(AA,[],1));
            AA(bsxfun(@lt, AA, max(AA, [], 1)*ratio)) = 0;
            [d, K] = size(AA);
            
            if K==2
                col = [4, 2];
            elseif K==3
                col = [4, 2, 1];
            else
                col = randi(6, 1, K);
            end
            if bg_white
                col = 7-col; 
            end
            col0 = col;
            img = zeros(d, 3);
            for m=3:-1:1
                img(:, m) = sum(bsxfun(@times, AA, mod(col, 2)), 2);
                col = floor(col/2);
            end
            temp = img/max(img(:))*(2^16);
            temp = uint16(temp);
            if bg_white 
                temp = 2^16-temp; 
            end
            temp = obj.reshape(temp, 3); 
            img = cell(1, d3); 
            for m=1:d3
               img{m} = squeeze(temp(:, :, m, :));  
            end
        end
        
        %% compare two neurons
        function compare_pairs(obj, ind_pair, save_figs)
            if ~exist('save_figs', 'var') || isempty(save_figs)
                save_figs = false; 
            end
            img_A = obj.overlap_neurons(ind_pair, 'A');
            obj.showImage(img_A);
            tmp_pos0 = get(gcf, 'position'); 
            tmp_pos0(1:2) = [0, 0]; 
            set(gcf, 'name', '2p footprints', 'position', tmp_pos0); 
            if save_figs
                saveas(gcf, sprintf('pair_%d_%d_A.pdf', ind_pair(1), ind_pair(2)));
            end
            
            img_A_corr = obj.overlap_neurons(ind_pair, 'A_corr');
            obj.showImage(img_A_corr);
            set(gcf, 'name', 'correlation images', 'position', tmp_pos0+[tmp_pos0(3), 0, 0, 0]); 
            if save_figs
                saveas(gcf, sprintf('pair_%d_%d_A_corr.pdf', ind_pair(1), ind_pair(2)));
            end
            img_A_mask = obj.overlap_neurons(ind_pair, 'A_mask');
            obj.showImage(img_A_mask);
            set(gcf, 'name', 'EM footprints', 'position', tmp_pos0+[2*tmp_pos0(3),0, 0 0]); 
            if save_figs
                saveas(gcf, sprintf('pair_%d_%d_A_mask.pdf', ind_pair(1), ind_pair(2)));
            end
            
            %% temporal
            figure('papersize', [14, 4], 'name', 'temporal traces');
            init_fig;
            T = size(obj.C, 2);
            pos = [0.07, 0.64, 0.85, 0.35];
            for m=1:2
                axes('position', pos);
                plot((1:T)/15, obj.C_raw(ind_pair(m), :));
                hold on;
                plot((1:T)/15, obj.C(ind_pair(m), :));
                axis  tight;
                if m~=length(ind_pair)
                    set(gca, 'xticklabel', []);
                end
                pos = pos + [0, -pos(4)-0.01, 0, 0];
                ylabel(sprintf('cell %d', ind_pair(m)));
            end
            xlabel('Time (sec.)');
            tmp_pos = get(gcf, 'position'); 
            tmp_pos(1:2) =[0, 0]; 
            set(gcf, 'position', tmp_pos); 
            if save_figs
                saveas(gcf, sprintf('pair_%d_%d_temporal.pdf', ind_pair(1), ind_pair(2)));
            end
        end
        
        %% normalize data 
        function normalize(obj, method)
            if ~exist('method', 'var') || isempty(method)
                method = 'c_noise'; 
            end
                
            if strcmpi(method, 'a_max')
                norm_A = max(obj.A, [], 1);
            elseif strcmpi(method, 'c_noise')
                norm_A = 1./GetSn(obj.C_raw)'; 
            elseif strcmpi(method, 'c_max')
                norm_A = 1./max(obj.C,[], 2)'; 
            else
                norm_A = 1./max(obj.C,[], 2)';
            end
            obj.A = bsxfun(@times, obj.A, 1./norm_A);
            obj.C = bsxfun(@times, obj.C, norm_A');
            obj.C_raw = bsxfun(@times, obj.C_raw, norm_A'); 
            obj.S = bsxfun(@times, obj.S, norm_A'); 
            
        end
        
        %% get EM ids 
        function em_ids = get_em_ids(obj, ind, em_info)
            ind = cell2mat(obj.match_status.em_ids(ind)); 
            em_ids = em_info(ind, 1); 
        end
    end
end














