classdef MF3D < handle
    %% properties
    properties
        % spatial
        A;          % spatial components of neurons
        A_em;       % spatial footprints of the corresponding EM components 
%         A_em;     % spatial mask for constraining the neuron footprints
        A_corr;     % correlation betwen the raw video and C
        % temporal
        C;          % temporal components of neurons
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
        white_list = []; % EM components that should be considered.
        black_list = []; % EM components that should be ignored
        %quality control
        ids;        % unique identifier for each neuron
        match_status = struct('status', [], 'em_ids', [], 'confidence', ...
            [], 'scores', []);  % the status of neuron matching. It's a struct variable with two fields:
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
        whitelist = [];  % add the neuron to white list (yes 1; no 0)
        labels = [];   % label the neuron as soma(1) or dendrite(2); negative value means that the neuron has been verified. 
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
        tuning_curve = {}; % tuning curve of each neuron
        
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
        function delete(obj, ind, add_to_black_list)
            % write the deletion into the log file
            if ~exist('ind', 'var') || isempty(ind)
                return;
            end
            if ~exist('add_to_black_list', 'var') || isempty(add_to_black_list)
                add_to_black_list = false;
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
            obj.A_em(:, ind) = [];
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
                % add these deleted neurons to a blacklist
                %                 temp = cell2mat(obj.match_status.em_ids(ind));
                %                 if add_to_black_list
                %                     obj.black_list = unique([obj.black_list; reshape(temp, [], 1)]);
                %                 end
                % delete
                obj.match_status.status(ind) = [];
                obj.match_status.em_ids(ind) = [];
                obj.match_status.confidence(ind) = [];
                obj.match_status.scores(ind) = [];
            end
            
            % save the log
            fprintf('%d neurons were deleted. \n', n_del);
        end
        
        
        %% view neurons
        showNeuron(obj, ind, orientation);
        
        %% determine spatial support
        Amask = determine_spatial_support(obj);
        
        %% update spatial components
        update_spatial(obj, Y, with_em, preprocess_Y);
        
        %% update temporal components
        [C_offset] = update_temporal(obj, Y, allow_deletion, wight_em, preprocess_Y);
        
        %% update background
        update_background(obj, Yr, preprocess_Y);
        
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
        function scores = calculate_matching_scores(obj, Aem)
            % check the spatial mask
            if iscell(Aem)
                Aem = obj.convert_matrix(Aem);
            end
            em_mask = (sum(Aem, 2)<=0);   % constrain to the area within the EM volume
            Aem(em_mask, :) = [];
            n = size(Aem, 1);
            Aem_sum = sum(Aem, 1);
            Aem_std = sqrt(sum(Aem.^2,1) - (Aem_sum.^2/n));
            
            %% compute matching score
            A_ = obj.A;
            A_(em_mask, :) = [];
            A_(A_<=0) = 0;
            P_ = obj.A; %obj.A_em .* tmpA_corr;
            P_(em_mask, :) = [];
            P_sum = sum(P_, 1);
            P_std = sqrt(sum(P_.^2, 1) - (P_sum.^2/n));
            
            temp1 = bsxfun(@times, A_'*(Aem>0), 1./sum(A_,1)'); % explained signal with different masks
            %             temp2 = (P_'*Aem-mean(P_,1)'*mean(Aem,1)) ./ ...
            %                 (std(P_, 0, 1)'*std(Aem, 0, 1));
            temp2 = P_'*Aem - P_sum'*(Aem_sum/n);
            temp2 = bsxfun(@times, temp2, 1./P_std');
            temp2 = bsxfun(@times, temp2, 1./Aem_std);
            
            temp = temp1 .* temp2;
            temp(isnan(temp)) = 0;
            scores = sparse(temp);
            obj.scores = scores;
            %
            %             K = size(obj.A, 2);
            %             Kem = size(scores, 2);
            %             ind = sub2ind([K, Kem], 1:K, cell2mat(obj.match_status.s
            %             obj.match_status.scores = scores(
        end
        
        %% evaluate matching performance
        confidence = evaluate_matching_confidence(obj, Yr, Aem, segment_ids);
        
        %% initialization given EM masks
        ind_voxels_em = initialize_em(obj, Aem, Y, options, black_list, white_list);
        %% initialization given EM masks
        [ai, ci, si, ci_raw] = initialize_one(obj, ai_em, Y);
        
        %% initialize multiple neurons given their EM masks
        initialize_multiple(obj, Aem, segment_ids, Y, options);
        
        %% deconvolve all temporal components
        C_ = deconvTemporal(obj, use_parallel, method_noise)
        
        %% preprocess Y
        function Y = preprocess(obj, Y)
            Y = obj.reshape(Y, 1);
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
        %% run HALS
        function hals(obj, Y, weight_em, with_EM_info, preprocess_Y)
            
            if ~exist('weight_em', 'var') || isempty(weight_em)
                weight_em = true;
            end
            if ~exist('with_EM_info', 'var') || isempty(with_EM_info)
                with_EM_info = true;
            else
                with_EM_info = false;
                weight_em = false;
            end
            if ~exist('preprocess_Y', 'var') || isempty(preprocess_Y)
                preprocess_Y = obj.options.pre_process_data;
            end
            
            %% preprocess Y
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
                preprocess_Y = false;
            end
            
            %% get the spatial range
            ind = obj.spatial_range;
            if ~isempty(ind)
                Y(~ind, :) = 0; % remove pixels outside of the EM volume
            end
            
            %% run HALS
            obj.update_background(Y, preprocess_Y);
            obj.update_temporal(Y, false, weight_em, preprocess_Y);
            obj.update_spatial(Y, with_EM_info, preprocess_Y);
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
        [srt] = orderROIs(obj, srt);
        
        %% normalize data
        function [Y, Y_sn] = normalize_data(obj, Y, update_var)
            Y = double(obj.reshape(Y,1));
            if isempty(obj.P.sn)
                T = min(size(Y, 2), 5000);
                Y_sn = GetSn(bsxfun(@minus, Y(:, 1:T), mean(Y(:, 1:T), 2)));
                obj.P.sn = Y_sn;
            else
                Y_sn = reshape(obj.P.sn, [], 1);
            end
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
            if obj.options.pre_process_data
                Y = obj.preprocess(Y);
            end
            
            Yres = bsxfun(@minus, obj.reshape(Y,1) - ...
                obj.b*obj.f - obj.A*obj.C, obj.b0);
        end
        
        %% compute RSS
        function rss = compute_rss(obj, Y, preprocess_Y)
            if ~exist('preprocess_Y', 'var') || isempty(preprocess_Y)
                preprocess_Y = true;
            end
            % preprocess Y
            if preprocess_Y
                Y = obj.preprocess(Y);
            end
            
            temp = bsxfun(@minus, obj.reshape(Y,1) - obj.b*obj.f-obj.A*obj.C, obj.b0);
            rss = sum(temp(:).^2);
        end
        %% decorrelate neuron traces
        C_ = decorrTemporal(obj, wd)
        
        %% show image
        function showImage(obj, ai, orientation, vlim, pixel_size, color_scalebar)
            if ~exist('pixel_size', 'var')
                pixel_size = [];
            end
            if ~exist('color_scalebar', 'var') || isempty(color_scalebar)
                color_scalebar = 'w';
            end
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
                if diff(vlim)==0
                    vlim(2) = vlim(2) + 1;
                end
            end
            
            if ~exist('orientation', 'var') || isempty(orientation)
                orientation = 'vertical';
            end
            if strcmpi(orientation, 'horizontal')
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
                axis tight; %equal off tight;
                box on;
                set(gca, 'xtick', [], 'ytick', []);
                pos_ax = pos_ax + dpos;
                
                if ~isempty(pixel_size) && (m==1) % draw scale bar
                    hold on;
                    plot(d2-3- [0, 20]/pixel_size, d1-[1,1]*3, 'color', color_scalebar, 'linewidth', 3);
                    %  text(d2-5-20/pixel_size, d1-7, '20um', 'color', 'w', ...
                    %     'fontsize', round(20/pixel_size/d2*75), 'fontweight', 'bold');
                end
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
            results.A_em = obj.reshape(obj.A_em(:, ind), 3);
            results.C_raw = obj.C_raw(ind, :);
            results.C = obj.C(ind, :);
            results.S = obj.S(ind, :);
            results.labels = obj.labels(ind, :);
            tmp_ids = cell2mat(obj.match_status.em_ids(ind));
            results.EM_IDs = em_ids(tmp_ids, 1);
            results.match_score = obj.match_status.scores;
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
            obj.A_em = sparse(obj.A_em);
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
        function [img, col0, AA] = overlapA(obj, ind, ratio, bg_white, use_em, zoomin)
            %merge all neurons' spatial components into one singal image
            if ~exist('use_em', 'var') || isempty(use_em)
                use_em = false;
            end
            if use_em
                A_ = obj.A_em;
            else
                A_ = obj.A;
            end
            if nargin<2 || isempty(ind)
                AA = A_;
                CC = obj.C;
                CC_raw = obj.C_raw;
                confidence = obj.match_status.confidence;
            else
                AA = A_(:, ind);
                CC = obj.C(ind,:);
                CC_raw = obj.C_raw(ind,:);
                confidence = obj.match_status.confidence(ind);
            end
            snr = var(CC, 0, 2) ./var(CC_raw-CC, 0, 2);
            snr(isnan(snr)) = 0;
            
            if ~exist('ratio', 'var')||isempty(ratio)
                ratio = 0.1;
            end
            if ~exist('bg_white', 'var') || isempty(bg_white)
                bg_white = false;
            end
            d3 = obj.options.d3;
            
            %% normalize the amplitude of each neuron
            v_norm = max(AA, [], 1);
            AA = bsxfun(@times, AA, 1./v_norm); %sqrt(sum(AA.^2,1)));
            AA(bsxfun(@lt, AA, max(AA, [], 1)*ratio)) = 0;
            %             AA = bsxfun(@times, AA, sqrt(snr)');
            %             AA = bsxfun(@times, AA, (snr.^0.5)');
            [d, K] = size(AA);
            
            col = prism(K);
            col = col(randperm(K),:);
            if bg_white
                col = 1-col;
            end
            
            col0 = col;
            img = zeros(d, 3);
            img = AA * col;
            temp = img/max(img(:))*(2^16);
            temp = uint16(temp);
            if bg_white
                temp = 2^16-temp;
            elseif ~isempty(obj.spatial_range)
                temp(~obj.spatial_range, :) = 0;
            end
            
            % rotate and zoomin
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
            img_A_em = obj.overlap_neurons(ind_pair, 'A_em');
            obj.showImage(img_A_em);
            set(gcf, 'name', 'EM footprints', 'position', tmp_pos0+[2*tmp_pos0(3),0, 0 0]);
            if save_figs
                saveas(gcf, sprintf('pair_%d_%d_A_em.pdf', ind_pair(1), ind_pair(2)));
            end
            
            %% temporal
            figure('papersize', [14, 4], 'name', 'temporal traces');
            tmp_pos0 = get(gcf, 'position');
            tmp_pos0(1:2) = [800, 0];
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
            figure(gcf); 
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
        function flag = check_existence(obj, target_ids)
            em_ids = cell2mat(obj.match_status.em_ids);
            flag = false(size(target_ids)); 
            for m=1:length(target_ids)
               flag(m) = any(em_ids==target_ids(m));  
            end
        end
        
        %% find EM boundary
        function em_ranges = get_em_boundary(obj)
            if isempty(neuron.spatial_range)
                return;
            end
            %% compute boundary
            ind = obj.reshape(obj.spatial_range, 3);
            em_ranges = cell(1, obj.options.d3);
            for m=1:3
                [y, x] = find(ind(:, :, m));
                k = convhull(x, y);
                em_ranges{m} = [x(k), y(k)];
            end
        end
        
        %% remove neurons in the blacklist
        function remove_black_list(obj, black_list)
            K = size(obj.A, 2);
            ind = false(K, 1);
            for m=1:K
                if obj.match_status.status(m)==1 && ...
                        any(black_list ==obj.match_status.em_ids{m})
                    ind(m) = true;
                end
            end
            obj.delete(ind, true);
        end
        
        %% select pairs to merge
        function [correlated_pairs, A_sim, C_sim] = find_correlated_pairs(obj, corr_thr)
            if ~exist('corr_thr', 'var') || isempty(corr_thr)
                corr_thr = [0.1, 0.3];
            elseif length(corr_thr)==1
                corr_thr = [0.1, corr_thr];
            end
            fprintf('minimum spatial correlation: %.2f\n', corr_thr(1)); 
            fprintf('minimum temporal correlation: %.2f\n', corr_thr(2)); 
            
            % corr_thr: [min_spatial_similarit, min_temporal_similarity]
            A_sim = cosine_similarity(obj.A);
            C_sim = cosine_similarity(obj.C');
            [ii, jj] = find((A_sim>=corr_thr(1)) .* (C_sim>=corr_thr(2)));
            ind = (ii<jj);
            correlated_pairs = [ii(ind), jj(ind)];
        end
        
        %% extract neuron boundaries
        function em_contours = find_em_contours(obj, ind)
            if ~exist('ind', 'var') || isempty(ind)
                A_ = obj.A_em;
            else
                A_ = obj.A_em(:, ind);
            end
            K = size(A_, 2);
            
            d3 = obj.options.d3;
            em_contours = cell(K, d3);
            for m=1:K
                % for each neuron
                ai_em = obj.reshape(A_(:, m), 3);
                for n=1:d3
                    % for each plane
                    img = ai_em(:, :, n);
                    B = bwboundaries(img);
                    for k=1:length(B)
                        if size(B{k},1)<6
                            break;
                        end
                        frame_length = 5;
                        
                        % for each connected components
                        smoothx = sgolayfilt(B{k}(:, 2), 2, frame_length);
                        smoothy = sgolayfilt(B{k}(:, 1) , 2, frame_length);
                        B{k} = [smoothy, smoothx];
                    end
                    em_contours{m, n} = B;
                end
            end
        end
        
        %% rotate image and zoom in
        function [img, rot_info] = rotate_zoomin(obj, img, rot_info)
            img = obj.reshape(img, 3);
            if exist('rot_info', 'var')
                rot_angle = rot_info.rot_angle;
                crange = rot_info.crange;
                rrange = rot_info.rrange;
            else
                ind_nnz = obj.reshape(obj.spatial_range, 3);
                [rot_angle, crange, rrange] = find_rotation(ind_nnz);
                rot_info = struct('rot_angle', rot_angle, ...
                    'crange', crange, ...
                    'rrange', rrange);
            end
            img = imrotate(img, rot_angle);
            img = img(rrange(1):rrange(2), crange(1):crange(2), :, :);
        end
        
        %% compute tuning curves for all neurons
        function compute_tuning_curve(obj, stimuli)%, sig)
            bins_nan = isnan(stimuli);
            ori = stimuli(~bins_nan);
            
            x = reshape(unique(ori), [], 1);
            n = hist(ori, x);
            y = bsxfun(@eq, x, reshape(ori, 1, [])) * obj.S(:, ~bins_nan)';
            y = bsxfun(@times, y, reshape(1./n, [], 1));
            y(:, sum(y)==0) = 1;
            y = bsxfun(@times, y, 1./sum(y, 1));
            
            obj.tuning_curve = struct('x', x, 'y', y);
        end
        
        %% bootstrap tuning curves
        function boostrap_tuning_curve(obj, shuffle)
            K = size(obj.A, 2);
            if  isempty(obj.tuning_curve) || K~=size(obj.tuning_curve.y, 2)
                error('computing tuning curve first!');
            else
                tc = obj.tuning_curve;
            end
            if ~exist('shuffle', 'var') || isempty(shuffle)
                shuffle = [];
            end
            tc.shuffle = shuffle;
            tc.yfit = zeros(size(tc.y));
            tc.pars = zeros(5, K);
            tc.rss = zeros(1, K);
            tc.pvals = zeros(1, K);
            for m=1:K
                fprintf('|');
            end
            fprintf('\n');
            for m=1:K
                if isempty(shuffle)
                    [tc.yfit(:, m), tc.pars(:, m), tc.rss(m)] =...
                        fit_von_mises2(tc.x, tc.y(:, m));
                else
                    [tc.yfit(:, m), tc.pars(:, m), tc.rss(m), tc.pvals(m)] =...
                        bootstrap_von_mises2(tc.x, tc.y(:, m), shuffle);
                end
                fprintf('.');
            end
            fprintf('\n');
            obj.tuning_curve = tc;
        end
        %% post process spatial shapes
        function post_process_spatial(obj, with_threshold)
            if exist('with_threshold', 'var') && with_threshold
                % compute residual
                sn = std(obj.compute_residual(), 0, 2);
                thresh = sn * (2./sqrt(sum(obj.C.^2, 2)'));
                obj.A(obj.A<thresh) = 0;
            end
            
            
            K = size(obj.A, 2);
            for m=1:K
                % keep pixels connecting to the spatial mask
                ai = obj.reshape(obj.A(:, m), 3);
                ind = obj.reshape(obj.A_em(:, m)>0, 3);
                for n=1:3
                    ind = (imdilate(ind, strel('square', 2)).*ai>0);
                end
                
                % remove isolated pixels
                ind = (imfilter(double(ind), ones(3,3)) > 1);
                obj.A(:, m) = ai(:).*ind(:);
            end
        end
        
        %% construct in & out area
        [ind_in, ind_out] = construct_in_out(obj, ai_em);
        
        %% match neurons with the EM components that has the largest matching score
        function n_corrections = rematch(obj, Aem, segment_ids, min_rank)
            if iscell(Aem)
                Aem = cell2mat(Aem);
                segment_ids = cell2mat(segment_ids);
            end
            if ~exist('min_rank', 'var') || isempty(min_rank)  % only rematch components whose EM match is not the top min_rank candidates
                min_rank = 2;
            end
            % find the mismatches
            confidence = obj.match_status.confidence;
            ind = find(confidence < 1);
            k_mismatch = length(ind);
            
            % match it with the best EM component
            n_corrections = 0;
            for k=1:k_mismatch
                cell_id = ind(k);
                [tmp_scores, idx_sort] = sort(obj.scores(cell_id,:), 'descend');
                if confidence(cell_id) >= tmp_scores(1) / tmp_scores(min_rank)
                    continue;
                end
                n_corrections = n_corrections + 1;
                obj.match_status.em_ids{cell_id} = segment_ids(idx_sort(1));
                obj.match_status.confidence(cell_id) = tmp_scores(1) / tmp_scores(2);
                obj.match_status.scores(cell_id) = tmp_scores(1);
                obj.A_em(:, cell_id) = Aem(:, idx_sort(1));
            end
        end
        %% merge repeated extraction
        function merge_repeats(obj)
            em_ids = cell2mat(obj.match_status.em_ids);
            em_ids = reshape(em_ids, [], 1);
            [unique_ids, ~, ic] = unique(em_ids);
            tmp_count = accumarray(ic, 1);
            idx = find(tmp_count>1);    % find EM components that have repeated extraction
            if ~isempty(idx)
                for m=1:length(idx)
                    tmp_id = unique_ids(idx(m));
                    ind = find(em_ids==tmp_id);
                    wi = obj.A_em(:, ind(1));
                    tmp_A = obj.A(:, ind);
                    tmp_Craw = obj.C_raw(ind, :);
                    temp = sum(tmp_A,1) .*sum(obj.C(ind,:), 2)'; % use the largest neuron as initialization
                    [~, tmp_ind] = max(temp);
                    ai = tmp_A(:, tmp_ind);
                    for miter=1:10
                        %% estimate ci
                        ci = ((wi.*ai)'*tmp_A*tmp_Craw) / ((wi.*ai)'*ai);
                        
                        % estimate ai
                        ai = max(0, tmp_A*(tmp_Craw*ci') / (ci*ci'));  % ai >= 0
                    end
                    obj.A(:, ind(1)) = ai;
                    obj.C_raw(ind(1), :) = ci;
                    obj.C(ind(1), :) = ci;
                    obj.match_status.confidence(ind(2:end)) = 0;
                end
                obj.delete(obj.match_status.confidence==0);
            end
            
        end
        
        %% construct weights for fixing boundaries effect 
        function weights_boundaries = fix_boundary_effects(obj)
            ind_in = obj.reshape(obj.spatial_range, 3); 
            ind_edge = imdilate(ind_in, strel('disk', 3)) - ind_in; 
            temp = imfilter(ind_edge, ones(3)); 
            weights_boundaries = obj.reshape(temp, 1);       
        end
        
        %% run EASE automatically
        at_ease(obj, configs);
        
        %% manually delete neurons 
        function vij = calculate_score(obj, ai, pj)
            % check the spatial footprint
            ai = obj.reshape(ai, 3); 
            pj = obj.reshape(pj, 3);
            ai(ai<0) = 0;
            tmp_img = image_local_corr(ai, pj, 7);
            vij = corr(tmp_img(:).*(pj(:)>0), pj(:));
            
        end
        
        %% detailed check EM matches 
        function [idx_top, v_top] = check_top_matches(obj, Aem)
            K = size(obj.A, 2); 
            ntop = 5; 
            idx_top = zeros(ntop, K); 
            v_top = zeros(ntop, K); 
            for m=1:K 
                [~, temp] = sort(obj.scores(m,:), 'descend'); 
                idx = temp(1:ntop); 
                idx_top(:, m) = idx; 
                tmp_ai = obj.reshape(obj.A_corr(:, m), 3);
                for n=1:ntop
                    ai_em = obj.reshape(Aem(:, idx(n)), 3);
                    tmp_img = image_local_corr(tmp_ai, ai_em, 7);
                    %             tmp_corr = corr(tmp_img(:).*tmp_ai(:), tmp_img(:).*ai(:));
                    tmp_corr = corr(tmp_img(:).*(ai_em(:)>0), ai_em(:));
                    v_top(n, m) = tmp_corr; 
                end
            end 
            
        end 
        
      
    end
end














