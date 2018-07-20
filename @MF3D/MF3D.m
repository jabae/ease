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
        %quality control
        ids;        % unique identifier for each neuron
        match_status = struct('status', [], 'em_ids', []);  % the status of neuron matching. It's a struct variable with two fields:
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
            
            obj.A(:, ind) = [];
            obj.A_mask(:, ind) = [];
            if ~isempty(obj.A_corr)
                obj.A_corr(:, ind) = [];
            end
            obj.C(ind, :) = [];
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
            end
            if size(obj.scores,1) == size(obj.A, 2)
                obj.scores(ind, :) = [];
            end
            % save the log
        end
        
        
        %% view neurons
        function showNeurons(obj, A_, ind)
            if ~exist('A_', 'var') || isempty(A_)
                A_ = obj.A;
            end
            if ~exist('ind', 'var') || isempty(ind)
                ind = 1:size(A_, 2);
            end
            
            figure;
            d3 = obj.options.d3;
            for m=1:length(ind)
                ai = obj.reshape(A_(:, ind(m)), 3);
                ai_max = max(ai(:));
                for n=d3:-1:1
                    subplot(d3, 1, n);
                    imagesc(ai(:, :, n), [0, ai_max*0.8]);
                    axis equal off tight;
                end
                title(sprintf('neuron %d', m)); %ind(m)));
                pause(.1);
            end
        end
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
        showDemixing(obj, Y, min_max, col_map, avi_nm, t_pause, ind_neuron);
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
            em_mask = (sum(Aem, 2)<=0);
            
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
        ind_voxels_em = initialize_em(obj, Aem, Y, Kmax, options);
        ind_voxels_em = initialize_all(obj, Aem, Y, Kmax, options);
        
        %% deconvolve all temporal components
        C_ = deconvTemporal(obj, use_parallel, method_noise)
        
        %% run HALS
        function hals(obj, Y)
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
            Y = obj.reshape(Y,1);
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
            Yres = bsxfun(@minus, obj.reshape(Y,1) - obj.b*obj.f-obj.A*obj.C, obj.b0);
        end
        
        function rss = compute_rss(obj, Y)
            temp = bsxfun(@minus, obj.reshape(Y,1) - obj.b*obj.f-obj.A*obj.C, obj.b0);
            rss = sum(temp(:).^2);
        end
        %% decorrelate neuron traces
        C_ = decorrTemporal(obj, wd)
        
        %% show image
        function showImage(obj, ai)
            if numel(ai) == 1
                ai = obj.reshape(obj.A(:, ai), 3);
            else
                ai = obj.reshape(ai, 3);
            end
            [d1, d2, d3, dc] = size(ai);
            ai_max = max(ai(:));
            ai_min = quantile(ai(:), 0.1);
            figure('papersize', [d2, d1*d3]*5/d2);
            init_fig;
            for m=1:d3
                subplot(d3, 1, m);
                if dc==3
                    imagesc(squeeze(ai(:, :, m, :)), [ai_min, ai_max + 1*(ai_max==ai_min)]);
                else
                    imagesc(ai(:, :, m), [ai_min, ai_max+1*(ai_max==ai_min)]);
                end
                axis equal off tight;
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
    end
    
    
end














