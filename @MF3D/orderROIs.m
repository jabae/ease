function [srt] = orderROIs(obj, srt)
% srt: sorting order
if nargin<2
    srt='srt';
end

if ischar(srt)
    switch lower(srt)
        case 'mean'
            if obj.options.deconv_flag
                temp = mean(obj.C,2)'.*sum(obj.A);
            else
                temp = mean(obj.C_raw.*(obj.C_raw>0),2)'.*sum(obj.A);
            end
            [~, srt] = sort(temp, 'descend');
        case 'sparsity_spatial'
            temp = sqrt(sum(obj.A.^2, 1))./sum(abs(obj.A), 1);
            [~, srt] = sort(temp);
        case 'sparsity_temporal'
            temp = sqrt(sum(obj.C_raw.^2, 2))./sum(abs(obj.C_raw), 2);
            [~, srt] = sort(temp, 'descend');
        case 'pnr'
            pnrs = max(obj.C, [], 2)./std(obj.C_raw-obj.C, 0, 2);
            [~, srt] = sort(pnrs, 'descend');
        case 'skewness'
            l3l2 = skewness(obj.C, 0, 2);
            [~, srt] = sort(l3l2, 'descend');
        case 'l3'
            C_ = bsxfun(@minus, obj.C_raw, mean(obj.C_raw,2));
            l3 = sum(obj.A.^3,1) .*sum(C_.^3, 2)';
            [~, srt] = sort(l3, 'descend');
        case 'temporal_cluster'
            fprintf('use the temporal correlations to order neurons\n');
            dd = pdist(obj.C_raw, 'cosine');
            tree = linkage(dd, 'complete');
            srt = optimalleaforder(tree, dd);
        case 'spatial_cluster'
            fprintf('use the spatial correlations to order neurons\n');
            A_ = bsxfun(@times, obj.A, 1./sqrt(sum(obj.A.^2, 1)));
            temp = 1-A_' * A_;
            dd = temp(tril(true(size(temp)), -1));
            dd = reshape(dd, 1, []);
            tree = linkage(dd, 'complete');
            srt = optimalleaforder(tree, dd);
        case 'confidence'
            fprintf('use the confidence score to order neurons\n');
            [~, srt] = sort(obj.match_status.confidence, 'descend');
             case 'labels'
            fprintf('use the labels to order neurons\n');
            [~, srt] = sort(obj.labels);
            
        case 'match_scores'
            fprintf('use the matching score to order neurons\n');
            [~, srt] = sort(obj.match_status.scores, 'descend');
        otherwise
            fprintf('use snr to order neurons\n');
            snrs = var(obj.C, 0, 2)./var(obj.C_raw-obj.C, 0, 2);
            [~, srt] = sort(snrs, 'descend');
    end
    
end
obj.A = obj.A(:, srt);
obj.C = obj.C(srt, :);

try
    obj.A_em = obj.A_em(:, srt);
    obj.C_raw = obj.C_raw(srt,:);
    obj.S = obj.S(srt, :);
    obj.ids = obj.ids(srt);
    obj.labels = obj.labels(srt);
    obj.match_status.status = obj.match_status.status(srt);
    obj.match_status.em_ids = obj.match_status.em_ids(srt);
    if ~isempty(obj.scores)
        obj.scores = obj.scores(srt, :);
    end
    if ~isempty(obj.match_status.confidence)
        obj.match_status.confidence = obj.match_status.confidence(srt);
    end
    if ~isempty(obj.match_status.confidence)
        obj.match_status.scores = obj.match_status.scores(srt);
    end
    obj.A_corr = obj.A_corr(:, srt);
    
end
end
