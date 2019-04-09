function confidence = evaluate_matching_confidence(obj, Yr, Aem, segment_ids)
%% evaluate the maching performance between EM segments and 2P neurons . 
% check the spatial mask
if iscell(Aem)
    Aem = cell2mat(Aem);
    segment_ids = cell2mat(segment_ids); 
end
em_mask = (sum(Aem, 2)<=0);   % constrain to the area within the EM volume
Aem(em_mask, :) = [];
n = size(Aem, 1);
Aem_sum = sum(Aem, 1);
Aem_std = sqrt(sum(Aem.^2,1) - (Aem_sum.^2/n));


%% get residual 
Yres = obj.compute_residual(Yr); 
Yres = bsxfun(@minus, Yres, mean(Yres, 2)); 
var_Yres = sum(Yres.^2, 2);
var_Yres(var_Yres==0) = inf; 

%% get unconstrained estimation of A  and estimate correlation between each ci & Y
C_ = obj.C; 
C_ = bsxfun(@minus, C_, mean(C_, 2)); 
A_ = obj.A + (Yres*C_')/sum(C_.^2, 2)'; 
% variance of each neuron 
var_neurons = var_Yres + bsxfun(@times, obj.A.^2, sum(C_.^2, 2)'); 

obj.A_corr = bsxfun(@times, A_./sqrt(var_neurons), sqrt(sum(C_.^2, 2)')); 

%% compute matching score
A_(em_mask, :) = [];

A_sum = sum(A_, 1);
A_std = sqrt(sum(A_.^2, 1) - (A_sum.^2/n));

% correlation with ai 
temp2 = A_'*Aem - A_sum'*(Aem_sum/n);  % a faster way of computing correlation coefficients
temp2 = bsxfun(@times, temp2, 1./A_std');
corr_aipj = bsxfun(@times, temp2, 1./Aem_std);

% explained signal 
A_(A_<0) = 0; 
explained_signal = bsxfun(@times, A_, 1./sum(A_,1))' * (Aem>0); 
% sparsity of pj 
% sparsity = 

temp = (explained_signal) .* corr_aipj;
% temp =  corr_aipj;
temp(isnan(temp)) = 0;
obj.scores = sparse(temp);

%% compute matching confidence
K = size(A_, 2);
confidence = zeros(1, K);
best_scores = zeros(1, K);
for m=1:K
    em_id = find(segment_ids==obj.match_status.em_ids{m}, 1, 'first');
    temp = obj.scores(m, :);
    v_select = temp(em_id);
    temp(em_id) = -inf;
    v_others = max(temp);
    confidence(m) = v_select / v_others;
    best_scores(m)= v_select;
end
obj.match_status.confidence = confidence;
obj.match_status.scores = best_scores;
end

