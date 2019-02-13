function confidence = evaluate_matching_confidence(obj, Aem, Yr)
%% evaluate the maching performance between EM segments and 2P neurons . 
d3 = obj.options.d3; 
% check the spatial mask
if iscell(Aem)
    Aem = obj.convert_matrix(Aem);
end
em_mask = (sum(Aem, 2)<=0);   % constrain to the area within the EM volume
Aem(em_mask, :) = [];
n = size(Aem, 1);
Aem_sum = sum(Aem, 1);
Aem_std = sqrt(sum(Aem.^2,1) - (Aem_sum.^2/n));

%% compute correlation
if ~exist('Yr', 'var') || isempty(Yr)
    if isempty(obj.frame_range)
        Yr = evalin('base', 'Y_cnmf');
    else
        temp = obj.frame_range;
        Yr = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
    end
end
Yres = obj.reshape(Yr, 1) - obj.A*obj.C - obj.b*obj.f;
Yres = bsxfun(@minus, Yres, mean(Yres, 2));
var_Yres = sum(Yres.^2, 2);
tmpA_corr = zeros(size(obj.A));
A_ = obj.A;

for m=1:size(obj.A, 2)
    ci = obj.C(m, :);
    ci = ci - mean(ci);
    ai = obj.A(:, m);
    temp = (Yres*ci' + ai*(ci*ci')); 
    tmpA_corr(:, m) = temp ./sqrt(var_Yres+(ai.^2)*(ci*ci'))/sqrt(ci*ci'); % approximate the variance
    
    % estimate a unconstrained ai
    ai = obj.reshape(obj.A_mask(:, m), 3);
    ai_mask = imdilate(repmat(sum(ai, 3)>0, [1, 1, d3]), strel('square', 3));
    %     ai_mask = imdilate(ai>0, strel('square', 3));
    A_(:, m) = temp.* reshape(ai_mask>0, [],1);
end
obj.A_corr = tmpA_corr;

%% compute matching score
A_(em_mask, :) = [];
A_p = A_; %positive elements 
A_p(A_p<0) = 0; 
% A_(A_<=0) = 0;
P_ = A_; %obj.A_mask .* tmpA_corr;
% P_(em_mask, :) = [];
P_sum = sum(P_, 1);
P_std = sqrt(sum(P_.^2, 1) - (P_sum.^2/n));

temp1 = bsxfun(@times, A_'*(Aem>0), 1./sum(A_p,1)'); % explained signal with different masks
temp2 = P_'*Aem - P_sum'*(Aem_sum/n);  % a faster way of computing correlation coefficients
temp2 = bsxfun(@times, temp2, 1./P_std');
temp2 = bsxfun(@times, temp2, 1./Aem_std);

temp = temp1 .* temp2;
% temp = temp2; 
temp(isnan(temp)) = 0;
obj.scores = sparse(temp);

%% compute matching confidence
K = size(A_, 2);
confidence = zeros(1, K);
best_scores = zeros(1, K);
for m=1:K
    em_id = obj.match_status.em_ids{m};
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