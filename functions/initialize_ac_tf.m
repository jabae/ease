function [ai, ci] = initialize_ac_tf(Yin, Yout, wi, ci, bg_rank)
% wi weights each pixel according to its EM components 
if ~exist('r', 'var') || isempty(r)
    bg_rank = 10;
end
[~, ~, V] = svdsecon(Yout, bg_rank);
T = length(ci);
Pvp = eye(T) - V*V';    % projection matrix for the null space of V
wi = reshape(wi, [], 1); 

% differential matrix
Dv = diff(V, 2, 1);

% project Yin to the null space of V
Y_p = Yin*Pvp;

%% run iterations
for miter=1:2
    % estimate ai
    ci_p = ci*Pvp;  %project ci to the null space of V
    ai = max(0, Y_p*ci_p' / (ci_p*ci_p'));  % ai >= 0
    
    %% estimate ci to minimize its trend
    yi_p = ((wi.*ai)'*Y_p) / ((wi.*ai)'*ai);
    
    % find a smooth ci
    Dyi = diff(yi_p, 2)';
    
    cvx_begin quiet
        variable w(bg_rank)
        minimize(norm(Dv*w+Dyi, 1))
    cvx_end
    
    ci = w'*V' + yi_p; 
end
ci = ci * norm(ai,2); 
ai = ai / norm(ai, 2); 