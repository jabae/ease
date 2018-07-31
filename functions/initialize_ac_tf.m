function [ai, ci] = initialize_ac_tf(Yin, Yout, wi, ci, rank_out)
%% initialize a neuron's spatial & temporal component given its spatial support 
%{
    This initialization algorithm uses the spatial mask of the neuron and
    remove the contaminating signals form its neighbors and background. 
%}

%% inputs: 
%{  
    Yin: d_in * T matrix, imaging data within the spaital support 
    Yout: d_out * T matrix, imaging data surrounding the spatial support 
    wi:  d_in * 1 vector, weights on each pixels in the regression problem 
    ci: 1*T vecotr, initialized temporal components 
    rank_out: integer, the truncated rank of Yout 
%}

%% outputs: 
%{
    ai: d_in * 1, spatial mask of the initialized neuron 
    ci: 1 * T, tmeporal trace of the initialized neuron 
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

%% code 
[din, ~] = size(Yin); 
if ~exist('wi', 'var') || isempty(wi)
    wi = ones(din, 1); 
else
    wi = reshape(wi, [], 1);
end
if ~exist('rank_out', 'var') || isempty(rank_out)
    rank_out = 20;
end
[~, ~, V] = svdsecon(Yout, rank_out); 
Vp = null(V');   % orthogonal basis for the null space of V 

% 2nd order difference of V and Vp 
Dv = diff(V, 2, 1);

% project Yin to the null space of V
Y_p = Yin - (Yin*V) *V';

%% iteratively update ai an dalpha 
ci_p = ci - (ci*V) * V';
for miter=1:4
    % estimate ai
    ai = max(0, Y_p*ci_p' / (ci_p*ci_p'));  % ai >= 0
    
    %% estimate ci_p
    ci_p = ((wi.*ai)'*Y_p) / ((wi.*ai)'*ai);
    
end

%% find ci by smoothing its trend 
Dci_p = diff(ci_p, 2)';

cvx_begin quiet
variable w(rank_out)
minimize(norm(Dv*w+Dci_p, 1))
cvx_end

ci = w'*V' + ci_p;

%% normalize results 
ci = ci * norm(ai,2); 
ai = ai / norm(ai, 2);


