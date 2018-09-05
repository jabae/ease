function [ai, ci] = initialize_ac_seminmf(Y, wi, ci, maxIter)
%% initialize a neuron's spatial & temporal component by running rank-1 NMF 
%{
    This initialization algorithm models the data Y as Y = ai*ci subject to
    ai >= 0;
%}

%% inputs: 
%{  
    Y: d_in * T matrix, imaging data within the spaital support 
    wi:  d_in * 1 vector, weights on each pixels in the regression problem 
    ci: 1*T vecotr, initialized temporal components 
    maxIter: integer, maximum number of iteration 
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
[din, ~] = size(Y); 
if ~exist('wi', 'var') || isempty(wi)
    wi = ones(din, 1); 
else
    wi = reshape(wi, [], 1);
end
if ~exist('maxIter', 'var') 
    maxIter = 10; 
end
%% iteratively update ai an dalpha 
for miter=1:maxIter
    % estimate ai
    ai = max(0, Y*ci' / (ci*ci'));  % ai >= 0
    
    %% estimate ci
    ci = ((wi.*ai)'*Y) / ((wi.*ai)'*ai);
end


%% normalize results 
ci = ci * norm(ai,2); 
ai = ai / norm(ai, 2);


