function [ai, ci] = initialize_ac_tf(Yin, Yout, wi, ci, rank_out, maxIter)
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
[din, ~] = size(Yin);
dout = size(Yout, 1);
if ~exist('wi', 'var') || isempty(wi)
    wi = ones(din, 1);
else
    wi = reshape(wi, [], 1);
end
if ~exist('rank_out', 'var') || isempty(rank_out)
    rank_out = 20;
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 10;
end
[~, ~, V] = svdsecon(randn(rank_out*5, dout) * Yout, rank_out);
% Vp = null(V');   % orthogonal basis for the null space of V

% 2nd order difference of V and Vp
Dv = diff(V, 2, 1);

% project Yin to the null space of V
Y_p = Yin - (Yin*V) *V';

%% iteratively update ai an dalpha
ci_p = ci - (ci*V) * V';
for miter=1:maxIter
    % estimate ai
    ai = max(0, Y_p*ci_p' / (ci_p*ci_p'));  % ai >= 0
    
    %% estimate ci_p
    ci_p = ((wi.*ai)'*Y_p) / ((wi.*ai)'*ai);
end

%% find ci by smoothing its trend
Dci_p = diff(ci_p, 2)';

% cvx_begin quiet
% variable w(rank_out)
% minimize(norm(Dv*w+Dci_p, 1))
% cvx_end

w = solve_w(Dv, Dci_p);
ci = w'*V' + ci_p;

%% normalize results
ci = ci * norm(ai,2);
ai = ai / norm(ai, 2);
end

function [w, iter]= solve_w(X, y)
% solve problem minimize norm(X*w+y, 1) given X & y
% use gradient descent algorithm and back tracking line search method to determine step size

p = size(X, 2);     % dimension of w
w = zeros(p, 1);

t = 1;              % initialize step size
beta = 0.9;
tol = 1e-10;
maxIter = 100;

for iter=1:maxIter
    z = X*w + y;
    g = X' * sign(z);   % gradient
    v0 = norm(z, 1);
    % back tracking line search
    while norm(z - t*X*g, 1) > v0 - t*norm(g, 2)^2/2
        t = beta*t;
    end
    
    if t*norm(g, 2)^2/v0 < tol
        break;
    end
    w = w - t *g;
end
end