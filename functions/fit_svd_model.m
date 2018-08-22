function [b, f, b0] = fit_svd_model(Y, nb, A, C, b_old, f_old, thresh_outlier, sn, ind_patch)
%% fit a patched data with SVD%% what does this function do 
%{
    details of this function 
%}

%% inputs: 
%{
%}

%% outputs: 
%{
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

%% code 

Ymean = mean(Y,2); 
if ~exist('A', 'var') || isempty(A)
    [d, T] = size(Y); 
    A = ones(d,1); 
    C = zeros(1, T); 
elseif issparse(A)
    A = full(A); 
end
Cmean = mean(C, 2); 
Y = bsxfun(@minus, double(Y), Ymean); 
C = bsxfun(@minus, C, Cmean);

Bf = Y - A*C;
if ~exist('ind_patch', 'var')
    ind_patch = true(size(A,1), 1); 
end
if exist('thresho_outlier', 'var')
    Bf_old = b_old*f_old;
    tmp_Bf = Bf(ind_patch, :);
    ind_outlier = bsxfun(@gt, tmp_Bf, bsxfun(@plus, Bf_old, thresh_outlier*reshape(sn, [], 1)));
    tmp_Bf(ind_outlier) = Bf_old(ind_outlier);
    Bf(ind_patch, :) = tmp_Bf;
    clear tmp_Bf Bf_old ind_outlier;
end

Bf_center = bsxfun(@minus, Bf, mean(Bf, 2)); 
tsub = max(1, floor(size(Bf_center, 2)/1000)); 
[u, ~, ~] = svdsecon(Bf_center(:, 1:tsub:end), nb); 
% [u, s, v] = rsvd(bsxfun(@minus, Bf, mean(Bf, 2)), nb); 
b = u(ind_patch, :); 
f = u'*Bf_center; 
b0 = Ymean(ind_patch) - A(ind_patch,:)*Cmean-b*mean(f, 2);

