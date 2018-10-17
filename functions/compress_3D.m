function [dims, ind_all, v_all, nnz_all]  = compress_3D(ijk, dims)
%% compress 3D matrix 
%% inputs: 
%   ijk: n*3 array, [i, j, k] of nonzero voxels 
%   dims: [d1, d2, d3], matrix dimension 

%% outputs: 
%   dims: [d1, d2, d3], array dimension
%   ind_all: k*1 vector, indices of nonzero voxels
%   v_all: k*1 vector, value of each 4*4 block
%   nnz_all: number of non-zero elements in each slice

%% author: Pengcheng Zhou, Columbia University, 2018
% zhoupc1988@gmail.com

%% code 
% input arguments 
ii = ijk(:,1); 
jj = ijk(:,2); 
kk = ijk(:,3); 

d1 = dims(1); 
d2 = dims(2); 
d3 = dims(3); 

d1_new = ceil(d1/4); 
d2_new = ceil(d2/4); 

% compute the value of each pixel in its belong block 
i_in = mod(ii-1, 4) + 1;  % index within each 4*4 block  
j_in = mod(jj-1, 4) + 1; 
temp = sub2ind([4, 4], i_in, j_in); 
v_all = 2.^(temp-1); 

% get the indices of those nonzero blocks 
i_block = (ii-i_in) / 4+1;  % block index 
j_block = (jj-j_in) / 4+1; 
ind_all = sub2ind([d1_new, d2_new], i_block, j_block); 

% order neurons based on [kvalue, ind]
temp = (d1_new*d2_new) * kk + ind_all; 
[~, idx]  = sort(temp); 
ind_all = ind_all(idx);
v_all = v_all(idx); 
kk = kk(idx); 

% save results 
N = length(kk); 
nnz_all = zeros(d3, 1); 

idx_z = 0;
idx_ind = nan; 
m_last = 0; 
for m=1:N 
    if kk(m)~=idx_z  % switch z planes
        idx_z = kk(m);
        idx_ind = nan; 
    end
    
    if ind_all(m) ~= idx_ind
        idx_ind = ind_all(m);
        m_last = m; 
        nnz_all(idx_z) = nnz_all(idx_z) + 1;
    else
        v_all(m_last) = v_all(m_last) + v_all(m);
        ind_all(m) = nan; 
    end
end 

ind = isnan(ind_all); 
ind_all(ind) = []; 
v_all(ind) = []; 

%% choose the right format for saving space 
v_all = uint16(v_all); 

temp = [8, 16, 32, 64]; 
idx = find(temp>log2(double(max(ind_all))), 1); 
eval(sprintf('ind_all = uint%d(ind_all); ', temp(idx))); 

idx = find(temp>log2(double(max(nnz_all))), 1); 

eval(sprintf('nnz_all = uint%d(nnz_all); ', temp(idx))); 
























