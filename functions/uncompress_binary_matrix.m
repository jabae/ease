function X = uncompress_binary_matrix(dims, ind, vv)
%% compress a binary matrix by replacing blocks of 4*4 with a uint16 number. 

%% input: 
%   dims: [d1, d2], matrix dimension 
%   ind: k*1 vector, indices of nonzero blocks 
%   v: k*1 vector, value of each 4*4 block 

%% output: 
%   X: d1*d2 matrix 

%% author: Pengcheng Zhou, Columbia University, 2018 
% zhoupc1988@gmail.com 

%% code 
d1 = dims(1); 
d2 = dims(2);  
d1_block = ceil(d1/4); 
d2_block = ceil(d2/4);
N = length(ind); 

% create a cell array for indices of non-zero values 
ii_all = cell(N, 1); 
jj_all = cell(N, 1); 
[ii, jj] = ind2sub([d1_block, d2_block], ind); 

for m=1:N 
    [tmp_ii, tmp_jj] = ind2sub([4, 4], find(de2bi(vv(m))')); 
    ii_all{m} = ii(m)*4+tmp_ii-4; 
    jj_all{m} = jj(m)*4+tmp_jj-4; 
end 
X = sparse(cell2mat(ii_all), cell2mat(jj_all), 1, d1_block*4, d2_block*4); 
X = X(1:d1, 1:d2); 