function [dims, ind, vv] = compress_binary_matrix(X)
%% compress a binary matrix by replacing blocks of 4*4 with a uint16 number. 
%% input: 
%   X: d1*d2 matrix 

%% output: 
%   dims: [d1, d2], matrix dimension 
%   ind: k*1 vector, indices of nonzero blocks 
%   v: k*1 vector, value of each 4*4 block 

%% author: Pengcheng Zhou, Columbia University, 2018 
% zhoupc1988@gmail.com 

%% code 
[d1, d2]  = size(X); 
d1_block = ceil(d1/4); 
d2_block = ceil(d2/4);

kernel = reshape(2.^(15:-1:0), 4, 4); 
X_conv = conv2(double(X), kernel); 
X = X_conv(4:4:(d1_block*4), 4:4:(d2_block*4)); 
% % pad the matrix to make its row number and column numbers 
% if (d1_block*4~=d1) || (d2_block*4~=d2)
%     X = padarray(X, d1_block*4-d1, d2_block*4-d2); 
% end 
% 
% % convert the matrix to array with elements as 4*4 matrix 
% X = mat2cell(X, ones(d1_block,1)*4, ones(d2_block, 1)*4); 
% 
% % convert each 4*4 matrix into a uint16 number 
% X = cellfun(@(M) bi2de(reshape(M, 1, [])), X); 
% 
[ii, jj, vv] = find(X); 
ind = sub2ind([d1_block, d2_block], ii, jj); 
d = d1_block*d2_block; 
if d<2^8
    ind = uint8(ind); 
elseif d<2^16
    ind = uint16(ind); 
elseif d<2^32
    ind = uint32(ind); 
else
    ind = uint64(ind); 
end

vv = uint16(vv); 
dims = [d1, d2]; 
end 