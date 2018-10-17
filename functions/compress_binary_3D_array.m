function [dims, ind_all, v_all, nnz_all] = compress_binary_3D_array(X, use_parallel)
%% compress a 3D binary array
%% input:
%   X: d1*d2*d3 array

%% output:
%   dims: [d1, d2, d3], array dimension
%   ind_all: k*1 vector, indices of nonzero voxels
%   v_all: k*1 vector, value of each 4*4 block
%   nnz_all: number of non-zero elements in each slice

%% author: Pengcheng Zhou, Columbia University, 2018
% zhoupc1988@gmail.com

%% code
[d1, d2, d3] = size(X);
dims = [d1, d2, d3];
ind_all = cell(d3, 1);
v_all = cell(d3, 1);

% empty outpty
if d1*d2<2^8
    ind0 = uint8([]);
elseif d1*d2<2^16
    ind0 = uint16([]);
elseif d1*d2<2^32
    ind0 = uint16([]);
else
    ind0 = uint64([]);
end

% compress each slice independently
if ~exist('use_parallel', 'var') || ~use_parallel
    nnz_all = zeros(d3,1);
    for m=1:d3
        tmpX = X(:, :, m);
        if sum(tmpX(:)) == 0
            ind_all{m} = ind0;
            v_all{m} = uint16([]);
            continue;
        else
            [~, ind_all{m}, v_all{m}] = compress_binary_matrix(tmpX);
            nnz_all(m) = length(v_all{m});
        end
    end
else
    X = mat2cell(reshape(X, [], d3), d1*d2, ones(d3, 1));
    nnz_all = cell(d3,1);
    
    parfor m=1:d3
        tmpX = X{m};
        if sum(tmpX) == 0
            ind_all{m} = ind0;
            v_all{m} = uint16([]);
            nnz_all{m} = 0;
            continue;
        else
            [~, ind_all{m}, v_all{m}] = compress_binary_matrix(reshape(tmpX, d1, d2));
            nnz_all{m} = length(v_all{m});
        end
    end
    nnz_all = cell2mat(nnz_all);
end
ind_all = cell2mat(ind_all);
v_all = cell2mat(v_all);

max_nnz = max(nnz_all); 
if max_nnz<2^8
    nnz_all = uint8(nnz_all); 
elseif max_nnz<2^16 
    nnz_all = uint16(nnz_all); 
elseif max_nnz<2^32
    nnz_all = uint32(nnz_all); 
else
    nnz_all = uint64(nnz_all); 
end
