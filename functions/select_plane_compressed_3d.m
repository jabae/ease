function [ind, vv, dims] = select_plane_compressed_3d(z, ind_all, v_all, nnz_all, dims)
%% select one plane from the compressed 3d data 

%% inputs: 
%   z: scalar or array, list of planes to be selected 
%   ind_all: indices of nonzero blocks 
%   v_all: values of those 4*4 blocks 
%   nnz_all: number of nnz values for each plane 

%% outputs: 
%   ind: indices of nonzero blocks 
%   vv: values of those 4*4 blocks 
%   dims: 2D dimensions 

%% code 

temp = cumsum(double(nnz_all)); 
if z==1
    idx = 1:temp(1); 
else
    idx = (temp(z-1)+1): temp(z);
end 

ind = ind_all(idx); 
vv = v_all(idx); 
dims = dims(1:2); 

