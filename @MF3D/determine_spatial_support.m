function A_mask = determine_spatial_support(obj)
%% determine the support to all spatial components 
[d, K] = size(obj.A);
A_ = mat2cell(obj.A, d, ones(1, K)); 
IND = cell(1, K); 
expandCore = strel(ones(3, 3, 3)); 
d1 = obj.options.d1; 
d2 = obj.options.d2; 
d3 = obj.options.d3; 

parfor m = 1:K
    A_temp = imdilate(reshape(full(A_{m}),d1,d2,d3),expandCore);
    IND{m} = A_temp(:)>0;
end

A_mask = cell2mat(IND); 
obj.A_mask = A_mask; 