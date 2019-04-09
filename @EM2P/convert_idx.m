function [idx_new, dims_new] = convert_idx(obj, idx) 
%% convert the voxel indices in the stack space into the indices in calcium imaging planes 
%{
	indices in the stack space is the indices in a 3D array of dimension = obj.dims_2p. 
    indices in the stack space is the indices in a 3D array of dimension =
    [obj.d1, obj.d2, num_scans * num_planes]
    The function has to deal with following things: 
    1. cropping (the new 3D array is a cropped area in x-y plane and
    subsampling in z axis)
    2. motion in xy planes are different for different z 
%}

%% inputs
%{
	obj: type; description
	idx: n*1 vector; indices in 3D array (dim = obj.dims_stack)
%}

%% outputs
%{
%}

%% Author
%{
	Pengcheng Zhou 
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License 
%}

%% stack space 
dims_2p = obj.dims_stack;      % stack dimension 
FOV_stack_ = obj.FOV_stack;     % FOV in 2p stack 
[idx_r, idx_c, idx_z] = ind2sub(dims_2p, idx); 
ssub = obj.ssub; 

% new z indices 
zvals = obj.video_zvals_updated; % zvalues of all planes 
zidx = bsxfun(@plus, (0:(obj.num_scans-1))'*obj.num_slices, 1:obj.num_slices); 
temp = zeros(dims_2p(3),1); 
temp(zvals(:)) = zidx(:);
new_idx_z = temp(idx_z); 

% new xy indices 
ii0 = zeros(dims_2p(3), 1);      % spatial shifts in y direction 
jj0 = zeros(dims_2p(3), 1);     % spatial shifts in x direction 
ii0(zvals(:)) = obj.em_shifts.ii(:) + FOV_stack_(1); 
jj0(zvals(:)) = obj.em_shifts.jj(:) + FOV_stack_(3); 
idx_r = idx_r - ii0(idx_z) + 1; 
idx_c = idx_c - jj0(idx_z) + 1;

new_idx_r = floor(idx_r / ssub);
new_idx_c = floor(idx_c / ssub); 

%% new space 
dims_new = [obj.d1, obj.d2, obj.num_scans*obj.num_slices];   % new dimension 
idx_new = sub2ind(dims_new, new_idx_r, new_idx_c, new_idx_z); 

