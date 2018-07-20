function [subs_2p, V] = mesh2volume_fast(vertices, faces, options)

%% convert triangular meshes into volume pixles.
%% inputs:
%   vertices: M*3 matrix, locations of vertices (x, y, z) 
%   faces: N*3 matrix, indices of vertices formulating a triangle 
%   options: 
%       voxel_em: em resolution for voxelization 
%       dims_2p: dimension of 2p imaging (d1, d2, d3) 
%       range_2p: spatial range of 2p imaging (l1, l2, l3) 
%       A, offset: 3*3 matrix & 1*3 vector, transform em coordinates to 2p
%           coordinates Y_em*A + offset
%       scale_factor: double, scale the EM coordinates because of some
%           error in the preprocessing step 
%% outputs:
%   subs_2p:  k*3 matrix, the indices for the voxels in 2p space 
%   V:        3D matrix of the voxelized component 

%% Author: Pengcheng Zhou, Columbia University, 2017

%% parameters 
voxel_em = options.voxel_em*4;    %em resolution
dims_2p = options.dims_2p; 
range_2p = options.range_2p; 
voxels_2p = range_2p./dims_2p;    %2p resolution
res_factor = 2;     % use a higher resolution before downsampling 
A = options.A; 
offset = options.offset; 
scale_factor = options.scale_factor; 

if isempty(faces)
    subs_2p = zeros(0, 3); 
    V = 0; 
    return; 
end 
%% voxelize mesh surfaces 
nz_voxels = cell(1, length(vertices));
parfor seg_id = 1:length(vertices)
    % convert the coordinates from EM space to 2p space 
    tmp_vert = bsxfun(@plus, vertices{seg_id} * scale_factor * A, offset);

    % the starting point 
    vert_0 = min(tmp_vert);
    vert_range = max(tmp_vert)-vert_0+1;
    
    % voxelize 
    FV = struct('vertices', bsxfun(@minus, tmp_vert, vert_0), ...
        'faces', faces{seg_id} + 1);
    tmp_V = polygon2voxel(FV, round(vert_range./voxel_em), 'auto', false);
    
    %% determine the xyz locations of nonzero voxels
    ind = find(tmp_V);
    [x, y, z] = ind2sub(size(tmp_V), ind);
    temp = bsxfun(@times, [x, y, z]-1, voxel_em);
    nz_voxels{seg_id} = bsxfun(@plus, temp', vert_0');
end

%% create a small volume for filling the empty holes
% xyz positions of all nonzero voxels 
subs = round(bsxfun(@times, cell2mat(nz_voxels), res_factor./voxels_2p')') + 1;
volume_size = range(subs)+1;
subs_0 = min(subs, [],1); 
idx_new = bsxfun(@minus, subs, subs_0)+1;
% get the indices of all nonzero voxels in the zoom-in space 
ind_new = unique(sub2ind(volume_size, idx_new(:,1), idx_new(:,2), ...
    idx_new(:,3)));
V = false(volume_size);
V(ind_new) = true;

% dilate and fill
se = strel('sphere',1);
V = imdilate(V,se);

% find the axes with fewer planes and fill holes in those planes
[n, idx] = min(size(V));
x = cell(n, 1); 
y = cell(n, 1); 
z = cell(n, 1); 
V_dilate = cell(n, 1); 
for m=1:n
   if idx==1
       V_dilate{m} = squeeze(V(m, :, :)); 
   elseif idx==2
       V_dilate{m} = squeeze(V(:, m, :)); 
   else
       V_dilate{m} = squeeze(V(:, :, m)); 
   end
end
parfor m=1:n
    if idx==1
        temp = imfill(V_dilate{m}, 'hole');
        [y{m}, z{m}, ~] = find(temp); 
        x{m} = m*ones(size(y{m})); 
    elseif idx==2
        temp = imfill(V_dilate{m}, 'hole');
        [x{m}, z{m}, ~] = find(temp); 
        y{m} = m*ones(size(x{m})); 
    else
        temp = imfill(V_dilate{m}, 'hole');
        [x{m}, y{m}, ~] = find(temp); 
        z{m} = m * ones(size(x{m})); 
    end
end
x = cell2mat(x); 
y = cell2mat(y); 
z = cell2mat(z); 

%% convert voxel locations to 2p locations
temp = unique(floor(([x, y, z] + subs_0 -1)/res_factor)+1, 'rows'); 
subs_2p= zeros(size(temp)); 
subs_2p(:,1) = dims_2p(2) - temp(:,2) + 1;    % the first dimension is y 
subs_2p(:,2) = temp(:,1);                   % second dimension is x 
subs_2p(:,3) = dims_2p(3) - temp(:,3) +1 ;  % the third dimension is z 


















