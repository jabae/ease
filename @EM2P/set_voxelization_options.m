function options = set_voxelization_options(obj, voxel_em)
% voxel size of EM segment
if ~exist('voxel_em', 'var') || isempty(voxel_em)
    voxel_em = [200, 200, 200]/1000;   % um
elseif length(voxel_em)==1
    voxel_em = ones(1,3) * voxel_em; 
end
options.voxel_em = voxel_em; % voxelize EM data with this resolution
options.scale_factor = obj.em_scale_factor;

% spatial transformation
[A_convert, offset] = obj.get_transformation();
options.A = A_convert;
options.offset = offset;

% voxel resolution of 2p stack data.
options.dims_2p = obj.dims_stack;
options.range_2p = obj.range_2p;

% parallel processing
options.use_parallel = true;
end