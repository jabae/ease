function options = set_projections_options(obj)

%% set the options for projecting EM masks onto video planes 
options.dims_2p = obj.dims_stack;
options.zblur = obj.em_zblur;
options.zvals = reshape(obj.video_zvals_updated', 1, []);
options.norm_z = 2*options.zblur^2;
options.min_voxels = 10; 
end