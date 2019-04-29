function options = set_projections_options(obj)

%% set the options for projecting EM masks onto video planes
options.dims_2p = obj.dims_stack;
options.zblur = obj.em_zblur;
options.zvals = reshape(obj.video_zvals_updated', 1, []);
options.norm_z = 2*options.zblur^2;
options.min_voxels = 10;

% save the option to database
key.zblur = options.zblur;
key.zvals = obj.video_zvals_updated;
key.hash = sum(prod(key.zvals,2)) * key.zblur;

rel = eval(sprintf('%s.Blurs', obj.dj_name));
[versions, tmp_hash] = rel.fetchn('version', 'hash');
if isempty(tmp_hash) 
    key.version = 1; 
    options.version = key.version;
    rel.insert(key);
elseif any(tmp_hash==key.hash)
    options.version = versions(tmp_hash==key.hash);
else
    key.version = max(versions) + 1;
    options.version = key.version;
    rel.insert(key);
end
obj.blur_version = options.version; 

end