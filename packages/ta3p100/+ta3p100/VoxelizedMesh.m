%{
# voxelized meshes
-> ta3p100.Segment
---
indices                     : longblob                      # indices of nonzero voxels
%}

classdef VoxelizedMesh < dj.Computed
	methods(Access=protected)

        function makeTuples(self, key)
            options = evalin('base', 'options'); 
            % computed the indices of nonzero voxels
            [vertices, faces]= fetch1((ta3p100.Mesh & key), 'vertices', 'triangles');
            [subs_2p, ~] = mesh2volume(vertices, faces, options);
            key.indices = sub2ind(options.dims_2p, subs_2p(:,1), subs_2p(:,2), subs_2p(:,3));
	
            self.insert(key)
		end
	end

end
