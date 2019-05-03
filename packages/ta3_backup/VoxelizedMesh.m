%{
# my newest table
-> ta3.Mesh
-----
n_fragments          : int                    # number of fragments
n_voxels             : int                    # number of voxels
n_vertices           : int                    # number of vertices
n_faces              : int                    # number of faces
indices              : longblob               # indices of nonzero voxels
%}

classdef VoxelizedMesh < dj.Computed
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            options = evalin('base', 'options');
            % computed the indices of nonzero voxels 
            [vertices, faces]= fetchn((ta3.MeshFragment & sprintf(...
                'segmentation=%d and segment_id=%d', key.segmentation, key.segment_id)),...
                'vertices', 'triangles');
            [subs_2p, ~] = mesh2volume(vertices, faces, options);
            key.indices = sub2ind(options.dims_2p, subs_2p(:,1), subs_2p(:,2), subs_2p(:,3));
            
            key.n_voxels = length(key.indices); 
            key.n_vertices = sum(cellfun(@length, vertices)); 
            key.n_faces = sum(cellfun(@length, faces)); 
            key.n_fragments = length(vertices); 
            self.insert(key)
        end
    end
    
end