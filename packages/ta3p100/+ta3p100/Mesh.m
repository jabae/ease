%{
# Pinky100 Mesh structures
-> ta3p100.Segment
-> ta3p100.Segment
---
n_vertices                  : bigint                        # number of vertices in this mesh
n_triangles                 : bigint                        # number of triangles in this mesh
vertices                    : longblob                      # x,y,z coordinates of vertices
triangles                   : longblob                      # triangles (triplets of vertices)
%}


classdef Mesh < dj.Imported

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			 self.insert(key)
		end
	end

end