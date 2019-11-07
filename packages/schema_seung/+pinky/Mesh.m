%{
# Segment meshes
-> pinky.Segment
---
n_vertices : int
n_triangles : int
vertices : longblob
triangles : longblob
%}

classdef Mesh < dj.Manual
end
