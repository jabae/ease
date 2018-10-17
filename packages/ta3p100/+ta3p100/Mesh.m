%{
# Mesh
-> ta3p100.Segment
---
n_vertices                  : int                           # number of vertices in this mesh
n_triangles                 : int                           # number of triangles in this mesh
vertices                    : longblob                      # x,y,z coordinates of vertices
triangles                   : longblob                      # triangles (triplets of vertices)
%}


classdef Mesh < dj.Manual
end
