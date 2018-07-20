%{
nda.Soma (manual) # somatic location
-> nda.EMCell
-> nda.Stack
---
soma_x                      : float                         # x position (um)
soma_y                      : float                         # y position (um)
soma_z                      : float                         # z position (um)
%}

classdef Soma < dj.Relvar
end
