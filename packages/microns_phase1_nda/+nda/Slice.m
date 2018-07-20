%{
nda.Slice (manual) # Z-axis slice in scan
-> nda.Scan
slice           : smallint 		# slice number
-----
z_offset 		: float 		# Z-offset from scan depth (microns)
%}

classdef Slice < dj.Relvar
end