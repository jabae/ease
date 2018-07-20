%{
nda.Trace (manual) # calcium traces
->nda.Mask
-----
trace 		: longblob 	# calcium trace (int16)
%}

classdef Trace < dj.Relvar
end