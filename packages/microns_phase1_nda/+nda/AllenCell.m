%{
nda.AllenCell (manual) # correct Allen cell ids in electron microscopy volume
allen_id    : smallint 		# Allen EM ID
em_id 		: smallint 		# functional mask ID
-----
%}

classdef AllenCell < dj.Relvar
end