%{
nda.Mask (manual) # cell masks
-> nda.Slice
-> nda.EMCell
-----
mask_pixels : blob          # vector of pixel indices in mask

%}

classdef Mask < dj.Relvar
end