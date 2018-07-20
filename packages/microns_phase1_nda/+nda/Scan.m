%{
nda.Scan (manual) # Two-photon functional scans
scan_idx        : smallint 		# scan ID
-----
depth           : int           # Scan depth from the surface (microns)
laser_power 	: int           # Laser power (mW)
wavelength      : int           # Laser wavelength (nm)
filename 		: varchar(255) 	# Scan base filename uploaded to S3
%}

classdef Scan < dj.Relvar
end