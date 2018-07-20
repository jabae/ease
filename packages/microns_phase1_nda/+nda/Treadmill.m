%{
# Treadmill activity synced to imaging
-> nda.Scan
---
treadmill_speed             : longblob                      # vector of treadmill velocities synchronized with slice 1 frame times (cm/s)
treadmill_vel               : longblob                      # velocity of treadmill (may be positive or negative)
%}

classdef Treadmill < dj.Manual
end
