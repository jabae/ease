%{
# Segment: a volumetric segmented object
-> ta3p100.Segmentation
-> ta3p100.Segmentation
segment_id                  : bigint                        # segment id unique within each Segmentation
%}


classdef Segment < dj.Manual
end