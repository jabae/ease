%{
# 
-> ta3p100.Segment
annotation_timestamp=CURRENT_TIMESTAMP: timestamp           # 
---
-> ta3p100.Proofreader
-> ta3p100.AnnotationLookup
annotation_comment          : varchar(4000)                 # 
%}


classdef Annotation < dj.Manual
end