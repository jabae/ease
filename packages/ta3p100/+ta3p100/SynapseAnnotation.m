%{
# 
-> ta3p100.Synapse
annotation_timestamp=CURRENT_TIMESTAMP: timestamp           # 
---
-> ta3p100.Proofreader
-> ta3p100.SynapseAnnotationLookup
annotation_comment          : varchar(4000)                 # 
%}


classdef SynapseAnnotation < dj.Manual
end