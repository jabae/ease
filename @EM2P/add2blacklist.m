function add2blacklist(obj, segment_ids, reason)
%% add a list of EM segments to the black list. 
% these neurons will be ignored by EASE 

key.segmentation = obj.em_segmentation; 
key.comment = reason; 
rel = eval(sprintf('%s.BlackList', obj.dj_name)); 

for m=1:length(segment_ids)
    key.segment_id = segment_ids(m); 
    rel.inserti(key); 
end 

fprintf('added %d segments to the black list\n', length(segment_ids)); 