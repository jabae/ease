function add2whitelist(obj, segment_ids, reason)
%% add a list of EM segments to the white list. 
%

%% avoid duplicate items 

% these neurons will be ignored by EASE 
if ~exist('reason', 'var')
    reason = ''; 
end
key.segmentation = obj.em_segmentation; 
key.comment = reason; 
rel = eval(sprintf('%s.WhiteList', obj.dj_name)); 

for m=1:length(segment_ids)
    key.segment_id = segment_ids(m); 
    rel.inserti(key); 
end 

fprintf('added %d segments to the white list\n', length(segment_ids)); 