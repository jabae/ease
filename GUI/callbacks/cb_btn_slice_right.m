%% move the next slice 
% update slice_id 
ease.slice_id = min(ease.num_slices, ease.slice_id+1);

% update the display 
set(ease.gui.edit_slice, 'string', num2str(ease.slice_id)); 