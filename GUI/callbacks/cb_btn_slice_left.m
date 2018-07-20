%% move the previous slice 
% update slice_id 
ease.slice_id = max(1, ease.slice_id-1);

% update the display 
set(ease.gui.edit_slice, 'string', num2str(ease.slice_id)); 