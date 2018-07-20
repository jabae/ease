%% set the slice_id to the typed value 

temp = round(str2double(get(ease.gui.edit_slice, 'string')));

if temp>=1 && temp<=ease.num_slices
    ease.slice_id = temp; 
else
    set(ease.gui.edit_slice, 'string', num2str(ease.slice_id)); 
end