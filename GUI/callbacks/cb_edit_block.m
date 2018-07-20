%% set the block_id to the typed value 

temp = round(str2double(get(ease.gui.edit_block, 'string')));

if temp>=1 && temp<=ease.num_blocks
    ease.block_id = temp; 
else
    set(ease.gui.edit_block, 'string', num2str(ease.block_id)); 
end