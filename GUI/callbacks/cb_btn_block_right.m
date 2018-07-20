%% move the next block 
% update block_id 
ease.block_id = min(ease.num_blocks, ease.block_id+1);

% update the display 
set(ease.gui.edit_block, 'string', num2str(ease.block_id)); 