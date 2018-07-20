%% move the previous block 
% update block_id 
ease.block_id = max(1, ease.block_id-1);

% update the display 
set(ease.gui.edit_block, 'string', num2str(ease.block_id)); 