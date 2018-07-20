%% move to the previous scan 
% update scan id 
ease.scan_id = max(1, ease.scan_id-1);

% update the display 
set(ease.gui.edit_scan, 'string', num2str(ease.scan_id)); 