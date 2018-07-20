%% move the next scan 
% update scan id
ease.scan_id = min(ease.num_scans, ease.scan_id+1);

% update the display 
set(ease.gui.edit_scan, 'string', num2str(ease.scan_id)); 