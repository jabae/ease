%% set the scan_id to the typed value 

temp = round(str2double(get(ease.gui.edit_scan, 'string')));

if temp>=1 && temp<=ease.num_scans
    ease.scan_id = temp; 
else
    set(ease.gui.edit_scan, 'string', num2str(ease.scan_id)); 
end