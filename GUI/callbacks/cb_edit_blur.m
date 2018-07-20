%% set the value of the blurring size 
temp = str2double(get(ease.gui.edit_blur, 'string'));
ease.em_zblur = max(0, round(temp)); 

set(ease.gui.edit_blur, 'string', ease.em_zblur); 