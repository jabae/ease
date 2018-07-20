%% descrease the blurring size 
ease.em_zblur = str2double(get(ease.gui.edit_blur, 'string')); 

ease.em_zblur = max(0, round(ease.em_zblur)-1);

set(ease.gui.edit_blur, 'string', num2str(ease.em_zblur));
