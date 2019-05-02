% ease_load_calcium_data; 
neuron = ease.get_MF3D(); 

Y = ease.load_Y(); 
summary_images = ease.summary_images; 
if ~isempty(ease.gui) && isvalid(ease.gui.fig_main)
    cb_btn_slice;
end 
% load EM 
[Aem, segment_ids] = ease.load_Aem(); 
