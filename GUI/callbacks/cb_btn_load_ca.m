%% load video data into the memory 
ease.load_video_mem(); 

%% calculate summary statistics for the selected data 
summary_images = ease.calculate_summary_images(); 

cb_btn_slice;

%% neuron 
neuron = ease.get_MF3D(ease.video_T); 
cell_id = 1; 
