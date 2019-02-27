%% load video data into the memory
erun.load_calcium_data; 
erun.load_em_projections; 

neuron.spatial_range = sparse(pixels_em); 
Y_in_use = neuron.reshape(Y_in_use, 1); 

% neuron.frame_range = [201, size(Y_in_use, 2)]; 
% 
% % normalize data
% [Y_cnmf, Y_sn] = neuron.normalize_data(Y_in_use); 
