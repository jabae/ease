if ease.block_id==0
    % load all blocks if block_id = 0;
    fprintf('You are going to load all %d blocks of scan %d\n', ...
        ease.num_blocks, ease.scan_id);
    neuron = ease.get_MF3D();

    T_all = 0; % total number of frames
    Y_all = cell(1, ease.num_slices);
    for mmblock=1:ease.num_blocks
        fprintf('loading block %d of scan %d...', mmblock, ease.scan_id);
        if ease.use_denoise
            tmpY =  Y_denoised{ease.scan_id, mmblock};
        else
            tmpY = Y_raw{ease.scan_id, mmblock};
        end
        if isempty(tmpY)
            tmp_neuron = ease.load_video_mem(ease.scan_id, mmblock);  % load data of one scan
            T_all = T_all + size(Y_in_use, 2);
            Y_all{mmblock} = tmp_neuron.reshape(Y_in_use, 1);
        else
            Y_all{mmblock} = neuron.reshape(tmpY, 1);
        end
        fprintf('Done\n');
    end
    summary_images = ease.calculate_summary_images();
    
    Y_in_use = cell2mat(Y_all);
    ease.block_id = 0;
    clear Y_all tmpY;
else
    %% load video data into the memory
    neuron = ease.load_video_mem();
    
    %% calculate summary statistics for the selected data
    summary_images = ease.calculate_summary_images();
end
