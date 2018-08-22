%% video shifts and FOV 
shifts_video_ii = ease.video_shifts.ii; 
shifts_video_jj = ease.video_shifts.jj; 
FOV_ = ease.FOV; 

%% loading data
for mscan = 1:ease.num_scans
    for mblock = 1:ease.num_blocks
        tmp_scan = cell(1, ease.num_slices); 
        for mslice = 1:ease.num_slices
            % data loader
            dl_video = ease.video_loader{mscan, mslice, mblock};

            % shifts in x & y direction 
            dy = round(shifts_video_ii(mscan, mslice));
            dx = round(shifts_video_jj(mscan, mslice));
            T = dl_video.num_frames;
            tmp_scan{mslice} = dl_video.load_tzrc([], [1,1], FOV_(1:2)+dy, FOV_(3:4)+dx);
        end 
        [d1, d2, T] = size(tmp_scan{1}); 
        data = zeros(d1, d2, ease.num_slices, T, 'like', tmp_scan{1}(1)); 
        for mslice=1:ease.num_slices
           data(:, :, mslice, :) = tmp_scan{mslice};  
        end
        
        save(fullfile(ease.raw_folder, sprintf('scan%d_block%d_complete.mat', mscan, mblock)), ...
            'data', '-v7.3'); 
        disp([mscan, mblock]); 
    end
end

