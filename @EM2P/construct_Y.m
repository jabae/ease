function [Y_raw, Y_denoised] = construct_Y(obj)
%% construct a container for storing video data. 

Y_raw = cell(obj.num_scans, obj.num_blocks); 
Y_denoised = cell(obj.num_scans, obj.num_blocks); 

for mscan=1:obj.num_scans
    for mblock=1:obj.num_blocks
        [dl_Yr, dl_Yd] = obj.create_dataloader(mscan, mblock, []);
        Y_raw{mscan, mblock} = dl_Yr; 
        Y_denoised{mscan, mblock} = dl_Yd; 
    end
end

assignin('base', 'Y_raw', Y_raw); 
assignin('base', 'Y_denoised', Y_denoised); 