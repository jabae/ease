function freespace_Y(obj)
%% free the space of ca imaging video by replacing it with a data loader

flag_in_use = false(obj.num_scans, obj.num_blocks);
if obj.block_id==0
    flag_in_use(obj.scan_id, :) = true;
else
    flag_in_use(obj.scan_id, obj.block_id) = true;
end

for mscan=1:obj.num_scans
    for mblock=1:obj.num_blocks
        if flag_in_use
            continue;
        else
            tmp_str = sprintf('[Y_raw{%d, %d}, Y_denoised{%d, %d}]=ease.create_dataloader(%d, %d,[]);', ...
                mscan, mblock, mscan, mblock, mscan, mblock);
            evalin('base', tmp_str);
        end
    end
end