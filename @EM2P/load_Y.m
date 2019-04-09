function Y = load_Y(obj, mscan, mblock)
%% load the video data for running source extraction
if ~exist('mscan', 'var') || isempty(mscan)
    mscan = obj.scan_id;
end
if ~exist('mblock', 'var')||isempty(mblock)
    mblock = obj.block_id;
end

% use the denoised data or not
if obj.use_denoise
    var_Y = 'Y_denoised';
else
    var_Y = 'Y_raw';
end

% get data dimension
tmp_idl = obj.create_dataloader(1,1);
dims = tmp_idl.dims;
dvoxel = prod(dims);

% check if (Y_raw, Y_denoised) exists
if ~exist_in_workspace('Y_raw', 'base')
    obj.construct_Y();
end
%% load data
if mblock>0
    fprintf('loading block %d of scan %d\n', ...
        mblock, mscan);
    % load one block only
    tmp_str = sprintf('%s{%d, %d}', var_Y, mscan, mblock);
    tmpY = evalin('base', tmp_str);
    
    if ~isnumeric(tmpY)
        % replace the dataloader with its actual data
        evalin('base', sprintf('%s{%d,%d}=%s{%d,%d}.load_tzrc();', var_Y, mscan, mblock, var_Y, mscan, mblock));
        tmpY = evalin('base', tmp_str);
    end
    Y = reshape(tmpY, dvoxel, []);
    fprintf('done\n');
else
    % load all blocks and then concatenate them
    fprintf('loading all %d blocks of scan %d\n', ...
        obj.num_blocks, mscan);
    Y = cell(1, obj.num_blocks);
    for m=1:obj.num_blocks
        tmp_str = sprintf('%s{%d, %d}', var_Y, mscan, m);
        tmpY = evalin('base', tmp_str);
        if ~isnumeric(tmpY)
            evalin('base', sprintf('%s{%d,%d}=%s{%d,%d}.load_tzrc();', var_Y, mscan, m, var_Y, mscan, m));
            tmpY = evalin('base', tmp_str);
        end
        Y{1, m} = reshape(tmpY, dvoxel, []);
        fprintf('block %d: done\n', m);
    end
    Y = cell2mat(Y);
    fprintf('done\n');
    
end

