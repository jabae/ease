function Y = get_Y(obj, mscan, mblock)
%% load the video data given the scan id and the block id.
%{
	a high level function for fetching the calcium imaging video in the
	given scan and block. it also supports concatenating multiple scans
	when mblock = 0. 
%}

%% inputs
%{
	obj: type; description
	mscan: integer; scan id
	mblock: integer; block id 
%}

%% outputs
%{
	Y: (d1*d2*d3)*T matrix; video data
%}

%% Author
%{
	Pengcheng Zhou 
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License 
%}

%% 
%% construct data loader
if ~exist('mscan', 'var') || isempty(mscan)
    mscan = obj.scan_id;
end
if ~exist('mblock', 'var') || isempty(mblock)
    mblock = obj.block_id;
end

if mblock>0
    % load a single block
    T = obj.video_T;
    if numel(T)>1
        T = T(mscan, mblock);
    end
    % create dataloader for fetching data. 
    [dl_Yr, dl_Yd] = obj.create_dataloader(mscan, mblock, T);
    
    % use the denoised video or not 
    if obj.use_denoise
        var_name = 'Y_denoised';
    else
        var_name = 'Y_raw';
    end
    Y = evalin('base', sprintf('%s{%d, %d}', var_name, mscan, mblock));
    
    if isempty(Y)
        fprintf('loading data for (scan %d, block %d)...\n', mscan, mblock);
        if obj.use_denoise
            Y = dl_Yd.load_tzrc();
        else
            Y = dl_Yr.load_tzrc();
        end
        assignin('base', 'tmpY', Y);
        evalin('base', sprintf('%s{%d, %d}=tmpY;', var_name, mscan, mblock), Y);
        assignin('base', 'tmpY', []);
    else
        fprintf('data loaded: (scan %d, block %d)...\n', mscan, mblock);
    end
    
else
    % load all blocks
    Y = cell(1, obj.num_blocks);
    d1 = obj.d1;
    d2 = obj.d2;
    d3 = obj.num_slices;
    for m=1:obj.num_blocks
        Y{m} = reshape(obj.get_Y(mscan, m), d1*d2*d3, []);
    end
    Y = cell2mat(Y);
end
