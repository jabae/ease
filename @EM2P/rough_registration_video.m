function rough_registration_video(obj)


%% code
if ~isempty(obj.video_zvals)
    fprintf('The videos have been registered once. \n Please set ease.video_zvals = [] if you want to align again.\n ');
    return;
else
    fprintf('Doing a rough registration betweeen the stack data and the video ')
end

% shifts for the stack data and video data
shifts_video_ii = zeros(obj.num_scans, obj.num_slices);
shifts_video_jj = shifts_video_ii;

% FOV info
ssub = obj.ssub;
dims_video = obj.dims_video;
FOV_ = [1, dims_video(1), 1, dims_video(2)];
d1_ = dims_video(1);
d2_ = dims_video(2);

%% compute the relative shift for 2p stack data

% normalize pixel values using pixels nearby
temp = (rand(15)>0.8);  % you can use temp = ones(15)/225;
tmp_kernel = temp / sum(temp(:));

if ~exist_in_workspace('stack_2p', 'base')
    % if stack data was not loaded, load it.
    obj.load_stack();
end
stack_2p = evalin('base', 'stack_2p');

% downsampling and normalizing
stack_2p = imresize(stack_2p, 1./ssub);
stack_2p = stack_2p ./imfilter(stack_2p, tmp_kernel, 'replicate');

% results
z_vals = zeros(obj.num_scans, obj.num_slices);

% show progress bar
fprintf('\n');
for m=1:obj.num_scans * obj.num_slices
    fprintf('|');
end
fprintf('\n');

for m=1:obj.num_scans
    for n=1:obj.num_slices
        dl_video = obj.video_loader{m, n, 1};
        
        % loading video data.
        T = dl_video.num_frames;
        tmpY = dl_video.load_tzrc(round([T*3/8, T*5/8]), [1,1], FOV_(1:2), FOV_(3:4));
        img_video = mean(tmpY, 3);
        img_video_norm = img_video ./ imfilter(img_video, tmp_kernel, 'replicate');
        % loading 2p stack data and find the best match
        z = 1;  % bottom candinate plane
        vmax = 0;
        while z<= obj.dims_stack(3)
            img_stack = stack_2p(:, :, z);
            temp = normxcorr2(img_stack, img_video_norm);
            tmp_max = max(temp(:));
            [ii, jj, ~] = find(temp==tmp_max);
            if tmp_max>vmax
                vmax = tmp_max;
                shifts_video_ii(m, n) = ii - d1_;
                shifts_video_jj(m, n) = jj - d2_;
                z_vals(m, n) = z;
            end
            z = z+1;
        end
        fprintf('.');
    end
end
fprintf('\n');

obj.video_shifts = struct('ii', shifts_video_ii,...
    'jj', shifts_video_jj);
obj.video_zvals = z_vals;
end