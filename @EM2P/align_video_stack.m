function align_video_stack(obj)
%% what does this function do
%{
    align the video data and the 2p stack data in the cropped FOV
%}

%% inputs:
%{
    2p stack data and video data within the cropped FOV
%}

%% outputs:
%{
    after FOV alginment, obj.video_zvals_updated will be updated and the
    shfits of 2p data in x-y plane are saved into obj.stack_shifts
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code
if ~isempty(obj.video_zvals_updated)
    fprintf('The data have been aligned once. \n Please set ease.video_zvals_updated = [] if you want to align again.\n');
    return;
else
    fprintf('Doing a fine registration between the stack data and the video data.\n');
end

%% shifts for the stack data and video data
shifts_stack_ii = zeros(obj.num_scans, obj.num_slices);
shifts_stack_jj = shifts_stack_ii;

shifts_video_ii = obj.video_shifts.ii;
shifts_video_jj = obj.video_shifts.jj;

% z position
z_vals = obj.video_zvals;
z_shift = obj.align_max_zshift;    % the best match will be searched within frames 4 planes above and below.

%% FOV info
ssub = obj.dims_stack(1) / obj.dims_video(1);

FOV_ = obj.FOV;
FOV_stack_ = obj.FOV_stack;
d1_ = obj.d1;
d2_ = obj.d2;
d1_stack = d1_ * ssub;
d2_stack = d2_ * ssub;

%% compute the relative shift for 2p stack data
if ~exist_in_workspace('stack_2p', 'base')
    % if stack data was not loaded, load it.
    obj.load_stack();
end
stack_2p = evalin('base', 'stack_2p');

% images showing the performances of motion correction
img_match_fov = cell(obj.num_scans, obj.num_slices);
temp = (rand(20)>0.5);
tmp_kernel = temp / sum(temp(:));   % used to normalize pixel values for balancing colors

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
        dy = round(shifts_video_ii(m, n));
        dx = round(shifts_video_jj(m, n));
        T = dl_video.num_frames;
        tmpY = dl_video.load_tzrc(round([T/4, T*3/4]), [1,1], FOV_(1:2)+dy, FOV_(3:4)+dx);
        img_video = imresize(mean(tmpY, 3), ssub);
        img_video_norm = img_video ./ imfilter(img_video, tmp_kernel, 'replicate');
        % loading 2p stack data and find the best match
        z = max(1, z_vals(m, n)-z_shift);   % bottom candinate plane
        vmax = 0;
        while z<= min(obj.dims_stack(3), z_vals(m,n)+z_shift)
            img_stack = stack_2p(FOV_stack_(1):FOV_stack_(2), FOV_stack_(3):FOV_stack_(4), z);
            img_stack_norm = img_stack ./ imfilter(img_stack, tmp_kernel, 'replicate');
            temp = normxcorr2(img_video_norm, img_stack_norm);
            tmp_max = max(temp(:));
            [ii, jj, ~] = find(temp==tmp_max);
            if tmp_max>vmax
                vmax = tmp_max;
                shifts_stack_ii(m, n) = ii - d1_stack;
                shifts_stack_jj(m, n) = jj - d2_stack;
                z_vals(m, n) = z;
            end
            z = z+1;
        end
        
        % overlap two images and check the performances of the MC
        img = zeros([d1_stack, d2_stack, 3]);
        img(:, :, 1) = img_video / quantile(img_video(:), 0.99);
        tmp_r = (FOV_stack_(1):FOV_stack_(2))+shifts_stack_ii(m, n);
        tmp_c = ( FOV_stack_(3):FOV_stack_(4))+shifts_stack_jj(m, n);
        ind_r = find(tmp_r>=1 & tmp_r<=obj.dims_stack(1));
        ind_c = find(tmp_c>=1 & tmp_c<=obj.dims_stack(2));
        img(ind_r, ind_c, 2) = stack_2p(tmp_r(ind_r), tmp_c(ind_c), ...
            z_vals(m, n)) / quantile(img_stack(:), 0.99);
        
        img_match_fov{m, n} = img;
        fprintf('.');
    end
end
fprintf('\n');

obj.FOV = FOV_;
obj.FOV_stack = FOV_stack_;
obj.stack_shifts = struct('ii', shifts_stack_ii, 'jj', shifts_stack_jj);
obj.em_shifts = obj.stack_shifts;
obj.video_zvals_updated = z_vals;
obj.aligned_images = img_match_fov;
end