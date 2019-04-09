%% show summary statistics of the video data 
% EM info 
if ~exist('EM_masks', 'var')
    ease.get_em_masks(); 
end

% get summary images 
if ~exist('summary_images', 'var') 
    summary_images = ease.summary_images; 
end 

ssub = ease.dims_stack(1) / ease.dims_video(1); 
mscan = ease.scan_id; 

for mslice = 1:ease.num_slices
    z = ease.video_zvals_updated(ease.scan_id, mslice);
    hold(ease.gui.ax_slice{mslice}, 'on')
    imagesc(eval(sprintf('summary_images.%s(:, :, mslice)', ease.nam_show)), 'parent', ease.gui.ax_slice{mslice}); 
    axis(ease.gui.ax_slice{mslice}, 'tight');
    set(ease.gui.ax_slice{mslice}, 'xtick', [], 'ytick', [], 'fontsize', 14); 
    ylabel(ease.gui.ax_slice{mslice}, sprintf('z=%d', z)); 
    
    tmp_em = ease.em_boundary{ease.scan_id}; 
    dx = ease.em_shifts.jj(mscan, mslice); 
    dy = ease.em_shifts.ii(mscan, mslice); 
    if ~isempty(tmp_em)
        xi = (tmp_em(:,2)-ease.FOV_stack(3)-dx)/ssub;
        yi = (tmp_em(:,1)-ease.FOV_stack(1)-dy)/ssub;
        plot(xi, yi, '-.r', 'linewidth', 2, 'parent', ease.gui.ax_slice{mslice}); 
    end
    set(ease.gui.ax_slice{mslice}, 'Ydir', 'reverse'); 
end