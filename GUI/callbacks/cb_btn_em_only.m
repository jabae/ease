%% merge
tmp_img = img_em_ca;
set(ease.gui.check_em_only, 'value', ease.show_em_only);

if ease.show_em_only
    tmp_img(:, :, 1, :) = 0;
    set(ease.gui.btn_em_only, 'backgroundcolor', ease.gui.color_pink);
else
    set(ease.gui.btn_em_only, 'backgroundcolor', ease.gui.color_gray);
end
for m=1:ease.num_slices
    hold(ease.gui.ax_em{m}, 'off');
    imagesc(squeeze(tmp_img(:, :, :, m)), 'parent', ease.gui.ax_em{m});
    set(ease.gui.ax_em{m}, 'xtick', [], 'ytick',[]);
    hold(ease.gui.ax_em{m}, 'on');
    
    z = ease.video_zvals_updated(ease.scan_id, mslice);
    dx = ease.em_shifts.jj(ease.scan_id, m);
    dy = ease.em_shifts.ii(ease.scan_id, m);
    tmp_em = ease.em_ranges{z};
    
    if ~isempty(tmp_em)
        xi = (tmp_em(:,2)-ease.FOV_stack(3)-dx)/ssub;
        yi = (tmp_em(:,1)-ease.FOV_stack(1)-dy)/ssub;
        plot(xi, yi, '-.r', 'linewidth', 2, 'parent', ease.gui.ax_em{m});
    end
end

