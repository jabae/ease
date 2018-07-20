%% use the specified slice 

tmp_id = find(z_ranges==str2double(get(text_align_2p, 'string')));
if isempty(tmp_id)
    tmp_id = 1;
else
    tmp_id = max(1, tmp_id);
    id_2p = min(tmp_id, length(z_ranges));
end
set(text_align_2p, 'string', num2str(z_range(id_2p)));

imagesc(imgs_2p(:, :, id_2p), 'parent', ax_align_2p);
set(ax_align_2p, 'xtick', [], 'ytick',[]);

temp = imgs_2p(:, :, id_2p);
max_2p = quantile(temp(:), 0.99);
img_2p_em(:, :, 1) = imgs_2p(:, :, id_2p)/max_2p;
imagesc(img_2p_em, 'parent', ax_align_merge);
set(ax_align_merge, 'xtick', [], 'ytick',[]);

