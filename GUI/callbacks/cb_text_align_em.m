%% use the specified slice 
tmp_id = find(z_ranges==str2double(get(text_align_em, 'string')));
if isempty(tmp_id)
    tmp_id = 1;
else
    tmp_id = max(1, tmp_id);
    id_em = min(tmp_id, length(z_ranges));
end
set(text_align_em, 'string', num2str(z_range(id_em)));

%% move the em slice down 
imagesc(imgs_em(:, :, id_em), 'parent', ax_align_em);
set(ax_align_em, 'xtick', [], 'ytick',[]); 

img_2p_em(:, :, 2) = imgs_em(:, :, id_em)/2;
imagesc(img_2p_em, 'parent', ax_align_merge);
set(ax_align_merge, 'xtick', [], 'ytick',[]); 