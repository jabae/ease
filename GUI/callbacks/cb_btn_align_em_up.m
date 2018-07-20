%% move the em slice up 
id_em = min(length(z_range), id_em+1);
set(text_align_em, 'string', num2str(z_range(id_em))); 

imagesc(imgs_em(:, :, id_em), 'parent', ax_align_em);
set(ax_align_em, 'xtick', [], 'ytick',[]); 

img_2p_em(:, :, 2) = imgs_em(:, :, id_em)/2;
imagesc(img_2p_em, 'parent', ax_align_merge);
set(ax_align_merge, 'xtick', [], 'ytick',[]); 
