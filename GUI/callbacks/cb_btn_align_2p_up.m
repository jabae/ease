%% move the 2p slice up 
id_2p = min(length(z_range), id_2p+1);
set(text_align_2p, 'string', num2str(z_range(id_2p))); 

imagesc(imgs_2p(:, :, id_2p), 'parent', ax_align_2p);
set(ax_align_2p, 'xtick', [], 'ytick',[]); 

temp = imgs_2p(:, :, id_2p); 
max_2p = quantile(temp(:), 0.99);
img_2p_em(:, :, 1) = imgs_2p(:, :, id_2p)/max_2p;
imagesc(img_2p_em, 'parent', ax_align_merge);
set(ax_align_merge, 'xtick', [], 'ytick',[]); 
