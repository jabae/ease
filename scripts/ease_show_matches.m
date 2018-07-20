%%
status = neuron.match_status.status;
ind = find(status==1);
match_ids_ca = ind;
match_ids_em = cell2mat(neuron.match_status.em_ids(ind));
ind = false(size(match_ids_ca));
A_corr = neuron.A_corr;
A_corr(~ind_em_volume, :) = 0;
% for m=1:length(ind)
%     if max(A_corr(:, m))<0.4
%         ind(m) = true;
%     end
% end
match_ids_ca(ind) = [];
match_ids_em(ind) = [];
d1 = ease.d1;
d2 = ease.d2;
ssub = ease.dims_stack(1) / ease.dims_video(1);
d1_2p = d1 * ssub;
d2_2p = d2 * ssub;
%%

for m_match = 1:length(match_ids_ca)
figure('papersize', [3*d2, 3*d1]/max(d1, d2)*5);
init_fig;
[ha, pos] = tight_subplot(4, 3);
    ha(10).Position(1) = ha(1).Position(1);
    ha(10).Position(3) = ha(3).Position(1) + ha(3).Position(3);
    
    cell_id = match_ids_ca(m_match);
    em_id = match_ids_em(m_match);
    ai_corr = neuron.reshape(neuron.A_corr(:, cell_id), 3);
    ai = neuron.reshape(neuron.A(:, cell_id), 3);
    ai_max = max(ai(:));
    ci_raw = neuron.C_raw(cell_id, :);
    ci = neuron.C(cell_id, :);
    % generate the EM masks
    ai_em = zeros(d1, d2, ease.num_slices);
    img_em_ca = zeros(d1, d2, 3, ease.num_slices);
    ai_em = neuron.reshape(Aem_ds(:, em_id), 3); 
%     for m=1:ease.num_slices
%         temp = Aem{m}(:, em_id);
%         temp = reshape(temp, d1_2p, d2_2p);
%         ai_em(:, :, m) = imresize(full(temp), [d1, d2], 'box');
%     end
    ai_em_max = max(ai_em(:));
    
    for m=1:ease.num_slices
        tmp_em = ease.em_ranges{ease.video_zvals_updated(ease.scan_id, m)};
        % get contours
        A_temp = ai(:, :, m);
        [tmp1, tmp2, ~] = find(A_temp);
        
        img_em_ca(:, :, 2, m) = ai(:, :, m) /ai_max;
        img_em_ca(:, :, 1, m) = (ai_em(:, :, m)>0)/ai_em_max;
        
        axes(ha(m)); cla;
        imagesc(ai_corr(:, :, m), [0, 1]);
        axis equal tight; hold on;
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        if ~isempty(tmp_em)
            xi = (tmp_em(:,2)-ease.FOV_stack(3))/ssub;
            yi = (tmp_em(:,1)-ease.FOV_stack(1))/ssub;
            plot(xi, yi, '-.r', 'linewidth', 2);
        end
        title(sprintf('slice %d', m));
        if m==1
            ylabel('corr(Y, c_i)');
        end
        axes(ha(m+6));
        imagesc(ai_em(:, :,  m), [0, ai_em_max]);
        axis equal  tight;
        
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        if ~isempty(tmp_em)
            hold on;
            xi = (tmp_em(:,2)-ease.FOV_stack(3))/ssub;
            yi = (tmp_em(:,1)-ease.FOV_stack(1))/ssub;
            plot(xi, yi, '-.r', 'linewidth', 2);
        end
        if m==1
            ylabel('EM');
        end
        axes(ha(m+3));
        imagesc(ai(:, :,  m), [0, ai_max]);
        axis equal  tight;
        
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        if ~isempty(tmp_em)
            hold on;
            xi = (tmp_em(:,2)-ease.FOV_stack(3))/ssub;
            yi = (tmp_em(:,1)-ease.FOV_stack(1))/ssub;
            plot(xi, yi, '-.r', 'linewidth', 2);
        end
        
        if m==1
            ylabel('a_i');
        end
    end
    
    
    axes(ha(10));
    cla;
    plot(ci_raw);
    hold on;
    plot(ci, 'r'); axis tight;
    legend('raw', 'denoised');
    saveas(gcf, fullfile(EASE_dir, 'matched_denoised', sprintf('match_%d.png', m_match)));
close; 
end