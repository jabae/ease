%%
if ~exist('cell_id', 'var')
    cell_id = 1;
end
K_ca = size(neuron.A, 2);

cell_id = min(cell_id, K_ca);

set(ease.gui.text_cell_id, 'string', num2str(cell_id));

%% show temporal activity
cla(ease.gui.ax_activity);
hold(ease.gui.ax_activity, 'on');
if ~isempty(neuron.C_raw)
    plot(neuron.C_raw(cell_id, :), 'parent', ease.gui.ax_activity);
    plot(neuron.C(cell_id, :), 'parent', ease.gui.ax_activity);
else
    plot(neuron.C(cell_id, :), 'parent', ease.gui.ax_activity);
end
axis(ease.gui.ax_activity, 'tight');

%% show spatial components
ssub = ease.ssub; 
ai = neuron.reshape(neuron.A(:, cell_id), 3);
if ~isempty(neuron.A_corr)
    ai_corr = neuron.reshape(neuron.A_corr(:, cell_id), 3);
else
    ai_corr = ai;
end
ai_max = max(ai(:))+0.1;
ai_corr_max = max(ai_corr(:))+0.1;
if isnan(ai_corr_max)
    ai_corr_max = 1; 
end
for m=1:ease.num_slices
    z = ease.video_zvals_updated(ease.scan_id, m);
    dx = ease.em_shifts.jj(ease.scan_id, m);
    dy = ease.em_shifts.ii(ease.scan_id, m);
    
    axes(ease.gui.ax_corr{m}); cla;
    imagesc(ai_corr(:, :, m), [0, ai_corr_max]);
    axis off tight;
    hold on;
    
    tmp_em = ease.em_boundary{ease.scan_id};
    
    if ~isempty(tmp_em)
        xi = (tmp_em(:,2)-ease.FOV_stack(3)-dx)/ssub;
        yi = (tmp_em(:,1)-ease.FOV_stack(1)-dy)/ssub;
        plot(xi, yi, '-.r', 'linewidth', 2);
    end
    %% estimate ai and compute the pairwise correlation
    axes(ease.gui.ax_em{m}); cla;
    imagesc(ai(:, :, m), [0, ai_max]);
    axis off tight;
end

%% matching status
if iscell(segment_ids)
    segment_ids = cell2mat(segment_ids); 
end
if iscell(Aem)
    Aem = cell2mat(Aem); 
end

tmp_status = neuron.match_status;
%show confidence score 
cb_btn_rate; 
cb_btn_label; 
if tmp_status.status(cell_id) == -1
    col_status = {ease.gui.color_gray, ease.gui.color_gray, ease.gui.color_gray, 'red'};
elseif tmp_status.status(cell_id) == 0
    col_status = {ease.gui.color_gray, [0.5, 0, 0], ease.gui.color_gray, ease.gui.color_gray};
else
    col_status = {[1, 0, 0], ease.gui.color_gray, ease.gui.color_gray, ease.gui.color_gray};
end
% 
% set(ease.gui.btn_em_rematch, 'backgroundcolor', col_status{1});
% set(ease.gui.btn_em_candidate, 'backgroundcolor', col_status{2});
% set(ease.gui.btn_em_zero, 'backgroundcolor', col_status{4});
% set(ease.gui.btn_em_ignore, 'backgroundcolor', col_status{3});

%% run matching neurons
cb_show_em_neuron;

