%% check how well NMf component matches the EM component
set(ease.gui.text_em_current, 'string', num2str(em_rank));
em_id = em_sort_id(em_rank);

%% show matching scores
cla(ease.gui.ax_score);
plot(em_sorted_scores, '-*', 'parent', ease.gui.ax_score);
hold(ease.gui.ax_score, 'on');
plot(em_rank, em_scores(em_id), 'or', 'parent', ease.gui.ax_score);
tmp_ids = tmp_status.em_ids{cell_id};
if tmp_status.status(cell_id)==0
    for m=1:length(tmp_ids)
        tmp_rank = find(em_sort_id==abs(tmp_ids(m)));
        if em_sort_id(tmp_rank)== tmp_ids(m)
            % candidate match
            plot(tmp_rank, em_sorted_scores(tmp_rank), 'sm', 'parent', ease.gui.ax_score);
        else
            % ignored
            plot(tmp_rank, em_sorted_scores(tmp_rank), 'sg', 'parent', ease.gui.ax_score);
        end
    end
elseif tmp_status.status(cell_id)==1
    plot(em_sorted_scores, '-*', 'color', ease.gui.color_gray*0.8, 'parent', ease.gui.ax_score);
    tmp_rank = find(em_sort_id==tmp_ids);
    plot(tmp_rank, em_sorted_scores(tmp_rank), 'sm', 'parent', ease.gui.ax_score);
end
set(ease.gui.ax_score, 'ylim',[0, max(em_scores*1.1)], 'xlim', [0, em_rank+20]);

%% generate the EM masks
% ai_em = zeros(ease.d1, ease.d2, ease.num_slices);
% ssub = ease.dims_stack(1) / ease.dims_video(1); 
% for m=1:ease.num_slices
%     temp = Aem{m}(:, em_id);
%     temp = reshape(temp, ease.d1*ssub, ease.d2*ssub);
%     ai_em(:, :, m) = imresize(full(temp), [ease.d1, ease.d2], 'box');
% end
ai_em = neuron.reshape(Aem(:, em_id), 3); 

%% merge
img_em_ca = zeros(ease.d1, ease.d2, ease.num_slices);
max_ai = max(ai(:));
max_ai_em = max(ai_em(:)); 
for m=1:ease.num_slices
    img_em_ca(:, :, 1, m) = ai(:, :, m)./ max_ai * 4;
    img_em_ca(:, :, 2, m) = ai_em(:, :, m) / max_ai_em *2;
end

cb_btn_em_only;

%%
switch neuron.match_status.status(cell_id)
    case 0
        if any(neuron.match_status.em_ids{cell_id}==-em_id)
            % ignore
            set(ease.gui.btn_em_candidate, 'backgroundcolor', ease.gui.color_gray);
            set(ease.gui.btn_em_ignore, 'backgroundcolor', 'red');
        elseif any(neuron.match_status.em_ids{cell_id}==em_id)
            % candidate
            set(ease.gui.btn_em_candidate, 'backgroundcolor', 'red');
            set(ease.gui.btn_em_ignore, 'backgroundcolor', ease.gui.color_gray);
        else
            % not sure
            set(ease.gui.btn_em_candidate, 'backgroundcolor', [0.5, 0., 0]);
            set(ease.gui.btn_em_ignore, 'backgroundcolor', ease.gui.color_gray);
        end
    case 1
        if em_id == neuron.match_status.em_ids{cell_id}
            set(ease.gui.btn_em_perfect, 'backgroundcolor', 'red');
        else
            set(ease.gui.btn_em_perfect, 'backgroundcolor', ease.gui.color_gray);
        end
    otherwise
        fprintf('this neuron has is has no match\n');
end