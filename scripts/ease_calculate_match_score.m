%% show matching scores
em_scores = neuron.scores(cell_id, :);
[em_sorted_scores, em_sort_id] = sort(em_scores, 'descend');
em_rank = 1;

tmp_status = neuron.match_status;
if tmp_status.status(cell_id) == -1
    col_status = {ease.gui.color_gray, ease.gui.color_gray, ease.gui.color_gray, 'red'};
    return;
elseif tmp_status.status(cell_id) == 0
    col_status = {ease.gui.color_gray, [0.5, 0, 0], ease.gui.color_gray, ease.gui.color_gray};
else
    col_status = {[1, 0, 0], ease.gui.color_gray, ease.gui.color_gray, ease.gui.color_gray};
    em_rank = find(em_sort_id==tmp_status.em_ids{cell_id});
end

if em_sorted_scores(1)>0
    ease_check_match;  
end













