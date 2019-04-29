%% show matching scores
em_scores = neuron.scores(cell_id, :);
[em_sorted_scores, em_sort_id] = sort(em_scores, 'descend');
em_rank = 1;

tmp_status = neuron.match_status;
switch tmp_status.status(cell_id)
    case 0
        col_status = {ease.gui.color_gray, [0.5, 0, 0], ease.gui.color_gray, ease.gui.color_gray};
        
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
        col_status = {[1, 0, 0], ease.gui.color_gray, ease.gui.color_gray, ease.gui.color_gray};
        temp = find(segment_ids==tmp_status.em_ids{cell_id}, 1, 'first');
        em_rank = find(em_sort_id==temp);
        set(ease.gui.btn_em_rematch, 'backgroundcolor', 'red');
    otherwise
        col_status = {ease.gui.color_gray, ease.gui.color_gray, ease.gui.color_gray, 'red'};
        fprintf('this neuron has is has no match\n');
end


if em_sorted_scores(1)>0
    %% check how well NMF component matches the EM component
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
end













