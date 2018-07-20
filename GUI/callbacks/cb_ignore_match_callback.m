neuron.match_status.status(cell_id) = 0;
temp = neuron.match_status.em_ids{cell_id};

if ~any(temp==-em_id)
    neuron.match_status.em_ids{cell_id} = [temp, -em_id];
    set(ease.gui.btn_em_ignore, 'backgroundcolor', 'red');
    set(ease.gui.btn_em_candidate, 'backgroundcolor', ease.gui.color_gray);
end
if any(temp==em_id)
    temp(temp==em_id) = [];
    neuron.match_status.em_ids{cell_id}(temp==em_id) = temp;
    set(ease.gui.btn_em_candidate, 'backgroundcolor', ease.gui.color_gray);
end