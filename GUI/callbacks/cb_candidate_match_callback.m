neuron.match_status.status(cell_id) = 0;
temp = neuron.match_status.em_ids{cell_id};

if ~any(temp==em_id);
    neuron.match_status.em_ids{cell_id} = [temp, em_id];
    set(ease.gui.btn_em_candidate, 'backgroundcolor', 'red'); 
end
