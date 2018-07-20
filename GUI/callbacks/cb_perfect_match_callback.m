neuron.match_status.status(cell_id) = 1; 
neuron.match_status.em_ids{cell_id} = em_id; 
set(gco, 'backgroundcolor', 'red'); 
set(ease.gui.btn_em_candidate, 'backgroundcolor', ease.gui.color_gray); 
set(ease.gui.btn_em_ignore, 'backgroundcolor', ease.gui.color_gray);
set(ease.gui.btn_em_zero, 'backgroundcolor', ease.gui.color_gray); 