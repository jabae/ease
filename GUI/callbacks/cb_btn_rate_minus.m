tmp_confidence = neuron.match_status.confidence(cell_id); 

tmp_confidence = max(round(tmp_confidence-1), 0); 

neuron.match_status.confidence(cell_id) = tmp_confidence; 

set(ease.gui.edit_confidence, 'string', int2str(tmp_confidence)); 