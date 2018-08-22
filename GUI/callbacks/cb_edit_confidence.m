tmp_confidence = round(str2double(get(gco, 'string')));

if isnumeric(tmp_confidence)
    tmp_confidence = max(0, min(tmp_confidence, 5));
else
    tmp_confidence = neuron.match_status.confidence(cell_id);
end
set(ease.gui.edit_confidence, 'string', int2str(tmp_confidence));

neuron.match_status.confidence(cell_id) = tmp_confidence;
