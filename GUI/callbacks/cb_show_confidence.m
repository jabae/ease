tmp_confidence = neuron.match_status.confidence;
if isempty(tmp_confidence) || length(tmp_confidence)~=size(neuron.A_corr, 2)
    fprintf('compute matching confidence...\n'); 
    neuron.evaluate_matching_confidence(Y, Aem, segment_ids);
    fprintf('done'); 
end

neuron.orderROIs('confidence');
cb_show_2p_neuron; 