ind_bad = find(neuron.match_status.confidence==0); 
neuron.delete(ind_bad);