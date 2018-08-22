%% remove bad components
neuron.scores = neuron.calculate_matching_scores(Aem_ds, 'sim');
ind_rank = neuron.get_match_rank();
neuron.delete(ind_rank>5);
neuron.hals(Y_cnmf);