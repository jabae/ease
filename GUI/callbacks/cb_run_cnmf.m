%% run HALS to update model variables. 
neuron.merge_repeats(); 
fprintf('EASE is running CNMF to update model variables\n...\n')
neuron.hals(Y); 
fprintf('Done\n\n'); 

