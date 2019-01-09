if ease.options_init.clear_results
   temp = input('do you really want to clear existing results?');  
end
if ~strcmpi(temp, 'y')
    ease.options_init.clear_results = false; 
end
fprintf('EASE is trying to add %d neurons from the top %d EM components\n...\n', ...
    ease.options_init.K_new, ease.options_init.K_candidate); 
neuron.initialize_em(Aem, [], ease.options_init);
fprintf('Done\n\n'); 

neuron.evaluate_matching_confidence(Aem); 