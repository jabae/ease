% to confirm whether you want to clear the old results 
if ease.options_init.clear_results
   temp = input('do you really want to clear existing results?');  
end
if ~strcmpi(temp, 'y')
    ease.options_init.clear_results = false; 
end

% do initialization 
fprintf('EASE is trying to add %d neurons from the top %d EM components\n...\n', ...
    ease.options_init.K_new, ease.options_init.K_candidate); 
neuron.ease_initialization(Y, Aem, segment_ids, ease.options_init, ...
    neuron.black_list, neuron.white_list);
fprintf('Done\n\n'); 

% update model variables and evaluate the matching scores 
neuron.hals(Y); 
neuron.evaluate_matching_confidence(Y, Aem, segment_ids); 