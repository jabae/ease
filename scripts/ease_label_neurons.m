%% start label neurons
K = size(neuron.A, 2);
if ~exist('cell_id', 'var')
    cell_Id = 1;
else
    cell_id = max(1, cell_id-1);
end
while true
    clc;
    %% guide
    fprintf('******************* Guide *******************\n');
    fprintf('input for cell ID: \n');
    fprintf('\t0: stop this iteration\n');
    fprintf('\t-k: back to the previous k-th neuron\n');
    fprintf('score the confidence\n');
    fprintf('\tscores span from 0 to 5 and 5 is the best.\n');
    fprintf('\tneurons with 0 scores will be deleted\n');
    fprintf('label somas and dendrites\n');
    fprintf('\ts: somas\n\td: dendrites\n');
    fprintf('*********************************************\n');
    
    % show neuron
    cell_id = min(K, cell_id+1);
    ease_show_2p_neuron;
    
    % choose cell id
    commandwindow;
    temp = input('Cell ID: ');
    if temp==0
        break;
    elseif temp<0
        cell_id = max(1, cell_id+temp);
        ease_show_2p_neuron;
    end
    
    
    % mark neuron confidences
    commandwindow;
    tmp_score = input('confidence: ');
    if ~isempty(tmp_score)
        neuron.match_status.confidence(cell_id)= min(max(0, round(tmp_score)), 5);
    end
    cb_btn_rate;
    
    % label neurons
    commandwindow;
    temp = input('soma (s) or dendrite (d): ', 's');
    if strcmpi(temp, 's')
        neuron.labels(cell_id) = 1; cb_btn_label;
    elseif strcmpi(temp, 'd')
        neuron.labels(cell_id) = 2; cb_btn_label;
    end
    
    if cell_id == K 
        break; 
    end
end