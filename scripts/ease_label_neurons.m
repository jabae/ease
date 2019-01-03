%% start label neurons
neuron.orderROIs('pnr'); 
ease.startGUI(); 
cell_id = 1; cb_btn_slice; ease_show_2p_neuron;
K = size(neuron.A, 2);
if ~exist('cell_id', 'var')
    cell_Id = 1;
else
    cell_id = max(1, cell_id-1);
end
tmp_white_list = false(K, 1); 
while true
    clc;
    %% guide
    fprintf('******************* Guide *******************\n');
    fprintf('input for cell ID: \n');
    fprintf('\t0: stop this iteration\n');
    fprintf('\t-k: back to the previous k-th neuron\n');
    fprintf('keep the neuron?\n');
    fprintf('\ty: yes\n\tn: no\n');
    fprintf('is it a soma?\n');
    fprintf('\ty: yes\n\tn: no\n');
    fprintf('add neurons to the white list\n');
    fprintf('\ty: yes\n\tn: no\n');
    fprintf('*********************************************\n');
    
    % show neuron
    cell_id = min(K, cell_id+1);
    ease_show_2p_neuron;
    drawnow; 
    
    % choose cell id
    commandwindow;
    temp = input('Cell ID: ');
    if temp==0
        break;
    elseif temp<0
        cell_id = max(1, cell_id+temp);
        ease_show_2p_neuron;
        drawnow; 
    end
    
    
    % mark neuron confidences
    commandwindow;
    temp = input('keep the neuron (y/n)?: ', 's');
    if strcmpi(temp, 'n')
        neuron.match_status.confidence(cell_id)= 0;
    end
    cb_btn_rate;
    
    % label neurons
    commandwindow;
    temp = input('is it a soma (y/n)?: ', 's');
    if strcmpi(temp, 'y')
        neuron.labels(cell_id) = 1; cb_btn_label;
    end
    
    % add to whitelist 
    commandwindow;
    temp = input('add it to white list (y/n)?: ', 's');
    if strcmpi(temp, 'y')
        tmp_white_list(cell_id) = true;
    elseif strcmpi(temp, 'n')
        tmp_white_list(cell_id) = false; 
    end
    if cell_id == K 
        break; 
    end
end

temp = cell2mat(neuron.match_status.em_ids(tmp_white_list)); 

if ~exist('white_list', 'var')
    white_list = temp; 
else
    white_list = union(white_list, temp); 
end