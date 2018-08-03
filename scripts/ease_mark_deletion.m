%% mark neuron to be deleted 
K = size(neuron.A, 2); 
if ~exist('ease_ind_delete', 'var') || (length(ease_ind_delete)~=K)
   ease_ind_delete = false(1, K);  
end

ease_ind_delete(cell_id) = true; 