snrs_all = zeros(ease.num_scans, ease.num_blocks, 1000); 
indices_all = zeros(ease.num_scans, ease.num_blocks, 1000); 
ids_all = zeros(1000, 1); 
k_neurons = 0; 
fprintf('Collecting neurons extracted in all scans\n'); 
for scan_id=1:5 
    for block_id=1:3 
        neuron = neurons_all{scan_id, block_id}; 
        
        em_ids = cell2mat(neuron.match_status.em_ids); 
        snrs = var(neuron.C, 0, 2) ./ var(neuron.C_raw-neuron.C, 0, 2);
        for m=1:length(em_ids)
            em_id = em_ids(m);
            idx = find(ids_all(1:k_neurons)==em_id);
            if isempty(idx) % add this neuron
                k_neurons = k_neurons + 1;
                ids_all(k_neurons) = em_id;
                idx = k_neurons;
            end
            snrs_all(scan_id, block_id, idx) = snrs(m);
            indices_all(scan_id, block_id, idx) = m; 
        end
    end
end
collections.ids_all = ids_all(1:k_neurons); 
collections.snrs_all = sparse(reshape(snrs_all(:, :, 1:k_neurons), [], k_neurons)); 
collections.indices_all = sparse(reshape(indices_all(:, :, 1:k_neurons), [], k_neurons)); 
fprintf('Done!\n'); 