%% load EM info from matlab file
if isempty(ease.em_data)
    ease.load_em();
end

%% create masks for EM volumes
if ~exist('EM_masks', 'var')
    ease.get_em_masks();
end

%% project all EM segments to the scanning planes
ease.get_Aem_scan();
pixels_em = neuron.reshape(sum(Aem, 2)>0, 3);
pixels_em = imerode(pixels_em, strel('disk', 2)); 
pixels_em = neuron.reshape(pixels_em>0, 1); 