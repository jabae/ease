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
pixels_em = (sum(Aem, 2)>0);
