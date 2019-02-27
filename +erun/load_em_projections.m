%% load EM info from matlab file
if isempty(ease.em_data)
    ease.load_em();
end

%% create masks for EM volumes
if ~exist('EM_masks', 'var')
    ease.get_em_masks();
end

%% project all EM segments to the scanning planes
ease.get_Aem_scan();        % load the projected EM footprints
d1 = diff(ease.FOV(1:2)) + 1; 
d2 = diff(ease.FOV(3:4)) + 1;
d3 = ease.num_slices; 

pixels_em = reshape(full(sum(Aem, 2)>0), d1, d2, d3);
pixels_em = imerode(pixels_em, strel('disk', 2)); 
pixels_em = sparse(reshape(pixels_em>0, [], 1)); 