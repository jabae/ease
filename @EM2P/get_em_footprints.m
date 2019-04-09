function [Aem, segment_ids, segment_del] = get_em_footprints(obj, em_ids, scan_ids)

if ~exist('scan_ids', 'var')||isempty(scan_ids)
    scan_ids = obj.scan_id;
end
if isnumeric(em_ids)
    fields = {'segmentation', 'segment_id'};
    em_ids = reshape(em_ids, [], 1);
    n = length(em_ids);
    vals = num2cell([ones(n, 1)*obj.em_segmentation, em_ids]); %, ones(n,1), [1,1]);
    
    conditions = cell2struct(vals, fields, 2);
else
    conditions = em_ids; 
end

[segment_ids, idx_values] = fetchn(obj.rel_footprints & conditions, 'segment_id', 'idx_value');

% remove empty element
k_per_segment = cellfun(@length, idx_values);
ind = (k_per_segment==0);
segment_del = segment_ids(ind);
segment_ids(ind) = [];
idx_values(ind) = [];
k_per_segment(ind) = [];

if isempty(segment_ids)
    fprintf('no EM components selected.\n');
    Aem = {};
    return;
end
%% convert idx_values to a sparse matrix to store all footprints
idx_values = cell2mat(idx_values);
[idx, dims_new] = obj.convert_idx(idx_values(:,1));
K = length(segment_ids);

cumk = cumsum(k_per_segment);
idx_id = zeros(size(idx));
idx_id(cumk(1:(end-1))+1) = 1;
idx_id(1) = 1;
idx_id = cumsum(idx_id);

tmp_Aem = sparse(idx, idx_id, idx_values(:,2), prod(dims_new),  K);
nscan = length(scan_ids);
Aem = cell(nscan, 1);
for m=1:nscan
    k = scan_ids(m);
    idx = false(dims_new);
    idx(:, :, (1:obj.num_slices)+(k-1)*obj.num_slices) = true;
    Aem{m} = tmp_Aem(idx(:), :);
end
end