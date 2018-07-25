function [Aem_blur, dims] = extract_Aem_mask(data, ind)
dims_Aem = data.dims_Aem;
dims = data.dims;

d3 = data.d3_2p;
Aem_blur = zeros(dims_Aem(1), length(ind), d3);
for m=1:d3
    temp = data.Aem_all{m};
    if ~isempty(temp)
        Aem_blur(:, :, m) = full(temp(:, ind));
    end
end

%% blur, change order and compress
Aem_blur = permute(Aem_blur, [1, 3, 2]);
Aem_blur = reshape(Aem_blur, [], size(Aem_blur, 3));
end