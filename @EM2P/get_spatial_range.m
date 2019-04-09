function spatial_range = get_spatial_range(obj)
%% function goal
%{
	details
%}

%% inputs
%{
	obj: type; description
%}

%% outputs
%{
%}

%% Author
%{
	Pengcheng Zhou
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	XXX License
%}

em_ranges = obj.em_data.em_ranges;
spatial_range = false(obj.d1, obj.d2, obj.num_scans*obj.num_slices);
for mscan=1:obj.num_scans
    for mslice = 1:obj.num_slices
        z = obj.video_zvals_updated(mscan, mslice);
    end
    em_range = em_ranges{z};
    if isempty(em_range)
        continue;
    end
    idx = sub2ind(obj.dims_stack, em_range(:,1), em_range(:,2), ones(size(em_range,1),1)*z);
    [idx_new, ~] = obj.convert_idx(idx);
    spatial_range(idx_new) = true;
end
