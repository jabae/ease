function visualize_em_footprints(obj, em_id, new_figure)
%% visualize the mesh version of 1 EM segment 
%{
%}

%% inputs
%{
	obj: type; description
	em_id: UINT64; EM ID
	new_figure: boolean; create a new figure (new_figure=true) or not
%}

%% outputs
%{
%}

%% Author
%{
	Pengcheng Zhou 
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License 
%}


%%
if ~exist('new_figure', 'var') || isempty(new_figure) || new_figure
    figure;
else
    hold on; 
end

% fetch EM voxels 
[idx_value] = fetch1(obj.rel_footprints & ...
    sprintf('segmentation=%d', obj.em_segmentation) & ...
    sprintf('segment_id=%d', em_id) & ...
    sprintf('version=%d', obj.blur_version), 'idx_value'); 
[idx_r, idx_c, idx_z] = ind2sub(obj.dims_stack, idx_value(:,1));
v = idx_value(:, 2); 
if isempty(v)
    warning('empty footprints.'); 
    return; 
end 
% resolutions 
resolution = obj.range_2p ./ obj.dims_stack; 
FOV_2p = obj.FOV_stack; 
zvals = unique(idx_z); 
d1 = diff(FOV_2p(1:2)) + 1; 
d2 = diff(FOV_2p(3:4)) + 1; 
[xx, yy] = meshgrid((FOV_2p(3):FOV_2p(4))*resolution(2), ...
    (FOV_2p(1):FOV_2p(2))*resolution(1));
zz = ones(d1, d2)*resolution(3); 
%% 

for m=1:length(zvals)
   ind = (idx_z==zvals(m)); 
   img = nan(d1, d2); 
   tmp_idx = sub2ind([d1, d2], idx_r(ind)-FOV_2p(1)+1, idx_c(ind)-FOV_2p(3)+1);  
   img(tmp_idx) = v(ind); 
   surf(xx, yy, zz*zvals(m), img, 'edgecolor', [1,1,1]*0.9, 'edgealpha', 0.05); 
   hold on; 
end
caxis([0, max(v(:))*0.8]);
