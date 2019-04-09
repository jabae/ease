function visualize_em_voxels(obj, em_id, new_figure)
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
if ~exist('new_figure', 'var') || isempty(new_figure)
    figure;
else
    hold on; 
end

% fetch EM voxels 
[indices] = fetch1(obj.rel_voxels & ...
    sprintf('segmentation=%d', obj.em_segmentation) & ...
    sprintf('segment_id=%d', em_id), 'indices'); 
[ind_r, ind_c, ind_z] = ind2sub(obj.dims_stack, indices);

% resolutions 
resolution = obj.range_2p ./ obj.dims_stack; 
x = (ind_c-1) * resolution(2); 
y = (ind_r-1) * resolution(1); 
z = (ind_z-1) * resolution(3); 

% show voxels 
plot3(x, y, z, '.k', 'markersize', 3); 