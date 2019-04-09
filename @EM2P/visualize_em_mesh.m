function visualize_em(obj, em_id, new_figure)
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

% fetch EM meshes
[vertices, faces] = fetchn(obj.rel_mesh & ...
    sprintf('segmentation=%d', obj.em_segmentation) & ...
    sprintf('segment_id=%d', em_id), 'vertices', 'triangles');

% get the transformation matrix between EM space and 2P space
[A_convert, offset] = obj.get_transformation();
scale_factor = obj.em_scale_factor;

% show meshes 
for m=1:length(faces)
    vert = bsxfun(@plus, vertices{m} * scale_factor * A_convert, offset);
    trisurf(faces{m}+1, vert(:,1),...
        obj.range_2p(2)-vert(:,2), ...
        obj.range_2p(3)-vert(:,3), ...
        'edgecolor', 'none', 'facecolor', [1, 0.7, 0]);       hold on;
end
xlabel('X (um)'); 
ylabel('Y (um)'); 
zlabel('Z (um)'); 
