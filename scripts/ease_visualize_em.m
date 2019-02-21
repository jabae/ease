% setup path and connect to a database
if ~exist('dj_connected', 'var') || ~dj_connected
    ease_connect_database;
end

%% 
if ~exist('new_figure', 'var')
    new_figure = true; 
end

%% get the transformation matrix between EM space and 2P space
[A_convert, offset] = ease.get_transformation(); 
scale_factor = ease.em_scale_factor; 

% create a struct variable storing parameters
if strcmpi(data_name, 'pinky40')  % convert the EM unit to um
    [vertices, faces] = fetchn(ta3.MeshFragment & 'segmentation=1' &...
        sprintf('segment_id=%d', em_id),...
        'vertices', 'triangles');
    for m=1:length(faces)
        vert = bsxfun(@plus, vertices{m} * scale_factor * A_convert, offset);
        k_vert = size(vert, 1);
        trisurf(faces{m}+1, vert(:,1),...
            ease.range_2p(2)-vert(:,2), ...
            ease.range_2p(3)-vert(:,3), ...
            'edgecolor', 'none', 'facecolor', [1, 0.7, 0]);       hold on;
    end
else
    [vertices, faces] = fetch1(ta3p100.Mesh & ...
        sprintf('segmentation=%d', ease.em_segmentation) & ...
        sprintf('segment_id=%d', em_id), 'vertices', 'triangles');
    if new_figure
        figure;
    end
    vert = bsxfun(@plus, vertices * scale_factor * A_convert, offset);
    k_vert = size(vert, 1);
    if ~exist('tmp_color', 'var')
        tmp_color = [1, 0.7, 0]; 
    end
      trisurf(faces+1, vert(:,1),...
        ease.range_2p(2)-vert(:,2), ...
        ease.range_2p(3)-vert(:,3), ...
        'edgecolor', 'none', 'facecolor', tmp_color);
end

%% collect EM information
