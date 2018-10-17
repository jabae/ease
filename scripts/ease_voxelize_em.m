% setup path and connect to a database 
addpath(fullfile(EASE_dir, 'packages'));
addpath(genpath(fullfile(EASE_dir, 'packages', 'polygon2voxel')));

temp = dj.conn;
if ~temp.isConnected
    ease_connect_database;
end
%% get the transformation matrix between EM space and 2P space 
if ~exist(fullfile(ease.data_folder, ease.matfile_transformation), 'file')
    ease_get_transformation; 
else
    load(fullfile(ease.data_folder, ease.matfile_transformation)); 
end

% create a struct variable storing parameters
options.voxel_em = [256, 256, 240]/1000; % voxelize EM data with this resolution 
options.dims_2p = ease.dims_stack;
options.range_2p = ease.range_2p;
options.A = A_convert;
options.offset = offset;
options.user_parallel = false; 
if strcmpi(data_name, 'pinky40')  % convert the EM unit to um
    options.scale_factor = 0.001*3.58/4;      
else
    options.scale_factor = 0.001;      
end

%% get all EM IDs and their number of vertices 
if strcmpi(data_name, 'pinky40')
    [segment_ids,n_vertices]=fetchn(ta3.Mesh,ta3.MeshFragment,...
        'segment_id','sum(n_vertices)->total_n_vertices');
else
    rel = ta3p100.Mesh();
    [segment_ids, n_vertices] = rel.fetchn('segment_id', 'n_vertices'); 
end

% order segment_ids according to the number of vertices
[n_vertices, idx] = sort(n_vertices, 'descend'); 
segment_ids = segment_ids(idx); 

%% create voxelized EM components and save them into database 
populate(ta3p100.VoxelizedMesh); 

%% select one neuron and check the visualization 
if strcmpi(data_name, 'pinky40')
    id = 24927722; 
    str_id = sprintf('segment_id=%d', id); 
    [vertices, faces] = fetchn(ta3.MeshFragment & str_id, 'vertices', 'triangles');
    
    [subs_2p, ~] = mesh2volume(vertices, faces, options);
    ind = sub2ind(options.dims_2p, subs_2p(:,1), subs_2p(:,2), subs_2p(:,3));
    V = zeros(options.dims_2p);
    V(ind) = 1;
else   
    id = segment_ids(40);    
    str_id = sprintf('segment_id=%d', id);
    indices = fetch1(ta3p100.VoxelizedMesh & str_id, 'indices'); 
    V = zeros(options.dims_2p); 
    V(indices) = 1; 
end


%% save the EM masks slice by slices
EM_masks = reshape(EM_masks, [dims_new, num_array]);
em_ranges = cell(dims_2p(3), 1);
for m=1:dims_new(3)
    temp = reshape(EM_masks(:, :, m, :),  dims_new(1) * dims_new(2), num_array);
    [ii, jj, v] = find(temp);
    if ~isempty(ii)
        [y, x] = ind2sub(dims_new(1:2), ii);
        % get the spatial ranges
        k = boundary(x, y, 0);
        em_ranges{m+sub0(3)-1} = [y(k)+sub0(1)-1, x(k)+sub0(2)-1];
        
        % save the nonzero voxels
        ii = sub2ind(dims_2p(1:2), y+sub0(1)-1, x+sub0(2));
        temp = [ii, jj, v]; % save a sparse matrix
        eval(sprintf('em_data.slice_%d = temp;', m+sub0(3)-1));    
    end
end

em_data.em_ranges = em_ranges;

%%
em_ranges = cell(dims_2p(3), 1);
for zz=1:dims_2p(3)
    tmp_str = sprintf('slice_%d', zz);
    if ismember(tmp_str, em_variables)
        temp = eval([em_nam, '.', tmp_str]);
        [y, x] = ind2sub(dims_2p(1:2), temp(:,1));
        k = boundary(x,y,0); 
        em_ranges{zz} = [y(k), x(k)]; 

    end
end

em_data.em_ranges = em_ranges;

