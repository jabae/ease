% add the packages used in this pipeline
addpath(fullfile(EASE_dir, 'packages'));
addpath(genpath(fullfile(EASE_dir, 'packages', 'microns_phase1_nda')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'pipeline', 'matlab')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'ta3')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'polygon2voxel')));


%% connect to the database
setenv('DJ_HOST','ninai.cluster-chjk7zcxhsgn.us-east-1.rds.amazonaws.com:3306')
setenv('DJ_USER','pczhou')
setenv('DJ_PASS','lilliput')
dj.conn()

% get all segments containing more than 15 fragments.
[segment_ids, n_fragments] = fetchn(aggr(ta3.Mesh,ta3.MeshFragment,'count(*)->n') & 'n>=10', ...
    'segment_id', 'n');
[segment_ids_all,n_vertices_all]=fetchn(ta3.Mesh,ta3.MeshFragment,'segment_id','sum(n_vertices)->total_n_vertices');
ind = find(n_vertices_all>6000);
segment_ids = segment_ids_all(ind);
n_vertices = n_vertices_all(ind);

%% fetch, voxelize and save
[vertices, faces] = fetchn(ta3.MeshFragment &  sprintf('segment_id=%d', id), 'vertices', 'triangles');

