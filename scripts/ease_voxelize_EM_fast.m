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

%% parameters
voxel_em = [25, 25, 20]/1000;   % the desired voxel size for voxelizing
% em meshes (not the actual voxel size)
dims_2p = [512, 512, 310];      % dimension of the stack data
dims_video = [256, 256];        % dimension of the functional imaging data
range_2p = [400, 400, 310];     % spatial range of 2p stack
voxels_2p = range_2p ./ dims_2p;  % voxel size of 2p stack data
scale_factor = 3.58/4;          % the factor for scaling EM coordinates.
% this is used because of some pre-processing
% error
% load the transformation
if exist(matfile_transformation, 'file')
    load(matfile_transformation);
else
    ease_get_transformation;
end

% create a struct variable storing parameters
options.voxel_em = voxel_em;
options.dims_2p = dims_2p;
options.range_2p = range_2p;
options.A = A_convert;
options.offset = offset;
options.scale_factor = scale_factor;
MAX_LENGTH = 64;

%% fetch all em segments and save them into a matfile
if ~exist(matfile_em, 'file')
    % open a mat file for saving information
    em_data = matfile(matfile_em, 'Writable', true);
    save(matfile_em, 'options', 'MAX_LENGTH', '-v7.3');
else
    em_data = matfile(matfile_em, 'Writable', true);
end

% get all segments containing more than 15 fragments.
% [segment_ids, n_fragments] = fetchn(aggr(ta3.Mesh,ta3.MeshFragment,'count(*)->n') & 'n>=10', ...
%     'segment_id', 'n');
[segment_ids_all,n_vertices_all]=fetchn(ta3.Mesh,ta3.MeshFragment,'segment_id','sum(n_vertices)->total_n_vertices');
ind = find(n_vertices_all>2500);
segment_ids = segment_ids_all(ind); 
n_vertices = n_vertices_all(ind); 

K = length(segment_ids);
nticks = min(K, 100);

num_array = ceil(K/MAX_LENGTH);
dims_old = dims_2p;
sub0 = [210, 70, 55];
dims_new = [135, 280, 170];
EM_masks = zeros(prod(dims_new), num_array, 'like', uint64(0));
EM_info = zeros(K, 4);

% extract the whole database 
if ~exist('EM_fragment.mat', 'file')
    [seg_all, vertices_all, faces_all] = fetchn(ta3.MeshFragment,'segment_id', 'vertices', 'triangles');
    save EM_fragments seg_all vertices_all faces_all -v7.3;
else
   load EM_fragments.mat;  
end
% voxelize neurons one by one 
for m=1:K
    % create a progress bar
    if m==1
        fprintf('loading EM data\n');
        for mm=1:nticks
            fprintf('|');
        end
        fprintf('\n');
    end
    id = segment_ids(m);
    
    % corrsponding bit
    ii = ceil(m/MAX_LENGTH);
    jj = 2^(mod(m-1, MAX_LENGTH));
    
    % fetch, voxelize and save
    ind = (seg_all==id); 
    vertices = vertices_all(ind); 
    faces = faces_all(ind); 
%     [vertices, faces] = fetchn(ta3.MeshFragment &  sprintf('segment_id=%d', id), 'vertices', 'triangles');
    
    [subs_2p, ~] = mesh2volume_fast(vertices, faces, options);
%     em_segment = struct('vertices', {vertices}, 'faces', {faces});
%     eval(sprintf('em_data.cell_%d = em_segment;', id));
    
    ind_2p = sub2ind(dims_new, subs_2p(:,1)-sub0(1)+1, ...
        subs_2p(:,2)-sub0(2)+1, subs_2p(:,3)-sub0(3)+1);
    
    EM_masks(ind_2p, ii) = EM_masks(ind_2p, ii) + jj; % the m-th bit is set as 1.
    EM_info(m, :) = [id, length(faces), length(ind_2p), n_vertices(m)];
    if round(m*nticks/K)~=round((m+1)*nticks/K)
        fprintf('.');
    end
end
fprintf('\n');
em_data.EM_info = EM_info;

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
    disp(m); 
end

em_data.em_ranges = em_ranges;

%%
% em_ranges = cell(dims_2p(3), 1);
% for zz=1:dims_2p(3)
%     tmp_str = sprintf('slice_%d', zz);
%     if ismember(tmp_str, em_variables)
%         temp = eval([em_nam, '.', tmp_str]);
%         [y, x] = ind2sub(dims_2p(1:2), temp(:,1));
%         k = boundary(x,y,0); 
%         em_ranges{zz} = [y(k), x(k)]; 
% 
%     end
% end
% 
% em_data.em_ranges = em_ranges;



