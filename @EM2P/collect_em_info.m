function collect_em_info(obj)
%% Collect EM information 
%{
	get all EM IDs and the spatial ranges of the EM volume 
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


%% collect EM information
tmp_file = fullfile(obj.data_folder, obj.matfile_em); % the matfile for storing results 
if ~exist(tmp_file, 'file')
    fprintf('Collecting EM info...\n'); 
    [segment_id, n_voxels] = fetchn(obj.rel_voxels & ...
        sprintf('version=%d', obj.em_segmentation),...
        'segment_id', 'n_voxels'); 
    EM_info = [segment_id, n_voxels];
    
    dims_2p = obj.dims_stack; 
    
    % determine EM ranges
    % in the future, this part should be done patch by patch for
    % scalability concerns 
    [indices] = fetchn(obj.rel_voxels & ...
        sprintf('version=%d', obj.em_segmentation),...
        'indices');
    [y, x, z] = ind2sub(dims_2p, unique(cell2mat(indices)));
    em_ranges = cell(dims_2p(3), 1);
    for zz=1:dims_2p(3)
        ind = (z==zz);
        if ~isempty(ind)
            tmp_y = y(ind);
            tmp_x = x(ind);
            k = boundary(tmp_x,tmp_y,0);
            em_ranges{zz} = [tmp_y(k), tmp_x(k)];
        end
        if mod(zz, 10)==0
            fprintf('%d...', zz);
        end
    end
    fprintf('\n'); 
    
    obj.em_info = EM_info; 
    obj.em_ranges = em_ranges; 
    save(fullfile(obj.data_folder, obj.matfile_em), 'EM_info', 'em_ranges',...
        '-v7.3');
    fprintf('Done.\n'); 
else
    fprintf('The EM info has been saved into \n %s.\n', tmp_file);
    fprintf('Delete the file first if you want to select data again\n');
    temp = load(tmp_file); 
    obj.em_info = temp.EM_info; 
    obj.em_ranges = temp.em_ranges; 
end
