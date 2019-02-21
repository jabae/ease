
%% collect EM information
tmp_file = fullfile(ease.data_folder, ease.matfile_em);
if ~exist(tmp_file, 'file')
    [segment_id, indices] = fetchn(rel_voxels & ...
        sprintf('segmentation=%d', ease.em_segmentation),...
        'segment_id', 'indices');
    EM_info = [segment_id, cellfun(@length, indices)];
    
    % determine EM ranges
    [y, x, z] = ind2sub(options.dims_2p, unique(cell2mat(indices)));
    em_ranges = cell(options.dims_2p(3), 1);
    for zz=1:options.dims_2p(3)
        ind = (z==zz);
        if ~isempty(ind)
            tmp_y = y(ind);
            tmp_x = x(ind);
            k = boundary(tmp_x,tmp_y,0);
            em_ranges{zz} = [tmp_y(k), tmp_x(k)];
        end
        if mod(zz, 10)==0
            disp(zz);
        end
    end
    
    
    EM_data.EM_info = EM_info;
    EM_data.options = options;
    EM_data.em_ranges = em_ranges;
    
    save(fullfile(ease.data_folder, ease.matfile_em), 'EM_info', 'options', 'em_ranges',...
        '-v7.3');
else
    fprintf('The EM segments has been selected and saved into \n\t%s.\n', tmp_file);
    fprintf('\tDelete the file first if you want to select data again\n');
end