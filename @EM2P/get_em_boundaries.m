function get_em_boundaries(obj)
for mscan=1:obj.num_scans
    zs = obj.video_zvals_updated(mscan, :);
    if iscell(zs)
        zs = cell2mat(zs); 
    end
    temp = obj.em_data.em_ranges;
    em_volume = cell2mat(temp(zs));
    if ~isempty(em_volume)
        k = convhull(em_volume(:, 2), em_volume(:, 1));
        obj.em_boundary{mscan} = em_volume(k, :);
    else
        obj.em_boundary{mscan} = [];
    end
end
end
