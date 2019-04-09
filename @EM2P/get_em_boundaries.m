function get_em_boundaries(obj)
%% function goal
%{
	details
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
em_ranges = obj.em_ranges;

for mscan=1:obj.num_scans
    zs = obj.video_zvals_updated(mscan, :);
    if iscell(zs)
        zs = cell2mat(zs); 
    end
    em_volume = cell2mat(em_ranges(zs));
    if ~isempty(em_volume)
        k = convhull(em_volume(:, 2), em_volume(:, 1));
        obj.em_boundary{mscan} = em_volume(k, :);
    else
        obj.em_boundary{mscan} = [];
    end
end
end
