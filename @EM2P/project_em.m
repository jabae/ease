function project_em(obj, use_parallel)
%% EM all EM voxels into the same scanning plane of the video data
%{%}

%% inputs
%{
	obj: class object
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


%% project all EM meshes
fprintf('project EM voxels onto each 2p scanning plane\n');
assignin('base', 'options', obj.set_projections_options()); % create a struct variable storing parameters

% check whether we should add this z value to the database 
rel = eval(sprintf('%s.Zblur;', obj.dj_name)); 
temp = fetchn(rel, 'zblur'); 
if ~any(temp==obj.em_zblur)
    key = struct('zblur', obj.em_zblur); 
    rel.insert(key); 
end 
if exist('use_parallel', 'var') && (use_parallel)
    parpopulate(obj.rel_footprints, sprintf('segmentation=%d and zblur=%d', obj.em_segmentation, obj.em_zblur));
else
    populate(obj.rel_footprints, sprintf('segmentation=%d and zblur=%d', obj.em_segmentation, obj.em_zblur));
end

fprintf('done!\n');
