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
if exist('use_parallel', 'var') && (~use_parallel)
    populate(obj.rel_footprints, sprintf('segmentation=%d', obj.em_segmentation));
else
    parpopulate(obj.rel_footprints, sprintf('segmentation=%d', obj.em_segmentation));
end

fprintf('done!\n');
