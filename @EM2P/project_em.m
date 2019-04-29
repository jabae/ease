function project_em(obj)
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
options = obj.set_projections_options();
assignin('base', 'options', options); % create a struct variable storing parameters

populate(obj.rel_footprints, sprintf('segmentation=%d and version=%d', obj.em_segmentation, options.version));

fprintf('done!\n');
