function voxelize_em(obj, use_parallel)
%% voxelize all EM meshes into voxels in the same space of 2p stack data 
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


% metainfo = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
% temp = metainfo.(obj.data_name);
% try  % voxelized meshes 
%     eval(sprintf('rel_voxels=%s.VoxelizedMesh;', temp.datajoint_name));
% catch 
%     error('errors in using %s.VoxelizedMesh\n', temp.datajoint_name); 
% end

%% voxelize all EM meshes
assignin('base', 'options', obj.set_voxelization_options()); % create a struct variable storing parameters
if exist('use_parallel', 'var') && (~use_parallel)
    populate(obj.rel_voxels, sprintf('version=%d', obj.em_segmentation));
else
    parpopulate(obj.rel_voxels, sprintf('version=%d', obj.em_segmentation));  
end 
% collect EM range information 
obj.collect_em_info(); 

