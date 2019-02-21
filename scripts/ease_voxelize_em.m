%% connect to a database 
if ~exist('dj_connected', 'var') || ~dj_connected
    ease_connect_database;
end

%% %% create a struct variable storing parameters
options = ease.set_voxelization_options(); 

%% voxelize all EM meshes 
parpopulate(rel_voxels);