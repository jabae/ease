function save_MF3D(obj)
% check the inputs
mscan = obj.scan_id; 
mblock = obj.block_id; 

var_name = sprintf('neuron_scan%d_block%d', mscan, mblock);

%% determine the matfile for saving the results
FOV_ = obj.FOV;
matfile_mf3d = fullfile(obj.output_folder, ...
    sprintf('neurons_%d_%d_%d_%d.mat',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));

%% save the results
flag_created(mscan, mblock) = true; %#ok<NASGU>
neuron = evalin('base', 'neuron'); 
eval(sprintf('%s = neuron;', var_name));
save(matfile_mf3d, var_name, '-append');
fprintf('The results of the current sources extraction were saved\n'); 
end