function save_MF3D(obj)
% check the inputs
mscan = obj.scan_id; 
mblock = obj.block_id;


%% determine the matfile for saving the results
FOV_ = obj.FOV;
matfile_mf3d = fullfile(obj.output_folder, ...
    sprintf('neurons_%d_%d_%d_%d.mat',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));

if mblock==0
    var_name = sprintf('neuron_scan%d_all_block', mscan);
else
    var_name = sprintf('neuron_scan%d_block%d', mscan, mblock);
end

%% save the results
if mblock>0
flag_created(mscan, mblock) = true; %#ok<NASGU>
end
neuron = evalin('base', 'neuron'); 
neuron.compress(); 
eval(sprintf('%s = neuron;', var_name));
save(matfile_mf3d, var_name, '-append');

fprintf('The results of the current sources extraction were saved\n'); 
end