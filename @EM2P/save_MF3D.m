function save_MF3D(obj)
% check the inputs
mscan = obj.scan_id;
mblock = obj.block_id;


%% determine the matfile for saving the results
FOV_ = obj.FOV;
matfile_mf3d = fullfile(obj.output_folder, ...
    sprintf('neurons_%d_%d_%d_%d.mat',...
    FOV_(1), FOV_(2), FOV_(3), FOV_(4)));
folder_MF3D = matfile_mf3d(1:(end-4));

if mblock==0
    var_name = sprintf('neuron_scan%d_all_blocks', mscan);
else
    var_name = sprintf('neuron_scan%d_block%d', mscan, mblock);
end

%% save the results
if mblock>0
    flag_created(mscan, mblock) = true; %#ok<NASGU>
end
neuron = evalin('base', 'neuron');
neuron.compress();
if ~exist(folder_MF3D, 'dir')
    mkdir(folder_MF3D); 
end
file_name = fullfile(folder_MF3D, sprintf('%s.mat', var_name)); 
eval(sprintf('%s = neuron;', var_name));
save(file_name, var_name, '-v7.3'); 

% save(matfile_mf3d, var_name, '-append');

fprintf('The results of the current sources extraction were saved\n');
end