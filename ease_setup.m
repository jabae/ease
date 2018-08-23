%% add path
EASE_dir = fileparts(which('ease_setup.m'));
addpath(EASE_dir);
addpath(fullfile(EASE_dir, 'scripts'));
addpath(fullfile(EASE_dir, 'packages'));
addpath(fullfile(EASE_dir, 'packages', 'mesh2volume'));
addpath(fullfile(EASE_dir, 'packages', 'reduce_poly'));
addpath(fullfile(EASE_dir, 'GUI'));
addpath(fullfile(EASE_dir, 'GUI', 'callbacks'));
addpath(fullfile(EASE_dir, 'functions'));

% addpath(genpath(fullfile(EASE_dir, 'packages', 'microns_phase1_nda')));
% addpath(genpath(fullfile(EASE_dir, 'packages', 'pipeline', 'matlab')));
% addpath(genpath(fullfile(EASE_dir, 'packages', 'ta3')));
% addpath(genpath(fullfile(EASE_dir, 'packages', 'polygon2voxel')));

%% setup imaging data loader 
if isempty(which('IDL.m'))
    try
        run(fullfile(EASE_dir, '..', 'ImagingDataLoader', 'idl_setup.m'));
    catch
        disp('install ImagingDataLoader https://github.com/zhoupc/ImagingDataLoader.git\n'); 
        % download imaging data loader and install it
        % TBD
    end
end

%% setup oasis_matlab 
if isempty(which('deconvolveCa.m'))
    try
        run(fullfile(EASE_dir, '..', 'OASIS_matlab', 'setup.m'));
    catch
        disp('install oasis-matlab https://github.com/zhoupc/OASIS_matlab.git\n'); 
        % download imaging data loader and install it
        % TBD
    end
end

%% install YAML-matlab 
addpath(fullfile(EASE_dir, '..', 'yamlmatlab')); 

%% warning off 
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle'); 