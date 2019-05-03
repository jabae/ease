%% load packages 
fi.usepkg({'yaml', ...  % load YAML files 
    'utils',...         % matlab util functions 
    'idl', ...          % easy access to scientific imaging data 
    'oasis' ...        % deconvolve and denoise temporal activity 
    'datajoint', ...    % use datajoint 
    }); 

%% add path
quiet = true; 
EASE_dir = fi.locate('ease', quiet);
addpath(EASE_dir);
addpath(fullfile(EASE_dir, 'scripts'));
addpath(fullfile(EASE_dir, 'packages'));
addpath(fullfile(EASE_dir, 'functions'));

addpath(genpath(fullfile(EASE_dir, 'packages', 'microns_phase1_nda')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'pipeline', 'matlab')));

addpath(fullfile(EASE_dir, 'GUI'));
addpath(fullfile(EASE_dir, 'GUI', 'callbacks'));

%% create a class object 
ease = EM2P(); 
addpath(fullfile(ease.dir_project, 'schemas')); 