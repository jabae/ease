%% load packages 
fi.usepkg({'yaml', ...  % load YAML files 
    'utils',...         % matlab util functions 
    'idl', ...          % easy access to scientific imaging data 
    'oasis' ...        % deconvolve and denoise temporal activity 
    }); 

%% add path
quiet = true; 
EASE_dir = fi.locate('ease', quiet);
addpath(EASE_dir);
addpath(fullfile(EASE_dir, 'scripts'));
addpath(fullfile(EASE_dir, 'packages'));

addpath(genpath(fullfile(EASE_dir, 'packages', 'microns_phase1_nda')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'pipeline', 'matlab')));

addpath(fullfile(EASE_dir, 'GUI'));
addpath(fullfile(EASE_dir, 'GUI', 'callbacks'));
addpath(fullfile(EASE_dir, 'functions'));

