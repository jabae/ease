%% load packages 
fi.usepkg({'yaml', ...  % load YAML files 
    'utils',...         % matlab util functions 
    'idl', ...          % easy access to scientific imaging data 
    'oasis' ...        % deconvolve and denoise temporal activity 
    }); 

%% add path
EASE_dir = fi.locate('ease');
addpath(EASE_dir);
addpath(fullfile(EASE_dir, 'scripts'));
addpath(fullfile(EASE_dir, 'packages'));
addpath(fullfile(EASE_dir, 'GUI'));
addpath(fullfile(EASE_dir, 'GUI', 'callbacks'));
addpath(fullfile(EASE_dir, 'functions'));
