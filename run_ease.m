if ~exist('ease', 'var')
    fi.usepkg('ease'); 
end 

%% determine project folder, dataset and database
ease.connect_database();
ease.select_data();

%% voxelize EM 
ease.voxelize_em(); 

%% choose FOV
show_fov = false;
ease.choose_FOV(show_fov);

%% crop and align the video data and the stack data 
% ease.load_stack();      % stack data
ease.load_video();      % video data

ease.rough_registration_video();  % rough registration 
ease.align_video_stack(); % fine registration within the cropped area 

% crop video
ease.crop_video;

% cell array for storing raw/denoised data and its summary statistics  
[Y_raw, Y_denoised] = ease.construct_Y();

%% load EM info
% ease.load_em(); 
% EM_info = ease.em_data.EM_info;
ease.get_em_boundaries(); 

%% project EM masks onto the scanning planes 
ease.project_em(); 
ease.construct_Aem(); 
ease.get_em_volume(); 

%% save the current configurations
ease.write_config();

%% GUI
if exist('with_GUI', 'var')  && with_GUI
    ease.startGUI();
end

%% create a flag indicating EASE is ready to be used.
flag_ease_running = true;
