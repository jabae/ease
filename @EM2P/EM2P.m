classdef EM2P < handle
    
    % This class is a wrapper for processing EM and 2P imaging data jointly
    
    % Author: Pengcheng Zhou, Columbia University, 2018
    % zhoupc1988@gmail.com
    
    %% properties
    properties
        % YAML file for storing the configurations
        yaml_path = '';
        
        % data folder 
        data_folder = ''; %folder for storing intermediate results 
        % output folder
        output_folder = '';
        
        % fig folder 
        fig_folder = ''; % folder for saving figures 
        
        % video_folder 
        video_folder = ''; %folder for savign videos 
        
        % video imaging data
        matfile_video = ''; % file path
        video_data = [];    % the mat file including information of video data
        video_loader = {};  % data loader for accesing video data
        video_shifts = struct('ii', [], 'jj', []);       % shifts in relative to 2p stack data
        video_zvals = [];  % z values for each slice
        video_zvals_updated = [];   % updated z values for different planes.
        video_Fs = 15;      % frame rate
        video_T = 8781;     % number of frames.
        denoised_folder = ''; % folder of the denoising results for the cropped area
        raw_folder = '';     % folder of the raw data of the cropped area
        use_denoise = false; % use the denoised data for running CNMF
        
        % mat file for storing all 2p stack data
        matfile_stack = '/data/lab/Dropbox/Tolias_PC/em_2p_pipeline/data/functional_data/stack_2p.mat';
        stack_shifts = {};
        
        % coordinate conversion between EM and 2p stack data
        registration_csv = ''; 
        matfile_transformation = '/data/lab/Dropbox/Tolias_PC/em_2p_pipeline/coor_convert.mat';
        
        % FOV information
        d1 = [];
        d2 = [];
        d3 = [];
        extra_margin = 5; % extra margins surrounding the EM volume
        FOV = [];
        FOV_stack = [];
        aligned_images = [];      % overlap of stack images and mean of video images to show the performances of alignment
        align_max_zshift = 8;     %
        ssub = 2;       % downsampling factor between stack data and video data 
        
        % data information
        num_scans = 8;          % number of total scans
        num_slices = 3;         % number of the imaged slices per scan
        num_blocks = 3;         % number of temporal blocks per slice
        dims_stack = [512, 512, 310];      % dimension of the stack data
        dims_video = [256, 256];        % dimension of the functional imaging data
        range_2p = [400, 400, 310];     % spatial range of 2p stack
        
        % default calcium imaging data to be analyzed.
        scan_id = 1;
        slice_id = 1;
        block_id = 1;
        
        % initialization 
        options_init = struct();    % options for initialization 
        init_method = 'regression';     % the default initialization method
        show_init = true;   % interactively running the initialization procedure
        pause_init = false; % show the intermediate results during the initialization
        maxIter_init = 1;  % nubmer of iterations in the initialization step
        
        show_em_only = false; % show EM component when you display the merge of the EM mask and CNMF mask
        
        
        match_mask = true;    % use the spatial mask (true) or the spatial footprint (false) in
        % the step of computing the matching score between EM & CNMF component
        
        % EM data
        matfile_em = ''; % matfile for storing all EM data
        em_segmentation = 1; 
        em_data = [];
        em_variables = [];
        em_info = [];
        em_ranges = [];
        K_em = 0;
        em_shifts = [];     % shifts for EM volumes
        % pre-load EM data for faster speed. (default: true)
        em_load_flag = true;
        em_zblur = 6;          % the projection of each EM semgent onto one slice takes the nearby (-blur_size:blur_size) planes
        PACK_SIZE = 64;     % pack size for compressing binary vectors.
        em_boundary = [];   % EM boundaries
        
        % method for computing the matching scores
        score_method = 'corr' ;
        
        % create new MF3D class object for each scan 
        create_new = false; 
        
        % GUI
        gui = [];
        nam_show = 'cn';    % show correlation image
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = EM2P(path_yaml_file)
            try
                EASE_folder = fileparts(which('ease_setup.m'));
                obj.read_config(path_yaml_file);
            catch
                obj.yaml_path = fullfile(EASE_folder, 'config_ease.yaml');
                obj.output_folder = '/data/lab/Dropbox/softwares/EASE/results';
                obj.matfile_video = '/data/paninski_lab/Tolias_lab/functional_data_analysis/functional_data.mat';
                obj.matfile_em = '/data/paninski_lab/Tolias_lab/EM/em.mat';
                obj.denoised_folder = '/data/paninski_lab/Tolias_lab/cropped_denoised/';
                obj.raw_folder = '/data/paninski_lab/Tolias_lab/cropped/';
            end
        end
        
        %% load and write configurations
        read_config(obj, path_file);
        write_config(obj, path_file);
        
        %% load 2p stack data
        load_stack(obj);
        
        %% load EM data
        load_em(obj);
        
        %% load functional video data
        load_video(obj);
        
        %% create data loader
        [dl_Yr, dl_Yd] = create_dataloader(obj, scan_id, block_id, T);
        
        %% choose FOV
        choose_FOV(obj, show_fov);
        
        %% align the video data and the 2p stack data
        align_video_stack(obj);
        
        %% extract EM volumes given the scan ID
        Aem_proj = extract_em_segments(obj, mscan, mslice, ssub);
        
        %% get EM masks in a specified scan ID
        Aem = get_Aem_scan(obj, mscan, ssub);
        
        %% create/load neuron objects for all video data
        neuron = get_MF3D(obj, T, create_new);
        
        %% get all calium imaging data
        summary_images = calculate_summary_images(obj, mscan, mblock)
        
        %% load video data into memory
        neuron = load_video_mem(obj, mscan, mblock, create_new);
        
        %% create masks for EM volumes
        get_em_masks(obj);
        
        %% save the current object 
        function save(obj)
            tmp_file = fullfile(obj.output_folder, 'ease.mat');
            save(tmp_file, 'obj');
        end
        
        %% save the current results of the sources extraction
        save_MF3D(obj);
        
        %% GUI
        startGUI(obj);
  
        %% load the projection of Aem onto the 2p stack 
        Aem_proj = project_em_to_stack(obj); 

        %% update initialization options 
        update_init_options(obj); 
        
        %% get EM boundaries
        get_em_boundaries(obj);  
    
    end
    
    
end

