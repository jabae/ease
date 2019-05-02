classdef EM2P < handle
    
    % This class is a wrapper for processing EM and 2P imaging data jointly
    
    % Author: Pengcheng Zhou, Columbia University, 2018
    % zhoupc1988@gmail.com
    
    %% properties
    properties
        % YAML file for storing the configurations
        yaml_path = '';
        
        % data name
        data_name = '';
        
        % project folder
        dir_project = '';
        
        % datasets
        datasets_list = {};
        
        % databases
        databases_list = {};
        
        % data folder
        data_folder = ''; %folder for storing intermediate results
        
        % output folder
        output_folder = '';
        
        % script folder
        script_folder = '';
        
        % fig folder
        fig_folder = ''; % folder for saving figures
        
        % video_folder
        video_folder = ''; %folder for savign videos
        
        % video imaging data
        matfile_video = 'functional_data.mat'; % file path
        video_data = [];    % the mat file including information of video data
        video_loader = {};  % data loader for accesing video data
        video_shifts = struct('ii', [], 'jj', []);       % shifts in relative to 2p stack data
        video_zvals = [];  % z values for each slice
        video_zvals_updated = [];   % updated z values for different planes.
        video_Fs = 15;      % frame rate
        video_T = 0;     % number of frames.
        denoised_folder = 'cropped_denoised_video'; % folder of the denoising results for the cropped area
        raw_folder = 'cropped_raw_video';     % folder of the raw data of the cropped area
        use_denoise = false; % use the denoised data for running CNMF
        normalize_data = true; % normalize data or not
        summary_images = [];
        
        % mat file for storing all 2p stack data
        matfile_stack = 'stack_2p.mat';
        stack_shifts = {};
        
        % coordinate conversion between EM and 2p stack data
        registration_csv = 'registration.csv';
        matfile_transformation = 'coor_convert.mat';
        transformation = struct('A_convert', [], 'offset', []);
        
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
        
        
        % EM data
        matfile_em = 'em.mat'; % matfile for storing all EM data
        em_segmentation = 1;
        em_data = [];
        em_variables = [];
        em_info = [];
        em_ranges = [];
        em_volume = [];     % binary variable to indicate whether a pixel is within the EM volume or not.
        K_em = 0;
        em_shifts = [];     % shifts for EM volumes
        % pre-load EM data for faster speed. (default: true)
        em_load_flag = true;
        em_zblur = 8;          % the projection of each EM semgent onto one slice takes the nearby (-blur_size:blur_size) planes
        em_boundary = [];   % EM boundaries
        em_scale_factor = 0.001;  % a factor to map EM unit to um
        show_em_only = false; % show EM component when you display the merge of the EM mask and CNMF mask
        
        % method for computing the matching scores
        score_method = 'corr' ;
        
        % create new MF3D class object for each scan
        create_new = false;
        
        % GUI
        gui = [];
        nam_show = 'cn';    % show correlation image
        
        % datajoint
        dj_connected = false;
        
        % datajoint
        dj_name = [];
        rel_mesh = [];
        rel_voxels = [];
        rel_footprints = [];
        
        % blur version
        blur_version = [];
        
        % stimuli 
        conditions = []; 
    end
    
    %% methods
    methods
        %% constructor and options setting
        function obj = EM2P(path_yaml_file, start_new_project)
            if exist('path_yaml_file', 'var') && exist(path_yaml_file, 'file')
                % start from a given yaml configuration
                obj.read_config(path_yaml_file);
            else
                % start from a default yaml configuration
                EASE_folder = fi.locate('ease', true);
                path_yaml_file = fullfile(EASE_folder, 'config.yaml');
                obj.read_config(path_yaml_file, true);
                obj.yaml_path = '';   % don't overwrite the default config
            end
            
            %
            if ~exist('start_new_project', 'var') || isempty(start_new_project)
                start_new_project = false;
            end
            
            % determine project folder
            if isempty(obj.dir_project)
                tmp_file = fullfile(fi.locate('ease', true), '.projects.mat');
                if exist(tmp_file, 'file')
                    load(tmp_file, 'projects');
                else
                    projects = {};
                end
                if isempty(projects)
                    obj.init_project();
                elseif length(projects)==1&& exist(projects{1}, 'dir') && ~start_new_project
                    obj.dir_project = projects{1};
                    fprintf('project folder: %s\n', projects{1});
                else
                    fprintf('there are following projects: \n');
                    
                    n_projects = length(projects);
                    for m=1:length(projects)
                        fprintf('--%d: %s\n', m, projects{m});
                    end
                    fprintf('\n***************** GUIDE *****************\n');
                    fprintf('type  0: create a new one\n');
                    fprintf('type  k: choose the k-th project shown above\n');
                    fprintf('type -k: delete the k-th project shown above\n');
                    fprintf('type k+1: select an existing project folder\n');
                    fprintf('*****************  END  *****************\n');
                    
                    while true
                        tmp_k = input('\nchoose: ');
                        if tmp_k==0
                            obj.init_project();
                            break;
                        elseif tmp_k > n_projects
                            [dir_name] = uigetdir();
                            projects{n_projects+1} = dir_name;
                            obj.dir_project = dir_name;
                            break;
                        elseif tmp_k>0  %select a project
                            obj.dir_project = projects{tmp_k};
                            break;
                        elseif abs(tmp_k)<= n_projects  % delete a project
                            projects{abs(tmp_k)} = [];
                        end
                    end
                    projects(cellfun(@isempty, projects)) = [];
                    save(tmp_file, 'projects');
                end
            end
            
            % load databases and datasets
            temp = yaml.ReadYaml(fullfile(obj.dir_project, 'metainfo.yaml'));
            if isempty(obj.datasets_list)
                obj.datasets_list = temp.datasets_list;
            else
                obj.datasets_list = union(obj.datasets_list, temp.datasets_list);
            end
            if isempty(obj.databases_list)
                obj.databases_list = temp.databases_list;
            else
                obj.databases_list = union(obj.databases_list, temp.databases_list);
            end
        end
        
        %% delete a project
        function del_project(obj)
            tmp_file = fullfile(fi.locate('ease', true), '.projects.mat');
            load(tmp_file, 'projects');
            if isempty(projects)
                return;
            else
                fprintf('there are following projects: \n');
                for m=1:length(projects)
                    fprintf('--%d: %s\n', m, projects{m});
                end
                tmp_k = input('\ndelete which one: ');
                try
                    projects(tmp_k) = [];
                    save(tmp_file, 'projects', '-append');
                catch
                    return;
                end
            end
        end
        
        %% load and write configurations
        read_config(obj, path_file, quiet);
        write_config(obj, path_file);
        
        %% select data
        data_name = select_data(obj, datasets);
        
        %% how-tos
        how_tos(obj);
        
        %% load 2p stack data
        load_stack(obj);
        
        %% load EM data
        load_em(obj);
        
        %% load functional video data
        load_video(obj);
        
        %% load the video data
        Y = load_Y(obj, mscan, mblock);
        
        %% create data loader
        [dl_Yr, dl_Yd] = create_dataloader(obj, scan_id, block_id, T);
        
        %% choose FOV
        choose_FOV(obj, show_fov);
        
        %% align the video data and the 2p stack data
        align_video_stack(obj);
        
        %% rough registration of the video data and stack data.
        [video_shifts, z_vals] = rough_registration_video(obj)
        
        %% extract EM volumes given the scan ID
        Aem_proj = extract_em_segments(obj, mscan, mslice, ssub);
        
        %% get EM masks in a specified scan ID
        Aem = get_Aem_scan(obj, mscan, ssub);
        
        %% new version of loading Aem
        [Aem, segment_ids] = load_Aem(obj, mscan);
        
        %% create/load neuron objects for all video data
        neuron = get_MF3D(obj, create_new);
        
        %% get all calium imaging data
        summary_images = calculate_summary_images(obj, mscan, mblock)
        
        %% create masks for EM volumes
        get_em_masks(obj);
        
        %% save the current object
        function save(obj)
            tmp_file = fullfile(obj.output_folder, 'obj.mat');
            save(tmp_file, 'obj');
        end
        
        %% voxelize EM meshes
        voxelize_em(obj, use_parallel);
        
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
        
        %% get transformation
        [A_convert, offset] = get_transformation(obj);
        
        %% set options for voxelization
        options = set_voxelization_options(obj, voxel_em);
        
        %% set options for projecting EM
        options = set_projections_options(obj);
        
        %% connect to database
        connect_database(obj, database_list, dj_username, dj_password)
        
        %% visualize the mesh of one EM segment
        visualize_em_mesh(obj, em_id, new_figure, mesh_color, z_zoomin);
        
        %% visualize the voxelized version of one EM segment
        visualize_em_voxels(obj, em_id, new_figure, voxel_color);
        
        %% collect EM information
        collect_em_info(obj);
        
        %% visualize EM footprints
        visualize_em_footprints(obj, em_id, new_figure)
        
        %% project EM
        project_em(obj);
        
        %% add a new dataset
        add_dataset(obj);
        
        %% convert indices from 2p stack to video
        [idx_new, dims_new] = convert_idx(obj, idx);
        
        %% extract EM footprints given a list of EM
        [Aem, segments_id, segments_del] = get_em_footprints(obj, em_ids, scan_id);
        
        %% extract video data
        Y = get_Y(obj, mscan, mblock);
        
        %% get the EM volume on the scanning planes
        spatial_range = get_spatial_range(obj);
        
        %% construct EM footprints on the calcium imaging planes
        construct_Aem(obj);
        
        %% add EM components to the white/black list
        add2whitelist(obj, segment_ids, reason)
        add2blacklist(obj, segment_ids, reason)
        
        %% construct a data container
        [Y_raw, Y_denoised] = construct_Y(obj);
        
        %% release space for storing the raw video
        freespace_Y(obj, mscan, mblock);
        
        %% close all windows except the main window
        function closeall(obj)
            %%
            set(obj.gui.fig_main, 'HandleVisibility', 'off');
            close all;
            set(obj.gui.fig_main, 'HandleVisibility', 'on');
        end
        
        %% create configurations for running EASE
        function configs = create_running_configs(obj, K_new, K_candidates, nb, min_rank, black_list, white_list)
            options = obj.options_init;
            niter = length(nb);
            configs = cell(niter, 1);
            if ~exist('black_list', 'var')
                black_list = []; 
            end 
            if ~exist('white_list', 'var') 
                white_list = []; 
            end
            % configurations for running EASE
            for m=1:niter
                options.K_candidate = K_candidates(m);
                options.K_new = K_new(m);
                if m>1
                    options.clear_results = false;
                end
                configs{m} = struct('options_init', options,...
                    'min_rank', min_rank(m), 'nb', nb(m), ...
                    'black_list', black_list, 'white_list', white_list);
            end
        end
        %% load em volume
        function get_em_volume(obj)
            
            file_name = fullfile(obj.output_folder, ...
                sprintf('segmentation_%d_zblur_%d', ...
                obj.em_segmentation, obj.em_zblur),  ...
                'em_volume.mat');
            if exist(file_name, 'file')
                temp = load(file_name);
                obj.em_volume = temp.em_volume;
            else
                obj.construct_Aem();
            end
        end
        
        %% load conditions 
        function conditions = load_stimuli(obj)
            file_stimuli = fullfile(obj.data_folder, 'conditions.mat');
            temp = load(file_stimuli);
            
            conditions = 2*pi/360*temp.conditions;
            obj.conditions = conditions; 
%             xbins = unique(conditions(~isnan(conditions)));
%             sig = pi/10;
        end 
        
        %% compute tuning curves 
        function compute_tuning_curves(obj, scan_list, shuffle)
            if ~exist('shuffle', 'var')
                shuffle = []; 
            end 
            n_list = length(scan_list); 
            if isempty(obj.conditions)
                obj.load_stimuli(); 
            end 
            for m=1:n_list 
                obj.scan_id = scan_list(m); 
                neuron = obj.get_MF3D(); 
                tmp_stimuli = obj.conditions(obj.scan_id, neuron.frame_range(1):neuron.frame_range(2)); 
                neuron.compute_tuning_curve(tmp_stimuli); 
                neuron.boostrap_tuning_curve(shuffle);
                obj.save_MF3D();
            end 
        end 
        %% process all scans given the configuration 
        function process_scans(obj, scan_list, configs, frame_range)
            n_list = length(scan_list); 
            for m=1:n_list
               obj.scan_id = scan_list(m); 
               neuron = obj.get_MF3D();
               if exist('frame_range', 'var') && ~isempty(frame_range)
                   neuron.frame_range = frame_range;
               end
               neuron.at_ease(configs); 
               obj.save_MF3D(); 
            end
        end 
        %% combine results
        function [results, neurons] = combine_results(obj, scan_list)
            if ~exist('scan_list', 'var') || isempty(scan_list)
                scan_list = 1:obj.num_scans;
            end
            n_list = length(scan_list);
            neurons = cell(n_list,1);
            em_ids = cell(n_list,1);
            snrs = cell(n_list, 1);
            pnrs = cell(n_list, 1);
            confidences = cell(n_list, 1);
            match_scores = cell(n_list, 1);
            labels_all = cell(n_list, 1);
            A_all = cell(n_list, 1);
            Acorr_all = cell(n_list, 1);
            Aem_all = cell(n_list, 1);
            C_all = cell(n_list, 1);
            Craw_all = cell(n_list, 1);
            S_all = cell(n_list, 1);
            tc_y = cell(n_list, 1); 
            tc_yfit = cell(n_list, 1); 
            tc_rss = cell(n_list, 1); 
            tc_pvals = cell(n_list, 1); 
            tc_pars = cell(n_list, 1); 
            % load all results
            for m=1:n_list
                obj.scan_id = scan_list(m);
                neuron = obj.get_MF3D(false);  % load results; don't create new
                neuron.normalize('c_noise'); 
                neurons{m} = neuron;
                
                em_ids{m} = cell2mat(neuron.match_status.em_ids)';
                confidences{m} = (neuron.match_status.confidence);
                match_scores{m} = (neuron.match_status.scores);
                labels_all{m} = neuron.labels;
                
                A_all{m} = neuron.A;
                Acorr_all{m} = neuron.A_corr;
                Aem_all{m} = neuron.A_em;
                C_all{m} = neuron.C;
                Craw_all{m} = neuron.C_raw;
                S_all{m} = neuron.S;
                snrs{m}= var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);
                pnrs{m} = max(neuron.C, [], 2) ./ std(neuron.C_raw-neuron.C, 0, 2);
                
                % tuning curves 
                tc = neuron.tuning_curve; 
                if ~isempty(tc)
                    tc_bins = tc.x; 
                    tc_y{m} = tc.y; 
                    tc_yfit{m} = tc.yfit; 
                    tc_pars{m} = tc.pars; 
                    tc_rss{m} = tc.rss; 
                    tc_pvals{m} = tc.pvals;
                    npars = size(tc.pars, 1); 
                end
            end
            
            % organize results by neurons
            unique_ids = unique(cell2mat(em_ids));
            K = length(unique_ids);
            d = size(A_all{1}, 1);
            T = size(C_all{1}, 2);
            
            results.em_ids = unique_ids;
            results.scan_list = scan_list;
            results.detected = false(n_list, K);
            results.confidences = zeros(n_list, K);
            results.match_scores = zeros(n_list, K);
            results.labels = zeros(n_list, K);
            results.snrs = zeros(n_list, K);
            results.pnrs = zeros(n_list, K);
            results.A = zeros(d, n_list, K);
            results.A_corr = zeros(d, n_list, K);
            results.A_em = zeros(d, n_list, K);
            results.C = zeros(T,n_list, K);
            results.Craw = zeros(T, n_list, K);
            results.S = zeros(T, n_list, K);
            if exist('tc_bins', 'var')
                nbins = length(tc_bins);
                results.bins = tc_bins; 
                results.y = zeros(nbins, n_list, K); 
                results.yfit = zeros(nbins, n_list, K); 
                results.pars = zeros(npars, n_list, K); 
                results.rss = zeros(n_list, K); 
                results.pvals = zeros(n_list, K); 
            end 
            for m=1:n_list
                [~, idx] = ismember(em_ids{m}, unique_ids);
                results.detected(m, idx) = true;
                results.confidences(m, idx) = confidences{m};
                results.match_scores(m, idx) = match_scores{m};
                results.labels(m, idx) = labels_all{m};
                results.snrs(m, idx) = snrs{m};
                results.pnrs(m, idx) = pnrs{m};
                results.A(:, m, idx) = A_all{m};
                results.A_corr(:, m, idx) = Acorr_all{m};
                results.A_em(:, m, idx) = Aem_all{m};
                results.C(:, m, idx) = C_all{m}';
                results.Craw(:, m, idx) = Craw_all{m}';
                results.S(:, m, idx) = S_all{m}';
                
                if exist('tc_bins', 'var')
                    results.y(:, m, idx) = tc_y{m}; 
                    results.yfit(:, m, idx) = tc_yfit{m}; 
                    results.pars(:, m, idx) = tc_pars{m}; 
                    results.rss(m, idx) = tc_rss{m}; 
                    results.pvals(m, idx) = tc_pvals{m}; 
                end 
            end
            
            % compress results 
            dims = [obj.d1, obj.d2, obj.d3, n_list]; 
            results.dims = dims; 
            d = prod(dims); 
            results.A = sparse(reshape(results.A, d, K)); 
            results.A_corr = sparse(reshape(results.A_corr, d, K)); 
            results.A_em = sparse(reshape(results.A_em, d, K)); 
            
        end
        
        %% backup results 
        function backup_results(obj)
            %% determine the matfile for saving the results
            FOV_ = obj.FOV;
            fov_str = sprintf('neurons_%d_%d_%d_%d', FOV_(1), FOV_(2), FOV_(3), FOV_(4)); 
            matfile_mf3d = fullfile(obj.output_folder, sprintf('%s.mat',fov_str));
            folder_MF3D = fullfile(obj.output_folder, fov_str);
            
            %% 
            tmp_str = input('backup name: ', 's'); 
            tmp_folder = fullfile(obj.output_folder, tmp_str); 
            if ~exist(tmp_folder, 'dir')
                mkdir(tmp_folder); 
                copyfile(matfile_mf3d, tmp_folder); 
                copyfile(folder_MF3D, fullfile(tmp_folder, fov_str)); 
                            fprintf('The results of the current sources extraction were backed up to the folder\n\t%s\n', tmp_folder);
            else
                fprintf('the backup name has been used.\n'); 
            end
           
            
        end
    end
end

