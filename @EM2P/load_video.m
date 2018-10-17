function load_video(obj)
%% what does this function do 
%{
    map calcium imaging data to the memory 
%}

%% inputs: 
%{
%}

%% outputs: 
%{
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

%% code 
video_mat_file = fullfile(obj.data_folder, obj.matfile_video); 
% load video data
if isempty(obj.video_data)
    % check the existance of video data
    if ~exist(video_mat_file, 'file')
        % TBD
        ease_distribute_functional_data;
    end
    
    % map data
    obj.video_data = matfile(video_mat_file, 'Writable', false);
    obj.video_loader = obj.video_data.dl_videos;
    obj.video_shifts.ii = obj.video_data.shifts_ca_ii;
    obj.video_shifts.jj = obj.video_data.shifts_ca_jj;
    obj.video_zvals = obj.video_data.z_vals;
end
fprintf('\nThe functional video data have been mapped to memory\n\n');

end