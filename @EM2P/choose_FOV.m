function choose_FOV(obj, show_fov)
%% what does this function do
%{
    choose FOV of CI data to include the whole EM volume
%}

%% inputs:
%{
    show_fov: boolean, show locations of the cropped FOV
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
% show FOV or not
if ~exist('show_fov', 'var')
    show_fov = 'false';
end
ssub = obj.dims_stack(1) / obj.dims_video(1);

% check whether the FOV has been selected.
if ~isempty(obj.FOV)
    % if selected already, skip
    fprintf('The FOV has been chosen. \n Please set ease.FOV = [] if you want to choose FOV again.\n\n ');
else
    
    % load EM info
    %     if isempty(obj.em_info)
    %         obj.load_em();
    %     end
    
    % spatial downsampling between 2p stack and video data
    
    %% crop a rectangle area to include all EM volumes
    em_ranges_all = cell2mat(obj.em_ranges);
    ii0 = max(round(min(em_ranges_all(:,1)/ssub))-obj.extra_margin, 1);
    ii1 = min(round(max(em_ranges_all(:,1)/ssub))+obj.extra_margin, obj.dims_video(1));
    jj0 = max(round(min(em_ranges_all(:,2)/ssub))-obj.extra_margin, 1);
    jj1 = min(round(max(em_ranges_all(:,2)/ssub))+obj.extra_margin, obj.dims_video(2));
    obj.FOV = [ii0, ii1, jj0, jj1];
    obj.FOV_stack = ssub * obj.FOV + [-1, 0, -1, 0];
    obj.d1 = diff(obj.FOV(1:2)) + 1;
    obj.d2 = diff(obj.FOV(3:4)) + 1;
    obj.d3 = obj.num_slices;
end


%% show the selected FOV
if show_fov
    if ~exist('em_ranges_all', 'var')
        em_ranges_all = cell2mat(obj.em_ranges);
    end
    d1_stack = obj.d1 * ssub;
    d2_stack = obj.d2 * ssub;
    FOV_stack_ = obj.FOV_stack;
    
    figure('papersize', [5, 5]);
    init_fig;
    set(gcf, 'defaultaxesfontsize', 15);
    axes('position', [0.01, 0.01, 0.98, 0.98]);
    fill([1, obj.dims_stack(1), obj.dims_stack(1), 1, 1], ...
        [1, 1, obj.dims_stack(2), obj.dims_stack(2), 1], ones(1,3)*0.8);
    hold on;
    k = convhull(em_ranges_all(:,2), em_ranges_all(:,1));
    
    fill([1, d2_stack, d2_stack,1, 1]+FOV_stack_(3)-1, ...
        [1, 1, d1_stack, d1_stack, 1]+FOV_stack_(1), 'r');
    fill(em_ranges_all(k, 2), em_ranges_all(k,1), 'g');
    set(gca, 'ydir', 'reverse');
    legend('whole FOV', 'cropped FOV', 'convex hull of the EM volume');
    axis equal off tight;
    saveas(gcf, fullfile(obj.output_folder, 'FOV.pdf'));
end
end