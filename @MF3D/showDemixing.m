function showDemixing(obj, Y, min_max, col_map, avi_nm, t_pause, ind_neuron, rot_info, fig_visible)
%% show results of running EASE
%{
    show videos of the following data:
    row 1: raw video
    row 2: background video
    row 3: background-subtracted video
    row 4: denoised video (A*C)
    row 5: residual
    row 6: demixing video (A*C). each neuron has its own pesudocolor
%}

%% inputs:
%{
    Y: video data
    col_mat: colormap to show video
    avi_nm: name of the saved avi file
    t_pause:
    ind_neuron: indices of neurons to be shown
    rot_info: rotation information
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
if ~exist('Y', 'var') || isempty(Y)
    if isempty(obj.frame_range)
        Y = evalin('base', 'Y_cnmf');
    else
        temp = obj.frame_range;
        Y = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
    end
end
Y = double(obj.reshape(Y, 1));
d1 = obj.options.d1;
d2 = obj.options.d2;
d3 = obj.options.d3;
nc = d3;
nr = 6;
if ~exist('ind_neuron', 'var') || isempty(ind_neuron)
    ind_neuron = 1:size(obj.A, 2);
end
if ~exist('rot_info', 'var')
    [rot_angle, rot_xlim, rot_ylim] = find_rotation(obj.reshape(obj.spatial_range, 3));
else
    rot_angle = rot_info.rot_angle;
    rot_xlim = rot_info.rot_xlim;
    rot_ylim = rot_info.rot_ylim;
end
if ~exist('fig_visible', 'var')
    fig_visible = 'on';
end
rot_flag = true;
img_width = diff(rot_xlim);
img_height = diff(rot_ylim);


% play movies
width = round(img_width*nc*3/(0.99-(nc-1)*0.005)); 
height = round(img_height*nr*3/(0.99-(nr-1)*0.01));
figure('papersize', 7*[width, height]/max([width, height]));
% width = round(img_width*nc/max(img_height*nr,img_width*nc)*1000);
% height =round(img_height*nr*1.2/max(img_height*nr, img_width*nc)*1000);

set(gcf, 'position', [100, 100, width, height], 'visible', fig_visible);
ha = tight_subplot(nr, nc, [0.01, 0.005], 0.005, 0.005);
ha = reshape(ha, nc, nr);
if ~exist('col_map', 'var') || isempty(col_map)
    col_map = jet;
end
set(gcf, 'colormap', col_map);
if exist('avi_nm', 'var') && ischar(avi_nm)
    avi_file = VideoWriter(avi_nm);
    if ~isnan(obj.Fs)
        avi_file.FrameRate = obj.Fs;
    end
    avi_file.Quality = 100; 
    avi_file.open();
    avi_flag = true;
else
    avi_flag = false;
end
if ismatrix(Y); Y=obj.reshape(Y, 2); end
T = size(Y, ndims(Y));
if (nargin<3) || (isempty(min_max))
    min_max = [1, 17];
end
if ~exist('t_pause', 'var'); t_pause=0.01; end

sort_frames = false;
noise_method = 'std_res';
if strcmpi(noise_method, 'std_res')
    Yres = obj.compute_residual(Y);
    sn = std(Yres, 0, 2);
else
    sn = GetSn(obj.reshape(Y, 1));
end
b_ = bsxfun(@times, obj.b, 1./sn);
b0_ = obj.b0./sn;
A_ = bsxfun(@times, obj.A, 1./sn);
Y = obj.reshape(Y,1);
Y = bsxfun(@times, Y, 1./sn);
Y = obj.reshape(Y, 3);

% get demixed video
Y_mixed = zeros(obj.options.d1*obj.options.d2*obj.options.d3, ...
    T, 3);
temp = prism;
K = size(obj.A(:, ind_neuron), 2);
col = temp(randi(64, K,1), :);

for m=1:3
    Y_mixed(:, :, m) = A_(:, ind_neuron)* (diag(col(:,m))*obj.C(ind_neuron, :));
end
Y_mixed = uint16(Y_mixed*2/min_max(2)*65536);

if sort_frames
    ind = (sum(A_,2)>0);
    x = sum(Yres(ind, :).^2, 1);
    [~, tt] = sort(x, 'descend');
else
    tt = 1:T;
end
kt = 3;
k_res = 1;
min_max_y = min_max;
fprintf('writing the demixing videos...\nProgress bar:\n');
for mframe=1:100
    fprintf('|'); 
end
fprintf('\n'); 
k0 = 0;
for mframe=1:kt:T
    t = tt(mframe);
    ybg = obj.reshape(b_ * obj.f(:, t)+b0_, 3);
    yac = obj.reshape(A_(:, ind_neuron) * obj.C(ind_neuron, t), 3);
    ymixed = obj.reshape(Y_mixed(:, t, :), 3);
    for z=d3:-1:1
        % show raw video data
        if rot_flag
            imagesc(imrotate(Y(:, :, z, t), rot_angle), 'parent', ha(z,1), ...
                min_max_y);
            colormap(col_map);
        else
            imagesc(Y(:, :, z, t), 'parent', ha(z,1), min_max_y); colormap(col_map);
        end
        set(ha(z,1), 'xlim', rot_xlim, 'xtick', [], 'ylim', rot_ylim, 'ytick', []);
        if z==1
            text(rot_xlim(1), rot_ylim(1)+3, sprintf('raw data ([%d, %d])', min_max_y(1), ...
                min_max_y(2)), 'fontsize', 12, 'fontweight', 'bold', 'color', 'w','parent', ha(z, 1))
        end
        
        % show background
        if rot_flag
            imagesc(imrotate(ybg(:, :, z), rot_angle), 'parent', ha(z,2), ...
                min_max_y);
        else
            imagesc(ybg(:, :, z), 'parent', ha(z,2), min_max_y); %#ok<*UNRCH>
        end
        set(ha(z,2), 'xlim', rot_xlim, 'xtick', [], 'ylim', rot_ylim, 'ytick', []);
        if z==1
           text(rot_xlim(1), rot_ylim(1)+3, ...
                sprintf('background ([%d, %d])', min_max_y(1), ...
                min_max_y(2)), 'fontsize', 12, 'fontweight', 'bold', ...
                'color', 'w', 'parent', ha(z,2));
        end
        
        % background subtracted video
        if rot_flag
            imagesc(imrotate(Y(:, :, z, t)-ybg(:, :, z), rot_angle),...
                'parent', ha(z, 3), min_max);
        else
            imagesc(Y(:, :, z, t)-ybg(:, :, z), 'parent', ha(z,3), ...
                min_max);
        end
        set(ha(z,3), 'xlim', rot_xlim, 'xtick', [], 'ylim', rot_ylim, 'ytick', []);
        
        if z==1
               text(rot_xlim(1), rot_ylim(1)+3, ...
                   sprintf('BG-subtracted ([%d, %d])',...
                   min_max(1), min_max(2)), 'fontsize', 12, 'fontweight', ...
                   'bold', 'color', 'w','parent', ha(z, 3))
        end
        
        % denoised video
        if rot_flag
            imagesc(imrotate(yac(:, :, z), rot_angle),...
                'parent', ha(z, 4), min_max);
        else
            imagesc(yac(:, :, z), 'parent', ha(z,4), min_max);
        end
        set(ha(z,4), 'xlim', rot_xlim, 'xtick', [], 'ylim', rot_ylim, 'ytick', []);
        if z==1
            text(rot_xlim(1), rot_ylim(1)+3, sprintf('denoised ([%d, %d])', min_max(1), ...
                min_max(2)), 'fontsize', 12, 'fontweight', 'bold',...
                'color', 'w','parent', ha(z, 4));
        end
        
        % residual
        if rot_flag
            imagesc(imrotate(Y(:, :, z, t)-ybg(:, :, z)-yac(:, :, z), rot_angle),...
                'parent', ha(z, 5), (min_max-mean(min_max))/k_res);
        else
            imagesc(Y(:, :, z, t)-ybg(:, :, z)-yac(:, :, z), ...
                'parent', ha(z, 5), (min_max-mean(min_max))/k_res);
        end
        set(ha(z,5), 'xlim', rot_xlim, 'xtick', [], 'ylim', rot_ylim, 'ytick', []);
        if z==1
            temp = (min_max - mean(min_max))/k_res;
            text(rot_xlim(1), rot_ylim(1)+3, ...
                sprintf('residual ([%d, %d])', temp(1), temp(2)),...
                'color', 'k', 'fontsize', 12, ...
                'fontweight', 'bold','parent', ha(z, 5));
        end
        
        % demixed video
        if rot_flag
            imagesc(imrotate(squeeze(ymixed(:, :, z, :)), rot_angle), ...
                'parent', ha(z, 6));
        else
            imagesc(squeeze(ymixed(:, :, z, :)), 'parent', ha(z, 6));
        end
        set(ha(z,6), 'xlim', rot_xlim, 'xtick', [], 'ylim', rot_ylim, 'ytick', []);
        
        if z==1
            text(rot_xlim(1), rot_ylim(1)+3,'demixed', 'color', 'w', ...
                'parent', ha(z, 6),'fontsize', 12, 'fontweight', 'bold');
%                 sprintf('demixed ([%d, %d])',...
%                 min_max(1), min_max(2)), 'color', 'w', 'parent', ha(z, 6),...
%                 'fontsize', 12, 'fontweight', 'bold'); 
        end
        if z==3
            if isnan(obj.Fs)
                text(rot_xlim(2)/2+1, rot_ylim(2)-3, sprintf('Frame %d', t),...
                    'fontsize', 12, 'color', 'w','parent', ha(z, 6));
            else
                text(rot_xlim(2)/2+1, rot_ylim(2)-3, sprintf('Time = %.2f sec', t/obj.Fs),...
                    'fontsize', 12, 'color', 'w', 'parent', ha(z, 6));
            end
        end
    end
    
    if ~strcmpi(fig_visible, 'off')
        pause(t_pause);
    end
    if avi_flag
        temp = getframe(gcf);
        temp.cdata = imresize(temp.cdata, [height, width]);
        avi_file.writeVideo(temp);
    end
    k1 = round(100*mframe/T); 
    if k1>k0
        for k=(k0+1):k1
            fprintf('.');
        end
        k0 = k1; 
    end
    
end
if avi_flag
    avi_file.close();
end

end