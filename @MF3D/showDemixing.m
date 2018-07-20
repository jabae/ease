function showDemixing(obj, Y, min_max, col_map, avi_nm, t_pause, ind_neuron)
Y = double(Y);
d1 = obj.options.d1;
d2 = obj.options.d2;
d3 = obj.options.d3;
nc = d3;
nr = 6;
if ~exist('ind_neuron', 'var') || isempty(ind_neuron)
    ind_neuron = 1:size(obj.A, 2);
end
rot_angle = 6.5;
rot_xlim = [13, 125];
rot_ylim = [32, 51];
rot_flag = true;
if rot_flag
    img_width = diff(rot_xlim);
    img_height = diff(rot_ylim);
else
    img_width = d2;
    img_height = d1;
end

% play movies
figure('papersize', [img_width*nc, img_height*nr*2]/...
    max(img_height*nr,img_width*nc)*16);
width = round(img_width*nc/max(img_height*nr,img_width*nc)*1500);
height =round(img_height*nr*1.4/max(img_height*nr, img_width*nc)*1500);
set(gcf, 'position', [100, 100, width, height]);
ha = tight_subplot(nr, nc);
ha = reshape(ha, nc, nr);
if ~exist('col_map', 'var') || isempty(col_map)
    col_map = jet;
end
if exist('avi_nm', 'var') && ischar(avi_nm)
    avi_file = VideoWriter(avi_nm);
    if ~isnan(obj.Fs)
        avi_file.FrameRate = obj.Fs;
    end
    avi_file.open();
    avi_flag = true;
else
    avi_flag = false;
end
if ismatrix(Y); Y=obj.reshape(Y, 2); end
T = size(Y, ndims(Y));
if (nargin<3) || (isempty(min_max))
    temp = Y(:, :, randi(T, min(100, T), 1));
    min_max = quantile(temp(:), [0.2, 0.9999]);
    min_max(1) = max(min_max(1), 0);
    min_max(2) = max(min_max(2), min_max(1)+0.1);
end
if ~exist('t_pause', 'var'); t_pause=0.01; end

sort_frames = true;
noise_method = 'std_res'; 
if strcmpi(noise_method, 'std_res')
    Yres = obj.compute_residual(Y);
    sn = std(Yres, 0, 2);
else
    sn = GetSn(obj.reshape(Y, 1));
end
b_ = obj.b./sn;
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
k_res = 2;
min_max_y = min_max;
for mframe=1:kt:T
    t = tt(mframe);
    ybg = obj.reshape(b_ * obj.f(:, t)+b0_, 3);
    yac = obj.reshape(A_(:, ind_neuron) * obj.C(ind_neuron, t), 3);
    ymixed = obj.reshape(Y_mixed(:, t, :), 3);
    for z=d3:-1:1
        axes(ha(z,1));
        if rot_flag
            imagesc(imrotate(Y(:, :, z, t), rot_angle), min_max_y); colormap(col_map);
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(Y(:, :, z, t), min_max_y); colormap(col_map);
            axis tight;
        end
        axis off;
        if z==ceil(d3/2)
            title(sprintf('raw data (1x, [%d, %d])', min_max_y(1), min_max_y(2)), 'color', 'm');
        end
        axes(ha(z, 2));
        if rot_flag
            imagesc(imrotate(ybg(:, :, z), rot_angle), min_max_y); colormap(col_map);
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(ybg(:, :, 3), min_max_y); colormap(col_map);
            axis tight;
        end
        axis off;
        if z==ceil(d3/2)
            title(sprintf('background (1x, [%d, %d])', min_max_y(1), min_max_y(2)), 'color', 'm');
        end
        
        axes(ha(z,3));
        if rot_flag
            imagesc(imrotate(Y(:, :, z, t)-ybg(:, :, z), rot_angle), min_max); colormap(col_map);
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(Y(:, :, z, t)-ybg(:, :, z), min_max); colormap(col_map);
            axis tight;
        end
        axis off;
        if z==ceil(d3/2)
            title(sprintf('background-subtracted (1x, [%d, %d])', min_max(1), min_max(2)), 'color', 'm');
        end
        
        axes(ha(z, 4));
        if rot_flag
            imagesc(imrotate(yac(:, :, z), rot_angle),min_max); colormap(col_map);
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(yac(:, :, z), min_max); colormap(col_map);
            axis tight;
        end
        axis off;
        if z==ceil(d3/2)
            title(sprintf('denoised data (1x, [%d, %d])', min_max(1), min_max(2)), 'color', 'm');
        end
        axes(ha(z, 5));
        if rot_flag
            imagesc(imrotate(Y(:, :, z, t)-ybg(:, :, z)-yac(:, :, z), rot_angle), (min_max-mean(min_max))/k_res); colormap(col_map);
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(Y(:, :, z, t)-ybg(:, :, z)-yac(:, :, z), (min_max-mean(min_max))/k_res); colormap(col_map);
            axis tight;
        end
        axis off;
        if z==ceil(d3/2)
            temp = (min_max - mean(min_max))/k_res;
            title(sprintf('residual (%dx, [%d, %d])', k_res, temp(1), temp(2)), 'color', 'm');
        end
        
        axes(ha(z, 6));
        if rot_flag
            imagesc(imrotate(squeeze(ymixed(:, :, z, :)), rot_angle));
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(squeeze(ymixed(:, :, z, :)));
            axis tight;
        end
        axis off;
        if z==ceil(d3/2)
            title(sprintf('demixed video (1x, [%d, %d])', min_max(1), min_max(2)), 'color', 'm');
            if isnan(obj.Fs)
                text(rot_xlim(1)+1, rot_ylim(1)+1, sprintf('Frame %d', t), 'fontsize', 15, 'color', 'w');
            else
                text(rot_xlim(1)+1, rot_ylim(1)+3, sprintf('Time = %.2f', t/obj.Fs), 'fontsize', 15, 'color', 'w');
            end
        end
    end
    
    pause(t_pause);
    if avi_flag
        temp = getframe(gcf);
        temp.cdata = imresize(temp.cdata, [height, width]);
        avi_file.writeVideo(temp);
    end
    
end
if avi_flag
    avi_file.close();
end
