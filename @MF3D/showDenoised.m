function showDenoised(obj, min_max, col_map, avi_nm, t_pause)
d1 = obj.options.d1;
d2 = obj.options.d2;
d3 = obj.options.d3;
nc = d3;
nr = 2;

rot_angle = 6.5;
rot_xlim = [11, 127];
rot_ylim = [30, 53];
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
    max(img_height*nr,img_width*nc)*8);
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
Y = obj.A * obj.C;
T = size(Y, ndims(Y));
if (nargin<3) || (isempty(min_max))
    d = numel(Y); 
    temp = Y(randi(d, min(1000, d)));
    min_max = quantile(temp(:), [0.2, 0.9999]);
    min_max(1) = max(min_max(1), 0);
    min_max(2) = max(min_max(2), min_max(1)+0.1);
end
if ~exist('t_pause', 'var'); t_pause=0.01; end

% get demixed video
Y_mixed = zeros(obj.options.d1*obj.options.d2*obj.options.d3, ...
    T, 3);
temp = prism;
K = size(obj.A, 2);
col = temp(randi(64, K,1), :);
for m=1:3
    Y_mixed(:, :, m) = obj.A* (diag(col(:,m))*obj.C);
end
Y_mixed = uint16(Y_mixed*2/min_max(2)*65536);

kt = 3;
k_res = 2;
for t=1:kt:T
    yac = obj.reshape(Y(:, t), 3);
    ymixed = obj.reshape(Y_mixed(:, t, :), 3);
    for z=d3:-1:1
        axes(ha(z,1));
        
        if rot_flag
            imagesc(imrotate(yac(:, :, z), rot_angle),min_max); colormap(col_map);
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(yac(:, :, z), min_max); colormap(col_map);
            axis tight;
        end
        set(ha(z, 1), 'xtick', [], 'ytick', []);
        if z==1
            ylabel('denoised (1x)');
        end
        
        
        axes(ha(z, 2));
        if rot_flag
            imagesc(imrotate(squeeze(ymixed(:, :, z, :)), rot_angle));
            xlim(rot_xlim);
            ylim(rot_ylim);
        else
            imagesc(squeeze(ymixed(:, :, z, :)));
            axis tight;
        end
        set(ha(z, 2), 'xtick', [], 'ytick', []);
        
        if z==1
            ylabel('demixed');
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
