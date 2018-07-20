function playMovie(obj, Y, min_max, col_map, avi_nm, t_pause)
Y = double(Y);
d1 = obj.options.d1;
d2 = obj.options.d2;
d3 = obj.options.d3;

nc = max(1, round(sqrt(d3*d1/d2)));
nr = d3/nc;
% play movies
figure('papersize', [d2*nc,d1*nr]/max(d1,d2)*5);
width = d2*nc/max(d1,d2)*500;
height =d1*nr/max(d1,d2)*500;
set(gcf, 'position', [100, 100, width, height]);
ha = tight_subplot(nr, nc);
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
[~, ~, T] = size(Y);
if (nargin<3) || (isempty(min_max))
    temp = Y(:, :, randi(T, min(100, T), 1));
    min_max = quantile(temp(:), [0.2, 0.9999]);
    min_max(1) = max(min_max(1), 0);
    min_max(2) = max(min_max(2), min_max(1)+0.1);
end
if ~exist('t_pause', 'var'); t_pause=0.01; end
for t=1:size(Y,4)
    for z=d3:-1:1
        axes(ha(z));
        imagesc(Y(:, :, z, t), min_max); colormap(col_map);
        axis equal; axis off tight;
    end
    if isnan(obj.Fs)
        title(sprintf('Frame %d', t));
    else
        text(1, 10, sprintf('Time = %.2f', t/obj.Fs), 'fontsize', 15, 'color', 'w');
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