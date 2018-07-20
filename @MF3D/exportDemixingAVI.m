function exportDemixingAVI(obj, Y, min_max, avi_nm, col_map)
% export matrix data to movie
%min_max: 1*2 vector, scale
%avi_nm: string, file name
%col_map: colormap

T = size(Y, ndims(Y));
Y = Y(:);
if ~exist('col_map', 'var') || isempty(col_map)
    col_map = jet;
end
if ~exist('avi_nm', 'var') || isempty(avi_nm)
    avi_nm = 'a_movie_with_no_name.avi';
end
if ~exist('min_max', 'var') || isempty(min_max)
    min_max = [min(Y(1:10:end)), max(Y(1:10:end))];
end

Y = uint8(64*(Y-min_max(1))/diff(min_max));
Y(Y<1) = 1;
Y(Y>64) = 64;
col_map = uint8(255*col_map);
Y = reshape(col_map(Y, :), obj.options.d1, [], T, 3);
Y = permute(Y, [1,2,4,3]);

avi_file = VideoWriter(avi_nm);
avi_file.open();
for m=1:T
    avi_file.writeVideo(squeeze(Y(:, :, :, m)));
end
avi_file.close();
end