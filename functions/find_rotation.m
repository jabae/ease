function [rot_angle, rot_xlim, rot_ylim] = find_rotation(spatial_range)
%% find an optimal rotation to show EM volume
img = double(sum(spatial_range, 3)>0);
[d1, d2] = size(img);
nnz_row = sum(sum(img,2)>0);
nnz_column = sum(sum(img,1)>0);

if nnz_row>nnz_column
    y = (1:d1);
    temp = bsxfun(@times, img, (1:d2));
    temp(temp==0) = nan;
    [v, x] = max(temp, [], 2);
    x(isnan(v)) = [];
    y(isnan(v)) = [];
else
    x = (1:d2);
    temp = bsxfun(@times, img, (1:d1)');
    temp(temp==0) = nan;
    [v, y] = max(temp, [], 1);
    x(isnan(v)) = [];
    y(isnan(v)) = [];
end

% remove the first and the last two pixels for removing outliers
x = x(3:(end-2));
y = y(3:(end-2));
k = polyfit(x, y, 1);
rot_angle = atan(k(1)) /pi *180;

img_new = imrotate(img, rot_angle);
img_new = imerode(img_new, strel('square', 4));
[a, b] = find(img_new);
rot_xlim = [min(b), max(b)];
rot_ylim = [min(a), max(a)];
end