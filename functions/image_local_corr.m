function img_corr = image_local_corr(img1, img2, sz)
%% compute local correlations between two images 

[d1, d2, d3] = size(img1); 
img_corr = zeros(d1, d2, d3); 
if ~exist('sz','var') || isempty(sz)
    sz = 8; 
end 
kernel = ones(sz); 
img12 = img1.*img2; 
npixels = imfilter(ones(d1,d2), kernel); 
img12_sum = imfilter(img12, kernel); 
img11_sum = imfilter(img1.^2, kernel); 
img22_sum = imfilter(img2.^2, kernel); 
img1_sum = imfilter(img1, kernel); 
img2_sum = imfilter(img2, kernel); 
for m=1:d3 
    xy_mean = img12_sum(:, :, m)./npixels; 
    x_mean = img1_sum(:, :, m)./ npixels; 
    xx_mean = img11_sum(:, :, m)./npixels; 
    y_mean = img2_sum(:, :, m)./npixels; 
    yy_mean = img22_sum(:, :, m)./npixels; 
    
    img_corr(:, :, m) = (xy_mean-x_mean.*y_mean) ./ (...
        sqrt(xx_mean-x_mean.^2) .* ...
        sqrt(yy_mean-y_mean.^2)); 
end 
img_corr(isnan(img_corr)) = 0; 