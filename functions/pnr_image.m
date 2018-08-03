function pnr = pnr_image(Y)
%% compute the peak-to-noise ratio (PNR) image given a video data  
%{
    (peak-median) / noise 
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

[d1, d2, T] = size(Y); 
Y = reshape(Y, d1*d2, T); 

Ymax = max(Y, [], 2); 
Ymedian = median(Y, 2); 
Ysn = GetSn(Y); 
Ysn(Ysn==0) = inf; 

pnr = reshape((Ymax-Ymedian)./Ysn, d1, d2); 
