function showNeuron(obj, ind, orientation, vlim)
%% show spatial & temporal components of the function neuron, and its EM footprints
%{
%}

%% inputs: 
%{
    ind: integer, neuron ID 
    orientation: {'vertical', 'horizental' (default)}
    vlim: [vmin, vmax], color range 
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

ai = obj.reshape(obj.A(:, ind), 3);

if ~exist('orientation', 'var') || isempty(orientation)
    orientation = 'vertical';
end
[d1, d2, d3] = size(ai);
img_max = max(ai(:))*0.8;
if ~exist('vlim', 'var') || isempty(vlim)
    vlim = [0, img_max];
end

if strcmpi(orientation, 'horizental')
    h = 3*(d1+2);
    w = d3*(d2+2); 
    
    pos_ax = [1-(d2+1)/w, 1-(d1+2)/h, d2/w, d1/h];
    dpos = [-(d2+2)/w, 0, 0, 0];
else
    h = d3*(d1+2);
    w = d2+2;
    pos_ax = [1/w, 1/h, 1-2/w, d1/h];
    dpos = [0, (d1+2)/h, 0, 0];
end
figure('papersize', [w, h] / max(w, h)*7);
set(gcf, 'color', 'w', ...
    'position', [200, 200, 100*get(gcf, 'papersize')], ...
    'paperposition', [0, 0, get(gcf, 'papersize')]);

for m=1:d3
    axes('position', pos_ax);
    imagesc(ai(:, :, d3-m+1), vlim);
    axis off; %equal off tight;
    pos_ax = pos_ax + dpos;
end
end