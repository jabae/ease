function showNeuron(obj, ind, orientation)
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
if ~exist('orientation', 'var') || isempty(orientation)
    orientation = 'horizental';
end
ai = obj.reshape(obj.A(:, ind), 3);
img_max = max(ai(:))*0.8;
vlim = [0, img_max];

[d1, d2, d3] = size(ai);
ai_em = obj.reshape(obj.A_mask(:, ind), 3);
img_max = max(ai_em(:))*0.8;
vlim_em = [0, img_max];

ai_corr = obj.reshape(obj.A_corr(:, ind), 3);
img_max = max(ai_corr(:));
vlim_corr = [0, img_max];
%% crop region
ai = ai(4:(end-4), 4:(end-4), :);
ai_em = ai_em(4:(end-4), 4:(end-4), :);
ai_corr = ai_corr(4:(end-4), 4:(end-4), :);
if strcmpi(orientation, 'horizental')
    h = 4*(d1+2) + d1/2;
    w = d3*(d2+2)+d2/5;
    
    pos_ax = [1-(d2+1)/w, 1-(d1+2)/h, d2/w, d1/h];
    dpos_h = [-(d2+2)/w, 0, 0, 0];
    dpos_v = [0, -(d1+2)/h, 0, 0];
else
    % we only use horizental display for now
    %     h = d3*(d1+2);
    %     w = d2+2;
    %     pos_ax = [1/w, 1/h, 1-2/w, d1/h];
    %     dpos = [0, (d1+2)/h, 0, 0];
end
figure('papersize', [w, h] / max(w, h)*8);
set(gcf, 'color', 'w', ...
    'position', [200, 200, 100*get(gcf, 'papersize')], ...
    'paperposition', [0, 0, get(gcf, 'papersize')]);

for m=1:d3
    axes('position', pos_ax);
    imagesc(ai_corr(:, :, d3-m+1), vlim_corr);
    pos_ax = pos_ax + dpos_h;
    if m==3
        ylabel('corr(Y, c_i)');
    end
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box on;
end

pos_ax = pos_ax - 3*dpos_h + dpos_v;

for m=1:d3
    axes('position', pos_ax);
    imagesc(ai(:, :, d3-m+1), vlim);
    pos_ax = pos_ax + dpos_h;
    if m==3
        ylabel('a_i');
    end
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box on;
end
pos_ax = pos_ax - 3*dpos_h + dpos_v;

for m=1:d3
    axes('position', pos_ax);
    imagesc(ai_em(:, :, d3-m+1), vlim_em);
    pos_ax = pos_ax + dpos_h;
    if m==3
        ylabel('p_i');
    end
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    box on;
end

pos_ax = pos_ax - dpos_h + dpos_v;
pos_ax(3) = (3*d2+4)/w;
axes('position', pos_ax);
T = size(obj.C, 2);
Fs = obj.Fs;

plot((1:T)/Fs, obj.C_raw(ind, :), 'b');
hold on;
plot((1:T)/Fs, obj.C(ind, :), 'r', 'linewidth', 2);
axis tight;
temp = get(gca, 'ylim');
set(gca, 'ylim', temp + [0, temp(2)*0.1]);
legend('raw', 'denoised', 'orientation', 'horizental');
xlabel('Time (sec.)');

end