function update_spatial(obj, Y)
%% udpate spatial components

%% inputs:
%   Y: d X T matrix, data
%   num: scalar. If method=='hals', then num is the number of iterations to
%       update A; If method=='nnls', then num is the maximum number of neurons
%       overlapping at one pixel
%   method: method for updating the spatial components {'hals', 'nnls'}.
%       default: 'nnls'

%% Author: Pengcheng Zhou, Carnegie Mellon University.
%% load data
if ~exist('Y', 'var') || isempty(Y)
    if isempty(obj.frame_range)
        Y = evalin('base', 'Y_cnmf');
    else
        temp = obj.frame_range;
        Y = evalin('base', sprintf('Y_cnmf(:, %d:%d)', temp(1), temp(2)));
    end
end
Y = obj.reshape(Y, 1);
%% get the spatial range
ind = obj.spatial_range;
if ~isempty(ind)
    Y(~ind, :) = 0; % remove pixels outside of the EM volume
end
Ysignal = Y - obj.reconstruct_background();
masks = zeros(size(obj.A_mask));
for m=1:size(masks, 2)
    ai = obj.reshape(obj.A_mask(:, m), 3);
    ai_mask = imdilate(repmat(sum(ai, 3)>0, [1, 1, 3]), strel('square', 3));
    %     ai_mask = imdilate(ai>0, strel('square', 3));
    masks(:, m) = ai_mask(:);
end
obj.A = HALS_spatial(Ysignal, obj.A, obj.C, masks, 20);
obj.b0 = mean(Y, 2) - obj.A*mean(obj.C, 2);
end
