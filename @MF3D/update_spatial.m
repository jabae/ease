function update_spatial(obj, Y, with_EM, preprocess_Y)
%% udpate spatial components

%% inputs:
%   Y: d X T matrix, data
%   num: scalar. If method=='hals', then num is the number of iterations to
%       update A; If method=='nnls', then num is the maximum number of neurons
%       overlapping at one pixel
%   method: method for updating the spatial components {'hals', 'nnls'}.
%       default: 'nnls'

%% Author: Pengcheng Zhou, Carnegie Mellon University.
%% code 

if ~exist('preprocess_Y', 'var') || isempty(preprocess_Y)
    preprocess_Y = true; 
end 

%% pre-process Y 
Y = obj.reshape(Y,1);   % reshape the data to a matrix

if preprocess_Y
    % select frames to be analyzed
    if isempty(obj.frame_range)
        Y = double(Y);
    else
        t0 = obj.frame_range(1);
        t1 = obj.frame_range(2);
        Y = double(Y(:, t0:t1));
    end
    
    % normalize data
    if obj.options.normalize_data
        sn = obj.reshape(obj.P.sn, 1);
        Y = bsxfun(@times, Y, 1./sn);
    end
    
    % remove all pixels outside of the EM volume
    if ~isempty(obj.spatial_range)
        Y(~obj.spatial_range, :) = 0; % remove pixels outside of the EM volume
    end
end
if ~exist('with_EM', 'var') || isempty(with_EM)
    with_EM = true; 
end

%% get the spatial range
Ysignal = Y - obj.reconstruct_background();
d3 = obj.options.d3; 
if with_EM
    masks = zeros(size(obj.A_mask));
    for m=1:size(masks, 2)
        ai = obj.reshape(obj.A_mask(:, m), 3);
        ai_mask = imdilate(repmat(sum(ai, 3)>0, [1, 1, d3]), strel('square', 3));
        %     ai_mask = imdilate(ai>0, strel('square', 3));
        masks(:, m) = ai_mask(:);
    end
else
    masks = zeros(size(obj.A));
    for m=1:size(masks, 2)
        ai = obj.reshape(obj.A(:, m), 3);
        ai_mask = imdilate(ai>0, strel('square', 3));
        %     ai_mask = imdilate(ai>0, strel('square', 3));
        masks(:, m) = ai_mask(:);
    end
end


if strcmpi(obj.options.spatial_algorithm, 'hals_thresh')
    sn = std(Ysignal-obj.A*obj.C, 0, 2); 
    obj.A = HALS_spatial_thresh(Ysignal, obj.A, obj.C, masks, 20, sn); 
else
    obj.A = HALS_spatial(Ysignal, obj.A, obj.C, masks, 20);
end 
obj.b0 = mean(Y, 2) - obj.A*mean(obj.C, 2);

%% post process spatial shapes 
obj.post_process_spatial(); 
end
