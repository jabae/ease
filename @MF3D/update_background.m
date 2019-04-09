function update_background(obj, Y, preprocess_Y)
%% update the background term of the model 
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

%% fit the model 
if strcmpi(obj.options.background_model, 'svd')
    nb = obj.options.nb;
    [obj.b, obj.f, obj.b0] = fit_svd_model(Y, nb, obj.A, obj.C);
end
end