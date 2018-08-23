function update_background(obj, Y)
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
if strcmpi(obj.options.background_model, 'svd')
    nb = obj.options.nb;
    [obj.b, obj.f, obj.b0] = fit_svd_model(Y, nb, obj.A, obj.C);
end
end