function update_background(obj, Yr)
if ~exist('Yr', 'var')
    Yr = obj.dataloader.load_tzrc();
end
Yr = obj.reshape(Yr, 1);
if strcmpi(obj.options.background_model, 'svd')
    nb = obj.options.nb;
    [obj.b, obj.f, obj.b0] = fit_svd_model(Yr, nb, obj.A, obj.C);
end
end