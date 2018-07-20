%% choose seed pixels and try to initialize a neuron
tmpY = Y{scan_id, slice_id, block_id};
while true
    axes(ax_slice{slice_id});
    [c, r] = ginput(1);
    c = round(c);
    r = round(r);
    if r<1 || r>d1 || c<1 || c>d2
        cb_btn_slice; 
        for m=1:3
            axes(ax_slice{m});
            hold on;
            active_seed = plot(c_old, r_old, '+r');
        end
        c_seed = c_old;
        r_seed = r_old;
        slice_seed = slice_id;
        break;
    end
    y = squeeze(tmpY(r, c, :));
    axes(ax_activity); hold off;
    plot(y);
    ylim(quantile(double(y), [0.001, 1]));
    title(ax_activity, 'activity');
    c_old = c;
    r_old = r;
end