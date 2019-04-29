%% populate the current axis limits to all panels
v_xlim = round(get(gca, 'xlim'));
v_ylim = round(get(gca, 'ylim'));

temp = {'slice', 'corr', 'em'};
for m=1:3
    for n=1:ease.num_slices
        tmp_h = eval(sprintf('ease.gui.ax_%s{%d}', temp{m}, n));
        set(tmp_h, 'xlim', v_xlim);
        set(tmp_h, 'ylim', v_ylim);
    end
end