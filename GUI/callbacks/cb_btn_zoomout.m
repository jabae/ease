%% zoom out 
v_xlim = round(get(gca, 'xlim'));
v_ylim = round(get(gca, 'ylim'));

temp = {'slice', 'corr', 'em'};
for m=1:3
    for n=1:num_slices
        axes(eval(sprintf('ax_%s{%d}', temp{m}, n)));
        axis tight; 
    end
end