switch neuron.labels(cell_id)
    case 1
        temp = 'soma'; 
    case 2 
        temp = 'dendrite'; 
    otherwise 
        temp = 'unknown'; 
end

color_gray = 0.94 * [1, 1, 1]; 
if strcmpi(temp, 'soma')
    set(ease.gui.btn_soma, 'backgroundcolor', 'r'); 
else
        set(ease.gui.btn_soma, 'backgroundcolor', color_gray); 
end

if strcmpi(temp, 'dendrite')
    set(ease.gui.btn_dendrite, 'backgroundcolor', 'r'); 
else
    set(ease.gui.btn_dendrite, 'backgroundcolor', color_gray); 
end