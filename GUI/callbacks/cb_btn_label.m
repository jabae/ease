temp = neuron.labels(cell_id); 
color_gray = 0.94 * [1, 1, 1]; 

if temp<0  % the neuron has been verified  
    set(ease.gui.btn_verified, 'backgroundcolor', 'r'); 
else 
    set(ease.gui.btn_verified, 'backgroundcolor', color_gray); 
end 

switch abs(temp)
    case 1
        temp = 'soma'; 
    case 2 
        temp = 'dendrite'; 
    otherwise 
        temp = 'unknown'; 
end

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