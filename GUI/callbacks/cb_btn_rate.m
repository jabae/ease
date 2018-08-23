for m=0:5 
    set(eval(sprintf('ease.gui.btn_rate_%d', m)), 'backgroundcolor', 0.93*[1, 1,1]); 
end 

set(eval(sprintf('ease.gui.btn_rate_%d', neuron.match_status.confidence(cell_id))), 'backgroundcolor', 'r'); 