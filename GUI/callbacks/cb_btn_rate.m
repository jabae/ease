% for m=0:5
%     set(eval(sprintf('ease.gui.btn_rate_%d', m)), 'backgroundcolor', 0.93*[1, 1,1]);
% end
% set(eval(sprintf('ease.gui.btn_rate_%d', round(neuron.match_status.confidence(cell_id)))), 'backgroundcolor', 'r');

temp = neuron.match_status.confidence(cell_id);
if temp==0
    set(ease.gui.text_rate, 'string', '0', 'foregroundcolor', 'red');
else
    set(ease.gui.text_rate, 'string', num2str(round(neuron.match_status.confidence(cell_id))), ...
        'foregroundcolor', 'blue');
end