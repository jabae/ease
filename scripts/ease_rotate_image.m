[d1, d2] = size(img_to_be_rotated); 
rot_info = struct('rot_angle', 0, 'rot_xlim', 0, 'rot_ylim', 0); 
close; 
h_rotate = figure(); 
ax_img = axes('parent', h_rotate, 'position', [0.3, 0.15, 0.68, 0.8], 'unit', ...
    'normalized'); 
imagesc(img_to_be_rotated, 'parent', ax_img); 
axis tight equal; 

pos_0 = [0.05, 0.8, 0.15, 0.09]; 
text_angle = uicontrol('parent', h_rotate, 'position', pos_0, 'unit', ...
    'normalized', 'style', 'pushbutton', 'string', 'rot. angle', 'fontweight', ...
    'bold'); 

%% to be done 