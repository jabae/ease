%% load coordinates of the registered voxels 
temp = csvread(fullfile(ease.data_folder, ease.registration_csv), 1); 
y_2p = temp(:, 4:6);    % locations in the 2p space (unit: um)
y_em = temp(:, 1:3)/1000;    % locations in the em space (unit: nm)
ref_ids = temp(:, end); 

%% compute the transformation: y_2p = y_em * A_convert + offset

X = [y_em, ones(size(y_em,1),1)]; 
tmpA = (X'*X)\(X'*y_2p); 

A_convert = tmpA(1:3, :); 
offset = tmpA(4, :); 
y_2p_hat = bsxfun(@plus, y_em*A_convert, offset);  

temp_matfile = fullfile(ease.data_folder, ease.matfile_transformation); 
save(temp_matfile, 'A_convert', 'offset', 'ref_ids', 'y_2p', 'y_em'); 
fprintf('The transformation matrix was saved to the mat file\n\t%s\n', temp_matfile); 

%% verify 
figure('papersize', [9, 3]); 
init_fig; 

subplot(131); 
plot(y_2p(:,1), y_2p(:,1), 'r'); hold on; 
plot(y_2p(:,1), y_2p_hat(:,1), '.'); hold on; 
title('x axis'); 

subplot(132);
plot(y_2p(:,2), y_2p(:,2), 'r'); hold on; 
plot(y_2p(:,2), y_2p_hat(:,2), '.'); hold on; 
title('y axis');

subplot(133); hold on; 
plot(y_2p(:,3), y_2p(:,3), 'r'); 
plot(y_2p(:,3), y_2p_hat(:,3), '.');  
title('z axis');