function [A_convert, offset] = get_transformation(obj)
%% infer the transformation function to map coordinates from EM space to 2P imaging space
%{
    y_2p = A*y_em + offset
%}

%% inputs:
%{
%}

%% outputs:
%{
    A_convert: 3*3 matrix 
    offset:    3*1 vector 
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code

if (~isempty(obj.transformation)) && (~isempty(obj.transformation.A_convert))
    % the transformation has been pre-computed
    A_convert = obj.transformation.A_convert;
    offset = obj.transformation.offset;
elseif exist(fullfile(obj.data_folder, obj.matfile_transformation), 'file')
    % the transformation has been computed but not loaded
    temp = matfile(fullfile(obj.data_folder, obj.matfile_transformation));
    A_convert = temp.A_convert;
    offset = temp.offset;
    obj.transformation = struct('A_convert', A_convert, 'offset', offset);
else
    %% compute the transformationa gain
    % load the coordinates of neuron centroids in two spaces
    temp = csvread(fullfile(obj.data_folder, obj.registration_csv), 1);
    y_2p = temp(:, 4:6);    % locations in the 2p space (unit: um)
    y_em = temp(:, 1:3)/1000;    % locations in the em space (unit: nm)
    ref_ids = temp(:, end);
    
    % compute the transformation: y_2p = y_em * A_convert + offset
    X = [y_em, ones(size(y_em,1),1)];
    tmpA = (X'*X)\(X'*y_2p);
    
    A_convert = tmpA(1:3, :);
    offset = tmpA(4, :);
    
    % save the results
    obj.transformation = struct('A_convert', temp.A_convert, 'offset', temp.offset);
    
    temp_matfile = fullfile(obj.data_folder, obj.matfile_transformation);
    save(temp_matfile, 'A_convert', 'offset', 'ref_ids', 'y_2p', 'y_em');
    fprintf('The transformation matrix was saved to the mat file\n\t%s\n', temp_matfile);
end
