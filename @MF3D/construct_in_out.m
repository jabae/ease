function [ind_in, ind_out] = construct_in_out(obj, ai_em)
%% construct in & out area given one EM segment 
%{%}

%% inputs:
%{
    ai_em: d_em * 1 matrix, spatial footprints given by EM segments
%}

%% outputs:
%{
%}

%% author:
%{
    Pengcheng Zhou
    Columbia University, 2018
    zhoupc1988@gmail.com
%}

%% code 
ai = obj.reshape(ai_em, 3);

% ai_mask = imdilate(ai>0, strel('square', 3)); 
ai_mask = imdilate(repmat(sum(ai, 3)>0, [1, 1, 3]), strel('disk', 3));

% choose pixels within the EM masks and pixels surrounding the EM masks
ind_in = obj.reshape(ai_mask, 3);
ind_out = xor(ind_in, imdilate(ai_mask, strel('square', 8)));
ind_in = reshape(ind_in, [], 1); 
ind_out = reshape(ind_out, [], 1); 
