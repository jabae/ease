function [ai, ci, si, ci_raw] = initialize_one(obj, ai_em, Y)
%% initialize one neuron given its EM presentation pi
%{%}

%% inputs:
%{
    ai_em: d_em * 1 matrix, spatial footprints given by EM segments
    Y: d*T matrix or struct data. if Y is a struct data, it has following
    fields: 
        1. ind_in, boolean array, spatial support of the EM component 
        2. ind_out, boolean array, pixels surrounding the EM support 
        3. Y_in, din * T, video within the spatial support 
        4. Y_out, dout * T, video surrounding the EM support 
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

%% start initialization

ai_em = obj.reshape(ai_em, 3);
ai_mask = imdilate(ai_em>0, strel('square', 4)); %imdilate(repmat(sum(ai, 3)>0, [1, 1, 3]), strel('square', 3));

if isstruct(Y)
    ind_in = Y.ind_in; 
    ind_out = Y.ind_out; 
    Yin = double(Y.Yin); 
    Yout = double(Y.Yout);
else
    % choose pixels within the EM masks and pixels surrounding the EM masks
    [ind_in, ind_out] = obj.construct_in_out(ai_em); 
    
    % select video data in & out of the selected EM mask
    Yin = double(Y(ind_in, :));
    Yout = double(Y(ind_out, :));
end 

if ~isempty(obj.P.sn)  % normalize data
    Yin = bsxfun(@times, Yin, obj.P.sn(ind_in,1));
    Yout = bsxfun(@times, Yout, obj.P.sn(ind_out,1));
end

if ~isempty(obj.b)  % subtract background 
    Yin = Yin - obj.b(ind_in, :)*obj.f;
    Yout = Yout - obj.b(ind_out, :)*obj.f;
end 
if ~isempty(obj.A) % subtract neural sigwnal 
    Yin = Yin - obj.A(ind_in,:)*obj.C;
    Yout = Yout - obj.A(ind_out,:)*obj.C;
end

% subtract mean 
Yin = bsxfun(@minus, Yin, mean(Yin,2)); 
Yout = bsxfun(@minus, Yout, mean(Yout, 2)); 
wi = reshape(ai_em(ind_in), [], 1); 
ci0 = wi' * Yin/(wi'*wi); 
ai = double(ind_in);

% initialize (ai, ci)
rank_out = ceil(sum(ind_out)/100); 
[ai(ind_in), ci_raw] = initialize_ac_tf(Yin, Yout, wi, ci0, rank_out);
[~, sn] = estimate_baseline_noise(ci_raw);

% deconvolve temporal traces 
[ci, si, tmp_options] = deconvolveCa(ci_raw, obj.options.deconv_options, 'sn', sn);
ci = reshape(ci, 1,[]);
ci_raw = ci_raw - tmp_options.b; 

% ai_in = Yin * (ci-mean(ci))'; 
% ai_out = Yout * (ci-mean(ci))'; 


