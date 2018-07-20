function update_temporal(obj, Y, allow_deletion, weight_em)
%% run HALS to update temporal components
% input:
%   Y:  d*T, fluorescence data
%   allow_deletion: boolean, allow deletion (default: true)
% output:
%   C_raw: K*T, temporal components without being deconvolved

% Author: Pengcheng Zhou, Columbia University

% options
maxIter = obj.options.maxIter;
deconv_options_0 = obj.options.deconv_options;

if ~exist('allow_deletion', 'var') || isempty(allow_deletion)
    allow_deletion = false;
end
if ~exist('weight_em', 'var') || isempty(weight_em)
    weight_em = 'false'; 
end

Y = obj.reshape(Y,1); 
Ymean = mean(Y, 2); 
Y = Y - obj.reconstruct_background(); 

%% initialization
A = obj.A;
if weight_em
    A_em = obj.A_mask;
else
    A_em = ones(size(A));
end
K = size(A, 2);     % number of components
C = obj.C;
C_raw = zeros(size(C));
S = zeros(size(C));
A = full(A);
U = (A.*A_em)'*Y;
V = (A.*A_em)'*A;
aa = diag(V);   % squares of l2 norm all all components
sn =  zeros(1, K);
smin = zeros(1,K);
% kernel = obj.kernel;
kernel_pars = cell(K,1);
%% updating
ind_del = false(K, 1);
for miter=1:maxIter
    for k=1:K
        if  aa(k)==0
            C_raw(k, :) = C_raw(k, :)*0;
            C(k,:) = C(k, :)*0;
            S(k, :) = S(k, :)*0;
            ind_del(k) = true;
            continue;
        end
        if ind_del
            continue;
        end
        temp = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);
        %remove baseline and estimate noise level
%         [b_hist, sn_hist] = estimate_baseline_noise(temp);
        b = mean(temp(temp<median(temp)));
        sn_psd = GetSn(temp);
%         if sn_psd<sn_hist
            tmp_sn = sn_psd;
%         else
%             tmp_sn = sn_hist;
%             b = b_hist;
%         end
        
        temp = temp -b;
        sn(k) = tmp_sn;
        
        % deconvolution
        if obj.options.deconv_flag
            [ck, sk, deconv_options]= deconvolveCa(temp, deconv_options_0, 'maxIter', 2, 'sn', tmp_sn);
            smin(k) = deconv_options.smin;
            kernel_pars{k} = reshape(deconv_options.pars, 1, []);
            temp = temp - deconv_options.b; 
        else
            if maxIter>2
                ck = threshold_ci(temp, 3);
            else
                tmp_c = trend_filter(temp);
                
                ck = tmp_c - median(tmp_c);
                temp = temp -median(tmp_c);
            end
        end
        % save convolution kernels and deconvolution results
        C(k, :) = ck;
        
        if sum(ck(2:end))==0
            ind_del(k) = true;
        end
        % save the spike count in the last iteration
        if miter==maxIter
            if obj.options.deconv_flag
                S(k, :) = sk;
            end
            C_raw(k, :) = temp;
        end
    end
end
obj.A = bsxfun(@times, A, sn);
obj.C = bsxfun(@times, C, 1./sn');
obj.C_raw = bsxfun(@times, C_raw, 1./sn');
obj.S = bsxfun(@times, S, 1./sn');
obj.P.kernel_pars =cell2mat( kernel_pars);
obj.P.smin = smin/sn;
obj.P.sn_neuron = sn;
obj.b0 = Ymean - obj.A * mean(obj.C, 2); 
if allow_deletion
    obj.delete(ind_del);
end
