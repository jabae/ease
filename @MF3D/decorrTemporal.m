function C_ = decorrTemporal(obj, wd)
if ~exist('wd', 'var') || isempty(wd)
    wd = 3; 
end
%% find neighbors of each neurons
A_ = obj.A;
S_ = full(obj.S);
neuron_sn = GetSn(obj.C_raw);
pvals = obj.P.kernel_pars;
S_ = bsxfun(@times, S_, 1./neuron_sn);
[~ ,ind] = max(obj.A, [], 1); 
d1 = obj.options.d1; 
d2 = obj.options.d2; 
d3 = obj.options.d3; 
[r, c, p] = ind2sub([d1, d2, d3], ind); 
temp = sqrt(bsxfun(@minus, r, r').^2 + bsxfun(@minus, c,c').^2 + ...
    bsxfun(@minus, p, p').^2); 
ind_neigh = (temp<10); 

[K, T] = size(S_);
C_ = obj.C;

num_per_row = 100;
for m=1:K
    fprintf('|');
    if mod(m, num_per_row)==0
        fprintf('\n');
    end
end
fprintf('\n');


%% get kernels
model = obj.options.deconv_options.type;

%% remove small spikes
Tk = 500; 
for k=1:K
    s = S_(k, :);
    tmpS = S_(ind_neigh(k, :), :);
    if isempty(tmpS)
        continue;
    else
        ind = (s~=max(tmpS, [], 1));
    end
    
    if wd>1
        ind = (conv(double(ind), ones(1, wd), 'same') > 0); 
    end
    s(ind) = 0;
    
    fprintf('.');
    
    % generate calcium trace
    if strcmpi(model, 'ar1')
        ht = (pvals(k)).^(0:(Tk-1));
    elseif strcmpi(model, 'ar2')
        ht = exp2kernel(ar2exp(pvals(k,:)), Tk);
    end
    c = conv(double(full(s')), ht/ht(1))*neuron_sn(k);
    C_(k, :) = c(1:T);
    if mod(k,num_per_row)==0
        fprintf('\n');
    end
end
fprintf('\n');


obj.C = C_;
end
