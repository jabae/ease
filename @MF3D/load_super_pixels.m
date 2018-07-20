function load_super_pixels(obj, sp)
%% load results of super-pixel method

A_sp = sp.fin_rlt.a;
C_sp = sp.fin_rlt.c_tf';
C_raw_sp = sp.fin_rlt.c';
if isempty(C_sp)
    C_sp = C_raw_sp;
end
K_sp = size(A_sp, 2);
b0 = sp.fin_rlt.b;
fb = sp.fin_rlt.fb;
ff = sp.fin_rlt.ff';

% neuron id
ids_sp = (1:size(A_sp, 2));
if ~isempty(obj.A)
    ids_sp = ids_sp + max(obj.ids);
end

% matching status

status_sp = [obj.match_status.status, zeros(1, K_sp)];
em_ids = [obj.match_status.em_ids, cell(1, K_sp)];

obj.A = [obj.A, A_sp];
obj.C_raw = [obj.C_raw; C_raw_sp];
obj.C =  [obj.C; C_sp];
obj.b0 = b0;
obj.f = ff;
obj.b = fb;

obj.ids = [obj.ids, ids_sp];
obj.match_status.status = status_sp;
obj.match_status.em_ids = em_ids;
