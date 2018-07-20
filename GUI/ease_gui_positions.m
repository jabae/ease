%% GUI objects 
w1_ = 85; 
w2_ = 25; 
w4_ = 50;
w3_ = d2*1.8; 
dw1_ = 3; 
dw2_ = 1; 
dw3_ = 10; 

h1_ = 20; 
dh1_ = 5; 
h2_ = 80; 
h3_ = d1*1.8; 

w0 = 50; 
h0 = 850; 
font_size = 11; 

%% load  data 
pos_btn_load_ca = [w0, h0, w1_, h1_]; 
pos_btn_load_em = [w0, h0-h1_-dh1_, w1_, h1_]; 

%% scan id
pos_text_scan = pos_btn_load_ca + [w1_+dw3_, 0, 0, 0]; 
pos_scan_left = pos_text_scan + [w1_, 0, w2_-w1_, 0]; 
pos_edit_scan = pos_scan_left + [w2_, 0, 0, 0]; 
pos_scan_right = pos_edit_scan + [w2_, 0, 0, 0]; 

%% plane id 
pos_text_plane = pos_scan_right + [w2_+dw3_, 0, w1_-w2_, 0]; 
pos_plane_left = pos_text_plane + [w1_, 0, w2_-w1_, 0]; 
pos_edit_plane = pos_plane_left + [w2_, 0, 0, 0]; 
pos_plane_right = pos_edit_plane + [w2_, 0, 0, 0]; 

%% block id 
pos_text_block = pos_plane_right + [w2_+dw3_, 0, w1_-w2_, 0]; 
pos_block_left = pos_text_block + [w1_, 0, w2_-w1_, 0]; 
pos_edit_block = pos_block_left + [w2_, 0, 0, 0]; 
pos_block_right = pos_edit_block + [w2_, 0, 0, 0]; 


%% load em 
pos_text_blur = pos_btn_load_em + [w1_+dw3_, 0, 0, 0]; 
pos_blur_left = pos_text_blur + [w1_, 0, w2_-w1_, 0]; 
pos_edit_blur = pos_blur_left + [w2_, 0, 0, 0]; 
pos_blur_right = pos_edit_blur + [w2_, 0, 0, 0]; 

%% run motion correction 


%% summary statistics 
pos_btn_pick = pos_btn_load_em + [0, -2*(h1_+dh1_), 0, 0]; 
pos_btn_init_seed = pos_btn_pick + [w1_+dw1_, 0, 0, 0];
pos_btn_match_em = pos_btn_init_seed + [w1_+dw1_, 0, 0, 0]; 
pos_btn_cn = pos_btn_match_em + [w1_+dw1_, 0, w4_-w1_, 0]; 
pos_btn_pnr = pos_btn_cn + [w4_+dw1_, 0, 0, 0]; 
pos_btn_max = pos_btn_pnr + [w4_+dw1_, 0, 0, 0]; 
pos_btn_mean = pos_btn_max + [w4_+dw1_, 0, 0, 0]; 
pos_btn_sn = pos_btn_mean + [w4_+dw1_, 0, 0, 0]; 

%% zoom-in & zoom-out 
pos_btn_zoomin = pos_btn_sn + [w4_+dw1_, 0, w1_-w4_, 0]; 
pos_btn_zoomout = pos_btn_zoomin + [w1_+dw1_, 0, 0, 0];

%% seed method 
pos_btn_cnmf = pos_btn_load_em + [0, -h1_-dh1_, 0,0]; 
pos_btn_auto_init = pos_btn_cnmf + [w1_+dw1_, 0, 0,0];
pos_btn_manual_init = pos_btn_auto_init + [w1_+dw1_, 0, 0,0];
pos_btn_load_init = pos_btn_manual_init + [w1_+dw1_, 0, 0,0];

pos_btn_ca_ahead = [pos_btn_load_init(1)+w1_, pos_btn_load_init(2), w2_, h1_];
pos_text_ca_current = pos_btn_ca_ahead + [w2_+dw1_, 0, w2_, 0];
pos_btn_ca_next = pos_text_ca_current + [2*w2_+dw1_, 0, -w2_, 0]; 
pos_accept_em= pos_btn_ca_next + [w2_+2*dw1_, 0, 0 , 0]; 

%% axes for showing spatial components  
pos_ax_slice = cell(num_slices,1); 
pos_ax_slice{1} = [w0, h0-4*h1_-h3_, w3_, h3_];
pos_ax_corr = cell(num_slices, 1); 
pos_ax_corr{1} = [w0+w3_+dw1_, h0-4*h1_-h3_, w3_, h3_];
pos_ax_em = cell(num_slices, 1); 
pos_ax_em{1} = [w0+2*(w3_+dw1_), h0-4*h1_-h3_, w3_, h3_];
for m=2:num_slices
    pos_ax_slice{m} = pos_ax_slice{m-1} - [0, h3_+dh1_, 0, 0]; 
    pos_ax_corr{m} = pos_ax_corr{m-1} - [0, h3_+dh1_, 0, 0]; 
    pos_ax_em{m} = pos_ax_em{m-1} - [0, h3_+dh1_, 0, 0];    
end

%% axes for showing activity 
tmp_x0 = pos_ax_slice{1}(1); 
tmp_dx = pos_ax_em{1}(1)+pos_ax_em{1}(3) - tmp_x0; 
pos_ax_activity = [tmp_x0, pos_ax_em{3}(2)-h2_-dh1_*4, tmp_dx, h2_]; 
    
%% axes for showing the match score 
pos_ax_score = pos_ax_activity + [0, -h2_-40, 0, 0]; 
pos_btn_em_ahead = [pos_ax_score(1), pos_ax_score(2)-2*h1_, w2_, h1_];
pos_text_em_current = pos_btn_em_ahead + [w2_+dw1_, 0, 0, 0];
pos_btn_em_next = pos_text_em_current + [w2_+dw1_, 0, 0, 0]; 

pos_btn_em_perfect = pos_btn_em_next + [w2_+dw1_*2, 0, w1_-w2_, 0]; 
pos_btn_em_candidate = pos_btn_em_perfect + [w1_*1.02, 0, 0, 0];
pos_btn_em_ignore = pos_btn_em_candidate + [w1_*1.02, 0, 0, 0];
pos_btn_em_zero = pos_btn_em_ignore + [w1_*1.02, 0, 0, 0]; 

pos_check_em= pos_btn_em_zero + [w2_+2*dw1_, 0, 0 , 0]; 
pos_btn_em_only = pos_check_em + [w2_*0.8, 0, w1_-w2_, 0]; 

%% the main figure
x0_fig = 50; 
y0_fig = 100; 
w_fig = max(pos_btn_load_em(1)+w1_, pos_ax_em{1}(1)+w3_)+3; 
h_fig = 900; 
pos_fig = [x0_fig, y0_fig, w_fig, h_fig]; 

%% figure for aligning em data and 2p stack data 
tmp_h = 3*(h3_+h1_); 
pos_fig_align = [pos_fig(1)+w_fig+dw1_, y0_fig+h_fig-tmp_h, w3_+w1_, tmp_h]; 


pos_ax_align_merge = [dw1_, dh1_, w3_, h3_]; 

pos_ax_align_em = pos_ax_align_merge + [0, h3_+h1_, 0, 0]; 
pos_btn_align_em_dn = pos_ax_align_em + [w3_+dw1_, 0, w4_-w3_, floor(h3_/3)-h3_];
pos_text_align_em = pos_ax_align_em + [w3_+dw1_, h3_/3, w4_-w3_, floor(h3_/3)-h3_];
pos_btn_align_em_up = pos_ax_align_em + [w3_+dw1_, h3_*2/3, w4_-w3_, floor(h3_/3)-h3_];

pos_ax_align_2p = pos_ax_align_em + [0, h3_+h1_, 0, 0]; 
pos_btn_align_2p_dn = pos_ax_align_2p + [w3_+dw1_, 0, w4_-w3_, floor(h3_/3)-h3_];
pos_text_align_2p = pos_ax_align_2p + [w3_+dw1_, h3_/3, w4_-w3_, floor(h3_/3)-h3_];
pos_btn_align_2p_up = pos_ax_align_2p + [w3_+dw1_, h3_*2/3, w4_-w3_, floor(h3_/3)-h3_];










