function startGUI(obj)
%% create a GUI for easy visualizing data 

%% GUI positions
% frequently used widths and heights
w1_ = 100;
w2_ = 25;
w4_ = 55;
w5_ = 80; 
w3_ = obj.d2*1.8;
dw1_ = 3;
% dw2_ = 1;
dw3_ = 10;

h1_ = 20;
dh1_ = 5;
h2_ = 80;
h3_ = obj.d1*1.8;

w0 = 50;
h0 = 850;
font_size = 11;

% data loading
pos_btn_load_video = [w0, h0, w1_, h1_];

% summary images 
pos_btn_summary_image = [w0, h0-(h1_+dh1_), w1_, h1_]; 

pos_btn_load_em = [w0, h0-2*(h1_+dh1_), w1_, h1_];

% choose scan id
pos_text_scan = pos_btn_load_video + [w1_+dw1_, 0, 0, 0];
pos_scan_left = pos_text_scan + [w1_, 0, w2_-w1_, 0];
pos_edit_scan = pos_scan_left + [w2_, 0, 0, 0];
pos_scan_right = pos_edit_scan + [w2_, 0, 0, 0];
pos_btn_cn = [w0+w1_+dw1_, h0-(h1_+dh1_), w4_, h1_];
pos_btn_pnr = pos_btn_cn + [w4_+dw1_, 0, 0, 0];
pos_btn_max = pos_btn_pnr + [w4_+dw1_, 0, 0, 0];
pos_btn_mean = pos_btn_max + [w4_+dw1_, 0, 0, 0];
pos_btn_sn = pos_btn_mean + [w4_+dw1_, 0, 0, 0];
% zoom-in & zoom-out
pos_btn_zoomin = pos_btn_sn + [w4_+dw1_, 0, w5_-w4_, 0];
pos_btn_zoomout = pos_btn_zoomin + [w5_+dw1_, 0, 0, 0];


% choose slice id
pos_text_plane = pos_scan_right + [w2_+dw1_, 0, w1_-w2_, 0];
pos_plane_left = pos_text_plane + [w1_, 0, w2_-w1_, 0];
pos_edit_slice = pos_plane_left + [w2_, 0, 0, 0];
pos_plane_right = pos_edit_slice + [w2_, 0, 0, 0];

% choose block id
pos_text_block = pos_plane_right + [w2_+dw3_, 0, w1_-w2_, 0];
pos_block_left = pos_text_block + [w1_, 0, w2_-w1_, 0];
pos_edit_block = pos_block_left + [w2_, 0, 0, 0];
pos_block_right = pos_edit_block + [w2_, 0, 0, 0];


% load em
pos_text_blur = pos_btn_load_em + [w1_+dw1_, 0, 0, 0];
pos_blur_left = pos_text_blur + [w1_, 0, w2_-w1_, 0];
pos_edit_blur = pos_blur_left + [w2_, 0, 0, 0];
pos_blur_right = pos_edit_blur + [w2_, 0, 0, 0];


% summary statistics
pos_btn_pick = pos_btn_load_em + [0, -2*(h1_+dh1_), 0, 0];
pos_btn_init_seed = pos_btn_pick + [w1_+dw1_, 0, 0, 0];
pos_btn_match_em = pos_btn_init_seed + [w1_+dw1_, 0, 0, 0];


% seed method
pos_btn_auto_init = pos_btn_load_em + [0, -h1_-dh1_, 0,0];
pos_btn_init_options = pos_btn_auto_init + [w1_+dw1_, 0, 0, 0]; 
pos_btn_cnmf = pos_btn_init_options + [w1_+dw1_, 0, 0,0];
pos_btn_auto_del = pos_btn_cnmf + [w1_+dw1_, 0, 0,0];
% pos_btn_load_init = pos_btn_auto_del + [w1_+dw1_, 0, 0,0];

pos_btn_visualize_neurons = pos_btn_load_em + [0, -2*(h1_+dh1_), 0, 0]; 
pos_btn_ca_ahead = pos_btn_visualize_neurons + [w1_+dw1_, 0, w2_-w1_, 0]; 
pos_text_cell_id = pos_btn_ca_ahead + [w2_+dw1_, 0, w2_, 0];
pos_btn_ca_next = pos_text_cell_id + [2*w2_+dw1_, 0, -w2_, 0];
pos_btn_confidence = pos_btn_ca_next + [w2_+dw1_, 0, 2*w2_, 0]; 
pos_btn_rate_5 = [pos_btn_confidence(1)+ 3*w2_, pos_btn_confidence(2), w2_, h1_];
pos_btn_rate_4 = pos_btn_rate_5 + [w2_, 0, 0, 0];
pos_btn_rate_3 = pos_btn_rate_4 + [w2_, 0, 0, 0];
pos_btn_rate_2 = pos_btn_rate_3 + [w2_, 0, 0, 0];
pos_btn_rate_1 = pos_btn_rate_2 + [w2_, 0, 0, 0];
pos_btn_rate_0 = pos_btn_rate_1 + [w2_, 0, 0, 0];
pos_btn_delete = pos_btn_rate_0 + [w2_, 0, w5_-w2_, 0]; 
pos_btn_soma = pos_btn_delete + [w5_+dw1_, 0, 0, 0]; 
pos_btn_dendrite = pos_btn_soma + [w5_, 0, 0, 0]; 

% axes for showing spatial components
pos_ax_slice = cell(obj.num_slices,1);
pos_ax_slice{1} = [w0, h0-6*h1_-h3_, w3_, h3_];
pos_ax_corr = cell(obj.num_slices, 1);
pos_ax_corr{1} = [w0+w3_+dw1_, h0-6*h1_-h3_, w3_, h3_];
pos_ax_em = cell(obj.num_slices, 1);
pos_ax_em{1} = [w0+2*(w3_+dw1_), h0-6*h1_-h3_, w3_, h3_];
for m=2:obj.num_slices
    pos_ax_slice{m} = pos_ax_slice{m-1} - [0, h3_+dh1_, 0, 0];
    pos_ax_corr{m} = pos_ax_corr{m-1} - [0, h3_+dh1_, 0, 0];
    pos_ax_em{m} = pos_ax_em{m-1} - [0, h3_+dh1_, 0, 0];
end

% axes for showing activity
tmp_x0 = pos_ax_slice{1}(1);
tmp_dx = pos_ax_em{1}(1)+pos_ax_em{1}(3) - tmp_x0;
pos_ax_activity = [tmp_x0, pos_ax_em{3}(2)-h2_-dh1_*4, tmp_dx, h2_];

% axes for showing the match score
pos_ax_score = pos_ax_activity + [0, -h2_-40, 0, 0];
pos_btn_em_ahead = [pos_ax_score(1), pos_ax_score(2)-2*h1_, w2_, h1_];
pos_text_em_current = pos_btn_em_ahead + [w2_+dw1_, 0, 0, 0];
pos_btn_em_next = pos_text_em_current + [w2_+dw1_, 0, 0, 0];

pos_btn_em_perfect = pos_btn_em_next + [w2_+dw1_*2, 0, w1_-w2_, 0];
pos_btn_em_candidate = pos_btn_em_perfect + [w1_*1.02, 0, 0, 0];
pos_btn_em_ignore = pos_btn_em_candidate + [w1_*1.02, 0, 0, 0];
pos_btn_em_zero = pos_btn_em_ignore + [w1_*1.02, 0, 0, 0];

pos_check_em= pos_btn_em_zero + [w1_*1.04, 0, w2_-w1_ , 0];
pos_btn_em_only = pos_check_em + [w2_*0.7, 0, w1_-w2_, 0];

% the main figure
x0_fig = 50;
y0_fig = 100;
w_fig = max(pos_btn_load_em(1)+w1_, pos_ax_em{1}(1)+w3_)+3;
h_fig = 900;
pos_fig = [x0_fig, y0_fig, w_fig, h_fig];

% figure for aligning em data and 2p stack data
% tmp_h = 3*(h3_+h1_);
% pos_fig_align = [pos_fig(1)+w_fig+dw1_, y0_fig+h_fig-tmp_h, w3_+w1_, tmp_h];


% pos_ax_align_merge = [dw1_, dh1_, w3_, h3_];

% pos_ax_align_em = pos_ax_align_merge + [0, h3_+h1_, 0, 0];
% pos_btn_align_em_dn = pos_ax_align_em + [w3_+dw1_, 0, w4_-w3_, floor(h3_/3)-h3_];
% pos_text_align_em = pos_ax_align_em + [w3_+dw1_, h3_/3, w4_-w3_, floor(h3_/3)-h3_];
% pos_btn_align_em_up = pos_ax_align_em + [w3_+dw1_, h3_*2/3, w4_-w3_, floor(h3_/3)-h3_];
% 
% pos_ax_align_2p = pos_ax_align_em + [0, h3_+h1_, 0, 0];
% pos_btn_align_2p_dn = pos_ax_align_2p + [w3_+dw1_, 0, w4_-w3_, floor(h3_/3)-h3_];
% pos_text_align_2p = pos_ax_align_2p + [w3_+dw1_, h3_/3, w4_-w3_, floor(h3_/3)-h3_];
% pos_btn_align_2p_up = pos_ax_align_2p + [w3_+dw1_, h3_*2/3, w4_-w3_, floor(h3_/3)-h3_];



%% GUI layouts
try % close old GUI 
    close(get(obj.gui.text_scan, 'parent')); 
end 
handles.fig_main = figure('position', pos_fig);
handles.color_gray = ones(1,3)*0.94;
handles.color_cyan = [0, 1, 1];
handles.color_pink = [255, 153, 255]/255;

%% scan IDs
handles.text_scan = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_text_scan, 'string', 'scan ID', 'fontsize', font_size, ...
    'fontweight', 'bold', 'backgroundcolor', handles.color_gray);
handles.edit_scan = uicontrol(handles.fig_main, 'style', 'edit', 'units', 'pixels', ...
    'position', pos_edit_scan, 'string', num2str(obj.scan_id), 'fontweight', 'bold',...
    'callback', 'cb_edit_scan;');
handles.btn_scan_left = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_scan_left, 'string', '<<', ...
    'callback', 'cb_btn_scan_left;');
handles.btn_scan_right = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_scan_right, 'string', '>>', ...
    'callback', 'cb_btn_scan_right;');


%% slice IDs
handles.text_plane = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_text_plane, 'string', 'slice ID', 'fontsize', font_size,...
    'fontweight', 'bold', 'backgroundcolor', handles.color_gray);
handles.edit_slice = uicontrol(handles.fig_main, 'style', 'edit', 'units', 'pixels', ...
    'position', pos_edit_slice, 'string', num2str(obj.slice_id), 'fontweight', 'bold',...
    'callback', 'cb_edit_slice;');
handles.btn_slice_left = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_plane_left, 'string', '<<', ...
    'callback', 'cb_btn_slice_left;');
handles.btn_slice_right = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_plane_right, 'string', '>>', ...
    'callback', 'cb_btn_slice_right;');

%% block IDs
handles.text_block = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_text_block, 'string', 'block ID', 'fontsize', font_size,...
    'fontweight', 'bold', 'backgroundcolor', handles.color_gray);
handles.edit_block = uicontrol(handles.fig_main, 'style', 'edit', 'units', 'pixels', ...
    'position', pos_edit_block, 'string', num2str(obj.block_id), 'fontweight', 'bold',...
    'callback', 'cb_edit_block;');
handles.btn_block_left = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_block_left, 'string', '<<', ...
    'callback', 'cb_btn_block_left;');
handles.btn_block_right = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_block_right, 'string', '>>', ...
    'callback', 'cb_btn_block_right;');

% axes for visualizing spatial shapes
ax_slice = cell(obj.num_slices, 1);
ax_corr = cell(obj.num_slices,1);
ax_em = cell(obj.num_slices, 1);
for m=1:obj.num_slices
    ax_slice{m} = axes('parent', handles.fig_main, 'units', 'pixels','box', 'on', ...
        'position', pos_ax_slice{m}, 'ytick', [], 'xtick', []);
    ax_corr{m} = axes('parent',handles.fig_main, 'units', 'pixels','box', 'on', ...
        'position', pos_ax_corr{m}, 'ytick', [], 'xtick', []);
    ax_em{m} = axes('parent',handles.fig_main, 'units', 'pixels', 'box', 'on', ...
        'position', pos_ax_em{m}, 'ytick', [], 'xtick', []);
end
handles.ax_slice = ax_slice; 
handles.ax_corr = ax_corr; 
handles.ax_em = ax_em; 

%% load button
handles.btn_load_ca = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_load_video, 'string', 'load Ca', 'backgroundcolor', handles.color_cyan,...
    'fontweight', 'bold', 'fontsize', font_size, 'callback', 'cb_btn_load_ca;');


%% EM blurring
handles.text_blur = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_text_blur, 'string', 'z blur', 'fontsize', font_size, ...
    'fontweight', 'bold', 'backgroundcolor', handles.color_gray);
handles.edit_blur = uicontrol(handles.fig_main, 'style', 'edit', 'units', 'pixels', ...
    'position', pos_edit_blur, 'string', num2str(obj.em_zblur), 'fontweight', 'bold',...
    'callback', 'cb_edit_blur;');
handles.btn_blur_left = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_blur_left, 'string', '<<', ...
    'callback', 'cb_btn_blur_left;');
handles.btn_blur_right = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_blur_right, 'string', '>>', ...
    'callback', 'cb_btn_blur_right;');

handles.btn_load_em = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_load_em, 'string', 'load EM', 'backgroundcolor', handles.color_cyan,...
    'fontweight', 'bold', 'fontsize', font_size, 'callback', 'cb_btn_load_em;');

%% summary statistics
% buttons
% handles.btn_pick = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
%     'position', pos_btn_pick, 'string', 'get seed', 'backgroundcolor', handles.color_cyan, ...
%     'fontsize', font_size, 'fontweight', 'bold','callback', 'cb_btn_pick');

handles.btn_cn = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_cn, 'string', 'CN', 'backgroundcolor', handles.color_pink, ...
    'fontsize', font_size, 'fontweight', 'bold','callback', 'ease.nam_show=''cn'';cb_btn_slice');

handles.btn_pnr = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_pnr, 'string', 'PNR', 'backgroundcolor', handles.color_pink, ...
    'fontsize', font_size, 'fontweight', 'bold','callback', 'ease.nam_show=''pnr'';cb_btn_slice');
%
handles.btn_max = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_max, 'string', 'MAX', 'backgroundcolor', handles.color_pink, ...
    'fontsize', font_size, 'fontweight', 'bold','callback', 'ease.nam_show=''max'';cb_btn_slice');
%
handles.btn_mean = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_mean, 'string', 'MEAN', 'backgroundcolor', handles.color_pink, ...
    'fontsize', font_size, 'fontweight', 'bold','callback', 'ease.nam_show=''mean'';cb_btn_slice');
%
handles.btn_noise = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_sn, 'string', 'SN', 'backgroundcolor', handles.color_pink, ...
    'fontsize', font_size, 'fontweight', 'bold','callback', 'ease.nam_show=''sn'';cb_btn_slice');

%% zoom-in & zoom-out
handles.btn_zoomin = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_zoomin, 'string', 'zoom in', 'backgroundcolor', handles.color_gray,...
    'fontsize', font_size, 'fontweight', 'bold', 'callback', 'cb_btn_zoomin;');
handles.btn_zoomout = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_zoomout, 'string', 'zoom out', 'backgroundcolor', handles.color_gray,...
    'fontsize', font_size, 'fontweight', 'bold', 'callback', 'cb_btn_zoomout;');

%% initialize CNMF
handles.btn_auto_init = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_auto_init, 'string', 'add cells', 'backgroundcolor', handles.color_cyan,...
    'fontsize', font_size, 'fontweight', 'bold', 'callback', 'cb_add_neurons;');

handles.btn_init_options = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_init_options, 'string', 'options', 'backgroundcolor', handles.color_gray,...
    'fontsize', font_size, 'fontweight', 'bold', 'callback', 'cb_btn_init_options;');

handles.btn_update_cnmf = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_cnmf, 'string', 'run CNMF', 'backgroundcolor', handles.color_cyan,...
    'fontsize', font_size, 'fontweight', 'bold', 'callback', 'ease_run_cnmf;');

handles.btn_auto_del = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_auto_del, 'string', 'auto del.', 'backgroundcolor', handles.color_cyan,...
    'fontsize', font_size, 'fontweight', 'bold', 'callback', 'ease_auto_del;');

% handles.btn_load_init = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
%     'position', pos_btn_load_init, 'string', 'load init.', 'backgroundcolor', handles.color_gray,...
%     'fontsize', font_size, 'fontweight', 'bold', 'callback', 'ease_import_results;');

%% visulize neurons 
handles.btn_visualize_neuron = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_visualize_neurons, 'string', 'check cell',...
    'fontweight', 'bold', 'fontsize', font_size, 'background', handles.color_cyan, ...
    'callback', 'cell_id = 1; cb_btn_slice; ease_show_2p_neuron;');

handles.btn_ca_ahead = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_ca_ahead, 'string', '<<',...
    'callback', 'cell_id = max(1, cell_id-1); ease_show_2p_neuron;');
handles.text_cell_id = uicontrol(handles.fig_main, 'style', 'edit', 'units', 'pixels', ...
    'position', pos_text_cell_id, 'string', '0', 'fontweight', 'bold', 'fontsize', font_size,...
    'callback', 'cell_id= str2double(get(gco, ''string''));ease_show_2p_neuron;');
handles.btn_ca_next = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_ca_next, 'string', '>>',...
    'callback', 'cell_id=cell_id+1; ease_show_2p_neuron;');

handles.btn_confidence = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_confidence, 'string', 'score', 'fontweight', 'bold',...
    'fontsize', font_size, 'callback', 'ease_mark_deletion;');

handles.btn_rate_5 = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_rate_5, 'string', '5', 'fontsize', font_size, ...
    'callback', 'neuron.match_status.confidence(cell_id)=5; cb_btn_rate;');
handles.btn_rate_4 = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_rate_4, 'string', '4','fontsize', font_size,  ...
    'callback', 'neuron.match_status.confidence(cell_id)=4; cb_btn_rate;');
handles.btn_rate_3 = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_rate_3, 'string', '3', 'fontsize', font_size, ...
    'callback', 'neuron.match_status.confidence(cell_id)=3; cb_btn_rate;');
handles.btn_rate_2 = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_rate_2, 'string', '2', 'fontsize', font_size, ...
    'callback', 'neuron.match_status.confidence(cell_id)=2; cb_btn_rate;');
handles.btn_rate_1 = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_rate_1, 'string', '1', 'fontsize', font_size, ...
    'callback', 'neuron.match_status.confidence(cell_id)=1; cb_btn_rate;');
handles.btn_rate_0 = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_rate_0, 'string', '0', 'fontsize', font_size, ...
    'callback', 'neuron.match_status.confidence(cell_id)=0; cb_btn_rate;');

handles.btn_delete = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_delete, 'string', 'delete', 'fontweight', 'bold',...
    'fontsize', font_size, 'callback', 'cb_btn_delete;');

handles.btn_soma = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_soma, 'string', 'soma', 'fontweight', 'bold',...
    'fontsize', font_size, 'callback', 'neuron.labels(cell_id) = 1; cb_btn_label;');

handles.btn_dendrite = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_dendrite, 'string', 'dendrite', 'fontweight', 'bold',...
    'fontsize', font_size, 'callback', 'neuron.labels(cell_id) = 2; cb_btn_label;');

%% initialize a component given a seed pixel
% handles.btn_init_seed = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
%     'position', pos_btn_init_seed, 'string', 'init. seed', 'backgroundcolor', handles.color_cyan,...
%     'fontsize', font_size, 'fontweight', 'bold', 'callback', 'ease_initialize_ac;');

% handles.btn_match_em = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
%     'position', pos_btn_match_em, 'string', 'match EM', 'backgroundcolor', handles.color_cyan,...
%     'fontsize', font_size, 'fontweight', 'bold', 'callback', 'ease_calculate_match_score;');

handles.btn_em_ahead = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_ahead, 'string', '<<',...
    'callback', 'em_rank = max(1, em_rank-1); ease_check_match;');
handles.text_em_current = uicontrol(handles.fig_main, 'style', 'edit', 'units', 'pixels', ...
    'position', pos_text_em_current, 'string', '1', 'fontweight', 'bold', 'fontsize', font_size,...
    'callback', 'em_rank= str2double(get(gco, ''string''));ease_check_match;');
handles.btn_em_next = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_next, 'string', '>>',...
    'callback', 'em_rank=em_rank+1; ease_check_match;');
handles.check_em_only = uicontrol(handles.fig_main, 'style', 'checkbox', 'units', 'pixels', ...
    'position', pos_check_em,...
    'callback', 'ease.show_em_only= logical(get(gco, ''value'')); cb_btn_em_only;');
handles.btn_em_only = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_only, 'string', 'EM only', 'backgroundcolor', handles.color_gray, 'fontsize', font_size,...
    'callback', 'ease.show_em_only = xor(true, ease.show_em_only); cb_btn_em_only');

handles.btn_em_perfect = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_perfect, 'string', 'perfect', 'fontsize', font_size,...
    'callback', 'cb_perfect_match_callback; ');

handles.btn_em_candidate = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_candidate, 'string', 'candidate', 'fontsize', font_size,...
    'callback', 'cb_candidate_match_callback');

handles.btn_em_ignore = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_ignore, 'string', 'ignore', 'fontsize', font_size,...
    'callback', 'cb_ignore_match_callback');

handles.btn_em_zero = uicontrol(handles.fig_main, 'style', 'pushbutton', 'units', 'pixels', ...
    'position', pos_btn_em_zero, 'string', 'zero', 'fontsize', font_size,...
    'callback', 'cb_zero_match_callback');

%% the axes for showing activity
handles.ax_activity = axes('parent',handles.fig_main, 'units', 'pixels',...
    'position', pos_ax_activity);
title(handles.ax_activity, 'activity');

handles.ax_score = axes('parent',handles.fig_main, 'units', 'pixels',...
    'position', pos_ax_score, 'title', 'match score');
title(handles.ax_score, 'match score');

obj.gui = handles; 
end
