function EM_masks = get_em_masks(obj)
%% generate binary masks for indicating EM volumes in each scanning planes 
%{
    
%}

%% inputs: 
%{
%}

%% outputs: 
%{
    EM_masks: num_scan X num_slice cell array, each cell element indicate the
    spatial mask of EM volume on that plane. 
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

EM_masks = cell(obj.num_scans, obj.num_slices);
FOV_2p = obj.FOV_stack;
blur_size = 2*obj.em_zblur;
ssub = obj.dims_stack(1) / obj.dims_video(1);

% overlay nearby planes
for mscan=1:obj.num_scans
    for mslice = 1:obj.num_slices
        z = obj.video_zvals_updated(mscan, mslice);
        ind = max(z-blur_size, 1):min(z+blur_size, obj.dims_stack(3)); 
        
        tmp_em = obj.em_ranges(ind);
        tmp_em = cell2mat(tmp_em);
        if isempty(tmp_em)
            EM_masks{mscan, mslice} = false(obj.d1, obj.d2);
        else
            k = convhull(tmp_em(:,2), tmp_em(:,1));
            xi = (tmp_em(k,2)-FOV_2p(3))/ssub;
            yi = (tmp_em(k,1)-FOV_2p(1))/ssub;
            EM_masks{mscan, mslice} = poly2mask(xi, yi, obj.d1, obj.d2);
        end
    end
end

assignin('base', 'EM_masks', EM_masks);
end