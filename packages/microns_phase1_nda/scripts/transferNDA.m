%% nda.Scan
key = fetch(experiment.Scan * experiment.ScanLaser & 'animal_id=8973',...
    'depth','power->laser_power','wavelength','filename')
insert(nda.Scan, rmfield(key,{'animal_id','session'}))
%% nda.ScanInfo
key = fetch(preprocess.PrepareGalvo & 'animal_id=8973',...
    'nframes','px_width','px_height','um_width','um_height','bidirectional','fps','zoom',...
    'nchannels','nslices','fill_fraction','raster_phase')
insert(nda.ScanInfo, rmfield(key,{'animal_id','session'}))
%% nda.Slice
key = fetch(nda.Scan)
for i=1:length(key)
    k=key(i);
    for j=1:3
        k.slice=j;
        k.z_offset=(j-1)*8;
        insert(nda.Slice,k)
    end
end
%% nda.Treadmill
key = fetch(nda.Scan)
for i=1:length(key)
    [bt,bv] = fetch1(preprocess.Treadmill & key(i) & 'animal_id=8973','treadmill_time','treadmill_vel');
    ft = fetch1(preprocess.BehaviorSync & key(i) & 'animal_id=8973','frame_times');
    bv = interp1(bt,abs(bv),ft(1:3:end));
    key(i).treadmill_speed = bv;
end
insert(nda.Treadmill,key)
%% nda.Pupil
key = fetch(nda.Scan)
for i=1:length(key)
    et = fetch1(preprocess.Eye & key(i) & 'animal_id=8973','eye_time');
    [xy,r] = fetchn(preprocess.EyeTrackingFrame * preprocess.Eye & key(i) & 'animal_id=8973',...
        'center','major_r'); % Specify eye quality?
    clear x y
    for j=1:length(xy)
        if isempty(xy{j})
            x(j)=nan; y(j)=nan;
        else
            x(j)=xy{j}(1); y(j)=xy{j}(2);
        end
    end
    ft = fetch1(preprocess.BehaviorSync & key(i) & 'animal_id=8973','frame_times');
    key(i).pupil_r = interp1(et,r,ft(1:3:end),'nearest');
    key(i).pupil_x = interp1(et,x,ft(1:3:end),'nearest');
    key(i).pupil_y = interp1(et,y,ft(1:3:end),'nearest');
end
insert(nda.Pupil,key)
%% nda.EMCell
clear key
for i=0:569
    key.em_id=i;
    insert(nda.EMCell,key);
end
%% nda.Stimulus
key = fetch(nda.Scan)
for i=1:length(key)
    k=key(i);
    k.movie = fetch1(preprocess.SyncedMonet & key(i) & 'animal_id=8973','movie');
    insert(nda.Stimulus, k);
end

%% nda.Mask
temp=load('C:\work\ssMaskFINAL.mat')
masks = temp.masks;
scans = fetchn(nda.Scan,'scan_idx');
for i=1:570
    for j=1:9
        for w=1:3
            if ~isempty(masks{i,j,w})
                key.scan_idx=scans(j);
                key.slice=w;
                key.em_id=fetch1(nda.AllenCell & ['em_id=' num2str(i)],'allen_id');
                key.mask_pixels=masks{i,j,w};
                insert(nda.Mask,key)
            end
        end
    end
end
%% nda.Trace
key=fetch(nda.Scan);
blockSize=500;
nframes=27300;
for i=1:length(key) 
    fprintf('\tScan %i\n', i);
    reader = nda.getReader(key(i));
    for islice = 1:3
        key(i).slice=islice;
        fixMotion = get_fix_motion_fun(preprocess.PrepareGalvoMotion & key(i) & 'animal_id=8973');
        fixRaster = get_fix_raster_fun(preprocess.PrepareGalvo & key(i) & 'animal_id=8973');
        
        maskKeys = fetch(nda.Mask & key(i));
        traceKeys = maskKeys;
        for imask = 1:length(traceKeys)
            traceKeys(imask).trace = zeros(1,nframes,'int16');
        end
        for start = 1:blockSize:nframes
            blockInd = start:min(nframes, start + blockSize - 1);
            fprintf('\tComputing traces for frames %i:%i\n', blockInd(1), blockInd(end));
            
            block=squeeze(reader(:,:,1,islice, blockInd));
            for w=1:length(blockInd)
                %block(:,:,w)=fixRaster(block(:,:,w));
                block(:,:,w)=fixMotion(fixRaster(block(:,:,w)),blockInd(w));
            end
            block = flip(block,1);
            a=mean(squeeze(block),3);
            block = reshape(block,256*256,length(blockInd));
            for imask = 1:length(maskKeys)
                pxls = fetch1(nda.Mask & maskKeys(imask),'mask_pixels');
                a(pxls)=a(pxls)*2;
                traceKeys(imask).trace(blockInd)=mean(block(pxls,:));
            end
            imagesc(a)
            drawnow
        end
        insert(nda.Trace,traceKeys);
    end
    clear reader
end