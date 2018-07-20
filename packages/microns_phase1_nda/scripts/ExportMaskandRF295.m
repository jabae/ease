cellKey = fetch(nda.EMCell & 'em_id in (295)'); scans=[3 4 5 6 9];
%cellKey = fetch(nda.EMCell & 'em_id in (514,541,496,70,206,295)')
%cellKey = fetch(nda.MeatCell & nda.RFScore);
%cellKey = fetch((nda.EMCell & nda.RFScore) - nda.MeatCell);
clear s
for i=1:length(cellKey)
    scores = fetchn(nda.RFScore & cellKey(i), 'score');
    s(i) = median(scores);
end
[~,ind]=sort(s,'descend');
cellKey=cellKey(ind);


for icell=1:length(cellKey)
    maskKey = fetch(nda.RF & cellKey(icell));
    clear meanRF

    for imask=1:length(maskKey)
        k=find(maskKey(imask).scan_idx==scans)+1;
        %%
       
        rf = fetch1(nda.RF & maskKey(imask),'rf');
        rf = squeeze(mean(rf(1:3,:,:)));
        meanRF(imask,:,:) = rf;
        imagesc(rf)
        axis image off
        set(gcf,'color','w')
        drawnow
        saveas(gcf,['RF_scan' num2str(k) '_slice_' num2str(maskKey(imask).slice) '.png']);
        
        
        %%
        
        g = flip(fetch1(preprocess.PrepareGalvoAverageFrame &...
            sprintf('animal_id=8973 and channel=1 and scan_idx=%d and slice=%d',maskKey(imask).scan_idx,maskKey(imask).slice),'frame'),1);
        r = flip(fetch1(preprocess.PrepareGalvoAverageFrame &...
            sprintf('animal_id=8973 and channel=2 and scan_idx=%d and slice=%d',maskKey(imask).scan_idx,maskKey(imask).slice),'frame'),1);
        
        g = g-min(g(:));
        g = g/quantile(g(:),.8);
        
        r = r-min(r(:));
        r = r/quantile(r(:),.9);
        
        g(g>1)=1; r(r>1)=1;
        g=g*.7;r=r*.6;
        
        g=imresize(g,2,'nearest');r=imresize(r,2,'nearest');
        b=ones(size(g));
        
        maskPixels = fetch1(nda.Mask & maskKey(imask),'mask_pixels');
        maskIm=zeros(size(g));
        im = zeros(256,256);
        im(maskPixels)=.1;
        im=imresize(im,2,'nearest');
        maskIm=maskIm+bwperim(im);
        
        g=g+bwperim(im);
        r=r+bwperim(im);
        g(g>1)=1; r(r>1)=1;
        
        %maskIm(maskIm>1)=1;
        image(cat(3,r,g,maskIm));
        axis image off
        maskCenter = [median(find(mean(maskIm,1))) median(find(mean(maskIm,2)))];
        xlim(maskCenter(1) + [-30 30])
        ylim(maskCenter(2) + [-30 30])
        %pos = get(gca,'position');
        %set(gca,'position',pos + [0 -.015 .028 .028])
        set(gcf,'color','w')
        drawnow
        saveas(gcf,['Mask_scan' num2str(k) '_slice_' num2str(maskKey(imask).slice) '.png']);
        clf
    end
    %set(gcf,'position',[100 -50 500 1000],'color','w')
    figure(200)
    imagesc(squeeze(mean(meanRF)));
    axis image off
    drawnow
    saveas(gcf,'MeanRF.png');
    %pause
    %figure(100);clf
    %figure(200);clf
end

