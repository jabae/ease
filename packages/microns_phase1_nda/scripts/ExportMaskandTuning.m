cellKey = fetch(nda.EMCell & 'em_id=362');
cellKey = fetch(nda.MeatCell * nda.AllenCell & nda.RFScore);
%cellKey = fetch(nda.EMCell * nda.AllenCell & nda.RFScore);
figure
for icell=1:length(cellKey)
    cellKey(icell).em_id=cellKey(icell).em_id-1; %% masks are incorrectly labeled with em_id-1
    maskKey = fetch(nda.Spike & nda.DriftResponse & cellKey(icell) );
    for imask=1:length(maskKey)
        %%
        subplot(length(maskKey),2,imask*2)
        s = fetch(nda.DriftResponse * nda.DriftTrial & maskKey(imask), 'direction','response');
        if isempty(s)
            continue
        end
        [responses, direction] = dj.struct.tabulate(s,'response', 'direction');
        angles = (0:size(responses,1)-1)/size(responses,1)*360;
        %subplot(nrows, ncols, icell)
        plot(angles, squeeze(responses)', 'k.')
        hold on
        m = nanmean(responses,2);
        s = nanstd(responses,[], 2)./sqrt(sum(~isnan(responses),2));
        errorbar(angles, m, s, 'r', 'LineWidth', 3)
        ylim([0 max(m)*4])
        xlim([0 360])
        box off
        keyTitle(maskKey(imask))
        
        %%
        subplot(length(maskKey),2,imask*2-1)
        g = flip(fetch1(preprocess.PrepareGalvoAverageFrame &...
            sprintf('animal_id=8973 and channel=1 and scan_idx=%d and slice=%d',maskKey(imask).scan_idx,maskKey(imask).slice),'frame'),1);
        r = flip(fetch1(preprocess.PrepareGalvoAverageFrame &...
            sprintf('animal_id=8973 and channel=2 and scan_idx=%d and slice=%d',maskKey(imask).scan_idx,maskKey(imask).slice),'frame'),1);
        
        g = g-min(g(:));
        g = g/quantile(g(:),.98);
        
        r = r-min(r(:));
        r = r/quantile(r(:),.99);
        
        g(g>1)=1; r(r>1)=1;
        g=g*.6;r=r*.6;
        
        g=imresize(g,2,'nearest');r=imresize(r,2,'nearest');
        b=ones(size(g));
        
        maskPixels = fetch1(nda.Mask & maskKey(imask),'mask_pixels');
        maskIm=zeros(size(g));
        im = zeros(256,256);
        im(maskPixels)=1;
        im=imresize(im,2,'nearest');
        maskIm=maskIm+bwperim(im);
        
        g=g+bwperim(im);
        r=r+bwperim(im);
        g(g>1)=1; r(r>1)=1;
        
        %maskIm(maskIm>1)=1;
        image(cat(3,r,g,maskIm));
        axis image off
        maskCenter = [median(find(mean(maskIm,1))) median(find(mean(maskIm,2)))];
        %xlim(maskCenter(1) + [-30 30])
        %ylim(maskCenter(2) + [-30 30])
        pos = get(gca,'position');
        set(gca,'position',pos + [0 -.015 .028 .028])
        
    end
    drawnow
    pause
    clf
    set(gcf,'position',[100 -50 500 1000],'color','w')
end

