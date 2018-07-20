scans = fetchn(experiment.Scan & 'animal_id=8973','scan_idx');
for i=1:9
    for j=1:3
        g = flip(fetch1(preprocess.PrepareGalvoAverageFrame &...
            sprintf('animal_id=8973 and channel=1 and scan_idx=%d and slice=%d',scans(i),j),'frame'),1);
        r = flip(fetch1(preprocess.PrepareGalvoAverageFrame &...
            sprintf('animal_id=8973 and channel=2 and scan_idx=%d and slice=%d',scans(i),j),'frame'),1);
        
        g = g-min(g(:));
        g = g/quantile(g(:),.98);
        
        r = r-min(r(:));
        r = r/quantile(r(:),.99);
        
        g(g>1)=1; r(r>1)=1;

        g=imresize(g,2,'nearest');r=imresize(r,2,'nearest');
        b=ones(size(g));
        
        
        maskIm=zeros(size(g));
        for w=1:570
            if ~isempty(masks{w,i,j})
                im = zeros(256,256);
                im(masks{w,i,j})=1;
                im=imresize(im,2,'nearest');
                maskIm=maskIm+bwperim(im);
            end
        end
        %maskIm(maskIm>1)=1;
        image(cat(3,r,g,maskIm));
        axis image off
        pause
        imwrite(cat(3,r,g,maskIm),sprintf('c:\\work\\exportFigs\\scan%d-%d.jpg',i,j))
        clf
    end
end
