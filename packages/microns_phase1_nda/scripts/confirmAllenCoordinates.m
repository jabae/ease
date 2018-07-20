%load('C:\work\microns\.mat files\pinky40lfinalV71oc.mat');
%load('C:\work\microns\.mat files\m8973stk.mat')
%stk=stk(:,:,:,310:-1:1);


for i=1:length(c)
    if count(nda.Trace & sprintf('em_id=%d',c(i).id))
        
        %% Display location in stack of Allen coordinates
        subplot(131)
        g = flip(squeeze(stk(:,:,1,round(c(i).xyz(3)+1))),1);
        r = flip(squeeze(stk(:,:,2,round(c(i).xyz(3)+1))),1);
        g = g-min(g(:)); g = g/quantile(g(:),.98);
        r = r-min(r(:)); r = r/quantile(r(:),.99);
        g(g>1)=1; r(r>1)=1; b=zeros(size(g));
        image(0:400,0:400,cat(3,r,g,b));
        hold on
        plot(c(i).xyz(1),c(i).xyz(2),'marker','o','markersize',8,'color','w','linewidth',3);
        axis image
        title(sprintf('Depth of cell %d', c(i).xyz(3)+1))
        
        %% Display location in functional slice for this cell
        subplot(132)
        mask = fetch(nda.Mask & sprintf('em_id=%d',c(i).id),'mask_pixels');
        
        % Choose mask with most pixels as soma
        ind = 1;
        maxPxls=0;
        for j=1:length(mask)
            if length(mask(j).mask_pixels)>maxPxls
                maxPxls=length(mask(j).mask_pixels);
                ind=j;
            end
        end
        mask = mask(ind);
        
        g = flip(fetch1(preprocess.PrepareGalvoAverageFrame & mask & 'animal_id=8973 and channel=1','frame'),1);
        r = flip(fetch1(preprocess.PrepareGalvoAverageFrame & mask & 'animal_id=8973 and channel=2','frame'),1);
        g = g-min(g(:)); g = g/quantile(g(:),.98);
        r = r-min(r(:)); r = r/quantile(r(:),.99);
        g(g>1)=1; r(r>1)=1;
        g=imresize(g,2,'nearest');r=imresize(r,2,'nearest');
        b=ones(size(g));
        
        maskIm=zeros(size(g));
        im = zeros(256,256);
        im(mask.mask_pixels)=1;
        im=imresize(im,2,'nearest');
        maskIm=maskIm+bwperim(im);
        image(cat(3,r,g,maskIm));
        axis image
        title('Functional scan')
        
        %% Display location in stack based on z-position of functional slice in stack
        subplot(133)
        depth = fetch1(nda.Scan & mask,'depth') + (mask.slice-1)*8;
        g = flip(squeeze(stk(:,:,1,depth)),1);
        r = flip(squeeze(stk(:,:,2,depth)),1);
        g = g-min(g(:)); g = g/quantile(g(:),.98);
        r = r-min(r(:)); r = r/quantile(r(:),.99);
        g(g>1)=1; r(r>1)=1; b=zeros(size(g));
        image(0:400,0:400,cat(3,r,g,b));
        hold on
        plot(c(i).xyz(1),c(i).xyz(2),'marker','o','markersize',8,'color','w','linewidth',3);
        axis image
        title(sprintf('Depth of slice %d', depth))
        
        pause
    end
end
        