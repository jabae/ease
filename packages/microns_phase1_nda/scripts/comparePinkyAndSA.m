pxPerMicron = 512/400;
clear x y z xloc yloc zloc
for i=1:570
    x(i)=c(i).xyz(1);
    y(i)=c(i).xyz(2);
    z(i)=c(i).xyz(3);
    xloc(i) = sa.loc(i).x;
    yloc(i) = sa.loc(i).y;
    zloc(i) = sa.loc(i).z;
end

%%
k=0
for i=1:length(c)
    if ~count(nda.MeatCell & ['em_id=' num2str(c(i).id)]) && c(i).inEM
        if count(nda.RF & ['em_id=' num2str(c(i).id)]) 
        k=k+1;
        end
    end
end
        
k


