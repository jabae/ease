% summarizes the quality of mice
%set(groot,'defaultAxesColorOrder', parula(4)*0.7, 'defaultAxesLineStyleOrder', '-|--|:')
key = fetch(nda.EMCell & nda.VonMises);
key = fetch(nda.MeatCell *nda.AllenCell& nda.VonMises);
key = fetch(nda.EMCell & (nda.VonMises * nda.Soma & 'stack_name="Pinky40"'));
clear osi
for i=1:length(key)
    %key(i).em_id=key(i).em_id-1; %% masks are incorrectly labeled with em_id-1
    s = fetch(nda.VonMises & key(i), 'von_base->w1','von_amp1->w2','von_amp2->w3','exp(-von_sharp)->r','von_pvalue->p');
    osi(i) = median(osifun(s));
    scores = fetchn(nda.RFScore & key(i), 'score');
    rfscore(i) = median(scores);
end
[~,ind]=sort(osi,'descend');
key=key(ind);

mx = quantile(osi, 0.99);
mn = 0.0;
bins = linspace(mn, mx, 100);
h = cumsum(hist(osi,bins));

plot(bins, h(:,end)-h, 'LineWidth', 1);
xlabel OSI
title 'number of cells above OSI'
grid on 
box off
xlim([mn mx])
hold on
plot([1 1]*0.4, ylim, 'r-')
%hold off
%%
%set(gcf,'PaperSize',[5 2.5], 'PaperPosition', [0 0 5 2.5])
%print('-dpdf',getLocalPath('~/Desktop/quality-tuning'))


plot(bins, 100-100*bsxfun(@rdivide, h, h(:,end))', 'LineWidth', 1);
xlabel OSI
title 'percent cells above OSI'
grid on 
box off
xlim([mn mx])
hold on
plot([1 1]*0.4, ylim, 'r-')
hold off

%set(gcf,'PaperSize',[5 2.5], 'PaperPosition', [0 0 5 2.5])
%print('-dpdf',getLocalPath('~/Desktop/quality-tuning-percent'))






    