% summarizes the quality of mice
%set(groot,'defaultAxesColorOrder', parula(4)*0.7, 'defaultAxesLineStyleOrder', '-|--|:')
%key = fetch(nda.EMCell & nda.RF);
%key = fetch((nda.MeatCell*nda.AllenCell) & nda.RF);
key = fetch(nda.EMCell & (nda.RF * nda.Soma & 'stack_name="AllenEM"'));
clear s
for i=1:length(key)
    %key(i).em_id=key(i).em_id-1; %% masks are incorrectly labeled with em_id-1
    scores = fetchn(nda.RFScore & key(i), 'score');
    s(i) = median(scores);
end
[~,ind]=sort(s,'descend');
key=key(ind);

mx = quantile(s, 0.99);
mn = 0.0;
bins = linspace(mn, mx, 100);
h = cumsum(hist(s,bins));

plot(bins, h(:,end)-h, 'LineWidth', 1);
xlabel Score
title 'number of cells above RF Score'
grid on 
box off
xlim([mn mx])
hold on
%plot([1 1]*0.4, ylim, 'r-')
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






    