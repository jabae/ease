%% nda.RFScore
clear all
key = fetch(nda.RFScore);
for i=1:length(key)
    tuple(i)=fetch(nda.RFScore & key(i),'*');
    key(i).em_id=key(i).em_id+1;
    tuple(i).em_id=fetch1(nda.AllenCell & key(i),'allen_id');
end
delQuick(nda.RFScore);
for i=1:length(tuple)
    fprintf('Inserting %d \n',i);
    insert(nda.RFScore,tuple(i));
end
%%
