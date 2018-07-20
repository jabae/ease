key.animal_id=8973;
key.session=1;

conditions=nan(12,27300);
k=1;
for i=[2 3 4 5 6 9 10 11 12]
    key.scan_idx=i;
    frameTimes = fetch1(preprocess.Sync & key, 'frame_times');
    frameTimes=frameTimes(1:3:end);
    
    trialOn = fetchn(nda.DriftTrial & key, 'onset');
    trialOff = fetchn(nda.DriftTrial & key, 'offset');
    dir = fetchn(nda.DriftTrial & key, 'direction');
    
    for j=1:length(trialOn)
        startInd = find(frameTimes>=trialOn(j),1,'first');
        stopInd = find(frameTimes<=trialOff(j),1,'last');
        test(k)=stopInd-startInd;
        conditions(i,startInd:stopInd) = rad2deg(dir(j));
        k=k+1;
    end
end
%%
%addAttribute(nda.Stimulus,'conditions : longblob # vector indicating direction of moving noise, or NaN for stationary synchronized with slice 1 frame times (degrees)')
key.animal_id=8973;
key.session=1;
for i=[2 3 4 5 6 9 10 11 12]
    key.scan_idx=i
    update(nda.Stimulus & key,'conditions',conditions(i,:));
end