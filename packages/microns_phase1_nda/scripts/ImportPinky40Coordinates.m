%uiopen('C:\Users\jreimer\Downloads\Baylor_Pinky40_xyz_um_flip_final_V7 (1).xlsx',1)
%import as cell array
dat = BaylorPinky40xyzumflipfinalV71;
%dat(1:2,:)=[];

clear c
for i=1:size(dat)
    c(i).id = -1*str2num(dat{i}(3:end));
    c(i).inEM = logical(dat{i,2});
    if c(i).inEM
        c(i).xyzEM = [dat{i,3} dat{i,4} dat{i,5}]; 
    else
        c(i).xyzEM = [nan nan nan];
    end
    c(i).xyz = [dat{i,6} dat{i,7} dat{i,8}]; 
end
    

%% add cells not in pinky40 to em coordinates
%uiopen('C:\Users\jreimer\Downloads\Baylor_Pinky40_xyz_um_flip_final_V7 (1).xlsx',1)
%import as cell array
dat = BaylorPinky40xyzumflipfinalV71;


clear c
for i=1:size(dat)
    c(i).id = -1*str2num(dat{i}(3:end));
    c(i).inEM = logical(dat{i,2});
    if c(i).inEM
        c(i).xyzEM = [dat{i,3} dat{i,4} dat{i,5}]; 
    else
        c(i).xyzEM = [nan nan nan];
    end
    c(i).xyz = [dat{i,6} dat{i,7} dat{i,8}]; 
end
    
    
    