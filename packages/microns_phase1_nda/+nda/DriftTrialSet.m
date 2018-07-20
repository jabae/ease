%{
nda.DriftTrialSet (computed) # all drift trials for this scan
-> preprocess.Sync
---
%}


classdef DriftTrialSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = preprocess.Sync & (psy.MovingNoise & 'speed>0') & 'animal_id=8973';
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            self.insert(key)
            iTrial = 0;
            for key = fetch(psy.Trial * preprocess.Sync & psy.MovingNoise & key)'
                iTrial = makeTuples(nda.DriftTrial, key, iTrial);
            end
        end
        
    end
end