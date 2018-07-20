%{
nda.DriftResponse (computed) # calcium responses to drift trials
-> nda.DriftTrial
-> nda.Spike
-----
response: float  # averaged response
%}

classdef DriftResponse < dj.Relvar & dj.AutoPopulate
    
    properties 
        popRel = nda.Spike;
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            frameTimes = fetch1(preprocess.Sync & key & 'animal_id=8973', 'frame_times');
            rate = fetch1(nda.Spike & key, 'rate');
            frameTimes = frameTimes(key.slice:3:end);
            
            latency = 0.03;
            trials = fetch(nda.DriftTrial & key, 'onset', 'offset', 'ORDER BY drift_trial');
            for trial = trials'
                if mod(trial.drift_trial, 20)==0 || trial.drift_trial==trials(end).drift_trial
                    fprintf('Trial %3d/%d\n', trial.drift_trial, trials(end).drift_trial)
                end
                tuple = dj.struct.join(key, rmfield(trial, {'onset', 'offset'}));
                tuple.response = mean(rate(frameTimes > trial.onset+latency & frameTimes < trial.offset+latency));
                self.insert(tuple)
            end
        end
    end
    
    
end
