%{
nda.DriftResponseRaw (computed) # calcium responses to drift trials
-> nda.DriftTrial
-> nda.Trace
-----
response: float  # averaged response
%}

classdef DriftResponseRaw < dj.Relvar & dj.AutoPopulate
    
    properties 
        popRel = nda.Trace;
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            frameTimes = fetch1(preprocess.Sync & key & 'animal_id=8973', 'frame_times');
            trace = double(fetch1(nda.Trace & key, 'trace'));
            frameTimes = frameTimes(key.slice:3:end);
            
            latency = 0.03;
            baselineWin = .5;
            trials = fetch(nda.DriftTrial & key, 'onset', 'offset', 'ORDER BY drift_trial');
            for trial = trials'
                tuple = dj.struct.join(key, rmfield(trial, {'onset', 'offset'}));
                baseline = median(trace(frameTimes > (trial.onset - baselineWin) & frameTimes < trial.onset));
                tuple.response = (median(trace(frameTimes > trial.onset+latency & frameTimes < trial.offset+latency))-baseline)./baseline;
                self.insert(tuple)
            end
        end
    end
    
    
end
