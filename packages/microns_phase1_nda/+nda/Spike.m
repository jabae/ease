%{
nda.Spike (computed) # inferred spike rates from traces
-> nda.Trace
-----
rate 		: longblob 	# inferred spike rate (float vector)
%}

classdef Spike < dj.Relvar & dj.AutoPopulate
    properties
        popRel = nda.Trace;
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            dt = 1/fetch1(nda.ScanInfo & key, 'fps');
            tr = fetch1(nda.Trace & key, 'trace');
            key.rate = fast_oopsi(double(tr),struct('dt',dt),struct('lambda',.3));
            self.insert(key)
        end
    end
end
