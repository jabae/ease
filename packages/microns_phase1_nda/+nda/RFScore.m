%{
nda.RFScore (computed) # quality score of receptive fields
-> nda.RF
-----
score  : float      # quality score
%}

classdef RFScore < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = nda.RF;
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            roi=false(90,160);
            roi(30:75,60:120)=true;
            rf = fetch1(nda.RF & key,'rf');
            rf=rf-min(rf(:));
            for w=1:5
                inside = squeeze(rf(w,:,:));
                inside(~roi)=nan;
                outside = squeeze(rf(w,:,:));
                outside(roi)=nan;
                
                score(w) = nanstd(inside(:))./nanstd(outside(:));
            end
            
            tuple = key;
            tuple.score = mean(score(1:3));
            self.insert(tuple);
            
        end
    end
end

