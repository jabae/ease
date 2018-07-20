%{
nda.RF (computed) # receptive fields from Monet stimulus
-> nda.Stimulus
-> nda.Spike
-----
rf  : longblob      # receptive field map
%}

classdef RF < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = nda.Stimulus;
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            nlags = 5;
            %key = fetch(nda.Spike & 'em_id=347 and scan_idx=3 and slice=3')
            [rate,spikeKey] = fetchn(nda.Spike & key,'rate');
            fprintf('Loading movie...\n');
            mov = fetch1(nda.Stimulus & key,'movie');
            sz=size(mov);
            mov = reshape(mov,[],27300);
            mov(isnan(mov))=0;
            for i=1:length(spikeKey)
                fprintf('Scan %d, cell %d of %d...\n',key.scan_idx,i,length(spikeKey))
                rf = conv2(mov',rate{i}(end:-1:nlags),'valid');
                rf = flipdim(reshape(rf,nlags,sz(1),sz(2)),1);
                tuple = spikeKey(i);
                tuple.rf=rf;
                self.insert(tuple);
                
                clf
                for i=1:nlags
                    subplot(1,5,i)
                    imagesc(squeeze(rf(i,:,:)));
                    axis image off
                end
                drawnow
            end
        end
    end
end
            
            
