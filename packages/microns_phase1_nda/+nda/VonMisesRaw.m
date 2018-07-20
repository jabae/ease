%{
nda.VonMisesRaw (computed) # directional tuning
-> nda.Trace
-----
von_r2     : float  # r-squared explaned by vonMises fit
von_pref   : float  #  preferred directions
von_base   : float  #  von mises base value
von_amp1   : float  #  amplitude of first peak
von_amp2   : float  #  amplitude of second peak
von_sharp  : float  #  sharpnesses
von_pvalue : float  # p-value by shuffling (nShuffles = 1e4)
%}

classdef VonMisesRaw < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = nda.Trace & nda.DriftResponseRaw;
    end
    
    
    methods
        function plot(self)
            cellKey = fetch(nda.EMCell & self);
            ncells = length(cellKey);
            ncols = ceil(sqrt(ncells));
            nrows = ceil(ncells/ncols);
            for icell=1:length(cellKey)
            figure
                maskKey = fetch(nda.Trace & nda.DriftResponseRaw & cellKey(icell) );
                for imask=1:length(maskKey)
                    subplot(1,length(maskKey),imask)
                    s = fetch(nda.DriftResponseRaw * nda.DriftTrial & maskKey(imask), 'direction','response');
                    if isempty(s)
                        continue
                    end
                    [responses, direction] = dj.struct.tabulate(s,'response', 'direction');
                    angles = (0:size(responses,1)-1)/size(responses,1)*360;
                    %subplot(nrows, ncols, icell)
                    plot(angles, squeeze(responses)', 'k.')
                    hold on
                    m = nanmean(responses,2);
                    s = nanstd(responses,[], 2)./sqrt(sum(~isnan(responses),2));
                    errorbar(angles, m, s, 'r', 'LineWidth', 3)
                    ylim([-.2 1])
                    xlim([0 360])
                    box off
                    keyTitle(maskKey(imask))
                end
            end
        end
    end
    
    
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            s = fetch(nda.DriftResponseRaw*nda.DriftTrial & key, 'direction','response');
            [r, direction] = dj.struct.tabulate(s,'response', 'direction');
            responses(1,:,:)=r;
            assert(all(diff(direction)>0), 'just a check that directions are monotonic')
            nShuffles = 1e4;
            [von, r2, p] = ne7.rf.VonMises2.computeSignificance(double(responses), nShuffles);
            r2(isnan(r2)) = 0;
            tuple = key;
            tuple.von_r2 = r2;
            tuple.von_base = von.w(1);
            tuple.von_amp1 = von.w(2);
            tuple.von_amp2 = von.w(3);
            tuple.von_sharp= von.w(4);
            tuple.von_pref = von.w(5);
            tuple.von_pvalue = p;
            self.insert(tuple)
        end
        
    end
end