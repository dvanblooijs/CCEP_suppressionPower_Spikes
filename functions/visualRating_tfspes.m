function dataBase = visualRating_tfspes(dataBase,stimps, chans)

dataBase.ERSP.checked = NaN(size(dataBase(1).cc_epoch_sorted,1),size(dataBase(1).cc_epoch_sorted,3) );
for stimp= stimps
        for chan = chans
            figure(1),
            a=dataBase(1).ERSP.allERSPboot{stimp,chan};
            a=flipud(a);
            imagesc(dataBase(1).ERSP.times,dataBase(1).ERSP.freqs,a,[-15,15])
            title(sprintf('Stimulated: %s-%s, response in: %s',dataBase(1).cc_stimchans{stimp,1},dataBase(1).cc_stimchans{stimp,2},dataBase(1).ch{chan}))
            xlabel('Time(ms)')
            ax = gca;
            ax.YTick = min(dataBase(1).ERSP.freqs):50:max(dataBase(1).ERSP.freqs);
            ax.YTickLabels = max(dataBase(1).ERSP.freqs):-50:min(dataBase(1).ERSP.freqs);
            ylabel('Frequency (Hz)')
            colormap jet
            x = input('BB? [y/n] ','s');
            if strcmp(x,'y')
                dataBase.ERSP.checked(chan,stimp) = 1 ;
            else
                dataBase.ERSP.checked(chan,stimp) = 0 ;
            end
        end
end


end
