
function dataBase = visualRating_ccep(dataBase,cfg)

for subj = 1:size(dataBase,2)
fs = dataBase(subj).ccep_header.Fs;

tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

for stimp = 1:size(dataBase(subj).cc_epoch_sorted_avg,2)
    visERs = [];
    for chan =1 :size(dataBase(subj).cc_epoch_sorted_avg,1)
        
        if ismember(chan,dataBase(subj).ERs(stimp).detERs)
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.11 0.77 0.8];
            
            subplot(1,2,1)
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:)),'r');
            hold on
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s-%s',dataBase(subj).ch{chan},dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2}))
            
            subplot(1,2,2)
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:)),'r');
            hold on
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
            hold off
            xlim([-0.5 1])
            ylim([-750 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
            x = input('ER? [y/n] ','s');
            if strcmp(x,'y')
                visERs = [visERs, chan] ;
            end
        end
        
    end
    
    dataBase(subj).ERs(stimp).checked = visERs;
end
end
