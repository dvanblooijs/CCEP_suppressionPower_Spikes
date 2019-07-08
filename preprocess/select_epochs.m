function dataBase = select_epochs(dataBase,cfg)

epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

for subj = 1:size(dataBase,2)
    tt = round(epoch_length*dataBase(subj).ccep_header.Fs);
    
    % allocation
    cc_epoch_sorted = NaN(size(dataBase(subj).data,1),dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets,1),tt);
    
    for elec = 1:size(dataBase(subj).data,1) % for all channels
        for ll = 1:length(dataBase(subj).cc_stimsets) % for all epochs with >4 stimuli
            if strcmp(cfg.dir,'no')
                eventnum1 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}]));
                eventnum2 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,2}, '-',dataBase(subj).cc_stimchans{ll,1}]));
                eventnum = [eventnum1;eventnum2];
            elseif strcmp(cfg.dir,'yes')
                eventnum = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}]));
            end
            
            if size(eventnum,1) >dataBase(subj).max_stim
                events = dataBase(subj).max_stim;
            else
                events = size(eventnum,1);
            end
            
            for n=1:events
                cc_epoch_sorted(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
            end
        end
    end
    
    cc_epoch_sorted_avg = squeeze(nanmean(cc_epoch_sorted,2));
    
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    
    fprintf('...%s has been epoched and averaged... \n',dataBase(subj).subj)
end

fprintf('Data of all subjects have been epoched and averaged\n')

end
