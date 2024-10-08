function dataBase = preprocess_ECoG_ccep(dataBase,cfg)

epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

if exist('minstim','var') == 0
    minstim = 5;
else
    minstim = cfg.minstim;
end

for nSubj = 1:size(dataBase,2)
    
    %% unique stimulation pairs
    stimpair = dataBase(nSubj).tb_events.electrical_stimulation_site(contains(dataBase(nSubj).tb_events.sub_type,'SPES') & ~contains(dataBase(nSubj).tb_events.electrical_stimulation_site,'n/a')) ;

    stimnum = NaN(size(stimpair,1),2);
    for stimp = 1:size(stimpair,1)
        stimchans = strsplit(stimpair{stimp},'-');
        for chan = 1:2
            stimnum(stimp,chan) = find(strcmp(stimchans{chan},dataBase(nSubj).ch)==1);
        end
    end
          
    stimelek = sort(stimnum,2);    
    [cc_stimsets,~,IC] = unique(stimelek,'rows');
    
    n = histcounts(IC,'BinMethod','integers');
    
    if any(diff(n) ~= 0)
        stimremove = find(n<minstim); % remove al stimulation pairs that are stimulated less than 5 times
        
        stimelek(IC==stimremove,:) = [];
        
        [cc_stimsets,~,IC] = unique(stimelek,'rows');
        n = histcounts(IC,'BinMethod','integers');
        if any(diff(n) ~= 0)
            fprintf('ERROR: %s some stimulation pairs are stimulated less/more than all others\n',dataBase(nSubj).sub_label)
        end        
    end
    
    cc_stimchans = cell(size(cc_stimsets,1),2);
    
    for stimp = 1:size(cc_stimsets,1)
        for chan =1:2
            cc_stimchans{stimp,chan} = dataBase(nSubj).ch{cc_stimsets(stimp,chan)};
        end
    end
    
    max_stim = median(n);
    
    dataBase(nSubj).cc_stimsets = cc_stimsets;
    dataBase(nSubj).cc_stimchans = cc_stimchans;
    dataBase(nSubj).max_stim = max_stim;
    
    stimdif = find(n ~= max_stim);
    for stimp =1:size(stimdif,2)
        [cc_stimchans{stimdif(stimp),1} '-' cc_stimchans{stimdif(stimp),2}] %#ok<NOPRT>
    end
    
    %% select epochs
    t = round(epoch_length*dataBase(nSubj).ccep_header.Fs);
    
    % allocation
    cc_epoch_sorted = NaN(size(dataBase(nSubj).data,1),dataBase(nSubj).max_stim,size(dataBase(nSubj).cc_stimsets,1),t);
    tt_epoch_sorted = NaN(dataBase(nSubj).max_stim,size(dataBase(nSubj).cc_stimsets,1),t); % samplenumbers for each epoch
    
    for elec = 1:size(dataBase(nSubj).data,1) % for all channels
        for ll = 1:length(dataBase(nSubj).cc_stimsets) % for all epochs with >4 stimuli
                eventnum1 = find(strcmp(dataBase(nSubj).tb_events.electrical_stimulation_site,[dataBase(nSubj).cc_stimchans{ll,1}, '-',dataBase(nSubj).cc_stimchans{ll,2}]));
                eventnum2 = find(strcmp(dataBase(nSubj).tb_events.electrical_stimulation_site,[dataBase(nSubj).cc_stimchans{ll,2}, '-',dataBase(nSubj).cc_stimchans{ll,1}]));
                eventnum = [eventnum1;eventnum2];
            
            if size(eventnum,1) >dataBase(nSubj).max_stim
                events = dataBase(nSubj).max_stim;
            else
                events = size(eventnum,1);
            end
            
            for n = 1:events
                cc_epoch_sorted(elec,n,ll,:) = dataBase(nSubj).data(elec,dataBase(nSubj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(nSubj).ccep_header.Fs)+1:dataBase(nSubj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(nSubj).ccep_header.Fs));
                tt_epoch_sorted(n,ll,:) = dataBase(nSubj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(nSubj).ccep_header.Fs)+1:dataBase(nSubj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(nSubj).ccep_header.Fs);
            end
        end
    end
    
    cc_epoch_sorted_avg = squeeze(mean(cc_epoch_sorted,2,'omitnan'));
    
    dataBase(nSubj).cc_epoch_sorted = cc_epoch_sorted;
    dataBase(nSubj).tt_epoch_sorted = tt_epoch_sorted;
    dataBase(nSubj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    
    fprintf('...%s has been epoched and averaged... \n',dataBase(nSubj).sub_label)
    
end