function dataBase = removeStimArt(dataBase)

for nSubj = 1:size(dataBase,2)
    
    % raw data
    data_raw = dataBase(nSubj).data;

    fs = dataBase(nSubj).ccep_header.Fs;
     
    % after stimulation, the signal is returning from saturated to
    % baseline. With subtracting a smooth signal, the slow trend is
    % removed, which is visible in the signal in the stimulation trial
    % after stimulating this specific electrode.
    data_smooth = NaN(size(data_raw));
    for nCh = 1:size(data_raw,1)
        data_smooth(nCh,:) = smooth(data_raw(nCh,:),fs,'moving');
    end

    % subtracting the smoothed signal
    data_noStimArt = dataBase(nSubj).data_reref-data_smooth;
    
    % make bad channels = 0
    badch = strcmp(dataBase(nSubj).tb_channels.status,'bad');
    data_noStimArt(badch,:) = NaN(sum(badch),size(data_noStimArt,2));
    
    % find value when signal is saturated (max and min))
    max_val = fix(max(data_raw(:))/10)*10 -1;
    min_val = fix(min(data_raw(:))/10)*10 +1;
    
    % extract SPES events
    tb_events_stim = dataBase(nSubj).tb_events(strcmp(dataBase(nSubj).tb_events.trial_type,'electrical_stimulation') & strcmp(dataBase(nSubj).tb_events.sub_type,'SPESclin'),:);
    
    for nEvent = 1:size(tb_events_stim,1)
                    
            % determine uV 10 ms pre-stim and 10ms post-stim for each channel (not stimulated)
            avg_volt = mean([data_noStimArt(:,tb_events_stim.sample_start(nEvent)-round(0.01*fs)),...
                data_noStimArt(:,tb_events_stim.sample_start(nEvent)+round(0.01*fs))],2);
            
            % remove stimulation artefact in each channel
            data_noStimArt(:,tb_events_stim.sample_start(nEvent)-round(0.01*fs):...
                tb_events_stim.sample_start(nEvent)+round(0.01*fs))=...
                repmat(avg_volt,1,size(-round(0.01*fs):round(0.01*fs),2));           
            
            %% remove saturated ECoG signal in stimulated channels and other channels with saturated signal
            % - determine stimulated channels
            stimchan = strsplit(tb_events_stim.electrical_stimulation_site{nEvent},'-');
            stimnum = NaN(1,2);
            for nCh=1:size(stimchan,2)
                stimnum(nCh) = find(strcmp(dataBase(nSubj).ch,stimchan{nCh})==1);
            end
            
            % define sample of next event
            if nEvent < size(tb_events_stim,1)
                samp_nextevent = tb_events_stim.sample_start(nEvent+1);
                
            elseif tb_events_stim.sample_start(nEvent) + round(fs *20) < size(data_noStimArt,2) % if data after last stimulus is longer than 20s
                samp_nextevent = tb_events_stim.sample_start(nEvent) + round(fs *20);
                
            else                
                samp_nextevent = size(data_noStimArt,2);
            end
            
            % find other channels with saturated values
            idx_sat = (data_raw(:,tb_events_stim.sample_start(nEvent):samp_nextevent)<=min_val | ...
                data_raw(:,tb_events_stim.sample_start(nEvent):samp_nextevent)>=max_val);
            satch = find(sum(idx_sat,2)>0)';
            
            % make saturated values = 0
            fixch = unique([satch, stimnum]);
            for nCh = 1:size(fixch,2)
                % find between sample and next event when signal is not saturated anymore
                samp_stop = find(data_raw(fixch(nCh),tb_events_stim.sample_start(nEvent):samp_nextevent)<= min_val | ...
                    data_raw(fixch(nCh),tb_events_stim.sample_start(nEvent):samp_nextevent)>= max_val,1,'last');
                
                if isempty(samp_stop)
                    data_noStimArt(fixch(nCh),tb_events_stim.sample_start(nEvent)-round(0.01*fs):...
                        tb_events_stim.sample_start(nEvent) + round(fs * 5))=...
                        zeros(size(fixch(nCh),2),size(round(-0.01*fs):round(fs * 5),2));
                else
                    % - remove data during 10ms pre-stim and post-stim in
                    % stimulated channels until signal is not saturated anymore
                    data_noStimArt(fixch(nCh),tb_events_stim.sample_start(nEvent)-round(0.01*fs):...
                        tb_events_stim.sample_start(nEvent) + samp_stop)=...
                        zeros(size(fixch(nCh),2),size(round(-0.01*fs):samp_stop,2));
                end
            end
    end
    
    dataBase(nSubj).data_rerefnoStimArt = data_noStimArt;
    
    fprintf('...Subject %s has been run...\n',dataBase(nSubj).sub_label)
    
end
