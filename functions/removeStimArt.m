function dataBase = removeStimArt(dataBase)

for subj = 1:size(dataBase,2)
    
    data = dataBase(subj).data;

    % pre-allocation
    data_noStimArt = dataBase(subj).data_reref;
    
    % make bad channels = 0
    badch = strcmp(dataBase(subj).tb_channels.status,'bad');
    data_noStimArt(badch,:) = NaN(sum(badch),size(data_noStimArt,2));

    % make signal 0 when saturated 
    max_val = fix(max(dataBase(subj).data(:))/10)*10 -1;
    min_val = fix(min(dataBase(subj).data(:))/10)*10 +1;
    
    % extract SPES events
    tb_events_stim = dataBase(subj).tb_events(strcmp(dataBase(subj).tb_events.trial_type,'electrical_stimulation') & strcmp(dataBase(subj).tb_events.sub_type,'SPESclin'),:);
    
    for eventnum = 1:size(tb_events_stim,1)
                    
            % determine uV 10 ms pre-stim and 10ms post-stim for each channel (not stimulated)
            avg_volt = mean([data_noStimArt(:,tb_events_stim.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs)),...
                data_noStimArt(:,tb_events_stim.sample_start(eventnum)+round(0.01*dataBase(subj).ccep_header.Fs))],2);
            
            % remove stimulation artefact in each channel
            data_noStimArt(:,tb_events_stim.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs):...
                tb_events_stim.sample_start(eventnum)+round(0.01*dataBase(subj).ccep_header.Fs))=...
                repmat(avg_volt,1,size(-round(0.01*dataBase(subj).ccep_header.Fs):round(0.01*dataBase(subj).ccep_header.Fs),2));           
            
            % remove saturated ECoG signal in stimulated channels and other
            % channels with saturated signal
            % - determine stimulated channels
            stimchan = strsplit(tb_events_stim.electrical_stimulation_site{eventnum},'-');
            stimnum = NaN(1,2);
            for i=1:size(stimchan,2)
                stimnum(i) = find(strcmp(dataBase(subj).ch,stimchan{i})==1);
            end
            
            % define sample of next event
            if eventnum < size(tb_events_stim,1)
                samp_nextevent = tb_events_stim.sample_start(eventnum+1);
                
            elseif tb_events_stim.sample_start(eventnum) + round(dataBase(subj).ccep_header.Fs *20) < size(data_noStimArt,2) % if data after last stimulus is longer than 20s
                samp_nextevent = tb_events_stim.sample_start(eventnum) + round(dataBase(subj).ccep_header.Fs *20);
                
            else                
                samp_nextevent = size(data_noStimArt,2);
            end
            
            % find other channels with saturated values
            idx_sat = (data(:,tb_events_stim.sample_start(eventnum):samp_nextevent)<=min_val | ...
                data(:,tb_events_stim.sample_start(eventnum):samp_nextevent)>=max_val);
            satch = find(sum(idx_sat,2)>0)';
            
            % make saturated values = 0
            fixch = unique([satch, stimnum]);
            for i=1:size(fixch,2)
                % find between sample and next event when signal is not saturated anymore
                samp_stop = find(data(fixch(i),tb_events_stim.sample_start(eventnum):samp_nextevent)<= min_val | ...
                    data(fixch(i),tb_events_stim.sample_start(eventnum):samp_nextevent)>= max_val,1,'last');
                
                if isempty(samp_stop)
                    data_noStimArt(fixch(i),tb_events_stim.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs):...
                        tb_events_stim.sample_start(eventnum) + round(dataBase(subj).ccep_header.Fs * 5))=...
                        zeros(size(fixch(i),2),size(round(-0.01*dataBase(subj).ccep_header.Fs):round(dataBase(subj).ccep_header.Fs * 5),2));
                else
                    % - remove data during 10ms pre-stim and post-stim in
                    % stimulated channels until signal is not saturated anymore
                    data_noStimArt(fixch(i),tb_events_stim.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs):...
                        tb_events_stim.sample_start(eventnum) + samp_stop)=...
                        zeros(size(fixch(i),2),size(round(-0.01*dataBase(subj).ccep_header.Fs):samp_stop,2));
                end
            end
    end
    
    dataBase(subj).data_rerefnoStimArt = data_noStimArt;
    
    fprintf('...Subject %s has been run...\n',dataBase(subj).sub_label)
    
end
