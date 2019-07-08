%% spike detection

%% determine channels with IEDs
sub_labels = { 'RESP0401', 'RESP0435', 'RESP0458', 'RESP0478', 'RESP0502',...
    'RESP0574', 'RESP0589', 'RESP0608', 'RESP0621', 'RESP0699'};

for subj = 1:size(dataBase,2)
    switch dataBase(subj).subj
        case 'RESP0401'
            IEDchan = {'soc4','soc5','soc6','sta3','sta4','stv5','stv6','stv7'};% de kanalen met spikes soc4-6, sta3,4 stv5-7
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0435'
            IEDchan = {'fh13','fh14','fh15','fh16','fh20','fh21','fh22','fh23',...
                'fh24','fh26','fh27','fh28','fh30','fh31','fl03','fl04','fl05','fl12','fl13'};% FH13-16, 20-24, 26-28,30,31, FL3-5,12,13
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0458'
            IEDchan = {'f08','f09','f14','f17','f18','f19'};% F8-9,14,17-19
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0478'
            IEDchan = {'ihh04','ihh05','ihl01','ihl02','c05','c06','c12','c20'};% IHH4,5, IHL1,2, C5,6,12,20
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0502'
            IEDchan = {'t19','t20','t21','t22','t23','t28','t29','t30','t31'};% T19-23, 28-31
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0574'
            IEDchan = {'dv1','dv2','dv3','da1','da2','da3','f04','f05','hf10',...
                'hf11','hf12','hf13','hf14'};% Dv1-3, Da1-3, F4,5, HF10-14  CHECK:F4,5 WEL ECHT?
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0589'
            IEDchan = {'stm5','stm6','sta2','sta3','sta4'}; % Stm5,6, Sta2-4 CHECK: STA2,3,4 ECHT?
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0608'
            IEDchan = {'d1','d2','d3','c06','c07','c08','c15','c16','c21','c22'};  % D1-3, C6-8, 15, 16, 21, 22
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0621'
            IEDchan = {'c43','c44','c45','c53'};  % C43,44,45,53   NOG CHECKEN MET FRANS
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'RESP0699'
            IEDchan = {'c13','c14','c20','c21','c22','c23','c28','c29'};  % C13, 14, 20-23, 28,29
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
    end
    
    dataBase(subj).IEDch = IEDch;
    dataBase(subj).IEDchan = IEDchan;
end

%% remove stimulation artefact

for subj = 1:size(dataBase,2)
    data_noStimArt  = dataBase(subj).data;
    
    for eventnum = 1:size(dataBase(subj).tb_events,1)
        if strcmp(dataBase(subj).tb_events.sub_type(eventnum),'SPES')
            % determine uV 10 ms pre-stim and 10ms post-stim for each channel
            av_volt = mean([data_noStimArt(:,dataBase(subj).tb_events.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs)),...
                data_noStimArt(:,dataBase(subj).tb_events.sample_start(eventnum)+round(0.01*dataBase(subj).ccep_header.Fs))],2);
            
            % remove stimulation artefact in each channel
            data_noStimArt(:,dataBase(subj).tb_events.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs):...
                dataBase(subj).tb_events.sample_start(eventnum)+round(0.01*dataBase(subj).ccep_header.Fs))=...
                repmat(av_volt,1,size(-round(0.01*dataBase(subj).ccep_header.Fs):round(0.01*dataBase(subj).ccep_header.Fs),2));
            
            % remove saturated ECoG signal in stimulated channels
            % - determine stimulated channels
            stimchan = strsplit(dataBase(subj).tb_events.electrical_stimulation_site{eventnum},'-');
            for i=1:size(stimchan,2)
                stimnum(i) = find(strcmp(dataBase(subj).ch,stimchan{i})==1);
            end
            
            % - remove data during 10ms pre-stim and 5s post-stim in stimulated channels
            data_noStimArt(stimnum,dataBase(subj).tb_events.sample_start(eventnum)-round(0.01*dataBase(subj).ccep_header.Fs):...
                dataBase(subj).tb_events.sample_start(eventnum)+round(5*dataBase(subj).ccep_header.Fs))=...
                zeros(size(stimnum,2),1,size(round(-0.01*dataBase(subj).ccep_header.Fs):round(5*dataBase(subj).ccep_header.Fs),2));
        end
    end
    
    dataBase(subj).data_noStimArt = data_noStimArt;
    
    fprintf('...Subject %s has been run...\n',sub_labels{subj})
    
end

fprintf('...All subjects has been run...\n')

%% figure with 15s with 1 channel

subj = 1;
chan = 1;
eventnum = find(strcmp(dataBase(subj).tb_events.sub_type,'SPES')==1,1);
eventsamp = dataBase(subj).tb_events.sample_start(eventnum);

% total signal in time(s)
t = 1/dataBase(subj).ccep_header.Fs:1/dataBase(subj).ccep_header.Fs:size(dataBase(subj).data_noStimArt,2)/ccep_header.Fs;
% 15 seconds aruond eventnum
tt = eventsamp/dataBase(subj).ccep_header.Fs-5:1/dataBase(subj).ccep_header.Fs:eventsamp/dataBase(subj).ccep_header.Fs+10;

figure(1),
plot(t,dataBase(subj).data(chan,:),'k')
hold on
plot(t,dataBase(subj).data_noStimArt(chan,:),'r')
hold off
xlim([round(tt(1)), round(tt(end))])
ylim([-3000 3000])
title(sprintf('%s: electrode %s, stimulating %s',sub_labels{subj},dataBase(subj).ch{chan},dataBase(subj).tb_events.electrical_stimulation_site{eventnum}))

%% detect spikes - only in IED channels

for subj = 1:size(dataBase,2)
    data_IED = dataBase(subj).data_noStimArt(dataBase(subj).IEDch,:);
    
    dataBase(subj).data_IED = data_IED;
    
    [M,Mrange]=findCortSpikes(dataBase(subj).data_IED,dataBase(subj).ccep_header);
    
    dataBase(subj).spikes = M;
    
end

% uiteindelijk heb ik alleen de spikes loc nodig (M.groups{6,:})


%% check spikes visually in the first 10 minutes

subj = 1;

spesstim = strcmp(dataBase(subj).tb_events.sub_type,'SPES');
samp_start = dataBase(subj).tb_events.sample_start(spesstim);

for elec = 1%:size(dataBase(subj).IEDch,2)
    
    yall = [];
    xall=[];
    
    for win = 1:10*4-1 % each window has 15s, so 40 windows to visualize 10min
        t = 1/dataBase(subj).ccep_header.Fs:1/dataBase(subj).ccep_header.Fs:size(dataBase(subj).data_IED,2)/dataBase(subj).ccep_header.Fs;
        
        % window showed in plot
        tt = [(round((win-1)*15)*dataBase(subj).ccep_header.Fs), round(((win)*15)*dataBase(subj).ccep_header.Fs)];
        
        eventnum = find(samp_start>tt(1) & samp_start<tt(2));
        
        figure(1),
        plot(t,dataBase(subj).data_IED(elec,:),'k')
        
        hold on
        % vertical lines when stimulated
        plot(repmat(samp_start(eventnum)'/ccep_header.Fs,size(-1000:1000,2),1),...
            repmat(-1000:1000,size(eventnum,1),1)','r')
        
        %         % plot detected spikes
        %         locSpikes = [M(elec).groups{6,:}];
        %         valSpikes = [M(elec).groups{7,:}];
        %
        %         locnum = find(locSpikes >tt(1) & locSpikes<tt(2));
        %
        %         plot(locSpikes(locnum)/ccep_header.Fs,...
        %             dataBase(subj).data_IED(elec,locSpikes(locnum)),'*b')
        hold off
        
        ylim([-2000 2000])
        xlim([(win-1)*15, win*15])
        title(sprintf('%s: electrode %s',sub_labels{subj},dataBase(subj).IEDchan{elec}))
        
        [x,y]=ginput;
        
        xall = [xall; round(x*ccep_header.Fs)];
        yall = [yall; y];
        
    end
    
    visSpikes(elec).x = xall;
    visSpikes(elec).y = yall;

end


%% determine performance by FP, TP, FN





%% 


