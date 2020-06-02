%% spike detection
% author: Dorien van Blooijs
% date: july 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath(genpath('git_rep/eeglab/'))
addpath('git_rep/fieldtrip/')
ft_defaults

%% settings
cfg.dataPath = '/Fridge/CCEP';
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', 'sub-RESP0502',...
    'sub-RESP0574', 'sub-RESP0589', 'sub-RESP0608', 'sub-RESP0621', 'sub-RESP0699'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031153','run-051138','run-011714','run-021549','run-031740',...
    'run-021358','run-021050','run-021057','run-021147','run-031717'};
cfg.ERpath = '/Fridge/users/dorien/derivatives/BB_article/CCEPderiv';

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(cfg);

%% preprocessing CCEP in ECoG

% % sort stimulation pairs
% cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
% cfg.amp = 'no'; % if you want to take stimulation current into account
%
% % select epochs and average
% cfg.epoch_length = 4; % in seconds, -2:2
% cfg.epoch_prestim = 2;
%
% dataBase = preprocess_ECoG_ccep(dataBase,cfg);
%
% disp('All ECoGs are preprocessed')

%% determine channels with IEDs

for subj = 1:size(dataBase,2)
    switch dataBase(subj).sub_label
        case 'sub-RESP0401'
            IEDchan = {'soc4','soc5','soc6','sta3','sta4','stv5','stv6','stv7'};% de kanalen met spikes soc4-6, sta3,4 stv5-7
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0435' % not really convincing spikes...
            IEDchan = {'fh13','fh14','fh15','fh16','fh20','fh21','fh22','fh23',...
                'fh24','fh26','fh27','fh28','fh30','fh31','fl03','fl04','fl05','fl12','fl13'};% FH13-16, 20-24, 26-28,30,31, FL3-5,12,13
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0458'
            IEDchan = {'f08','f09','f14','f17','f18','f19'};% F8-9,14,17-19
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0478'
            IEDchan = {'ihh04','ihh05','ihl01','ihl02','c05','c06','c12','c20'};% IHH4,5, IHL1,2, C5,6,12,20
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0502' %% difficult to differentiate from mu
            IEDchan = {'t19','t20','t21','t22','t23','t28','t29','t30','t31'};% T19-23, 28-31
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0574'
            IEDchan = {'dv1','dv2','dv3','da1','da2','da3','f04','f05','hf10',...
                'hf11','hf12','hf13','hf14'};% Dv1-3, Da1-3, F4,5, HF10-14  CHECK:F4,5 WEL ECHT?
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0589'
            IEDchan = {'stm5','stm6','sta2','sta3','sta4'}; % Stm5,6, Sta2-4 CHECK: STA2,3,4 ECHT?
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0608'
            IEDchan = {'d1','d2','d3','c06','c07','c08','c15','c16','c21','c22'};  % D1-3, C6-8, 15, 16, 21, 22
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0621'
            IEDchan = {'c43','c44','c45','c53'};  % C43,44,45,53   NOG CHECKEN MET FRANS
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0699'
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
            % determine uV 10 ms pre-stim and 10ms post-stim for each channel (not stimulated)
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
    
    fprintf('...Subject %s has been run...\n',dataBase(subj).sub_label)
    
end

fprintf('...All subjects has been run...\n')

%% figure with 15s with 1 channel

subj = 1;
chan = 1;
eventnum = find(strcmp(dataBase(subj).tb_events.sub_type,'SPES')==1,1);
eventsamp = dataBase(subj).tb_events.sample_start(eventnum);

% total signal in time(s)
t = 1/dataBase(subj).ccep_header.Fs:1/dataBase(subj).ccep_header.Fs:size(dataBase(subj).data_noStimArt,2)/dataBase(subj).ccep_header.Fs;
% 15 seconds aruond eventnum
tt = eventsamp/dataBase(subj).ccep_header.Fs-5:1/dataBase(subj).ccep_header.Fs:eventsamp/dataBase(subj).ccep_header.Fs+10;

figure(1),
plot(t,dataBase(subj).data(chan,:),'k')
hold on
plot(t,dataBase(subj).data_noStimArt(chan,:),'r')
hold off
xlim([round(tt(1)), round(tt(end))])
ylim([-3000 3000])
title(sprintf('%s: electrode %s, stimulating %s',dataBase(subj).sub_label,dataBase(subj).ch{chan},dataBase(subj).tb_events.electrical_stimulation_site{eventnum}))

%% visually score spikes in first 10 minutes
subj=1;
fs = dataBase(subj).ccep_header.Fs;

win_start = 0:15:10*60-15;
win_stop = 15:15:10*60;

for t=1:size(win_start,2)
    figure(1)%('units','normalized','position',[0.01 0.01 0.9 0.9]);
    
    hold on
    for chan = 1%:size(dataBase(subj).IEDch,2)
        plot(1/fs:1/fs:size(dataBase(subj).data_noStimArt,2)/fs,dataBase(subj).data_noStimArt(dataBase(subj).IEDch(chan),:))
    end
    hold off
    xlim([win_start(t),win_stop(t)])
    
    set(gcf,'units','normalized','OuterPosition',[0.01 0.01 0.9 0.9],...
        'InnerPosition',[0.05 0.05 0.85 0.85]);
    pause
    
end

%% detect spikes - only in IED channels

for subj = 1:size(dataBase,2)
    data_IED = dataBase(subj).data_noStimArt(dataBase(subj).IEDch,:);
    
    dataBase(subj).data_IED = data_IED;
    
    [M,Pharmat, Pharmatall]=findCortSpikes(dataBase(subj).data_IED,dataBase(subj).ccep_header);
    
    dataBase(subj).spikes.M = M;
    dataBase(subj).spikes.Pharmat = Pharmat;
    dataBase(subj).spikes.Pharmatall = Pharmatall;
    
    filename = ['/Fridge/users/dorien/derivatives/BB_article/CCEPderiv/' dataBase(subj).sub_label, '/',dataBase(subj).ses_label, '/',...
        dataBase(subj).sub_label,'_', dataBase(subj).ses_label, '_', dataBase(subj).task_label,'_', dataBase(subj).run_label,'_Mvalues.mat'];
    
    spikespat.sub_label = dataBase(subj).sub_label;
    spikespat.ses_label = dataBase(subj).ses_label;
    spikespat.run_label = dataBase(subj).run_label;
    spikespat.task_label = dataBase(subj).task_label;
    spikespat.spikesdet = dataBase(subj).spikes;
    spikespat.IEDchan = dataBase(subj).IEDchan;
    spikespat.IEDch = dataBase(subj).IEDch;
    
    save(filename,'spikespat')
    
    
end


%% check spikes visually in the first 10 minutes

subj = 1;

spesstim = strcmp(dataBase(subj).tb_events.sub_type,'SPES');
samp_start = dataBase(subj).tb_events.sample_start(spesstim);
M = dataBase(subj).spikes;
dataBase(subj).data_IED = dataBase(subj).data(dataBase(subj).IEDch,:);

for elec = 1:size(dataBase(subj).IEDch,2)
    
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
        
        % plot detected spikes
        %         locSpikes = [M(elec).groups{6,:}];
        %         valSpikes = [M(elec).groups{7,:}];
        %
        %         locnum = find(locSpikes >tt(1) & locSpikes<tt(2));
        %
        %         plot(locSpikes(locnum)/ccep_header.Fs,...
        %             dataBase(subj).data_IED(elec,locSpikes(locnum)),'*b')
        %         hold off
        
        ax = gca;
        h = gcf;
        h.PaperPosition = [-5 -5 20 10];
        ax.Position = [0.02 0.02 0.95 0.95];
        ylim([-7000 7000])
        xlim([(win-1)*15, win*15])
        title(sprintf('%s: electrode %s',sub_labels{subj},dataBase(subj).IEDchan{elec}))
        grid minor
        
        [x,y]=ginput;
        xall = [xall; round(x*ccep_header.Fs)];
        yall = [yall; y];
        
        %         str = input('Select other point in figure[y/n]: ','s');
        %         n=1;
        %         while strcmp(str,'n')
        %             str = input('Select other point in figure[y/n]: ','s');
        %             dcm = datacursormode;
        %             dcm.DisplayStyle = 'datatip';
        %             dcm.SnaptoDataVertex = 'on';
        %             s(n) = getCursorInfo(dcm);
        %             n=n+1;
        %         end
        
    end
    
    visSpikes(elec).x = xall;
    visSpikes(elec).y = yall;
    
end


%% determine performance by FP, TP, FN





%% put location spikes to stimulation pair and save






