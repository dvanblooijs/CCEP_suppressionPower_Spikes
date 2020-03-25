%% visually score spikes
% This script is used to annotate interictal epileptic discharges in ECoG
% data

%% config

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath(genpath('git_rep/eeglab/'))
addpath('git_rep/fieldtrip/')
ft_defaults

%% patient characteristics

cfg.ERpath = '/Fridge/users/dorien/derivatives/BB_article/CCEPderiv';
cfg.dataPath = '/Fridge/CCEP';
cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX or REC2StimXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPESclin';
cfg.run_label = {['run-' input('Run [daydayhhminmin]: ','s')]};
cfg.IEDch = input('Channels showing interictal epileptic discharges: ','s');

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(cfg);

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

%% find IED channels

dataBase.IEDchan = strsplit(cfg.IEDch,',');
dataBase.IED = NaN(size(dataBase.IEDchan));
for i=1:size(dataBase.IEDchan,2)
    charIED = regexp(lower(dataBase.IEDchan{i}),'[a-z]');
    numIED = regexp(dataBase.IEDchan{i},'[0-9]');
    
    type1 = [dataBase.IEDchan{i}(charIED),dataBase.IEDchan{i}(numIED)];
    type2 = [dataBase.IEDchan{i}(charIED),'0',dataBase.IEDchan{i}(numIED)];
    if sum(strcmpi(dataBase.ch,type1))
        dataBase.IED(i) = find(strcmpi(dataBase.ch,type1)==1);
    elseif sum(strcmpi(dataBase.ch,type2))
        dataBase.IED(i) = find(strcmpi(dataBase.ch,type2)==1);
    else
        error('Channel %s or %s is not found',type1,type2)
    end
end

%% visually score IEDs
close all
vis_spikes = struct;

runs = ceil(size(dataBase.IED,2)/9);

for run=1:runs
    
    if runs==1
        numplots = size(dataBase(subj).IED,2);
        IEDchan = dataBase(subj).IEDchan;
        IED = dataBase(subj).IED;
    elseif runs>1 && run~=runs
        numplots = 7;
        IEDchan = dataBase(subj).IEDchan(numplots*(run-1)+1:numplots*run);
        IED = dataBase(subj).IED(numplots*(run-1)+1:numplots*run);
    elseif runs>1 && run==runs
        numplots=7;
        IEDchan = dataBase(subj).IEDchan(numplots*(run-1)+1:end);
        IED = dataBase(subj).IED(numplots*(run-1)+1:end);
        numplots = size(IED,2);
    end
   
    htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
    startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;
    win_start = 0:15:10*60-15;
    win_stop = 15:15:10*60;
    fs = dataBase(subj).ccep_header.Fs;
    time = 1/fs:1/fs:size(dataBase(subj).data_noStimArt,2)/fs;
    
    yscale = 2500; %uV/cm, same as how you would view in Micromed
    
    % downsample signal to speed up visualisation
    downsamp = 16;
    data_down_temp = downsample(dataBase(subj).data_noStimArt',downsamp);
    data_down = data_down_temp';
    fs_down = round(fs/downsamp);
    
    time_down = 1/fs_down:1/fs_down:size(data_down,2)/fs_down;
    % vis_spikes = struct;
    x_all = [];
    y_all = [];
    
    figure%('units','normalized','position',[0.01 0.01 0.9 0.9]);
    disp('Maximize figure and continue by clicking any button')
    pause
    
    h = gcf;
    h.Units = 'centimeters';
    h.Position(4); % height of figure
    sizesubplot = h.Position(4)*(htsub-0.02) ; %height of each subplot
    
    for t=1 :size(win_start,2)
        
        yaxis = yscale* sizesubplot;
        ymin = -1*yaxis/2;
        ymax = yaxis/2;
        
        for chan=1:numplots
            subplot(numplots,1,chan)
            plot(time_down,data_down(IED(chan),:))
            hold on
            plot(0.1*fs*ones(round(ymax)-round(ymin),1),round(ymin):round(ymax)-1,'k:')
            xlim([win_start(t),win_stop(t)])
            ylim([ymin,ymax])
            ylabel(IEDchan{chan});
            ax1(chan) = gca;
            ax1(chan).Position = [0.05 startsub(chan) 0.93 htsub-0.02];%[0.05 startsub(chan) 0.9 htsub];
            
            if exist('vis_spikes_all','var')
%                 hold on
                for num = 1:size(vis_spikes_all(chan).xall,2)
                    sampnum = find(time_down>=vis_spikes_all(chan).xall(num),1,'first');
                    plot(time_down(sampnum),data_down(IED(chan),sampnum),'ro')
                end
            end
            
        end
        ax1(1).Title.String = (sprintf('IED channels of %s, scale = %d %sV/cm',dataBase(subj).sub_label,yscale,char(181)));
        ax1(numplots).XLabel.String= 'Time (s)';
        
        currkey = 0;
        n=1;
        x=zeros(size(IED,2),1); y=zeros(size(IED,2),1);
        while ~strcmp(currkey,'c')
            w = waitforbuttonpress;
            disp('Click on figure or press "c" to continue')
            if w == 0     %0 for mouse click, 1 for button press
                
                cp = get(gca,'CurrentPoint');
                axcur = gca;
                % find in which subplot spikes were pointed
                elec = find(startsub == axcur.Position(2));
                %             elec=1;
                % find sample number closest to the selected point
                [~,sampnum] = min(abs(time_down-cp(1,1)));
                %             x(elec,n) = time_down(sampnum);
                %             y(elec,n) = data_down(dataBase(subj).IEDch(elec),sampnum);
                
                % find nearby peak
                [~,locs] = findpeaks(data_down(IED(elec),...
                    sampnum-round(0.1*fs_down):sampnum+round(0.1*fs_down)),'NPeaks',1,'SortStr','descend');
                % find x-position of nearby peak
                locsamp = sampnum-round(0.1*fs_down)+locs-1;
                
                %             if exist('vis_spikes_all','var')
                if ismembertol(time_down(locsamp),x(elec,:),0.1,'DataScale',1)
                    [~,locB] = ismembertol(time_down(locsamp),x(elec,:),0.1,'DataScale',1)
                    hold on
                    plot(x(elec,locB),y(elec,locB),'wo'); drawnow;
                    x(elec,locB) = 0;
                    y(elec,locB) = 0;
                    
                    disp('selected spike twice')
                else
                    x(elec,n) = time_down(locsamp);
                    y(elec,n) = data_down(IED(elec),locsamp);
                    hold on
                    plot(x(elec,n),y(elec,n),'ro'); drawnow;
                    n=n+1;
                end
                %             else
                %                 x(elec,n) = time_down(locsamp);
                %                 y(elec,n) = data_down(dataBase(subj).IEDch(elec),locsamp);
                %                 hold on
                %                 plot(x(elec,n),y(elec,n),'ro'); drawnow;
                %                 n=n+1;
                %             end
            elseif w==1
                currkey = get(gcf,'CurrentCharacter');
                
            end
            hold off
        end
        x_all = [x_all,x];
        y_all = [y_all,y];
        %     vis_spikes(t).x = x;
        %     vis_spikes(t).y = y;
        
    end
    vis_spikes(run).x_all = x_all;
    vis_spikes(run).y_all = y_all;

end

disp('Annotation finished!')

%% reshap vis_spikes

if runs == 3
x_all = [vis_spikes(1).x_all, zeros(size(vis_spikes(1).x_all,1),size(vis_spikes(2).x_all,2)), zeros(size(vis_spikes(1).x_all,1),size(vis_spikes(3).x_all,2));...
    zeros(size(vis_spikes(2).x_all,1),size(vis_spikes(1).x_all,2)), vis_spikes(2).x_all,  zeros(size(vis_spikes(2).x_all,1),size(vis_spikes(3).x_all,2));...
    zeros(size(vis_spikes(3).x_all,1),size(vis_spikes(1).x_all,2)), zeros(size(vis_spikes(3).x_all,1),size(vis_spikes(2).x_all,2)), vis_spikes(3).x_all];

y_all = [vis_spikes(1).y_all, zeros(size(vis_spikes(1).y_all,1),size(vis_spikes(2).y_all,2)), zeros(size(vis_spikes(1).y_all,1),size(vis_spikes(3).y_all,2));...
    zeros(size(vis_spikes(2).y_all,1),size(vis_spikes(1).y_all,2)), vis_spikes(2).y_all,  zeros(size(vis_spikes(2).y_all,1),size(vis_spikes(3).y_all,2));...
    zeros(size(vis_spikes(3).y_all,1),size(vis_spikes(1).y_all,2)), zeros(size(vis_spikes(3).y_all,1),size(vis_spikes(2).y_all,2)), vis_spikes(3).y_all];

elseif runs ==2
    x_all = [vis_spikes(1).x_all, zeros(size(vis_spikes(1).x_all,1),size(vis_spikes(2).x_all,2));...
        zeros(size(vis_spikes(2).x_all,1),size(vis_spikes(1).x_all,2)), vis_spikes(2).x_all];
    
    y_all = [vis_spikes(1).y_all, zeros(size(vis_spikes(1).y_all,1),size(vis_spikes(2).y_all,2));...
        zeros(size(vis_spikes(2).y_all,1),size(vis_spikes(1).y_all,2)), vis_spikes(2).y_all];
    
end
    

%% save all annotated spikes

folderName = '/Fridge/users/dorien/derivatives/BB_article/IEDs/';
fileName = [dataBase(subj).sub_label, '_', dataBase(subj).ses_label, '_',...
    dataBase(subj).task_label, '_', dataBase(subj).run_label, '_visIEDs.mat'];

sub_label = dataBase(subj).sub_label;
ses_label = dataBase(subj).ses_label;
task_label = dataBase(subj).task_label;
run_label = dataBase(subj).run_label;
IED = dataBase(subj).IED;
IEDchan = dataBase(subj).IEDchan;

save([folderName, fileName],'sub_label','ses_label','task_label','run_label',...
    'IED','IEDchan','x_all','y_all')

disp('saved')

% %% put all annotated spikes into one array for each electrode
% 
% vis_spikes_all=struct;
% % for each electrode in IEDchan
% for elec = 1:size(dataBase(subj).IEDch,2)
%     xall_indiv = [];
%     % for each time window
%     for t = 1:size(vis_spikes,2)
%         xindiv = vis_spikes(t).x(elec,vis_spikes(t).x(elec,:)~=0);
%         xall_indiv = [xall_indiv, xindiv];
%     end
%     vis_spikes_all(elec).xall = xall_indiv;
% end
% 
% %% find peak in 0.2s around the selected point - not necessary since I added this in spike rating
% 
% % for elec = 1:size(dataBase(subj).IEDch,2)
% %    for num = 1:size(vis_spikes_all(elec).xall,2 )
% %
% %        samp = find(time>=vis_spikes_all(elec).xall(num),1,'first');
% %
% %        if  samp < round(0.1*fs)
% %            start_samp = 1;
% %            stop_samp = samp + round(0.1*fs);
% %        else
% %            start_samp = samp - round(0.1*fs);
% %            stop_samp = samp + round(0.1*fs);
% %        end
% %
% %        % find largest peak in a 0.2s time window
% %        [pks,locs] = findpeaks(dataBase(subj).data_noStimArt(dataBase(subj).IEDch(elec),...
% %            start_samp:stop_samp),'NPeaks',1,'SortStr','descend');
% %
% %        loc_time = locs +start_samp -1;
% %
% %        vis_spikes_all(elec).xallspike(num) = time(loc_time);
% %    end
% % end
% 
% %% plot back these peaks and see whether all spikes are correct
% close all
% 
% numplots = size(dataBase(subj).IEDch,2);
% htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
% startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;
% win_start = 0:15:10*60-15;
% win_stop = 15:15:10*60;
% fs = dataBase(subj).ccep_header.Fs;
% time = 1/fs:1/fs:size(dataBase(subj).data_noStimArt,2)/fs;
% 
% yscale = 1000; %uV/cm, same as how you would view in Micromed
% 
% % downsample signal to speed up visualisation
% downsamp = 16;
% data_down_temp = downsample(dataBase(subj).data_noStimArt',downsamp);
% data_down = data_down_temp';
% fs_down = round(fs/downsamp);
% 
% time_down = 1/fs_down:1/fs_down:size(data_down,2)/fs_down;
% vis_spikes = struct;
% 
% figure%('units','normalized','position',[0.01 0.01 0.9 0.9]);
% disp('Maximize figure and continue by clicking any button')
% pause
% 
% h = gcf;
% h.Units = 'centimeters';
% h.Position(4); % height of figure
% sizesubplot = h.Position(4)*(htsub-0.02) ; %height of each subplot
% 
% for t=1 :3%:size(win_start,2)
%     
%     yaxis = yscale* sizesubplot;
%     ymin = -1*yaxis/2;
%     ymax = yaxis/2;
%     
%     for chan=1:numplots
%         subplot(numplots,1,chan)
%         plot(time_down,data_down(dataBase(subj).IEDch(chan),:)),
%         hold on
%         for num = 1:size(vis_spikes_all(chan).xallspike,2)
%             sampnum = find(time_down>=vis_spikes_all(chan).xallspike(num),1,'first');
%             plot(time_down(sampnum),data_down(dataBase(subj).IEDch(chan),sampnum),'ro')
%         end
%         
%         xlim([win_start(t),win_stop(t)])
%         ylim([ymin,ymax])
%         ylabel(dataBase(subj).IEDchan{chan});
%         ax1(chan) = gca;
%         ax1(chan).Position = [0.05 startsub(chan) 0.93 htsub-0.02];%[0.05 startsub(chan) 0.9 htsub];
%     end
%     ax1(1).Title.String = (sprintf('IED channels of %s, scale = %d %sV/cm',dataBase(subj).sub_label,yscale,char(181)));
%     ax1(numplots).XLabel.String= 'Time (s)';
%     
%     currkey = 0;
%     n=1;
%     x=zeros(size(dataBase(subj).IEDch,2),1); y=zeros(size(dataBase(subj).IEDch,2),1);
%     while ~strcmp(currkey,'c')
%         w = waitforbuttonpress;
%         disp('Click on figure or press "c" to continue')
%         if w == 0     %0 for mouse click, 1 for button press
%             
%             cp = get(gca,'CurrentPoint')
%             axcur = gca;
%             % find in which subplot spikes were pointed
%             elec = find(startsub == axcur.Position(2));
%             %             elec=1;
%             % find sample number closest to the selected point
%             [~,sampnum] = min(abs(time_down-cp(1,1)));
%             x(elec,n) = time_down(sampnum);
%             y(elec,n) = data_down(dataBase(subj).IEDch(elec),sampnum);
%             
%             hold on
%             plot(x(elec,n),y(elec,n),'ro'); drawnow;
%             n=n+1;
%         elseif w==1
%             currkey = get(gcf,'CurrentCharacter');
%             
%         end
%         hold off
%     end
%     vis_spikes(t).x = x;
%     vis_spikes(t).y = y;
%     
% end
% 
% 
% %% test, en werkt!! :D
% 
% numplots = 3;
% htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
% startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;
% 
% 
% fs = 10;
% t = 1/fs:1/fs:10;
% sign(1,:) = rand(size(t,2),1);
% sign(2,:) = rand(size(t,2),1);
% sign(3,:) = rand(size(t,2),1);
% 
% figure,
% ax(1)=subplot(3,1,1);
% ax(1).Position = [0.05 startsub(1) 0.93 htsub-0.04];
% plot(t,sign(1,:));
% ax(2)=subplot(3,1,2);
% ax(2).Position = [0.05 startsub(2) 0.93 htsub-0.04];
% plot(t,sign(2,:));
% ax(3)=subplot(3,1,3);
% ax(3).Position = [0.05 startsub(3) 0.93 htsub-0.04];
% plot(t,sign(3,:));
% currkey = 0;
% n=1;
% x=[]; y=[];
% while ~strcmp(currkey,'c')
%     w = waitforbuttonpress;
%     if w == 0     %0 for mouse click, 1 for button press
%         disp('Button click')
%         
%         cp = get(gca,'CurrentPoint');
%         axcur = gca;
%         % find in which subplot spikes were pointed
%         elec = find(startsub == axcur.Position(2));
%         
%         [~,sampnum] = min(abs(t-cp(1,1)));
%         x(elec,n) = t(sampnum)
%         y(elec,n) = sign(elec,sampnum)
%         
%         hold on
%         plot(x(elec,n),y(elec,n),'ro'); drawnow;
%         n=n+1;
%     elseif w==1
%         currkey = get(gcf,'CurrentCharacter');
%         
%     end
%     hold off
% end