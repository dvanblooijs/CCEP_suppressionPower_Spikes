function vis_spikes = visualRateSpikes(dataBase)

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
    
    yscale = 1600; %uV/cm, same as how you would view in Micromed
    
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
    
    for t = 1:size(win_start,2)
        
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
            ax1(chan) = gca;  %#ok<AGROW>
            ax1(chan).Position = [0.05 startsub(chan) 0.93 htsub-0.02];%#ok<AGROW> %[0.05 startsub(chan) 0.9 htsub];
            
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
                    [~,locB] = ismembertol(time_down(locsamp),x(elec,:),0.1,'DataScale',1) %#ok<NOPRT>
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
            elseif w==1
                currkey = get(gcf,'CurrentCharacter');
                
            end
            hold off
        end
        x_all = [x_all,x]; %#ok<AGROW>
        y_all = [y_all,y]; %#ok<AGROW>
        
    end
    vis_spikes(run).x_all = x_all;
    vis_spikes(run).y_all = y_all;

end

disp('Annotation finished!')
