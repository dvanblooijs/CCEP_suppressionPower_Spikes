%% visually score spikes
close all

numplots = size(dataBase(subj).IEDch,2);
htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;
win_start = 0:15:10*60-15;
win_stop = 15:15:10*60;
fs = dataBase(subj).ccep_header.Fs;
time = 1/fs:1/fs:size(dataBase(subj).data_noStimArt,2)/fs;

yscale = 1000; %uV/cm, same as how you would view in Micromed

% downsample signal to speed up visualisation
downsamp = 16;
data_down_temp = downsample(dataBase(subj).data_noStimArt',downsamp);
data_down = data_down_temp';
fs_down = round(fs/downsamp);

time_down = 1/fs_down:1/fs_down:size(data_down,2)/fs_down;
vis_spikes = struct;
x_all = [];
y_all = [];

figure%('units','normalized','position',[0.01 0.01 0.9 0.9]);
disp('Maximize figure and continue by clicking any button')
pause

h = gcf;
h.Units = 'centimeters';
h.Position(4); % height of figure
sizesubplot = h.Position(4)*(htsub-0.02) ; %height of each subplot

for t=1 :3%:size(win_start,2)
    
    yaxis = yscale* sizesubplot;
    ymin = -1*yaxis/2;
    ymax = yaxis/2;
    
    for chan=1:numplots
        subplot(numplots,1,chan)
        plot(time_down,data_down(dataBase(subj).IEDch(chan),:))
        xlim([win_start(t),win_stop(t)])
        ylim([ymin,ymax])
        ylabel(dataBase(subj).IEDchan{chan});
        ax1(chan) = gca;
        ax1(chan).Position = [0.05 startsub(chan) 0.93 htsub-0.02];%[0.05 startsub(chan) 0.9 htsub];
    
        if exist('vis_spikes_all','var')
            hold on
            for num = 1:size(vis_spikes_all(chan).xall,2)
                sampnum = find(time_down>=vis_spikes_all(chan).xall(num),1,'first');
                plot(time_down(sampnum),data_down(dataBase(subj).IEDch(chan),sampnum),'ro')
            end
        end
    
    end
    ax1(1).Title.String = (sprintf('IED channels of %s, scale = %d %sV/cm',dataBase(subj).sub_label,yscale,char(181)));
    ax1(numplots).XLabel.String= 'Time (s)';
    
    currkey = 0;
    n=1;
    x=zeros(size(dataBase(subj).IEDch,2),1); y=zeros(size(dataBase(subj).IEDch,2),1);
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
            [~,locs] = findpeaks(data_down(dataBase(subj).IEDch(elec),...
                sampnum-round(0.1*fs_down):sampnum+round(0.1*fs_down)),'NPeaks',1,'SortStr','descend');
            % find x-position of nearby peak
            locsamp = sampnum-round(0.1*fs_down)+locs-1;
            
            if exist('vis_spikes_all','var')
                if ismembertol(time_down(locsamp),vis_spikes_all(elec).xall,0.1,'DataScale',1)
                    [~,locB] = ismembertol(time_down(locsamp),vis_spikes_all(elec).xall,0.1,'DataScale',1)
                    plot(x(elec,locB),y(elec,locB),'r*'); drawnow;
                    x(elec,locB) = 0;
                    y(elec,locB) = 0;
                    
                    disp('selected spike twice')
                else
                    x(elec,n) = time_down(locsamp);
                    y(elec,n) = data_down(dataBase(subj).IEDch(elec),locsamp);
                    hold on
                    plot(x(elec,n),y(elec,n),'ro'); drawnow;
                    n=n+1;
                end
            else
                x(elec,n) = time_down(locsamp);
                y(elec,n) = data_down(dataBase(subj).IEDch(elec),locsamp);
                hold on
                plot(x(elec,n),y(elec,n),'ro'); drawnow;
                n=n+1;
            end
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

%% put all annotated spikes into one array for each electrode

vis_spikes_all=struct;
% for each electrode in IEDchan
for elec = 1:size(dataBase(subj).IEDch,2)
    xall_indiv = [];
    % for each time window
    for t = 1:size(vis_spikes,2)
        xindiv = vis_spikes(t).x(elec,vis_spikes(t).x(elec,:)~=0);
        xall_indiv = [xall_indiv, xindiv];
    end
    vis_spikes_all(elec).xall = xall_indiv;
end

%% find peak in 0.2s around the selected point - not necessary since I added this in spike rating

% for elec = 1:size(dataBase(subj).IEDch,2)
%    for num = 1:size(vis_spikes_all(elec).xall,2 )
%        
%        samp = find(time>=vis_spikes_all(elec).xall(num),1,'first');
%        
%        if  samp < round(0.1*fs)
%            start_samp = 1;
%            stop_samp = samp + round(0.1*fs);
%        else
%            start_samp = samp - round(0.1*fs);
%            stop_samp = samp + round(0.1*fs);
%        end
%        
%        % find largest peak in a 0.2s time window
%        [pks,locs] = findpeaks(dataBase(subj).data_noStimArt(dataBase(subj).IEDch(elec),...
%            start_samp:stop_samp),'NPeaks',1,'SortStr','descend');
%        
%        loc_time = locs +start_samp -1;
%        
%        vis_spikes_all(elec).xallspike(num) = time(loc_time);
%    end    
% end

%% plot back these peaks and see whether all spikes are correct
close all

numplots = size(dataBase(subj).IEDch,2);
htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;
win_start = 0:15:10*60-15;
win_stop = 15:15:10*60;
fs = dataBase(subj).ccep_header.Fs;
time = 1/fs:1/fs:size(dataBase(subj).data_noStimArt,2)/fs;

yscale = 1000; %uV/cm, same as how you would view in Micromed

% downsample signal to speed up visualisation
downsamp = 16;
data_down_temp = downsample(dataBase(subj).data_noStimArt',downsamp);
data_down = data_down_temp';
fs_down = round(fs/downsamp);

time_down = 1/fs_down:1/fs_down:size(data_down,2)/fs_down;
vis_spikes = struct;

figure%('units','normalized','position',[0.01 0.01 0.9 0.9]);
disp('Maximize figure and continue by clicking any button')
pause

h = gcf;
h.Units = 'centimeters';
h.Position(4); % height of figure
sizesubplot = h.Position(4)*(htsub-0.02) ; %height of each subplot

for t=1 :3%:size(win_start,2)
    
    yaxis = yscale* sizesubplot;
    ymin = -1*yaxis/2;
    ymax = yaxis/2;
    
    for chan=1:numplots
        subplot(numplots,1,chan)
        plot(time_down,data_down(dataBase(subj).IEDch(chan),:)),
        hold on
        for num = 1:size(vis_spikes_all(chan).xallspike,2)
            sampnum = find(time_down>=vis_spikes_all(chan).xallspike(num),1,'first');
            plot(time_down(sampnum),data_down(dataBase(subj).IEDch(chan),sampnum),'ro')
        end
        
        xlim([win_start(t),win_stop(t)])
        ylim([ymin,ymax])
        ylabel(dataBase(subj).IEDchan{chan});
        ax1(chan) = gca;
        ax1(chan).Position = [0.05 startsub(chan) 0.93 htsub-0.02];%[0.05 startsub(chan) 0.9 htsub];
    end
    ax1(1).Title.String = (sprintf('IED channels of %s, scale = %d %sV/cm',dataBase(subj).sub_label,yscale,char(181)));
    ax1(numplots).XLabel.String= 'Time (s)';
    
    currkey = 0;
    n=1;
    x=zeros(size(dataBase(subj).IEDch,2),1); y=zeros(size(dataBase(subj).IEDch,2),1);
    while ~strcmp(currkey,'c')
        w = waitforbuttonpress;
        disp('Click on figure or press "c" to continue')
        if w == 0     %0 for mouse click, 1 for button press
            
            cp = get(gca,'CurrentPoint')
            axcur = gca;
            % find in which subplot spikes were pointed
            elec = find(startsub == axcur.Position(2));
%             elec=1;
            % find sample number closest to the selected point
            [~,sampnum] = min(abs(time_down-cp(1,1)));
            x(elec,n) = time_down(sampnum);
            y(elec,n) = data_down(dataBase(subj).IEDch(elec),sampnum);
            
            hold on
            plot(x(elec,n),y(elec,n),'ro'); drawnow;
            n=n+1;
        elseif w==1
            currkey = get(gcf,'CurrentCharacter');
            
        end
        hold off 
    end
    vis_spikes(t).x = x;
    vis_spikes(t).y = y;
    
end


%% test, en werkt!! :D

numplots = 3;
htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;


fs = 10;
t = 1/fs:1/fs:10;
sign(1,:) = rand(size(t,2),1);
sign(2,:) = rand(size(t,2),1);
sign(3,:) = rand(size(t,2),1);

figure,
ax(1)=subplot(3,1,1);
ax(1).Position = [0.05 startsub(1) 0.93 htsub-0.04];
plot(t,sign(1,:));
ax(2)=subplot(3,1,2);
ax(2).Position = [0.05 startsub(2) 0.93 htsub-0.04];
plot(t,sign(2,:));
ax(3)=subplot(3,1,3);
ax(3).Position = [0.05 startsub(3) 0.93 htsub-0.04];
plot(t,sign(3,:));
currkey = 0;
n=1;
x=[]; y=[];
while ~strcmp(currkey,'c')
    w = waitforbuttonpress;
    if w == 0     %0 for mouse click, 1 for button press
        disp('Button click')
        
        cp = get(gca,'CurrentPoint');
        axcur = gca;
        % find in which subplot spikes were pointed
        elec = find(startsub == axcur.Position(2));
        
        [~,sampnum] = min(abs(t-cp(1,1)));
        x(elec,n) = t(sampnum)
        y(elec,n) = sign(elec,sampnum)
        
        hold on
        plot(x(elec,n),y(elec,n),'ro'); drawnow;
        n=n+1;
    elseif w==1
        currkey = get(gcf,'CurrentCharacter');
        
    end
    hold off
end