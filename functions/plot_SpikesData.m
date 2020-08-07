
function plot_SpikesData(dataBase,subj,IEDch, IED)

string = [repmat('%s, ',1,size(IEDch,2)-1), '%s'];
visChan = input(sprintf(['Write down which channels you would like to visualize [ch1, ch2]? Choose from ',string,': \n'],IEDch{:}),'s');
if ~isempty(visChan)
    visChansplit = strsplit(visChan,{', ',','});
    channum = NaN(1,size(visChansplit,2));
    for i=1:size(visChansplit,2)
        dig = regexp(visChansplit{i},'[0-9]');
        chars = regexp(lower(visChansplit{i}),'[a-z]');
        test1 = [visChansplit{i}(chars),visChansplit{i}(dig)];
        test2 = [visChansplit{i}(chars),'0', visChansplit{i}(dig)];
        if sum(strcmpi(test1,IEDch))
            channum(i) = find(strcmpi(IEDch,test1)==1);
        elseif sum(strcmpi(test2,IEDch))
            channum(i) = find(strcmpi(IEDch,test2)==1);
        end
    end
    
else
    channum = 1:size(IEDch,2);
end

numIED = size(channum,2);
runs = ceil(numIED/9);

for run=1:runs
    
    if runs==1
        numplots = size(channum,2);
        IEDchan = IEDch(channum);
        IEDnum = IED(channum);
    elseif runs>1 && run~=runs
        numplots = 9;
        IEDchan = IEDch(numplots*(run-1)+1:numplots*run);
        IEDnum = IED(channum(numplots*(run-1)+1:numplots*run));
    elseif runs>1 && run==runs
        numplots=9;
        IEDchan = IEDch(numplots*(run-1)+1:end);
        IEDnum = IED(subj).visIEDs.IED(channum(numplots*(run-1)+1:end));
        numplots = size(IEDnum,2);
    end
    
    win_start = 0:15:10*60-15;
    win_stop = 15:15:10*60;
    htsub = 0.9/numplots; %0.11%0.9/numplots;%0.07;%0.9/numplots;
    startsub = 0.98-htsub:-htsub:0.05;%0.95-htsub:-htsub:0.05; %0.8492:-0.1055:0.1103;%0.05:0.11:0.95;%0.05:htsub:0.95;%0.95-htsub:-htsub:0.05;
        
    figure
    h = gcf;
    h.Units = 'centimeters';
    h.Position(4); % height of figure

    sizesubplot = h.Position(4)*(htsub-0.02) ; %height of each subplot

    fs = dataBase(subj).ccep_header.Fs;   
    yscale = 1600 * sizesubplot; %uV/cm, same as how you would view in Micromed
    
    % downsample signal to speed up visualisation
    downsamp = 16;
    data_down_temp = downsample(dataBase(subj).data_rerefnoStimArt',downsamp);
    data_down = data_down_temp';
    fs_down = round(fs/downsamp);
    
    time_down = 1/fs_down:1/fs_down:size(data_down,2)/fs_down;
    
    sampnum = cell(1,size(channum,2));
    for chan = channum
        for num = 1:size(dataBase(subj).detIED{chan},2)
            if ~isnan(dataBase(subj).detIED{chan}(num))
                sampnum{chan}(num) = find(time_down>=dataBase(subj).detIED{chan}(num)/fs,1,'first');
            else
                sampnum{chan}(num) = NaN;
            end
            
        end
        sampnum{chan} = sampnum{chan}(~isnan(sampnum{chan}));
    end
    
    for t = 1 :size(win_start,2)
        
        yaxis = yscale;
        ymin = -1*yaxis/2;
        ymax = yaxis/2;
        
        for n = 1:numplots
            subplot(numplots,1,n)
            plot(time_down,data_down(IEDnum(n),:))
            hold on
            plot(time_down(sampnum{channum(n)}),data_down(IEDnum(n),sampnum{channum(n)}),'ro')
            
            xlim([win_start(t),win_stop(t)])
            ylim([ymin,ymax])
            ylabel(IEDchan{n});
            ax1(n) = gca;  %#ok<AGROW>
            ax1(n).Position = [0.05 startsub(n) 0.93 htsub-0.02];  %#ok<AGROW>
        end
        ax1(1).Title.String = (sprintf('IED channels of %s, scale = %d %sV/cm',dataBase(subj).sub_label,yscale,char(181)));
        ax1(numplots).XLabel.String= 'Time (s)';
        pause
    end
end