% Figure 5: spikes in epochs (not) connected and (no) suppressed power
% we visualize the number of spikes in each epoch of 200ms 2s prior and 2s
% post-stimulation for several conditions: in response electrodes (not)
% connected to the stimulus pair and in which (no) suppressed power is
% observed. 

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

%% general parameters

pFDR = 0.05;
cmap = parula(101); % 0-1 with steps of 0.01

%% combine ERSPs, CCEPs (only include the IED channels), 
% and spikes of all subjects

% pre-allocation
all_ERSPmat = [];
all_CCEPmat = [];
all_spikessamp = [];
all_metaFile = [];

for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)

        ERSPmat = dataBase(subj).ERSPmat(:,dataBase(subj).IEDs.IEDch);
        all_ERSPmat = [all_ERSPmat; ERSPmat(:)]; % with ERSPmat(:), all columns will be put below each other, so first column 1, then below this column 2 etc.
        
        CCEPmat =  dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        all_CCEPmat = [all_CCEPmat; CCEPmat(:)]; 
        
        % spikessamp had a size [stimp x chan x trials], and we want to
        % reshape it into [stimp*chan x trials]. 
        spikessamp = dataBase(subj).IEDs.totEpochspikessamp;
        [m,n,o] = size(spikessamp);
        spikessamp_reshape = reshape(spikessamp,m*n,o); % for channel 1, all stimps are mentioned, then channel 2 etc. 
        all_spikessamp = [all_spikessamp; spikessamp_reshape];
        
        % this metaFile contains 4 columns: 
        % 1: subject number, 2&3: stimulation pair, 4: channel
        metaFile = [subj*ones(size(ERSPmat,1)*size(ERSPmat,2),1),...
            repmat(dataBase(subj).cc_stimsets,size(ERSPmat,2),1),...
            reshape(repmat(dataBase(subj).IEDs.IEDch',size(ERSPmat,1),1),m*n,1)];
        all_metaFile = [all_metaFile; metaFile]; 
    end
end

% housekeeping
clear CCEPmat ERSPmat m n o spikessamp spikessamp_reshape subj metaFile

%% one overall analysis
spikessampCCEPERSP = []; % CCEP and ERSP
spikessampCCEPnERSP = []; % CCEP and no ERSP
spikessampnCCEPERSP = []; % no CCEP but ERSP
spikessampnCCEPnERSP = []; % no CCEP and no ERSP

for n = 1:size(all_CCEPmat,1)
    if all_CCEPmat(n) == 1 && all_ERSPmat(n) == 1
        spikessampCCEPERSP = [spikessampCCEPERSP; all_spikessamp(n,:)]; 

    elseif all_CCEPmat(n) ==1 && all_ERSPmat(n) == 0
        spikessampCCEPnERSP = [spikessampCCEPnERSP; all_spikessamp(n,:)]; 

    elseif all_CCEPmat(n) == 0 && all_ERSPmat(n) == 1
        spikessampnCCEPERSP = [spikessampnCCEPERSP; all_spikessamp(n,:)]; 

    elseif all_CCEPmat(n) == 0 && all_ERSPmat(n) == 0
        spikessampnCCEPnERSP = [spikessampnCCEPnERSP; all_spikessamp(n,:)]; 

    end
end

%% compare number of spikes per time window for CCEPERSP, CCEPnERSP, nCCEPERSP and nCCEPnERSP channels
tstep = 0.2;

spikessamp_all{1} = spikessampCCEPERSP; % CCEP & ERSP
spikessamp_all{2} = spikessampCCEPnERSP; % CCEP & no ERSP
spikessamp_all{3} = spikessampnCCEPERSP; % no CCEP & ERSP
spikessamp_all{4} = spikessampnCCEPnERSP; % no CCEP & no ERSP
spikessamp_all{5} = [spikessampCCEPERSP; spikessampCCEPnERSP]; % CCEP
spikessamp_all{6} = [spikessampnCCEPERSP; spikessampnCCEPnERSP]; % no CCEP
spikessamp_all{7} = [spikessampCCEPERSP; spikessampnCCEPERSP]; % ERSP
spikessamp_all{8} = [spikessampCCEPnERSP; spikessampnCCEPnERSP]; % no ERSP
spikessamp_all{9} = [spikessampCCEPERSP; spikessampCCEPnERSP; spikessampnCCEPERSP; spikessampnCCEPnERSP]; % all

count_all = cell(size(spikessamp_all));
for m = 1:size(spikessamp_all,2)

    spikessamp = spikessamp_all{m};

    % pre-allocation: count of spikes
    count = NaN(size(spikessamp,1),cfg.epoch_length/0.2);

    % count number of spikes per time window (tstep)
    for num = 1:size(spikessamp,1)
        [n,edges ] = histcounts(vertcat(spikessamp{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

        count(num,1:size(n,2)) = n/10;
    end

count_all{m} = count;

end

%% statistical analysis
% convert edges to mean between two edges
edges_mid =  NaN(1,size(edges,2)-1);
for nn = 1:size(edges,2)-1
    edges_mid(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); 
end

n_poststim = find(round(t(edges_mid),1)>tstep,1,'first'); 
p = NaN(size(spikessamp_all,2),size(count_all{1},2)-n_poststim+1);

for m = 1:size(spikessamp_all,2)
    for n = n_poststim:size(count_all{1},2)

        value_pre = mean(count_all{m}(:,round(t(edges_mid),1)<-tstep),2);
        [~, p(m,n-n_poststim+1)] = ttest(count_all{m}(:,n),value_pre(:));
    end
end

% FDR correction 
m = length(p(:));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*pFDR;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

%% figure

% use edges and plot counts as block signal
edges_disc = NaN(2,size(edges,2)-1); 

for m = 1:size(spikessamp_all,2)

    count = count_all{m};
    
    % use edges and plot counts as block signal
    edges_disc = NaN(2,size(edges,2)-1);
    count_disc = NaN([size(count,1),2,size(edges,2)-1]);

    for nn = 1:size(edges,2)-1
        edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1];
        count_disc(:,:,nn) = repmat(count(:,nn),1,2);
    end

    edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1));
    count_disc = reshape(count_disc, ...
        size(count_disc,1),2*size(count,2));

    % number of IED
    avNumIED_temp = mean(count);
    avNumIED_pre = mean(avNumIED_temp(round(t(edges_mid),1)<-tstep));
    avNumIED = mean(count_disc);
    sterrNumIED = std(count_disc)./sqrt(size(count_disc,1));
    lower_err = avNumIED - sterrNumIED;
    upper_err = avNumIED + sterrNumIED;

    % plot figures
    ymin = 0.1;
    ymax = 0.3;

    figure(m),
    fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))], ...
        [upper_err(t(edges_disc)<-tstep), flip(lower_err(t(edges_disc)<-tstep))], ...
        cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
    hold on
    fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
        [upper_err(t(edges_disc)>tstep) flip(lower_err(t(edges_disc)>tstep))], ...
        cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
    plot(t(edges_disc),avNumIED_pre*ones(size(edges_disc)),'-.','Color',[200/256, 200/256, 200/256])
    plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIED(t(edges_disc)<-tstep),'Color',cmap(10,:),'LineWidth',2);
    plot(t(edges_disc(t(edges_disc)>tstep)),avNumIED(t(edges_disc)>tstep),'Color',cmap(10,:),'LineWidth',2)

    ylim([ymin, ymax])

    % add p<0.05 with FDR correction
    y = 0.95*ymax;
    for n = 1:size(p_sig,2) % number of time windows (9)
        if p_sig(m,n) == 1 % significant
            if p(m,n)<0.001
                text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(10,:));
            elseif p(m,n) <0.01
                text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(10,:));
            elseif p(m,n) < 0.05
                text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(10,:));
            end
        end
    end

    if m == 1
        title('Number of Spikes in electrodes with CCEP and with SP')
        str = 'CCEPERSP';
    elseif m == 2
        title('Number of Spikes in electrodes with CCEP and without SP')
        str = 'CCEPnERSP';
    elseif m == 3
        title('Number of Spikes in electrodes without CCEP and with SP')
        str = 'nCCEPERSP';
    elseif m == 4
        title('Number of Spikes in electrodes without CCEP and without SP')
        str = 'nCCEPnERSP';
    elseif m == 5
        title('Number of Spikes in electrodes with CCEP')
        str = 'CCEP';
    elseif m == 6
        title('Number of Spikes in electrodes without CCEP')
        str = 'nCCEP';
    elseif m == 7
        title('Number of Spikes in electrodes with SP')
        str = 'ERSP';
    elseif m == 8
        title('Number of Spikes in electrodes without SP')
        str = 'nERSP';
    elseif m == 9
        title('Number of Spikes in all electrodes')
        str = 'all';
    end

    xlabel('Time (s)')
    ylabel('Mean spike rate')

    % save figure
    figureName = sprintf('%s/fig6_IEDs%s',...
        myDataPath.figures,str);

    set(gcf,'PaperPositionMode','auto','Renderer','painters')
    print('-dpng','-r300',figureName)
    print('-vector','-depsc',figureName)

    fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)
end

%% end of script