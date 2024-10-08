% Figure 6: spikes in epochs (not) connected and (no) suppressed power

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

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
        all_metaFile = [all_metaFile; metaFile]; %#ok<AGROW>
    end
end

% housekeeping
% clear CCEPmat ERSPmat m n o spikessamp spikessamp_reshape subj metaFile

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

% convert edges to mean between two edges
edges_mid =  NaN(1,size(edges,2)-1);
for nn = 1:size(edges,2)-1
    edges_mid(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); 
end

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
        myDataPath.Figures,str);

    set(gcf,'PaperPositionMode','auto','Renderer','painters')
    print('-dpng','-r300',figureName)
    print('-vector','-depsc',figureName)

    fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)
end

% %% add all spike samples into two categories: nCCEP or CCEP
% 
% spikessampCCEP = [];
% spikessampnCCEP = [];
% metaFileCCEP = []; metaFilenCCEP = [];
% 
% for n = 1:size(all_CCEPmat,1)
% 
%     if all_CCEPmat(n) == 1
%         spikessampCCEP = [spikessampCCEP; all_spikessamp(n,:)]; 
%         metaFileCCEP = [metaFileCCEP; all_metaFile(n,:)]; 
% 
%     elseif all_CCEPmat(n) == 0
%         spikessampnCCEP = [spikessampnCCEP; all_spikessamp(n,:)]; 
%         metaFilenCCEP = [metaFilenCCEP; all_metaFile(n,:)];
%     end
% end
% 
% % housekeeping
% clear I n
% 
% % compare number of spikes per time window for CCEP and nCCEP channels
% tstep = 0.2; % this is the time window
% 
% % pre-allocation: count of spikes
% countCCEP = NaN(size(spikessampCCEP,1),cfg.epoch_length/tstep); % [[responses with CCEP x epochs of 0.2s]
% countnCCEP = NaN(size(spikessampnCCEP,1),cfg.epoch_length/tstep); % [responses without CCEP x epochs of 0.2s
% 
% % count number of spikes per time window (tstep)
% % in electrodes with CCEPs
% for num = 1:size(spikessampCCEP,1)
%     [n,edges ] = histcounts(vertcat(spikessampCCEP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countCCEP(num,1:size(n,2)) = n/10; % divided by ten, because ten trials
% end
% 
% % count number of spikes per time window (tstep)
% % in electrodes without CCEPs
% for num = 1:size(spikessampnCCEP,1)
%     [n,edges ] = histcounts(vertcat(spikessampnCCEP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countnCCEP(num,1:size(n,2)) = n/10;
% end
% 
% % convert edges to mean between the two edges: 0-410, 410-819,..., 7782-8192
% edges_mid = NaN(1,size(edges,2)-1);
% for nn = 1:size(edges,2)-1
%     edges_mid(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); 
% end
% 
% % use edges and plot counts as block signal
% % pre-allocation
% edges_disc = NaN(2,size(edges,2)-1); countnCCEP_disc = NaN([size(countnCCEP,1),2,size(edges,2)-1]); countCCEP_disc = NaN([size(countCCEP,1),2,size(edges,2)-1]);
% for nn = 1:size(edges,2)-1
%     edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1]; 
%     countnCCEP_disc(:,:,nn) = repmat(countnCCEP(:,nn),1,2); 
%     countCCEP_disc(:,:,nn) = repmat(countCCEP(:,nn),1,2); 
% end
% 
% edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1)); % was [2x20] and becomes [1x40] 1-409-411-818-820 etc. Each time window (1-409, 411-818) becomes a horizontal line
% countnCCEP_disc = reshape(countnCCEP_disc,size(countnCCEP_disc,1),2*size(countnCCEP,2)); 
% countCCEP_disc = reshape(countCCEP_disc,size(countCCEP_disc,1),2*size(countCCEP,2));
% 
% % number of IED when CCEP
% avNumIEDCCEP_temp = mean(countCCEP); % average all counts per time window
% avNumIEDCCEP_pre = mean(avNumIEDCCEP_temp(round(t(edges_mid),1)<-tstep)); % we want the step around 0 to be excluded due to the artefact
% avNumIEDCCEP = mean(countCCEP_disc);
% sterrNumIEDCCEP = std(countCCEP_disc)./sqrt(size(countCCEP_disc,1)); % standard error of the mean
% lower_errCCEP = avNumIEDCCEP - sterrNumIEDCCEP;
% upper_errCCEP = avNumIEDCCEP + sterrNumIEDCCEP;
% 
% % number of IED when no CCEP
% avNumIEDnCCEP_temp = mean(countnCCEP);
% avNumIEDnCCEP_pre = mean(avNumIEDnCCEP_temp(round(t(edges_mid),1)<-tstep)); % we want the step around 0 to be excluded due to the artefact
% avNumIEDnCCEP = mean(countnCCEP_disc);
% sterrNumIEDnCCEP = std(countnCCEP_disc)./sqrt(size(countnCCEP_disc,1));
% lower_errnCCEP = avNumIEDnCCEP - sterrNumIEDnCCEP;
% upper_errnCCEP = avNumIEDnCCEP + sterrNumIEDnCCEP;
% 
% % statistical analysis: mann whitney u test
% n_poststim = find(round(t(edges_mid),1) > tstep,1,'first');
% p = NaN(2,size(countCCEP,2)-n_poststim+1);
% for n = n_poststim:size(countCCEP,2)
%         CCEP_pre = mean(countCCEP(:,round(t(edges_mid),1)<-tstep),2);
%         nCCEP_pre = mean(countnCCEP(:,round(t(edges_mid),1)<-tstep),2);
%         [~, p(1,n-n_poststim+1)] = ttest(countCCEP(:,n),CCEP_pre(:));
%         [~, p(2,n-n_poststim+1)] = ttest(countnCCEP(:,n),nCCEP_pre(:));
% %     CCEP_pre = median(countCCEP(:,round(t(edges_mid),1)<-tstep),2);
% %     nCCEP_pre = median(countnCCEP(:,round(t(edges_mid),1)<-tstep),2);
% %     p(1,n-n_poststim+1) = signrank(countCCEP(:,n),CCEP_pre(:));
% %     p(2,n-n_poststim+1) = signrank(countnCCEP(:,n),nCCEP_pre(:));
% end
% 
% % FDR correction
% m = length(p(:));
% [p_sort, p_ind] = sort(p(:));
% thisVal = NaN(size(p_sort));
% for kk = 1:length(p_sort)
%     thisVal(kk) = (kk/m)*pFDR;
% end
% 
% p_sig = p;
% p_sig(p_ind) = p_sort<thisVal;
% 
% % plot figures
% ymin = round(0.9*min([avNumIEDnCCEP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDCCEP(t(edges_disc)<-tstep| t(edges_disc)>tstep)]),2,'significant');
% ymax = round(1.1*max([avNumIEDnCCEP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDCCEP(t(edges_disc)<-tstep| t(edges_disc)>tstep)]),2,'significant');
% 
% figure,
% subplot(2,1,1);
% % filled box for standard error of number of spikes in responses WITH ccep
% % pre-stim
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))],...
%     [upper_errCCEP(t(edges_disc)<-tstep) flip(lower_errCCEP(t(edges_disc)<-tstep))],...
%     cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
% hold on
% % filled box for standard error of number of spikes in responses WITH ccep
% % post-stim
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errCCEP(t(edges_disc)>tstep) flip(lower_errCCEP(t(edges_disc)>tstep))], ...
%     cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
% % plot average pre, plot pre spikes, post spikes
% plot(t(edges_disc),avNumIEDCCEP_pre*ones(size(edges_disc)),'-.','Color',[200/256, 200/256, 200/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDCCEP(t(edges_disc)<-tstep),'Color',cmap(10,:),'LineWidth',2)
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDCCEP(t(edges_disc)>tstep),'Color',cmap(10,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin,ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n = 1:size(p_sig,2) % number of time windows (9)
%     m = 1; %CCEP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when CCEP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% subplot(2,1,2);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))],...
%     [upper_errnCCEP(t(edges_disc)<-tstep) flip(lower_errnCCEP(t(edges_disc)<-tstep))],...
%     cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errnCCEP(t(edges_disc)>tstep) flip(lower_errnCCEP(t(edges_disc)>tstep))], ...
%     cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
% plot(t(edges_disc),avNumIEDnCCEP_pre*ones(size(edges_disc)),'-.','Color',[180/256, 180/256, 180/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDnCCEP(t(edges_disc)<-tstep),'Color',cmap(20,:),'LineWidth',2)
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDnCCEP(t(edges_disc)>tstep),'Color',cmap(20,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin,ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 2; %nCCEP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% 
% title('Number of Spikes when no CCEP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% h = gcf;
% h.Units = 'normalized';
% h.Position = [0.1 0.4 0.4 0.5];
% 
% figureName = sprintf('%s/fig6_IEDsCCEP',...
%     myDataPath.Figures);
% 
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-vector','-depsc',figureName)
% 
% fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)
% 
% % % housekeeping
% % clear avNumIEDCCEP avNumIEDCCEP_pre avNumIEDCCEP_temp CCEP_pre nCCEP_pre nn
% % clear avNumIEDnCCEP avNumIEDnCCEP_pre avNumIEDnCCEP_temp num p p_ind p_sort p_sig s
% % clear countCCEP countCCEP_disc countCCEP_disc n n_poststim
% % clear countnCCEP countnCCEP_disc countnCCEP_disc m metaFileCCEP metaFilenCCEP
% % clear edges edges_mid edges_disc edges_disc h kk lower_errCCEP lower_errnCCEP
% % clear spikessampCCEP spikessampnCCEP sterrNumIEDCCEP sterrNumIEDnCCEP thisVal upper_errnCCEP
% % clear upper_errCCEP y ymax ymin h figureName ans
% 
% %% add all spike samples into two categories: nERSP or ERSP
% 
% spikessampERSP = []; % ERSP
% spikessampnERSP = []; % no ERSP
% metaFileERSP = []; metaFilenERSP = [];
% 
% for n = 1:size(all_ERSPmat,1)
%     if all_ERSPmat(n) == 1
%         spikessampERSP = [spikessampERSP; all_spikessamp(n,:)]; 
%         metaFileERSP = [metaFileERSP; all_metaFile(n,:)]; 
% 
%     elseif all_ERSPmat(n) == 0
%         spikessampnERSP = [spikessampnERSP; all_spikessamp(n,:)];
%         metaFilenERSP = [metaFilenERSP; all_metaFile(n,:)]; 
%     end
% end
% 
% % housekeeping
% clear I n
% 
% % compare number of spikes per time window for ERSP and nERSP channels
% tstep = 0.2;
% 
% % pre-allocation: count of spikes
% countERSP = NaN(size(spikessampERSP,1),cfg.epoch_length/0.2);
% countnERSP = NaN(size(spikessampnERSP,1),cfg.epoch_length/0.2);
% 
% % count number of spikes per time window (tstep)
% % in electrodes with ERSP
% for num = 1:size(spikessampERSP,1)
%     [n,edges ] = histcounts(vertcat(spikessampERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countERSP(num,1:size(n,2)) = n/10;
% end
% 
% % count number of spikes per time window (tstep)
% % in electrodes with ERSPs
% for num = 1:size(spikessampnERSP,1)
%     [n,edges ] = histcounts(vertcat(spikessampnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countnERSP(num,1:size(n,2)) = n/10;
% end
% 
% % convert edges to mean between the two edges
% edges_mid = NaN(1,size(edges,2)-1);
% for nn = 1:size(edges,2)-1
%     edges_mid(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); 
% end
% 
% % use edges and plot counts as block signal
% edges_disc = NaN(2,size(edges,2)-1); countnERSP_disc = NaN([size(countnERSP,1),2,size(edges,2)-1]); countERSP_disc = NaN([size(countERSP,1),2,size(edges,2)-1]);
% for nn = 1:size(edges,2)-1
%     edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1];
%     countnERSP_disc(:,:,nn) = repmat(countnERSP(:,nn),1,2); 
%     countERSP_disc(:,:,nn) = repmat(countERSP(:,nn),1,2); 
% end
% 
% edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1));
% countnERSP_disc = reshape(countnERSP_disc,size(countnERSP_disc,1),2*size(countnERSP,2));
% countERSP_disc = reshape(countERSP_disc,size(countERSP_disc,1),2*size(countERSP,2));
% 
% % number of IED when ERSP
% % MEAN
% avNumIEDERSP_temp = mean(countERSP);
% avNumIEDERSP_pre = mean(avNumIEDERSP_temp(round(t(edges_mid),1)<-tstep));
% avNumIEDERSP = mean(countERSP_disc);
% sterrNumIEDERSP = std(countERSP_disc)./sqrt(size(countERSP_disc,1));
% lower_errERSP = avNumIEDERSP - sterrNumIEDERSP;
% upper_errERSP = avNumIEDERSP + sterrNumIEDERSP;
% % MEDIAN 
% % avNumIEDERSP_temp = median(countERSP);
% % avNumIEDERSP_pre = median(avNumIEDERSP_temp(round(t(edges_mid),1)<-tstep));
% % avNumIEDERSP = median(countERSP_disc);
% % sterrNumIEDERSP = quantile(countERSP_disc,[.25 .75]);
% % lower_errERSP = avNumIEDERSP - sterrNumIEDERSP(1,:);
% % upper_errERSP = avNumIEDERSP + sterrNumIEDERSP(2,:);
% 
% % number of IED when no ERSP
% % MEAN
% avNumIEDnERSP_temp = mean(countnERSP);
% avNumIEDnERSP_pre = mean(avNumIEDnERSP_temp(round(t(edges_mid),1)<-tstep));
% avNumIEDnERSP = mean(countnERSP_disc);
% sterrNumIEDnERSP = std(countnERSP_disc)./sqrt(size(countnERSP_disc,1));
% lower_errnERSP = avNumIEDnERSP - sterrNumIEDnERSP;
% upper_errnERSP = avNumIEDnERSP + sterrNumIEDnERSP;
% % MEDIAN
% % avNumIEDnERSP_temp = median(countnERSP);
% % avNumIEDnERSP_pre = median(avNumIEDnERSP_temp(round(t(edges_mid),1)<-tstep));
% % avNumIEDnERSP = median(countnERSP_disc);
% % sterrNumIEDnERSP = quantile(countnERSP_disc,[.25 .75]);
% % lower_errnERSP = avNumIEDnERSP - sterrNumIEDnERSP(1,:);
% % upper_errnERSP = avNumIEDnERSP + sterrNumIEDnERSP(2,:);
% 
% % statistical analysis: mann whitney u test
% n_poststim = find(round(t(edges_mid),1)>tstep,1,'first');
% p = NaN(2,size(countERSP,2)-n_poststim+1);
% for n = n_poststim:size(countERSP,2)
%     ERSP_pre = mean(countERSP(:,round(t(edges_mid),1)<-tstep),2);
%     nERSP_pre = mean(countnERSP(:,round(t(edges_mid),1)<-tstep),2);
%    [~, p(1,n-n_poststim+1)] = ttest(countERSP(:,n),ERSP_pre(:)); 
%    [~,p(2,n-n_poststim+1)] = ttest(countnERSP(:,n),nERSP_pre(:));
% %    ERSP_pre = median(countERSP(:,round(t(edges_mid),1)<-tstep),2);
% %     nERSP_pre = median(countnERSP(:,round(t(edges_mid),1)<-tstep),2);
% %     p(1,n-n_poststim+1) = signrank(countERSP(:,n),ERSP_pre(:)); 
% %     p(2,n-n_poststim+1) = signrank(countnERSP(:,n),nERSP_pre(:));
% 
% end
% 
% % FDR correction
% m = length(p(:));
% [p_sort,p_ind] = sort(p(:));
% thisVal = NaN(size(p_sort));
% for kk = 1:length(p_sort)
%     thisVal(kk) = (kk/m)*pFDR;
% end
% 
% p_sig = p;
% p_sig(p_ind) = p_sort<thisVal;
% 
% % plot figures
% ymin = round(0.9*min([avNumIEDnERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep)]),2,'significant');
% ymax = round(1.1*max([avNumIEDnERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep)]),2,'significant');
% 
% figure,
% subplot(2,1,1);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))],...
%     [upper_errERSP(t(edges_disc)<-tstep) flip(lower_errERSP(t(edges_disc)<-tstep))],...
%     cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errERSP(t(edges_disc)>tstep) flip(lower_errERSP(t(edges_disc)>tstep))], ...
%     cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
% plot(t(edges_disc),avNumIEDERSP_pre*ones(size(edges_disc)),'-.','Color',[200/256, 200/256, 200/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDERSP(t(edges_disc)<-tstep),'Color',cmap(10,:),'LineWidth',2)
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDERSP(t(edges_disc)>tstep),'Color',cmap(10,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin,ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 1; %ERSP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when ERSP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% subplot(2,1,2);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))],...
%     [upper_errnERSP(t(edges_disc)<-tstep) flip(lower_errnERSP(t(edges_disc)<-tstep))],...
%     cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errnERSP(t(edges_disc)>tstep) flip(lower_errnERSP(t(edges_disc)>tstep))], ...
%     cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
% plot(t(edges_disc),avNumIEDnERSP_pre*ones(size(edges_disc)),'-.','Color',[180/256, 180/256, 180/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDnERSP(t(edges_disc)<-tstep),'Color',cmap(20,:),'LineWidth',2)
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDnERSP(t(edges_disc)>tstep),'Color',cmap(20,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin,ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 2; %nERSP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when no ERSP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% h = gcf;
% h.Units = 'normalized';
% h.Position = [0.1 -0.2 0.4 0.5];
% 
% figureName = sprintf('%s/fig6_IEDsERSP',...
%     myDataPath.Figures);
% 
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-vector','-depsc',figureName)
% 
% fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)
% 
% % housekeeping
% % clear avNumIEDERSP avNumIEDERSP_pre avNumIEDERSP_temp ERSP_pre nERSP_pre nn
% % clear avNumIEDnERSP avNumIEDnERSP_pre avNumIEDnERSP_temp num p p_ind p_sort p_sig s
% % clear countERSP countERSP_disc countERSP n n_poststim
% % clear countnERSP countnERSP_disc countnERSP m metaFileERSP metaFilenERSP
% % clear edges edges_mid edges_disc edges_disc h kk lower_errERSP lower_errnERSP
% % clear spikessampERSP spikessampnERSP sterrNumIEDERSP sterrNumIEDnERSP thisVal upper_errnERSP
% % clear upper_errERSP y ymax ymin figureName
% 
% %% add all spike samples into four categories: CCEPERSP, CCEPnERSP, nCCEPERSP or nCCEPnERSP
% 
% spikessampCCEPERSP = []; % CCEP and ERSP
% spikessampCCEPnERSP = []; % CCEP and no ERSP
% spikessampnCCEPERSP = []; % no CCEP but ERSP
% spikessampnCCEPnERSP = []; % no CCEP and no ERSP
% metaFileCCEPERSP = []; metaFilenCCEPERSP = [];
% metaFileCCEPnERSP = []; metaFilenCCEPnERSP = [];
% 
% for n = 1:size(all_CCEPmat,1)
%     if all_CCEPmat(n) == 1 && all_ERSPmat(n) == 1
%         spikessampCCEPERSP = [spikessampCCEPERSP; all_spikessamp(n,:)]; 
%         metaFileCCEPERSP = [metaFileCCEPERSP; all_metaFile(n,:)]; 
% 
%     elseif all_CCEPmat(n) ==1 && all_ERSPmat(n) == 0
%         spikessampCCEPnERSP = [spikessampCCEPnERSP; all_spikessamp(n,:)]; 
%         metaFileCCEPnERSP = [metaFileCCEPnERSP; all_metaFile(n,:)]; 
% 
%     elseif all_CCEPmat(n) == 0 && all_ERSPmat(n) == 1
%         spikessampnCCEPERSP = [spikessampnCCEPERSP; all_spikessamp(n,:)]; 
%         metaFilenCCEPERSP = [metaFilenCCEPERSP; all_metaFile(n,:)]; 
% 
%     elseif all_CCEPmat(n) == 0 && all_ERSPmat(n) == 0
%         spikessampnCCEPnERSP = [spikessampnCCEPnERSP; all_spikessamp(n,:)]; 
%         metaFilenCCEPnERSP = [metaFilenCCEPnERSP; all_metaFile(n,:)]; 
% 
%     end
% end
% 
% % housekeeping
% clear I n
% 
% % compare number of spikes per time window for CCEPERSP, CCEPnERSP, nCCEPERSP and nCCEPnERSP channels
% tstep = 0.2;
% 
% % pre-allocation: count of spikes
% countCCEPERSP = NaN(size(spikessampCCEPERSP,1),cfg.epoch_length/0.2);
% countnCCEPERSP = NaN(size(spikessampnCCEPERSP,1),cfg.epoch_length/0.2);
% countCCEPnERSP = NaN(size(spikessampCCEPnERSP,1),cfg.epoch_length/0.2);
% countnCCEPnERSP = NaN(size(spikessampnCCEPnERSP,1),cfg.epoch_length/0.2);
% 
% % count number of spikes per time window (tstep)
% for num = 1:size(spikessampCCEPERSP,1)
%     [n,edges ] = histcounts(vertcat(spikessampCCEPERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countCCEPERSP(num,1:size(n,2)) = n/10;
% end
% 
% for num = 1:size(spikessampnCCEPERSP,1)
%     [n,edges ] = histcounts(vertcat(spikessampnCCEPERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countnCCEPERSP(num,1:size(n,2)) = n/10;
% end
% 
% for num = 1:size(spikessampCCEPnERSP,1)
%     [n,edges ] = histcounts(vertcat(spikessampCCEPnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countCCEPnERSP(num,1:size(n,2)) = n/10;
% end
% 
% for num = 1:size(spikessampnCCEPnERSP,1)
%     [n,edges ] = histcounts(vertcat(spikessampnCCEPnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);
% 
%     countnCCEPnERSP(num,1:size(n,2)) = n/10;
% end
% 
% % convert edges to mean between two edges
% edges_mid =  NaN(1,size(edges,2)-1);
% for nn = 1:size(edges,2)-1
%     edges_mid(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); 
% end
% 
% % use edges and plot counts as block signal
% edges_disc = NaN(2,size(edges,2)-1); 
% countnCCEPERSP_disc = NaN([size(countnCCEPERSP,1),2,size(edges,2)-1]); 
% countCCEPERSP_disc = NaN([size(countCCEPERSP,1),2,size(edges,2)-1]);  
% countnCCEPnERSP_disc = NaN([size(countnCCEPnERSP,1),2,size(edges,2)-1]); 
% countCCEPnERSP_disc = NaN([size(countCCEPnERSP,1),2,size(edges,2)-1]); 
% 
% for nn = 1:size(edges,2)-1
%     edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1]; 
%     countnCCEPERSP_disc(:,:,nn) = repmat(countnCCEPERSP(:,nn),1,2);
%     countnCCEPnERSP_disc(:,:,nn) = repmat(countnCCEPnERSP(:,nn),1,2);
%     countCCEPERSP_disc(:,:,nn) = repmat(countCCEPERSP(:,nn),1,2); 
%     countCCEPnERSP_disc(:,:,nn) = repmat(countCCEPnERSP(:,nn),1,2);
% end
% 
% edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1));
% countnCCEPERSP_disc = reshape(countnCCEPERSP_disc, ...
%     size(countnCCEPERSP_disc,1),2*size(countnCCEPERSP,2));
% countCCEPERSP_disc = reshape(countCCEPERSP_disc, ...
%     size(countCCEPERSP_disc,1),2*size(countCCEPERSP,2));
% countnCCEPnERSP_disc = reshape(countnCCEPnERSP_disc, ...
%     size(countnCCEPnERSP_disc,1),2*size(countnCCEPnERSP,2));
% countCCEPnERSP_disc = reshape(countCCEPnERSP_disc, ...
%     size(countCCEPnERSP_disc,1),2*size(countCCEPnERSP,2));
% 
% % number of IED when CCEP & ERSP
% avNumIEDCCEPERSP_temp = mean(countCCEPERSP);
% avNumIEDCCEPERSP_pre = mean(avNumIEDCCEPERSP_temp(round(t(edges_mid),1)<-tstep));
% avNumIEDCCEPERSP = mean(countCCEPERSP_disc);
% sterrNumIEDCCEPERSP = std(countCCEPERSP_disc)./sqrt(size(countCCEPERSP_disc,1));
% lower_errCCEPERSP = avNumIEDCCEPERSP - sterrNumIEDCCEPERSP;
% upper_errCCEPERSP = avNumIEDCCEPERSP + sterrNumIEDCCEPERSP;
% 
% % number of IED when no CCEP but ERSP
% avNumIEDnCCEPERSP_temp = mean(countnCCEPERSP);
% avNumIEDnCCEPERSP_pre = mean(avNumIEDnCCEPERSP_temp(round(t(edges_mid),1)<-tstep));
% avNumIEDnCCEPERSP = mean(countnCCEPERSP_disc);
% sterrNumIEDnCCEPERSP = std(countnCCEPERSP_disc)./sqrt(size(countnCCEPERSP_disc,1));
% lower_errnCCEPERSP = avNumIEDnCCEPERSP - sterrNumIEDnCCEPERSP;
% upper_errnCCEPERSP = avNumIEDnCCEPERSP + sterrNumIEDnCCEPERSP;
% 
% % number of IED when CCEP but no ERSP
% avNumIEDCCEPnERSP_temp = mean(countCCEPnERSP);
% avNumIEDCCEPnERSP_pre = mean(avNumIEDCCEPnERSP_temp(round(t(edges_mid),1)<-tstep));
% avNumIEDCCEPnERSP = mean(countCCEPnERSP_disc);
% sterrNumIEDCCEPnERSP = std(countCCEPnERSP_disc)./sqrt(size(countCCEPnERSP_disc,1));
% lower_errCCEPnERSP = avNumIEDCCEPnERSP - sterrNumIEDCCEPnERSP;
% upper_errCCEPnERSP = avNumIEDCCEPnERSP + sterrNumIEDCCEPnERSP;
% 
% % number of IED when no CCEP & no ERSP
% avNumIEDnCCEPnERSP_temp = mean(countnCCEPnERSP);
% avNumIEDnCCEPnERSP_pre = mean(avNumIEDnCCEPnERSP_temp(round(t(edges_mid),1)<-tstep));
% avNumIEDnCCEPnERSP = mean(countnCCEPnERSP_disc);
% sterrNumIEDnCCEPnERSP = std(countnCCEPnERSP_disc)./sqrt(size(countnCCEPnERSP_disc,1));
% lower_errnCCEPnERSP = avNumIEDnCCEPnERSP - sterrNumIEDnCCEPnERSP;
% upper_errnCCEPnERSP = avNumIEDnCCEPnERSP + sterrNumIEDnCCEPnERSP;
% 
% % statistical analysis: mann whitney u test
% n_poststim = find(round(t(edges_mid),1)>tstep,1,'first'); 
% p = NaN(4,size(countCCEPERSP,2)-n_poststim+1);
% for n = n_poststim:size(countCCEPERSP,2)
%     CCEPERSP_pre = mean(countCCEPERSP(:,round(t(edges_mid),1)<-tstep),2);
%     nCCEPERSP_pre = mean(countnCCEPERSP(:,round(t(edges_mid),1)<-tstep),2);
%     CCEPnERSP_pre = mean(countCCEPnERSP(:,round(t(edges_mid),1)<-tstep),2);
%     nCCEPnERSP_pre = mean(countnCCEPnERSP(:,round(t(edges_mid),1)<-tstep),2);
%     [~, p(1,n-n_poststim+1)] = ttest(countCCEPERSP(:,n),CCEPERSP_pre(:)); 
%     [~, p(2,n-n_poststim+1)] = ttest(countnCCEPERSP(:,n),nCCEPERSP_pre(:));
%     [~, p(3,n-n_poststim+1)] = ttest(countCCEPnERSP(:,n),CCEPnERSP_pre(:));
%     [~, p(4,n-n_poststim+1)] = ttest(countnCCEPnERSP(:,n),nCCEPnERSP_pre(:));
% end
% 
% m = length(p(:));
% [p_sort,p_ind] = sort(p(:));
% thisVal = NaN(size(p_sort));
% for kk = 1:length(p_sort)
%     thisVal(kk) = (kk/m)*pFDR;
% end
% 
% p_sig = p;
% p_sig(p_ind) = p_sort<thisVal;
% 
% % plot figures
% ymin = round(0.9*min([avNumIEDnCCEPERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDCCEPERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDnCCEPnERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDnCCEPnERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep)]),2,'significant');
% ymax = round(1.1*max([avNumIEDnCCEPERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDCCEPERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDnCCEPnERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep), ...
%     avNumIEDCCEPnERSP(t(edges_disc)<-tstep| t(edges_disc)>tstep)]),2,'significant');
% 
% figure,
% subplot(4,1,1);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))], ...
%     [upper_errCCEPERSP(t(edges_disc)<-tstep), flip(lower_errCCEPERSP(t(edges_disc)<-tstep))], ...
%     cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errCCEPERSP(t(edges_disc)>tstep) flip(lower_errCCEPERSP(t(edges_disc)>tstep))], ...
%     cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% plot(t(edges_disc),avNumIEDCCEPERSP_pre*ones(size(edges_disc)),'-.','Color',[200/256, 200/256, 200/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDCCEPERSP(t(edges_disc)<-tstep),'Color',cmap(10,:),'LineWidth',2);
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDCCEPERSP(t(edges_disc)>tstep),'Color',cmap(10,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin, ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 1; %CCEPERSP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when CCEP and ERSP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% subplot(4,1,2);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))], ...
%     [upper_errnCCEPERSP(t(edges_disc)<-tstep) flip(lower_errnCCEPERSP(t(edges_disc)<-tstep))], ...
%     cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errnCCEPERSP(t(edges_disc)>tstep) flip(lower_errnCCEPERSP(t(edges_disc)>tstep))], ...
%     cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% plot(t(edges_disc),avNumIEDnCCEPERSP_pre*ones(size(edges_disc)),'-.','Color',[180/256, 180/256, 180/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDnCCEPERSP(t(edges_disc)<-tstep),'Color',cmap(20,:),'LineWidth',2);
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDnCCEPERSP(t(edges_disc)>tstep),'Color',cmap(20,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin, ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 2; %nCCEPERSP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when no CCEP but ERSP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% subplot(4,1,3);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))], ...
%     [upper_errCCEPnERSP(t(edges_disc)<-tstep) flip(lower_errCCEPnERSP(t(edges_disc)<-tstep))], ...
%     cmap(30,:),'EdgeColor',cmap(30,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errCCEPnERSP(t(edges_disc)>tstep) flip(lower_errCCEPnERSP(t(edges_disc)>tstep))], ...
%     cmap(30,:),'EdgeColor',cmap(30,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% plot(t(edges_disc),avNumIEDCCEPnERSP_pre*ones(size(edges_disc)),'-.','Color',[180/256, 180/256, 180/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDCCEPnERSP(t(edges_disc)<-tstep),'Color',cmap(30,:),'LineWidth',2);
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDCCEPnERSP(t(edges_disc)>tstep),'Color',cmap(30,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin,ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 3; %CCEPnERSP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when CCEP but no ERSP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% subplot(4,1,4);
% fill([t(edges_disc(t(edges_disc)<-tstep)), flip(t(edges_disc(t(edges_disc)<-tstep)))], ...
%     [upper_errnCCEPnERSP(t(edges_disc)<-tstep) flip(lower_errnCCEPnERSP(t(edges_disc)<-tstep))], ...
%     cmap(40,:),'EdgeColor',cmap(40,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% hold on
% fill([t(edges_disc(t(edges_disc)>tstep)), flip(t(edges_disc(t(edges_disc)>tstep)))], ...
%     [upper_errnCCEPnERSP(t(edges_disc)>tstep) flip(lower_errnCCEPnERSP(t(edges_disc)>tstep))], ...
%     cmap(40,:),'EdgeColor',cmap(40,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
% plot(t(edges_disc),avNumIEDnCCEPnERSP_pre*ones(size(edges_disc)),'-.','Color',[180/256, 180/256, 180/256])
% plot(t(edges_disc(t(edges_disc)<-tstep)),avNumIEDnCCEPnERSP(t(edges_disc)<-tstep),'Color',cmap(40,:),'LineWidth',2);
% plot(t(edges_disc(t(edges_disc)>tstep)),avNumIEDnCCEPnERSP(t(edges_disc)>tstep),'Color',cmap(40,:),'LineWidth',2)
% 
% ylim([0.1 0.3]) %ylim([ymin,ymax])
% 
% % add p<0.05 with FDR correction
% y = 0.95*ymax;
% for n=1:size(p_sig,2) % number of time windows (9)
%     m = 4; %nCCEPnERSP
%     if p_sig(m,n) == 1 % significant
%         if p(m,n)<0.001
%             text(t(edges_mid(n_poststim+n-1)),y,'***','Color',cmap(m*10,:));
%         elseif p(m,n) <0.01
%             text(t(edges_mid(n_poststim+n-1)),y,'**','Color',cmap(m*10,:));
%         elseif p(m,n) < 0.05
%             text(t(edges_mid(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
%         end
%     end
% end
% 
% title('Number of Spikes when no CCEP and no ERSP')
% xlabel('Time (s)')
% ylabel('Mean number of spikes')
% 
% h = gcf;
% h.Units = 'normalized';
% h.Position = [0.5 -0.2 0.4 1.1];
% 
% % save figure
% figureName = sprintf('%s/fig6_IEDsCCEPERSP',...
%     myDataPath.Figures);
% 
% set(gcf,'PaperPositionMode','auto','Renderer','painters')
% print('-dpng','-r300',figureName)
% print('-vector','-depsc',figureName)
% 
% fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)
% 
% % housekeeping
% % clear avNumIEDCCEPERSP avNumIEDCCEPERSP_pre avNumIEDCCEPERSP_temp CCEPERSP_pre nCCEPERSP_pre nn
% % clear avNumIEDCCEPnERSP avNumIEDCCEPnERSP_pre avNumIEDCCEPnERSP_temp CCEPnERSP_pre nCCEPnERSP_pre 
% % clear avNumIEDnCCEPERSP avNumIEDnCCEPERSP_pre avNumIEDnCCEPERSP_temp num p p_ind p_sort p_sig s
% % clear avNumIEDnCCEPnERSP avNumIEDnCCEPnERSP_pre avNumIEDnCCEPnERSP_temp 
% % clear countCCEPERSP countCCEPERSP_disc countCCEPERSP_sel n n_poststim
% % clear countCCEPnERSP countCCEPnERSP_disc countCCEPnERSP_sel 
% % clear countnCCEPERSP countnCCEPERSP_disc countnCCEPERSP_sel m metaFileCCEPERSP metaFilenCCEPERSP
% % clear countnCCEPnERSP countnCCEPnERSP_disc countnCCEPnERSP_sel metaFileCCEPnERSP metaFilenCCEPnERSP
% % clear edges edges_mid edges_disc edges_disc h kk lower_errCCEPERSP lower_errnCCEPERSP lower_errCCEPnERSP lower_errnCCEPnERSP
% % clear spikessampCCEPERSP spikessampnCCEPERSP sterrNumIEDCCEPERSP sterrNumIEDnCCEPERSP thisVal upper_errnCCEPERSP
% % clear spikessampCCEPnERSP spikessampnCCEPnERSP sterrNumIEDCCEPnERSP sterrNumIEDnCCEPnERSP upper_errnCCEPnERSP
% % clear upper_errCCEPERSP y ymax ymin 
% % clear upper_errCCEPnERSP figureName
% 
% 
%% end