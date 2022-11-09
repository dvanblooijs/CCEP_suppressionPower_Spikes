%% Figure 3: odds ratio CCEPs vs Power Suppression
% in all separate subjects
% in all subjects combined

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

%% combine all matrices with CCEPs and ERSPs of all subjects
all_CCEPmat = [];
all_ERSPmat = [];

for subj = 1:size(dataBase,2)

    all_CCEPmat = [all_CCEPmat; dataBase(subj).CCEPmat(:)]; %#ok<AGROW>
    all_ERSPmat = [all_ERSPmat; dataBase(subj).ERSPmat(:)];  %#ok<AGROW>
end

% housekeeping
clear subj

%% calculate odds ratio in all separate subjects

% determine co-occurrence of ccep and ERSP suppression
clc

for subj = 1:size(dataBase,2)
    tbl = crosstab(dataBase(subj).CCEPmat(:),dataBase(subj).ERSPmat(:));

    [~,p,stats] = fishertest(tbl);

    dataBase(subj).stat_ccep_ERSP.OR = stats.OddsRatio; %#ok<SAGROW> 
    dataBase(subj).stat_ccep_ERSP.CI = stats.ConfidenceInterval; %#ok<SAGROW> 
    dataBase(subj).stat_ccep_ERSP.tbl = tbl; %#ok<SAGROW> 
    dataBase(subj).stat_ccep_ERSP.p = p;  %#ok<SAGROW> 

    fprintf('--- %s: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
        dataBase(subj).sub_label, stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)

end

% housekeeping
clear stats h p subj tbl

%% calculate odds ratio in all patients combined

% determine co-occurrence of ccep and ERSP suppression
clc

cmap = parula(101); % 0-1 with steps of 0.01
tbl = crosstab(all_CCEPmat,all_ERSPmat);

[~,p,stats] = fishertest(tbl);

stat_ccep_ERSP.OR = stats.OddsRatio;
stat_ccep_ERSP.CI = stats.ConfidenceInterval;
stat_ccep_ERSP.tbl = tbl;
stat_ccep_ERSP.p = p;

fprintf('--- all subjects: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
    stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)

% housekeeping
clear stats h p subj tbl

%% apply FDR correction: p<0.001
pFDR = 0.05;

pVals = NaN(size(dataBase));
for subj = 1:size(dataBase,2)
    pVals(subj) = dataBase(subj).stat_ccep_ERSP.p;
end
pVals(subj+1) = stat_ccep_ERSP.p;

[pSort,pInd] = sort(pVals(:));

m = length(pVals);
thisVal = NaN(size(pSort));
for kk = 1:length(pSort)
    thisVal(kk) = (kk/m)*pFDR;
end

pSig = pVals;
pSig(pInd) = pSort < thisVal;

%% make forest plot cceps vs ERSP
close all

maxCI = NaN(size(dataBase,2),1);
h = figure;
hold on
for subj = 1:size(dataBase,2)
    xline = dataBase(subj).stat_ccep_ERSP.CI(1):0.1:dataBase(subj).stat_ccep_ERSP.CI(2);

    maxCI(subj) = dataBase(subj).stat_ccep_ERSP.CI(2);
    plot(xline,-subj*ones(size(xline,2),1),'k')

    plot(dataBase(subj).stat_ccep_ERSP.OR,-subj,'o','MarkerFaceColor',cmap(50,:),'MarkerEdgeColor',cmap(50,:))

%     if dataBase(subj).stat_ccep_ERSP.p <0.001
%         p_str = '***';
%     elseif dataBase(subj).stat_ccep_ERSP.p <0.01 && dataBase(subj).stat_ccep_ERSP.p > 0.001
%         p_str = '**';
%     elseif dataBase(subj).stat_ccep_ERSP.p <0.05 && dataBase(subj).stat_ccep_ERSP.p > 0.01
%         p_str = '*';
%     end

% FDR corrected p
    if pSig(subj) == 1
        p_str = '***';
    end

    text(dataBase(subj).stat_ccep_ERSP.CI(2)+0.5,-subj,p_str)%,'FontSize',8)
end

% all patients combined
xline_all = stat_ccep_ERSP.CI(1):0.1:stat_ccep_ERSP.CI(2);

maxCI_all = stat_ccep_ERSP.CI(2);
plot(xline_all,-11*ones(size(xline_all,2),1),'k')

plot(stat_ccep_ERSP.OR,-11,'o','MarkerFaceColor',cmap(30,:),'MarkerEdgeColor',cmap(30,:))

% if stat_ccep_ERSP.p <0.001
%     p_str = '***';
% elseif stat_ccep_ERSP.p <0.01 && stat_ccep_ERSP.p > 0.001
%     p_str = '**';
% elseif stat_ccep_ERSP.p <0.05 && stat_ccep_ERSP.p > 0.01
%     p_str = '*';
% end

if pSig(subj+1)==1
    p_str = '***';
end

text(stat_ccep_ERSP.CI(2)+0.5,-11,p_str)%,'FontSize',8)

plot(ones(size(dataBase,2)+3,1),0:size(dataBase,2)+2,'k:')
hold off

ylim([-size(dataBase,2)-2 0])
xlim([0 ceil(max([maxCI; maxCI_all]))+2])

ax = gca;
ax.YTick = -12:0;
ax.YTickLabel = flip([{' '},replace({dataBase(:).sub_label},'sub-',''),{'all subjects combined'},{' '}]);
ylabel('Subjects')
xlabel('Odds ratio')
title('Occurrence of power suppression and CCEP')

h.Units = 'normalized';
h.Position = [0.6 0.37 0.5 0.5];

figureName = sprintf('%s/fig3_OR',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear all_CCEPmat all_ERSPmat ax figureName maxCI maxCI_all p_str 
clear h subj stat_ccep_ERSP xline xline_all




%% end of script