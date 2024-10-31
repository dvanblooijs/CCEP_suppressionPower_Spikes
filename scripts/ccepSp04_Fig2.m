%% Figure 2: odds ratio CCEPs vs Power Suppression
% in all separate subjects
% in all subjects combined

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

%% general parameters
pFDR = 0.05;
cmap = parula(101); % 0-1 with steps of 0.01

%% combine all matrices with CCEPs and ERSPs of all subjects
all_CCEPmat = [];
all_ERSPmat = [];

for nSubj = 1:size(dataBase,2)

    all_CCEPmat = [all_CCEPmat; dataBase(nSubj).CCEPmat(:)]; 
    all_ERSPmat = [all_ERSPmat; dataBase(nSubj).ERSPmat(:)]; 
end

% housekeeping
clear nSubj

%% calculate odds ratio in all separate subjects

% determine co-occurrence of ccep and ERSP suppression
clc

for nSubj = 1:size(dataBase,2)
    tbl = crosstab(dataBase(nSubj).CCEPmat(:),dataBase(nSubj).ERSPmat(:));

    [~,p,stats] = fishertest(tbl);

    dataBase(nSubj).stat_ccep_ERSP.OR = stats.OddsRatio; 
    dataBase(nSubj).stat_ccep_ERSP.CI = stats.ConfidenceInterval; 
    dataBase(nSubj).stat_ccep_ERSP.tbl = tbl;
    dataBase(nSubj).stat_ccep_ERSP.p = p;  

    fprintf('--- %s: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
        dataBase(nSubj).sub_label, stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)

end

% housekeeping
clear stats h p nSubj tbl

%% calculate odds ratio in all patients combined

% determine co-occurrence of ccep and ERSP suppression
clc

tbl = crosstab(all_CCEPmat,all_ERSPmat);

[~,p,stats] = fishertest(tbl);

stat_ccep_ERSP.OR = stats.OddsRatio;
stat_ccep_ERSP.CI = stats.ConfidenceInterval;
stat_ccep_ERSP.tbl = tbl;
stat_ccep_ERSP.p = p;

fprintf('--- all subjects: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
    stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)

% housekeeping
clear stats h p nSubj tbl

%% apply FDR correction: p<0.001

pVals = NaN(size(dataBase));
for nSubj = 1:size(dataBase,2)
    pVals(nSubj) = dataBase(nSubj).stat_ccep_ERSP.p;
end
pVals(nSubj+1) = stat_ccep_ERSP.p;

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
for nSubj = 1:size(dataBase,2)
    xline = dataBase(nSubj).stat_ccep_ERSP.CI(1):0.1:dataBase(nSubj).stat_ccep_ERSP.CI(2);

    maxCI(nSubj) = dataBase(nSubj).stat_ccep_ERSP.CI(2);
    plot(xline,-nSubj*ones(size(xline,2),1),'k')

    plot(dataBase(nSubj).stat_ccep_ERSP.OR,-nSubj,'o','MarkerFaceColor',cmap(50,:),'MarkerEdgeColor',cmap(50,:))

% FDR corrected p
    if pSig(nSubj) == 1
        p_str = '***';
    end

    text(dataBase(nSubj).stat_ccep_ERSP.CI(2)+0.5,-nSubj,p_str)%,'FontSize',8)
end

% all patients combined
xline_all = stat_ccep_ERSP.CI(1):0.1:stat_ccep_ERSP.CI(2);

maxCI_all = stat_ccep_ERSP.CI(2);
plot(xline_all,-11*ones(size(xline_all,2),1),'k')

plot(stat_ccep_ERSP.OR,-11,'o','MarkerFaceColor',cmap(30,:),'MarkerEdgeColor',cmap(30,:))

% FDR corrected [
if pSig(nSubj+1)==1
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
    myDataPath.figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear all_CCEPmat all_ERSPmat ax figureName maxCI maxCI_all p_str 
clear h nSubj stat_ccep_ERSP xline xline_all

%% end of script