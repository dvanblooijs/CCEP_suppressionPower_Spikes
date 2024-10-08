% Figure 6: ERSPs and CCEP vs spike ratio

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

%% combine spike ratio and CCEPs/ERSPs of all subjects

all_spikeratio = [];
all_ERSPmat = [];
all_CCEPmat = [];

for subj = 1:size(dataBase,2)
    if any(dataBase(subj).IEDmat)
        all_spikeratio = [all_spikeratio; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>
        ERSPmat = dataBase(subj).ERSPmat(:,dataBase(subj).IEDs.IEDch);
        all_ERSPmat = [all_ERSPmat; ERSPmat(:)]; %#ok<AGROW>
        CCEPmat = dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        all_CCEPmat = [all_CCEPmat; CCEPmat(:)]; %#ok<AGROW>

    end
end

names = cell(size(all_ERSPmat));
[names{all_ERSPmat == 1 & all_CCEPmat == 1}] = deal('cc');
[names{all_ERSPmat == 1 & all_CCEPmat == 0}] = deal('cn');
[names{all_ERSPmat == 0 & all_CCEPmat == 1}] = deal('nc');
[names{all_ERSPmat == 0 & all_CCEPmat == 0}] = deal('nn');
[names{isnan(all_ERSPmat) | isnan(all_CCEPmat)}] = deal('z');

% housekeeping
% clear subj CCEPmat

%% statistics

%% logarithmic spikes
all_spikeratio_log = log(all_spikeratio);

% exclude the infinite values, because otherwise violinplot cannot be made
% value Inf is when the ratio=0, which means that there is no difference
% between pre and post-stimulation
all_spikeratio_log(isinf(all_spikeratio_log)) = NaN;

fprintf('Median (normal SR) powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_log(all_ERSPmat == 1 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (normal SR) no powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_log(all_ERSPmat == 0 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (normal SR) powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_log(all_ERSPmat == 1 & all_CCEPmat == 0),'omitnan'))
fprintf('Median (normal SR) no powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_log(all_ERSPmat == 0 & all_CCEPmat == 0),'omitnan'))

fprintf('Mean (normal SR) powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_log(all_ERSPmat ==1 & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (normal SR) no powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_log(all_ERSPmat ==0 & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (normal SR) powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_log(all_ERSPmat ==1 & all_CCEPmat == 0),'omitnan'))
fprintf('Mean (normal SR) no powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_log(all_ERSPmat ==0 & all_CCEPmat == 0),'omitnan'))

fprintf('--- No significance tested ---\n\n')

%% absolute values of spikes

all_spikeratio_abs = abs(log(all_spikeratio));

[pAbs, tab, stats] = kruskalwallis(all_spikeratio_abs,names,"off");
[a,b,c,d] = multcompare(stats)

% [pAbs,~] = ranksum(all_spikeratio_abs(all_ERSPmat ==1), all_spikeratio_abs(all_ERSPmat == 0));

% exclude the infinite values, because otherwise violinplot cannot be made
% value Inf is when the ratio=0, which means that there is no difference
% between pre and post-stimulation
all_spikeratio_abs(isinf(all_spikeratio_abs)) = NaN;

fprintf('Median (absolute SR) powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_abs(all_ERSPmat ==1 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (absolute SR) no powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_abs(all_ERSPmat ==0  & all_CCEPmat == 1),'omitnan'))
fprintf('Median (absolute SR) powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_abs(all_ERSPmat ==1  & all_CCEPmat == 0),'omitnan'))
fprintf('Median (absolute SR) no powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_abs(all_ERSPmat ==0  & all_CCEPmat == 0),'omitnan'))

fprintf('Mean (absolute SR) powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_abs(all_ERSPmat ==1  & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (absolute SR) no powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_abs(all_ERSPmat ==0  & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (absolute SR) powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_abs(all_ERSPmat ==1  & all_CCEPmat == 0),'omitnan'))
fprintf('Mean (absolute SR) no powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_abs(all_ERSPmat ==0  & all_CCEPmat == 0),'omitnan'))

fprintf('p = %f \n \n',pAbs)

%% negative spikes

all_spikeratio_neg = log(all_spikeratio);
all_spikeratio_neg(all_spikeratio_neg > 0) = NaN;

[pNeg,anovatab,stats] = kruskalwallis(all_spikeratio_neg,names,"off");

% [pNeg,~] = ranksum(all_spikeratio_neg(all_ERSPmat ==1), all_spikeratio_neg(all_ERSPmat == 0));

% exclude the infinite values, because otherwise violinplot cannot be made
% value Inf is when the ratio=0, which means that there is no difference
% between pre and post-stimulation
all_spikeratio_neg(isinf(all_spikeratio_neg)) = NaN;

fprintf('Median (negative SR) powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_neg(all_ERSPmat ==1 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (negative SR) no powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_neg(all_ERSPmat ==0 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (negative SR) powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_neg(all_ERSPmat ==1 & all_CCEPmat == 0),'omitnan'))
fprintf('Median (negative SR) no powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_neg(all_ERSPmat ==0 & all_CCEPmat == 0),'omitnan'))

fprintf('Mean (negative SR) powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_neg(all_ERSPmat ==1 & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (negative SR) no powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_neg(all_ERSPmat ==0 & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (negative SR) powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_neg(all_ERSPmat ==1 & all_CCEPmat == 0),'omitnan'))
fprintf('Mean (negative SR) no powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_neg(all_ERSPmat ==0 & all_CCEPmat == 0),'omitnan'))

fprintf('p = %f \n \n',pNeg)

%% positive spikes

all_spikeratio_pos = log(all_spikeratio);
all_spikeratio_pos(all_spikeratio_pos<0) = NaN;

% [pPos,~] = ranksum(all_spikeratio_pos(all_ERSPmat ==1), all_spikeratio_pos(all_ERSPmat == 0));
[pPos, anovatab,stats] = kruskalwallis(all_spikeratio_pos,names,"off");

% exclude the infinite values, because otherwise violinplot cannot be made
% value Inf is when the ratio=0, which means that there is no difference
% between pre and post-stimulation
all_spikeratio_pos(isinf(all_spikeratio)) = NaN;

fprintf('Median (positive SR) powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_pos(all_ERSPmat ==1 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (positive SR) no powerSup + CCEP = %1.2f\n',...
    median(all_spikeratio_pos(all_ERSPmat ==0 & all_CCEPmat == 1),'omitnan'))
fprintf('Median (positive SR) powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_pos(all_ERSPmat ==1 & all_CCEPmat == 0),'omitnan'))
fprintf('Median (positive SR) no powerSup + no CCEP = %1.2f\n',...
    median(all_spikeratio_pos(all_ERSPmat ==0 & all_CCEPmat == 0),'omitnan'))

fprintf('Mean (positive SR) powerSup + CCEP = %1.2f\n',...
    mean(all_spikeratio_pos(all_ERSPmat ==1 & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (positive SR) no powerSup + CCEP= %1.2f\n',...
    mean(all_spikeratio_pos(all_ERSPmat ==0 & all_CCEPmat == 1),'omitnan'))
fprintf('Mean (positive SR) powerSup + no CCEP = %1.2f\n',...
    mean(all_spikeratio_pos(all_ERSPmat ==1 & all_CCEPmat == 0),'omitnan'))
fprintf('Mean (positive SR) no powerSup + no CCEP= %1.2f\n',...
    mean(all_spikeratio_pos(all_ERSPmat ==0 & all_CCEPmat == 0),'omitnan'))

fprintf('p = %f \n \n',pPos)

%% apply FDR correction: p<0.05

pVals = [pAbs, pNeg, pPos];

[pSort,pInd] = sort(pVals(:));

m = length(pVals);
thisVal = NaN(size(pSort));
for kk = 1:length(pSort)
    thisVal(kk) = (kk/m)*0.05;
end

pSig = pVals;
pSig(pInd) = pSort < thisVal;

%% figures: violin plots

%% violin plot with logarithmic spike ratio
% in all subjects combined

ymin = min(all_spikeratio_log);
ymax = max(all_spikeratio_log);

h = figure(5);
violinplot(all_spikeratio_log,names,'Width',0.3);

h.Units = 'normalized';
h.Position = [0.1 0.1 0.21 0.8];

ylim([ymin ymax])
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 60;
ax.XTickLabel = {'PS + CCEP','PS + no CCEP','no PS + CCEP','no PS + no CCEP',''};
ax.YLabel.String = 'Logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig6_SR_powerSup_CCEP_normal',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear figureName h ymax ymin ans ax

%% violin plot with absolute logarithmic spike ratio
% in all subjects combined
% statistical test with mann whitney u test

ymax = max(all_spikeratio_abs);

h = figure(6);
violinplot(all_spikeratio_abs,names,'Width',0.3);
hold on

% if pAbs < 0.1 && pAbs > 0.05
%     text(1.5,ymax-0.4,'~')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% elseif pAbs < 0.05 && pAbs > 0.01
%     text(1.5,ymax-0.4,'*')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% elseif pAbs < 0.01 && pAbs > 0.001
%     text(1.5,ymax-0.4,'**')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% elseif pAbs < 0.001
%     text(1.45,ymax-0.4,'***')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% end

% FDR corrected p
if pSig(1) == 1
    text(1.45,ymax-0.4,'***')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
end

hold off

h.Units = 'normalized';
h.Position = [0.3 0.1 0.21 0.8];

ylim([-0.1 ymax])
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 60;
ax.XTickLabel = {'PS + CCEP','PS + no CCEP','no PS + CCEP','no PS + no CCEP',''};
ax.YLabel.String = 'Logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig6_SR_powerSup_CCEP_absolute',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear ans ax figureName h ymax

%% violin plot with only negative logarithmic spike ratios
% which means a decrease in spikes after stimulation
% statistical test with mann whitney u test

ymin = min(all_spikeratio_neg);
ymax = 0.1;

h = figure(7);
violinplot(all_spikeratio_neg,names,'Width',0.3);
hold on

% if pNeg < 0.1 && pNeg > 0.05
%     text(1.5,ymin+0.4,'~')
%     plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
% elseif pNeg < 0.05 && pNeg > 0.01
%     text(1.5,ymin+0.4,'*')
%     plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
% elseif pNeg < 0.01 && pNeg > 0.001
%     text(1.5,ymin+0.4,'**')
%     plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
% elseif pNeg < 0.001
%     text(1.45,ymin+0.4,'***')
%     plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
% end

% FDR corrected p
if pSig(2) == 1
    text(1.45,ymin+0.4,'***')
    plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
end

hold off

h.Units = 'normalized';
h.Position = [0.5 0.1 0.21 0.8];

ylim([ymin ymax])
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 60;
ax.XTickLabel = {'PS + CCEP','PS + no CCEP','no PS + CCEP','no PS + no CCEP',''};
ax.YLabel.String = 'Negative logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig6_SR_powerSup_CCEP_negative',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear ans ax figureName h p ymax ymin

%% violin plot with only positive logarithmic spike ratios
% which means an increase in spikes after stimulation
% statistical analysis with mann whitney u test

ymin = -0.1;
ymax = max(all_spikeratio_pos);

h = figure(8);
violinplot(all_spikeratio_pos,names,'Width',0.3);
hold on

% if pPos < 0.1 && pPos > 0.05
%     text(1.5,ymax-0.4,'~')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% elseif pPos < 0.05 && pPos > 0.01
%     text(1.5,ymax-0.4,'*')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% elseif pPos < 0.01 && pPos > 0.001
%     text(1.5,ymax-0.4,'**')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% elseif pPos < 0.001
%     text(1.45,ymax-0.4,'***')
%     plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
% end

% FDR corrected p
if pSig(3) == 1
    text(1.45,ymax-0.4,'***')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
end

hold off

h.Units = 'normalized';
h.Position = [0.7 0.1 0.21 0.8];

ylim([ymin ymax])
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 60;
ax.XTickLabel = {'PS + CCEP','PS + no CCEP','no PS + CCEP','no PS + no CCEP',''};
ax.YLabel.String = 'Logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig6_SR_powerSup_CCEP_positive',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear ans ax figureName h p ymax ymin

%% end

% housekeeping
clear all_ERSPmat all_spikeratio all_spikeratio_abs all_spikeratio_log
clear all_spikeratio_neg all_spikeratio_pos names
