% Figure 3: CCEPs vs spike ratio
close all
clc

%% violin plot with logarithmic spike ratio
% in all patients combined

spikeratioAll = [];
CCEPmatAll = [];
for subj = 1:size(dataBase,2)
    if any(dataBase(subj).IEDmat)
        spikeratioAll = [spikeratioAll; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>
        CCEPmat = dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        CCEPmatAll = [CCEPmatAll; CCEPmat(:)]; %#ok<AGROW>
    end
end

spikeratioAll = log(spikeratioAll);

names = cell(size(CCEPmatAll));
[names{CCEPmatAll == 1}] = deal('c');
[names{CCEPmatAll == 0}] = deal('nc');
[names{isnan(CCEPmatAll)}] = deal('z');

spikeratioAll(isinf(spikeratioAll)) = NaN;

fprintf('Median (normal SR) connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Median (normal SR) not connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==0),'omitnan'))

fprintf('Mean (normal SR) connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Mean (normal SR) not connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==0),'omitnan'))

ymin = min(spikeratioAll);
ymax = max(spikeratioAll);

h = figure(1);
violinplot(spikeratioAll,names,'Width',0.3);

h.Units = 'normalized';
h.Position = [0.1 0.1 0.21 0.8];

ylim([ymin ymax])
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 60;
ax.XTickLabel = {'connected','not connected',' '};
ax.YLabel.String = 'Logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig4_SR_CCEP_normal',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

%% violin plot with absolute logarithmic spike ratio
% in all patients combined
% statistical test with mann whitney u test

spikeratioAll = [];
CCEPmatAll = [];
for subj = 1:size(dataBase,2)
    if any(dataBase(subj).IEDmat)
        spikeratioAll = [spikeratioAll; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>
        CCEPmat = dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        CCEPmatAll = [CCEPmatAll; CCEPmat(:)]; %#ok<AGROW>
    end
end

spikeratioAll = abs(log(spikeratioAll));

[p,~] = ranksum(spikeratioAll(CCEPmatAll ==1), spikeratioAll(CCEPmatAll == 0));

names = cell(size(CCEPmatAll));
[names{CCEPmatAll == 1}] = deal('c');
[names{CCEPmatAll == 0}] = deal('nc');
[names{isnan(CCEPmatAll)}] = deal('z');

spikeratioAll(isinf(spikeratioAll)) = NaN;

fprintf('Median (absolute SR) connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Median (absolute SR) not connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==0),'omitnan'))

fprintf('Mean (absolute SR) connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Mean (absolute SR) not connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==0),'omitnan'))

ymax = max(spikeratioAll);

h = figure(2);
violinplot(spikeratioAll,names,'Width',0.3);
hold on

if p < 0.1 && p > 0.05
    text(1.5,ymax-0.4,'~')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
elseif p < 0.05 && p > 0.01
    text(1.5,ymax-0.4,'*')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
elseif p < 0.01 && p > 0.001
    text(1.5,ymax-0.4,'**')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
elseif p < 0.001
    text(1.45,ymax-0.4,'***')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
end

hold off

h.Units = 'normalized';
h.Position = [0.3 0.1 0.21 0.8];

ylim([0 ymax])
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 60;
ax.XTickLabel = {'connected','not connected',' '};
ax.YLabel.String = 'Absolute logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig4_SR_CCEP_absolute',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

%% violin plot with only negative logarithmic spike ratios
% which means a decrease in spikes after stimulation
% statistical test with mann whitney u test

spikeratioAll = [];
CCEPmatAll = [];
for subj = 1:size(dataBase,2)
    if any(dataBase(subj).IEDmat)
        spikeratioAll = [spikeratioAll; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>
        CCEPmat = dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        CCEPmatAll = [CCEPmatAll; CCEPmat(:)]; %#ok<AGROW>
    end
end

spikeratioAll = log(spikeratioAll);
spikeratioAll(spikeratioAll > 0) = NaN;

[p,~] = ranksum(spikeratioAll(CCEPmatAll ==1), spikeratioAll(CCEPmatAll == 0));

names = cell(size(CCEPmatAll));
[names{CCEPmatAll == 1}] = deal('c');
[names{CCEPmatAll == 0}] = deal('nc');
[names{isnan(CCEPmatAll)}] = deal('z');

spikeratioAll(isinf(spikeratioAll)) = NaN;

fprintf('Median (negative SR) connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Median (negative SR) not connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==0),'omitnan'))

fprintf('Mean (negative SR) connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Mean (negative SR) not connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==0),'omitnan'))

ymin = min(spikeratioAll);
ymax = 0.1;

h = figure(3);
violinplot(spikeratioAll,names,'Width',0.3);
hold on

if p < 0.1 && p > 0.05
    text(1.5,ymin+0.4,'~')
    plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
elseif p < 0.05 && p > 0.01
    text(1.5,ymin+0.4,'*')
    plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
elseif p < 0.01 && p > 0.001
    text(1.5,ymin+0.4,'**')
    plot(1.1:0.1:1.9,ymin+0.37*ones(9,1),'k')
elseif p < 0.001
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
ax.XTickLabel = {'connected','not connected',' '};
ax.YLabel.String = 'Negative logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig4_SR_CCEP_negative',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

%% violin plot with only positive logarithmic spike ratios
% which means an increase in spikes after stimulation
% statistical analysis with mann whitney u test

spikeratioAll = [];
CCEPmatAll = [];
for subj = 1:size(dataBase,2)
    if any(dataBase(subj).IEDmat)
        spikeratioAll = [spikeratioAll; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>
        CCEPmat = dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        CCEPmatAll = [CCEPmatAll; CCEPmat(:)]; %#ok<AGROW>
    end
end

spikeratioAll = log(spikeratioAll);

spikeratioAll(spikeratioAll<0) = NaN;

[p,~] = ranksum(spikeratioAll(CCEPmatAll ==1), spikeratioAll(CCEPmatAll == 0));

names = cell(size(CCEPmatAll));
[names{CCEPmatAll == 1}] = deal('c');
[names{CCEPmatAll == 0}] = deal('nc');
[names{isnan(CCEPmatAll)}] = deal('z');

spikeratioAll(isinf(spikeratioAll)) = NaN;

fprintf('Median (positive SR) connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Median (positive SR) not connected = %1.2f\n',...
    median(spikeratioAll(CCEPmatAll ==0),'omitnan'))

fprintf('Mean (positive SR) connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==1),'omitnan'))
fprintf('Mean (positive SR) not connected = %1.2f\n',...
    mean(spikeratioAll(CCEPmatAll ==0),'omitnan'))

ymin = -0.1;
ymax = max(spikeratioAll);

h = figure(4);
violinplot(spikeratioAll,names,'Width',0.3);
hold on

if p < 0.1 && p > 0.05
    text(1.5,ymax-0.4,'~')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
elseif p < 0.05 && p > 0.01
    text(1.5,ymax-0.4,'*')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
elseif p < 0.01 && p > 0.001
    text(1.5,ymax-0.4,'**')
    plot(1.1:0.1:1.9,ymax-0.43*ones(9,1),'k')
elseif p < 0.001
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
ax.XTickLabel = {'connected','not connected',' '};
ax.YLabel.String = 'Positive logarithmic spike ratio';
ax.XLabel.String = ' ';
ax.Title.String = 'Spike ratio';

figureName = sprintf('%s/fig4_SR_CCEP_positive',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

%% end