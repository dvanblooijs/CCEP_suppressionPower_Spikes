% Figure 6: spikes in epochs (not) connected and (no) suppressed power

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

%% combine ERSPs, CCEPs (only include the IED channels), 
% and spikes, spike ratios of all subjects

all_ERSPmat = [];
all_CCEPmat = [];
all_spikessamp = [];
all_spikeratio = [];
all_metaFile = [];

for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)

        ERSPmat = dataBase(subj).ERSPmat(:,dataBase(subj).IEDs.IEDch);
        all_ERSPmat = [all_ERSPmat; ERSPmat(:)]; %#ok<AGROW> % with ERSPmat(:), all columns will be put below each other, so first column 1, then below this column 2 etc.
        
        CCEPmat =  dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        all_CCEPmat = [all_CCEPmat; CCEPmat(:)]; %#ok<AGROW>
        
        % spikessamp had a size [stimp x chan x trials], and we want to
        % reshape it into [stimp*chan x trials]. 
        spikessamp = dataBase(subj).IEDs.totEpochspikessamp;
        [m,n,o] = size(spikessamp);
        spikessamp_reshape = reshape(spikessamp,m*n,o); % for channel 1, all stimps are mentioned, then channel 2 etc. 
        all_spikessamp = [all_spikessamp; spikessamp_reshape]; %#ok<AGROW>
        
        all_spikeratio = [all_spikeratio; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>

        % this metaFile contains 4 columns: 
        % 1: subject number, 2&3: stimulation pair, 4: channel
        metaFile = [subj*ones(size(ERSPmat,1)*size(ERSPmat,2),1),...
            repmat(dataBase(subj).cc_stimsets,size(ERSPmat,2),1),...
            reshape(repmat(dataBase(subj).IEDs.IEDch',size(ERSPmat,1),1),m*n,1)];
        all_metaFile = [all_metaFile; metaFile]; %#ok<AGROW>
    end
end

% housekeeping
clear CCEPmat ERSPmat m n o spikessamp spikessamp_reshape subj metaFile

%% add all spike samples into two categories: nCCEP or CCEP

% first sort in order of ascending spike ratio
[~,I] = sort(log(all_spikeratio),'ascend');

spikessampCCEP = [];
spikessampnCCEP = [];
metaFileCCEP = []; metaFilenCCEP = [];

for n = 1:size(I,1)
    num = I(n);
    if all_CCEPmat(num) == 1
        spikessampCCEP = [spikessampCCEP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFileCCEP = [metaFileCCEP; all_metaFile(num,:)]; %#ok<AGROW>

    elseif all_CCEPmat(num) == 0
        spikessampnCCEP = [spikessampnCCEP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFilenCCEP = [metaFilenCCEP; all_metaFile(num,:)]; %#ok<AGROW>
    end
end

% housekeeping
clear I num n

% compare number of spikes per time window for CCEP and nCCEP channels

tstep = 0.2;

% pre-allocation
countCCEP_cont = NaN(size(spikessampCCEP,1),cfg.epoch_length/0.2);
countnCCEP_cont = NaN(size(spikessampnCCEP,1),cfg.epoch_length/0.2);

% count number of spikes per time window (tstep)
for num = 1:size(spikessampCCEP,1)
    [n,edges ] = histcounts(vertcat(spikessampCCEP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countCCEP_cont(num,1:size(n,2)) = n/10; % divided by ten, because ten trials
end

for num = 1:size(spikessampnCCEP,1)
    [n,edges ] = histcounts(vertcat(spikessampnCCEP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnCCEP_cont(num,1:size(n,2)) = n/10;
end

% convert edges to mean between the two edges
edges_cont = [];
for nn = 1:size(edges,2)-1
    edges_cont(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); %#ok<SAGROW>
end

% use edges and plot counts as block signal
edges_disc = []; countnCCEP_disc = []; countCCEP_disc = [];
for nn = 1:size(edges,2)-1
    edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1]; %#ok<SAGROW>
    countnCCEP_disc(:,:,nn) = repmat(countnCCEP_cont(:,nn),1,2); %#ok<SAGROW>
    countCCEP_disc(:,:,nn) = repmat(countCCEP_cont(:,nn),1,2); %#ok<SAGROW>
end

edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1));
countnCCEP_disc = reshape(countnCCEP_disc,size(countnCCEP_disc,1),2*size(countnCCEP_cont,2));
countCCEP_disc = reshape(countCCEP_disc,size(countCCEP_disc,1),2*size(countCCEP_cont,2));

s = input('Do you want to display a block signal or a continuous signal? [block/continuous]: ','s');

if strcmpi(s,'continuous')
    edges_sel = edges_cont;
    countnCCEP_sel = countnCCEP_cont;
    countCCEP_sel = countCCEP_cont;

    % OR ...
elseif strcmpi(s,'block')
    edges_sel = edges_disc;
    countnCCEP_sel = countnCCEP_disc;
    countCCEP_sel = countCCEP_disc;

else
    error('The answer to the previous question was not "block" or "continuous".')
end

% number of IED when CCEP
avNumIEDCCEP_temp = mean(countCCEP_sel);
avNumIEDCCEP_pre = mean(avNumIEDCCEP_temp(t(edges_sel)<-tstep));
avNumIEDCCEP = avNumIEDCCEP_temp ;%- avNumIEDCCEP_pre;
sterrNumIEDCCEP = std(countCCEP_sel)./sqrt(size(countCCEP_sel,1));
lower_errCCEP = avNumIEDCCEP - sterrNumIEDCCEP;
upper_errCCEP = avNumIEDCCEP + sterrNumIEDCCEP;

% number of IED when no CCEP
avNumIEDnCCEP_temp = mean(countnCCEP_sel);
avNumIEDnCCEP_pre = mean(avNumIEDnCCEP_temp(t(edges_sel)<-tstep));
avNumIEDnCCEP = avNumIEDnCCEP_temp ;%- avNumIEDnCCEP_pre;
sterrNumIEDnCCEP = std(countnCCEP_sel)./sqrt(size(countnCCEP_sel,1));
lower_errnCCEP = avNumIEDnCCEP - sterrNumIEDnCCEP;
upper_errnCCEP = avNumIEDnCCEP + sterrNumIEDnCCEP;

% FDR correction
n_poststim = find(t(edges_cont)>tstep,1,'first');
p = [];
for n = n_poststim:size(countCCEP_cont,2)
    CCEP_pre = countCCEP_cont(:,t(edges_cont)<-tstep);
    nCCEP_pre = countnCCEP_cont(:,t(edges_cont)<-tstep);
    p(1,n-n_poststim+1) = ranksum(countCCEP_cont(:,n),CCEP_pre(:)); %#ok<SAGROW> 
    p(2,n-n_poststim+1) = ranksum(countnCCEP_cont(:,n),nCCEP_pre(:)); %#ok<SAGROW> 
end

m = length(p(:));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

% plot figures
ymin = round(0.9*min([avNumIEDnCCEP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDCCEP(t(edges_sel)<-tstep| t(edges_sel)>tstep)]),2,'significant');
ymax = round(1.1*max([avNumIEDnCCEP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDCCEP(t(edges_sel)<-tstep| t(edges_sel)>tstep)]),2,'significant');

figure,
subplot(2,1,1);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))],...
    [upper_errCCEP(t(edges_sel)<-tstep) flip(lower_errCCEP(t(edges_sel)<-tstep))],...
    cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errCCEP(t(edges_sel)>tstep) flip(lower_errCCEP(t(edges_sel)>tstep))], ...
    cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
plot(t(edges_sel),avNumIEDCCEP_pre*ones(size(edges_sel)),'-.','Color',[200/256, 200/256, 200/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDCCEP(t(edges_sel)<-tstep),'Color',cmap(10,:),'LineWidth',2)
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDCCEP(t(edges_sel)>tstep),'Color',cmap(10,:),'LineWidth',2)

ylim([ymin,ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 1; %CCEP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
            %         elseif p(n) <0.01
            %             text(tt(n),0.2,'**')
            %         elseif p(n) < 0.05
            %             text(tt(n),0.2,'*')
        end
    end
end

title('Number of Spikes when CCEP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

subplot(2,1,2);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))],...
    [upper_errnCCEP(t(edges_sel)<-tstep) flip(lower_errnCCEP(t(edges_sel)<-tstep))],...
    cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errnCCEP(t(edges_sel)>tstep) flip(lower_errnCCEP(t(edges_sel)>tstep))], ...
    cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
plot(t(edges_sel),avNumIEDnCCEP_pre*ones(size(edges_sel)),'-.','Color',[180/256, 180/256, 180/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDnCCEP(t(edges_sel)<-tstep),'Color',cmap(20,:),'LineWidth',2)
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDnCCEP(t(edges_sel)>tstep),'Color',cmap(20,:),'LineWidth',2)

ylim([ymin,ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 2; %nCCEP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
        end
    end
end

title('Number of Spikes when no CCEP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

h = gcf;
h.Units = 'normalized';
h.Position = [0.1 0.4 0.4 0.5];

figureName = sprintf('%s/fig6_IEDsCCEP_%s',...
    myDataPath.Figures,s);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear avNumIEDCCEP avNumIEDCCEP_pre avNumIEDCCEP_temp CCEP_pre nCCEP_pre nn
clear avNumIEDnCCEP avNumIEDnCCEP_pre avNumIEDnCCEP_temp num p p_ind p_sort p_sig s
clear countCCEP_cont countCCEP_disc countCCEP_sel n n_poststim
clear countnCCEP_cont countnCCEP_disc countnCCEP_sel m metaFileCCEP metaFilenCCEP
clear edges edges_cont edges_disc edges_sel h kk lower_errCCEP lower_errnCCEP
clear spikessampCCEP spikessampnCCEP sterrNumIEDCCEP sterrNumIEDnCCEP thisVal upper_errnCCEP
clear upper_errCCEP y ymax ymin h figureName ans

%% add all spike samples into two categories: nERSP or ERSP

% first sort in order of ascending spike ratio
[~,I] = sort(log(all_spikeratio),'ascend');

spikessampERSP = []; % ERSP
spikessampnERSP = []; % no ERSP
metaFileERSP = []; metaFilenERSP = [];

for n = 1:size(all_ERSPmat,1)
    num = I(n);
    if all_ERSPmat(num) == 1
        spikessampERSP = [spikessampERSP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFileERSP = [metaFileERSP; all_metaFile(num,:)]; %#ok<AGROW>

    elseif all_ERSPmat(num) == 0
        spikessampnERSP = [spikessampnERSP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFilenERSP = [metaFilenERSP; all_metaFile(num,:)]; %#ok<AGROW>
    end
end

% housekeeping
clear I num n

% compare number of spikes per time window for ERSP and nERSP channels

tstep = 0.2;

% pre-allocation
countERSP_cont = NaN(size(spikessampERSP,1),cfg.epoch_length/0.2);
countnERSP_cont = NaN(size(spikessampnERSP,1),cfg.epoch_length/0.2);

for num = 1:size(spikessampERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countERSP_cont(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnERSP_cont(num,1:size(n,2)) = n/10;
end

% convert edges to mean between the two edges
edges_cont = [];
for nn = 1:size(edges,2)-1
    edges_cont(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); %#ok<SAGROW>
end

% use edges and plot counts as block signal
edges_disc = []; countnERSP_disc = []; countERSP_disc = [];
for nn = 1:size(edges,2)-1
    edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1]; %#ok<SAGROW>
    countnERSP_disc(:,:,nn) = repmat(countnERSP_cont(:,nn),1,2); %#ok<SAGROW>
    countERSP_disc(:,:,nn) = repmat(countERSP_cont(:,nn),1,2); %#ok<SAGROW>
end

edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1));
countnERSP_disc = reshape(countnERSP_disc,size(countnERSP_disc,1),2*size(countnERSP_cont,2));
countERSP_disc = reshape(countERSP_disc,size(countERSP_disc,1),2*size(countERSP_cont,2));

s = input('Do you want to display a block signal or a continuous signal? [block/continuous]: ','s');

if strcmpi(s,'continuous')
    edges_sel = edges_cont;
    countnERSP_sel = countnERSP_cont;
    countERSP_sel = countERSP_cont;

    % OR ...
elseif strcmpi(s,'block')
    edges_sel = edges_disc;
    countnERSP_sel = countnERSP_disc;
    countERSP_sel = countERSP_disc;

else
    error('The answer to the previous question was not "block" or "continuous".')
end

% number of IED when ERSP
% MEAN
avNumIEDERSP_temp = mean(countERSP_sel);
avNumIEDERSP_pre = mean(avNumIEDERSP_temp(t(edges_sel)<-tstep));
avNumIEDERSP = avNumIEDERSP_temp ;%- avNumIEDERSP_pre;
sterrNumIEDERSP = std(countERSP_sel)./sqrt(size(countERSP_sel,1));
lower_errERSP = avNumIEDERSP - sterrNumIEDERSP;
upper_errERSP = avNumIEDERSP + sterrNumIEDERSP;
% MEDIAN 
% avNumIEDERSP_temp = median(countERSP_sel);
% avNumIEDERSP_pre = median(avNumIEDERSP_temp(t(edges_sel)<-tstep));
% avNumIEDERSP = avNumIEDERSP_temp ;%- avNumIEDERSP_pre;
% sterrNumIEDERSP = quantile(countERSP_sel,[.25 .75]);
% lower_errERSP = avNumIEDERSP - sterrNumIEDERSP(1,:);
% upper_errERSP = avNumIEDERSP + sterrNumIEDERSP(2,:);

% number of IED when no ERSP
% MEAN
avNumIEDnERSP_temp = mean(countnERSP_sel);
avNumIEDnERSP_pre = mean(avNumIEDnERSP_temp(t(edges_sel)<-tstep));
avNumIEDnERSP = avNumIEDnERSP_temp ;%- avNumIEDERSP_pre;
sterrNumIEDnERSP = std(countnERSP_sel)./sqrt(size(countnERSP_sel,1));
lower_errnERSP = avNumIEDnERSP - sterrNumIEDnERSP;
upper_errnERSP = avNumIEDnERSP + sterrNumIEDnERSP;
% MEDIAN
% avNumIEDnERSP_temp = median(countnERSP_sel);
% avNumIEDnERSP_pre = median(avNumIEDnERSP_temp(t(edges_sel)<-tstep));
% avNumIEDnERSP = avNumIEDnERSP_temp ;%- avNumIEDERSP_pre;
% sterrNumIEDnERSP = quantile(countnERSP_sel,[.25 .75]);
% lower_errnERSP = avNumIEDnERSP - sterrNumIEDnERSP(1,:);
% upper_errnERSP = avNumIEDnERSP + sterrNumIEDnERSP(2,:);

% FDR correction
n_poststim = find(t(edges_cont)>tstep,1,'first');
p = [];
for n = n_poststim:size(countERSP_cont,2)
    ERSP_pre = countERSP_cont(:,t(edges_cont)<-tstep);
    nERSP_pre = countnERSP_cont(:,t(edges_cont)<-tstep);
    p(1,n-n_poststim+1) = ranksum(countERSP_cont(:,n),ERSP_pre(:)); %#ok<SAGROW> 
    p(2,n-n_poststim+1) = ranksum(countnERSP_cont(:,n),nERSP_pre(:)); %#ok<SAGROW> 
%    [~, p(1,n-n_poststim+1)] = ttest2(countERSP_cont(:,n),ERSP_pre(:)); %#ok<SAGROW> 
%    [~,p(2,n-n_poststim+1)] = ttest2(countnERSP_cont(:,n),nERSP_pre(:)); %#ok<SAGROW>
end

m = length(p(:));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

% plot figures
ymin = round(0.9*min([avNumIEDnERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep)]),2,'significant');
ymax = round(1.1*max([avNumIEDnERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep)]),2,'significant');

figure,
subplot(2,1,1);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))],...
    [upper_errERSP(t(edges_sel)<-tstep) flip(lower_errERSP(t(edges_sel)<-tstep))],...
    cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errERSP(t(edges_sel)>tstep) flip(lower_errERSP(t(edges_sel)>tstep))], ...
    cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5)
plot(t(edges_sel),avNumIEDERSP_pre*ones(size(edges_sel)),'-.','Color',[200/256, 200/256, 200/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDERSP(t(edges_sel)<-tstep),'Color',cmap(10,:),'LineWidth',2)
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDERSP(t(edges_sel)>tstep),'Color',cmap(10,:),'LineWidth',2)

ylim([ymin,ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 1; %ERSP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
            %         elseif p(n) <0.01
            %             text(tt(n),0.2,'**')
            %         elseif p(n) < 0.05
            %             text(tt(n),0.2,'*')
        end
    end
end

title('Number of Spikes when ERSP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

subplot(2,1,2);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))],...
    [upper_errnERSP(t(edges_sel)<-tstep) flip(lower_errnERSP(t(edges_sel)<-tstep))],...
    cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errnERSP(t(edges_sel)>tstep) flip(lower_errnERSP(t(edges_sel)>tstep))], ...
    cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5)
plot(t(edges_sel),avNumIEDnERSP_pre*ones(size(edges_sel)),'-.','Color',[180/256, 180/256, 180/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDnERSP(t(edges_sel)<-tstep),'Color',cmap(20,:),'LineWidth',2)
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDnERSP(t(edges_sel)>tstep),'Color',cmap(20,:),'LineWidth',2)

ylim([ymin,ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 2; %nERSP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
        end
    end
end

title('Number of Spikes when no ERSP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

h = gcf;
h.Units = 'normalized';
h.Position = [0.1 -0.2 0.4 0.5];

figureName = sprintf('%s/fig6_IEDsERSP_%s',...
    myDataPath.Figures,s);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear avNumIEDERSP avNumIEDERSP_pre avNumIEDERSP_temp ERSP_pre nERSP_pre nn
clear avNumIEDnERSP avNumIEDnERSP_pre avNumIEDnERSP_temp num p p_ind p_sort p_sig s
clear countERSP_cont countERSP_disc countERSP_sel n n_poststim
clear countnERSP_cont countnERSP_disc countnERSP_sel m metaFileERSP metaFilenERSP
clear edges edges_cont edges_disc edges_sel h kk lower_errERSP lower_errnERSP
clear spikessampERSP spikessampnERSP sterrNumIEDERSP sterrNumIEDnERSP thisVal upper_errnERSP
clear upper_errERSP y ymax ymin figureName

%% add all spike samples into four categories: CCEPERSP, CCEPnERSP, nCCEPERSP or nCCEPnERSP

% first sort in order of ascending spike ratio
[~,I] = sort(log(all_spikeratio),'ascend');

spikessampCCEPERSP = []; % CCEP and ERSP
spikessampCCEPnERSP = []; % CCEP and no ERSP
spikessampnCCEPERSP = []; % no CCEP but ERSP
spikessampnCCEPnERSP = []; % no CCEP and no ERSP
metaFileCCEPERSP = []; metaFilenCCEPERSP = [];
metaFileCCEPnERSP = []; metaFilenCCEPnERSP = [];

for n = 1:size(all_CCEPmat,1)
    num = I(n);
    if all_CCEPmat(num) == 1 && all_ERSPmat(num) == 1
        spikessampCCEPERSP = [spikessampCCEPERSP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFileCCEPERSP = [metaFileCCEPERSP; all_metaFile(num,:)]; %#ok<AGROW>

    elseif all_CCEPmat(num) ==1 && all_ERSPmat(num) == 0
        spikessampCCEPnERSP = [spikessampCCEPnERSP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFileCCEPnERSP = [metaFileCCEPnERSP; all_metaFile(num,:)]; %#ok<AGROW>

    elseif all_CCEPmat(num) == 0 && all_ERSPmat(num) == 1
        spikessampnCCEPERSP = [spikessampnCCEPERSP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFilenCCEPERSP = [metaFilenCCEPERSP; all_metaFile(num,:)]; %#ok<AGROW>

    elseif all_CCEPmat(num) == 0 && all_ERSPmat(num) == 0
        spikessampnCCEPnERSP = [spikessampnCCEPnERSP; all_spikessamp(num,:)]; %#ok<AGROW>
        metaFilenCCEPnERSP = [metaFilenCCEPnERSP; all_metaFile(num,:)]; %#ok<AGROW>

    end
end

% housekeeping
clear I num n

% compare number of spikes per time window for CCEPERSP, CCEPnERSP, nCCEPERSP and nCCEPnERSP channels

tstep = 0.2;

% pre-allocation
countCCEPERSP_cont = NaN(size(spikessampCCEPERSP,1),cfg.epoch_length/0.2);
countnCCEPERSP_cont = NaN(size(spikessampnCCEPERSP,1),cfg.epoch_length/0.2);
countCCEPnERSP_cont = NaN(size(spikessampCCEPnERSP,1),cfg.epoch_length/0.2);
countnCCEPnERSP_cont = NaN(size(spikessampnCCEPnERSP,1),cfg.epoch_length/0.2);

% count number of spikes per time window (tstep)
for num = 1:size(spikessampCCEPERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampCCEPERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countCCEPERSP_cont(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnCCEPERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampnCCEPERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnCCEPERSP_cont(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampCCEPnERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampCCEPnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countCCEPnERSP_cont(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnCCEPnERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampnCCEPnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnCCEPnERSP_cont(num,1:size(n,2)) = n/10;
end

% convert edges to mean between two edges
edges_cont = [];
for nn = 1:size(edges,2)-1
    edges_cont(nn)= round(edges(nn)+((edges(nn+1)-edges(nn))/2)); %#ok<SAGROW>
end

% use edges and plot counts as block signal
edges_disc = []; countnCCEPERSP_disc = []; countCCEPERSP_disc = []; 
countnCCEPnERSP_disc = []; countCCEPnERSP_disc = [];
for nn = 1:size(edges,2)-1
    edges_disc(:,nn) = [round(edges(nn))+1,round(edges(nn+1))-1]; %#ok<SAGROW>
    countnCCEPERSP_disc(:,:,nn) = repmat(countnCCEPERSP_cont(:,nn),1,2); %#ok<SAGROW>
    countnCCEPnERSP_disc(:,:,nn) = repmat(countnCCEPnERSP_cont(:,nn),1,2); %#ok<SAGROW>
    countCCEPERSP_disc(:,:,nn) = repmat(countCCEPERSP_cont(:,nn),1,2); %#ok<SAGROW>
    countCCEPnERSP_disc(:,:,nn) = repmat(countCCEPnERSP_cont(:,nn),1,2); %#ok<SAGROW>
end

edges_disc = reshape(edges_disc,1,2*(size(edges,2)-1));
countnCCEPERSP_disc = reshape(countnCCEPERSP_disc, ...
    size(countnCCEPERSP_disc,1),2*size(countnCCEPERSP_cont,2));
countCCEPERSP_disc = reshape(countCCEPERSP_disc, ...
    size(countCCEPERSP_disc,1),2*size(countCCEPERSP_cont,2));
countnCCEPnERSP_disc = reshape(countnCCEPnERSP_disc, ...
    size(countnCCEPnERSP_disc,1),2*size(countnCCEPnERSP_cont,2));
countCCEPnERSP_disc = reshape(countCCEPnERSP_disc, ...
    size(countCCEPnERSP_disc,1),2*size(countCCEPnERSP_cont,2));

s = input('Do you want to display a block signal or a continuous signal? [block/continuous]: ','s');

if strcmpi(s,'continuous')
    edges_sel = edges_cont;
    countnCCEPERSP_sel = countnCCEPERSP_cont;
    countCCEPERSP_sel = countCCEPERSP_cont;
    countnCCEPnERSP_sel = countnCCEPnERSP_cont;
    countCCEPnERSP_sel = countCCEPnERSP_cont;

    % OR ...
elseif strcmpi(s,'block')
    edges_sel = edges_disc;
    countnCCEPERSP_sel = countnCCEPERSP_disc;
    countCCEPERSP_sel = countCCEPERSP_disc;
    countnCCEPnERSP_sel = countnCCEPnERSP_disc;
    countCCEPnERSP_sel = countCCEPnERSP_disc;

else
    error('The answer to the previous question was not "block" or "continuous".')
end

% number of IED when CCEP & ERSP
avNumIEDCCEPERSP_temp = mean(countCCEPERSP_sel);
avNumIEDCCEPERSP_pre = mean(avNumIEDCCEPERSP_temp(t(edges_sel)<-tstep));
avNumIEDCCEPERSP = avNumIEDCCEPERSP_temp;% - avNumIEDCCEPERSP_pre;
sterrNumIEDCCEPERSP = std(countCCEPERSP_sel)./sqrt(size(countCCEPERSP_sel,1));
lower_errCCEPERSP = avNumIEDCCEPERSP - sterrNumIEDCCEPERSP;
upper_errCCEPERSP = avNumIEDCCEPERSP + sterrNumIEDCCEPERSP;

% number of IED when no CCEP but ERSP
avNumIEDnCCEPERSP_temp = mean(countnCCEPERSP_sel);
avNumIEDnCCEPERSP_pre = mean(avNumIEDnCCEPERSP_temp(t(edges_sel)<-tstep));
avNumIEDnCCEPERSP = avNumIEDnCCEPERSP_temp;% - avNumIEDnCCEPERSP_pre;
sterrNumIEDnCCEPERSP = std(countnCCEPERSP_sel)./sqrt(size(countnCCEPERSP_sel,1));
lower_errnCCEPERSP = avNumIEDnCCEPERSP - sterrNumIEDnCCEPERSP;
upper_errnCCEPERSP = avNumIEDnCCEPERSP + sterrNumIEDnCCEPERSP;

% number of IED when CCEP but no ERSP
avNumIEDCCEPnERSP_temp = mean(countCCEPnERSP_sel);
avNumIEDCCEPnERSP_pre = mean(avNumIEDCCEPnERSP_temp(t(edges_sel)<-tstep));
avNumIEDCCEPnERSP = avNumIEDCCEPnERSP_temp;% - avNumIEDCCEPnERSP_pre;
sterrNumIEDCCEPnERSP = std(countCCEPnERSP_sel)./sqrt(size(countCCEPnERSP_sel,1));
lower_errCCEPnERSP = avNumIEDCCEPnERSP - sterrNumIEDCCEPnERSP;
upper_errCCEPnERSP = avNumIEDCCEPnERSP + sterrNumIEDCCEPnERSP;

% number of IED when no CCEP & no ERSP
avNumIEDnCCEPnERSP_temp = mean(countnCCEPnERSP_sel);
avNumIEDnCCEPnERSP_pre = mean(avNumIEDnCCEPnERSP_temp(t(edges_sel)<-tstep));
avNumIEDnCCEPnERSP = avNumIEDnCCEPnERSP_temp;% - avNumIEDnCCEPnERSP_pre;
sterrNumIEDnCCEPnERSP = std(countnCCEPnERSP_sel)./sqrt(size(countnCCEPnERSP_sel,1));
lower_errnCCEPnERSP = avNumIEDnCCEPnERSP - sterrNumIEDnCCEPnERSP;
upper_errnCCEPnERSP = avNumIEDnCCEPnERSP + sterrNumIEDnCCEPnERSP;

% FDR correction
n_poststim = find(t(edges_cont)>tstep,1,'first'); p = [];
for n = n_poststim:size(countCCEPERSP_cont,2)
    CCEPERSP_pre = countCCEPERSP_cont(:,t(edges_cont)<-tstep);
    nCCEPERSP_pre = countnCCEPERSP_cont(:,t(edges_cont)<-tstep);
    CCEPnERSP_pre = countCCEPnERSP_cont(:,t(edges_cont)<-tstep);
    nCCEPnERSP_pre = countnCCEPnERSP_cont(:,t(edges_cont)<-tstep);
    p(1,n-n_poststim+1) = ranksum(countCCEPERSP_cont(:,n),CCEPERSP_pre(:)); %#ok<SAGROW> 
    p(2,n-n_poststim+1) = ranksum(countnCCEPERSP_cont(:,n),nCCEPERSP_pre(:)); %#ok<SAGROW> 
    p(3,n-n_poststim+1) = ranksum(countCCEPnERSP_cont(:,n),CCEPnERSP_pre(:)); %#ok<SAGROW> 
    p(4,n-n_poststim+1) = ranksum(countnCCEPnERSP_cont(:,n),nCCEPnERSP_pre(:)); %#ok<SAGROW> 
end

m = length(p(:));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

% plot figures
ymin = round(0.9*min([avNumIEDnCCEPERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDCCEPERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDnCCEPnERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDnCCEPnERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep)]),2,'significant');
ymax = round(1.1*max([avNumIEDnCCEPERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDCCEPERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDnCCEPnERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep), ...
    avNumIEDCCEPnERSP(t(edges_sel)<-tstep| t(edges_sel)>tstep)]),2,'significant');

figure,
subplot(4,1,1);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))], ...
    [upper_errCCEPERSP(t(edges_sel)<-tstep), flip(lower_errCCEPERSP(t(edges_sel)<-tstep))], ...
    cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errCCEPERSP(t(edges_sel)>tstep) flip(lower_errCCEPERSP(t(edges_sel)>tstep))], ...
    cmap(10,:),'EdgeColor',cmap(10,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
plot(t(edges_sel),avNumIEDCCEPERSP_pre*ones(size(edges_sel)),'-.','Color',[200/256, 200/256, 200/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDCCEPERSP(t(edges_sel)<-tstep),'Color',cmap(10,:),'LineWidth',2);
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDCCEPERSP(t(edges_sel)>tstep),'Color',cmap(10,:),'LineWidth',2)

ylim([ymin, ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 1; %CCEPERSP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
            %         elseif p(n) <0.01
            %             text(tt(n),0.2,'**')
            %         elseif p(n) < 0.05
            %             text(tt(n),0.2,'*')
        end
    end
end

title('Number of Spikes when CCEP and ERSP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

subplot(4,1,2);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))], ...
    [upper_errnCCEPERSP(t(edges_sel)<-tstep) flip(lower_errnCCEPERSP(t(edges_sel)<-tstep))], ...
    cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errnCCEPERSP(t(edges_sel)>tstep) flip(lower_errnCCEPERSP(t(edges_sel)>tstep))], ...
    cmap(20,:),'EdgeColor',cmap(20,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
plot(t(edges_sel),avNumIEDnCCEPERSP_pre*ones(size(edges_sel)),'-.','Color',[180/256, 180/256, 180/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDnCCEPERSP(t(edges_sel)<-tstep),'Color',cmap(20,:),'LineWidth',2);
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDnCCEPERSP(t(edges_sel)>tstep),'Color',cmap(20,:),'LineWidth',2)

ylim([ymin, ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 2; %nCCEPERSP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
        end
    end
end

title('Number of Spikes when no CCEP but ERSP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

subplot(4,1,3);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))], ...
    [upper_errCCEPnERSP(t(edges_sel)<-tstep) flip(lower_errCCEPnERSP(t(edges_sel)<-tstep))], ...
    cmap(30,:),'EdgeColor',cmap(30,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errCCEPnERSP(t(edges_sel)>tstep) flip(lower_errCCEPnERSP(t(edges_sel)>tstep))], ...
    cmap(30,:),'EdgeColor',cmap(30,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
plot(t(edges_sel),avNumIEDCCEPnERSP_pre*ones(size(edges_sel)),'-.','Color',[180/256, 180/256, 180/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDCCEPnERSP(t(edges_sel)<-tstep),'Color',cmap(30,:),'LineWidth',2);
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDCCEPnERSP(t(edges_sel)>tstep),'Color',cmap(30,:),'LineWidth',2)

ylim([ymin,ymax])

% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 3; %CCEPnERSP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
        end
    end
end

title('Number of Spikes when CCEP but no ERSP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

subplot(4,1,4);
fill([t(edges_sel(t(edges_sel)<-tstep)), flip(t(edges_sel(t(edges_sel)<-tstep)))], ...
    [upper_errnCCEPnERSP(t(edges_sel)<-tstep) flip(lower_errnCCEPnERSP(t(edges_sel)<-tstep))], ...
    cmap(40,:),'EdgeColor',cmap(40,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
hold on
fill([t(edges_sel(t(edges_sel)>tstep)), flip(t(edges_sel(t(edges_sel)>tstep)))], ...
    [upper_errnCCEPnERSP(t(edges_sel)>tstep) flip(lower_errnCCEPnERSP(t(edges_sel)>tstep))], ...
    cmap(40,:),'EdgeColor',cmap(40,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
plot(t(edges_sel),avNumIEDnCCEPnERSP_pre*ones(size(edges_sel)),'-.','Color',[180/256, 180/256, 180/256])
plot(t(edges_sel(t(edges_sel)<-tstep)),avNumIEDnCCEPnERSP(t(edges_sel)<-tstep),'Color',cmap(40,:),'LineWidth',2);
plot(t(edges_sel(t(edges_sel)>tstep)),avNumIEDnCCEPnERSP(t(edges_sel)>tstep),'Color',cmap(40,:),'LineWidth',2)

ylim([ymin,ymax])
% add p<0.05 with FDR correction
y = 0.95*ymax;
for n=1:size(p_sig,2) % number of time windows (9)
    m = 4; %nCCEPnERSP
    if p_sig(m,n) == 1 % significant
        if p(m,n)<0.05
            text(t(edges_cont(n_poststim+n-1)),y,'*','Color',cmap(m*10,:));
        end
    end
end

title('Number of Spikes when no CCEP and no ERSP')
xlabel('Time (s)')
ylabel('Mean number of spikes')

h = gcf;
h.Units = 'normalized';
h.Position = [0.5 -0.2 0.4 1.1];

% save figure
figureName = sprintf('%s/fig6_IEDsCCEPERSP_%s',...
    myDataPath.Figures,s);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-vector','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)

% housekeeping
clear avNumIEDCCEPERSP avNumIEDCCEPERSP_pre avNumIEDCCEPERSP_temp CCEPERSP_pre nCCEPERSP_pre nn
clear avNumIEDCCEPnERSP avNumIEDCCEPnERSP_pre avNumIEDCCEPnERSP_temp CCEPnERSP_pre nCCEPnERSP_pre 
clear avNumIEDnCCEPERSP avNumIEDnCCEPERSP_pre avNumIEDnCCEPERSP_temp num p p_ind p_sort p_sig s
clear avNumIEDnCCEPnERSP avNumIEDnCCEPnERSP_pre avNumIEDnCCEPnERSP_temp 
clear countCCEPERSP_cont countCCEPERSP_disc countCCEPERSP_sel n n_poststim
clear countCCEPnERSP_cont countCCEPnERSP_disc countCCEPnERSP_sel 
clear countnCCEPERSP_cont countnCCEPERSP_disc countnCCEPERSP_sel m metaFileCCEPERSP metaFilenCCEPERSP
clear countnCCEPnERSP_cont countnCCEPnERSP_disc countnCCEPnERSP_sel metaFileCCEPnERSP metaFilenCCEPnERSP
clear edges edges_cont edges_disc edges_sel h kk lower_errCCEPERSP lower_errnCCEPERSP lower_errCCEPnERSP lower_errnCCEPnERSP
clear spikessampCCEPERSP spikessampnCCEPERSP sterrNumIEDCCEPERSP sterrNumIEDnCCEPERSP thisVal upper_errnCCEPERSP
clear spikessampCCEPnERSP spikessampnCCEPnERSP sterrNumIEDCCEPnERSP sterrNumIEDnCCEPnERSP upper_errnCCEPnERSP
clear upper_errCCEPERSP y ymax ymin 
clear upper_errCCEPnERSP figureName



%% figure with all IEDs as dots

countCCEP = 1;

figure(1),
for n = 1:size(spikessampCCEP,1)
    disp(n)
    %     for m = 1:size(spikessampCCEP,2)
    hold on
    samps = vertcat(spikessampCCEP{n,:});
    if ~isempty(samps)
        y = ones(size(samps));
        plot(t(samps),countCCEP*y,'.','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',1)
        countCCEP = countCCEP+1;
    end
    %     end
end

title('Spikes pre and post stimulation when CCEP')

countCCEP = 1;

figure(2),
for n = 1:size(spikessampnCCEP,1)
    disp(n)
    %     for m = 1:size(spikessampnCCEP,2)
    hold on
    samps = vertcat(spikessampnCCEP{n,:});
    if ~isempty(samps)
        y = ones(size(samps));
        plot(t(samps),countCCEP*y,'.','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',1)
        countCCEP = countCCEP+1;
    end
    %     end
end

title('Spikes pre and post stimulation when no CCEP')

%% a check to see whether the spikes make sense

num = 1; %size(metaFileCCEP,1);
subj = metaFileCCEP(num,1);
idx_stim = find(dataBase(subj).cc_stimsets(:,1)==metaFileCCEP(num,2) & dataBase(subj).cc_stimsets(:,2)==metaFileCCEP(num,3));
ch = metaFileCCEP(num,4);

figure,
subplot(2,1,1),
for n=1:10
    plot(-1*squeeze(dataBase(subj).cc_epoch_sorted(ch,n,idx_stim,:))+n*1000,'b')
    hold on
    plot(spikessampCCEP{num,n},...
        -1*squeeze(dataBase(subj).cc_epoch_sorted(ch,n,idx_stim,spikessampCCEP{num,n}))+n*1000,'k*')
end
xlim([0,size(dataBase(subj).cc_epoch_sorted,4)])
subplot(2,1,2),
samps = vertcat(spikessampCCEP{num,:});
plot(samps,ones(size(samps)),'*')
xlim([0,size(dataBase(subj).cc_epoch_sorted,4)])


%% end