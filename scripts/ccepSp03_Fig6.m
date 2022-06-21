% Figure 5: spikes in epochs (not) connected and (no) suppressed power

close all
clc

cfg = dataBase(1).ERSP.cfg;
fs = dataBase(1).ccep_header.Fs;

%% combine all

ERSPmatAll = [];
CCEPmatAll = [];
spikessampAll = [];
spikeratioAll = [];
metaFileAll = [];
for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)
        ERSPmat = dataBase(subj).ERSPmat(:,dataBase(subj).IEDs.IEDch);
        ERSPmatAll = [ERSPmatAll; ERSPmat(:)]; %#ok<AGROW> % with ERSPmat(:), all columns will be put below each other, so first column 1, then below this column 2 etc.
        CCEPmat =  dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);
        CCEPmatAll = [CCEPmatAll; CCEPmat(:)]; %#ok<AGROW>
        spikessamp = dataBase(subj).IEDs.totEpochspikessamp;
        [m,n,o]=size(spikessamp);
        spikessamp_reshape = reshape(spikessamp,m*n,o);
        spikessampAll = [spikessampAll; spikessamp_reshape]; %#ok<AGROW>
        spikeratioAll = [spikeratioAll; dataBase(subj).IEDs.spikesratio(:)]; %#ok<AGROW>

        metaFile = [subj*ones(size(ERSPmat,1)*size(ERSPmat,2),1),...
            repmat(dataBase(subj).cc_stimsets,size(ERSPmat,2),1),...
            reshape(repmat(dataBase(subj).IEDs.IEDch',size(ERSPmat,1),1),m*n,1)];
        metaFileAll = [metaFileAll; metaFile]; %#ok<AGROW>
    end
end

%%

[~,I] = sort(log(spikeratioAll),'ascend');

spikessampCCEP = [];
spikessampnCCEP = [];
metaFileCCEP = []; metaFilenCCEP = [];

for n = 1:size(CCEPmatAll,1)
    num = I(n);
    if CCEPmatAll(num) == 1
        spikessampCCEP = [spikessampCCEP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFileCCEP = [metaFileCCEP; metaFileAll(num,:)]; %#ok<AGROW>

    elseif CCEPmatAll(num) == 0
        spikessampnCCEP = [spikessampnCCEP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFilenCCEP = [metaFilenCCEP; metaFileAll(num,:)]; %#ok<AGROW>
    end
end

% %% plot with all spikes as dots
% fs = dataBase(subj).ccep_header.Fs;
% t = (1/fs:1/fs:cfg.epoch_length)-cfg.epoch_prestim;
% countCCEP = 1;
%
% figure(1),
% for n = 1:size(spikessampCCEP,1)
%     n
%     %     for m = 1:size(spikessampCCEP,2)
%     hold on
%     samps = vertcat(spikessampCCEP{n,:});
%     if ~isempty(samps)
%         y = ones(size(samps));
%         plot(t(samps),countCCEP*y,'.','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',1)
%         countCCEP = countCCEP+1;
%     end
%     %     end
% end
%
% title('Spikes pre and post stimulation when CCEP')
%
% countCCEP = 1;
%
% figure(2),
% for n = 1:size(spikessampnCCEP,1)
%     n
%     %     for m = 1:size(spikessampnCCEP,2)
%     hold on
%     samps = vertcat(spikessampnCCEP{n,:});
%     if ~isempty(samps)
%         y = ones(size(samps));
%         plot(t(samps),countCCEP*y,'.','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',1)
%         countCCEP = countCCEP+1;
%     end
%     %     end
% end
%
% title('Spikes pre and post stimulation when no CCEP')

%% a check to see whether the spikes make sense

% %num = 1; %size(metaFileCCEP,1);
% subj = metaFileCCEP(num,1);
% idx_stim = find(dataBase(subj).cc_stimsets(:,1)==metaFileCCEP(num,2) & dataBase(subj).cc_stimsets(:,2)==metaFileCCEP(num,3));
% ch = metaFileCCEP(num,4);
%
% figure,
% subplot(2,1,1),
% for n=1:10
%     plot(-1*squeeze(dataBase(subj).cc_epoch_sorted(ch,n,idx_stim,:))+n*1000,'b')
%     hold on
%     plot(spikessampCCEP{num,n},...
%         -1*squeeze(dataBase(subj).cc_epoch_sorted(ch,n,idx_stim,spikessampCCEP{num,n}))+n*1000,'k*')
% end
% xlim([0,size(dataBase(subj).cc_epoch_sorted,4)])
% subplot(2,1,2),
% samps = vertcat(spikessampCCEP{num,:});
% plot(samps,ones(size(samps)),'*')
% xlim([0,size(dataBase(subj).cc_epoch_sorted,4)])
%
%% summarize all spikes into a nice plot

tstep = 0.2;
tt = (tstep/2:tstep:cfg.epoch_length) - cfg.epoch_prestim;

countCCEP = NaN(size(spikessampCCEP,1),cfg.epoch_length/0.2);
countnCCEP = NaN(size(spikessampnCCEP,1),cfg.epoch_length/0.2);

for num = 1:size(spikessampCCEP,1)
    [n,edges ] = histcounts(vertcat(spikessampCCEP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countCCEP(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnCCEP,1)
    [n,edges ] = histcounts(vertcat(spikessampnCCEP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnCCEP(num,1:size(n,2)) = n/10;
end

avNumIEDCCEP_temp = mean(countCCEP);
avNumIEDCCEP_pre = mean(avNumIEDCCEP_temp(tt<-0.15));
avNumIEDCCEP = avNumIEDCCEP_temp - avNumIEDCCEP_pre;
sterrNumIEDCCEP = std(countCCEP)./sqrt(size(countCCEP,1));
lower_errCCEP = avNumIEDCCEP - sterrNumIEDCCEP;
upper_errCCEP = avNumIEDCCEP + sterrNumIEDCCEP;

avNumIEDnCCEP_temp = mean(countnCCEP);
avNumIEDnCCEP_pre = mean(avNumIEDnCCEP_temp(tt<-0.15));
avNumIEDnCCEP = avNumIEDnCCEP_temp - avNumIEDnCCEP_pre;
sterrNumIEDnCCEP = std(countnCCEP)./sqrt(size(countnCCEP,1));
lower_errnCCEP = avNumIEDnCCEP - sterrNumIEDnCCEP;
upper_errnCCEP = avNumIEDnCCEP + sterrNumIEDnCCEP;

cmap = parula(4);

figure,
fill([tt(tt<-0.15), flip(tt(tt<-0.15))],...
    [upper_errCCEP(tt<-0.15) flip(lower_errCCEP(tt<-0.15))],...
    cmap(1,:),'EdgeColor',cmap(1,:),'FaceAlpha',0.5)
hold on
fill([tt(tt>0.15), flip(tt(tt>0.15))], ...
    [upper_errCCEP(tt>0.15) flip(lower_errCCEP(tt>0.15))], ...
    cmap(1,:),'EdgeColor',cmap(1,:),'FaceAlpha',0.5)
plot(tt(tt<-0.15),avNumIEDCCEP(tt<-0.15),'Color',cmap(1,:),'LineWidth',2)
plot(tt(tt>0.15),avNumIEDCCEP(tt>0.15),'Color',cmap(1,:),'LineWidth',2)

fill([tt(tt<-0.15), flip(tt(tt<-0.15))], ...
    [upper_errnCCEP(tt<-0.15) flip(lower_errnCCEP(tt<-0.15))], ...
    cmap(2,:),'EdgeColor',cmap(2,:),'FaceAlpha',0.5)
fill([tt(tt>0.15), flip(tt(tt>0.15))], ...
    [upper_errnCCEP(tt>0.15) flip(lower_errnCCEP(tt>0.15))], ...
    cmap(2,:),'EdgeColor',cmap(2,:),'FaceAlpha',0.5)
plot(tt(tt<-0.15),avNumIEDnCCEP(tt<-0.15),'Color',cmap(2,:),'LineWidth',2)
plot(tt(tt>0.15),avNumIEDnCCEP(tt>0.15),'Color',cmap(2,:),'LineWidth',2)
% ylim([0 0.25])

% figure,
% fill([tt, flip(tt)],[upper_errCCEP flip(lower_errCCEP)],'b','EdgeColor','b')
% hold on
% fill([tt, flip(tt)],[upper_errCCEP flip(lower_errCCEP)],'b','EdgeColor','b')
% plot(tt,avNumIEDCCEP,'w')
% plot(tt,avNumIEDCCEP,'w')
%
% fill([tt, flip(tt)],[upper_errnCCEP flip(lower_errnCCEP)],'g','EdgeColor','g')
% hold on
% fill([tt, flip(tt)],[upper_errnCCEP flip(lower_errnCCEP)],'g','EdgeColor','g')
% plot(tt,avNumIEDnCCEP,'w')
% plot(tt,avNumIEDnCCEP,'w')
%
% FDR correction
n_poststim = find(tt>0.15,1,'first');
for n=n_poststim:size(countCCEP,2)
    CCEP_pre = countCCEP(:,tt<-0.15);
    nCCEP_pre = countnCCEP(:,tt<-0.15);
    p(1,n-n_poststim+1) = ranksum(countCCEP(:,n),CCEP_pre(:));
    p(2,n-n_poststim+1) = ranksum(countnCCEP(:,n),nCCEP_pre(:));
end

% FDR correction
m = length(p(:));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

y = round(1.3*max([avNumIEDnCCEP, avNumIEDCCEP]),3,'significant');
for n=1:size(p_sig,2)
    for m = 1:size(p_sig,1)

        if p_sig(m,n) == 1 % significant
            if p(m,n)<0.01
                text(tt(n_poststim+n-1),y,'*','Color',cmap(m,:))
                %         elseif p(n) <0.01
                %             text(tt(n),0.2,'**')
                %         elseif p(n) < 0.05
                %             text(tt(n),0.2,'*')
            end
        end
    end
end

%% taking into account both ERSP and CCEP
%%

[~,I] = sort(log(spikeratioAll),'ascend');

spikessampCCEPERSP = []; % CCEP and ERSP
spikessampCCEPnERSP = []; % CCEP and no ERSP
spikessampnCCEPERSP = []; % no CCEP but ERSP
spikessampnCCEPnERSP = []; % no CCEP and no ERSP
metaFileCCEPERSP = []; metaFilenCCEPERSP = [];
metaFileCCEPnERSP = []; metaFilenCCEPnERSP = [];

for n = 1:size(CCEPmatAll,1)
    num = I(n);
    if CCEPmatAll(num) == 1 && ERSPmatAll(num) == 1
        spikessampCCEPERSP = [spikessampCCEPERSP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFileCCEPERSP = [metaFileCCEPERSP; metaFileAll(num,:)]; %#ok<AGROW>

    elseif CCEPmatAll(num) ==1 && ERSPmatAll(num) == 0
        spikessampCCEPnERSP = [spikessampCCEPnERSP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFileCCEPnERSP = [metaFileCCEPnERSP; metaFileAll(num,:)]; %#ok<AGROW>

    elseif CCEPmatAll(num) == 0 && ERSPmatAll(num) == 1
        spikessampnCCEPERSP = [spikessampnCCEPERSP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFilenCCEPERSP = [metaFilenCCEPERSP; metaFileAll(num,:)]; %#ok<AGROW>

    elseif CCEPmatAll(num) == 0 && ERSPmatAll(num) == 0
        spikessampnCCEPnERSP = [spikessampnCCEPnERSP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFilenCCEPnERSP = [metaFilenCCEPnERSP; metaFileAll(num,:)]; %#ok<AGROW>

    end
end

%% summarize all spikes into a nice plot

tstep = 0.2;
countCCEPERSP = NaN(size(spikessampCCEPERSP,1),cfg.epoch_length/0.2);
countnCCEPERSP = NaN(size(spikessampnCCEPERSP,1),cfg.epoch_length/0.2);
countCCEPnERSP = NaN(size(spikessampCCEPnERSP,1),cfg.epoch_length/0.2);
countnCCEPnERSP = NaN(size(spikessampnCCEPnERSP,1),cfg.epoch_length/0.2);

for num = 1:size(spikessampCCEPERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampCCEPERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countCCEPERSP(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnCCEPERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampnCCEPERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnCCEPERSP(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampCCEPnERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampCCEPnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countCCEPnERSP(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnCCEPnERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampnCCEPnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnCCEPnERSP(num,1:size(n,2)) = n/10;
end

avNumIEDCCEPERSP_temp = mean(countCCEPERSP);
avNumIEDCCEPERSP_pre = mean(avNumIEDCCEPERSP_temp(tt<-0.15));
avNumIEDCCEPERSP = avNumIEDCCEPERSP_temp - avNumIEDCCEPERSP_pre;
sterrNumIEDCCEPERSP = std(countCCEPERSP)./sqrt(size(countCCEPERSP,1));
lower_errCCEPERSP = avNumIEDCCEPERSP - sterrNumIEDCCEPERSP;
upper_errCCEPERSP = avNumIEDCCEPERSP + sterrNumIEDCCEPERSP;

avNumIEDnCCEPERSP_temp = mean(countnCCEPERSP);
avNumIEDnCCEPERSP_pre = mean(avNumIEDnCCEPERSP_temp(tt<-0.15));
avNumIEDnCCEPERSP = avNumIEDnCCEPERSP_temp - avNumIEDnCCEPERSP_pre;
sterrNumIEDnCCEPERSP = std(countnCCEPERSP)./sqrt(size(countnCCEPERSP,1));
lower_errnCCEPERSP = avNumIEDnCCEPERSP - sterrNumIEDnCCEPERSP;
upper_errnCCEPERSP = avNumIEDnCCEPERSP + sterrNumIEDnCCEPERSP;

avNumIEDCCEPnERSP_temp = mean(countCCEPnERSP);
avNumIEDCCEPnERSP_pre = mean(avNumIEDCCEPnERSP_temp(tt<-0.15));
avNumIEDCCEPnERSP = avNumIEDCCEPnERSP_temp - avNumIEDCCEPnERSP_pre;
sterrNumIEDCCEPnERSP = std(countCCEPnERSP)./sqrt(size(countCCEPnERSP,1));
lower_errCCEPnERSP = avNumIEDCCEPnERSP - sterrNumIEDCCEPnERSP;
upper_errCCEPnERSP = avNumIEDCCEPnERSP + sterrNumIEDCCEPnERSP;

avNumIEDnCCEPnERSP_temp = mean(countnCCEPnERSP);
avNumIEDnCCEPnERSP_pre = mean(avNumIEDnCCEPnERSP_temp(tt<-0.15));
avNumIEDnCCEPnERSP = avNumIEDnCCEPnERSP_temp - avNumIEDnCCEPnERSP_pre;
sterrNumIEDnCCEPnERSP = std(countnCCEPnERSP)./sqrt(size(countnCCEPnERSP,1));
lower_errnCCEPnERSP = avNumIEDnCCEPnERSP - sterrNumIEDnCCEPnERSP;
upper_errnCCEPnERSP = avNumIEDnCCEPnERSP + sterrNumIEDnCCEPnERSP;

tt = (tstep/2:tstep:cfg.epoch_length) - cfg.epoch_prestim;

cmap = parula(5);
figure,
fill([tt(tt<-0.15), flip(tt(tt<-0.15))], ...
    [upper_errCCEPERSP(tt<-0.15) flip(lower_errCCEPERSP(tt<-0.15))], ...
    cmap(1,:),'EdgeColor',cmap(1,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
hold on
fill([tt(tt>0.15), flip(tt(tt>0.15))], ...
    [upper_errCCEPERSP(tt>0.15) flip(lower_errCCEPERSP(tt>0.15))], ...
    cmap(1,:),'EdgeColor',cmap(1,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
h1 = plot(tt(tt<-0.15),avNumIEDCCEPERSP(tt<-0.15),'Color',cmap(1,:),'LineWidth',2);
plot(tt(tt>0.15),avNumIEDCCEPERSP(tt>0.15),'Color',cmap(1,:),'LineWidth',2)

fill([tt(tt<-0.15), flip(tt(tt<-0.15))], ...
    [upper_errnCCEPERSP(tt<-0.15) flip(lower_errnCCEPERSP(tt<-0.15))], ...
    cmap(2,:),'EdgeColor',cmap(2,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
fill([tt(tt>0.15), flip(tt(tt>0.15))], ...
    [upper_errnCCEPERSP(tt>0.15) flip(lower_errnCCEPERSP(tt>0.15))], ...
    cmap(2,:),'EdgeColor',cmap(2,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
h2 = plot(tt(tt<-0.15),avNumIEDnCCEPERSP(tt<-0.15),'Color',cmap(2,:),'LineWidth',2);
plot(tt(tt>0.15),avNumIEDnCCEPERSP(tt>0.15),'Color',cmap(2,:),'LineWidth',2)

fill([tt(tt<-0.15), flip(tt(tt<-0.15))], ...
    [upper_errCCEPnERSP(tt<-0.15) flip(lower_errCCEPnERSP(tt<-0.15))], ...
    cmap(3,:),'EdgeColor',cmap(3,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
fill([tt(tt>0.15), flip(tt(tt>0.15))], ...
    [upper_errCCEPnERSP(tt>0.15) flip(lower_errCCEPnERSP(tt>0.15))], ...
    cmap(3,:),'EdgeColor',cmap(3,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
h3 = plot(tt(tt<-0.15),avNumIEDCCEPnERSP(tt<-0.15),'Color',cmap(3,:),'LineWidth',2);
plot(tt(tt>0.15),avNumIEDCCEPnERSP(tt>0.15),'Color',cmap(3,:),'LineWidth',2)

fill([tt(tt<-0.15), flip(tt(tt<-0.15))], ...
    [upper_errnCCEPnERSP(tt<-0.15) flip(lower_errnCCEPnERSP(tt<-0.15))], ...
    cmap(4,:),'EdgeColor',cmap(4,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
fill([tt(tt>0.15), flip(tt(tt>0.15))], ...
    [upper_errnCCEPnERSP(tt>0.15) flip(lower_errnCCEPnERSP(tt>0.15))], ...
    cmap(4,:),'EdgeColor',cmap(4,:),'FaceAlpha',0.5,'EdgeAlpha',0.5)
h4 = plot(tt(tt<-0.15),avNumIEDnCCEPnERSP(tt<-0.15),'Color',cmap(4,:),'LineWidth',2);
plot(tt(tt>0.15),avNumIEDnCCEPnERSP(tt>0.15),'Color',cmap(4,:),'LineWidth',2)

legend([h1,h2,h3,h4],'CCEP & ERSP','nCCEP & ERSP','CCEP & nERSP','nCCEP & nERSP')

% FDR correction
n_poststim = find(tt>0.15,1,'first');
for n=n_poststim:size(countCCEPERSP,2)
    CCEPERSP_pre = countCCEPERSP(:,tt<-0.15);
    nCCEPERSP_pre = countnCCEPERSP(:,tt<-0.15);
    CCEPnERSP_pre = countCCEPnERSP(:,tt<-0.15);
    nCCEPnERSP_pre = countnCCEPnERSP(:,tt<-0.15);
    p(1,n-n_poststim+1) = ranksum(countCCEPERSP(:,n),CCEPERSP_pre(:));
    p(2,n-n_poststim+1) = ranksum(countnCCEPERSP(:,n),nCCEPERSP_pre(:));
    p(3,n-n_poststim+1) = ranksum(countCCEPnERSP(:,n),CCEPnERSP_pre(:));
    p(4,n-n_poststim+1) = ranksum(countnCCEPnERSP(:,n),nCCEPnERSP_pre(:));
end

% FDR correction
m = length(p(:));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

y = round(1.5*max([avNumIEDnCCEP, avNumIEDCCEP]),2);
for n=1:size(p_sig,2)
    for m = 1:size(p_sig,1)

        if p_sig(m,n) == 1 % significant
            if p(m,n)<0.05
                text(tt(n_poststim+n-1),y+m*0.001,'*','Color',cmap(m,:))
                %         elseif p(n) <0.01
                %             text(tt(n),0.2,'**')
                %         elseif p(n) < 0.05
                %             text(tt(n),0.2,'*')
            end
        end
    end
end


%% looking into only ERSP
%%

[~,I] = sort(log(spikeratioAll),'ascend');

spikessampERSP = []; % CCEP and ERSP
spikessampnERSP = []; % CCEP and no ERSP
metaFileERSP = []; metaFilenERSP = [];

for n = 1:size(ERSPmatAll,1)
    num = I(n);
    if ERSPmatAll(num) == 1
        spikessampERSP = [spikessampERSP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFileERSP = [metaFileERSP; metaFileAll(num,:)]; %#ok<AGROW>

    elseif ERSPmatAll(num) == 0
        spikessampnERSP = [spikessampnERSP; spikessampAll(num,:)]; %#ok<AGROW>
        metaFilenERSP = [metaFilenERSP; metaFileAll(num,:)]; %#ok<AGROW>
    end
end

%% summarize all spikes into a nice plot

tstep = 0.2;
countERSP = NaN(size(spikessampERSP,1),cfg.epoch_length/0.2);
countnERSP = NaN(size(spikessampnERSP,1),cfg.epoch_length/0.2);

for num = 1:size(spikessampERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countERSP(num,1:size(n,2)) = n/10;
end

for num = 1:size(spikessampnERSP,1)
    [n,edges ] = histcounts(vertcat(spikessampnERSP{num,:}),'BinWidth',tstep*fs,'BinLimits',[0, cfg.epoch_length*fs]);

    countnERSP(num,1:size(n,2)) = n/10;
end

avNumIEDERSP = mean(countERSP);
sterrNumIEDERSP = std(countERSP)./sqrt(size(countERSP,1));
lower_errERSP = avNumIEDERSP - sterrNumIEDERSP;
upper_errERSP = avNumIEDERSP + sterrNumIEDERSP;

avNumIEDnERSP = mean(countnERSP);
sterrNumIEDnERSP = std(countnERSP)./sqrt(size(countnERSP,1));
lower_errnERSP = avNumIEDnERSP - sterrNumIEDnERSP;
upper_errnERSP = avNumIEDnERSP + sterrNumIEDnERSP;

tt = (tstep/2:tstep:cfg.epoch_length) - cfg.epoch_prestim;

figure,
fill([tt(tt<-0.15), flip(tt(tt<-0.15))],[upper_errERSP(tt<-0.15) flip(lower_errERSP(tt<-0.15))],'b','EdgeColor','b')
hold on
fill([tt(tt>0.15), flip(tt(tt>0.15))],[upper_errERSP(tt>0.15) flip(lower_errERSP(tt>0.15))],'b','EdgeColor','b')
plot(tt(tt<-0.15),avNumIEDERSP(tt<-0.15),'w')
plot(tt(tt>0.15),avNumIEDERSP(tt>0.15),'w')

fill([tt(tt<-0.15), flip(tt(tt<-0.15))],[upper_errnERSP(tt<-0.15) flip(lower_errnERSP(tt<-0.15))],'g','EdgeColor','g')
fill([tt(tt>0.15), flip(tt(tt>0.15))],[upper_errnERSP(tt>0.15) flip(lower_errnERSP(tt>0.15))],'g','EdgeColor','g')
plot(tt(tt<-0.15),avNumIEDnERSP(tt<-0.15),'w')
plot(tt(tt>0.15),avNumIEDnERSP(tt>0.15),'w')

% FDR correction
for n=1:size(countERSP,2)
    p(n) = ranksum(countERSP(:,n),countnCCEPERSP(:,n));

end

p(tt>=-0.15&tt<=0.15) = NaN;

% FDR correction
m = length(p(~isnan(p)));
[p_sort,p_ind] = sort(p(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end

p_sig = p;
p_sig(p_ind) = p_sort<thisVal;

for n=1:size(p_sig,2)
    if p_sig(n) == 1 % significant
        if p(n)<0.001
            text(tt(n),0.2,'*')
            %         elseif p(n) <0.01
            %             text(tt(n),0.2,'**')
            %         elseif p(n) < 0.05
            %             text(tt(n),0.2,'*')
        end
    end
end



%% end