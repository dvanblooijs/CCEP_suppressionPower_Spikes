% CCEP_powerSuppression_ERs_spikes
% this code analyzes the relationship between ERs, spikes and
% powerSuppression.
% author: D van Blooijs
% date: April 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath('git_rep/SPES_SOZ/detectERs')
addpath(genpath('git_rep/eeglab/'))
addpath('git_rep/fieldtrip/')
ft_defaults


%% process the data to enable construction of TFSPES plots

preprocess_mainfile

%% sort ERs_BB to cc_stimsets

for subj = 1:size(dataBase,2)
    stimchans = reshape([dataBase(subj).ERs_BB(:).stimchans],2, size(dataBase(subj).ERs_BB,2))';
    if isequal(dataBase(subj).cc_stimsets,stimchans)
        for stimp = 1:size(dataBase(subj).cc_stimsets,1)
            
            dataBase(subj).cc_ERs(stimp).vis = dataBase(subj).ERs_BB(stimp).ERsvis;
            dataBase(subj).cc_BBs(stimp).det = dataBase(subj).ERs_BB(stimp).BB;
        end
        
    else
        for stimp = 1:size(dataBase(subj).cc_stimsets,1)
            stimnum = find(stimchans(:,1) == dataBase(subj).cc_stimsets(stimp,1) & stimchans(:,2) == dataBase(subj).cc_stimsets(stimp,2));
            
            dataBase(subj).cc_ERs(stimp).vis = dataBase(subj).ERs_BB(stimnum).ERsvis;
            dataBase(subj).cc_BBs(stimp).det = dataBase(subj).ERs_BB(stimnum).BB;
            
        end
    end
end

%% preprocessing ERs & BBs for crosstab

for subj = 1:size(dataBase,2)
    ERmat = zeros(size(dataBase(subj).ch,1),size(dataBase(subj).cc_ERs,2));
    BBmat = zeros(size(dataBase(subj).ch,1),size(dataBase(subj).cc_BBs,2));
    
    for stimp =1:size(dataBase(subj).cc_ERs,2)
        % the stimulated electrodes must be NaN
        ERmat(dataBase(subj).cc_stimsets(stimp,:),stimp) = NaN;
        % the electrodes with ER must be 1
        ERmat(dataBase(subj).cc_ERs(stimp).vis,stimp) = 1;
    end
    
    for stimp =1:size(dataBase(subj).cc_BBs,2)
        % the stimulated electrodes must be NaN
        BBmat(dataBase(subj).cc_stimsets(stimp,:),stimp) = NaN;
        % the electrodes with BB must be 1
        BBmat(dataBase(subj).cc_BBs(stimp).det,stimp) = 1;
    end
    
    dataBase(subj).ERmat = ERmat;
    dataBase(subj).BBmat = BBmat;
    
end


%% load spikes

cfg.spikesinput = '/Fridge/users/dorien/derivatives/BB_article/IEDs/';

for subj = [1:4,6,8,10] %1:size(dataBase,2)
    load([cfg.spikesinput, dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label, '_detIED.mat']);
    
    dataBase(subj).IEDchan = spikespat.IEDchan;
    dataBase(subj).IEDch = spikespat.IEDch;
    dataBase(subj).spikes = spikespat.spikesdet;
    
end

%% locate spikes per stimulus

for subj = 1:size(dataBase,2)
    
    % pre-allocation
    spikelocs = cell(size(dataBase(subj).spikes,2),size(dataBase(subj).tt_epoch_sorted,2),size(dataBase(subj).tt_epoch_sorted,1));
    spks1slocspre = cell(size(dataBase(subj).spikes,2),size(dataBase(subj).tt_epoch_sorted,2),size(dataBase(subj).tt_epoch_sorted,1));
    spks1slocspost = cell(size(dataBase(subj).spikes,2),size(dataBase(subj).tt_epoch_sorted,2),size(dataBase(subj).tt_epoch_sorted,1));
    spks1sprenum = NaN(size(dataBase(subj).spikes,2),size(dataBase(subj).tt_epoch_sorted,2),size(dataBase(subj).tt_epoch_sorted,1));
    spks1spostnum = NaN(size(dataBase(subj).spikes,2),size(dataBase(subj).tt_epoch_sorted,2),size(dataBase(subj).tt_epoch_sorted,1));
    
    for IEDchan = 1:size(dataBase(subj).spikes,2)
        allspikes = [dataBase(subj).spikes{IEDchan}];
        
        for stimp=1:size(dataBase(subj).tt_epoch_sorted,2)
            for n=1:size(dataBase(subj).tt_epoch_sorted,1)
                
                % find locations of spikes in cc_epoch_sorted
                locs = allspikes(ismember(allspikes,dataBase(subj).tt_epoch_sorted(n,stimp,:)));
                
                % epoch is now 4s: 2s pre and 2s post stimulation --> start of stimulation
                startstim = round((dataBase(subj).tt_epoch_sorted(n,stimp,1)+dataBase(subj).tt_epoch_sorted(n,stimp,end))/2)-1;
                
                % 10 ms pre and post stimulation, no spikes were detected,
                % so determine spikes 1.01 s pre and post stimulation
                start1spre = startstim - round(1.01*dataBase(subj).ccep_header.Fs);
                stop1spre = startstim - round(0.01*dataBase(subj).ccep_header.Fs);
                
                start1spost = startstim + round(0.01*dataBase(subj).ccep_header.Fs);
                stop1spost = startstim + round(1.01*dataBase(subj).ccep_header.Fs);
                
                spikes1spre = locs(ismember(locs,start1spre:stop1spre));
                spikes1spost = locs(ismember(locs,start1spost:stop1spost));
                
                
                spikelocs(IEDchan,stimp,n) = {locs};
                
                spks1slocspre(IEDchan,stimp,n) = {spikes1spre};
                spks1slocspost(IEDchan,stimp,n) = {spikes1spost};
                
                spks1sprenum(IEDchan,stimp,n) = numel(spikes1spre);
                spks1spostnum(IEDchan,stimp,n) = numel(spikes1spost);
                
            end
        end
    end
    
    dataBase(subj).spikelocs        = spikelocs;
    dataBase(subj).spks1slocspre    = spks1slocspre;
    dataBase(subj).spks1slocspost   = spks1slocspost;
    dataBase(subj).spks1sprenum     = spks1sprenum;
    dataBase(subj).spks1spostnum    = spks1spostnum;
end


%% make ratio spike change

for subj = 1:size(dataBase,2)
    spikeratio = NaN(size(dataBase(subj).IEDchan,2),size(dataBase(subj).tt_epoch_sorted,2));
    
    for IEDchan = 1:size(dataBase(subj).IEDchan,2)
        for stimp = 1:size(dataBase(subj).tt_epoch_sorted,2)
            
            % if no spikes occurred before, and no spikes occurred after,
            % than both are 0 and this leads to a NaN, although, in fact,
            % this should be 0 (no spike change)
            if all(squeeze(dataBase(subj).spks1spostnum(IEDchan,stimp,:)) == 0 & squeeze(dataBase(subj).spks1sprenum(IEDchan,stimp,:)) == 0)
                spikeratio(IEDchan,stimp)= 0;
            else
                spikeratio(IEDchan,stimp) = nanmedian((squeeze(dataBase(subj).spks1spostnum(IEDchan,stimp,:)) - squeeze(dataBase(subj).spks1sprenum(IEDchan,stimp,:)))./...
                    (squeeze(dataBase(subj).spks1spostnum(IEDchan,stimp,:)) + squeeze(dataBase(subj).spks1sprenum(IEDchan,stimp,:))));
            end
        end
    end
    
    dataBase(subj).spikeratio = spikeratio;
end

%% categorize into : 0=no change, 1=decrease, 2=increase
spikeratioall = [];

for subj = 1:size(dataBase,2)
   spikeratioall = [spikeratioall; dataBase(subj).spikeratio(:)];
end

threshall = quantile(spikeratioall,[0.25,0.5,0.75]);

thresh = [-1/7 0 1/7];

for subj = 1: size(dataBase,2)
%     thresh = quantile(dataBase(subj).spikeratio(:),[0.25,0.5,0.75]);
    SR_change = dataBase(subj).spikeratio;
    SR_inc = dataBase(subj).spikeratio;
    SR_dec = dataBase(subj).spikeratio;

    SR_change(dataBase(subj).spikeratio>thresh(1) & dataBase(subj).spikeratio<thresh(3)) = 0;
    SR_change(dataBase(subj).spikeratio<=thresh(1) | dataBase(subj).spikeratio>=thresh(3)) = 1;
    
    % increase in spike ratio
    SR_inc(dataBase(subj).spikeratio>thresh(1) & dataBase(subj).spikeratio<thresh(3)) = 0;
    SR_inc(dataBase(subj).spikeratio<=thresh(1) ) = NaN;
    SR_inc(dataBase(subj).spikeratio>=thresh(3)) = 1;
    
    % decrease in spike ratio
    SR_dec(dataBase(subj).spikeratio>thresh(1) & dataBase(subj).spikeratio<thresh(3)) = 0;
    SR_dec(dataBase(subj).spikeratio<=thresh(1) ) = 1;
    SR_dec(dataBase(subj).spikeratio>=thresh(3)) = NaN;
    
%     dataBase(subj).thresh = thresh;
    dataBase(subj).SR_change = SR_change;
    dataBase(subj).SR_inc = SR_inc;
    dataBase(subj).SR_dec = SR_dec;
end

disp('Calculated spike ratios')
% --> too much variation in distribution of spike ratios! Therefore,
% combine all spike ratios and then determine thresholds!

thresh = [0];

for subj = 1: size(dataBase,2)
    SR_change = dataBase(subj).spikeratio;

    SR_change(dataBase(subj).spikeratio>thresh) = 1;
    SR_change(dataBase(subj).spikeratio<=thresh) = -1;
    
    % increase in spike ratio
%     SR_inc(dataBase(subj).spikeratio>thresh(1) & dataBase(subj).spikeratio<thresh(3)) = 0;
%     SR_inc(dataBase(subj).spikeratio<=thresh(1) ) = NaN;
%     SR_inc(dataBase(subj).spikeratio>=thresh(3)) = 1;
%     
%     % decrease in spike ratio
%     SR_dec(dataBase(subj).spikeratio>thresh(1) & dataBase(subj).spikeratio<thresh(3)) = 0;
%     SR_dec(dataBase(subj).spikeratio<=thresh(1) ) = 1;
%     SR_dec(dataBase(subj).spikeratio>=thresh(3)) = NaN;
%     
% %     dataBase(subj).thresh = thresh;
    dataBase(subj).SR_change = SR_change;
%     dataBase(subj).SR_inc = SR_inc;
%     dataBase(subj).SR_dec = SR_dec;
end

disp('Calculated spike ratios')
% --> too much variation in distribution of spike ratios! Therefore,
% combine all spike ratios and then determine thresholds!



%% analysis ERs_PS

% pre-allocation
ERs_PS_pat = NaN(size(dataBase,2),4);

for subj = 1:size(dataBase,2)
    
    % chi squared test
    [tbl,~,pchi] = crosstab(dataBase(subj).BBmat(:),dataBase(subj).ERmat(:));
    
    % odds ratio
    OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
    
    ERs_PS_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    
    dataBase(subj).OR_ERBB = OR;
    dataBase(subj).p_ERBB = pchi;
end

%% figure stacked bar plot

figure(1), bar(ERs_PS_pat)
legend({'ER & PS','ER & no PS','No ER & PS','No ER & no PS'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis ERs_spikes && BB_spikes --> only in part of the data

% only part of ERs is used
for subj = 1:size(dataBase,2)
    dataBase(subj).ERmatIED = dataBase(subj).ERmat(dataBase(subj).IEDch,:);
    dataBase(subj).BBmatIED = dataBase(subj).BBmat(dataBase(subj).IEDch,:);
    
end

%% analysis ERs_spikes
% pre-allocation
SRchange_ER_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = [1:4,6,8,10]%1:size(dataBase,2)
    
    % chi squared test
    [tbl,~,pchi] = crosstab(dataBase(subj).ERmatIED(:),dataBase(subj).SR_change(:));
    
    % odds ratio
    OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
    
    SRchange_ER_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    
    dataBase(subj).OR_ER_SRchange = OR;
    dataBase(subj).p_ER_SRchange = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRchange_ER_pat)
legend({'ER & spike change','ER & no spike change','No ER & spike change','No ER & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis ERs_spikes increase
% pre-allocation
SRinc_ER_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)
    
    if dataBase(subj).p_ER_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).ERmatIED(:),dataBase(subj).SR_inc(:));
        
        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
        
        SRinc_ER_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRinc_ER_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end
    
    dataBase(subj).OR_ER_SRinc = OR;
    dataBase(subj).p_ER_SRinc = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRinc_ER_pat)
legend({'ER & spike increase','ER & no spike change','No ER & spike increase','No ER & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis ERs_spikes decrease
% pre-allocation
SRdec_ER_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)
    
    if dataBase(subj).p_ER_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).ERmatIED(:),dataBase(subj).SR_dec(:));
        
        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
        
        SRdec_ER_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRdec_ER_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end
    
    dataBase(subj).OR_ER_SRdec = OR;
    dataBase(subj).p_ER_SRdec = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRdec_ER_pat)
legend({'ER & spike decrease','ER & no spike change','No ER & spike decrease','No ER & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])


%% analysis PSs_spikes
% pre-allocation
SRchange_BB_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = [1:4,6,8,10]%1:size(dataBase,2)
    
    % chi squared test
    [tbl,~,pchi] = crosstab(dataBase(subj).BBmatIED(:),dataBase(subj).SR_change(:));
    
    % odds ratio
    OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
    
    SRchange_BB_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    
    dataBase(subj).OR_BB_SRchange = OR;
    dataBase(subj).p_BB_SRchange = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRchange_BB_pat)
legend({'SP & spike change','SP & no spike change','No SP & spike change','No SP & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis BBs_spikes increase
% pre-allocation
SRinc_BB_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)
    
    if dataBase(subj).p_BB_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).BBmatIED(:),dataBase(subj).SR_inc(:));
        
        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
        
        SRinc_BB_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRinc_BB_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end
    
    dataBase(subj).OR_BB_SRinc = OR;
    dataBase(subj).p_BB_SRinc = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRinc_BB_pat)
legend({'SP & spike increase','SP & no spike change','No SP & spike increase','No SP & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis SP_spikes decrease
% pre-allocation
SRdec_BB_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)
    
    if dataBase(subj).p_BB_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).BBmatIED(:),dataBase(subj).SR_dec(:));
        
        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));
        
        SRdec_BB_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRdec_BB_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end
    
    dataBase(subj).OR_BB_SRdec = OR;
    dataBase(subj).p_BB_SRdec = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRdec_BB_pat)
legend({'SP & spike decrease','SP & no spike change','No SP & spike decrease','No SP & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])


%% figure

figure, 

subj=1

for i=1%:size(dataBase(subj).spikes,2)
    plot(dataBase(subj).spikes{i},i*ones(1,size(dataBase(subj).spikes{i},2)),'.')

end

