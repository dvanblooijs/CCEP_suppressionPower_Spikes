% CCEP_powerSuppression_cceps_spikes
% this code analyzes the relationship between cceps, spikes and
% powerSuppression.
% author: D van Blooijs
% date: April 2019

clc
clear
myDataPath = setLocalDataPath(1);

%% settings

files = dir(myDataPath.dataPath);
idx_subj = contains({files(:).name},'sub-');
files_subj = files(idx_subj);
cfg = struct([]);

for subj = 1:size(files_subj,1)

    cfg(subj).sub_labels = files_subj(subj).name;

    files = dir(fullfile(files_subj(subj).folder,files_subj(subj).name));
    idx_ses = contains({files(:).name},'ses-');
    files_ses = files(idx_ses);

    cfg(subj).ses_label = files_ses(1).name;

    cfg(subj).task_label = 'task-SPESclin';

    files = dir(fullfile(files_ses(1).folder,files_ses(1).name,'ieeg'));
    idx_eeg = contains({files(:).name},'.eeg');
    files_eeg = files(idx_eeg);

    for run = 1:size(files_eeg,1)
        runTemp = extractBetween(files_eeg(run).name,'run-','_ieeg');
        cfg(subj).run_label{run} = ['run-', runTemp{1}];
    end
end

%% load ECoGs with SPES from X patients

dataBase = load_ECoGdata(myDataPath,cfg);

%% preprocessing CCEP in ECoG
cfg = [];

% sort stimulation pairs
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% load CCEPs

files = (dir(myDataPath.CCEPpath));

for subj = 1:size(dataBase,2)

    idx_sub = contains({files(:).name}, dataBase(subj).sub_label) & [files(:).isdir] == 1;

    ses = dir(fullfile(files(idx_sub).folder,files(idx_sub).name));
    idx_ses = contains({ses(:).name}, dataBase(subj).ses_label);

    file = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name));
    idx_ccepfile = contains({file(:).name},'N1s');

    if sum(idx_ccepfile)>0
        dataBase(subj).ccep = load(fullfile(file(idx_ccepfile).folder,file(idx_ccepfile).name));
    end

end

disp('CCEPs are loaded')

clear idx_ERSPfile idx_ccepfile file files idx_ses idx_sub ses subj ans

%% load ERSP

files = (dir(myDataPath.ERSPoutput));

for subj = 1:size(dataBase,2)

    idx_sub = contains({files(:).name}, dataBase(subj).sub_label) & [files(:).isdir] == 1;

    ses = dir(fullfile(files(idx_sub).folder,files(idx_sub).name));
    idx_ses = contains({ses(:).name}, dataBase(subj).ses_label);

    run = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name));
    idx_run = contains({run(:).name},dataBase(subj).run_label);

    file = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name,run(idx_run).name));
    idx_ERSPfile = contains({file(:).name},'ERSP.mat');

    if sum(idx_ERSPfile)>0
        dataBase(subj).ERSP = load(fullfile(file(idx_ERSPfile).folder,file(idx_ERSPfile).name));
    end
end

disp('ERSPs are loaded')

clear idx_ERSPfile idx_ccepfile file files idx_ses idx_sub ses subj ans run idx_run

%% load spikes

files = dir(myDataPath.visIEDpath);

for subj = 1:size(dataBase,2)

    idx_sub = contains({files(:).name}, dataBase(subj).sub_label);

    if any(idx_sub)
        ses = dir(fullfile(files(idx_sub).folder,files(idx_sub).name));
        idx_ses = contains({ses(:).name}, dataBase(subj).ses_label);

        file = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name));
        idx_detIED = contains({file(:).name},'detIEDs');

        if sum(idx_detIED)>0
            load(fullfile(file(idx_detIED).folder,file(idx_detIED).name));
            dataBase(subj).IEDs = spikespat;
        end
    end
end

disp('Spikes are loaded')

clear subj ses spikespat idx_detIED idx_ses idx_sub file  files

%% make both the ccep and ersp matrix equal regarding stimulation pairs and response channels

for subj = 1:size(dataBase,2)

    % duplicate original cc_stimchans, because you will adapt those later
    % in this section of the script
    dataBase(subj).cc_stimchans_orig = dataBase(subj).cc_stimchans;
    dataBase(subj).cc_stimsets_orig = dataBase(subj).cc_stimsets;
    dataBase(subj).tt_epoch_sorted_orig = dataBase(subj).tt_epoch_sorted;

    CCEPmat = dataBase(subj).ccep.checked;
    ERSPmat = dataBase(subj).ERSP.checked;
    CCEPmat_comb = [];
    ERSPmat_comb = [];
    exclude_ccep = [];
    exclude_ersp = [];

    if size(CCEPmat,1) > size(CCEPmat,2) % [channels x stim pairs]
        CCEPmat = CCEPmat';
    end

    if size(ERSPmat,1) > size(ERSPmat,2) % [channels x stim pairs]
        ERSPmat = ERSPmat';
    end

    % only continue if the order of channels is equal in ccep and ersp
    % preprocessing
    if isequal(dataBase(subj).ccep.ch,dataBase(subj).ERSP.ch)

        % find total number of stimulation pairs and whether this one is
        % larger in ccep or in ersp preprocessing
        [stim_size,I] = max([size(CCEPmat,1), size(ERSPmat,1)]);

        % for each stimulation pair
        for stimp = 1:stim_size

            % when more stims in CCEP than in ERSP
            if I == 1

                % find where stimulus of CCEP can be found in ERSP
                idx = find(strcmp(dataBase(subj).ERSP.cc_stimchans(:,1),dataBase(subj).ccep.cc_stimchans(stimp,1)) & ...
                    strcmp(dataBase(subj).ERSP.cc_stimchans(:,2),dataBase(subj).ccep.cc_stimchans(stimp,2)), 1);

                % if this stimulus is found in both ERSP and CCEP data
                if ~isempty(idx)

                    CCEPmat_comb(stimp,:) = CCEPmat(stimp,:); %#ok<SAGROW>
                    ERSPmat_comb(stimp,:) = ERSPmat(idx,:);%#ok<SAGROW>

                    % the stimulated electrodes must be NaN
                    CCEPmat_comb(stimp, dataBase(subj).ccep.cc_stimsets(stimp,:)) = NaN; %#ok<SAGROW>
                    ERSPmat_comb(stimp, dataBase(subj).ERSP.cc_stimsets(idx,:)) = NaN; %#ok<SAGROW>

                else
                    % exclude this row in the ccep-matrix
                    exclude_ccep = [exclude_ccep, stimp];  %#ok<AGROW>

                end

            elseif I == 2 % when more stims in ERSP than in CCEP

                idx = find(strcmp(dataBase(subj).ccep.cc_stimchans(:,1),dataBase(subj).ERSP.cc_stimchans(stimp,1)) & ...
                    strcmp(dataBase(subj).ccep.cc_stimchans(:,2),dataBase(subj).ERSP.cc_stimchans(stimp,2)), 1);

                % if this stimulus is found in both ERSP and CCEP data
                if ~isempty(idx)

                    CCEPmat_comb(stimp,:) = CCEPmat(idx,:); %#ok<SAGROW>
                    ERSPmat_comb(stimp,:) = ERSPmat(stimp,:); %#ok<SAGROW>

                    % the stimulated electrodes must be NaN
                    CCEPmat_comb(stimp, dataBase(subj).ccep.cc_stimsets(idx,:)) = NaN; %#ok<SAGROW>
                    ERSPmat_comb(stimp, dataBase(subj).ERSP.cc_stimsets(stimp,:)) = NaN; %#ok<SAGROW>

                else
                    % exclude this row in the ersp-matrix
                    exclude_ersp = [exclude_ersp, stimp]; %#ok<AGROW>

                end
            end
        end

        % if some rows need to be excluded
        if ~isempty(exclude_ersp)

            % if last row needs to be excluded, this one is not in the
            % matrix at all
            if exclude_ersp(end) > size(ERSPmat_comb,1)
                ERSPmat_comb(exclude_ersp(1:end-1),:) = []; %#ok<SAGROW>

            else
                ERSPmat_comb(exclude_ersp,:) = []; %#ok<SAGROW>

            end
        end

        % if some rows need to be excluded
        if ~isempty(exclude_ccep)

            % if last row needs to be excluded, this one is not in the
            % matrix at all
            if exclude_ccep(end) > size(CCEPmat_comb,1)

                CCEPmat_comb(exclude_ccep(1:end-1),:) = []; %#ok<SAGROW>
            else
                CCEPmat_comb(exclude_ccep,:) = []; %#ok<SAGROW>

            end
        end

        if size(ERSPmat_comb,1) == size(CCEPmat_comb,1) && size(ERSPmat_comb,2) == size(CCEPmat_comb,2)

        else
            error('Sizes of ERSP matrix and ccep matrix in subject %s are inequal!',dataBase(subj).sub_label)
        end

        cc_stimchans = dataBase(subj).cc_stimchans_orig;
        cc_stimsets = dataBase(subj).cc_stimsets_orig;
        tt_epoch_sorted = dataBase(subj).tt_epoch_sorted_orig;
        cc_stimchans([exclude_ersp, exclude_ccep],:) = [];
        cc_stimsets([exclude_ersp, exclude_ccep],:) = [];
        tt_epoch_sorted(:,[exclude_ersp, exclude_ccep],:) = [];
        dataBase(subj).cc_stimchans = cc_stimchans;
        dataBase(subj).cc_stimsets = cc_stimsets;
        dataBase(subj).tt_epoch_sorted = tt_epoch_sorted;
        dataBase(subj).CCEPmat = CCEPmat_comb;
        dataBase(subj).ERSPmat = ERSPmat_comb;

    else
        error('Channels are not equal in ERSP and CCEP analysis')

    end
end

clear ERSPmat CCEPmat ERSPmat_comb CCEPmat_comb stimp exclude_ccep exclude_ersp I idx stim_size subj cc_stimchans cc_stimsets tt_epoch_sorted

%% find spikes in each stimulus pair

selDur = 1.1; % seconds pre and post stimulation to find spikes
t = -cfg.epoch_prestim+1/2048:1/2048 : (cfg.epoch_length - cfg.epoch_prestim);

for subj = 1:size(dataBase,2)

    if ~isempty(dataBase(subj).IEDs)

        % pre-allocation
        totEpochspikessamp = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2),10); % [stimpairs x IED channels x trials]
        selEpochspikessamp = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2),10); % [stimpairs x IED channels x trials]
        spikespre = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        spikespost = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        spikesratio = NaN(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        p_spikes = NaN(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]

        for stim = 1:size(dataBase(subj).tt_epoch_sorted,2) % for each stim pair
            for ch = 1:size(dataBase(subj).IEDs.IEDch,1)

                if strcmpi(dataBase(subj).cc_stimchans{stim,1},dataBase(subj).IEDs.IEDchan{ch}) || ...
                        strcmpi(dataBase(subj).cc_stimchans{stim,2},dataBase(subj).IEDs.IEDchan{ch})

                    % when channel with IEDs is stimulated, the spike ratio
                    % should be NaN
                    spikespre{stim,ch}(n) = NaN;
                    spikespost{stim,ch}(n) = NaN;
                    spikesratio(stim,ch) = NaN;
                    p_spikes(stim,ch) = NaN;

                else
                    for n = 1:size(dataBase(subj).cc_epoch_sorted,2) % for each trial in each stim pair (so usually 10)

                        % start and end of epoch for each stimulus pair
                        total_epoch = squeeze(dataBase(subj).tt_epoch_sorted(n,stim,:));

                        % 200 ms around stimulus artefact should be not included in
                        % spike counts
                        samp_stim = round((cfg.epoch_prestim-0.1) * 2048) : round((cfg.epoch_prestim+0.1) * 2048);
                        total_epoch(samp_stim) = NaN; % which is the complete duration of the period that is used for analysis (usually 4s: 2s pre-stim, 2s post-stim)

                        % selected epoch
                        samp_sel = round((cfg.epoch_prestim-selDur) * 2048) : round((cfg.epoch_prestim+selDur) * 2048);
                        selEpoch = total_epoch;
                        selEpoch(setdiff(1:size(total_epoch,1),samp_sel)) = NaN;

                        % find all detected spikes in this epoch of
                        % cfg.epoch_length
                        spikessamp = intersect(total_epoch,dataBase(subj).IEDs.spikesdet{ch});

                        % find sample number and calculate spikes prior and post stim
                        epochsamp = NaN(size(spikessamp)); removeCel = [];
                        preCount = 0; postCount = 0;
                        for num = 1:size(spikessamp,1)
                            epochsamp(num) = find(total_epoch == spikessamp(num));

                            if t(epochsamp(num)) < 0 && t(epochsamp(num)) > -selDur
                                preCount = preCount+1;
                            elseif t(epochsamp(num)) > 0  && t(epochsamp(num)) < selDur
                                postCount = postCount+1;
                            else
                                removeCel = [removeCel, num]; %#ok<AGROW>
                            end
                        end

                        selEpochsamp = epochsamp;
                        selEpochsamp(removeCel) = [];

                        %                     figure(1),
                        %                     plot(total_epoch,ones(size(total_epoch))),
                        %                     hold on,
                        %                     plot(selEpoch,1.5*ones(size(selEpoch))),
                        %                     plot(dataBase(subj).IEDs.spikesdet{ch},2*ones(1,size(dataBase(subj).IEDs.spikesdet{ch},2)),'*'),
                        %                     hold off,
                        %                     xlim([min(total_epoch),max(total_epoch)]),
                        %                     ylim([0 3])
                        %                     title(sprintf('%s-%s, %s, trial %d',...
                        %                         dataBase(subj).cc_stimchans{stim,:},...
                        %                         dataBase(subj).ch{dataBase(subj).IEDs.IEDch(ch)}, ...
                        %                         n))

                        %pause

                        spikespre{stim,ch}(n) = preCount;
                        spikespost{stim,ch}(n) = postCount;

                        totEpochspikessamp{stim,ch,n} = epochsamp;
                        selEpochspikessamp{stim,ch,n} = selEpochsamp;

                    end


                    spikesratio(stim,ch) = sum(spikespost{stim,ch})/sum(spikespre{stim,ch});
                    p = signtest(spikespre{stim,ch},spikespost{stim,ch});
                    p_spikes(stim,ch) = p;
                end
            end
        end

        dataBase(subj).IEDs.totEpochspikessamp  = totEpochspikessamp;
        dataBase(subj).IEDs.selEpochspikessamp  = selEpochspikessamp;
        dataBase(subj).IEDs.spikespre   = spikespre;
        dataBase(subj).IEDs.spikespost  = spikespost;
        dataBase(subj).IEDs.spikesratio = spikesratio;
        dataBase(subj).IEDs.spikes_p    = p_spikes;

    else
        dataBase(subj).IEDtest = 0;
    end
end

clear subj ch i samp_stim spikessamp startepoch stim stopepoch time_stim total_epoch

%% make figure with spikes per channel in all stimulus pairs
% close all
%
% for subj = 1%:size(dataBase,2)
%
%     if ~isempty(dataBase(subj).IEDs)
%         spikessamp = dataBase(subj).IEDs.spikessamp;
%         spikes_p = dataBase(subj).IEDs.spikes_p;
%
%         %         ccep_stimchans = vertcat(dataBase(subj).ccep.cceps(:).stimchans);
%
%         time_stim = cfg.epoch_prestim :5:50;
%         time_stimpre = time_stim - dur;
%         time_stimpost = time_stim + dur;
%
%         stimlabels = cell(size(dataBase(subj).cc_stimchans,1),1);
%         for i=1:size(dataBase(subj).cc_stimchans,1)
%             stimlabels{i} = sprintf('%s-%s',dataBase(subj).cc_stimchans{i,1},dataBase(subj).cc_stimchans{i,2});
%         end
%
%         ymin = 0;
%         ymax = size(spikessamp,1)+1;
%
%         for ch = 1:size(dataBase(subj).IEDs.IEDch,2)
%             h=figure;
%
%             % plot moment of stimulation
%             plot(time_stim.*ones(ymax+1,1),ymin:ymax,'k')
%
%             hold on
%             plot(time_stimpre.*ones(ymax+1,1),ymin:ymax,':k')
%             plot(time_stimpost.*ones(ymax+1,1),ymin:ymax,':k')
%
%             for stim = 1:size(spikessamp,1)
%                 % if channel is not part of stimulus pair
%                 if ~strcmpi(dataBase(subj).cc_stimchans{stim,1},dataBase(subj).IEDs.IEDchan{ch}) &&...
%                         ~strcmpi(dataBase(subj).cc_stimchans{stim,2},dataBase(subj).IEDs.IEDchan{ch})
%
%
%                     %                     ccepstimnum = dataBase(subj).cc_stimsets(stim,1) == ccep_stimchans(:,1) &...
%                     %                         dataBase(subj).cc_stimsets(stim,2) == ccep_stimchans(:,2);
%                     %
%                     % if there is a connection
%                     %                     if ismember(dataBase(subj).IEDs.IEDch(ch),dataBase(subj).ccep.cceps(ccepstimnum).ccepsvis)
%                     %                         plot(spikessamp{stim,ch}/2048,stim*ones(size(spikessamp{stim,ch})),'.','Color',[0.2 0.2 0.2]) % darker gray
%                     %                     else
%                     %                         plot(spikessamp{stim,ch}/2048,stim*ones(size(spikessamp{stim,ch})),'.','Color',[0.8 0.8 0.8]) % lighter gray
%                     %                     end
%                     %                     % if there is a significant modulation
%                     if spikes_p(stim,ch) <0.05
%                         plot(spikessamp{stim,ch}/2048,stim*ones(size(spikessamp{stim,ch})),'.','Color',[0.2 0.2 0.2]) % darker gray
%                     else
%                         plot(spikessamp{stim,ch}/2048,stim*ones(size(spikessamp{stim,ch})),'.','Color',[0.8 0.8 0.8]) % lighter gray
%                     end
%                 end
%             end
%
%             hold off
%             h.Units = 'normalized';
%             h.Position = [0.14 0.07 0.6 0.8];
%             ax = gca;
%             ax.YTick = ymin:ymax;
%             ax.YTickLabel = [{' '},stimlabels(:)',{' '}];
%             ylim([ymin ymax])
%             xlim([0 50])
%             title(sprintf('Detected spikes per stimulus pair in %s, %s, channel %s',dataBase(subj).sub_label,dataBase(subj).run_label,dataBase(subj).IEDs.IEDchan{ch}))
%         end
%     end
% end
%
% clear ax ccep_stimchans ccepstimnum ch h i spikessamp stim stimlabels subj time_stim ymax ymin






%% OLD
%% spikes vs cceps, based on whether p<0.05 was found with sign-test pre- and post stimulation

for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)
        spikes = zeros(size(dataBase(subj).IEDs.spikes_p));

        spikes(dataBase(subj).IEDs.spikes_p<0.05) = 1;
        CCEPmat = dataBase(subj).CCEPmat(:,dataBase(subj).IEDs.IEDch);

        [tbl] = crosstab(CCEPmat(:),spikes(:));

        if size(tbl,1) == 2 && size(tbl,2) == 2

            [~,p,stats] = fishertest(tbl);
        else
            p = NaN;
            stats.OddsRatio = NaN;
            stats.ConfidenceInterval = [-Inf Inf];
        end

        dataBase(subj).stat_ccep_spikes.OR = stats.OddsRatio;
        dataBase(subj).stat_ccep_spikes.CI = stats.ConfidenceInterval;
        dataBase(subj).stat_ccep_spikes.tbl = tbl;
        dataBase(subj).stat_ccep_spikes.p = p;

        fprintf('--- %s: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
            dataBase(subj).sub_label, stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)
    end
end

%% make forest plot cceps vs spikes, based on whether p<0.05 was found with sign-test pre- and post stimulation
close all

maxCI = NaN(size(dataBase,2),1);
h=figure(1);
hold on
for subj =1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)

        xline = dataBase(subj).stat_ccep_spikes.CI(1):0.1:dataBase(subj).stat_ccep_spikes.CI(2);

        maxCI(subj) = dataBase(subj).stat_ccep_spikes.CI(2);
        plot(xline,subj*ones(size(xline,2),1),'k')

        plot(dataBase(subj).stat_ccep_spikes.OR,subj,'o','MarkerFaceColor','r','MarkerEdgeColor','r')

        if dataBase(subj).stat_ccep_spikes.p <0.001
            p_str = '***';
        elseif dataBase(subj).stat_ccep_spikes.p <0.01 && dataBase(subj).stat_ccep_spikes.p > 0.001
            p_str = '**';
        elseif dataBase(subj).stat_ccep_spikes.p <0.05 && dataBase(subj).stat_ccep_spikes.p > 0.01
            p_str = '*';
        else
            p_str = ' ';
        end

        text(dataBase(subj).stat_ccep_spikes.CI(2)+0.5,subj,p_str)%,'FontSize',8)
    end
end
plot(ones(size(dataBase,2)+2,1),0:size(dataBase,2)+1,'k:')
hold off

ylim([0 size(dataBase,2)+1])
xlim([0 ceil(max(maxCI))+2])

h.Units = 'normalized';
h.Position = [0.3 0.4 0.5 0.5];

ax = gca;
ax.YTickLabel = [{' '},{dataBase(:).sub_label},{' '}];
xlabel('Odds ratio')
title('Occurrence of spikes modulation and CCEP')

%% spikes vs power suppression, based on whether p<0.05 was found with sign-test pre- and post stimulation

for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)
        spikes = zeros(size(dataBase(subj).IEDs.spikes_p));

        spikes(dataBase(subj).IEDs.spikes_p<0.05) = 1;
        ERSPmat = dataBase(subj).ERSPmat(:,dataBase(subj).IEDs.IEDch);

        tbl = crosstab(ERSPmat(:),spikes(:));

        if size(tbl,1) == 2 && size(tbl,2) == 2
            [~,p,stats] = fishertest(tbl);
        else
            p = NaN;
            stats.OddsRatio = NaN;
            stats.ConfidenceInterval = [-Inf Inf];
        end

        dataBase(subj).stat_ERSP_spikes.OR = stats.OddsRatio;
        dataBase(subj).stat_ERSP_spikes.CI = stats.ConfidenceInterval;
        dataBase(subj).stat_ERSP_spikes.tbl = tbl;
        dataBase(subj).stat_ERSP_spikes.p = p;

        fprintf('--- %s: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
            dataBase(subj).sub_label, stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)
    end
end

%% make forest plot power suppression vs spikes, based on whether p<0.05 was found with sign-test pre- and post stimulation
close all

maxCI = NaN(size(dataBase,2),1);
h=figure(1);
hold on
for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)

        xline = dataBase(subj).stat_ERSP_spikes.CI(1):0.1:dataBase(subj).stat_ERSP_spikes.CI(2);

        maxCI(subj) = dataBase(subj).stat_ERSP_spikes.CI(2);
        plot(xline,subj*ones(size(xline,2),1),'k')

        plot(dataBase(subj).stat_ERSP_spikes.OR,subj,'o','MarkerFaceColor','r','MarkerEdgeColor','r')

        if dataBase(subj).stat_ERSP_spikes.p <0.001
            p_str = '***';
        elseif dataBase(subj).stat_ERSP_spikes.p <0.01 && dataBase(subj).stat_ERSP_spikes.p > 0.001
            p_str = '**';
        elseif dataBase(subj).stat_ERSP_spikes.p <0.05 && dataBase(subj).stat_ERSP_spikes.p > 0.01
            p_str = '*';
        else
            p_str = ' ';
        end

        text(dataBase(subj).stat_ERSP_spikes.CI(2)+0.5,subj,p_str)%,'FontSize',8)
    end
end
plot(ones(size(dataBase,2)+2,1),0:size(dataBase,2)+1,'k:')
hold off

ylim([0 size(dataBase,2)+1])
xlim([0 ceil(max(maxCI))+2])

h.Units = 'normalized';
h.Position = [0.3 0.4 0.5 0.5];

ax = gca;
ax.YTickLabel = [{' '},{dataBase(:).sub_label},{' '}];
xlabel('Odds ratio')
title('Occurrence of spikes modulation and power suppression')


%% make scatter plot latency en spikes ratio
clc
for subj = [1:5,7:size(dataBase,2)]
    if ~isempty(dataBase(subj).IEDs)

        CCEPmat = dataBase(subj).ccep.n1_peak_sample;
        CCEPmat(dataBase(subj).ccep.checked == 0) = NaN;
        if size(CCEPmat,1) > size(CCEPmat,2)
            CCEPmat = CCEPmat';
        end

        CCEPmat = (CCEPmat(:,dataBase(subj).IEDs.IEDch)- 2*2048)/2048*1000;
        spikesratio = dataBase(subj).IEDs.spikesratio ;

        spikesratio(isnan(CCEPmat) | isnan(spikesratio)) = NaN;
        CCEPmat(isnan(CCEPmat) | isnan(spikesratio)) = NaN;

        [rho,p] = corr(CCEPmat(:),abs(log(spikesratio(:))),'rows','pairwise','type','spearman');

        fprintf('--- %s: rho = %1.3f, p = %1.10f --- \n',...
            dataBase(subj).sub_label, rho, p)

        figure(subj),
        scatter(CCEPmat(:),abs(log(spikesratio(:))))
    end
end

%% spikes vs power suppression, based on whether p<0.05 was found with sign-test pre- and post stimulation

for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)
        spikes = zeros(size(dataBase(subj).IEDs.spikes_p));

        spikes(dataBase(subj).IEDs.spikes_p<0.05) = 1;
        ERSPmat = dataBase(subj).ERSPmat(:,dataBase(subj).IEDs.IEDch);

        tbl = crosstab(ERSPmat(:),spikes(:));

        [~,p,stats] = fishertest(tbl);

        dataBase(subj).stat_ERSP_spikes.OR = stats.OddsRatio;
        dataBase(subj).stat_ERSP_spikes.CI = stats.ConfidenceInterval;
        dataBase(subj).stat_ERSP_spikes.tbl = tbl;
        dataBase(subj).stat_ERSP_spikes.p = p;

        fprintf('--- %s: OR = %2.1f, CI = %2.1f-%2.1f, p = %1.10f --- \n',...
            dataBase(subj).sub_label, stats.OddsRatio, stats.ConfidenceInterval(1), stats.ConfidenceInterval(2) , p)
    end
end

%% make forest plot ERSPs vs spikes based on sign test
close all

maxCI = NaN(size(dataBase,2),1);
h=figure(1);
hold on
for subj =1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDs)

        xline = dataBase(subj).stat_ERSP_spikes.CI(1):0.1:dataBase(subj).stat_ERSP_spikes.CI(2);

        maxCI(subj) = dataBase(subj).stat_ERSP_spikes.CI(2);
        plot(xline,subj*ones(size(xline,2),1),'k')

        plot(dataBase(subj).stat_ERSP_spikes.OR,subj,'o','MarkerFaceColor','r','MarkerEdgeColor','r')

        if dataBase(subj).stat_ERSP_spikes.p <0.001
            p_str = '***';
        elseif dataBase(subj).stat_ERSP_spikes.p <0.01 && dataBase(subj).stat_ERSP_spikes.p > 0.001
            p_str = '**';
        elseif dataBase(subj).stat_ERSP_spikes.p <0.05 && dataBase(subj).stat_ERSP_spikes.p > 0.01
            p_str = '*';
        else
            p_str = ' ';
        end

        text(dataBase(subj).stat_ERSP_spikes.CI(2)+0.5,subj,p_str)%,'FontSize',8)
    end
end
plot(ones(size(dataBase,2)+2,1),0:size(dataBase,2)+1,'k:')
hold off

ylim([0 size(dataBase,2)+1])
xlim([0 ceil(max(maxCI))+2])

h.Units = 'normalized';
h.Position = [0.3 0.4 0.5 0.5];

ax = gca;
ax.YTickLabel = [{' '},{dataBase(:).sub_label},{' '}];
xlabel('Odds ratio')
title('Occurrence of spikes modulation and power suppression')


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
    spikeratio = NaN(size(dataBase(subj).IEDs.IEDch,2),size(dataBase(subj).tt_epoch_sorted,2));

    for IEDchan = 1:size(dataBase(subj).IEDs.IEDch,2)
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



%% analysis cceps_PS

% pre-allocation
cceps_PS_pat = NaN(size(dataBase,2),4);

for subj = 1:size(dataBase,2)

    % chi squared test
    [tbl,~,pchi] = crosstab(dataBase(subj).ERSPmat(:),dataBase(subj).ccepmat(:));

    % odds ratio
    OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

    cceps_PS_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));

    dataBase(subj).OR_ccepERSP = OR;
    dataBase(subj).p_ccepERSP = pchi;
end

%% figure stacked bar plot

figure(1), bar(cceps_PS_pat)
legend({'ccep & PS','ccep & no PS','No ccep & PS','No ccep & no PS'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis cceps_spikes && ERSP_spikes --> only in part of the data

% only part of cceps is used
for subj = 1:size(dataBase,2)
    dataBase(subj).ccepmatIED = dataBase(subj).ccepmat(dataBase(subj).IEDch,:);
    dataBase(subj).ERSPmatIED = dataBase(subj).ERSPmat(dataBase(subj).IEDch,:);

end

%% analysis cceps_spikes
% pre-allocation
SRchange_ccep_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = [1:4,6,8,10]%1:size(dataBase,2)

    % chi squared test
    [tbl,~,pchi] = crosstab(dataBase(subj).ccepmatIED(:),dataBase(subj).SR_change(:));

    % odds ratio
    OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

    SRchange_ccep_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));

    dataBase(subj).OR_ccep_SRchange = OR;
    dataBase(subj).p_ccep_SRchange = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRchange_ccep_pat)
legend({'ccep & spike change','ccep & no spike change','No ccep & spike change','No ccep & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis cceps_spikes increase
% pre-allocation
SRinc_ccep_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)

    if dataBase(subj).p_ccep_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).ccepmatIED(:),dataBase(subj).SR_inc(:));

        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

        SRinc_ccep_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRinc_ccep_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end

    dataBase(subj).OR_ccep_SRinc = OR;
    dataBase(subj).p_ccep_SRinc = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRinc_ccep_pat)
legend({'ccep & spike increase','ccep & no spike change','No ccep & spike increase','No ccep & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis cceps_spikes decrease
% pre-allocation
SRdec_ccep_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)

    if dataBase(subj).p_ccep_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).ccepmatIED(:),dataBase(subj).SR_dec(:));

        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

        SRdec_ccep_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRdec_ccep_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end

    dataBase(subj).OR_ccep_SRdec = OR;
    dataBase(subj).p_ccep_SRdec = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRdec_ccep_pat)
legend({'ccep & spike decrease','ccep & no spike change','No ccep & spike decrease','No ccep & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])


%% analysis PSs_spikes
% pre-allocation
SRchange_ERSP_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = [1:4,6,8,10]%1:size(dataBase,2)

    % chi squared test
    [tbl,~,pchi] = crosstab(dataBase(subj).ERSPmatIED(:),dataBase(subj).SR_change(:));

    % odds ratio
    OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

    SRchange_ERSP_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));

    dataBase(subj).OR_ERSP_SRchange = OR;
    dataBase(subj).p_ERSP_SRchange = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRchange_ERSP_pat)
legend({'SP & spike change','SP & no spike change','No SP & spike change','No SP & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis ERSPs_spikes increase
% pre-allocation
SRinc_ERSP_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)

    if dataBase(subj).p_ERSP_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).ERSPmatIED(:),dataBase(subj).SR_inc(:));

        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

        SRinc_ERSP_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRinc_ERSP_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end

    dataBase(subj).OR_ERSP_SRinc = OR;
    dataBase(subj).p_ERSP_SRinc = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRinc_ERSP_pat)
legend({'SP & spike increase','SP & no spike change','No SP & spike increase','No SP & no spike change'})
xlabel('Patient #')
ylabel('Percentage of total number of connections')
ylim([0 1])

%% analysis SP_spikes decrease
% pre-allocation
SRdec_ERSP_pat = NaN(size(dataBase,2),4);

% chi squared test
for subj = 1:size(dataBase,2)

    if dataBase(subj).p_ERSP_SRchange <= 0.05
        % chi squared test
        [tbl,~,pchi] = crosstab(dataBase(subj).ERSPmatIED(:),dataBase(subj).SR_dec(:));

        % odds ratio
        OR = (tbl(4)/tbl(3))/(tbl(2)/tbl(1));

        SRdec_ERSP_pat(subj,:) = [tbl(4), tbl(3), tbl(2), tbl(1)]/sum(tbl(:));
    else
        OR = [];
        SRdec_ERSP_pat(subj,:) = [NaN, NaN, NaN, NaN];
        pchi = [];
    end

    dataBase(subj).OR_ERSP_SRdec = OR;
    dataBase(subj).p_ERSP_SRdec = pchi;
end

%% figure stacked bar plot

figure(1), bar(SRdec_ERSP_pat)
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

