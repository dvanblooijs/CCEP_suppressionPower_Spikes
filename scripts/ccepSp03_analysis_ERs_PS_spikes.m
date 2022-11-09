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

% housekeeping
clear files files_eeg files_ses files_subj idx_eeg idx_ses idx_subj run runTemp subj

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

% housekeeping
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

% housekeeping
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

%% find spikes in a specific time window for each response electrode 
% after stimulating each stimulus pair
% and calculate the spikesratio

% define period in which spikes are found for spikesratio
selDur = 1; % seconds pre and post stimulation to find spikes
fs = dataBase(1).ccep_header.Fs;
t = -cfg.epoch_prestim+1/fs:1/fs : (cfg.epoch_length - cfg.epoch_prestim);

for subj = 1:size(dataBase,2)

    if ~isempty(dataBase(subj).IEDs)

        % pre-allocation
        totEpochspikessamp = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2),10); % [stimpairs x IED channels x trials]
        selEpochspikessamp = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2),10); % [stimpairs x IED channels x trials]
        spikespre = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        spikespost = cell(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        spikesratio = NaN(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]
%         p_spikes = NaN(size(dataBase(subj).cc_stimchans,1),size(dataBase(subj).IEDs.IEDch,2)); % [stimpairs x IED channels]

        for stim = 1:size(dataBase(subj).tt_epoch_sorted,2) % for each stim pair
            for ch = 1:size(dataBase(subj).IEDs.IEDch,1)
            
                % if channel is not stimulated in this specific stim trial
                if strcmpi(dataBase(subj).cc_stimchans{stim,1},dataBase(subj).IEDs.IEDchan{ch}) || ...
                        strcmpi(dataBase(subj).cc_stimchans{stim,2},dataBase(subj).IEDs.IEDchan{ch})

                    % when channel with IEDs is stimulated, the spike ratio
                    % should be NaN
                    spikespre{stim,ch}(1:size(dataBase(subj).cc_epoch_sorted,2)) = NaN;
                    spikespost{stim,ch}(1:size(dataBase(subj).cc_epoch_sorted,2)) = NaN;
                    spikesratio(stim,ch) = NaN;
%                     p_spikes(stim,ch) = NaN;

                else
                    for n = 1:size(dataBase(subj).cc_epoch_sorted,2) % for each trial in each stim pair (so usually 10)

                        % time window in samples of epoch (usually 4s: 2s pre-stim, 2s post-stim)
                        total_epoch = squeeze(dataBase(subj).tt_epoch_sorted(n,stim,:));

                        % 200 ms around stimulus artefact should be not included in
                        % spike counts due to stimulus artefact (100 ms pre
                        % and 100 ms post stimulus)
                        samp_stim = round((cfg.epoch_prestim-0.1) * fs) : round((cfg.epoch_prestim+0.1) * fs);
                        total_epoch(samp_stim) = NaN; 

                        % find all detected spikes in this epoch (excluding
                        % those in window around stimulus artefact)
                        spikessamp = intersect(total_epoch,dataBase(subj).IEDs.spikesdet{ch});

                        % find sample number and calculate spikes pre and post stim
                        epochsamp = NaN(size(spikessamp)); removeSpike = [];
                        preCount = 0; postCount = 0;
                        for num = 1:size(spikessamp,1)
                            epochsamp(num) = find(total_epoch == spikessamp(num));

                            if t(epochsamp(num)) < 0 && t(epochsamp(num)) > -selDur
                                preCount = preCount+1;
                            elseif t(epochsamp(num)) > 0  && t(epochsamp(num)) < selDur
                                postCount = postCount+1;
                            else
                                % remove spike if it is not in the selected
                                % epoch (selDur around stimulus, excluding 
                                % the window around the stimulus artefact)
                                removeSpike = [removeSpike, num]; %#ok<AGROW>
                            end
                        end

                        selEpochsamp = epochsamp;
                        selEpochsamp(removeSpike) = [];

                        spikespre{stim,ch}(n) = preCount;
                        spikespost{stim,ch}(n) = postCount;

                        totEpochspikessamp{stim,ch,n} = epochsamp;
                        selEpochspikessamp{stim,ch,n} = selEpochsamp;

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

                    end

                    % if no spikes are found pre or post stim, one of the
                    % values == 0, which leads to Infinite values when
                    % determining the logarithmic value and difficulties in
                    % visualisation. So the 0-value is changed to 0.1.
                    if sum(spikespost{stim,ch}) == 0 && sum(spikespre{stim,ch}) ~= 0

                        spikesratio(stim,ch) = 0.1/sum(spikespre{stim,ch});

                    elseif sum(spikespost{stim,ch}) ~= 0 && sum(spikespre{stim,ch}) == 0

                        spikesratio(stim,ch) = sum(spikespost{stim,ch})/0.1;

                    elseif sum(spikespost{stim,ch}) == 0 && sum(spikespre{stim,ch}) == 0

                        spikesratio(stim,ch) = 0.1/0.1;

                    else
                        spikesratio(stim,ch) = sum(spikespost{stim,ch})/sum(spikespre{stim,ch});
%                     p_temp = signtest(spikespre{stim,ch},spikespost{stim,ch});
%                     p_spikes(stim,ch) = p_temp;
                    end
                end
            end
        end

        dataBase(subj).IEDs.totEpochspikessamp  = totEpochspikessamp; % complete time window (usually 4s: 2s pre and 2 s post-stim, excluding the window around the stimulus artefact)
        dataBase(subj).IEDs.selEpochspikessamp  = selEpochspikessamp; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
        dataBase(subj).IEDs.spikespre   = spikespre; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
        dataBase(subj).IEDs.spikespost  = spikespost; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
        dataBase(subj).IEDs.spikesratio = spikesratio; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
%         dataBase(subj).IEDs.spikes_p    = p_spikes;

    else
        % do nothing, because no spikes are annotated in this subject
    end
end

disp('Spike ratio is calculated.')

% housekeeping
clear ans ch epochsamp n num postCount preCount removeSpike samp_sel
clear samp_stim selDur selEpoch selEpochsamp selEpochspikessamp spikespost 
clear spikespre spikesratio spikessamp stim subj total_epoch totEpochspikessamp

%% end of script
% make figures used in article by running the scripts ccepSp03_FigX.m
