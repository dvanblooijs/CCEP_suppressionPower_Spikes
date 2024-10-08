% CCEP_powerSuppression_cceps_spikes
% this code analyzes the relationship between cceps, spikes and
% powerSuppression.

% author: D van Blooijs
% date: April 2019

%% set paths
clear
close all
clc

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

% set other paths and get paths where data is collected and where
% derivatives will be saved
myDataPath = ccepSp_setLocalDataPath(1);

% housekeeping
clear matlabFolder RepoPath rootPath 

%% settings

files = dir(myDataPath.proj_dirinput);
idx_subj = contains({files(:).name},'sub-');
files_subj = files(idx_subj);
cfg = struct([]);

for nSubj = 1:size(files_subj,1)

    cfg(nSubj).sub_labels = files_subj(nSubj).name;

    files = dir(fullfile(files_subj(nSubj).folder,files_subj(nSubj).name));
    idx_ses = contains({files(:).name},'ses-');
    files_ses = files(idx_ses);

    cfg(nSubj).ses_label = files_ses(1).name;

    cfg(nSubj).task_label = 'task-SPESclin';

    files = dir(fullfile(files_ses(1).folder,files_ses(1).name,'ieeg'));
    idx_eeg = contains({files(:).name},'.eeg');
    files_eeg = files(idx_eeg);

    for run = 1:size(files_eeg,1)
        runTemp = extractBetween(files_eeg(run).name,'run-','_ieeg');
        cfg(nSubj).run_label{run} = ['run-', runTemp{1}];
    end
end

% housekeeping
clear files files_eeg files_ses files_subj idx_eeg idx_ses idx_subj run runTemp nSubj

%% load ECoGs with SPES from X subjects

dataBase = load_ECoGdata(myDataPath,cfg);

%% preprocessing CCEP in ECoG
cfg = [];

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% load CCEPs

files = dir(fullfile(myDataPath.proj_dirinput,'derivatives','CCEPs'));

for nSubj = 1:size(dataBase,2)

    idx_sub = contains({files(:).name}, dataBase(nSubj).sub_label) & [files(:).isdir] == 1;

    ses = dir(fullfile(files(idx_sub).folder,files(idx_sub).name));
    idx_ses = contains({ses(:).name}, dataBase(nSubj).ses_label);

    file = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name));
    idx_ccepfile = contains({file(:).name},'N1s');

    if sum(idx_ccepfile)>0
        dataBase(nSubj).ccep = load(fullfile(file(idx_ccepfile).folder,file(idx_ccepfile).name));
    end

end

disp('CCEPs are loaded')

% housekeeping
clear idx_ERSPfile idx_ccepfile file files idx_ses idx_sub ses nSubj ans

%% load ERSP

files = dir(fullfile(myDataPath.proj_dirinput,'derivatives','ERSP'));

for nSubj = 1:size(dataBase,2)

    idx_sub = contains({files(:).name}, dataBase(nSubj).sub_label) & [files(:).isdir] == 1;

    ses = dir(fullfile(files(idx_sub).folder,files(idx_sub).name));
    idx_ses = contains({ses(:).name}, dataBase(nSubj).ses_label);

    run = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name));
    idx_run = contains({run(:).name},dataBase(nSubj).run_label);

    file = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name,run(idx_run).name));
    idx_ERSPfile = contains({file(:).name},'ERSP.mat');

    if sum(idx_ERSPfile)>0
        dataBase(nSubj).ERSP = load(fullfile(file(idx_ERSPfile).folder,file(idx_ERSPfile).name));
    end
end

disp('ERSPs are loaded')

% housekeeping
clear idx_ERSPfile idx_ccepfile file files idx_ses idx_sub ses nSubj ans run idx_run

%% load spikes

files = dir(fullfile(myDataPath.proj_dirinput,'derivatives','IEDs'));

for nSubj = 1:size(dataBase,2)

    idx_sub = contains({files(:).name}, dataBase(nSubj).sub_label);

    if any(idx_sub)
        ses = dir(fullfile(files(idx_sub).folder,files(idx_sub).name));
        idx_ses = contains({ses(:).name}, dataBase(nSubj).ses_label);

        file = dir(fullfile(ses(idx_ses).folder,ses(idx_ses).name));
        idx_detIED = contains({file(:).name},'detIEDs');

        if sum(idx_detIED)>0
            load(fullfile(file(idx_detIED).folder,file(idx_detIED).name));
            dataBase(nSubj).IEDs = spikespat;
        end
    end
end

disp('Spikes are loaded')

clear nSubj ses spikespat idx_detIED idx_ses idx_sub file  files

%% make both the ccep and ersp matrix equal regarding stimulation pairs and response channels

for nSubj = 1:size(dataBase,2)

    % duplicate original cc_stimchans, because you will adapt those later
    % in this section of the script
    dataBase(nSubj).cc_stimchans_orig = dataBase(nSubj).cc_stimchans;
    dataBase(nSubj).cc_stimsets_orig = dataBase(nSubj).cc_stimsets;
    dataBase(nSubj).tt_epoch_sorted_orig = dataBase(nSubj).tt_epoch_sorted;

    CCEPmat = dataBase(nSubj).ccep.checked;
    ERSPmat = dataBase(nSubj).ERSP.checked;

    if size(CCEPmat,1) > size(CCEPmat,2) % [channels x stim pairs]
        CCEPmat = CCEPmat';
    end

    if size(ERSPmat,1) > size(ERSPmat,2) % [channels x stim pairs]
        ERSPmat = ERSPmat';
    end

    % only continue if the order of channels is equal in ccep and ersp
    % preprocessing
    if isequal(dataBase(nSubj).ccep.ch,dataBase(nSubj).ERSP.ch)

        % if the order of stimchans is also equal in ccep and ersp, than
        % those matrices can be used without further processing
        if isequal(dataBase(nSubj).ccep.cc_stimchans,dataBase(nSubj).ERSP.cc_stimchans)
            ERSPmat_comb = ERSPmat;
            CCEPmat_comb = CCEPmat;
            exclude = [];
        else
            % find total number of stimulation pairs and whether this one is
            % larger in ccep or in ersp preprocessing
            [stim_size,I] = max([size(CCEPmat,1), size(ERSPmat,1)]);

            % pre-allocation
            CCEPmat_comb = NaN(stim_size,size(dataBase(nSubj).ch,1)); % [stimulus pairs x channels]
            ERSPmat_comb = NaN(stim_size,size(dataBase(nSubj).ch,1));

            % when more (or equal) stims in CCEP than in ERSP
            if I == 1
            
                % for each stimulation pair
                for nStimp = 1:stim_size

                    % find where stimulus of CCEP can be found in ERSP
                    idx = find(strcmp(dataBase(nSubj).ERSP.cc_stimchans(:,1),dataBase(nSubj).ccep.cc_stimchans(nStimp,1)) & ...
                        strcmp(dataBase(nSubj).ERSP.cc_stimchans(:,2),dataBase(nSubj).ccep.cc_stimchans(nStimp,2)), 1);

                    % if this stimulus is found in both ERSP and CCEP data
                    if ~isempty(idx)

                        CCEPmat_comb(nStimp,:) = CCEPmat(nStimp,:);
                        ERSPmat_comb(nStimp,:) = ERSPmat(idx,:);

                        % the stimulated electrodes must be NaN
                        CCEPmat_comb(nStimp, dataBase(nSubj).ccep.cc_stimsets(nStimp,:)) = NaN;
                        ERSPmat_comb(nStimp, dataBase(nSubj).ERSP.cc_stimsets(idx,:)) = NaN;

                    else
                        % do nothing, this stimulus pair should not be
                        % included in the final matrices

                    end
                end
            elseif I == 2 % when more stims in ERSP than in CCEP

                % for each stimulation pair
                for nStimp = 1:stim_size

                    idx = find(strcmp(dataBase(nSubj).ccep.cc_stimchans(:,1),dataBase(nSubj).ERSP.cc_stimchans(nStimp,1)) & ...
                        strcmp(dataBase(nSubj).ccep.cc_stimchans(:,2),dataBase(nSubj).ERSP.cc_stimchans(nStimp,2)), 1);

                    % if this stimulus is found in both ERSP and CCEP data
                    if ~isempty(idx)

                        CCEPmat_comb(nStimp,:) = CCEPmat(idx,:);
                        ERSPmat_comb(nStimp,:) = ERSPmat(nStimp,:);

                        % the stimulated electrodes must be NaN
                        CCEPmat_comb(nStimp, dataBase(nSubj).ccep.cc_stimsets(idx,:)) = NaN;
                        ERSPmat_comb(nStimp, dataBase(nSubj).ERSP.cc_stimsets(nStimp,:)) = NaN;

                    else
                         % do nothing, this stimulus pair should not be
                        % included in the final matrices

                    end
                end
            end
           
            % exclude NaN-rows that exists both in CCEPmat_comb and
            % ERSPmat_comb
            exclude = all(isnan(CCEPmat_comb),2) & all(isnan(ERSPmat_comb),2);
            CCEPmat_comb(exclude,:) = [];
            ERSPmat_comb(exclude,:) = [];

            if size(ERSPmat_comb,1) == size(CCEPmat_comb,1) && size(ERSPmat_comb,2) == size(CCEPmat_comb,2)

            else
                error('Sizes of ERSP matrix and ccep matrix in subject %s are inequal!',dataBase(nSubj).sub_label)
            end           
        end

        cc_stimchans = dataBase(nSubj).cc_stimchans_orig;
        cc_stimsets = dataBase(nSubj).cc_stimsets_orig;
        tt_epoch_sorted = dataBase(nSubj).tt_epoch_sorted_orig;
        cc_stimchans(exclude,:) = [];
        cc_stimsets(exclude,:) = [];
        tt_epoch_sorted(:,exclude,:) = [];
        dataBase(nSubj).cc_stimchans = cc_stimchans;
        dataBase(nSubj).cc_stimsets = cc_stimsets;
        dataBase(nSubj).tt_epoch_sorted = tt_epoch_sorted;
        dataBase(nSubj).CCEPmat = CCEPmat_comb;
        dataBase(nSubj).ERSPmat = ERSPmat_comb;

    else
        error('Channels are not equal in ERSP and CCEP analysis')

    end
end

disp('ERSPmat and CCEPmat are equal in size')

% housekeeping
clear ERSPmat exclude CCEPmat ERSPmat_comb CCEPmat_comb nStimp exclude_ccep exclude_ersp I idx stim_size nSubj cc_stimchans cc_stimsets tt_epoch_sorted

%% find spikes in a specific time window for each response electrode 
% after stimulating each stimulus pair and calculate the spikesratio

% define period in which spikes are found for spikesratio
selectDur = 1; % seconds pre and post stimulation to find spikes
fs = dataBase(1).ccep_header.Fs;
t = -cfg.epoch_prestim+1/fs:1/fs : (cfg.epoch_length - cfg.epoch_prestim);

for nSubj = 1:size(dataBase,2)

    if ~isempty(dataBase(nSubj).IEDs)

        % pre-allocation
        totEpochspikessamp = cell(size(dataBase(nSubj).cc_stimchans,1),size(dataBase(nSubj).IEDs.IEDch,2),10); % [stimpairs x IED channels x trials]
        selectEpochspikessamp = cell(size(dataBase(nSubj).cc_stimchans,1),size(dataBase(nSubj).IEDs.IEDch,2),10); % [stimpairs x IED channels x trials]
        spikespre = cell(size(dataBase(nSubj).cc_stimchans,1),size(dataBase(nSubj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        spikespost = cell(size(dataBase(nSubj).cc_stimchans,1),size(dataBase(nSubj).IEDs.IEDch,2)); % [stimpairs x IED channels]
        spikesratio = NaN(size(dataBase(nSubj).cc_stimchans,1),size(dataBase(nSubj).IEDs.IEDch,2)); % [stimpairs x IED channels]

        for nStimp = 1:size(dataBase(nSubj).tt_epoch_sorted,2) % for each stim pair
            for nCh = 1:size(dataBase(nSubj).IEDs.IEDch,1) % for each IED channel
            
                % if channel is stimulated in this specific stim trial
                if strcmpi(dataBase(nSubj).cc_stimchans{nStimp,1},dataBase(nSubj).IEDs.IEDchan{nCh}) || ...
                        strcmpi(dataBase(nSubj).cc_stimchans{nStimp,2},dataBase(nSubj).IEDs.IEDchan{nCh})

                    % when channel with IEDs is stimulated, the spike ratio
                    % should be NaN
                    spikespre{nStimp,nCh}(1:size(dataBase(nSubj).cc_epoch_sorted,2)) = NaN;
                    spikespost{nStimp,nCh}(1:size(dataBase(nSubj).cc_epoch_sorted,2)) = NaN;
                    spikesratio(nStimp,nCh) = NaN;

                else
                    for nTrial = 1:size(dataBase(nSubj).cc_epoch_sorted,2) % for each trial in each stim pair (so usually 10)

                        % time window in samples of epoch (usually 4s: 2s pre-stim, 2s post-stim)
                        total_epoch = squeeze(dataBase(nSubj).tt_epoch_sorted(nTrial,nStimp,:));

                        % 200 ms around stimulus artefact should be not included in
                        % spike counts due to stimulus artefact (100 ms pre
                        % and 100 ms post stimulus)
                        samp_stim = round((cfg.epoch_prestim-0.1) * fs) : round((cfg.epoch_prestim+0.1) * fs);
                        total_epoch(samp_stim) = NaN; 

                        % find all detected spikes in this epoch (excluding
                        % those in window around stimulus artefact)
                        spikessamp = intersect(total_epoch,dataBase(nSubj).IEDs.spikesdet{nCh});

                        % find sample number and calculate spikes pre and post stim
                        % pre-allocation
                        epochsamp = NaN(size(spikessamp)); removeSpike = [];
                        preCount = 0; postCount = 0;

                        for nSpikes = 1:size(spikessamp,1)
                            epochsamp(nSpikes) = find(total_epoch == spikessamp(nSpikes));

                            if t(epochsamp(nSpikes)) < 0 && t(epochsamp(nSpikes)) > -selectDur
                                preCount = preCount+1;
                            elseif t(epochsamp(nSpikes)) > 0  && t(epochsamp(nSpikes)) < selectDur
                                postCount = postCount+1;
                            else
                                % remove spike if it is not in the selected
                                % epoch (selectDur around stimulus, excluding 
                                % the window around the stimulus artefact)
                                removeSpike = [removeSpike, nSpikes]; 
                            end
                        end

                        selectEpochsamp = epochsamp;
                        selectEpochsamp(removeSpike) = [];

                        spikespre{nStimp,nCh}(nTrial) = preCount;
                        spikespost{nStimp,nCh}(nTrial) = postCount;

                        totEpochspikessamp{nStimp,nCh,nTrial} = epochsamp;
                        selectEpochspikessamp{nStimp,nCh,nTrial} = selectEpochsamp;

                                            % figure(1),
                                            % plot(total_epoch,ones(size(total_epoch))),
                                            % hold on,
                                            % plot(dataBase(nSubj).IEDs.spikesdet{nCh},2*ones(1,size(dataBase(nSubj).IEDs.spikesdet{nCh},2)),'*'),
                                            % hold off,
                                            % xlim([min(total_epoch),max(total_epoch)]),
                                            % ylim([0 3])
                                            % title(sprintf('%s-%s, %s, trial %d',...
                                            %     dataBase(nSubj).cc_stimchans{nStimp,:},...
                                            %     dataBase(nSubj).ch{dataBase(nSubj).IEDs.IEDch(nCh)}, ...
                                            %     nTrial))

                        %pause

                    end

                    % if no spikes are found pre or post stim, one of the
                    % values == 0, which leads to Infinite values when
                    % determining the logarithmic value and difficulties in
                    % visualisation. So the 0-value is changed to 0.1.
                    if sum(spikespost{nStimp,nCh}) == 0 && sum(spikespre{nStimp,nCh}) ~= 0

                        spikesratio(nStimp,nCh) = 0.1/sum(spikespre{nStimp,nCh});

                    elseif sum(spikespost{nStimp,nCh}) ~= 0 && sum(spikespre{nStimp,nCh}) == 0

                        spikesratio(nStimp,nCh) = sum(spikespost{nStimp,nCh})/0.1;

                    elseif sum(spikespost{nStimp,nCh}) == 0 && sum(spikespre{nStimp,nCh}) == 0

                        spikesratio(nStimp,nCh) = 0.1/0.1;

                    else
                        spikesratio(nStimp,nCh) = sum(spikespost{nStimp,nCh})/sum(spikespre{nStimp,nCh});
                    end
                end
            end
        end

        dataBase(nSubj).IEDs.totEpochspikessamp  = totEpochspikessamp; % complete time window (usually 4s: 2s pre and 2 s post-stim, excluding the window around the stimulus artefact)
        dataBase(nSubj).IEDs.selEpochspikessamp  = selectEpochspikessamp; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
        dataBase(nSubj).IEDs.spikespre   = spikespre; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
        dataBase(nSubj).IEDs.spikespost  = spikespost; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)
        dataBase(nSubj).IEDs.spikesratio = spikesratio; % in the selected period (selDur around stimulus, excluding the window of the stimulus artefact)

    else
        % do nothing, because no spikes are annotated in this subject
    end
end

disp('Spike ratio is calculated.')

% housekeeping
clear ans nCh epochsamp nTrial nSpikes postCount preCount removeSpike samp_sel
clear samp_stim selectDur selEpoch selectEpochsamp selectEpochspikessamp spikespost 
clear spikespre spikesratio spikessamp nStimp nSubj total_epoch totEpochspikessamp

%% end of script
% make figures used in article by running the scripts ccepSp04_FigX.m
