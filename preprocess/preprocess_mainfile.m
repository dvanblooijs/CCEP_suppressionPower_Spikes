%% mainfile CCEP_powerSuppression_ERs_spikes
% this code prepocesses BIDS data and analyzes the relationship between ERs, spikes and powerSuppression.
% author: D van Blooijs
% date: April 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath('git_rep/SPES_SOZ/detectERs')
addpath(genpath('git_rep/eeglab/'))     
addpath('git_rep/fieldtrip/')
ft_defaults

%% settings
cfg.dataPath = '/Fridge/CCEP';
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', 'sub-RESP0502',...
    'sub-RESP0574', 'sub-RESP0589', 'sub-RESP0608', 'sub-RESP0621', 'sub-RESP0699'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031153','run-051138','run-011714','run-021549','run-031740',...
    'run-021358','run-021050','run-021057','run-021147','run-031717'};
cfg.ERpath = '/Fridge/users/dorien/derivatives/BB_article/CCEPderiv';


%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(cfg);

%% load power suppression (PS) and ERs
for subj = 1:size(dataBase,2)
    load([cfg.ERpath, '/', cfg.sub_labels{subj}, '_ERs_BB.mat'])
    
    dataBase(subj).ERs_BB = stimp;
    
    if max([dataBase(subj).ERs_BB.stimchans]) > size(dataBase(subj).ch,1)
        fprintf('Warning: number of channels and stimulus pairs does not match! subject:%s, #channels:%1.0f, #stimpairs: %1.0f',sub_labels{subj},size(dataBase(subj).ch,1),max([dataBase(subj).ERs_BB.stimchans]))
    end
end

%% preprocessing CCEP in ECoG

% sort stimulation pairs
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

%% plot avg epoch
subj = 4;
trial = 2;
elec = 6;
fs = dataBase(subj).ccep_header.Fs;
tt =-round(cfg.epoch_prestim)+1/fs:1/fs:round((cfg.epoch_length-cfg.epoch_prestim));

figure(1),
plot(tt,squeeze(dataBase(subj).cc_epoch_sorted(elec,:,trial,:)));
hold on
plot(tt,squeeze(dataBase(subj).cc_epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
hold off
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('%s: Electrode %s, stimulating %s-%s',dataBase(subj).sub_label,dataBase(subj).ch{elec},dataBase(subj).cc_stimchans{trial,1},dataBase(subj).cc_stimchans{trial,2}))





