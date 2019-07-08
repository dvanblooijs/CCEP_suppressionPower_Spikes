%% mainfile CCEP_powerSuppression_ERs_spikes
% this code prepocesses BIDS data and analyzes the relationship between ERs, spikes and powerSuppression.
% author: D van Blooijs
% date: April 2019

addpath('Desktop/git_rep/CCEP_suppressionPower_Spikes')
addpath('Desktop/git_rep/CCEP_suppressionPower_Spikes/analysis_ERs_PS_spikes')
addpath('Desktop/git_rep/CCEP_suppressionPower_Spikes/makeSpikeDet')
addpath('Desktop/git_rep/CCEP_suppressionPower_Spikes/makeSVM_powSup')
addpath('Desktop/git_rep/CCEP_suppressionPower_Spikes/makeTFSPES')
addpath('Desktop/git_rep/CCEP_suppressionPower_Spikes/preprocess')
addpath('Desktop/git_rep/SPES_SOZ/detectERs')
addpath(genpath('Desktop/git_rep/eeglab/'))     
addpath('Desktop/git_rep/fieldtrip/')
ft_defaults

%% settings
cfg.dataPath = '/Fridge/CCEP';
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
cfg.sub_labels = { 'RESP0401', 'RESP0435', 'RESP0458', 'RESP0478', 'RESP0502',...
    'RESP0574', 'RESP0589', 'RESP0608', 'RESP0621', 'RESP0699'};
cfg.ses_label = '1';
cfg.task_label = 'SPESclin';
cfg.ERpath = '/home/dorien/local_drives/CCEPderiv';


%% load ECoGs with SPES from 10 patients

dataBase = load_data(cfg);

%% load power suppression (PS) and ERs
for subj = 1:size(dataBase,2)
    
    load([cfg.ERpath, '/', cfg.sub_labels{subj}, '_ERs_BB.mat'])
    
    dataBase(subj).ERs_BB = stimp;
    
    if max([dataBase(subj).ERs_BB.stimchans]) > size(dataBase(subj).ch,1)
        fprintf('Warning: number of channels and stimulus pairs does not match! subject:%s, #channels:%1.0f, #stimpairs: %1.0f',sub_labels{subj},size(dataBase(subj).ch,1),max([dataBase(subj).ERs_BB.stimchans]))
    end
end

%% sort stimulation pairs

% if you want to take negative/positive stimulation into account
cfg.dir = 'no';
% if you want to take stimulation current into account
cfg.amp = 'no';

dataBase = unique_stimpairs(dataBase,cfg);

%% select epochs and average

cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = select_epochs(dataBase,cfg);

%% plot avg epoch
subj = 4;
trial = 2;
elec = 6;
tt =-round(cfg.epoch_prestim*dataBase(subj).ccep_header.Fs)+1:round((cfg.epoch_length-cfg.epoch_prestim)*dataBase(subj).ccep_header.Fs);

figure(1),
plot(tt,squeeze(dataBase(subj).cc_epoch_sorted(elec,:,trial,:)));
hold on
plot(tt,squeeze(dataBase(subj).cc_epoch_sorted_avg(elec,trial,:)),'k','linewidth',2);
hold off
xlabel('time(s)')
ylabel('amplitude(uV)')
title(sprintf('%s: Electrode %s, stimulating %s-%s',dataBase(subj).subj,dataBase(subj).ch{elec},dataBase(subj).cc_stimchans{trial,1},dataBase(subj).cc_stimchans{trial,2}))





