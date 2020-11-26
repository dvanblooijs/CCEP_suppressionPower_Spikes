%% mainfile_TFSPES
% author: Dorien van Blooijs
% date: july 2019

clc
clear
cfg = setLocalDataPath(1);

%% patient settings

cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPESclin';
cfg.run_label = {['run-' input('Run [daydayhhminmin]: ','s')]};

%% load ECoGs with SPES from X patients

dataBase = load_ECoGdata(cfg);

%% preprocessing CCEP in ECoG

% sort stimulation pairs
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% make TF-SPES Event-Related - Stimulus - Perturbations

close all
cfg.saveERSP = 'yes';

dataBase = makeTFSPES(dataBase,cfg);



