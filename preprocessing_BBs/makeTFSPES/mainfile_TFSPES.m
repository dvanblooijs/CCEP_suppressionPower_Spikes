%% mainfile_TFSPES
% author: Dorien van Blooijs
% date: july 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
% addpath('git_rep/SPES_SOZ/detectERs')
addpath(genpath('git_rep/eeglab/'))     
addpath('git_rep/fieldtrip/')
ft_defaults

%% patient settings

cfg.dataPath = '/Fridge/CCEP';
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
% cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', 'sub-RESP0502',...
%     'sub-RESP0574', 'sub-RESP0589', 'sub-RESP0608', 'sub-RESP0621', 'sub-RESP0699'};
% train SVM detection for BB
% old database: PAT119, PAT126, PAT130, PAT135
cfg.sub_labels = { 'sub-RESP0607'};%, 'sub-RESP0638', 'sub-RESP0676', 'sub-RESP0690'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031211'};% {'run-021120','run-031049','run-021423','run-041139'};
cfg.ERpath = '/Fridge/users/dorien/derivatives/BB_article/CCEPderiv';

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

cfg.output = '/Fridge/users/dorien/derivatives/TFSPES/';
cfg.saveERSP = 'yes';

dataBase = makeTFSPES(dataBase,cfg);



