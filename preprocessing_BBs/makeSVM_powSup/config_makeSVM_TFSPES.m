%% config makeSVM_TFSPES

%% set paths

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath('git_rep/SPES_SOZ/detectERs')
addpath(genpath('git_rep/eeglab/'))
addpath('git_rep/fieldtrip/')
ft_defaults

%% set configuration
clear

cfg.dataPath = '/Fridge/CCEP';
% old database: PAT119,  PAT126, PAT130, PAT135
cfg.sub_labels = { 'sub-RESP0607', 'sub-RESP0638','sub-RESP0676','sub-RESP0690'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031211','run-031049','run-021423','run-041139'};
cfg.train = 1:3;
cfg.test = 4;

% Directory ERSP TFSPES figures
cfg.dir_ERSP = '/Fridge/users/dorien/derivatives/TFSPES/';

% Directory visual ratings
cfg.dir_visrate = '/Fridge/users/dorien/derivatives/BB_article/BB_visrating';
