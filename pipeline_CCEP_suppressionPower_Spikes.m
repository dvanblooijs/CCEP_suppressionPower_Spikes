%% pipeline CCEP_suppressionPower_Spikes
% author: Dorien van Blooijs
% date: September 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath(genpath('git_rep/eeglab/'))
addpath(genpath('git_rep/BasicCode_ECoG_DvB/'))
addpath(genpath('git_rep/REC2Stim/'))
addpath('git_rep/fieldtrip/')
ft_defaults

%% patient settings
cfg.dataPath = '/Fridge/CCEP';
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', 'sub-RESP0502',...
    'sub-RESP0574', 'sub-RESP0589', 'sub-RESP0608', 'sub-RESP0621', 'sub-RESP0699'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.ERpath = '/Fridge/users/dorien/derivatives/BB_article/CCEPderiv';


%% load detected ERs 
for subj = 1:size(dataBase,2)
    if exist(fullfile(cfg.ERpath,dataBase(subj).sub_label,dataBase(subj).ses_label,...
            [dataBase(subj).sub_label,'_', dataBase(subj).ses_label,'_',dataBase(subj).task_label, '_',dataBase(subj).run_label, '_ERs.mat']),'file')
        
        load(fullfile(cfg.ERpath,dataBase(subj).sub_label,dataBase(subj).ses_label,...
            [dataBase(subj).sub_label,'_', dataBase(subj).ses_label,'_',dataBase(subj).task_label, '_',dataBase(subj).run_label, '_ERs.mat']));
        
        dataBase(subj).ERs = ERs;
    else
        fprintf('Run pipeline_CCEP for patient %s',dataBase(subj).sub_label)
    end
end


%% load detected BB
for subj = 1:size(dataBase,2)
    if exist(fullfile(cfg.ERpath,dataBase(subj).sub_label,dataBase(subj).ses_label,...
            [dataBase(subj).sub_label,'_', dataBase(subj).ses_label,'_',dataBase(subj).task_label, '_',dataBase(subj).run_label, '_BBs.mat']),'file')
        
        load(fullfile(cfg.ERpath,dataBase(subj).sub_label,dataBase(subj).ses_label,...
            [dataBase(subj).sub_label,'_', dataBase(subj).ses_label,'_',dataBase(subj).task_label, '_',dataBase(subj).run_label, '_BBs.mat']));
        
        dataBase(subj).BBs = BBs;
    else
        fprintf('Run pipeline_BB for patient %s',dataBase(subj).sub_label)
    end
end


%% load detected spikes --> nog maken

for subj = 1:size(dataBase,2)
    if exist(fullfile(cfg.ERpath,...
            [dataBase(subj).sub_label, '_detSpikes.mat']),'file')
        
        load(fullfile(cfg.ERpath,...
            [dataBase(subj).sub_label, '_detSpikes.mat']));
        
        dataBase(subj).ERs_BB = stimp;
    else
        fprintf('Run pipeline_spikes for patient %s\n',dataBase(subj).sub_label)
    end
end

%% merge ERs, BBs and detected spikes so that the order of stimulus pairs is equal for all situations


%% analysis