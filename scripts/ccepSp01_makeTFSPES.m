%% ccepSP01_makeTFSPES
% author: Dorien van Blooijs
% date: july 2019

% with this script, time-frequency plots are constructed around single
% pulse electrical stimuli

% 1. set paths
% 2. select subject
% 3. load ecog data (BIDS)
% 4. epoch the data in [-2,2] around each single pulse.
% 5. re-reference the data
% 6. calculate event-related spectral perturbation (ERSP) plots
% 7. save indiviual figures and a matrix in ERSP.mat

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

%% patient settings

subFolders = dir(myDataPath.proj_dirinput);
idx_subjects = contains({subFolders(:).name},'sub-');
subjects = {subFolders(idx_subjects).name};
string = [repmat('%s, ',1,size(subjects,2)-1),'%s'];

cfg.sub_labels = input(sprintf(['Select one of these subjects: [',string,']: \n'],subjects{:}),'s');

sesFiles = dir(fullfile(myDataPath.proj_dirinput,cfg.sub_labels));
idx_ses = contains({sesFiles(:).name},'ses-');
cfg.ses_label = sesFiles(idx_ses).name;

cfg.task_label = 'task-SPESclin';

runFiles = dir(fullfile(myDataPath.proj_dirinput,cfg.sub_labels, ...
    cfg.ses_label,'ieeg'));
idx_run = contains({runFiles(:).name},'.eeg');
temp_run_label = extractBetween(runFiles(idx_run).name,'run-','_ieeg');
cfg.run_label = {['run-' temp_run_label{1}]};

%% load ECoGs with SPES from selected subject(s)

dataBase = load_ECoGdata(myDataPath,cfg);

%% preprocessing CCEP in ECoG

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% rereference data
cfg.reref = 1; % (1 = re-reference, 0 = no re-reference)

for nSubj = 1:size(dataBase,2)
    
    if cfg.reref == 1
        for nStimp = 1:size(dataBase(nSubj).cc_epoch_sorted,3) % for each stimulus pair
            
            for nTrial = 1:size(dataBase(nSubj).cc_epoch_sorted,2) % for each trial of each stimulus pair
                
                % find 10 signals with lowest variance and not being a bad channel or part stimulus pair
                variance = var(squeeze(dataBase(nSubj).cc_epoch_sorted(:,nTrial,nStimp,:)),1,2);
                [~,idx_var] = sort(variance,'ascend');
                
                idx_var = setdiff(idx_var,[find(strcmp(dataBase(nSubj).tb_channels.status,'bad'));dataBase(nSubj).cc_stimsets(nStimp,:)'],'stable');
                
                ref = median(squeeze(dataBase(nSubj).cc_epoch_sorted(idx_var(1:10),nTrial,nStimp,:)));
                
                dataBase(nSubj).cc_epoch_sorted_reref(:,nTrial,nStimp,:) = squeeze(dataBase(nSubj).cc_epoch_sorted(:,nTrial,nStimp,:)) - ref;
                
            end
        end
    else
        
        dataBase(nSubj).cc_epoch_sorted_reref = dataBase(nSubj).cc_epoch_sorted;        
    end
    
    dataBase(nSubj).cc_epoch_sorted_avg = squeeze(mean(dataBase(nSubj).cc_epoch_sorted_reref,2,'omitnan'));
end

%% make TF-SPES Event-Related - Stimulus - Perturbations

close all
cfg.saveERSP = 'yes';

dataBase = makeTFSPES(dataBase,myDataPath,cfg);



