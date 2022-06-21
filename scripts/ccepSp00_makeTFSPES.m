%% mainfile_TFSPES
% author: Dorien van Blooijs
% date: july 2019

clc
clear
myDataPath = setLocalDataPath(1);

%% patient settings

subFolders = dir(myDataPath.dataPath);
idx_subjects = contains({subFolders(:).name},'sub-');
subjects = {subFolders(idx_subjects).name};
string = [repmat('%s, ',1,size(subjects,2)-1),'%s'];

cfg.sub_labels = {input(sprintf(['Select one of these files [',string,']: \n'],subjects{:}),'s')};

sesFiles = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{:}));
idx_ses = contains({sesFiles(:).name},'ses-');
cfg.ses_label = sesFiles(idx_ses).name;

cfg.task_label = 'task-SPESclin';

runFiles = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{:}, ...
    cfg.ses_label,'ieeg'));
idx_run = contains({runFiles(:).name},'.eeg');
temp_run_label = extractBetween(runFiles(idx_run).name,'run-','_ieeg');
cfg.run_label = {['run-' temp_run_label{1}]};

%% load ECoGs with SPES from X patients

dataBase = load_ECoGdata(myDataPath,cfg);

%% preprocessing CCEP in ECoG

% sort stimulation pairs
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% rereference data
cfg.reref = 1; % (1 = re-reference, 0 = no re-reference)

for subj=1:size(dataBase,2)
    
    if cfg.reref == 1
        for stimp = 1:size(dataBase(subj).cc_epoch_sorted,3)
            
            for numstim =1:size(dataBase(subj).cc_epoch_sorted,2)
                
                % find 10 signals with lowest variance and not being a bad channel or part stimulus pair
                variance = var(squeeze(dataBase(subj).cc_epoch_sorted(:,numstim,stimp,:)),1,2);
                [~,idx_var] = sort(variance,'ascend');
                
                idx_var = setdiff(idx_var,[find(strcmp(dataBase(subj).tb_channels.status,'bad'));dataBase(subj).cc_stimsets(stimp,:)'],'stable');
                
                ref = median(squeeze(dataBase(subj).cc_epoch_sorted(idx_var(1:10),numstim,stimp,:)));
                
                dataBase(subj).cc_epoch_sorted_reref(:,numstim,stimp,:) = squeeze(dataBase(subj).cc_epoch_sorted(:,numstim,stimp,:)) - ref;
                
            end
        end
    else
        
        dataBase(subj).cc_epoch_sorted_reref = dataBase(subj).cc_epoch_sorted;        
    end
    
    dataBase(subj).cc_epoch_sorted_avg = squeeze(mean(dataBase(subj).cc_epoch_sorted_reref,2,'omitnan'));
end

%% make TF-SPES Event-Related - Stimulus - Perturbations

close all
cfg.saveERSP = 'yes';

dataBase = makeTFSPES(dataBase,myDataPath,cfg);



