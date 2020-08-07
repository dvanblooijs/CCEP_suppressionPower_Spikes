%% script to use svm (constructed with makeSVM_TFSPES.m) to detect and visually power suppression in TFSPES plots

%% settings for using SVM
config_useSVM_powSup

%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg.sub_labels,2)
    
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    
    dataBase(subj).ERSP = load(fullfile(cfg.dir_ERSP, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
end

disp('All ERSPs are loaded')

%% detect power suppression



%% visually check TFSPES with detected power suppression


