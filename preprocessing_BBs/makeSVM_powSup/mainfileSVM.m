%% mainfile SVM_TFSPES

%% settings for constructing SVM
config_makeSVM_TFSPES

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

%% load visual ratings

dataBase = load_visual_BBs(dataBase,cfg);
disp('All visual scores are loaded')

%% compare visual ratings from 2 scorers

% something with kappa?

%% make SVM

% makeSVM_TFSPES --> needs to be adapted to be BIDS compatible!!

%% detectPowSup

cfg.inputERSP = '/Fridge/CCEP/derivatives/TFSPES/';
cfg.outputdetPowSup = '/Fridge/CCEP/derivatives/TFSPES/';
dataBase = detectPowSup(dataBase,cfg);


%% check detected Power Suppressions visually

