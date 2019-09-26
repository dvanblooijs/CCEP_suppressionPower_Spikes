%% mainfile SVM_TFSPES

%% process the data to enable construction of TFSPES plots

preprocess_mainfile

%% compare visual ratings from 2 scorers



%% make SVM

% makeSVM_TFSPES --> needs to be adapted to be BIDS compatible!!

%% detectPowSup

cfg.inputERSP = '/Fridge/CCEP/derivatives/TFSPES/';
cfg.outputdetPowSup = '/Fridge/CCEP/derivatives/TFSPES/';
dataBase = detectPowSup(dataBase,cfg);


%% check detected Power Suppressions visually

