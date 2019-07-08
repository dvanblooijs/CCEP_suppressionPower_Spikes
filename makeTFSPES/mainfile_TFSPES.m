%% mainfile_TFSPES

%% process the data to enable construction of TFSPES plots

preprocess_mainfile

%% make TF-SPES figures

cfg.output = '/Fridge/CCEP/derivatives/TFSPES/';
cfg.saveERSP = 'no';

dataBase = makeTFSPES(dataBase,cfg);

