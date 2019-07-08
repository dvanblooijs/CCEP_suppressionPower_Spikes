%% mainfile_TFSPES
% author: Dorien van Blooijs
% date: july 2019

%% process the data to enable construction of TFSPES plots

preprocess_mainfile

%% make TF-SPES figures

cfg.output = '/Fridge/CCEP/derivatives/TFSPES/';
cfg.saveERSP = 'no';

dataBase = makeTFSPES(dataBase,cfg);

