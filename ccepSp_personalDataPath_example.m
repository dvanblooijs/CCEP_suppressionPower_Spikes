%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own datapaths

function localDataPath = ccepSp_personalDataPath_example()

% function that contains local data path, is ignored in .gitignore

localDataPath.proj_dirinput = '/datapath/to/data/on/openneuro/';
localDataPath.proj_diroutput = '/where/you/would/like/to/save/derivatives/';
localDataPath.figures = '/where/you/would/like/to/save/figures/';

addpath(genpath('/blabla/git_rep/CCEP_suppressionPower_Spikes'))
addpath(genpath('/blabla/git_rep/eeglab/'))
addpath('/blabla/git_rep/fieldtrip/')
ft_defaults

end