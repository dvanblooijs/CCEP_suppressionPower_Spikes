%% pipeline_CCEP
% author: Dorien van Blooijs
% date: September 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath(genpath('git_rep/eeglab/'))
% addpath(genpath('git_rep/BasicCode_ECoG_DvB/'))
% addpath(genpath('git_rep/REC2Stim/'))
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

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG

% sort stimulation pairs
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% detect ERs

dataBase = detectERs(dataBase,cfg);

disp('All ERs are detected')


%% visually check detected cceps

dataBase = visualRating_ccep(dataBase,cfg);

%% save ccep

for subj=1:size(dataBase,2)
    targetFolder = fullfile(cfg.ERpath, dataBase(subj).sub_label,dataBase(subj).ses_label);
    
    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    start_filename = strfind(dataBase(subj).dataName,'/');
    stop_filename = strfind(dataBase(subj).dataName,'_ieeg');
    
    fileName=[dataBase(subj).dataName(start_filename(end)+1:stop_filename-1),'_ERs.mat'];
    
    ERs = dataBase(subj).ERs;
    dataName = dataBase(subj).dataName;
    ch = dataBase(subj).ch;
    % detection parameters must be described when known!!
    detpar.thresh = NaN;
    detpar.minSD = NaN;
    detpar.sel = NaN;
    
    save(fullfile(targetFolder,fileName), 'ERs','dataName','ch','detpar');
end



