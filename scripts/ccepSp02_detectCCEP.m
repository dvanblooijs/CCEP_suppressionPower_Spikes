%% pipeline_CCEP
% author: Dorien van Blooijs
% date: September 2019

clc
clear
cfg = setLocalDataPath(1);

%% patient settings
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', 'sub-RESP0502',...
    'sub-RESP0574', 'sub-RESP0589', 'sub-RESP0608', 'sub-RESP0621', 'sub-RESP0699'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031153','run-051138','run-011714','run-021549','run-031740'...
    'run-021358','run-021050','run-021057','run-021147','run-031717'};

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

cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

dataBase = preprocess_ECoG_ccep(dataBase,cfg);

disp('All ECoGs are preprocessed')

%% detect ERs

dataBase = detect_n1peak_ECoG_ccep(dataBase, cfg);

disp('All ERs are detected')


%% visually check detected ccep

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
    detpar.thresh = cfg.thresh;
    detpar.minSD = cfg.minSD;

    save(fullfile(targetFolder,fileName), 'ERs','dataName','ch','detpar');
end



