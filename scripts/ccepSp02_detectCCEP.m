%% pipeline_CCEP
% author: Dorien van Blooijs
% date: September 2019

clc
clear
myDataPath = setLocalDataPath(1);

%% patient settings

files = dir(myDataPath.dataPath);
idx_subj = contains({files(:).name},'sub-');
files_subj = files(idx_subj);
cfg = struct([]);

for subj = 1:size(files_subj,1)

    cfg(subj).sub_labels = files_subj(subj).name;

    files = dir(fullfile(files_subj(subj).folder,files_subj(subj).name));
    idx_ses = contains({files(:).name},'ses-');
    files_ses = files(idx_ses);

    cfg(subj).ses_label = files_ses(1).name;

    cfg(subj).task_label = 'task-SPESclin';

    files = dir(fullfile(files_ses(1).folder,files_ses(1).name,'ieeg'));
    idx_eeg = contains({files(:).name},'.eeg');
    files_eeg = files(idx_eeg);

    for run = 1:size(files_eeg,1)
        runTemp = extractBetween(files_eeg(run).name,'run-','_ieeg');
        cfg(subj).run_label{run} = ['run-', runTemp{1}];
    end
end

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG
cfg = [];

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

%% detect CCEPs

dataBase = detect_n1peak_ECoG_ccep(dataBase, cfg);

disp('All CCEPs are detected')

%% visually check detected ccep
close all

% select subject
subs = {dataBase(:).sub_label};
string = [repmat('%s, ',1,size(subs,2)-1), '%s'];
substring = input(sprintf(['Choose subject: ',string,'\n'],subs{:}),'s');
subj = find(contains({dataBase(:).sub_label},substring));
    
% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.CCEPpath, dataBase(subj).sub_label,dataBase(subj).ses_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_N1sChecked.mat']),'file')
   dataBase(subj).ccep = load(fullfile(myDataPath.CCEPpath, dataBase(subj).sub_label,dataBase(subj).ses_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_N1sChecked.mat']));   
end

% continue with the stimulation pair after the last saved stimulation pair
if any(contains(fieldnames(dataBase(subj).ccep),'checkUntilStimp'))
    endstimp = dataBase(subj).ccep.checkUntilStimp;
else
    endstimp = 0;
end

dataBase = visualRating_ccep(myDataPath,dataBase,subj,cfg,endstimp);
    
%% save ccep

for subj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.CCEPpath, dataBase(subj).sub_label,dataBase(subj).ses_label);
    
    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    start_filename = strfind(dataBase(subj).dataName,'/');
    stop_filename = strfind(dataBase(subj).dataName,'_ieeg');
    
    fileName = [dataBase(subj).dataName(start_filename(end)+1:stop_filename-1),'_N1sChecked.mat'];
    
    ccep = dataBase(subj).ccep;
    [~,eegFile,eegExt] = fileparts(dataBase(subj).dataName);
    ccep.dataName = [eegFile, eegExt];
    ccep.ch = dataBase(subj).ch;
    ccep.cc_stimchans = dataBase(subj).cc_stimchans;
    ccep.cc_stimsets = dataBase(subj).cc_stimsets;
    
    % detection parameters must be described when known!!
    ccep.detpar.amplitude_thresh = cfg.amplitude_thresh;
    ccep.detpar.n1_peak_range = cfg.n1_peak_range;

    save(fullfile(targetFolder,fileName), '-struct','ccep');
end



