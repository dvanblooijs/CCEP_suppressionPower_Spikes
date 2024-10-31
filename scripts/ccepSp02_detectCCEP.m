%% pipeline_CCEP
% author: Dorien van Blooijs
% date: September 2019

% This script is used to detect CCEPs. 
% 1. set paths
% 2. select subjects
% 3. load ecog data (BIDS)
% 4. epoch the data in [-2,2] around each single pulse.
% 5. Detect CCEPs in each averaged epoch
% 6. Visually check the detected CCEPs. 
% 7. Save this in N1sChecked.mat

%% set paths

clear
close all
clc

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

% set other paths and get paths where data is collected and where
% derivatives will be saved
myDataPath = ccepSp_setLocalDataPath(1);

% housekeeping
clear matlabFolder RepoPath rootPath 

%% patient settings

files = dir(myDataPath.proj_dirinput);
idx_subj = contains({files(:).name},'sub-');
files_subj = files(idx_subj);
cfg = struct([]);

for nSubj = 1:size(files_subj,1)

    cfg(nSubj).sub_labels = files_subj(nSubj).name;

    files = dir(fullfile(files_subj(nSubj).folder,files_subj(nSubj).name));
    idx_ses = contains({files(:).name},'ses-');
    files_ses = files(idx_ses);

    cfg(nSubj).ses_label = files_ses(1).name;

    cfg(nSubj).task_label = 'task-SPESclin';

    files = dir(fullfile(files_ses(1).folder,files_ses(1).name,'ieeg'));
    idx_eeg = contains({files(:).name},'.eeg');
    files_eeg = files(idx_eeg);

    for nRun = 1:size(files_eeg,1)
        runTemp = extractBetween(files_eeg(nRun).name,'run-','_ieeg');
        cfg(nSubj).run_label{nRun} = ['run-', runTemp{1}];
    end
end

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG
cfg = [];

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
nSubj = find(contains({dataBase(:).sub_label},substring));
    
% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label,...
        [dataBase(nSubj).sub_label, '_', dataBase(nSubj).ses_label,'_',dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label,'_N1sChecked.mat']),'file')

   dataBase(nSubj).ccep = load(fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label,...
        [dataBase(nSubj).sub_label, '_', dataBase(nSubj).ses_label,'_',dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label,'_N1sChecked.mat']));   
end

% continue with the stimulation pair after the last saved stimulation pair
if any(contains(fieldnames(dataBase(nSubj).ccep),'checkUntilStimp'))
    endstimp = dataBase(nSubj).ccep.checkUntilStimp;
else
    endstimp = 0;
end

dataBase = visualRating_ccep(myDataPath,dataBase,nSubj,cfg,endstimp);
    
%% save ccep

for nSubj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.proj_diroutput, dataBase(nSubj).sub_label);
    
    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    start_filename = strfind(dataBase(nSubj).dataName,'/');
    stop_filename = strfind(dataBase(nSubj).dataName,'_ieeg');
    
    fileName = [dataBase(nSubj).dataName(start_filename(end)+1:stop_filename-1),'_N1sChecked.mat'];
    
    ccep = dataBase(nSubj).ccep;
    [~,eegFile,eegExt] = fileparts(dataBase(nSubj).dataName);
    ccep.dataName = [eegFile, eegExt];
    ccep.ch = dataBase(nSubj).ch;
    ccep.cc_stimchans = dataBase(nSubj).cc_stimchans;
    ccep.cc_stimsets = dataBase(nSubj).cc_stimsets;
    
    % detection parameters must be described when known!!
    ccep.detpar.amplitude_thresh = cfg.amplitude_thresh;
    ccep.detpar.n1_peak_range = cfg.n1_peak_range;

    save(fullfile(targetFolder,fileName), '-struct','ccep');
end



