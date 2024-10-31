%% ccepSp02_detectSpikes
% author: Dorien van Blooijs
% date: july 2019

% this script is used to detect interictal epileptiform discharges (spikes)

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

    for run = 1:size(files_eeg,1)
        runTemp = extractBetween(files_eeg(run).name,'run-','_ieeg');
        cfg(nSubj).run_label{run} = ['run-', runTemp{1}];
    end
end

% housekeeping
clear files files_eeg files_ses files_subj idx_eeg idx_ses idx_subj nSubj run runTemp

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% determine channels with IEDs

dataBase = findIEDchannels(dataBase);

%% re-reference data to reduce artefacts

for nSubj = 1:size(dataBase,2)

    IED = dataBase(nSubj).IEDch;
    dataBase = rerefData(dataBase,IED,nSubj);
    
end

%% remove stimulation artefact

dataBase = removeStimArt(dataBase);

fprintf('...All subjects has been run...\n')

%% plot avg and rereferenced signal
% negative is negative, so signal is upside down compared to signal in
% Micromed

for nSubj = 1:size(dataBase,2)
    if ~isempty(dataBase(nSubj).IEDch)
        IED = dataBase(nSubj).IEDch(1);
    else
        IED = 1;
    end

    plot_DataRerefStimfree(dataBase,nSubj,IED)
end

%% detect spikes - only in IED channels

% those thresholds are determined in ccepSp01_trainTestSpikeDetection
thresh_t_dif = 0.1;
thresh_SD = 6;

for nSubj = 1:size(dataBase,2)
    
    if ~isempty(dataBase(nSubj).IEDch)
        data_IED = dataBase(nSubj).data_rerefnoStimArt(dataBase(nSubj).IEDch,:);
        dataBase(nSubj).data_IED = data_IED;
        fs = dataBase(nSubj).ccep_header.Fs;

        [M,norm_M, Pharmat_norm, Pharmatall_norm] = findMahalanobisDist(data_IED,fs);

        dataBase(nSubj).spikes.M = M;
        dataBase(nSubj).spikes.norm_M = norm_M;
        dataBase(nSubj).spikes.Pharmat_norm = Pharmat_norm;
        dataBase(nSubj).spikes.Pharmatall_norm = Pharmatall_norm;

        dataBase(nSubj).detIED = findCortSpikes(data_IED,fs,norm_M,Pharmat_norm,thresh_t_dif,thresh_SD);

        spikespat.sub_label = dataBase(nSubj).sub_label;
        spikespat.ses_label = dataBase(nSubj).ses_label;
        spikespat.run_label = dataBase(nSubj).run_label;
        spikespat.task_label = dataBase(nSubj).task_label;
        spikespat.spikesdet = dataBase(nSubj).detIED;
        spikespat.IEDchan = dataBase(nSubj).IEDchan;
        spikespat.IEDch = dataBase(nSubj).IEDch;
        spikespat.Pharmat_norm = Pharmat_norm;
        spikespat.thresh_t_dif = thresh_t_dif;
        spikespat.thresh_SD = thresh_SD;

        foldername = fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label);
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end

        filename = [dataBase(nSubj).sub_label,'_', dataBase(nSubj).ses_label, '_', dataBase(nSubj).task_label,'_', dataBase(nSubj).run_label,'_detIEDs.mat'];

        save(fullfile(foldername,filename),'spikespat')
    end
end


%% check spikes visually in the first 10 minutes

nSubj = 1;
IEDch = dataBase(nSubj).IEDchan';
IED = dataBase(nSubj).IEDch';
plot_SpikesData(dataBase,nSubj,IEDch, IED)







