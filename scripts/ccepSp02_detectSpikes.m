%% spike detection
% author: Dorien van Blooijs
% date: july 2019

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

%% determine channels with IEDs

dataBase = findIEDchannels(dataBase);

%% re-reference data to reduce artefacts

for subj = 1:size(dataBase,2)

    IED = dataBase(subj).IEDch;
    dataBase = rerefData(dataBase,IED,subj);
    
end

%% remove stimulation artefact

dataBase = removeStimArt(dataBase);

fprintf('...All subjects has been run...\n')

%% plot avg and rereferenced signal
% negative is negative, so signal is upside down compared to signal in
% Micromed

for subj = 1:size(dataBase,2)
    if ~isempty(dataBase(subj).IEDch)
        IED = dataBase(subj).IEDch(1);
    else
        IED = 1;
    end

    plot_DataRerefStimfree(dataBase,subj,IED)
end

%% detect spikes - only in IED channels

% those thresholds are determined in ccepSp01_trainTestSpikeDetection
thresh_t_dif = 0.1;
thresh_SD = 6;

for subj = 1:size(dataBase,2)
    
    if ~isempty(dataBase(subj).IEDch)
        data_IED = dataBase(subj).data_rerefnoStimArt(dataBase(subj).IEDch,:);
        dataBase(subj).data_IED = data_IED;
        fs = dataBase(subj).ccep_header.Fs;

        [M,norm_M, Pharmat_norm, Pharmatall_norm] = findMahalanobisDist(data_IED,fs);

        dataBase(subj).spikes.M = M;
        dataBase(subj).spikes.norm_M = norm_M;
        dataBase(subj).spikes.Pharmat_norm = Pharmat_norm;
        dataBase(subj).spikes.Pharmatall_norm = Pharmatall_norm;

        dataBase(subj).detIED = findCortSpikes(data_IED,fs,norm_M,Pharmat_norm,thresh_t_dif,thresh_SD);

        spikespat.sub_label = dataBase(subj).sub_label;
        spikespat.ses_label = dataBase(subj).ses_label;
        spikespat.run_label = dataBase(subj).run_label;
        spikespat.task_label = dataBase(subj).task_label;
        spikespat.spikesdet = dataBase(subj).detIED;
        spikespat.IEDchan = dataBase(subj).IEDchan;
        spikespat.IEDch = dataBase(subj).IEDch;
        spikespat.Pharmat_norm = Pharmat_norm;
        spikespat.thresh_t_dif = thresh_t_dif;
        spikespat.thresh_SD = thresh_SD;

        foldername = fullfile(myDataPath.visIEDpath,dataBase(subj).sub_label, dataBase(subj).ses_label);
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end

        filename = [dataBase(subj).sub_label,'_', dataBase(subj).ses_label, '_', dataBase(subj).task_label,'_', dataBase(subj).run_label,'_detIEDs.mat'];

        save([foldername,'/',filename],'spikespat')
    end
end


%% check spikes visually in the first 10 minutes

subj = 1;
IEDch = dataBase(subj).IEDchan';
IED = dataBase(subj).IEDch';
plot_SpikesData(dataBase,subj,IEDch, IED)







