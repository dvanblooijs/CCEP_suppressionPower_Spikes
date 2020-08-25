%% spike detection
% author: Dorien van Blooijs
% date: july 2019

clc
clear
cfg = setLocalDataPath(1);

%% settings
cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', ...
    'sub-RESP0574',  'sub-RESP0608',  'sub-RESP0699'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031153','run-051138','run-011714','run-021549',...
    'run-021358','run-021057','run-031717'};

%% load ECoGs with SPES from 10 patients

dataBase = load_ECoGdata(cfg);

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

for subj = 1:size(dataBase,2)
    IED = dataBase(subj).IEDch(1);
    
    plot_DataRerefStimfree(dataBase,subj,IED)
end

%% detect spikes - only in IED channels

thresh_t_dif = 0.1;
thresh_SD = 6;

for subj = 1:size(dataBase,2)
    
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
    
    foldername = fullfile(cfg.visIEDpath,dataBase(subj).sub_label, dataBase(subj).ses_label);
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end
    
    filename = [dataBase(subj).sub_label,'_', dataBase(subj).ses_label, '_', dataBase(subj).task_label,'_', dataBase(subj).run_label,'_detIEDs.mat'];

    save([foldername,filename],'spikespat')

end


%% check spikes visually in the first 10 minutes

subj = 5;
IEDch = dataBase(subj).IEDchan;
IED = dataBase(subj).IEDch;
plot_SpikesData(dataBase,subj,IEDch, IED)







