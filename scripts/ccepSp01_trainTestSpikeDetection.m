
clc
clear
cfg = setLocalDataPath(1);

%% settings

cfg.sub_labels = { 'sub-RESP0638', 'sub-RESP0639', 'sub-RESP0677', 'sub-RESP0698'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031049','run-021222','run-030954','run-021707'};

%% load ECoGs with SPES from 4 patients in whom spikes are annotated during 10 minutes

dataBase = load_ECoGdata(cfg);

%% load scored IEDs

for subj = 1:size(dataBase,2)
    
    fileName = [cfg.visIEDpath, dataBase(subj).sub_label, '_', dataBase(subj).ses_label, '_',...
        dataBase(subj).task_label, '_', dataBase(subj).run_label, '_visIEDs.mat'];
    
    if exist(fileName,'file')
        dataBase(subj).visIEDs = load(fileName);
    end
end

disp('Visually scored spikes are loaded')


%% re-reference data to reduces artefacts
for subj = 1:size(dataBase,2)

    IED = dataBase.IEDch; % CHECK! dit klopt wellicht niet

    dataBase = rerefData(dataBase,subj,IED);
end

%% remove stimulation artefact

dataBase = removeStimArt(dataBase);

fprintf('...All subjects has been run...\n')

%% plot avg and rereferenced signal

for subj = 1:size(dataBase,2)
    IED = dataBase(subj).IEDch(1); % CHECK! dit klopt wellicht niet
    
    plot_DataRerefStimfree(dataBase,subj,IED)
end

%% detect Mvalues

for subj = 1:size(dataBase,2)
    
    data_IED = dataBase(subj).data_rerefnoStimArt(dataBase(subj).visIEDs.IED,:);
    
    dataBase(subj).data_IED = data_IED;
    fs = dataBase(subj).ccep_header.Fs;

    [M,norm_M, Pharmat_norm, Pharmatall_norm] = findMahalanobisDist(data_IED,fs);
    
    dataBase(subj).spikes.M = M;
    dataBase(subj).spikes.norm_M = norm_M;
    dataBase(subj).spikes.Pharmat_norm = Pharmat_norm;
    dataBase(subj).spikes.Pharmatall_norm = Pharmatall_norm;
    
end

%% train spike detector on 5 minutes of scored IEDs

train_thresh = 1:50;
thresh_t_dif = 0.1:0.1:1;
train_threshold = struct;

for subj = 1:size(dataBase,2)
    
    [sens,prec,F] = trainSpikeDetector(dataBase,subj,train_thresh,thresh_t_dif);
    
    train_threshold(subj).F = F;
    train_threshold(subj).sens = sens;
    train_threshold(subj).prec = prec;

end

disp('All thresholds are varied')

%% plot showing the differences in relative and absolute threshold values

plot_threshSpikeDetector(dataBase,train_threshold,train_thresh,thresh_t_dif)

%% combine all train thresholds for calculation of overall performance

mode = {'sens','prec','F'};

for i=1:size(mode,2)
    
    % pre-allocation
    performance = NaN(size(dataBase,2),size(train_thresh,2),size(thresh_t_dif,2),2); %[subj, SD thresh, t thresh, chan/pat]
    
    for subj = 1:size(dataBase,2)
        for chanpat=1:2 % 1= per channel, 2 = per patient
            performance(subj,:,:,chanpat) = train_threshold(subj).(mode{i})(:,:,chanpat);
        end
    end
    train_thresholds_all.(mode{i}) = performance;
    
    % median of performance per patient
    train_thresholds_all.([mode{i},'med']) = squeeze(median(performance,1));
end

%% plot heatmaps per patient and patients combined

plot_HeatmapThreshSpikeDetector(dataBase,train_thresholds_all,thresh_t_dif,train_thresh)

%% values with max performance

display_optimalValSpikeDetector(dataBase,train_thresholds_all,train_thresh,thresh_t_dif)

%% test spike detector on 5 minutes of scored IEDs
thresh_SD_opt = 6; 
thresh_t_dif_opt = 0.1;
test_threshold = struct;

for subj = 1:size(dataBase,2)
    
    [sens,prec,F] = testSpikeDetector(dataBase,subj,thresh_SD_opt,thresh_t_dif_opt);
    
    test_threshold(subj).F = F;
    test_threshold(subj).sens = sens;
    test_threshold(subj).prec = prec;
   
end

for subj = 1:size(dataBase,2)
    
    fprintf('%s: F = %0.2f, sens = %0.2f, prec = %0.2f \n',...
    dataBase(subj).sub_label, test_threshold(subj).F, test_threshold(subj).sens, test_threshold(subj).prec)

end

%% detect spikes

thresh_SD = thresh_SD_opt;
thresh_t_dif = thresh_t_dif_opt;

for subj = 1:size(dataBase,2)
    fs = dataBase(subj).ccep_header.Fs;
    
    Mvalues = dataBase(subj).spikes.norm_M;
    Pharmat = dataBase(subj).spikes.Pharmat_norm;
    data_stimfree = dataBase(subj).data_rerefnoStimArt(dataBase(subj).visIEDs.IED,:);

    dataBase(subj).detIED = findCortSpikes(data_stimfree,fs,Mvalues,Pharmat,thresh_t_dif,thresh_SD);
    
end

disp('All spikes are detected')

%% visualize detected spikes

subj = 1;
IEDch = dataBase(subj).IEDchan; % CHECK: dit klopt wellicht niet...
IED = dataBase(subj).IEDch;
plot_SpikesData(dataBase,subj,IEDch, IED)

