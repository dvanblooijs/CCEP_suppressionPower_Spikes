%% script to use svm (constructed with ccepSp01_trainTestSpDetection.m) 
% to detect and visually power suppression in ERSP plots

%% settings for using SVM
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

%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg,2)
    
    dataBase(subj).sub_label = cfg(subj).sub_labels;
    dataBase(subj).ses_label = cfg(subj).ses_label;
    dataBase(subj).task_label = cfg(subj).task_label;
    dataBase(subj).run_label = cfg(subj).run_label{1};
    
    dataBase(subj).ERSP = load(fullfile(myDataPath.ERSPoutput, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
    
   % load channels to exclude bad channels in a later stage
    channelsName = fullfile(myDataPath.dataPath,dataBase(subj).sub_label,dataBase(subj).ses_label,'ieeg',...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_channels.tsv']);
    tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');

    % ECoG electrodes
    idx_ch_incl = strcmp(tb_channels.type,'ECOG') ;
    tb_channels = tb_channels(idx_ch_incl,:);
    dataBase(subj).tb_channels = tb_channels;
    ch_ecog =  tb_channels.name;
    
    if ~isequal(ch_ecog,dataBase(subj).ERSP.ch)
        error('%s has a mismatch in ECoG electrodes and electrodes in ERSP loaded',dataBase(subj).sub_label)
    end
    
    % bad (noisy) ECoG electrodes
    idx_ch_bad = strcmp(tb_channels.status,'bad');    
    dataBase(subj).idx_ch_bad = idx_ch_bad; 
    
end

disp('All ERSPs are loaded')

%% load SVMmodel

pathname = myDataPath.SVMpath;

files = dir(pathname);
idx_SVM = contains({files(:).name},'SVM');

if sum(idx_SVM) >1
   error('More SVMmodels found, make sure the correct one is selected!')
end

SVM = load(fullfile(files(idx_SVM).folder, files(idx_SVM).name));
SVM.SVMfilename = extractBefore(files(idx_SVM).name,'.mat');

fprintf('Loaded %s \n',SVM.SVMfilename)

%% detect power suppression

for subj = 1:size(dataBase,2)
    
    [detPar,dataBase(subj).ERSPdet.statsPost] = getfeaturesTrain(dataBase(subj));

    % vector [3*allERSPimages x 2]
    tStart_all = vertcat(detPar.tStart_conc);
    tWidth_all = vertcat(detPar.tWidth_conc);
    fStart_all = vertcat(detPar.fStart_conc);
    fWidth_all = vertcat(detPar.fWidth_conc);
    Area_all = vertcat(detPar.Area_conc);
    
    % X-values in Support vector machine
    X = [];
    X(:,1:2) = Area_all;    % store area in X
    X(:,3:4) = tStart_all;  % time of start of suppression
    X(:,5:6) = fStart_all;  % minimal frequency of suppression
    X(:,7:8) = tWidth_all;  % duration of suppression
    X(:,9:10) = fWidth_all; % frequency range of suppression
    
    label_raw = predict(SVM.SVMModel,X);
    
    idxnan = isnan(X);
    idxnanrows = sum(idxnan,2);
    
    % put all detected values which were in stimulus pairs (NaNs) to 0
    label_raw(idxnanrows>0) = 0;
    label = reshape(label_raw,size(dataBase(subj).ERSP.allERSPboot));
    
    % put all detected values in bad channels to 0
    label(:,dataBase(subj).idx_ch_bad) = 0;
    
    dataBase(subj).ERSPdet.sub_label = dataBase(subj).sub_label;
    dataBase(subj).ERSPdet.ses_label = dataBase(subj).ses_label;
    dataBase(subj).ERSPdet.task_label = dataBase(subj).task_label;
    dataBase(subj).ERSPdet.run_label = dataBase(subj).run_label;
    dataBase(subj).ERSPdet.cc_stimchans = dataBase(subj).ERSP.cc_stimchans;
    dataBase(subj).ERSPdet.cc_stimsets = dataBase(subj).ERSP.cc_stimsets;
    dataBase(subj).ERSPdet.ch = dataBase(subj).ERSP.ch;

    dataBase(subj).ERSPdet.detected = label;
    dataBase(subj).ERSPdet.Area{1} = reshape(Area_all(:,1),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.Area{2} = reshape(Area_all(:,2),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.tStart{1} = reshape(tStart_all(:,1),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.tStart{2} = reshape(tStart_all(:,2),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.tWidth{1} = reshape(tWidth_all(:,1),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.tWidth{2} = reshape(tWidth_all(:,2),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.fStart{1} = reshape(fStart_all(:,1),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.fStart{2} = reshape(fStart_all(:,2),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.fWidth{1} = reshape(fWidth_all(:,1),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.fWidth{2} = reshape(fWidth_all(:,2),size(dataBase(subj).ERSP.allERSPboot));
    dataBase(subj).ERSPdet.SVM = SVM;
    
    fprintf('Detection of ERSPs in %s is completed\n',dataBase(subj).sub_label)
        
end

disp('Detection is completed.')


%% visually check ERSP with detected power suppression
close all

% select subject
subs = {dataBase(:).sub_label};
string = [repmat('%s, ',1,size(subs,2)-1), '%s'];
substring = input(sprintf(['Choose subject: ',string,'\n'],subs{:}),'s');
subj = find(contains({dataBase(:).sub_label},substring));

if isempty(subj)
   error('No present subject was selected') 
end

if exist(fullfile(myDataPath.ERSPoutput, dataBase(subj).sub_label, dataBase(subj).ses_label, dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPdet.mat']),'file')
    dataBase(subj).ERSPdet = load(fullfile(myDataPath.ERSPoutput, dataBase(subj).sub_label,dataBase(subj).ses_label, dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPdet.mat']));

    fprintf('Loading %s  \n',...
        fullfile(myDataPath.ERSPoutput, dataBase(subj).sub_label,dataBase(subj).ses_label, dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPdet.mat']))
    
end

% continue with the stimulation pair after the last saved stimulation pair
if any(contains(fieldnames(dataBase(subj).ERSPdet),'checkUntilStimp'))
    endstimp = dataBase(subj).ERSPdet.checkUntilStimp;
else
    endstimp = 0;
end

lookPic = dataBase(subj).ERSPdet.detected;

dataBase = visualRating_tfspes(dataBase,subj,lookPic,myDataPath,endstimp);

%% save detected and checked to the ERSP output (with actual ERSPs)
checked = dataBase(subj).ERSPdet.checked;
detected = dataBase(subj).ERSPdet.detected;
SVM = dataBase(subj).ERSPdet.SVM;

output = fullfile(myDataPath.ERSPoutput,dataBase(subj).sub_label,dataBase(subj).ses_label,...
    dataBase(subj).run_label);

filename = ['/' dataBase(subj).sub_label,'_' dataBase(subj).ses_label,...
    '_', dataBase(subj).task_label,'_',dataBase(subj).run_label '_ERSP.mat'];

save(fullfile(output,filename),'checked','detected','SVM','-append')
fprintf('%s is saved \n',[output,filename])

