%% ccepSp02_detectSp
% script to use a trained svm to detect and visually check power 
% suppression in ERSP plots

% 1. set paths
% 2. select subjects
% 3. load ERSP data 
% 4. load SVM model
% 5. detect power suppression
% 6. visually check power supperessions
% 7. save matrix with detected and checked power suppression additionally in ERSP.mat

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

%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for nSubj = 1:size(cfg,2)
    
    dataBase(nSubj).sub_label = cfg(nSubj).sub_labels;
    dataBase(nSubj).ses_label = cfg(nSubj).ses_label;
    dataBase(nSubj).task_label = cfg(nSubj).task_label;
    dataBase(nSubj).run_label = cfg(nSubj).run_label{1};
    
    dataBase(nSubj).ERSP = load(fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label,...
        [dataBase(nSubj).sub_label, '_',dataBase(nSubj).ses_label,'_' dataBase(nSubj).task_label,...
        '_', dataBase(nSubj).run_label,'_ERSP.mat']));
    
   % load channels to exclude bad channels in a later stage
    channelsName = fullfile(myDataPath.proj_dirinput,dataBase(nSubj).sub_label,dataBase(nSubj).ses_label,'ieeg',...
        [dataBase(nSubj).sub_label, '_',dataBase(nSubj).ses_label,'_' dataBase(nSubj).task_label,...
        '_', dataBase(nSubj).run_label,'_channels.tsv']);
    tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');

    % ECoG electrodes
    idx_ch_incl = strcmp(tb_channels.type,'ECOG') ;
    tb_channels = tb_channels(idx_ch_incl,:);
    dataBase(nSubj).tb_channels = tb_channels;
    ch_ecog =  tb_channels.name;
    
    if ~isequal(ch_ecog,dataBase(nSubj).ERSP.ch)
        error('%s has a mismatch in ECoG electrodes and electrodes in ERSP loaded',dataBase(nSubj).sub_label)
    end
    
    % bad (noisy) ECoG electrodes
    idx_ch_bad = strcmp(tb_channels.status,'bad');    
    dataBase(nSubj).idx_ch_bad = idx_ch_bad; 
    
end

disp('All ERSPs are loaded')

% housekeeping
clear ch_ecog channelsName files files_eeg files_ses files_subj idx_ch_bad idx_ch_incl idx_eeg idx_ses idx_subj run runTemp nSubj tb_channels

%% load SVMmodel

pathname = fullfile(myDataPath.proj_diroutput,'SVM');

files = dir(pathname);
idx_SVM = contains({files(:).name},'SVM');

if sum(idx_SVM) >1
   error('More SVMmodels found, make sure the correct one is selected!')
end

SVM = load(fullfile(pathname,'SVMmodel_trained_BB_20210331.mat'));
SVM.SVMfilename = extractBefore(files(idx_SVM).name,'.mat');

fprintf('Loaded %s \n',SVM.SVMfilename)

%% detect power suppression

for nSubj = 1:size(dataBase,2)
    
    [detPar,dataBase(nSubj).ERSPdet.statsPost] = getfeaturesTrain(dataBase(nSubj));

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
    label = reshape(label_raw,size(dataBase(nSubj).ERSP.allERSPboot));
    
    % put all detected values in bad channels to 0
    label(:,dataBase(nSubj).idx_ch_bad) = 0;
    
    dataBase(nSubj).ERSPdet.sub_label = dataBase(nSubj).sub_label;
    dataBase(nSubj).ERSPdet.ses_label = dataBase(nSubj).ses_label;
    dataBase(nSubj).ERSPdet.task_label = dataBase(nSubj).task_label;
    dataBase(nSubj).ERSPdet.run_label = dataBase(nSubj).run_label;
    dataBase(nSubj).ERSPdet.cc_stimchans = dataBase(nSubj).ERSP.cc_stimchans;
    dataBase(nSubj).ERSPdet.cc_stimsets = dataBase(nSubj).ERSP.cc_stimsets;
    dataBase(nSubj).ERSPdet.ch = dataBase(nSubj).ERSP.ch;

    dataBase(nSubj).ERSPdet.detected = label;
    dataBase(nSubj).ERSPdet.Area{1} = reshape(Area_all(:,1),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.Area{2} = reshape(Area_all(:,2),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.tStart{1} = reshape(tStart_all(:,1),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.tStart{2} = reshape(tStart_all(:,2),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.tWidth{1} = reshape(tWidth_all(:,1),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.tWidth{2} = reshape(tWidth_all(:,2),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.fStart{1} = reshape(fStart_all(:,1),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.fStart{2} = reshape(fStart_all(:,2),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.fWidth{1} = reshape(fWidth_all(:,1),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.fWidth{2} = reshape(fWidth_all(:,2),size(dataBase(nSubj).ERSP.allERSPboot));
    dataBase(nSubj).ERSPdet.SVM = SVM;
    
    fprintf('Detection of ERSPs in %s is completed\n',dataBase(nSubj).sub_label)
        
end

disp('Detection is completed.')


%% visually check ERSP with detected power suppression
close all

% select subject
subs = {dataBase(:).sub_label};
string = [repmat('%s, ',1,size(subs,2)-1), '%s'];
substring = input(sprintf(['Choose subject: ',string,'\n'],subs{:}),'s');
nSubj = find(contains({dataBase(:).sub_label},substring));

if isempty(nSubj)
   error('No present subject was selected') 
end

if exist(fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label,...
        [dataBase(nSubj).sub_label, '_', dataBase(nSubj).ses_label,'_',dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label,'_ERSPdet.mat']),'file')
    
    dataBase(nSubj).ERSPdet = load(fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label,...
        [dataBase(nSubj).sub_label, '_', dataBase(nSubj).ses_label,'_',dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label,'_ERSPdet.mat']));

    fprintf('Loading %s  \n',...
        fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label,...
        [dataBase(nSubj).sub_label, '_', dataBase(nSubj).ses_label,'_',dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label,'_ERSPdet.mat']))
    
end

% continue with the stimulation pair after the last saved stimulation pair
if any(contains(fieldnames(dataBase(nSubj).ERSPdet),'checkUntilStimp'))
    endstimp = dataBase(nSubj).ERSPdet.checkUntilStimp;
else
    endstimp = 0;
end

lookPic = dataBase(nSubj).ERSPdet.detected;

dataBase = visualRating_tfspes(dataBase,nSubj,lookPic,myDataPath,endstimp);

%% save detected and checked to the ERSP output (with actual ERSPs)
checked = dataBase(nSubj).ERSPdet.checked;
detected = dataBase(nSubj).ERSPdet.detected;
SVM = dataBase(nSubj).ERSPdet.SVM;

output = fullfile(myDataPath.proj_diroutput,dataBase(nSubj).sub_label);

filename = [dataBase(nSubj).sub_label,'_' dataBase(nSubj).ses_label,...
    '_', dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label '_ERSP.mat'];

save(fullfile(output,filename),'checked','detected','SVM','-append')
fprintf('%s is saved \n',fullfile(output,filename))

