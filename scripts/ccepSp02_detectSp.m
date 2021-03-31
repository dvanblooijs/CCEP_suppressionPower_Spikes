%% script to use svm (constructed with makeSVM_TFSPES.m) to detect and visually power suppression in TFSPES plots

%% settings for using SVM
clc
clear
myDataPath = setLocalDataPath(1);

%% patient settings
% old database: PAT54, PAT78, PAT88, PAT97, PAT99, PAT114, PAT115, PAT120, PAT123, PAT137
cfg.sub_labels = { 'sub-RESP0401', 'sub-RESP0435', 'sub-RESP0458', 'sub-RESP0478', 'sub-RESP0502',...
    'sub-RESP0574', 'sub-RESP0589', 'sub-RESP0608', 'sub-RESP0621', 'sub-RESP0699'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031153','run-051138','run-011714','run-021549','run-031740'...
    'run-021358','run-021050','run-021057','run-021147','run-031717'};

%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg.sub_labels,2)
    
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    dataBase(subj).ERSP = load(fullfile(myDataPath.TFSPESoutput, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
    
   % load channels to exclude bad channels in a later stage
    channelsName = fullfile(myDataPath.dataPath,cfg.sub_labels{subj},cfg.ses_label,'ieeg',...
        [cfg.sub_labels{subj}, '_', cfg.ses_label, '_', cfg.task_label, '_', cfg.run_label{subj},'_channels.tsv']);
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

load(fullfile(files(idx_SVM).folder, files(idx_SVM).name));

%% detect power suppression

for subj = 1:size(dataBase,2)
    
    [detPar,dataBase(subj).ERSP.statsPost] = getfeaturesTrain(dataBase(subj));

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
    
    label_raw = predict(SVMModel,X);
    label_raw = str2double(label_raw);
    
    idxnan = isnan(X);
    idxnanrows = sum(idxnan,2);
    
    % put all detected values which were in stimulus pairs (NaNs) to 0
    label_raw(idxnanrows>0) = 0;

    dataBase(subj).ERSPdet.sub_label = dataBase(subj).sub_label;
    dataBase(subj).ERSPdet.ses_label = dataBase(subj).ses_label;
    dataBase(subj).ERSPdet.task_label = dataBase(subj).task_label;
    dataBase(subj).ERSPdet.run_label = dataBase(subj).run_label;
    dataBase(subj).ERSPdet.detected = reshape(label_raw,size(dataBase(subj).ERSP.allERSPboot));
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
    
    fprintf('Detection of ERSPs in %s is completed\n',dataBase(subj).sub_label)
        
end

disp('Detection is completed.')


%% visually check TFSPES with detected power suppression
close all

% select subject
subs = {dataBase(:).sub_label};
string = [repmat('%s, ',1,size(subs,2)-1), '%s'];
substring = input(sprintf(['Choose subject: ',string,'\n'],subs{:}),'s');
subj = find(contains({dataBase(:).sub_label},substring));

if isempty(subj)
   error('No present subject was selected') 
end

if exist(fullfile(myDataPath.dir_visrate, dataBase(subj).sub_label, dataBase(subj).ses_label, dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPsChecked.mat']),'file')
   dataBase(subj).ERSP = load(fullfile(myDataPath.dir_visrate, dataBase(subj).sub_label,dataBase(subj).ses_label, dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPsChecked.mat']));   
end

% continue with the stimulation pair after the last saved stimulation pair
if any(contains(fieldnames(dataBase(subj).ERSPdet),'checkUntilStimp'))
    endstimp = dataBase(subj).ERSPdet.checkUntilStimp;
else
    endstimp = 0;
end

lookPic = dataBase(subj).ERSPdet.detected;

dataBase = visualRating_tfspes(dataBase,subj,lookPic,myDataPath,endstimp);

% save detected and checked to the TFSPES output
checked = dataBase(subj).ERSPdet.checked;
detected = dataBase(subj).ERSPdet.detected;

output = fullfile(myDataPath.TFSPESoutput,dataBase(subj).sub_label,dataBase(subj).ses_label,...
    dataBase(subj).run_label);

filename = ['/' dataBase(subj).sub_label,'_' dataBase(subj).ses_label,...
    '_', dataBase(subj).task_label,'_',dataBase(subj).run_label '_ERSP.mat'];

save(fullfile(output,filename),'checked','detected','-append')
fprintf('%s is saved \n',[output,filename])

