%% script to construct svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

%% settings for constructing SVM
clc
clear
cfg = setLocalDataPath(1);

%%

% old database: PAT119,  PAT126, PAT130, PAT135
cfg.sub_labels = { 'sub-RESP0607', 'sub-RESP0638','sub-RESP0676','sub-RESP0690'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031211','run-031049','run-021423','run-041139'};
cfg.train = 1:3;
cfg.test = 4;

%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg.sub_labels,2)
    
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    dataBase(subj).ERSP = load(fullfile(cfg.TFSPESoutput, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
end

disp('All ERSPs are loaded')

%% load visual ratings

% dataBase = load_visual_BBs(dataBase,cfg);
% disp('All visual scores are loaded')

for subj =1:size(dataBase,2)
    folderName = fullfile(cfg.dir_visrate,[cfg.sub_labels{subj},'_',cfg.ses_label,'_',cfg.run_label{subj},'_visBB_combined.mat']);
    
    dataBase(subj).visBB = load(folderName);

end

%% determine true labels for fitcsvm: y 

[Y_alltrain,Y_alltest, Y_conc, Y] = determine_Target_SVM(dataBase,cfg);

%% step 1: train thresholds hysteresis : this takes a few hours
% we vary the thresholds ThU and ThL to find the best thresholds for
% detecting BBs in each patient. We then make a heatmap to find a threshold
% that suits all three patients.

trainPar.C = 2.^(-4:2:8); % range of BoxConstraint in svm model

[train_threshold,train_threshold_all] = trainSVM_PowSup(cfg,dataBase,trainPar,Y_conc);

disp('Training completed')

%% step 1b: make heatmaps for each patient and for patients combined

mode = input('Visualize Loss or sensitivity or specificity or precision or F-value [L/sens/spec/prec/F]: ','s');
par1 = input('Visualize parameter1 (ThU or ThL or C) [ThU/ThL/C]: ','s');
par2 = input('Visualize parameter2 (ThU or ThL or C) [ThU/ThL/C]: ','s');

heatmap_trainSVM(cfg,dataBase,train_threshold,trainPar,mode,par1,par2)
heatmap_trainSVM(cfg,dataBase,train_threshold_all,trainPar,mode,par1,par2)

%% find optimal values for each individual patient and all patients in train dataset combined

clc
trainPar.mode = input('Find optimal for each patient [optimal], or define an option yourself [option]: ','s');

if strcmp(trainPar.mode,'option')
    trainPar.ThU_opt = str2double(input('Value of ThU: ','s'));
    trainPar.ThL_opt = str2double(input('Value of ThL: ','s'));
    trainPar.C_opt = str2double(input('Value of C: ','s'));
    
    optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold,subj)
    optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold_all,1)
    
else
    trainPar.ThU_opt = [];
    trainPar.ThL_opt = [];
    trainPar.C_opt = [];
    
    for n = 1:size(cfg.train,2)
        subj = cfg.train(n);
        optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold,subj)
    end
    
    optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold_all,1)
end

%% step 2: fit SVM with optimal thresholds for hysteresis
% we now train an SVMmodel with optimal thresholds as determined in step 1

% optimal value upper en lower threshold (determined in step1)
ThU_opt = 0.7;
ThL_opt = 0.4;
C_opt = 1;

% all images of patients in train set
n = NaN(size(cfg.train,2),1);
for i=1:size(cfg.train,2)
    n(i) = numel([dataBase(i).ERSP.allERSPboot]);
end
N = sum(n);

% get features
tStart = cell(size(cfg.train,2),1); tWidth = cell(size(cfg.train,2),1);
fStart = cell(size(cfg.train,2),1); fWidth = cell(size(cfg.train,2),1);
Area = cell(size(cfg.train,2),1);

for subj=1:size(cfg.train,2)
    
    [Area{subj}, tStart{subj}, fStart{subj}, tWidth{subj}, fWidth{subj},...
        dataBase(subj).ERSP.Area, dataBase(subj).ERSP.tStart,  ...
        dataBase(subj).ERSP.fStart, dataBase(subj).ERSP.tWidth, ...
        dataBase(subj).ERSP.fWidth] = getfeaturesTrain(dataBase(cfg.train(subj)),ThL_opt, ThU_opt);
    
end

% vector [3*allERSPimages x 2]
tStart_all = vertcat(tStart{:});
tWidth_all = vertcat(tWidth{:});
fStart_all = vertcat(fStart{:});
fWidth_all = vertcat(fWidth{:});
Area_all = vertcat(Area{:});

% X-values in Support vector machine
X = [];
X(:,1:2) = Area_all;                    % store area and duration in X
X(:,3:4) = tStart_all;
X(:,5:6) = fStart_all;
X(:,7:8) = tWidth_all;
X(:,9:10) = fWidth_all;
            
X_train = X;
% optimal SVM model
SVMModel = fitcsvm(X_train, Y_alltrain, ...
    'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
    'Cost',[0 1;4 0],'BoxConstraint', C_opt);
%      'Cost',[0 1;4 0],'BoxConstraint', C_opt);
%     'Prior','empirical','BoxConstraint', C_opt);

resubLoss(SVMModel)

CVModel = crossval(SVMModel);
kfoldLoss(CVModel)

% save SVMmodel
pathname = cfg.SVMpath;
trainPar.ThU_opt = ThU_opt;
trainPar.ThL_opt = ThL_opt;
traomPar.C_opt = C_opt;
filename = sprintf('SVMmodel_trained_BB_%1.1f_%1.1f_%1.2f_%s.mat',ThU_opt,ThL_opt,C_opt,datestr(now,'yyyymmdd'));
save(fullfile(pathname,filename),'SVMModel','trainPar')
fprintf('SVMmodel is saved in %s\n',filename)

%% test SVM model with one patient: RESP0690

% detectPowSup(dataBase,cfg)
subj = cfg.test;

[Area, tStart, fStart, tWidth, fWidth,...
    dataBase(subj).ERSP.Area, dataBase(subj).ERSP.tStart,  ...
        dataBase(subj).ERSP.fStart, dataBase(subj).ERSP.tWidth, ...
        dataBase(subj).ERSP.fWidth] = getfeaturesTrain(dataBase(cfg.test),ThL_opt, ThU_opt);

% X-values in Support vector machine
X = [];
X(:,1:2) = Area;                    % store area and duration in X
X(:,3:4) = tStart;
X(:,5:6) = fStart;
X(:,7:8) = tWidth;
X(:,9:10) = fWidth;

valData = X;
valTarg = Y_alltest;

label_raw = predict(SVMModel,valData);

label = str2num(cell2mat(label_raw)); %#ok<ST2NM>
TP = numel(find(label == 1 & valTarg == 1));
FP = numel(find(label == 1 & valTarg == 0));
TN = numel(find(label == 0 & valTarg == 0));
FN = numel(find(label == 0 & valTarg == 1));
if sum([TP;FP;TN;FN]) ~= size(valTarg,1)
    error('Number of TP, FP, TN, FN does not add up to all labels')
else
    sens = TP/(TP+FN);
    spec = TN/(TN+FP);
    prec = TP/(TP+FP);
    F = 2 * (prec * sens)/(prec + sens);
end

%% check FN or FP

mode = input('What would you like to check? [FP/FN/TP/TN]: ','s');

if strcmp(mode,'FP')
    lookPic = reshape((label==1 & valTarg == 0),size(dataBase(subj).ERSP.allERSP));
elseif strcmp(mode,'FN')
    lookPic = reshape((label==0 & valTarg == 1),size(dataBase(subj).ERSP.allERSP));
elseif strcmp(mode,'TP')
    lookPic = reshape((label==1 & valTarg == 1),size(dataBase(subj).ERSP.allERSP));   
elseif strcmp(mode,'TN')
     lookPic = reshape((label==0 & valTarg == 0),size(dataBase(subj).ERSP.allERSP));   
end

visualRating_tfspes(dataBase,subj,lookPic);

% %% detect ERSPs
% 
% for subj = 1:size(dataBase,2)
%     
%     [Area_conc, tStart_conc, fStart_conc, tWidth_conc, fWidth_conc, ...
%         dataBase(subj).ERSP.Area, dataBase(subj).ERSP.tStart,  ...
%         dataBase(subj).ERSP.fStart, dataBase(subj).ERSP.tWidth, ...
%         dataBase(subj).ERSP.fWidth] = getfeaturesTrain(dataBase(subj),ThL_opt, ThU_opt);
%     
%     % store area and duration in X
%     X = [];
%     X(:,1:2) = Area_conc;                    
%     X(:,3:4) = tStart_conc;
%     X(:,5:6) = fStart_conc;
%     X(:,7:8) = tWidth_conc;
%     X(:,9:10) = fWidth_conc;
%     
%     label_raw = predict(SVMModel,X);
%     label_raw = str2double(label_raw);
%     
%     dataBase(subj).ERSP.detected = reshape(label_raw,size(dataBase(subj).ERSP.allERSPboot));
%    fprintf('Detection of ERSPs in %s is completed\n',dataBase(subj).sub_label) 
%    
%    
% end
% 
% disp('Detection is completed.')
% 
% %% check detected ERSPs - optional
% 
% dataBase = visualRating_tfspes(dataBase,subj);
% checked = dataBase(subj).ERSP.checked;
% 
% output = fullfile(cfg.TFSPESoutput,dataBase(subj).sub_label,dataBase(subj).ses_label,...
%     [dataBase(subj).task_label,'_',dataBase(subj).run_label]);
% 
% filename = ['/' dataBase(subj).sub_label,'_' dataBase(subj).ses_label,...
%     '_', dataBase(subj).task_label,'_',dataBase(subj).run_label '_ERSP.mat'];
% 
% save([output,filename],'checked','-append')
