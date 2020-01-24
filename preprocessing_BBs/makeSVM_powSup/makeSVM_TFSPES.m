%% script to construct svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019


%% determine true labels for fitcsvm: y

[Y_alltrain,Y_alltest, Y_conc, Y] = determine_Target_SVM(dataBase,cfg);

%% determine range of ERSPvalues

% distribution_ERSPval(cfg,dataBase)

%% step 1a: train thresholds hysteresis
% we vary the thresholds ThU and ThL to find the best thresholds for
% detecting BBs in each patient. We then make a heatmap to find a threshold
% that suits all three patients.

trainPar.ThU = 0.1:0.1:0.9;  % range upper threshold hysteresis
trainPar.ThL = 0.1:0.1:0.9;  % range lower threshold hysteresis
trainPar.C = 2.^(-4:2:8); % range of BoxConstraint in svm model
trainPar.boot = input('Use bootstrapped ERSP [yes] or ERSP without bootstrapping [no]:', 's');

[train_threshold,train_threshold_all] = trainSVM_PowSup(cfg,dataBase,trainPar,Y_conc);

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
ThL_opt = 0.2;
C_opt = 0.25;

% all images of patients in train set
n = NaN(size(cfg.train,2),1);
for i=1:size(cfg.train,2)
    n(i) = numel([dataBase(i).ERSP.allERSPboot]);
end
N = sum(n);

% get features
D = cell(size(cfg.train,2),1); A = cell(size(cfg.train,2),1);
for i=1:size(cfg.train,2)
    
    [D{i},A{i}] = getfeaturesTrain(dataBase(cfg.train(i)),ThL_opt, ThU_opt);
    
end

% vector [3*allERSPimages x 2]
D_all = vertcat(D{:});
A_all = vertcat(A{:});

X = [];
% X-values in Support vector machine
X(:,1:2) = A_all;                    % store area and duration in X
X(:,3:4) = D_all;

% optimal SVM model
SVMModel = fitcsvm(X, Y_alltrain, ...
    'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
    'Prior','empirical','Cost',[0 1;4 0],'BoxConstraint', C_opt,'CrossVal','on');


% save SVMmodel
pathname = '/Fridge/users/dorien/derivatives/BB_article/';
filename = sprintf('SVMmodel_trained_BB_%1.1f_%1.1f_%1.2f_%s.mat',ThU_opt,ThL_opt,C_opt,datestr(now,'yyyymmdd'));
save(fullfile(pathname,filename),'SVMModel','trainPar')

%% test SVM model with one patient: RESP0690

% detectPowSup(dataBase,cfg)

[D_all,A_all] = getfeaturesTrain(dataBase(cfg.test),ThL_opt, ThU_opt);

% X-values in Support vector machine
X = [];
X(:,1:2) = A_all;                    % store area and duration in X
X(:,3:4) = D_all;

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
    sens(k) = TP/(TP+FN);
    spec(k) = TN/(TN+FP);
    prec(k) = TP/(TP+FP);
    F(k) = 2 * (prec(k) * sens(k))/(prec(k) + sens(k));
end



%% check detected ERSPs

subj = 1;
stimps = 1:size(dataBase(subj).cc_epoch_sorted,3);
chans = dataBase(subj).soz;%1:size(dataBase(1).cc_epoch_sorted,1);

dataBase = visualRating_tfspes(dataBase,stimps, chans);
checked = dataBase(1).ERSP.checked;

output = fullfile(cfg.TFSPESpath,dataBase(subj).sub_label,dataBase(subj).ses_label,...
    ['task-',dataBase(subj).task_label,'_',dataBase(subj).run_label]);

filename = ['/sub-' dataBase(subj).sub_label,'_' dataBase(subj).ses_label,...
    '_task-', dataBase(subj).task_label,'_',dataBase(subj).run_label '_ERSP.mat'];

save([output,filename],'checked','-append')
