%% script to construct svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

%% set paths

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath('git_rep/SPES_SOZ/detectERs')
addpath(genpath('git_rep/eeglab/'))
addpath('git_rep/fieldtrip/')
ft_defaults

%% set configuration
clear

cfg.dataPath = '/Fridge/CCEP';
% old database: PAT119,  PAT126, PAT130, PAT135
cfg.sub_labels = { 'sub-RESP0607', 'sub-RESP0638','sub-RESP0676','sub-RESP0690'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031211','run-031049','run-021423','run-041139'}; 
cfg.train = 1:3;
cfg.test = 4;

% Directory ERSP TFSPES figures
cfg.dir_ERSP = '/Fridge/users/dorien/derivatives/TFSPES/';

% Directory visual ratings
cfg.dir_visrate = '/Fridge/users/dorien/derivatives/BB_article/BB_visrating';


%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg.sub_labels,2)
    
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    
    dataBase(subj).ERSP = load(fullfile(cfg.dir_ERSP, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
end

disp('All ERSPs are loaded')


%% load visual ratings

dataBase = load_visual_BBs(dataBase,cfg);
disp('All visual scores are loaded')

%% determine true labels for fitcsvm: y

[Y_all, Y_conc,Y] = determine_Target_SVM(dataBase,cfg);

%% determine range of ERSPvalues

for subj = 1:3
    nBB = [];
    BB = [];
    
    allERSP = dataBase(subj).ERSP.allERSPboot;
    max_ERSP_BB = NaN(size(allERSP));
    max_ERSP_nBB = NaN(size(allERSP));
    for stimpair=1:size(allERSP,1)                            % for each stimulation pair
        for chan=1:size(allERSP,2)                        % for each recording electrode
            
            % hysteresis3d assumes non-negative image. We want to detect "blue
            % blobs" (negative values). Therefor, I multiply the ERSP with -1
            % and remove all values <0.
            ERSP = allERSP{stimpair,chan};
            ERSP2 = -1* ERSP;
            ERSP2(ERSP2<0) = 0;
            
            if Y{subj}(stimpair,chan) == 0 && ~ismember(chan,dataBase(subj).ERSP.cc_stimsets(stimpair,:)) % no visual scored BB
                nBB = [nBB; ERSP2(:)];          
                max_ERSP_nBB(stimpair,chan) = max(ERSP2(:));
                max_ERSP_BB(stimpair,chan) = NaN;

            elseif Y{subj}(stimpair,chan) == 1 && ~ismember(chan,dataBase(subj).ERSP.cc_stimsets(stimpair,:)) % visual scored BB 
                BB = [BB; ERSP2(:)];
                max_ERSP_BB(stimpair,chan) = max(ERSP2(:));
                max_ERSP_nBB(stimpair,chan) = NaN;
            end
        end
    end
    
    distribution(subj).nBB = nBB;
    distribution(subj).BB = BB;
    distribution(subj).max_ERSP_BB = max_ERSP_BB(:);
    distribution(subj).max_ERSP_nBB = max_ERSP_nBB(:);
end

%% make violin plots


%% step 1: train thresholds hysteresis 
% we vary the thresholds ThU and ThL to find the best thresholds for
% detecting BBs in each patient. We then make a heatmap to find a threshold
% that suits all three patients. 

ThU = 5:0.5:8;  % range upper threshold hysteresis
ThL = 2:0.5:4;  % range lower threshold hysteresis

for i = 1:size(cfg.train,2)
    
    subj = cfg.train(i);
    Y_indiv = Y_conc{subj};

    % all images of one patient in train set
    N = numel(dataBase(subj).ERSP.allERSPboot);
    
    % randomly assign which images are left out in each fold (for cross validation)
    kFolds = 10;     % number of folds
    kIdx = crossvalind('Kfold',N, kFolds); % N is the total number of figures
    
    L_min = NaN(length(ThU),length(ThL));
    
    for p = 1:length(ThU)
        for q = 1:length(ThL)
            
            [D_all,A_all] = getfeaturesTrain(dataBase(subj),ThL(q), ThU(p));
            
            % X-values in Support vector machine
            X = [];
            X(:,1:2) = A_all;                    % store area and duration in X
            X(:,3:4) = D_all;
            % pre-allocation
            gridC = 2.^(-4:2:8);
            L_med = NaN(size(gridC));
            
            for m = 1:size(gridC,2)
                
                C = gridC(m);
                
                L = NaN(kFolds,1);
                for k = 1:kFolds
                    trainData = X(kIdx~=k, :);      % data of 1 of the k folds
                    trainTarg = Y_indiv(kIdx~=k);
                    valData = X(kIdx==k,:);
                    valTarg = Y_indiv(kIdx==k);
                    
                    anSVMModel = fitcsvm(trainData, trainTarg, ...
                        'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
                        'Prior','empirical','Cost',[0 1;4 0],'BoxConstraint', C);
                    
                    L(k) = loss(anSVMModel,valData, valTarg);
                    
                end
                
                % median value for all leave-one-out options
                L_med(m) = median(L);
                
            end
            
            % minimal error with optimizing SVM for this specific set of
            % ThL and ThU
            L_min(p,q) = min(L_med);
            
        end
    end
    
    train_threshold(i).subj = dataBase(subj).sub_label;
    train_threshold(i).L = L_min;

end

%% make heatmaps

minval = min(min([train_threshold(:).L]));
maxval = max(max([train_threshold(:).L]));
cmap = colormap(autumn);
cmap = flipud(cmap);

figure(1), 
for subj=1:3
    subplot(1,3,subj)
    a= train_threshold(subj).L;
%     a = flipud(a);
    imagesc(a,[minval, maxval])
    colormap(cmap)
    xticklabels(ThL)
    yticklabels(ThU)
    title(dataBase(subj).sub_label)
    xlabel('ThL threshold')
    ylabel('ThU threshold')
end

%% 




%% step 2: fit SVM with optimal thresholds for hysteresis
% we now train an SVMmodel with optimal thresholds as determined in step 1
% TODO: nog aanpassen!!

% optimal value upper en lower threshold (determined in step1)
ThU_opt = []; 
ThL_opt = []; 

% all images of patients in train set
n = NaN(size(cfg.train,2),1);
for i=1:size(cfg.train,2)
    n(i) = numel([dataBase(i).ERSP.allERSPboot]);
end
N = sum(n);

% randomly assign which images are left out in each fold (for cross validation)
kFolds = 10;     % number of folds
kIdx = crossvalind('Kfold',N, kFolds); % N is the total number of figures

% get features
D = cell(size(cfg.train,2),1); A = cell(size(cfg.train,2),1);
for i=1:size(cfg.train,2)
    
    [D{i},A{i}] = getfeaturesTrain(dataBase(cfg.train(i)),ThL_opt, ThU_opt);
    
end
        
% vector [3*allERSPimages x 2]
D_all = vertcat(D{:});
A_all = vertcat(A{:});

% X-values in Support vector machine
X(:,1:2) = A_all;                    % store area and duration in X
X(:,3:4) = D_all;
       
gridC = 2.^(-4:2:8);
for C = gridC
    for k = 1:kFolds
        
        % divide all data in train and verification data
        trainData = X(kIdx~=k, :);      % data of 1 of the k folds
        trainTarg = Y(kIdx~=k);
        valData = X(kIdx==k,:);
        valTarg = Y(kIdx==k);
        
        anSVMModel = fitcsvm(trainData, trainTarg, ...
            'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
            'Prior','empirical','Cost',[0 1;4 0],'BoxConstraint', C);
        
        L(k) = loss(anSVMModel,valData, valTarg);
        
    end
    
    L_med = median(L); %--> this should help in determining the optimal SVM and C-value
    
end
  
%% step 3: test SVM on test patient


%% train SVM model with ERSPs from three patients: RESP0607, RESP0638, RESP0676

ThU = 5:0.5:8;  % range upper threshold hysteresis
ThL = 2:0.5:4;  % range lower threshold hysteresis

kFolds = 10;     % number of folds

% all images of patients in train set
n = NaN(size(cfg.train,2),1);
for i=1:size(cfg.train,2)
    n(i) = numel([dataBase(i).ERSP.allERSPboot]);
end
N = sum(n);

bestmodel.SVM = struct('SVMModel', NaN, ...     % this is to store the best SVM
    'C', NaN, 'ThL', NaN, 'ThU', NaN, 'Score', Inf);

% randomly assign which images are left out in each fold (for cross validation)
kIdx = crossvalind('Kfold',N, kFolds); % N is the total number of figures
tic

for k = 1:kFolds
    % forward feature selection starts
    
    bestmodel.ThUScore = inf;
    bestmodel.ThUCombo = struct('SVM', NaN,'ThU', NaN, 'ThL', NaN, 'C', NaN);
    for p = 1:length(ThU)
        bestmodel.ThLScore = inf;
        bestmodel.ThLCombo = struct('SVM', NaN, 'ThL', NaN, 'C', NaN);
        for q = 1:length(ThL)
            clear D A D_all A_all trainData trainTarg testData testTarg
            
            D = cell(size(cfg.train,2),1); A = cell(size(cfg.train,2),1);
            for i=1:size(cfg.train,2)
                
                [D{i},A{i}] = getfeaturesTrain(dataBase(cfg.train(i)),ThL(q), ThU(p));
                
            end
            
            % vector [3*allERSPimages x 2]
            D_all = vertcat(D{:});
            A_all = vertcat(A{:});
                        
            %%Support vector machine
            X(:,1:2) = A_all;                    % store area and duration in X
            X(:,3:4) = D_all;
            
            trainData = X(kIdx~=k, :);      % data of 1 of the k folds
            trainTarg = Y(kIdx~=k);
            valData = X(kIdx==k, :);
            valTarg = Y(kIdx==k);
            
            % this is the grid search for the BoxConstraint
            bestmodel.CScore = inf;
            bestmodel.C = NaN;
            bestmodel.CSVM = NaN;
            gridC = 2.^(-4:2:8);
            for C = gridC
                anSVMModel = fitcsvm(trainData, trainTarg, ...
                    'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
                    'Prior','empirical','Cost',[0 1;4 0],'BoxConstraint', C);
                L = loss(anSVMModel,valData, valTarg);
                if L < bestmodel.CScore        % saving best SVM on parameter
                    bestmodel.CScore = L;      % selection
                    bestmodel.C = C;
                    bestmodel.CSVM = anSVMModel;
                end
            end
            
            % saving the best SVM on lower threshold
            if (bestmodel.CScore <= bestmodel.ThLScore)
                bestmodel.ThLScore = bestmodel.CScore;
                bestmodel.ThLCombo.SVM = bestmodel.CSVM;
                bestmodel.ThLCombo.ThL = ThL(q);
                bestmodel.ThLCombo.C = bestmodel.C;
            end
            
        end
        
        % saving the best SVM on upper threshold
        if (bestmodel.ThLScore <= bestmodel.ThUScore)
            bestmodel.ThUScore = bestmodel.ThLScore;
            bestmodel.ThUCombo.SVM = bestmodel.ThLCombo.SVM;
            bestmodel.ThUCombo.ThL = ThL(q);
            bestmodel.ThUCombo.ThU = ThU(p);
            bestmodel.ThUCombo.C = bestmodel.C;
        end
        
    end
    
    % saving the best SVM over all folds
    if bestmodel.ThUScore <= bestmodel.SVM.Score
        bestmodel.SVM.SVMModel = bestmodel.ThUCombo.SVM;
        bestmodel.SVM.C = bestmodel.ThUCombo.C;
        bestmodel.SVM.ThL = bestmodel.ThUCombo.ThL;
        bestmodel.SVM.ThU = bestmodel.ThUCombo.ThU;
        bestmodel.SVM.Score = bestmodel.ThUScore;
    end
end

toc

%% test SVM model with one patient: RESP0690

detectPowSup(dataBase,cfg)


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
