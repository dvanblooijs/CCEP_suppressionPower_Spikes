%% script to construct svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

addpath(genpath('git_rep/CCEP_suppressionPower_Spikes'))
addpath('git_rep/SPES_SOZ/detectERs')
addpath(genpath('git_rep/eeglab/'))     
addpath('git_rep/fieldtrip/')
ft_defaults

cfg.dataPath = '/Fridge/CCEP';
% old database: PAT119, PAT130, PAT126
cfg.sub_labels = { 'sub-RESP0607', 'sub-RESP0638','sub-RESP0676','sub-RESP0690'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';


%% load visual ratings
cfg.dir_visrate = '/Fridge/users/dorien/derivatives/BB_article/BB_visrating';

D = dir(cfg.dir_visrate);
filenames = {D.name};
dataBase = struct;

for subj = 1 : size(cfg.sub_labels,2)
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    
    clear BS_visscores stimorder_visscores
    
    files = find(contains(filenames,cfg.sub_labels{subj}));
    for i=1:size(files,2)
        % load visual scoring from excel-file
        [num,txt,raw] = xlsread(fullfile(cfg.dir_visrate,filenames{files(i)}));
        % remove first row with stimulation pairs (see txt)
        raw(1,:) = [];
        % remove first column with numbers of electrodes
        raw(:,1) = [];
        % remove NaN-values
        raw(all(cellfun(@(x) any(isnan(x)),raw),2),:) = [];
        raw = celltomat(raw);
        BS_visscores(i,:,:) = raw;
        stimorder_visscores(i,:,:) = txt;
    end
    
    if isequal(stimorder_visscores(1,:), stimorder_visscores(2,:))
    
    dataBase(subj).stimorder_visscores =  stimorder_visscores(1,:);
    else
        disp('ERROR: order of stimulus pairs in visual scores differs')
    end
    
    dataBase(subj).BS_visscores = BS_visscores;
    
end

disp('All visual scores are loaded')

%% Load ERSP-data

cfg.dir_ERSP = '/Fridge/users/dorien/derivatives/TFSPES/';

for subj = 1:size(cfg.sub_labels,2)
    
    D = dir(fullfile(cfg.dir_ERSP, dataBase(subj).sub_label,...
        dataBase(subj).ses_label))
    
    dataBase(subj).run_label  = D(contains({D(:).name},'run')).name;
    
    dataBase(subj).ERSP = load(fullfile(cfg.dir_ERSP, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
end

% load Output_PAT_119
% load times
% load BS_score_PAT119_Dorien.mat
% load BS_score_PAT119_Michelle.mat
% 
% % Load data from mstimefspes_test PAT130
% load Output_PAT_130
% load BS_score_PAT130_Dorien.mat
% load BS_score_PAT130_Michelle.mat
% 
% % Load data from mstimefspes_test PAT126
% load Output_PAT_126
% load BS_score_PAT126_Dorien.mat
% load BS_score_PAT126_Michelle.mat
% allERSP2_126(1:22,:)=[];                    % file contains a SPES session that was not finished, so remove
% stimpchan_126(1:22,:)=[];

%%
ThU = 5:0.5:8;  % range upper threshold hysteresis
ThL = 2:0.5:4;  % range lower threshold hysteresis

kFolds = 10;     % number of folds
N = numel([dataBase(5).ERSP(1).allERSPboot]); % CHECK of het werkt!!!

bestSVM = struct('SVMModel', NaN, ...     % this is to store the best SVM
    'C', NaN, 'ThL', NaN, 'ThU', NaN, 'Score', Inf);

kIdx = crossvalind('Kfold',N, kFolds); % N is the total number of figures
tic

for k = 1:kFolds
    % forward feature selection starts   
    
    bestThUScore = inf;
    bestThUCombo = struct('SVM', NaN,'ThU', NaN, 'ThL', NaN, 'C', NaN);
    for p = 1:length(ThU)
        bestThLScore = inf;
        bestThLCombo = struct('SVM', NaN, 'ThL', NaN, 'C', NaN);
        for q = 1:length(ThL)
            clear S1 D1 A1 S2 D2 A2 S3 D3 A3 A D S trainData trainTarg testData testTarg
            [S1,D1,A1] = getfeaturesTrain(times, allERSP2_119, ThL(q), ThU(p), BS_score_PAT119_Michelle, BS_score_PAT119_Dorien, stimpchan_119);
            [S2,D2,A2] = getfeaturesTrain(times, allERSP2_130, ThL(q), ThU(p), BS_score_PAT130_Michelle, BS_score_PAT130_Dorien, stimpchan_130);
            [S3,D3,A3] = getfeaturesTrain(times, allERSP2_126, ThL(q), ThU(p), BS_score_PAT126_Michelle, BS_score_PAT126_Dorien, stimpchan_126);
            
            A=[A1 A2 A3];
            D=[D1 D2 D3];
            S=[S1 S2 S3];
            
%feature scaling; mean normalization % dit hoeft niet want zit al in
%fitcsvm
%             A = (A-mean(A))/std(A);           
%             D = (D-mean(D))/std(D);
            
            %%Support vector machine
            X(:,1) = A';                    % store area and duration in X
            X(:,2) = D';
            Y=S';                           % store true label in Y
            
            trainData = X(kIdx~=k, :);      % data of 1 of the k folds
            trainTarg = Y(kIdx~=k);
            testData = X(kIdx==k, :);
            testTarg = Y(kIdx==k);
            
            % this is the grid search for the BoxConstraint
            bestCScore = inf;
            bestC = NaN;
            bestCSVM = NaN;
            gridC = 2.^(-4:2:8);
            for C = gridC
                anSVMModel = fitcsvm(trainData, trainTarg, ...
                    'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
                    'Prior','empirical','Cost',[0 1;4 0],'BoxConstraint', C);
                L = loss(anSVMModel,testData, testTarg);
                if L < bestCScore        % saving best SVM on parameter
                    bestCScore = L;      % selection
                    bestC = C;
                    bestCSVM = anSVMModel;
                end
            end
            
            % saving the best SVM on lower threshold
            if (bestCScore < bestThLScore) || ...
                    (bestCScore == bestThLScore)
                bestThLScore = bestCScore;
                bestThLCombo.SVM = bestCSVM;
                bestThLCombo.ThL = ThL(q);
                bestThLCombo.C = bestC;
            end
            
            
        end
        
        % saving the best SVM on upper threshold
        if (bestThLScore < bestThUScore) || ...
                (bestThLScore == bestThUScore)
            bestThUScore = bestThLScore;
            bestThUCombo.SVM = bestThLCombo.SVM;
            bestThUCombo.ThL = ThL(q);
            bestThUCombo.ThU = ThU(p);
            bestThUCombo.C = bestC;
        end
           
    end
    
    % saving the best SVM over all folds
    if bestThUScore < bestSVM.Score
        bestSVM.SVMModel = bestThUCombo.SVM;
        bestSVM.C = bestThUCombo.C;
        bestSVM.ThL = bestThUCombo.ThL;
        bestSVM.ThU = bestThUCombo.ThU;
        bestSVM.Score = bestThUScore
    end
end

toc