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
cfg.run_label = {'run-021120','run-031049','run-021423','run-041139'};
cfg.train = 1:3;
cfg.test = 4;

%% load visual ratings
cfg.dir_visrate = '/Fridge/users/dorien/derivatives/BB_article/BB_visrating';

D = dir(cfg.dir_visrate);
filenames = {D.name};
% dataBase = struct;

for subj = 1 : size(cfg.sub_labels,2)
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    clear BS_visscores stimorder_visscores chan_visscores stimp_visscores stimpnum_visscores
    
    files = find(contains(filenames,cfg.sub_labels{subj}));
    if ~isempty(files)
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
            stimorder_visscores(i,:,:) = txt(2:end,1);
            chan_visscores(i,:,:) = txt(1,2:end)';
        end
        
        if isequal(stimorder_visscores(1,:), stimorder_visscores(2,:)) && isequal(chan_visscores(1,:),chan_visscores(2,:))
            
            splitcell = cellfun( @(x) strsplit(x,'-'),stimorder_visscores(1,:),'UniformOutput',false);
            
            for n=1:size(splitcell,2)
                stimp_visscores{n,1} = splitcell{n}{1};
                stimp_visscores{n,2} = splitcell{n}{2};
                
                stimpnum_visscores(n,1) = find(strcmpi(splitcell{n}{1},squeeze(chan_visscores(1,:))));
                stimpnum_visscores(n,2) = find(strcmpi(splitcell{n}{2},squeeze(chan_visscores(1,:))));
            
            end
                
            dataBase(subj).stimorder_visscores = stimorder_visscores(1,:)';
            dataBase(subj).stimpnum_visscores = stimpnum_visscores;
            dataBase(subj).stimp_visscores = stimp_visscores;
            dataBase(subj).chan_visscores = chan_visscores(1,:)';
        else
            disp('ERROR: order of stimulus pairs in visual scores differs')
        end
        
        dataBase(subj).BS_visscores = BS_visscores;
    else
        disp('No visual scores for this patient')
    end
end

disp('All visual scores are loaded')

%% Load ERSP-data

cfg.dir_ERSP = '/Fridge/users/dorien/derivatives/TFSPES/';

for subj = 1:size(cfg.sub_labels,2)
        
    dataBase(subj).ERSP = load(fullfile(cfg.dir_ERSP, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
end

disp('All ERSPs are loaded')

%% train SVM model with ERSPs from three patients: RESP0607, RESP0638, RESP0676

ThU = 5:0.5:8;  % range upper threshold hysteresis
ThL = 2:0.5:4;  % range lower threshold hysteresis

kFolds = 10;     % number of folds

% all images of patients in train set
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
            clear S D A S_all D_all A_all trainData trainTarg testData testTarg
            
            for i=1:size(cfg.train,2)
                        
                [S{i},D{i},A{i}] = getfeaturesTrain(dataBase(cfg.train(i)),ThL(q), ThU(p));
                
            end
            
            % vector [3*allERSPimages x 2]
            S_all = vertcat(S{:});
            D_all = vertcat(D{:});
            A_all = vertcat(A{:});
                        
%feature scaling; mean normalization % dit hoeft niet want zit al in
%fitcsvm
%             A = (A-mean(A))/std(A);           
%             D = (D-mean(D))/std(D);
            
            %%Support vector machine
            X(:,1:2) = A_all;                    % store area and duration in X
            X(:,3:4) = D_all;
            Y=S_all;                           % store true label in Y
            
            trainData = X(kIdx~=k, :);      % data of 1 of the k folds
            trainTarg = Y(kIdx~=k);
            testData = X(kIdx==k, :);
            testTarg = Y(kIdx==k);
            
            % this is the grid search for the BoxConstraint
            bestmodel.CScore = inf;
            bestmodel.C = NaN;
            bestmodel.CSVM = NaN;
            gridC = 2.^(-4:2:8);
            for C = gridC
                anSVMModel = fitcsvm(trainData, trainTarg, ...
                    'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
                    'Prior','empirical','Cost',[0 1;4 0],'BoxConstraint', C);
                L = loss(anSVMModel,testData, testTarg);
                if L < bestmodel.CScore        % saving best SVM on parameter
                    bestmodel.CScore = L;      % selection
                    bestmodel.C = C;
                    bestmodel.CSVM = anSVMModel;
                end
            end
            
            % saving the best SVM on lower threshold
            if (bestmodel.CScore < bestmodel.ThLScore) || ...
                    (bestmodel.CScore == bestmodel.ThLScore)
                bestmodel.ThLScore = bestmodel.CScore;
                bestmodel.ThLCombo.SVM = bestmodel.CSVM;
                bestmodel.ThLCombo.ThL = ThL(q);
                bestmodel.ThLCombo.C = bestmodel.C;
            end
            
            
        end
        
        % saving the best SVM on upper threshold
        if (bestmodel.ThLScore < bestmodel.ThUScore) || ...
                (bestmodel.ThLScore == bestmodel.ThUScore)
            bestmodel.ThUScore = bestmodel.ThLScore;
            bestmodel.ThUCombo.SVM = bestmodel.ThLCombo.SVM;
            bestmodel.ThUCombo.ThL = ThL(q);
            bestmodel.ThUCombo.ThU = ThU(p);
            bestmodel.ThUCombo.C = bestmodel.C;
        end
           
    end
    
    % saving the best SVM over all folds
    if bestmodel.ThUScore < bestmodel.SVM.Score
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
