% This function trains a support vector machine for detecting power
% suppression in ERSP-figures. 
% PARAMETERS INPUT:
% - cfg.train                   : which subjects are in the train dataset
% - trainPar.ThU                : range of values for the upperThreshold in hysteresis function (for delineating area of power suppression)
% - trainPar.ThL                : range of values for the lowerThreshold in hysteresis function (for delineating area of power suppression)
% - trainPar.C                  : range of values for the boxconstraint for fitting the SVM
% - dataBase.sub_label          : names of each subject {name}
% - dataBase.ERSP.allERSPboot   : ERSP-arrays [timesxfreqs] for each stimuluspair-electrode combination
% - dataBase.ERSP.cc_stimsets   : containt the stimsets in this dataset [stimpairsx2]

% PARAMETERS OUTPUT:
% - train_threshold.subj        : name of the subject
% - train_threshold.L           : Loss for each combination of values in trainPar
% - train_threshold.sens        : sensitivity for each combination of values in trainPar
% - train_threshold.spec        : specifiticy for each combination of values in trainPar
% - train_threshold.prec        : precision for each combination of values in trainPar
% - train_threshold.F           : F-value (2* ((prec*sens)/(prec+sens))) for each combination of values in trainPar
% - train_threshold_all         : the performance of all subjects in the train dataset combined

function [train_threshold,train_threshold_all] = trainSVM_PowSup(cfg,dataBase,trainPar,Y_conc)

ThU = trainPar.ThU;  % range upper threshold hysteresis
ThL = trainPar.ThL;  % range lower threshold hysteresis
gridC = trainPar.C; % range of BoxConstraint in svm model

% pre-allocation
train_threshold = struct; 
L_comb      = NaN(size(cfg.train,3),size(ThU,2),size(ThL,2),size(gridC,2));
sens_comb   = NaN(size(cfg.train,3),size(ThU,2),size(ThL,2),size(gridC,2));
spec_comb   = NaN(size(cfg.train,3),size(ThU,2),size(ThL,2),size(gridC,2));
prec_comb   = NaN(size(cfg.train,3),size(ThU,2),size(ThL,2),size(gridC,2));
F_comb      = NaN(size(cfg.train,3),size(ThU,2),size(ThL,2),size(gridC,2));

for subj = cfg.train
    
    Y_indiv = Y_conc{subj};
    
    % all images of one patient in train set
    N = numel(dataBase(subj).ERSP.allERSPboot);
    
    % randomly assign which images are left out in each fold (for cross validation)
    kFolds = 10;     % number of folds
    kIdx = crossvalind('Kfold',N, kFolds); % N is the total number of figures
    
    L_all = NaN(length(ThU),length(ThL),length(gridC));
    sens_all = NaN(length(ThU),length(ThL),length(gridC));
    prec_all = NaN(length(ThU),length(ThL),length(gridC));
    spec_all = NaN(length(ThU),length(ThL),length(gridC));
    F_all = NaN(length(ThU),length(ThL),length(gridC));
    
    for p = 1:length(ThU)
        for q = 1:length(ThL)
            
            [Area_conc, tStart_conc, fStart_conc, tWidth_conc, fWidth_conc] = getfeaturesTrain(dataBase(subj),ThL(q), ThU(p));
            
            % X-values in Support vector machine
            X = [];
            X(:,1:2)  = Area_conc;                    % store area, onset time, onset frequency, duration and frequency band in X
            X(:,3:4)  = tStart_conc;
            X(:,5:6)  = fStart_conc;
            X(:,7:8)  = tWidth_conc;
            X(:,9:10) = fWidth_conc;
            
            for m = 1:size(gridC,2)
                
                C = gridC(m);
                
                L = NaN(kFolds,1); sens = NaN(kFolds,1); spec = NaN(kFolds,1); prec = NaN(kFolds,1); F = NaN(kFolds,1);
                for k = 1:kFolds
                    trainData = X(kIdx~=k, :);      % data of 1 of the k folds
                    trainTarg = Y_indiv(kIdx~=k);
                    valData = X(kIdx==k,:);
                    valTarg = Y_indiv(kIdx==k);
                    
                    anSVMModel = fitcsvm(trainData, trainTarg, ...
                        'Standardize',true,'ClassNames',{'0','1'},'KernelFunction', 'RBF', 'KernelScale', 'auto', ...
                        'Prior','empirical','BoxConstraint', C);
                    
                    L(k) = loss(anSVMModel,valData, valTarg);
                    label_raw = predict(anSVMModel,valData);
                    
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
                end
                
                % median value for all leave-one-out options
                L_med = median(L);
                sens_med = median(sens);
                spec_med = median(spec);
                prec_med = median(prec);
                F_med = median(F);
                
                % error with optimizing SVM for this specific set of ThL and ThU
                L_all(p,q,m) = L_med;
                sens_all(p,q,m) = sens_med;
                spec_all(p,q,m) = spec_med;
                prec_all(p,q,m) = prec_med;
                F_all(p,q,m) = F_med;
                
                fprintf('---- subject = %s, ThL = %g, ThU = %g, C = %g ----\n',dataBase(subj).sub_label, ThL(q),ThU(p),C)
                
            end
        end
    end
    
    train_threshold(subj).subj  = dataBase(subj).sub_label;
    train_threshold(subj).L     = L_all;
    train_threshold(subj).sens  = sens_all;
    train_threshold(subj).spec  = spec_all;
    train_threshold(subj).prec  = prec_all;
    train_threshold(subj).F     = F_all;
    
    L_comb(subj,:,:,:)        = L_all;
    sens_comb(subj,:,:,:)     = sens_all;
    spec_comb(subj,:,:,:)     = spec_all;
    prec_comb(subj,:,:,:)     = prec_all;
    F_comb(subj,:,:,:)        = F_all;
end

% combine all patients
train_threshold_all.L       = squeeze(median(L_comb,1));
train_threshold_all.sens    = squeeze(median(sens_comb,1));
train_threshold_all.spec    = squeeze(median(spec_comb,1));
train_threshold_all.prec    = squeeze(median(prec_comb,1));
train_threshold_all.F       = squeeze(median(F_comb,1));

end