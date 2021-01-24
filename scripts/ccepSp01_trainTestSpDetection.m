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
    
    % load ERSPs
    dataBase(subj).ERSP = load(fullfile(cfg.TFSPESoutput, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));

    % load channels to exclude bad channels in a later stage
    channelsName = fullfile(cfg.dataPath,cfg.sub_labels{subj},cfg.ses_label,'ieeg',...
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

%% load visual ratings

% dataBase = load_visual_BBs(dataBase,cfg);
% disp('All visual scores are loaded')

for subj =1:size(dataBase,2)
    folderName = fullfile(cfg.dir_visrate,[cfg.sub_labels{subj},'_',cfg.ses_label,'_',cfg.run_label{subj},'_visBB_combined.mat']);
    
    dataBase(subj).visBB = load(folderName);

end

disp('All visual ratings are loaded')

%% determine true labels for fitcsvm: y 

[Y_alltrain,Y_alltest, Y_conc, Y] = determine_Target_SVM(dataBase,cfg);

%% step 1: train thresholds hysteresis : this takes a few hours
% we vary the thresholds ThU and ThL to find the best thresholds for
% detecting BBs in each patient. We then make a heatmap to find a threshold
% that suits all three patients.

% trainPar.C = 2.^(-4:2:8); % range of BoxConstraint in svm model
% 
% [train_threshold,train_threshold_all] = trainSVM_PowSup(cfg,dataBase,trainPar,Y_conc);
% 
% disp('Training completed')

%% step 1b: make heatmaps for each patient and for patients combined

% plotPerformance_trainSVM(train_threshold,trainPar)

%% find optimal values for each individual patient and all patients in train dataset combined

% trainPar.mode = input('Find optimal for each patient [optimal], or define an option yourself [option]: ','s');
% 
% if strcmp(trainPar.mode,'option')
%     trainPar.C_opt = str2double(input('Value of C: ','s'));
%     
%     optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold,subj)
%     optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold_all,1)
%     
% else
%     trainPar.C_opt = [];
%     
%     for n = 1:size(cfg.train,2)
%         subj = cfg.train(n);
%         optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold,subj)
%     end
%     
%     optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold_all,1)
% end

%% step 2: optimize SVM 
% we now train an SVMmodel with optimal thresholds as determined in step 1

%% get features
tStart = cell(size(cfg.train,2),1); tWidth = cell(size(cfg.train,2),1);
fStart = cell(size(cfg.train,2),1); fWidth = cell(size(cfg.train,2),1);
Area = cell(size(cfg.train,2),1);

for subj = 1:size(cfg.train,2)
    
    [Area{subj}, tStart{subj}, fStart{subj}, tWidth{subj}, fWidth{subj},...
        dataBase(subj).ERSP.Area, dataBase(subj).ERSP.tStart,  ...
        dataBase(subj).ERSP.fStart, dataBase(subj).ERSP.tWidth, ...
        dataBase(subj).ERSP.fWidth] = getfeaturesTrain(dataBase(cfg.train(subj)));
    
    fprintf('---- features are calculated in %s ----- \n',dataBase(subj).sub_label)
    
end

%% fit SVM
% vector [3*allERSPimages x 2]
tStart_all = vertcat(tStart{:});
tWidth_all = vertcat(tWidth{:});
fStart_all = vertcat(fStart{:});
fWidth_all = vertcat(fWidth{:});
Area_all = vertcat(Area{:});

% X-values in Support vector machine
X = [];
X(:,1:2) = Area_all;    % store area in X
X(:,3:4) = tStart_all;  % time of start of suppression
X(:,5:6) = fStart_all;  % minimal frequency of suppression
X(:,7:8) = tWidth_all;  % duration of suppression
X(:,9:10) = fWidth_all; % frequency range of suppression
     
% remove NaNs from X and Y
X_train = X;
idxnan = isnan(X_train);
idxnanrows = sum(idxnan,2);
X_train(idxnanrows>0,:) = [];
Y_alltrain(idxnanrows>0,:) = [];

% crossvalidaton
c = cvpartition(size(X_train,1),'KFold',10);
% optimalization variables
opts = struct('Optimizer','bayesopt','ShowPlots',true,'CVPartition',c,...
    'AcquisitionFunctionName','expected-improvement-plus');

% varied Costs
costvar(:,:,1) = [0 1; 1 0];
costvar(:,:,2) = [0 1; 2 0];
costvar(:,:,3) = [0 1; 3 0];
costvar(:,:,4) = [0 1; 4 0];
costvar(:,:,5) = [0 1; 5 0];

% pre-allocation
lossnew = NaN(1,size(costvar,3));
FN = NaN(1,size(costvar,3));
FP = NaN(1,size(costvar,3));
TP = NaN(1,size(costvar,3));
TN = NaN(1,size(costvar,3));
sens = NaN(1,size(costvar,3));
spec = NaN(1,size(costvar,3));
prec = NaN(1,size(costvar,3));
F = NaN(1,size(costvar,3));
MCC = NaN(1,size(costvar,3));

for i = 1:size(costvar,3)
    
    % train and crossvalidate SVM
    svmmod = fitcsvm(X_train,Y_alltrain,...
        'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',opts,'Cost',costvar(:,:,i));

    % predict visual rating
    label(:,i) = predict(svmmod, X_train); %#ok<SAGROW>
    
    % calculate loss
    if strcmp(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),'polynomial')
        lossnew(i) = kfoldLoss(fitcsvm(X_train,Y_alltrain,'CVPartition',c,...
            'KernelFunction',char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),...
            'PolynomialOrder',svmmod.HyperparameterOptimizationResults.XAtMinObjective.PolynomialOrder,...
            'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
            'Standardize',strcmpi(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.Standardize),'true')));
        
    elseif strcmp(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),'gaussian')
        
        lossnew(i) = kfoldLoss(fitcsvm(X_train,Y_alltrain,'CVPartition',c,...
            'KernelFunction',char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),...
            'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
            'KernelScale',svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelScale,...
            'Standardize',strcmpi(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.Standardize),'true')));
    end
    
    % calculate preformance (sensitivity, specificity, precision, F-score,
    % MMC) based on predicted visual rating
    TP(i) = numel(find(label(:,i) == 1 & Y_alltrain == 1));
    FP(i) = numel(find(label(:,i) == 1 & Y_alltrain == 0));
    TN(i) = numel(find(label(:,i) == 0 & Y_alltrain == 0));
    FN(i) = numel(find(label(:,i) == 0 & Y_alltrain == 1));
    
    if sum([TP(i);FP(i);TN(i);FN(i)]) ~= size(Y_alltrain,1)
        error('Number of TP, FP, TN, FN does not add up to all labels')
        
    else
        sens(i) = TP(i)/(TP(i)+FN(i));
        spec(i) = TN(i)/(TN(i)+FP(i));
        prec(i) = TP(i)/(TP(i)+FP(i));
        F(i) = 2 * (prec(i) * sens(i))/(prec(i) + sens(i));
        MCC(i) = ((TP(i)*TN(i)) - (FP(i)*FN(i)))/sqrt((TP(i)+FP(i))*(TP(i)+FN(i))*(TN(i)+FP(i))*(TN(i)+FN(i)));
    end
    
    svmmodsall(i).svmmod = svmmod; %#ok<SAGROW>
    svmmodsall(i).cost = costvar(:,:,i); %#ok<SAGROW>
    
    fprintf('----- Calculated SVM with Costs [%d %d; %d %d] ----- \n',costvar(1,1,i), costvar(1,2,i), costvar(2,1,i), costvar(2,2,i))
end

% %% calculate performance
% label(:,1) = predict(svmmod11,X_train);
% label(:,2) = predict(svmmod12,X_train);
% label(:,3) = predict(svmmod13,X_train);
% label(:,4) = predict(svmmod14,X_train);
% 
% for i=1:size(label,2)
%     TP(i) = numel(find(label(:,i) == 1 & Y_alltrain == 1));
%     FP(i) = numel(find(label(:,i) == 1 & Y_alltrain == 0));
%     TN(i) = numel(find(label(:,i) == 0 & Y_alltrain == 0));
%     FN(i) = numel(find(label(:,i) == 0 & Y_alltrain == 1));
%     if sum([TP(i);FP(i);TN(i);FN(i)]) ~= size(Y_alltrain,1)
%         error('Number of TP, FP, TN, FN does not add up to all labels')
%     else
%         sens(i) = TP(i)/(TP(i)+FN(i));
%         spec(i) = TN(i)/(TN(i)+FP(i));
%         prec(i) = TP(i)/(TP(i)+FP(i));
%         F(i) = 2 * (prec(i) * sens(i))/(prec(i) + sens(i));
%         MCC(i) = ((TP(i)*TN(i)) - (FP(i)*FN(i)))/sqrt((TP(i)+FP(i))*(TP(i)+FN(i))*(TN(i)+FP(i))*(TN(i)+FN(i)));
%     end
% end

%% combine all figures, and variables

allERSP = [];
for i=1:3
    allERSP = [allERSP; dataBase(i).ERSP.allERSPboot(:)]; %#ok<AGROW>
end

exclude = zeros(size(allERSP));
for i=1:size(allERSP,1)
    if any(any(isnan(allERSP{i} )))
        exclude(i) = 1;
    end
end

allERSPnew = allERSP(idxnanrows==0);
Area_allnew = Area_all(idxnanrows==0,:);
fStart_allnew = fStart_all(idxnanrows==0,:);
fWidth_allnew = fWidth_all(idxnanrows==0,:);
tStart_allnew = tStart_all(idxnanrows==0,:);
tWidth_allnew = tWidth_all(idxnanrows==0,:);

%% check FN or FP

% mode = input('What would you like to check? [FP/FN/TP/TN]: ','s');
t(1) = find(times<0,1,'first');
t(2) = find(times>0,1,'first');

for i=randperm(size(label,1))
    
    if any(label(i,[2,3]) ~= Y_alltrain(i))
        
        figure(1),
        a = allERSPnew{i};
        imagesc([times(1) times(end)],[freqs(1) freqs(end)],a,[-15,15])
        colormap jet
        ax = gca;
        ax.YDir = 'normal';
        pixelWidth = (times(end)-times(1))/(size(times,2)-1);
        ax.XLim = [times(1)-0.5*pixelWidth times(end)+0.5*pixelWidth];
        xlabel('Time(ms)')
        ylabel('Frequency (Hz)')
        
        hold on
        
        C = [0.7 0.7 0.7]; % light grey
        
        for win=1:2
            if Area_allnew(i,win) >0
                pixelWidth = (times(end)-times(1))/(size(times,2)-1);
                pixelHeigth = (freqs(end)-freqs(1))/(size(freqs,2)-1);
                
                % plot bounding box on largest power suppression
                if tStart_allnew(i,win) > 0.5
                    x1 = (times(floor(tStart_allnew(i,win)+t(win)-1)) + times(ceil(tStart_allnew(i,win)+t(win)-1)))/2;
                else
                    x1 = times(1)-0.5*pixelWidth;
                end
                
                if fStart_allnew(i,win) > 0.5
                    y1 = (freqs(floor(fStart_allnew(i,win))) + freqs(ceil(fStart_allnew(i,win))))/2;
                else
                    y1 = freqs(1)-0.5*pixelHeigth;
                end
                w = tWidth_allnew(i,win)*pixelWidth;
                h = fWidth_allnew(i,win)*pixelHeigth;
                
                rectangle('Position',[x1 y1 w h],'EdgeColor',C),
            end
        end
        
        str = cell(1,4);
        for n=1:4
            if Y_alltrain(i)==1 && label(i,n) == 1
                str{1,n} = 'TP';
            elseif Y_alltrain(i)==1 && label(i,n) == 0
                str{1,n} = 'FN';
            elseif Y_alltrain(i)==0 && label(i,n) == 1
                str{1,n} = 'FP';
            elseif Y_alltrain(i)== 0 && label(i,n) == 0
                str{1,n} = 'TN';
            end
        end
        
        title(sprintf('mod11 = %s, mod12 = %s, mod13 = %s, mod14 = %s',str{:}))
        pause
    end
end


%% save SVMmodel

SVMModel = svmmod14;
pathname = cfg.SVMpath;
filename = sprintf('SVMmodel_trained_BB_%s.mat',datestr(now,'yyyymmdd'));
save(fullfile(pathname,filename),'SVMModel')
fprintf('SVMmodel is saved in %s\n',fullfile(pathname,filename))

%% test SVM model with one patient: RESP0690

% detectPowSup(dataBase,cfg)
subj = cfg.test;

[Area, tStart, fStart, tWidth, fWidth,...
    dataBase(subj).ERSP.Area, dataBase(subj).ERSP.tStart,  ...
        dataBase(subj).ERSP.fStart, dataBase(subj).ERSP.tWidth, ...
        dataBase(subj).ERSP.fWidth] = getfeaturesTrain(dataBase(cfg.test));

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

% label = str2num(cell2mat(label_raw)); %#ok<ST2NM>
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
