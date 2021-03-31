%% script to construct svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

%% settings for constructing SVM
clc
clear
myDataPath = setLocalDataPath(1);

%%

% old database: PAT119,  PAT126, PAT130, PAT135
cfg.sub_labels = { 'sub-RESP0607', 'sub-RESP0638','sub-RESP0676','sub-RESP0690'};
cfg.ses_label = 'ses-1';
cfg.task_label = 'task-SPESclin';
cfg.run_label = {'run-031211','run-031049','run-021423','run-041139'};
cfg.train = 1:3;
cfg.test = 4;

%% Load ERSP-data
% these ERSPs are made and saved in ccepSp00_makeTFSPES

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg.sub_labels,2)
    
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    % load ERSPs
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

clear channelsName tb_channels idx_ch_incl ch_ecog idx_ch_bad

disp('All ERSPs are loaded')

%% load visual ratings
% These visual ratings are performed in ccepSP00_visualRateSp twice (two
% different observers). All ERSPs where rating was not in agreement, were
% checked again (ccepSp00_visualRateSpCombined), and these consent ratings
% are loaded in this part

for subj = 1:size(dataBase,2)
    
    folderName = fullfile(myDataPath.dir_visrate,[cfg.sub_labels{subj},...
        '_',cfg.ses_label,'_',cfg.run_label{subj},'_visBB_combined.mat']);
    
    dataBase(subj).visBB = load(folderName);

end

clear folderName

disp('All visual ratings are loaded')

%% determine true labels for fitcsvm: y 
% the visual ratings loaded in the previous section, are combined with the
% ERSPs and converted to labels: 
% 1 = power suppression was observed, 
% 0 = no power suppression was observed.

[Y_alltrain,Y_alltest] = determine_Target_SVM(dataBase,cfg);

%% get features
% this function extracts features of each ERSP.
% These features are: 
%   - Area_conc: size of the area with the largest power suppression - pre
%   and post stimulus [concatenated ERSPs x 2], first column is
%   pre-stimulus, second column is post-stimulus
%   - tStart_conc: start in time of this largest power suppression - pre
%   and post stimulus [concatenated ERSPs x 2], first column is
%   pre-stimulus, second column is post-stimulus
%   - fStart_conc: start in frequency of this largest power suppression -
%   pre and post stimulus [concatenated ERSPs x 2], first column is
%   pre-stimulus, second column is post-stimulus
%   - tWidth_conc: the duration of this largest power suppression - pre and
%   post stimulus [concatenated ERSPs x 2], first column is pre-stimulus,
%   second column is post-stimulus
%   -  fWidth_conc: the range in frequencies of this largest power
%   suppression - pre and post stimulus [concatenated ERSPs x 2], first
%   column is pre-stimulus, second column is post-stimulus

for subj = 1:size(cfg.train,2)
    
    detPar_train(subj) = getfeaturesTrain(dataBase(cfg.train(subj))); %#ok<SAGROW>
    
    fprintf('---- features are calculated in %s ----- \n',dataBase(subj).sub_label)
    
end

%% fit SVM - this takes a lot of time! 

% concatenate all features into five vectors [3*allERSPimages x 2]
tStart_train = vertcat(detPar_train(:).tStart_conc);
tWidth_train = vertcat(detPar_train(:).tWidth_conc);
fStart_train = vertcat(detPar_train(:).fStart_conc);
fWidth_train = vertcat(detPar_train(:).fWidth_conc);
Area_train = vertcat(detPar_train(:).Area_conc);

% store these features as X-values to use for fitting of the Support vector
% machine
X = [];
X(:,1:2) = Area_train;    % store area in X
X(:,3:4) = tStart_train;  % time of start of suppression
X(:,5:6) = fStart_train;  % minimal frequency of suppression
X(:,7:8) = tWidth_train;  % duration of suppression
X(:,9:10) = fWidth_train; % frequency range of suppression
     
% remove NaNs from X and Y --> these NaNs are when an ERSP is calculated in
% an electrode which is part of the current stimulus pair
X_train = X;
idxnan = isnan(X_train);
idxnanrows = sum(idxnan,2);
X_train(idxnanrows>0,:) = [];
Y_train = Y_alltrain;
Y_train(idxnanrows>0,:) = [];

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

% pre-allocation of the performance calculations for varying Costs
loss_train = NaN(1,size(costvar,3));
FN_train = NaN(1,size(costvar,3));
FP_train = NaN(1,size(costvar,3));
TP_train = NaN(1,size(costvar,3));
TN_train = NaN(1,size(costvar,3));
sens_train = NaN(1,size(costvar,3));
spec_train = NaN(1,size(costvar,3));
prec_train = NaN(1,size(costvar,3));
F_train = NaN(1,size(costvar,3));
MCC_train = NaN(1,size(costvar,3));

% for each Cost, optimize the SVM
for i = 1:size(costvar,3)
    
    % train and crossvalidate SVM
    svmmod = fitcsvm(X_train,Y_train,...
        'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',opts,'Cost',costvar(:,:,i));

    % predict visual rating
    label_train(:,i) = predict(svmmod, X_train); %#ok<SAGROW>
    
    % calculate loss
    if strcmp(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),'polynomial')
        loss_train(i) = kfoldLoss(fitcsvm(X_train,Y_train,'CVPartition',c,...
            'KernelFunction',char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),...
            'PolynomialOrder',svmmod.HyperparameterOptimizationResults.XAtMinObjective.PolynomialOrder,...
            'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
            'Standardize',strcmpi(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.Standardize),'true')));
        
    elseif strcmp(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),'gaussian')
        
        loss_train(i) = kfoldLoss(fitcsvm(X_train,Y_train,'CVPartition',c,...
            'KernelFunction',char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),...
            'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
            'KernelScale',svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelScale,...
            'Standardize',strcmpi(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.Standardize),'true')));
    
    elseif strcmp(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),'linear')
        
        loss_train(i) = kfoldLoss(fitcsvm(X_train,Y_train,'CVPartition',c,...
            'KernelFunction',char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelFunction),...
            'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
            'Standardize',strcmpi(char(svmmod.HyperparameterOptimizationResults.XAtMinObjective.Standardize),'true')));
    end
    
    % calculate preformance (sensitivity, specificity, precision, F-score,
    % MMC) based on predicted visual rating
    TP_train(i) = numel(find(label_train(:,i) == 1 & Y_train == 1));
    FP_train(i) = numel(find(label_train(:,i) == 1 & Y_train == 0));
    TN_train(i) = numel(find(label_train(:,i) == 0 & Y_train == 0));
    FN_train(i) = numel(find(label_train(:,i) == 0 & Y_train == 1));
    
    if sum([TP_train(i);FP_train(i);TN_train(i);FN_train(i)]) ~= size(Y_train,1)
        error('Number of TP, FP, TN, FN does not add up to all labels')
        
    else
        sens_train(i) = TP_train(i)/(TP_train(i)+FN_train(i));
        spec_train(i) = TN_train(i)/(TN_train(i)+FP_train(i));
        prec_train(i) = TP_train(i)/(TP_train(i)+FP_train(i));
        F_train(i) = 2 * (prec_train(i) * sens_train(i))/(prec_train(i) + sens_train(i));
        MCC_train(i) = ((TP_train(i)*TN_train(i)) - (FP_train(i)*FN_train(i)))/sqrt((TP_train(i)+FP_train(i))*(TP_train(i)+FN_train(i))*(TN_train(i)+FP_train(i))*(TN_train(i)+FN_train(i)));
    end
    
    svmmodstrainall(i).svmmod = svmmod; %#ok<SAGROW>
    svmmodstrainall(i).cost = costvar(:,:,i); %#ok<SAGROW>
    
    fprintf('----- Calculated SVM with Costs [%d %d; %d %d] ----- \n',costvar(1,1,i), costvar(1,2,i), costvar(2,1,i), costvar(2,2,i))
end

clear idxnan idxnanrows X svmmod

%% combine all figures, and  variables and check which SVM build with 
% specific Cost, would lead to which decision (TP/FP/FN/TN)

% combine all figures, and variables
allERSP_train = [];
for i = cfg.train
    allERSP_train = [allERSP_train; dataBase(i).ERSP.allERSPboot(:)]; %#ok<AGROW>
end

% these NaNs are when calculated in a channel part of the current stimulus
% pair
idxnan = isnan(Area_train);
idxnanrows = sum(idxnan,2);

allERSPtrain_rmnan = allERSP_train(idxnanrows == 0);
Area_train_rmnan = Area_train(idxnanrows ==0,:);
fStart_train_rmnan = fStart_train(idxnanrows ==0,:);
fWidth_train_rmnan = fWidth_train(idxnanrows ==0,:);
tStart_train_rmnan = tStart_train(idxnanrows ==0,:);
tWidth_train_rmnan = tWidth_train(idxnanrows ==0,:);

if size(label_train,1) ~= size(allERSPtrain_rmnan,1)
   error('Number of ERSPs is not equal to number of ERSPs in which power suppressions was detected!') 
end

% check FN or FP
times = dataBase(1).ERSP.times;
freqs = dataBase(1).ERSP.freqs;

t(1) = find(times<0,1,'first');
t(2) = find(times>0,1,'first');

for i = randperm(size(label_train,1)) % randomly show ERSPs from all used patients
    
    if any(label_train(i,:) ~= Y_train(i))
        
        figure(1),
        a = allERSPtrain_rmnan{i};
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
            if Area_train_rmnan(i,win) >0
                pixelWidth = (times(end)-times(1))/(size(times,2)-1);
                pixelHeigth = (freqs(end)-freqs(1))/(size(freqs,2)-1);
                
                % plot bounding box on largest power suppression
                if tStart_train_rmnan(i,win) > 0.5
                    x1 = (times(floor(tStart_train_rmnan(i,win)+t(win)-1)) + times(ceil(tStart_train_rmnan(i,win)+t(win)-1)))/2;
                else
                    x1 = times(1)-0.5*pixelWidth;
                end
                
                if fStart_train_rmnan(i,win) > 0.5
                    y1 = (freqs(floor(fStart_train_rmnan(i,win))) + freqs(ceil(fStart_train_rmnan(i,win))))/2;
                else
                    y1 = freqs(1)-0.5*pixelHeigth;
                end
                w = tWidth_train_rmnan(i,win)*pixelWidth;
                h = fWidth_train_rmnan(i,win)*pixelHeigth;
                
                rectangle('Position',[x1 y1 w h],'EdgeColor',C),
            end
        end
        
        str = cell(1,size(label_train,2));
        for n=1:size(label_train,2)
            if Y_train(i)==1 && label_train(i,n) == 1
                str{1,n} = 'TP';
            elseif Y_train(i)==1 && label_train(i,n) == 0
                str{1,n} = 'FN';
            elseif Y_train(i)==0 && label_train(i,n) == 1
                str{1,n} = 'FP';
            elseif Y_train(i)== 0 && label_train(i,n) == 0
                str{1,n} = 'TN';
            end
        end
        
        title(sprintf('mod11 = %s, mod12 = %s, mod13 = %s, mod14 = %s, mod15 = %s',str{:}))
        pause
    end
end

% small cleanup
clear a ax C freqs h i idxnan idxnanrows n pixelHeigth pixelWidth str t times w win x1 y1

%% save SVMmodel

% choose svm with sensitivity >0.9 and largest F-score
idx_svmopt = sens_train>0.9;
idx_svm = (F_train == max(F_train(idx_svmopt)));
SVMModel = svmmodstrainall(idx_svm).svmmod;
cost = svmmodstrainall(idx_svm).cost;

% save SVM
pathname = myDataPath.SVMpath;
filename = sprintf('SVMmodel_trained_BB_%s.mat',datestr(now,'yyyymmdd'));
save(fullfile(pathname,filename),'SVMModel','cost')
fprintf('SVMmodel is saved in %s\n',fullfile(pathname,filename))

% small cleanup
clear idx_svmopt idx_svm cost pathname filename

%% find features in each ERSP

detPar_test = getfeaturesTrain(dataBase(cfg.test)); 

%% test SVM model with one patient: RESP0690

% vector [3*allERSPimages x 2]
tStart_test = vertcat(detPar_test(:).tStart_conc);
tWidth_test = vertcat(detPar_test(:).tWidth_conc);
fStart_test = vertcat(detPar_test(:).fStart_conc);
fWidth_test = vertcat(detPar_test(:).fWidth_conc);
Area_test = vertcat(detPar_test(:).Area_conc);

% X-values in Support vector machine
X = [];
X(:,1:2) = Area_test;    % store area in X
X(:,3:4) = tStart_test;  % time of start of suppression
X(:,5:6) = fStart_test;  % minimal frequency of suppression
X(:,7:8) = tWidth_test;  % duration of suppression
X(:,9:10) = fWidth_test; % frequency range of suppression

% remove NaNs from X and Y
X_test_rmnan = X;
idxnan = isnan(X_test_rmnan);
idxnanrows = sum(idxnan,2);
X_test_rmnan(idxnanrows>0,:) = [];

Y_test_rmnan = Y_alltest;
Y_test_rmnan(idxnanrows>0) = [];

label_test_rmnan = predict(SVMModel,X_test_rmnan);

% put all detected values which were in stimulus pairs (NaNs) to 0 
% label_test(idxnanrows>0) = [];

% label = str2num(cell2mat(label_raw)); %#ok<ST2NM>
TP_test = numel(find(label_test_rmnan == 1 & Y_test_rmnan == 1));
FP_test = numel(find(label_test_rmnan == 1 & Y_test_rmnan == 0));
TN_test = numel(find(label_test_rmnan == 0 & Y_test_rmnan == 0));
FN_test = numel(find(label_test_rmnan == 0 & Y_test_rmnan == 1));

if sum([TP_test;FP_test;TN_test;FN_test]) ~= size(Y_test_rmnan,1)
    error('Number of TP, FP, TN, FN does not add up to all labels')
else
    sens_test = TP_test/(TP_test+FN_test);
    spec_test = TN_test/(TN_test+FP_test);
    prec_test = TP_test/(TP_test+FP_test);
    F_test = 2 * (prec_test * sens_test)/(prec_test + sens_test);
    MCC_test = ((TP_test*TN_test) - (FP_test*FN_test))/sqrt((TP_test+FP_test)*(TP_test+FN_test)*(TN_test+FP_test)*(TN_test+FN_test));
end

%% check FN or FP

% vector [3*allERSPimages x 2]
tStart_test = vertcat(detPar_test.tStart_conc);
tWidth_test = vertcat(detPar_test.tWidth_conc);
fStart_test = vertcat(detPar_test.fStart_conc);
fWidth_test = vertcat(detPar_test.fWidth_conc);
Area_test = vertcat(detPar_test.Area_conc);

% X-values in Support vector machine
X = [];
X(:,1:2) = Area_test;    % store area in X
X(:,3:4) = tStart_test;  % time of start of suppression
X(:,5:6) = fStart_test;  % minimal frequency of suppression
X(:,7:8) = tWidth_test;  % duration of suppression
X(:,9:10) = fWidth_test; % frequency range of suppression

label_test = predict(SVMModel,X);

idxnan = isnan(X);
idxnanrows = sum(idxnan,2);

% put all detected values which were in stimulus pairs (NaNs) to 0
label_test(idxnanrows>0) = 0;

dataBase(subj).ERSPdet.sub_label = dataBase(subj).sub_label;
dataBase(subj).ERSPdet.ses_label = dataBase(subj).ses_label;
dataBase(subj).ERSPdet.task_label = dataBase(subj).task_label;
dataBase(subj).ERSPdet.run_label = dataBase(subj).run_label;
dataBase(subj).ERSPdet.detected = reshape(label_test,size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.Area{1} = reshape(Area_test(:,1),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.Area{2} = reshape(Area_test(:,2),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.tStart{1} = reshape(tStart_test(:,1),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.tStart{2} = reshape(tStart_test(:,2),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.tWidth{1} = reshape(tWidth_test(:,1),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.tWidth{2} = reshape(tWidth_test(:,2),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.fStart{1} = reshape(fStart_test(:,1),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.fStart{2} = reshape(fStart_test(:,2),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.fWidth{1} = reshape(fWidth_test(:,1),size(dataBase(subj).ERSP.allERSPboot));
dataBase(subj).ERSPdet.fWidth{2} = reshape(fWidth_test(:,2),size(dataBase(subj).ERSP.allERSPboot));

mode = input('What would you like to check? [FP/FN/TP/TN]: ','s');
subj = cfg.test;

Y_test = Y_alltest;

if strcmp(mode,'FP')
    lookPic = reshape((label_test==1 & Y_test == 0),size(dataBase(subj).ERSP.allERSP));
elseif strcmp(mode,'FN')
    lookPic = reshape((label_test==0 & Y_test == 1),size(dataBase(subj).ERSP.allERSP));
elseif strcmp(mode,'TP')
    lookPic = reshape((label_test==1 & Y_test == 1),size(dataBase(subj).ERSP.allERSP));   
elseif strcmp(mode,'TN')
     lookPic = reshape((label_test==0 & Y_test == 0),size(dataBase(subj).ERSP.allERSP));   
end

endstimp = 0;

visualRating_tfspes(dataBase,subj,lookPic,[],endstimp);


