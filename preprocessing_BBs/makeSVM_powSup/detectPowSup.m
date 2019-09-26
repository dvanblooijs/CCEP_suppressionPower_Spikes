%% use svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

function detectPowSup(dataBase,cfg)

load SVMmodel9                          % get trained SVM model

for subj = 1:size(dataBase,2)
    % addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/hysteresis3d'))
    % addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/output/OUTPUT_PAT123'))
    % addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/output/SVM'))
    % addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/Data'))
    
    %change data to patient you want to detect SP for
    load([cfg.inputERSP, dataBase(subj).subj, '_ERSP'])
    
    [D,A] = getfeatures(times, freqs, ERSP, 4, 5, dataBase(subj).cc_stimsets);        % Use 4 and 5 as lower and upper threshold hysteresis (optimum decided via cross-validation)
    
    %%features in X
    X(:,1) = A';                            % Put area in matrix X
    X(:,2) = D';                            % Put duration in matrix X
    
    [labels,~] = predict(SVMModel,X);  % Predict scores by SVM with features area and duration
    
    labels = cellfun(@str2num, labels);
    Label_final=vec2mat(labels,size(ERSP,2));   % Scored by SVM
    
    %save data
    save([cfg.outputdetPowSup,dataBase(subj).subj,'_detectedPowSup'], 'Label_final');       % Save file in targetfolder
    
end
%%Manually check Label_final with figures. Adjust matrix and save as:
%%Label_final_adjusted in the same folder