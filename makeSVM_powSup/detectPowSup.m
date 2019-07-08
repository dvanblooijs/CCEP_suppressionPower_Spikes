%% use svm to detect power suppression in TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

clear all
close all

addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/hysteresis3d'))
addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/output/OUTPUT_PAT123'))
addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/output/SVM'))
addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/Data'))

%change data to patient you want to detect SP for

load SVMmodel9                          % get trained SVM model 
load times                              % for time resolution
load freqs
load ('Output_PAT_123','allERSP2_123')  % get ERSP matrix to calculate features
load ('PAT_123_allstimp','stimpchan')

[D,A] = getfeatures(times, freqs, allERSP2_123, 4, 5, stimpchan);        % Use 4 and 5 as lower and upper threshold hysteresis (optimum decided via cross-validation)

%%features in X
X(:,1) = A';                            % Put area in matrix X
X(:,2) = D';                            % Put duration in matrix X

[labels,scores] = predict(SVMModel,X);  % Predict scores by SVM with features area and duration

labels = cellfun(@str2num, labels);     
Label_final=vec2mat(labels,size(allERSP2_123,2));   % Scored by SVM 

%save data
fileName='Label_final';
targetFolder = '/Users/michelle/Google Drive/M3/Onderzoek/Matlab/output/OUTPUT_PAT123/';
save([targetFolder,fileName], 'Label_final');       % Save file in targetfolder 


%%Manually check Label_final with figures. Adjust matrix and save as:
%%Label_final_adjusted in the same folder