%% compare visualratings of suppressed power


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

%% Load ERSP-data

%pre-allocation
dataBase = struct([]);

for subj = 1:size(cfg.sub_labels,2)
    
    dataBase(subj).sub_label = cfg.sub_labels{subj};
    dataBase(subj).ses_label = cfg.ses_label;
    dataBase(subj).task_label = cfg.task_label;
    dataBase(subj).run_label = cfg.run_label{subj};
    
    dataBase(subj).ERSP = load(fullfile(myDataPath.TFSPESoutput, dataBase(subj).sub_label,...
        dataBase(subj).ses_label,dataBase(subj).run_label,...
        [dataBase(subj).sub_label, '_',dataBase(subj).ses_label,'_' dataBase(subj).task_label,...
        '_', dataBase(subj).run_label,'_ERSP.mat']));
end

disp('All ERSPs are loaded')

%% load visual rating

for subj = 1:size(dataBase,2)
    folderName1 = fullfile(myDataPath.dir_visrate,[cfg.sub_labels{subj},'_',cfg.ses_label,'_',cfg.run_label{subj},'_visBB_DvB.mat']);
    folderName2 = fullfile(myDataPath.dir_visrate,[cfg.sub_labels{subj},'_',cfg.ses_label,'_',cfg.run_label{subj},'_visBB_MvdS.mat']);

    dataBase(subj).visBB1 = load(folderName1);
    dataBase(subj).visBB2 = load(folderName2);
   
end

%% calculate kappa



%% check visual rating of suppressed power where there is no consensus
% - There should be no blue region before the stimulus.
% - There should be a clear blue region after the stimulus. 

subj = input('Which subject would you like to annotate?: '); 
vis_BB_combined = zeros(size(dataBase(subj).ERSP.allERSPboot));
n=1;

for stimp = 1:size(dataBase(subj).ERSP.allERSPboot,1)
    for chan = 1:size(dataBase(subj).ERSP.allERSPboot,2)
        
        if dataBase(subj).visBB1.vis_BB(stimp,chan) ~= dataBase(subj).visBB2.vis_BB(stimp,chan)
            fig = 1;
            plot_ERSP(dataBase(subj).ERSP,stimp,chan,fig)
            
            perc = n / size(dataBase(subj).ERSP.allERSPboot(:),1) *100;
            
            s = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this an ERSP? [y/n]: ',perc,dataBase(subj).ERSP.cc_stimchans{stimp,:},dataBase(subj).ERSP.ch{chan}),'s');
            
            if strcmp(s,'y')
                vis_BB_combined(stimp,chan) = 1;
            end
        else
            vis_BB_combined(stimp,chan) = dataBase(subj).visBB1.vis_BB(stimp,chan);
            
        end
        n=n+1;
        
    end
end

disp('Visual rating in completed')

%% save combined visual rating

folderName = fullfile(myDataPath.dir_visrate,[cfg.sub_labels{subj},'_',cfg.ses_label,'_',cfg.run_label{subj},'_visBB_combined.mat']);

sub_label = dataBase(subj).sub_label;
ses_label = dataBase(subj).ses_label;
task_label = dataBase(subj).task_label;
run_label = dataBase(subj).run_label;
cc_stimchans = dataBase(subj).ERSP.cc_stimchans;
cc_stimsets = dataBase(subj).ERSP.cc_stimsets;
ch = dataBase(subj).ERSP.ch;

save(folderName,'sub_label','ses_label','task_label','run_label',...
    'vis_BB_combined','ch','cc_stimsets','cc_stimchans')

fprintf('saved vis_BB in %s \n',folderName)




