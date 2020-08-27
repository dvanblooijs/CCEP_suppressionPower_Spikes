%% visually score spikes
% This script is used to annotate interictal epileptic discharges in ECoG
% data

%% config

clc
clear
cfg = setLocalDataPath(1);

%% patient characteristics

cfg.sub_labels = {['sub-' input('Patient number (RESPXXXX): ','s')]};
cfg.ses_label = input('Session number (ses-X): ','s');
cfg.task_label = 'task-SPESclin';
cfg.run_label = {['run-' input('Run [daydayhhminmin]: ','s')]};
cfg.IEDch = input('Channels showing interictal epileptic discharges: ','s');

%% load ECoGs with SPES from 1 patient

dataBase = load_ECoGdata(cfg);

%% remove stimulation artefact

dataBase = removeStimArt(dataBase);

fprintf('...All subjects has been run...\n')

%% find IED channels

dataBase.IEDchan = strsplit(cfg.IEDch,',');
dataBase.IED = NaN(size(dataBase.IEDchan));
for i=1:size(dataBase.IEDchan,2)
    charIED = regexp(lower(dataBase.IEDchan{i}),'[a-z]');
    numIED = regexp(dataBase.IEDchan{i},'[0-9]');
    
    type1 = [dataBase.IEDchan{i}(charIED),dataBase.IEDchan{i}(numIED)];
    type2 = [dataBase.IEDchan{i}(charIED),'0',dataBase.IEDchan{i}(numIED)];
    if sum(strcmpi(dataBase.ch,type1))
        dataBase.IED(i) = find(strcmpi(dataBase.ch,type1)==1);
    elseif sum(strcmpi(dataBase.ch,type2))
        dataBase.IED(i) = find(strcmpi(dataBase.ch,type2)==1);
    else
        error('Channel %s or %s is not found',type1,type2)
    end
end

%% visually score IEDs

vis_spikes = visualRateSpikes(dataBase);

%% reshape vis_spikes

runs = ceil(size(dataBase.IED,2)/9);

if runs == 3
x_all = [vis_spikes(1).x_all, zeros(size(vis_spikes(1).x_all,1),size(vis_spikes(2).x_all,2)), zeros(size(vis_spikes(1).x_all,1),size(vis_spikes(3).x_all,2));...
    zeros(size(vis_spikes(2).x_all,1),size(vis_spikes(1).x_all,2)), vis_spikes(2).x_all,  zeros(size(vis_spikes(2).x_all,1),size(vis_spikes(3).x_all,2));...
    zeros(size(vis_spikes(3).x_all,1),size(vis_spikes(1).x_all,2)), zeros(size(vis_spikes(3).x_all,1),size(vis_spikes(2).x_all,2)), vis_spikes(3).x_all];

y_all = [vis_spikes(1).y_all, zeros(size(vis_spikes(1).y_all,1),size(vis_spikes(2).y_all,2)), zeros(size(vis_spikes(1).y_all,1),size(vis_spikes(3).y_all,2));...
    zeros(size(vis_spikes(2).y_all,1),size(vis_spikes(1).y_all,2)), vis_spikes(2).y_all,  zeros(size(vis_spikes(2).y_all,1),size(vis_spikes(3).y_all,2));...
    zeros(size(vis_spikes(3).y_all,1),size(vis_spikes(1).y_all,2)), zeros(size(vis_spikes(3).y_all,1),size(vis_spikes(2).y_all,2)), vis_spikes(3).y_all];

elseif runs ==2
    x_all = [vis_spikes(1).x_all, zeros(size(vis_spikes(1).x_all,1),size(vis_spikes(2).x_all,2));...
        zeros(size(vis_spikes(2).x_all,1),size(vis_spikes(1).x_all,2)), vis_spikes(2).x_all];
    
    y_all = [vis_spikes(1).y_all, zeros(size(vis_spikes(1).y_all,1),size(vis_spikes(2).y_all,2));...
        zeros(size(vis_spikes(2).y_all,1),size(vis_spikes(1).y_all,2)), vis_spikes(2).y_all];
    
end
    

%% save all annotated spikes

folderName = cfg.visIEDpath;
fileName = [dataBase(subj).sub_label, '_', dataBase(subj).ses_label, '_',...
    dataBase(subj).task_label, '_', dataBase(subj).run_label, '_visIEDs.mat'];

sub_label = dataBase(subj).sub_label;
ses_label = dataBase(subj).ses_label;
task_label = dataBase(subj).task_label;
run_label = dataBase(subj).run_label;
IED = dataBase(subj).IED;
IEDchan = dataBase(subj).IEDchan;

save([folderName, fileName],'sub_label','ses_label','task_label','run_label',...
    'IED','IEDchan','x_all','y_all')

disp('saved')


