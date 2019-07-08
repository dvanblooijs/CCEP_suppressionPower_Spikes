% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_data(cfg)

dataPath = cfg.dataPath;
sub_labels = cfg.sub_labels;
ses_label = cfg.ses_label;
task_label = cfg.task_label;

dataBase = struct([]);
for i=1:size(sub_labels,2)
    D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
        ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_ieeg.eeg']));
    
    dataName = fullfile(D(1).folder, D(1).name);
    
    ccep_data = ft_read_data(dataName,'dataformat','brainvision_eeg');
    ccep_header = ft_read_header(dataName);
    
    % load channels
    D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
        ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_events.tsv']));
    
    eventsName = fullfile(D(1).folder, D(1).name);
    
    tb_events = readtable(eventsName,'FileType','text','Delimiter','\t');
    
    % load events
    D = dir(fullfile(dataPath,['sub-' sub_labels{i}],['ses-' ses_label],'ieeg',...
        ['sub-' sub_labels{i} '_ses-' ses_label '_task-' task_label ,'_run-*_channels.tsv']));
    
    channelsName = fullfile(D(1).folder, D(1).name);
    
    tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
    log_ch_incl = strcmp(tb_channels.type,'ECOG');
    
    ch = tb_channels.name;
    ch_incl = {ch{log_ch_incl}}';
    
    data = -1*ccep_data(log_ch_incl,:);
    
    dataBase(i).subj = sub_labels{i};
    dataBase(i).dataName = dataName;
    dataBase(i).ccep_header = ccep_header;
    dataBase(i).tb_events = tb_events;
    dataBase(i).tb_channels = tb_channels;
    dataBase(i).ch = ch_incl;
    dataBase(i).data = data;
    fprintf('...Subject %s has been run...\n',sub_labels{i})
end

disp('All subjects are loaded')
