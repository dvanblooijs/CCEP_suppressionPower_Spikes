% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_ECoGdata(myDataPath,cfg)

dataPath = myDataPath.proj_dirinput;
dataBase = struct([]);

for nSubj = 1:size(cfg,2)
    sub_label = cfg(nSubj).sub_labels;
    ses_label = cfg(nSubj).ses_label;
    task_label = cfg(nSubj).task_label;
    run_label = cfg(nSubj).run_label{1};

    D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
        [sub_label '_' ses_label '_' task_label ,'_',run_label, '_ieeg.eeg']));
   
    dataName = fullfile(D(1).folder, D(1).name);

    ccep_data = ft_read_data(dataName,'dataformat','brainvision_eeg');
    ccep_header = ft_read_header(dataName);

    % load events
    D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
        [sub_label '_' ses_label '_' task_label ,'_',run_label,'_events.tsv']));

    eventsName = fullfile(D(1).folder, D(1).name);

    tb_events = readtable(eventsName,'FileType','text','Delimiter','\t');

    % load electrodes
    D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
        [sub_label '_' ses_label ,'_electrodes.tsv']));

    elecsName = fullfile(D(1).folder, D(1).name);

    tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
    idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
    tb_electrodes = tb_electrodes(idx_elec_incl,:);

    % load channels
    D = dir(fullfile(dataPath,sub_label, ses_label,'ieeg',...
        [sub_label '_' ses_label '_' task_label ,'_',run_label,'_channels.tsv']));

    channelsName = fullfile(D(1).folder, D(1).name);

    tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
    idx_ch_incl = strcmp(tb_channels.type,'ECOG')|strcmp(tb_channels.type,'SEEG');

    tb_channels = tb_channels(idx_ch_incl,:);
    ch_incl = tb_channels.name;

    data = ccep_data(idx_ch_incl,:);

    dataBase(nSubj).sub_label = sub_label;
    dataBase(nSubj).ses_label = ses_label;
    dataBase(nSubj).task_label = task_label;
    dataBase(nSubj).run_label = run_label;
    dataBase(nSubj).dataName = dataName;
    dataBase(nSubj).ccep_header = ccep_header;
    dataBase(nSubj).tb_events = tb_events;
    dataBase(nSubj).tb_channels = tb_channels;
    dataBase(nSubj).tb_electrodes = tb_electrodes;
    dataBase(nSubj).ch = ch_incl;
    dataBase(nSubj).data = data;
    fprintf('...Subject %s has been run...\n',sub_label)
end

disp('All subjects are loaded')
