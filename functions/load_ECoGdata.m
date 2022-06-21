% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_ECoGdata(myDataPath,cfg)

dataPath = myDataPath.dataPath;
dataBase = struct([]);

for subj = 1:size(cfg,2)
    sub_label = cfg(subj).sub_labels;
    ses_label = cfg(subj).ses_label;
    task_label = cfg(subj).task_label;

    if isfield(cfg,'run_label')
        if size(cfg(subj).run_label{1},2)>4 % if label is more than run-
            run_label = cfg(subj).run_label{1};
        else
            run_label = {'run-*'};
        end
    else
        run_label = {'run-*'};
    end


    D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
        [sub_label '_' ses_label '_' task_label ,'_',run_label, '_ieeg.eeg']));

    if size(D,1) == 0
        error('%s does not exist',fullfile(dataPath,sub_label, ses_label,'ieeg',...
            [sub_label '_' ses_label '_' task_label ,'_',run_label, '_ieeg.eeg']))
    end

    % determine run_label
    if ~isfield(cfg,'run_label') || size(cfg(subj).run_label{1},2)<=4
        if size(D,1) == 1
            run_label = D(1).name(strfind(D(1).name,'run-'):strfind(D(1).name,'_ieeg')-1);
        else
            run_label = D(1).name(strfind(D(1).name,'run-'):strfind(D(1).name,'_ieeg')-1);
            fprintf('WARNING: More runs were available for %s_%s_%s, so determine run_label! \n',sub_label,ses_label,task_label)
        end

        dataName = fullfile(D(1).folder, D(1).name);
    else
        dataName = fullfile(D(1).folder, D(1).name);
    end

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

    dataBase(subj).sub_label = sub_label;
    dataBase(subj).ses_label = ses_label;
    dataBase(subj).task_label = task_label;
    dataBase(subj).run_label = run_label;
    dataBase(subj).dataName = dataName;
    dataBase(subj).ccep_header = ccep_header;
    dataBase(subj).tb_events = tb_events;
    dataBase(subj).tb_channels = tb_channels;
    dataBase(subj).tb_electrodes = tb_electrodes;
    dataBase(subj).ch = ch_incl;
    dataBase(subj).data = data;
    fprintf('...Subject %s has been run...\n',sub_label)
end

disp('All subjects are loaded')
