function dataBase = unique_stimpairs(dataBase,setting)

for subj = 1:size(dataBase,2)
    
stimpair = dataBase(subj).tb_events.electrical_stimulation_site(strcmp(dataBase(subj).tb_events.sub_type,'SPES'));

stimnum = NaN(size(stimpair,1),2);
for stimp = 1:size(stimpair,1)
    stimchans = strsplit(stimpair{stimp},'-');
    for chan = 1:2
        stimnum(stimp,chan) = find(strcmp(stimchans{chan},dataBase(subj).ch)==1);
    end
end

stimcur = str2double(dataBase(subj).tb_events.electrical_stimulation_current(strcmp(dataBase(subj).tb_events.sub_type,'SPES')));

if strcmp(setting.dir,'yes') && strcmp(setting.amp,'yes')
    stimelek = [stimnum stimcur];
elseif strcmp(setting.dir,'yes') && strcmp(setting.amp,'no')
    stimelek = stimnum;
elseif strcmp(setting.dir,'no') && strcmp(setting.amp,'yes')
    stimelek = [sort(stimnum,2) stimcur];
elseif strcmp(setting.dir,'no') && strcmp(setting.amp,'no')
    stimelek = sort(stimnum,2);
end

[cc_stimsets,~,IC] = unique(stimelek,'rows');

n = histcounts(IC,'BinMethod','integers');

if any(diff(n) ~= 0)
    stimremove = find(n<5); % remove al stimulation pairs that are stimulated less than 5 times
    
    stimelek(IC==stimremove,:) = [];
    
    [cc_stimsets,~,IC] = unique(stimelek,'rows');
    n = histcounts(IC,'BinMethod','integers');
    if any(diff(n) ~= 0)
        fprintf('ERROR: %s some stimulation pairs are stimulated less/more than all others\n',dataBase(subj).sub_label)
    end
    
end

cc_stimchans = cell(size(cc_stimsets,1),2);

for stimp = 1:size(cc_stimsets,1)
   for chan =1:2
       cc_stimchans{stimp,chan} = dataBase(subj).ch{cc_stimsets(stimp,chan)};
   end
    
end

max_stim = median(n);

dataBase(subj).cc_stimsets = cc_stimsets;
dataBase(subj).cc_stimchans = cc_stimchans;
dataBase(subj).max_stim = max_stim;

stimdif = find(n ~= max_stim);
for stimp =1:size(stimdif,2)
    [cc_stimchans{stimdif(stimp),1} '-' cc_stimchans{stimdif(stimp),2}]
end

end