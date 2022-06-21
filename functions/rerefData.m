function dataBase = rerefData(dataBase,IED,subj)

nvisChan = setdiff(1:size(dataBase(subj).ch,1),[IED',... 
    find(contains(dataBase(subj).tb_channels.status,'bad')==1)']);
dataBase(subj).data_avg = mean(dataBase(subj).data(nvisChan,:),1);

dataBase(subj).data_reref = dataBase(subj).data - dataBase(subj).data_avg;

end


