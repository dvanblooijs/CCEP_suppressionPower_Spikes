function dataBase = findIEDchannels(dataBase)

for nSubj = 1:size(dataBase,2)

    if any(contains(fieldnames(dataBase(nSubj).tb_electrodes),'ied'))
        idx_ied = contains(dataBase(nSubj).tb_electrodes.ied,'yes');
        
    else
        
        idx_ied = false(size(dataBase(nSubj).tb_electrodes,1),1);
    end

    IEDch = find(idx_ied ==1);
    IEDchan = dataBase(nSubj).tb_electrodes.name(idx_ied);

    dataBase(nSubj).IEDch = IEDch;
    dataBase(nSubj).IEDchan = IEDchan;
end