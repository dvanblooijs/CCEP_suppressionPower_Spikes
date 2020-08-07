function dataBase = findIEDchannels(dataBase)


for subj = 1:size(dataBase,2)
    switch dataBase(subj).sub_label
        case 'sub-RESP0401'
            IEDchan = {'soc4','soc5','soc6','sta3','sta4','stv5','stv6','stv7'};% de kanalen met spikes soc4-6, sta3,4 stv5-7
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0435' % not really convincing spikes...
            IEDchan = {'fh13','fh14','fh15','fh16','fh20','fh21','fh22','fh23',...
                'fh24','fh26','fh27','fh28','fh30','fh31','fl03','fl04','fl05','fl12','fl13'};% FH13-16, 20-24, 26-28,30,31, FL3-5,12,13
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0458'
            IEDchan = {'f08','f09','f14','f17','f18','f19'};% F8-9,14,17-19
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0478'
            IEDchan = {'ihh04','ihh05','ihl01','ihl02','c05','c06','c12','c20'};% IHH4,5, IHL1,2, C5,6,12,20
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0574'
            IEDchan = {'dv1','dv2','dv3','da1','da2','da3','f04','f05','hf10',...
                'hf11','hf12','hf13','hf14'};% Dv1-3, Da1-3, F4,5, HF10-14  
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0608'
            IEDchan = {'d1','d2','d3','c06','c07','c08','c15','c16','c21','c22'};  % D1-3, C6-8, 15, 16, 21, 22
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
            
        case 'sub-RESP0699'
            IEDchan = {'c13','c14','c20','c21','c22','c23','c28','c29'};  % C13, 14, 20-23, 28,29
            IEDch = NaN(size(IEDchan));
            for chan =1:size(IEDchan,2)
                IEDch(chan) = find(strcmpi(dataBase(subj).ch,IEDchan{chan})==1);
            end
    end
    
    dataBase(subj).IEDch = IEDch;
    dataBase(subj).IEDchan = IEDchan;
end