
function [sens,prec,F]=testSpikeDetector(dataBase,subj,thresh_SD_opt,thresh_t_dif_opt)
clc

% pre-allocation
TP = NaN(size(dataBase(subj).visIEDs.IED)); FN = NaN(size(dataBase(subj).visIEDs.IED)); FP = NaN(size(dataBase(subj).visIEDs.IED));

fs = dataBase(subj).ccep_header.Fs;

samp_test = round([5*60*fs+1, fs*10*60]);

% include Mvalues
Mvalues = dataBase(subj).spikes.norm_M;
Mtest = Mvalues(:,round(samp_test(1):samp_test(2)));
Pharmat = dataBase(subj).spikes.Pharmat_norm;

data_stimfree = dataBase(subj).data_rerefnoStimArt(dataBase(subj).visIEDs.IED,:);

% vary threshold (X*SD)
dif = round(thresh_t_dif_opt * fs);

for chan = 1:size(Mvalues,1)
    % include visually scored IEDs
    visIEDs = round(dataBase(subj).visIEDs.x_all(chan,:)*fs);
    visIEDs = visIEDs(visIEDs~=0);
    visIEDs = unique(visIEDs(visIEDs>samp_test(1) & visIEDs<samp_test(end)) - samp_test(1));
    
    Mthresh = (Pharmat(chan,3) + thresh_SD_opt*Pharmat(chan,2));
    
    % per kanaal, locatie van M's en de waardes van M's die binnen grenzen zitten
    locs = find(Mtest(chan,:)>Mthresh);
    
    if ~isempty(locs)
        % dit is iedere keer de laatste sample van een groepje
        lastsampgroups = [locs(diff(locs)>dif) locs(end)];
        % en een eerste sample van een groepje
        firstsampgroups = [locs(1) locs(find(diff(locs)>dif)+1)];
        
        sampnums = cell(1,size(lastsampgroups,2)); uV_sampnums = cell(1,size(lastsampgroups,2));
        sampnum_maxuV = NaN(1,size(lastsampgroups,2)); maxuV = NaN(1,size(lastsampgroups,2));
        for i=1:size(lastsampgroups,2)
            % samplenummers met Mrange
            sampnums{1,i} = firstsampgroups(i):lastsampgroups(i);
            % uV values bij samplenummers
            uV_sampnums{1,i} = data_stimfree(chan,firstsampgroups(i):lastsampgroups(i));
            % location with max uV value
            sampnum_maxuV(1,i) = sampnums{1,i}(find(uV_sampnums{1,i} == max(uV_sampnums{1,i}),1,'first'));
            % max uV value
            maxuV(1,i) = max(uV_sampnums{1,i});
            
            % remove all saturated results (these are mainly stimulation effect)
            % and spikes with an amplitude < 0 uV
            if maxuV(1,i) >2900
                sampnum_maxuV(1,i) = NaN;
                maxuV(1,i)  = NaN;
            end
        end
        
        TPchan = NaN(size(sampnum_maxuV)); FPchan = NaN(size(sampnum_maxuV));
        for m = 1:size(sampnum_maxuV,2)
            if ismembertol(sampnum_maxuV(1,m),visIEDs,round(0.05*fs),'DataScale',1)
                TPchan(m) = 1;
            elseif ~ismembertol(sampnum_maxuV(1,m),visIEDs,round(0.05*fs),'DataScale',1)
                FPchan(m) = 1;
            end
        end
        
        FNchan = NaN(size(visIEDs)); TPcheck = NaN(size(visIEDs));
        for m = 1:size(visIEDs,2)
            if ~ismembertol(visIEDs(m),sampnum_maxuV(1,:),round(0.05*fs),'DataScale',1)
                FNchan(m) = 1;
            elseif ismembertol(visIEDs(m),sampnum_maxuV(1,:),round(0.05*fs),'DataScale',1)
                TPcheck(m) = 1;
            end
        end
        
        if sum(~isnan(FNchan)) + sum(~isnan(TPchan)) ~= size(visIEDs,2)
            error('FN + TP does not add up to number of visually scored spikes')
        end
        
        TP(chan) = sum(TPchan(~isnan(TPchan)));
        FP(chan) = sum(FPchan(~isnan(FPchan)));
        FN(chan) = sum(FNchan(~isnan(FNchan)));
    else
        
    end
    
end

sens = sum(TP)/(sum(TP) + sum(FN));
prec = sum(TP)/(sum(TP) + sum(FP));
F = 2*((sens*prec)/(sens+prec));

disp('All performance is calculated')

    
end