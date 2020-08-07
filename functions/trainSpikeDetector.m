function [sens, prec, F] = trainSpikeDetector(dataBase,subj,train_thresh,train_t_dif)

fs = dataBase(subj).ccep_header.Fs;

samp_train = round([1, fs*5*60]);

% include Mvalues
Mvalues = dataBase(subj).spikes.norm_M;
Mtrain = Mvalues(:,samp_train(1):samp_train(2));
Pharmat = dataBase(subj).spikes.Pharmat_norm;
Pharmatall = dataBase(subj).spikes.Pharmatall_norm;

data_stimfree = dataBase(subj).data_rerefnoStimArt(dataBase(subj).visIEDs.IED,:);

% pre-allocation
sens = NaN(size(train_thresh,2),size(train_t_dif,2),2); %[SD thresh] x [t_dif] x [threshold per channel/patient]
prec = NaN(size(train_thresh,2),size(train_t_dif,2),2);
F = NaN(size(train_thresh,2),size(train_t_dif,2),2);
TP = NaN(size(dataBase(subj).visIEDs.IEDchan,2));
FP = NaN(size(dataBase(subj).visIEDs.IEDchan,2));
FN = NaN(size(dataBase(subj).visIEDs.IEDchan,2));

% vary threshold (X*SD)
for thresh_SD = 1:size(train_thresh,2)
    for thresh_t = 1:size(train_t_dif,2)
        dif = round(train_t_dif(thresh_t) * fs);
        
        for chan = 1:size(Mvalues,1)
            
            % include visually scored IEDs
            visIEDs = round(dataBase(subj).visIEDs.x_all(chan,:)*fs);
            visIEDs = visIEDs(visIEDs~=0);
            visIEDs = visIEDs(visIEDs<samp_train(end));
            
            % first threshold is per channel, second threshold is in
            % total for this patient, so with all channels used
            Mthresh_both = [(Pharmat(chan,3) + train_thresh(thresh_SD) * Pharmat(chan,2));...
                (Pharmatall(3) + train_thresh(thresh_SD) * Pharmatall(2))];
            
            for n=1:2 % for threshold per channel and for threshold determined in all channels at same time
                Mthresh = Mthresh_both(n);
                
                % per kanaal, locatie van M's en de waardes van M's die binnen grenzen zitten
                locs = find(Mtrain(chan,:)>Mthresh);
                
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
                    
                    FNchan = NaN(size(visIEDs));
                    for m = 1:size(visIEDs,2)
                        if ~ismembertol(visIEDs(m),sampnum_maxuV(1,:),round(0.05*fs),'DataScale',1)
                            FNchan(m) = 1;
                        end
                    end
                    
                    if sum(~isnan(FNchan)) + sum(~isnan(TPchan)) ~= size(visIEDs,2)
                        error('FN + TP does not add up to number of visually scored spikes')
                    end
                    
                    TP(chan,n) = sum(TPchan(~isnan(TPchan)));
                    FP(chan,n) = sum(FPchan(~isnan(FPchan)));
                    FN(chan,n) = sum(FNchan(~isnan(FNchan)));
                else
                    
                end
                
                sens(thresh_SD,thresh_t,n) = sum(TP(:,n))/(sum(TP(:,n)) + sum(FN(:,n)));
                prec(thresh_SD,thresh_t,n) = sum(TP(:,n))/(sum(TP(:,n)) + sum(FP(:,n)));
                F(thresh_SD,thresh_t,n) = 2*((sens(thresh_SD,thresh_t,n)*prec(thresh_SD,thresh_t,n))/(sens(thresh_SD,thresh_t,n)+prec(thresh_SD,thresh_t,n)));
                
                fprintf('--- subj = %s, Thresh = %g, T = %g, n = %g, chan = %s ---\n',dataBase(subj).sub_label,train_thresh(thresh_SD),train_t_dif(thresh_t),n,dataBase(subj).visIEDs.IEDchan{chan})
            end
            
        end
    end
end


