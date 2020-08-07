function display_optimalValSpikeDetector(dataBase,train_thresholds_all,train_thresh,train_t_dif)

for subj = 1:size(dataBase,2)
    clc
    for chanpat = 1:2
        
        maxF = max(max(max(train_thresholds_all.F(subj,:,:,chanpat))));
        
        [SD,t] = find(squeeze(train_thresholds_all.F(subj,:,:,chanpat)) == maxF,1,'first');
        
        if chanpat==1
            string = 'channel-specific';
        elseif chanpat ==2
            string = ' patient-specific';
        end
        
        sens_maxF = train_thresholds_all.sens(subj,SD,t,chanpat);
        prec_maxF = train_thresholds_all.prec(subj,SD,t,chanpat);
        
        fprintf('--- For %s, F = %0.2f, sens = %0.2f, prec = %0.2f threshSD = %g, t_dif = %g, %s --- \n',...
            dataBase(subj).sub_label,maxF,sens_maxF,prec_maxF,train_thresh(SD),train_t_dif(t),string)
        
        for subj2 = 1:size(dataBase,2)
            
            fprintf('For %s, F = %0.2f, sens = %0.2f, prec = %0.2f threshSD = %g, t_dif = %g, %s \n',...
                dataBase(subj2).sub_label,train_thresholds_all.F(subj2,SD,t,chanpat),...
                train_thresholds_all.sens(subj2,SD,t,chanpat),...
                train_thresholds_all.prec(subj2,SD,t,chanpat),train_thresh(SD),train_t_dif(t), string)
            
        end
        
        fprintf('For combined patients, F = %0.2f, sens = %0.2f, prec = %0.2f, threshSD = %g, t_dif = %g \n\n',...
            train_thresholds_all(1).Fmed(SD,t,chanpat),...
            train_thresholds_all(1).sensmed(SD,t,chanpat),...
            train_thresholds_all(1).precmed(SD,t,chanpat),train_thresh(SD),train_t_dif(t))
        
    end
    pause
end
clc

% max value for patients combined
for chanpat = 1:2
    maxF = max(max(train_thresholds_all(1).Fmed(:,:,chanpat)));
    
    [SD,t] = find(train_thresholds_all.Fmed(:,:,chanpat) == maxF,1,'first');
    
    sens_maxF = train_thresholds_all.sensmed(SD,t,chanpat);
    prec_maxF = train_thresholds_all.precmed(SD,t,chanpat);
    
    if chanpat==1
        string = 'channel-specific';
    elseif chanpat ==2
        string = ' patient-specific';
    end
    
    fprintf('\nFor combined patients, F = %0.2f, sens = %0.2f, prec = %0.2f, threshSD = %g, t_dif = %g, %s\n',...
        maxF,sens_maxF,prec_maxF,train_thresh(SD),train_t_dif(t),string)
    
    for subj2 = 1:size(dataBase,2)
        
        fprintf('For %s, F = %0.2f, sens = %0.2f, prec = %0.2f threshSD = %g, t_dif = %g, %s \n',...
            dataBase(subj2).sub_label,train_thresholds_all.F(subj2,SD,t,chanpat),...
            train_thresholds_all.sens(subj2,SD,t,chanpat),...
            train_thresholds_all.prec(subj2,SD,t,chanpat),train_thresh(SD),train_t_dif(t), string)
        
    end
end