function optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold,subj)

if isempty(trainPar.ThU_opt)
    minL = min(min(min(train_threshold(subj).L)));
    
    [r,~] = find(train_threshold(subj).L == minL,1,'first'); %ThU
    [r2,c] = find(squeeze(train_threshold(subj).L(r,:,:)) == minL,1,'first'); % ThL and C
    
    ThUopt = trainPar.ThU(r) ;
    ThLopt = trainPar.ThL(r2) ;
    Copt = trainPar.C(c);
    Lopt = train_threshold(subj).L(r,r2,c);
    sensopt = train_threshold(subj).sens(r,r2,c);
    specopt = train_threshold(subj).spec(r,r2,c);
    precopt = train_threshold(subj).prec(r,r2,c);
    Fopt = train_threshold(subj).F(r,r2,c);
    
    if size(train_threshold,2) == 1
        fprintf('\nOptimal thresholds for all patients combined are ThU = %1.1f, ThL = %1.1f, C = %1.4f \n',ThUopt,ThLopt,Copt)
        fprintf('Resulting in a performance: L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',Lopt, sensopt, specopt,precopt,Fopt)
        
    elseif size(train_threshold,2) >1
        fprintf('\nOptimal thresholds for %s are ThU = %1.1f, ThL = %1.1f, C = %1.4f \n',dataBase(subj).sub_label,ThUopt,ThLopt,Copt)
        fprintf('Resulting in a performance: L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',Lopt, sensopt, specopt,precopt,Fopt)
        
        other_subj = setdiff(cfg.train,subj);
        for n = 1:size(other_subj,2)
            osubj = other_subj(n);
            dataBase(osubj).sub_label;
            Lopt = train_threshold(osubj).L(r,r2,c);
            sensopt = train_threshold(osubj).sens(r,r2,c);
            specopt = train_threshold(osubj).spec(r,r2,c);
            precopt = train_threshold(osubj).prec(r,r2,c);
            Fopt = train_threshold(osubj).F(r,r2,c);
            
            fprintf('\n These thresholds for %s result in a performance in %s of: \n L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',...
                dataBase(subj).sub_label, dataBase(osubj).sub_label,Lopt, sensopt, specopt,precopt,Fopt)
            
        end
        
    end
else
    r = find(trainPar.ThU == trainPar.ThU_opt);
    r2 = find(trainPar.ThL == trainPar.ThL_opt);
    c = find(trainPar.C == trainPar.C_opt);
    
    ThUopt = trainPar.ThU(r) ;
    ThLopt = trainPar.ThL(r2) ;
    Copt = trainPar.C(c);
    
    if size(train_threshold,2) == 1
        
        Lopt = train_threshold(1).L(r,r2,c);
        sensopt = train_threshold(1).sens(r,r2,c);
        specopt = train_threshold(1).spec(r,r2,c);
        precopt = train_threshold(1).prec(r,r2,c);
        Fopt = train_threshold(1).F(r,r2,c);

        fprintf('\nThe optional thresholds are ThU = %1.1f, ThL = %1.1f, C = %1.4f \n',ThUopt,ThLopt,Copt)
        fprintf('\n These thresholds result in a performance in all patients combined of: \n L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',...
            Lopt, sensopt, specopt,precopt,Fopt)

    elseif size(train_threshold,2) >1
        fprintf('\nThe optional thresholds are ThU = %1.1f, ThL = %1.1f, C = %1.4f \n',ThUopt,ThLopt,Copt)
        
        for n = 1:size(cfg.train,2)
            subj = cfg.train(n);
            Lopt = train_threshold(subj).L(r,r2,c);
            sensopt = train_threshold(subj).sens(r,r2,c);
            specopt = train_threshold(subj).spec(r,r2,c);
            precopt = train_threshold(subj).prec(r,r2,c);
            Fopt = train_threshold(subj).F(r,r2,c);
            
            fprintf('\n These thresholds result in a performance in %s of: \n L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',...
                dataBase(subj).sub_label, Lopt, sensopt, specopt,precopt,Fopt)
            
        end
        
    end
end



