function optTrainPar_SVM(cfg,dataBase,trainPar,train_threshold,subj)

if isempty(trainPar.C_opt)
    minL = min(train_threshold(subj).L);
    
    r = find(train_threshold(subj).L == minL,1,'first'); 
    Copt = trainPar.C(r);
    Lopt = train_threshold(subj).L(r);
    sensopt = train_threshold(subj).sens(r);
    specopt = train_threshold(subj).spec(r);
    precopt = train_threshold(subj).prec(r);
    Fopt = train_threshold(subj).F(r);
    
    if size(train_threshold,2) == 1
        fprintf('\nOptimal thresholds for all patients combined are C = %1.4f \n',Copt)
        fprintf('Resulting in a performance: L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',Lopt, sensopt, specopt,precopt,Fopt)
        
    elseif size(train_threshold,2) >1
        fprintf('\nOptimal thresholds for %s are C = %1.4f \n',dataBase(subj).sub_label,Copt)
        fprintf('Resulting in a performance: L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',Lopt, sensopt, specopt,precopt,Fopt)
        
        other_subj = setdiff(cfg.train,subj);
        for n = 1:size(other_subj,2)
            osubj = other_subj(n);
            dataBase(osubj).sub_label;
            Lopt = train_threshold(osubj).L(r);
            sensopt = train_threshold(osubj).sens(r);
            specopt = train_threshold(osubj).spec(r);
            precopt = train_threshold(osubj).prec(r);
            Fopt = train_threshold(osubj).F(r);
            
            fprintf('\n These thresholds for %s result in a performance in %s of: \n L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',...
                dataBase(subj).sub_label, dataBase(osubj).sub_label,Lopt, sensopt, specopt,precopt,Fopt)
            
        end
        
    end
else
    c = find(trainPar.C == trainPar.C_opt);
    Copt = trainPar.C(c);
    
    if size(train_threshold,2) == 1
        
        Lopt = train_threshold(1).L(c);
        sensopt = train_threshold(1).sens(c);
        specopt = train_threshold(1).spec(c);
        precopt = train_threshold(1).prec(c);
        Fopt = train_threshold(1).F(c);

        fprintf('\nThe optional thresholds are C = %1.4f \n',Copt)
        fprintf('\n These thresholds result in a performance in all patients combined of: \n L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',...
            Lopt, sensopt, specopt,precopt,Fopt)

    elseif size(train_threshold,2) >1
        fprintf('\nThe optional thresholds are C = %1.4f \n',Copt)
        
        for n = 1:size(cfg.train,2)
            subj = cfg.train(n);
            Lopt = train_threshold(subj).L(c);
            sensopt = train_threshold(subj).sens(c);
            specopt = train_threshold(subj).spec(c);
            precopt = train_threshold(subj).prec(c);
            Fopt = train_threshold(subj).F(c);
            
            fprintf('\n These thresholds result in a performance in %s of: \n L = %1.2f, sensitivity = %1.2f, specificity = %1.2f, precision = %1.2f, F = %1.2f \n',...
                dataBase(subj).sub_label, Lopt, sensopt, specopt,precopt,Fopt)
            
        end
        
    end
end



