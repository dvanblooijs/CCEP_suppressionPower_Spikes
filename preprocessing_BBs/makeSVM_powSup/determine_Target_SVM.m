% this function uses the visual scores to make TargetValues (Y).
% 1 = power suppression was scored
% 0 = no power suppression was scored

function [Y_alltrain,Y_alltest, Y_conc, Y] = determine_Target_SVM(dataBase,cfg)

% pre-allocation
S = cell(size(cfg.sub_labels,2),1);
Y_conc = cell(size(cfg.sub_labels,2),1);

for subj=1:size(cfg.sub_labels,2)
    
    for stimpair = 1:size(dataBase(subj).ERSP.allERSPboot,1)
        for chan = 1:size(dataBase(subj).ERSP.allERSPboot,2)
            
            %Get stimulus pair in visual scoring
            if isequal(dataBase(subj).ERSP.ch,dataBase(subj).vis_scores.chan)
                if isequal(dataBase(subj).ERSP.cc_stimsets,dataBase(subj).vis_scores.stimpnum)
                    %  fprintf('%s: Channels and stimsets are equal\n',dataBase(subj).sub_label)
                    stimp = stimpair;
                    channum = chan;
                else
                    if stimpair == 1
                        fprintf('%s: Channels are equal, but stimsets are not, but the correct location was found \n',dataBase(subj).sub_label)
                    end
                    stimp = find(dataBase(subj).ERSP.cc_stimsets(stimpair,1)==dataBase(subj).vis_scores.stimpnum(:,1) ...
                        & dataBase(subj).ERSP.cc_stimsets(stimpair,2)== dataBase(subj).vis_scores.stimpnum(:,2));
                    channum = chan;
                end
            else
                if chan == 1 && stimpair == 1 
                    fprintf('%s: Channels are not equal, but the correct location was found \n',dataBase(subj).sub_label)
                end
                channum = find(strcmpi(dataBase(subj).ERSP.ch{chan},dataBase(subj).vis_scores.chan));
                stimp = find(strcmpi([dataBase(subj).ERSP.cc_stimchans{stimpair,1},'-',dataBase(subj).ERSP.cc_stimchans{stimpair,2}],dataBase(subj).vis_scores.stimorder));
            end
            
            %Get indexes of similar scoring
            if dataBase(subj).vis_scores.BS(1,stimp,channum)==dataBase(subj).vis_scores.BS(2,stimp,channum)
                S{subj}(stimpair,chan) = dataBase(subj).vis_scores.BS(1,stimp,channum);
            else
                S{subj}(stimpair,chan) = 0;
            end
        end
    end
    Y{subj} = S{subj};
    Y_conc{subj} = reshape(S{subj},numel(S{subj}),1);
end

Y_alltrain = vertcat(Y_conc{cfg.train});                           % store true label in Y
Y_alltest = vertcat(Y_conc{cfg.test});                           % store true label in Y

end