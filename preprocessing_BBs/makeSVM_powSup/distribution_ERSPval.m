% this function determines the distribution of ERSPvalues and makes a
% histogram.
% PARAMETERS:
% - cfg.train                   : which subjects are in the train dataset
% - dataBase.sub_label          : names of each subject {name}
% - dataBase.ERSP.allERSPboot   : ERSP-arrays [timesxfreqs] for each stimuluspair-electrode combination
% - dataBase.ERSP.cc_stimsets   : containt the stimsets in this dataset [stimpairsx2]


function distribution_ERSPval(cfg,dataBase)

distribution = struct;
for subj = cfg.train
    nBB = []; %NaN(sum(Y{subj}(:)==0),size(dataBase(subj).ERSP.ch,1)-2,size(dataBase(subj).ERSP.cc_stimchans,1)); % FIX: pre-allocation (and line 86,91)
    BB = []; %NaN(sum(Y{subj}(:)==1),size(dataBase(subj).ERSP.ch,1)-2,size(dataBase(subj).ERSP.cc_stimchans,1));
    
    allERSP = dataBase(subj).ERSP.allERSPboot;
    max_ERSP_BB = NaN(size(allERSP));
    max_ERSP_nBB = NaN(size(allERSP));
    for stimpair=1:size(allERSP,1)                            % for each stimulation pair
        for chan=1:size(allERSP,2)                        % for each recording electrode
            
            % hysteresis3d assumes non-negative image. We want to detect "blue
            % blobs" (negative values). Therefor, I multiply the ERSP with -1
            % and remove all values <0.
            ERSP = allERSP{stimpair,chan};
            ERSP2 = -1* ERSP;
            ERSP2(ERSP2<0) = 0;
            
            if Y{subj}(stimpair,chan) == 0 && ~ismember(chan,dataBase(subj).ERSP.cc_stimsets(stimpair,:)) % no visual scored BB and no member of stimpair
                nBB = [nBB; ERSP2(:)];
                max_ERSP_nBB(stimpair,chan) = max(ERSP2(:));
                max_ERSP_BB(stimpair,chan) = NaN;
                
            elseif Y{subj}(stimpair,chan) == 1 && ~ismember(chan,dataBase(subj).ERSP.cc_stimsets(stimpair,:)) % visual scored BB and no member of stimpair
                BB = [BB; ERSP2(:)];
                max_ERSP_BB(stimpair,chan) = max(ERSP2(:));
                max_ERSP_nBB(stimpair,chan) = NaN;
            end
        end
    end
    
    distribution(subj).nBB = nBB;
    distribution(subj).BB = BB;
    distribution(subj).max_ERSP_BB = max_ERSP_BB(:);
    distribution(subj).max_ERSP_nBB = max_ERSP_nBB(:);
    fprintf('%s is completed\n',dataBase(subj).sub_label)
end

% make histograms of ERSP values for each individual patient

figure,

for subj = cfg.train
    
    nBB = distribution(subj).nBB;
    BB = distribution(subj).BB;
    
    n = histcounts([nBB;BB],'BinMethod','integers');
    n_sort = sort(n,'descend');
    ymax = round(n_sort(2),1,'significant');
    ymin = 0;
    
    subplot(3,1,subj)
    histogram(nBB,'BinMethod','integers','FaceColor','g')
    hold on
    histogram(BB,'BinMethod','integers','FaceColor','b')
    hold off
    
    ylim([ymin,ymax])
    xlabel('ERSP value')
    ylabel('Counts')
    title(sprintf('Patient %s',dataBase(subj).sub_label))
    legend('nBB','BB')
end

end