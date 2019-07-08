% CCEP_powerSuppression_ERs_spikes
% this code analyzes the relationship between ERs, spikes and
% powerSuppression.
% author: D van Blooijs
% date: April 2019


% stappen:
% 1: eerste 15 minuten van de file pieken scoren in 1 kanaal (triggers in
%       ECoG Systemplus)
% 2: SPES van deze 10 patienten in BIDS zetten --> check! :)
% 3: power suppression van Michelle in derivatives BIDS zetten --> check :)
% 4: ERs in derivatives BIDS zetten --> check
% 5: per patient goede spike detectie maken adhv eerst gescoorde pieken in
%       deze patienten --> 10 min trainen, 5 min testen
% 6: spikes detecteren in rest van ECoG
% 7: latency tussen spikes bepalen voor en na stimulatie :)
% 8: statistical analysis
% 9: make figures
% 10: check SVM-scripts Michelle and make them BIDS compatible


%% sort ERs_BB to cc_stimsets

for subj = 1:size(dataBase,2)
   stimchans = reshape([dataBase(subj).ERs_BB(:).stimchans],2, size(dataBase(subj).ERs_BB,2))';
    if isequal(dataBase(subj).cc_stimsets,stimchans)
        for stimp = 1:size(dataBase(subj).cc_stimsets,1)
            
            dataBase(subj).cc_ERs(stimp).vis = dataBase(subj).ERs_BB(stimp).ERsvis;
            dataBase(subj).cc_BBs(stimp).det = dataBase(subj).ERs_BB(stimp).BB;
            
        end
        
    else
        for stimp = 1:size(dataBase(subj).cc_stimsets,1)
            stimnum = find(stimchans(:,1) == dataBase(subj).cc_stimsets(stimp,1)) & find(stimchans(:,2) == dataBase(subj).cc_stimsets(stimp,2));
            
            dataBase(subj).cc_ERs(stimp).vis = dataBase(subj).ERs_BB(stimnum).ERsvis;
            dataBase(subj).cc_BBs(stimp).det = dataBase(subj).ERs_BB(stimnum).BB;
            
        end
    end    
end


%% determine spike detection parameters and performance

par = setParCortSpikes; % nog maken

%% detect spikes

detCortSpikes(ECoG,par); % nog maken/aanpassen van findCortSpikes

%% make ratio spike change



%% categorize into : 0=no change, 1=decrease, 2=increase

%% analysis ERs_PS
% chi squared test
% [~,~,pchi] = crosstab(round(pmat),ERmat);

% figure stacked bar plot

%% analysis ERs_spikes
% chi squared test
% figure stacked bar plot

%% analysis PS_spikes
% chi squared test
% figure stacked bar plot
