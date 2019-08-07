# CCEP_suppressionPower_Spikes
This is the code used to analyse the effect of CCEP on suppression of power and epileptic spikes. 

1. Use preprocess_mainfile to load the data into matlab and epoch around stimuli (-2s:2s). This works very well with data in BIDS structure. 
2. Use mainfile_TFSPES to make ERSP-figures in whom power suppression (blue blob) can be observed. 
3. Use rate_visTFSPESscoring to load individual scorings (by DvB and MvdS). 
4. Use mainfileSVM to first make the SVM detector and then detect power suppressions in other patients. 
5. Use makeSpikeDet to detect epileptic spikes in the data. 
6. Use analysis_ERs_spikes to analyze the associations between power suppression, occurrence of CCEP, and spike suppression. 
