%% script to generate TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

function dataBase = makeTFSPES(dataBase,myDataPath,cfg)

for subj = 1:size(dataBase,2)
    fs = dataBase(subj).ccep_header.Fs;
    
    % used rereferenced or not rereferenced data
    if cfg.reref == 1
        cc_epoch_sorted = dataBase(subj).cc_epoch_sorted_reref;
    else
        cc_epoch_sorted = dataBase(subj).cc_epoch_sorted;
    end

    % pre-allocation
    allERSP = cell(size(cc_epoch_sorted,3),size(cc_epoch_sorted,1));
    allERSPboot = cell(size(cc_epoch_sorted,3),size(cc_epoch_sorted,1));
    del_stimp = zeros(size(cc_epoch_sorted,3),1);
    
    %% Choose stimulation from specific electrode
    for stimp=1:size(cc_epoch_sorted,3) % number for stim pair
        for chan = 1:size(cc_epoch_sorted,1) % number for recording electrode
            
            tmpsig = squeeze(cc_epoch_sorted(chan,:,stimp,:));
            tmpsig=tmpsig(:,round(fs):round(3*fs)-1);                                             % get 2 seconds (1sec before stim,1 sec after)
            
            if ~any(isnan(tmpsig(:))) % only run ERSP if there are 10 stimuli, no less
                EEG.pnts=size(tmpsig,2);                                                % Number of samples
                EEG.srate=fs;                                                           % sample frequency
                EEG.xmin=-1;                                                            % x axis limits
                EEG.xmax=1;
                tlimits = [EEG.xmin, EEG.xmax]*1000;                                    % time limit
                tmpsig = reshape(tmpsig', [1, size(tmpsig,2)*10]);                       % Get all 10 stim in 1 row
                
                %%Define parameters for EEGlab
                cfg.cycles = [3 0.8];
                frange = [10 250];                                                      % frequency range to plot/calculate spectrogram
                cfg.alpha = 0.05;
                
                pointrange1 = round(max((tlimits(1)/1000-EEG.xmin)*EEG.srate, 1));
                pointrange2 = round(min((tlimits(2)/1000-EEG.xmin)*EEG.srate, EEG.pnts));
                pointrange = pointrange1:pointrange2;
                
                chlabel = sprintf('TF-SPES stim %s-%s response %s',...
                    dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2},...
                    dataBase(subj).ch{chan});
                
                % Michelle used gjnewtimef but the results are equal so I choose to use newtimef
                %             [ERSPgj,~,powbasegj,timesgj,freqsgj,erspbootgj,~]=gjnewtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cycles,'type','coher','title',chlabel, 'padratio',4,...
                %                 'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15','alpha',0.05,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
                [ERSP,~,~,times,freqs,erspboot,~]=newtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cfg.cycles,'type','coher','title',chlabel, 'padratio',4,...
                    'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15,'alpha',cfg.alpha,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
                
                figsize
                
                outlabel=sprintf('Stimpair%s-%s_Response%s.jpg',...
                    dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2},dataBase(subj).ch{chan});
                
                if strcmp(cfg.saveERSP,'yes')
                    % Create a name for a subfolder within output
                    output = fullfile(myDataPath.ERSPoutput,dataBase(subj).sub_label,dataBase(subj).ses_label,...
                        dataBase(subj).run_label);
                    
                    newSubFolder = sprintf('%s/Stimpair%s-%s/', output,...
                        dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2});
                    % Create the folder if it doesn't exist already.
                    if ~exist(newSubFolder, 'dir')
                        mkdir(newSubFolder);
                    end
                    
                    saveas(gcf,[newSubFolder,outlabel],'jpg')
                end
                
                close(gcf);
                                
                ERSP_boot = ERSP;
                ERSP_boot(ERSP>erspboot(:,1) & ERSP<erspboot(:,2)) = 0;
                
                allERSP{stimp,chan} = ERSP;
                allERSPboot{stimp,chan} = ERSP_boot;
                
                clear ERSP
                
            else
                del_stimp(stimp) = 1;
                
                allERSPboot{stimp,chan} = [];                
                allERSP{stimp,chan} = [];                
            end
        end
        
        if strcmp(cfg.saveERSP,'yes')
            
            allERSPboot = allERSPboot(~del_stimp,:);
            allERSP = allERSP(~del_stimp,:);
            
            targetFolder = output;
            fileName=['/' dataBase(subj).sub_label,'_' dataBase(subj).ses_label,...
                '_', dataBase(subj).task_label,'_',dataBase(subj).run_label '_ERSP.mat'];
            cc_stimchans = dataBase(subj).cc_stimchans(~del_stimp,:);
            cc_stimsets = dataBase(subj).cc_stimsets(~del_stimp,:);
            ch = dataBase(subj).ch;
            
            save([targetFolder,fileName], '-v7.3','allERSP', 'allERSPboot','times','freqs','cc_stimchans','cc_stimsets','ch','cfg');
        end
        
        dataBase(subj).allERSP = allERSP;
        dataBase(subj).allERSPboot = allERSPboot;
        dataBase(subj).times = times;
        dataBase(subj).freqs = freqs;
    end
end
