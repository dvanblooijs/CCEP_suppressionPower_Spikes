%% script to generate TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

% this function uses the epoched SPES-data and calculates the Event-Related
% Spectral Perturbation plot.

function dataBase = makeTFSPES(dataBase,myDataPath,cfg)

for nSubj = 1:size(dataBase,2)

    fs = dataBase(nSubj).ccep_header.Fs;
    
    % used rereferenced or not rereferenced data
    if cfg.reref == 1
        cc_epoch_sorted = dataBase(nSubj).cc_epoch_sorted_reref;
    else
        cc_epoch_sorted = dataBase(nSubj).cc_epoch_sorted;
    end

    % pre-allocation
    allERSP = cell(size(cc_epoch_sorted,3),size(cc_epoch_sorted,1));
    allERSPboot = cell(size(cc_epoch_sorted,3),size(cc_epoch_sorted,1));
    del_stimp = zeros(size(cc_epoch_sorted,3),1);
    
    %% Choose stimulation from specific electrode
    for nStimp=1:size(cc_epoch_sorted,3) % number for stim pair
        for nChan = 1:size(cc_epoch_sorted,1) % number for recording electrode
            
            tmpsig = squeeze(cc_epoch_sorted(nChan,:,nStimp,:));
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
                    dataBase(nSubj).cc_stimchans{nStimp,1},dataBase(nSubj).cc_stimchans{nStimp,2},...
                    dataBase(nSubj).ch{nChan});
                
                [ERSP,~,~,times,freqs,erspboot,~] = newtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cfg.cycles,'type','coher','title',chlabel, 'padratio',4,...
                    'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15,'alpha',cfg.alpha,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
                
                figsize
                
                outlabel = sprintf('Stimpair%s-%s_Response%s.jpg',...
                    dataBase(nSubj).cc_stimchans{nStimp,1},dataBase(nSubj).cc_stimchans{nStimp,2},dataBase(nSubj).ch{nChan});
                
                if strcmp(cfg.saveERSP,'yes')
                    % Create a name for a subfolder within output
                    output = fullfile(myDataPath.proj_diroutput,'ERSP',dataBase(nSubj).sub_label,dataBase(nSubj).ses_label,...
                        dataBase(nSubj).run_label);
                    
                    newSubFolder = sprintf('%s/Stimpair%s-%s/', output,...
                        dataBase(nSubj).cc_stimchans{nStimp,1},dataBase(nSubj).cc_stimchans{nStimp,2});
                    
                    % Create the folder if it doesn't exist already.
                    if ~exist(newSubFolder, 'dir')
                        mkdir(newSubFolder);
                    end
                    
                    saveas(gcf,[newSubFolder,outlabel],'jpg')
                end
                
                close(gcf);
                                
                ERSP_boot = ERSP;
                ERSP_boot(ERSP>erspboot(:,1) & ERSP<erspboot(:,2)) = 0;
                
                allERSP{nStimp,nChan} = ERSP;
                allERSPboot{nStimp,nChan} = ERSP_boot;
                
                clear ERSP
                
            else
                del_stimp(nStimp) = 1;
                
                allERSPboot{nStimp,nChan} = [];                
                allERSP{nStimp,nChan} = [];                
            end
        end
        
        if strcmp(cfg.saveERSP,'yes')
            
            allERSPboot = allERSPboot(~del_stimp,:);
            allERSP = allERSP(~del_stimp,:);
            
            targetFolder = output;
            fileName = [dataBase(nSubj).sub_label,'_' dataBase(nSubj).ses_label,...
                '_', dataBase(nSubj).task_label,'_',dataBase(nSubj).run_label '_ERSP.mat'];
            cc_stimchans = dataBase(nSubj).cc_stimchans(~del_stimp,:);
            cc_stimsets = dataBase(nSubj).cc_stimsets(~del_stimp,:);
            ch = dataBase(nSubj).ch;
            
            save(fullfile(targetFolder,fileName), '-v7.3','allERSP', 'allERSPboot','times','freqs','cc_stimchans','cc_stimsets','ch','cfg');
        end
        
        dataBase(nSubj).allERSP = allERSP;
        dataBase(nSubj).allERSPboot = allERSPboot;
        dataBase(nSubj).times = times;
        dataBase(nSubj).freqs = freqs;
    end
end
