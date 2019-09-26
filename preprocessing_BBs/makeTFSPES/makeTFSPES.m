%% script to generate TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

function dataBase = makeTFSPES(dataBase,cfg)

for subj = 1:size(dataBase,2)
    fs = dataBase(subj).ccep_header.Fs;
    
    % pre-allocation
    allERSP = cell(size(dataBase(subj).cc_epoch_sorted,3),size(dataBase(subj).cc_epoch_sorted,1));
    allERSP2 = cell(size(dataBase(subj).cc_epoch_sorted,3),size(dataBase(subj).cc_epoch_sorted,1));
    allERSPboot = cell(size(dataBase(subj).cc_epoch_sorted,3),size(dataBase(subj).cc_epoch_sorted,1));
    allERSPboot2 = cell(size(dataBase(subj).cc_epoch_sorted,3),size(dataBase(subj).cc_epoch_sorted,1));
    
    %% Choose stimulation from specific electrode
    for stimp=1:size(dataBase(subj).cc_epoch_sorted,3) % number for stim pair
        for chan = 1:size(dataBase(subj).cc_epoch_sorted,1) % number for recording electrode
            
            tmpsig = squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:));
            tmpsig=tmpsig(:,round(fs):round(3*fs)-1);                                             % get 2 seconds (1sec before stim,1 sec after)
            
            if ~any(isnan(tmpsig)) % only run ERSP if there are 10 stimuli, no less
                EEG.pnts=size(tmpsig,2);                                                % Number of samples
                EEG.srate=fs;                                                           % sample frequency
                EEG.xmin=-1;                                                            % x axis limits
                EEG.xmax=1;
                tlimits = [EEG.xmin, EEG.xmax]*1000;                                    % time limit
                tmpsig = reshape(tmpsig', [1, size(tmpsig,2)*10]);                       % Get all 10 stim in 1 row
                
                %%Define parameters for EEGlab
                cycles = [3 0.8];
                frange = [10 250];                                                      % frequency range to plot/calculate spectrogram
                
                pointrange1 = round(max((tlimits(1)/1000-EEG.xmin)*EEG.srate, 1));
                pointrange2 = round(min((tlimits(2)/1000-EEG.xmin)*EEG.srate, EEG.pnts));
                pointrange = pointrange1:pointrange2;
                
                chlabel = sprintf('TF-SPES stim %s-%s response %s',...
                    dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2},...
                    dataBase(subj).ch{chan});
                
                % Michelle used gjnewtimef but the results are equal so I choose to use newtimef
                %             [ERSPgj,~,powbasegj,timesgj,freqsgj,erspbootgj,~]=gjnewtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cycles,'type','coher','title',chlabel, 'padratio',4,...
                %                 'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15','alpha',0.05,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
                [ERSP,~,~,times,freqs,erspboot,~]=newtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cycles,'type','coher','title',chlabel, 'padratio',4,...
                    'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15','alpha',0.05,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
                
                figsize
                
                outlabel=sprintf('Stimpair%s-%s_Response%s.jpg',...
                    dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2},dataBase(subj).ch{chan});
                
                if strcmp(cfg.saveERSP,'yes')
                    % Create a name for a subfolder within output
                    output = fullfile(cfg.output,dataBase(1).sub_label,dataBase(1).ses_label,...
                        dataBase(1).run_label);
                        
                    newSubFolder = sprintf('%s/Stimpair%s-%s/', output,...
                        dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2});
                    % Create the folder if it doesn't exist already.
                    if ~exist(newSubFolder, 'dir')
                        mkdir(newSubFolder);
                    end
                    
                    saveas(gcf,[newSubFolder,outlabel],'jpg')
                end
                
                close(gcf);
                
                %% Hysteresis
                t_start = length(times)/2+1;                        % t_start after stim
                ERSP2 = ERSP * -1;                                  % Hysteresis takes highest threshold (not absolute)
                ERSP2 = ERSP2(:,t_start:end);                       % Get part after stim
                
                ERSP_boot = ERSP;
                ERSP_boot(ERSP>erspboot(:,1) & ERSP<erspboot(:,2)) = 0;
                ERSP_boot2 = ERSP_boot * -1;                                  % Hysteresis takes highest threshold (not absolute)
                ERSP_boot2 = ERSP_boot2(:,t_start:end);                       % Get part after stim
                
                allERSP{stimp,chan} = ERSP;
                allERSPboot{stimp,chan} = ERSP_boot;
                allERSP2{stimp,chan} = ERSP2;
                allERSPboot2{stimp,chan} = ERSP_boot2;
                
                clear ERSP ERSP2
                
            else
                allERSP{stimp,chan} = [];
                allERSPboot{stimp,chan} = [];
                allERSP2{stimp,chan} = [];
                allERSPboot2{stimp,chan} = [];

            end
        end
        
        if strcmp(cfg.saveERSP,'yes')
            
            targetFolder = output;
            fileName=['/' dataBase(1).sub_label,'_' dataBase(1).ses_label,...
                        '_', dataBase(1).task_label,'_',dataBase(1).run_label '_ERSP.mat'];
            save([targetFolder,fileName], 'allERSP', 'allERSP2', 'allERSPboot','allERSPboot2','times','freqs');
        end
        
        dataBase(subj).ERSP = allERSP;
        dataBase(subj).ERSPboot = allERSPboot;
        
    end
end


%Plot image manually
% figure(1),
% a=allERSP{stimp,chan};
% a=flipud(a);
% imagesc(a,[-15,15])
% colormap jet