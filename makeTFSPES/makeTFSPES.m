%% script to generate TFSPES plots
% author: Michelle van der Stoel
% date: Sep2017-Sep2018
% made BIDS compatible by: Dorien van Blooijs
% date: July 2019

function dataBase = makeTFSPES(dataBase,cfg)

for subj = 1:size(dataBase,2)
    fs = dataBase(subj).ccep_header.Fs;
    
    %% Choose stimulation from specific electrode
    for n=1:size(dataBase(subj).cc_epoch_sorted,3) % number for stim pair
        for m = 1:size(dataBase(subj).cc_epoch_sorted,1) % number for recording electrode

            tmpsig = squeeze(dataBase(subj).cc_epoch_sorted(n,:,m,:));
            tmpsig=tmpsig(:,fs:3*fs-1)';                                             % get 2 seconds (1sec before stim,1 sec after)
            
            EEG.pnts=size(tmpsig,1);                                                % Number of samples
            EEG.srate=fs;                                                           % sample frequency
            EEG.xmin=-1;                                                            % x axis limits
            EEG.xmax=1;
            tlimits = [EEG.xmin, EEG.xmax]*1000;                                    % time limit
            tmpsig = reshape(tmpsig, [1, size(tmpsig,1)*10]);                       % Get all 10 stim in 1 row
            
            %%Define parameters for EEGlab
            cycles=[3 0.8];
            frange = [10 250];                                                      % frequency range to plot/calculate spectrogram
            
            pointrange1 = round(max((tlimits(1)/1000-EEG.xmin)*EEG.srate, 1));
            pointrange2 = round(min((tlimits(2)/1000-EEG.xmin)*EEG.srate, EEG.pnts));
            pointrange = pointrange1:pointrange2;
            
            chlabel=sprintf('TF-SPES stim %s-%s response %s',...
                dataBase(subj).cc_stimchans{n,1},dataBase(subj).cc_stimchans{n,2},...
                dataBase(subj).ch{m});
           
            % Michelle used gjnewtimef but the results are equal so I choose to use newtimef
%             [ERSPgj,~,powbasegj,timesgj,freqsgj,erspbootgj,~]=gjnewtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cycles,'type','coher','title',chlabel, 'padratio',4,...
%                 'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15','alpha',0.05,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
            [ERSP,~,powbase,times,freqs,erspboot,~]=newtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cycles,'type','coher','title',chlabel, 'padratio',4,...
                'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15','alpha',0.05,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
            
            figsize
            
            outlabel=sprintf('Stimpair%s-%s_Response%s.jpg',...
                dataBase(subj).cc_stimchans{n,1},dataBase(subj).cc_stimchans{n,2},dataBase(subj).ch{m});

            if strcmp(cfg.saveERSP,'yes')
                % Create a name for a subfolder within output
                output = [cfg.output dataBase(subj).subj];
                newSubFolder = sprintf('%s/Stimpair%s-%s/', output,...
                    dataBase(subj).cc_stimchans{n,1},dataBase(subj).cc_stimchans{n,2});
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
            
            allERSP{n,m} = ERSP;
            allERSP2{n,m} = ERSP2;
            
            clear hys tri ERSP ERSP2 stats
        end
    end
    
    if strcmp(cfg.saveERSP,'yes')
        
        targetFolder = output;
        fileName='/ERSP';
        save([targetFolder,fileName], 'allERSP', 'allERSP2');
    end
    
    dataBase(subj).ERSP = allERSP2;
end


% %Plot image manually
% a=allERSP{1,2};
% a=flipud(a);
% imagesc(a,[-15,15])
% colormap jet