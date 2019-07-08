%% script to generate TFSPES plots
% author: Michelle van der Stoel
% date: jan-sept 2018

clear all 
close all

%% Add needed folders
addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/eeglab14_1_1b'));       % with genpath add all subfolders 
addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/Data'));       
addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/SPES_GJ'));

%% Load data
load PAT_99_allstimp                                                       % ctrl f --> change patient number (replace all)

pat=99;
%% Choose stimulation from specific electrode
for n = 1:size(allstims,1)                                                  % number for stim pair
    for m = 1:size(allstims,2)                                               % number for recording electrode
        tmpsig=allstims{n,m};                                                   % specific stim pair + rec electrode
        tmpsig=tmpsig(fs:3*fs-1,:);                                             % get 2 seconds (1sec before stim,1 sec after)
        
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
        
        stimpair=stimpchan(n,:);
        stimName=deblank(char(channels(stimpair)));
        rec=deblank(char(channels(m)));       
        chlabel=sprintf('TF-SPES stim %s%s (%d %d) response %s (%d)',stimName(1,:),stimName(2,:),stimpair(1),stimpair(2),rec,m);
        
       [ERSP,~,powbase,times,freqs,erspboot,~]=gjnewtimef(tmpsig(:,:),length(pointrange),[tlimits(1) tlimits(2)],EEG.srate,cycles,'type','coher','title',chlabel, 'padratio',4,...
           'plotphase','off','freqs',frange,'plotitc','off','newfig','off','erspmax',15','alpha',0.05,'baseline',[-1000 -100]); % change alpha number to NaN to turn bootstrapping off
       
        figsize
        
        outlabel=sprintf('Stimpair%d-%dResponse%d.jpg',stimpair(1),stimpair(2),m);
%         output = '/Users/michelle/Google Drive/M3/LaTeX/images/';
%         print('-depsc','-r300',[output,outlabel])
%         
        % Create a name for a subfolder within output
        output = '/Users/michelle/Google Drive/M3/Onderzoek/Matlab/output/OUTPUT_PAT99/';
        newSubFolder = sprintf('%s/Stimpair%d-%d/', output,stimpair(1),stimpair(2));
        % Finally, create the folder if it doesn't exist already.
        if ~exist(newSubFolder, 'dir')
            mkdir(newSubFolder);
        end
        
        saveas(gcf,[newSubFolder,outlabel],'jpg')
        close(gcf);
        
        %%Get name channels in image name      
        %outlabel=sprintf('%s%s_%d-%d%s%d.jpg',stimName(1,:),stimName(2,:), stimpair(1),stimpair(2),rec,m);
        
        %option write ERSP matrix 
        %outlabel=sprintf('%schannel%s.mat',char(stimName), char(rec));
        %save(outlabel,'ERSP');
       
       
        %% Hysteresis
        % addpath(genpath('/Users/michelle/Google Drive/M3/Onderzoek/Matlab/hysteresis3d'))
        t_start = length(times)/2+1;                        % t_start after stim
        ERSP2 = ERSP * -1;                                  % Hysteresis takes highest threshold (not absolute)
        ERSP2 = ERSP2(:,t_start:end);                       % Get part after stim
        
        allERSP{n,m} = ERSP;
        allERSP2_99{n,m} = ERSP2; 
        
        clear hys tri ERSP ERSP2 stats
    end
end

targetFolder = output; 
fileName='Output_PAT_99';
save([targetFolder,fileName], 'allERSP', 'allERSP2_99','stimpchan');

% %Plot image manually
% a=allERSP{1,2};
% a=flipud(a);
% imagesc(a,[-15,15])
% colormap jet