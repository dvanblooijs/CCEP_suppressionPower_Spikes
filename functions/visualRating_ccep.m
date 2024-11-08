function dataBase = visualRating_ccep(myDataPath,dataBase,subj,cfg,endstimp)

filefolder = fullfile(myDataPath.proj_diroutput,dataBase(subj).sub_label);
if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

ccep = dataBase(subj).ccep;
fs = dataBase(subj).ccep_header.Fs;
tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

% if n1 peaks are not detected yet
if sum(strcmp(fieldnames(dataBase),'ccep'))==0
    ccep.n1_peak_sample = zeros(size(dataBase(subj).cc_epoch_sorted_avg,1),size(dataBase(subj).cc_epoch_sorted_avg,2));
elseif  sum(strcmp(fieldnames(dataBase(subj).ccep),'n1_peak_amplitude'))==0
    ccep.n1_peak_sample = zeros(size(dataBase(subj).cc_epoch_sorted_avg,1),size(dataBase(subj).cc_epoch_sorted_avg,2));
end

% if no n1 peaks are checked yet
if sum(strcmp(fieldnames(dataBase),'ccep'))==0
    ccep.checked = zeros(size(dataBase(subj).cc_epoch_sorted_avg,1),size(dataBase(subj).cc_epoch_sorted_avg,2));
elseif  sum(strcmp(fieldnames(dataBase(subj).ccep),'checked'))==0
    ccep.checked = zeros(size(dataBase(subj).cc_epoch_sorted_avg,1),size(dataBase(subj).cc_epoch_sorted_avg,2));
end

n=numel(ccep.checked(:,1:endstimp))+1;
for stimp = endstimp+1:size(dataBase(subj).cc_epoch_sorted_avg,2)
    for chan = 1 :size(dataBase(subj).cc_epoch_sorted_avg,1)
        
        if ~isnan(ccep.n1_peak_amplitude(chan,stimp))
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.11 0.77 0.8];
            
            low_ci = quantile(squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:)),.16,1); % find "mean - SD"
            high_ci = quantile(squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:)),.84,1); % find "mean + SD"
            
            tt0 = find(tt>= -0.1,1,'first');
            tt1 = find(tt>= 0.5,1,'first');
            
            subplot(1,2,1)
            fill([tt(tt0:tt1) tt(tt1:-1:tt0)],[low_ci(tt0:tt1) high_ci(tt1:-1:tt0)],'r','LineStyle','none','FaceAlpha',0.3)
            hold on
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:)),'r:');
            plot(tt,low_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
            plot(tt,high_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
            plot(tt(ccep.n1_peak_sample(chan,stimp)),ccep.n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s-%s',dataBase(subj).ch{chan},dataBase(subj).cc_stimchans{stimp,1},dataBase(subj).cc_stimchans{stimp,2}))
            
            subplot(1,2,2)
            fill([tt(tt0:tt1) tt(tt1:-1:tt0)],[low_ci(tt0:tt1) high_ci(tt1:-1:tt0)],'r','LineStyle','none','FaceAlpha',0.3)
            hold on
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted(chan,:,stimp,:)),'r:');
            plot(tt,low_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
            plot(tt,high_ci,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
            plot(tt,squeeze(dataBase(subj).cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
            plot(tt(ccep.n1_peak_sample(chan,stimp)),ccep.n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',6)
            hold off
            
            xlim([-0.1 0.5])
            ylim([-750 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
            
            perc = n / size(ccep.n1_peak_amplitude(:),1) *100;
            
            x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this an ER (y/n) and correct N1 detection (y/yn)? [y/yn/n]: ',perc,dataBase(subj).cc_stimchans{stimp,:},dataBase(subj).ch{chan}),'s');
            
            if strcmp(x,'yn') % so it is an ER, but N1 is not correctly detected
                currkey = 0;
                while ~strcmp(currkey,'c') % currkey == 0 %
                    w = waitforbuttonpress;
                    
                    if w == 0     %0 for mouse click, 1 for button press
                        cp = get(gca,'CurrentPoint');
                        % find sample number closest to the selected point
                        [~,sampnum] = min(abs(tt-cp(1,1)));
                        
                        % find nearby peak
                        [~,locs] = findpeaks(-1*squeeze(dataBase(subj).cc_epoch_sorted_avg(chan,stimp,...
                            sampnum-round(0.01*fs):sampnum+round(0.01*fs))),'NPeaks',1,'SortStr','ascend');
                        % find x-position of nearby peak
                        locsamp = sampnum-round(0.01*fs)+locs-1;
                        
                        hold on
                        plot(tt(locsamp),dataBase(subj).cc_epoch_sorted_avg(chan,stimp,locsamp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6); drawnow;
                        hold off
                        
                        fprintf('sample: %1.4f, amplitude: %2.1f, Press "c" when correct \n',tt(locsamp),dataBase(subj).cc_epoch_sorted_avg(chan,stimp,locsamp))
                        
                    elseif w==1
                        currkey = get(gcf,'CurrentCharacter');
                    end
                end
                
                ccep.n1_peak_amplitude(chan,stimp) = dataBase(subj).cc_epoch_sorted_avg(chan,stimp,locsamp);
                ccep.n1_peak_sample(chan,stimp) = locsamp;
                ccep.checked(chan,stimp) = 1;
                
            elseif strcmp(x,'y')
                ccep.checked(chan,stimp) = 1 ;
            else
                ccep.checked(chan,stimp) = 0 ;
            end
        end
        n=n+1;
    end
    
    % save also till which stimpair visual N1s are checked.
    ccep.checkUntilStimp = stimp;
    
    filename = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_N1sChecked.mat'];
    
    % save file during scoring in case of error
    save(fullfile(filefolder,filename),'-struct','ccep');
    
end

dataBase(subj).ccep = ccep;
fprintf('Visual rating of %s completed\n',dataBase(subj).sub_label)

end
