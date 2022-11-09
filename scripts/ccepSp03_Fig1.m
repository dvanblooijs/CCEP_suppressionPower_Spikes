% Figure 1: example of ten responses to spes stimulation, an CCEP, and an
% ERSP when there is a large change in spike ratio

%% first run ccepSp03_analysis_ERs_PS_spikes.m

%% plot a figure for each subject showing the spikes before and after stimulation
% in ten epochs, in whom power suppression and a CCEP was observed

close all

stimart_spikes = t>-0.1 & t < 0.1; % 100 ms around stimulus artefact should be removed, because a CCEP 100ms after stimulus could be wrongly detected as spike
% stimart_CCEPs = t>=0 & t < 0.01; % 10ms after stimulus artefact should be removed, because no signal is recorded, only stimulus artefact

for subj = 1:size(dataBase,2)

    if any(contains(fieldnames(dataBase(subj).tb_electrodes),'ied')) % patient with observed IEDs
        ERSPmat = dataBase(subj).ERSPmat;
        CCEPmat = dataBase(subj).CCEPmat;

        % selection of channels in which IEDs were observed
        ERSPmatIED = ERSPmat(:,dataBase(subj).IEDs.IEDch);
        CCEPmatIED = CCEPmat(:,dataBase(subj).IEDs.IEDch);

        IEDmat = dataBase(subj).IEDs.spikesratio;

        bothCCEP_ERSP = find(CCEPmatIED == 1 & ERSPmatIED == 1); % responses with both CCEP and ERSP

        [~,I] = sort(log(IEDmat(bothCCEP_ERSP)),'ascend');

        selResp = bothCCEP_ERSP(I(3)); % select the one stim-resp combination with the most negative spike ratio (spike decrease)

        [r,c] = ind2sub(size(IEDmat),selResp);

        selChan = dataBase(subj).IEDs.IEDch(c);
        selStim = r;
        yminSpikes = -2000;
        ymaxSpikes = 16000;

        f = figure(subj);
        fill([0 10 10 0],[yminSpikes yminSpikes ymaxSpikes ymaxSpikes],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for counting spikes
%         fill([0 0.01 0.01 0],[yminSpikes yminSpikes ymaxSpikes ymaxSpikes],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for counting spikes
        hold on

%         plot([0.1 0.1],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) % epoch CCEP
%         plot([-0.1 -0.1],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) % epoch spikes
        plot([100 100],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) % epoch CCEP
        plot([-100 -100],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) % epoch spikes
        plot([1000 1000],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) % epoch CCEP
        plot([-1000 -1000],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) % epoch spikes
%         plot([dataBase(subj).ERSP.times(1)/1000 dataBase(subj).ERSP.times(1)/1000],[yminSpikes ymaxSpikes], 'Color',[220/256 220/256 255/256]) % epoch ERSP
%         plot([dataBase(subj).ERSP.times(end)/1000 dataBase(subj).ERSP.times(end)/1000],[yminSpikes ymaxSpikes], 'Color',[220/256 220/256 255/256]) % epoch ERSP
        plot([dataBase(subj).ERSP.times(1) dataBase(subj).ERSP.times(1)],[yminSpikes ymaxSpikes], 'Color',[220/256 220/256 255/256]) % epoch ERSP
        plot([dataBase(subj).ERSP.times(end) dataBase(subj).ERSP.times(end)],[yminSpikes ymaxSpikes], 'Color',[220/256 220/256 255/256]) % epoch ERSP
        
        for n = 1:10

            signal_orig = -1*squeeze(dataBase(subj).cc_epoch_sorted(selChan,n,selStim,:))+n*1400;
            signal_spikes = signal_orig;
            signal_spikes(stimart_spikes) = NaN;
%             signal_nonart = signal_orig;
%             signal_nonart(stimart_CCEPs) = NaN;
         hold on

            scatter(t(dataBase(subj).IEDs.totEpochspikessamp{r,c,n})*1000,...
                signal_spikes(dataBase(subj).IEDs.totEpochspikessamp{r,c,n}),'o',...
                'MarkerEdgeColor','k','MarkerFaceColor',[220/256 220/256 255/256])
            plot(t*1000,signal_orig,'b','Color','k') % original signal
%             plot(t,signal_orig,'b','Color',[220/256 220/256 255/256]) % original signal
%             plot(t,signal_nonart,'b','Color',[45/256 45/256 255/256]) % signal without artefact
%             plot(t,signal_spikes,'b','Color',[155/256 155/256 255/256]) % signal for spikes


        end
        title(sprintf('%s: %s-%s, %s', ...
            dataBase(subj).sub_label,dataBase(subj).cc_stimchans{selStim,:},dataBase(subj).ch{selChan}))
        xlabel('Time (ms)')
        ylim([yminSpikes ymaxSpikes])
        f.Units = 'normalized';
        f.Position = [0.35 0.5 0.5 0.4];

        % plot CCEP
        signal_avg_orig = squeeze(dataBase(subj).cc_epoch_sorted_avg(selChan,selStim,:));
%         signal_avg_orig2 = mean(squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)));
%         signal_sterr = std(squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)))./sqrt(10);
%         lower_err = -1*signal_avg_orig2 - signal_sterr;
%         upper_err = -1*signal_avg_orig2 + signal_sterr;
%         signal_avg = signal_avg_orig;
%         signal_avg(stimart_CCEPs) = NaN;

        figure(subj+10)
        fill([0 10 10 0],[yminSpikes yminSpikes ymaxSpikes ymaxSpikes],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for counting spikes
        hold on
%         fill([t*1000 t(end:-1:1)*1000],[upper_err lower_err(end:-1:1)],[150/256 150/256 255/256],'EdgeColor',[1 1 1])
        %         plot(t,-1*squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)),'Color',[150/256 150/256 255/256])
        plot(t*1000, -1*squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)),':','Color','k','LineWidth',0.5)
        plot(t*1000, -1*signal_avg_orig,'Color','k','LineWidth',2)
%         plot(t*1000, -1*signal_avg_orig,'Color',[220/256 220/256 255/256],'LineWidth',2)
 %         plot(t*1000, -1*signal_avg,'Color',[45/256 45/256 255/256],'LineWidth',2)

        title(sprintf('%s: Averaged CCEP: %s-%s, %s',...
            dataBase(subj).sub_label,dataBase(subj).cc_stimchans{selStim,:},dataBase(subj).ch{selChan}))
        ylim([-1500 3500])
        xlim([-100 100])
        xlabel('Time (ms)')

        n = subj+20;
        ERSPall = dataBase(subj).ERSP;
        [~, axes1] = plot_ERSP(ERSPall,selStim,selChan,n);
        colorbar(axes1)

        figure(subj+30)
        minVal = -1; maxVal = 1;
        randNumbers = minVal+(maxVal-minVal)*rand(size(IEDmat));
        scatter(randNumbers(:),log(IEDmat(:)),20,[220/256 220/256 255/256],'filled')
        hold on
        scatter(randNumbers(selStim,c),log(IEDmat(selStim,c)),30,'k','filled')
        hold off
        xlim([-5 5])
        title(sprintf('Spike ratio of %s',dataBase(subj).sub_label))
        ylabel('Logarithmic spike ratio')
    end

end

% housekeeping
clear axes1 bothCCEP_ERSP c CCEPmat CCEPmatIED ERSPall ERSPmat ERSPmatIED
clear fig1 I IEDmat lower_err n r selChan selResp selStim signal_avg signal_avg_orig
clear signal_avg_orig2 signal_nonart signal_orig signal_spikes signal_sterr
clear stimart_CCEPs stimart_spikes subj upper_err

%% save figures for manuscript

subj = 8;
print(figure(subj),'-vector','-depsc',fullfile(myDataPath.Figures,['fig1_IED_' dataBase(subj).sub_label]))
print(figure(subj+10),'-vector','-depsc',fullfile(myDataPath.Figures,['fig1_CCEP_' dataBase(subj).sub_label]))
print(figure(subj+20),'-vector','-depsc',fullfile(myDataPath.Figures,['fig1_ERSP_' dataBase(subj).sub_label]))
print(figure(subj+30),'-vector','-depsc',fullfile(myDataPath.Figures,['fig1_logSR_' dataBase(subj).sub_label]))

fprintf('Figures for subject %d are saved. \n',subj)

%housekeeping
clear subj

% %% figure displaying in what time window what is analyzed
% 
% figure,
% fill([-2 -0.1 -0.1 -2],[0 0 1 1],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for counting spikes
% hold on
% fill([0.1 2 2 0.1],[0 0 1 1],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for counting spikes
% fill([dataBase(1).ERSP.times(1)/1000 dataBase(1).ERSP.times(end)/1000 dataBase(1).ERSP.times(end)/1000 dataBase(1).ERSP.times(1)/1000], ...
%     [1 1 2 2],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window ERSP
% fill([0.01 0.1 0.1 0.01],[2 2 3 3],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for CCEPs
% plot([0 0],[-1 4],'k')
% hold off
% ylim([-1 4])

%% end of script