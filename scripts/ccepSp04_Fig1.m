% Figure 1: 
% B) example of ten responses to spes stimulation, 
% C) a CCEP
% D) an ERSP when there is a large change in spike ratio
% E) the spike ratio of this subject

%% first run ccepSp03_analysis_ERs_PS_spikes.m

%% plot a figure for each subject showing the spikes before and after stimulation
% in ten epochs, in whom power suppression and a CCEP was observed

close all

% 100 ms around stimulus artefact should be removed, because a CCEP 100ms after stimulus could be wrongly detected as spike
stimart_spikes = t>-0.1 & t < 0.1; 

for nSubj = 1:size(dataBase,2)

    if any(contains(fieldnames(dataBase(nSubj).tb_electrodes),'ied')) % patient with observed IEDs
        ERSPmat = dataBase(nSubj).ERSPmat;
        CCEPmat = dataBase(nSubj).CCEPmat;

        % selection ERSPmat/CCEPmat of channels in which IEDs were observed
        ERSPmatIED = ERSPmat(:,dataBase(nSubj).IEDs.IEDch);
        CCEPmatIED = CCEPmat(:,dataBase(nSubj).IEDs.IEDch);

        IEDmat = dataBase(nSubj).IEDs.spikesratio;

        % responses with both CCEP and ERSP
        bothCCEP_ERSP = find(CCEPmatIED == 1 & ERSPmatIED == 1); 

        % sort spikeratios with both CCEP and ERSP
        [~,I] = sort(log(IEDmat(bothCCEP_ERSP)),'ascend');

        selResp = bothCCEP_ERSP(I(3)); % select the one stim-resp combination with one of the most negative spike ratio (spike decrease)

        [row,column] = ind2sub(size(IEDmat),selResp);

        selectChan = dataBase(nSubj).IEDs.IEDch(column);
        selectStim = row;
        yminSpikes = -2000;
        ymaxSpikes = 16000;

        % FIGURE of spikes in all ten SPES-epochs
        f = figure(nSubj);
        % window of stim artefact
        fill([0 10 10 0],[yminSpikes yminSpikes ymaxSpikes ymaxSpikes],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) 
        hold on

        % epoch CCEP
        plot([100 100],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) 
        plot([1000 1000],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) 
        % epoch IEDs
        plot([-100 -100],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) 
        plot([-1000 -1000],[yminSpikes ymaxSpikes],'Color',[220/256 220/256 255/256]) 
        % epoch ERSP
        plot([dataBase(nSubj).ERSP.times(1) dataBase(nSubj).ERSP.times(1)],[yminSpikes ymaxSpikes], 'Color',[220/256 220/256 255/256]) % epoch ERSP
        plot([dataBase(nSubj).ERSP.times(end) dataBase(nSubj).ERSP.times(end)],[yminSpikes ymaxSpikes], 'Color',[220/256 220/256 255/256]) % epoch ERSP
        
        % plot each of the ten SPES-epochs
        for nSignal = 1:10

            signal_orig = -1*squeeze(dataBase(nSubj).cc_epoch_sorted(selectChan,nSignal,selectStim,:))+nSignal*1400;
            signal_spikes = signal_orig;
            signal_spikes(stimart_spikes) = NaN;
         hold on

         % plot the detected IEDs
            scatter(t(dataBase(nSubj).IEDs.totEpochspikessamp{row,column,nSignal})*1000,...
                signal_spikes(dataBase(nSubj).IEDs.totEpochspikessamp{row,column,nSignal}),'o',...
                'MarkerEdgeColor','k','MarkerFaceColor',[220/256 220/256 255/256])
            plot(t*1000,signal_orig,'b','Color','k') % original signal
        end

        title(sprintf('%s: %s-%s, %s', ...
            dataBase(nSubj).sub_label,dataBase(nSubj).cc_stimchans{selectStim,:},dataBase(nSubj).ch{selectChan}))
        xlabel('Time (ms)')
        ylim([yminSpikes ymaxSpikes])
        f.Units = 'normalized';
        f.Position = [0.35 0.5 0.5 0.4];

        % FIGURE CCEP
        signal_avg_orig = squeeze(dataBase(nSubj).cc_epoch_sorted_avg(selectChan,selectStim,:));

        figure(nSubj+10)
        % window of stim artefact
        fill([0 10 10 0],[yminSpikes yminSpikes ymaxSpikes ymaxSpikes],[220/256 220/256 255/256],'EdgeColor',[220/256 220/256 255/256]) % window for counting spikes
        hold on
        % plot each CCEP
        plot(t*1000, -1*squeeze(dataBase(nSubj).cc_epoch_sorted(selectChan,:,selectStim,:)),':','Color','k','LineWidth',0.5)
        % plot averaged CCEP
        plot(t*1000, -1*signal_avg_orig,'Color','k','LineWidth',2)

        title(sprintf('%s: Averaged CCEP: %s-%s, %s',...
            dataBase(nSubj).sub_label,dataBase(nSubj).cc_stimchans{selectStim,:},dataBase(nSubj).ch{selectChan}))
        ylim([-1500 3500])
        xlim([-100 100])
        xlabel('Time (ms)')

        % FIGURE ERSP
        nFig = nSubj+20;
        ERSPall = dataBase(nSubj).ERSP;
        [~, axes1] = plot_ERSP(ERSPall,selectStim,selectChan,nFig);
        colorbar(axes1)

        % FIGURE spike ratio
        figure(nSubj+30)
        minVal = -1; maxVal = 1;
        randNumbers = minVal+(maxVal-minVal)*rand(size(IEDmat));
        scatter(randNumbers(:),log(IEDmat(:)),20,[220/256 220/256 255/256],'filled')
        hold on
        scatter(randNumbers(selectStim,column),log(IEDmat(selectStim,column)),30,'k','filled')
        hold off
        xlim([-5 5])
        title(sprintf('Spike ratio of %s',dataBase(nSubj).sub_label))
        ylabel('Logarithmic spike ratio')
    end

end

% housekeeping
clear axes1 bothCCEP_ERSP column CCEPmat CCEPmatIED ERSPall ERSPmat ERSPmatIED
clear fig1 I IEDmat lower_err nSignal row selectChan selResp selectStim signal_avg signal_avg_orig
clear signal_avg_orig2 signal_nonart signal_orig signal_spikes signal_sterr
clear stimart_CCEPs stimart_spikes nSubj upper_err f maxVal minVal nFig randNumbers ymaxSpikes yminSpikes

%% save figures for manuscript

nSubj = 8;
print(figure(nSubj),'-vector','-depsc',fullfile(myDataPath.figures,['fig1_IED_' dataBase(nSubj).sub_label]))
print(figure(nSubj+10),'-vector','-depsc',fullfile(myDataPath.figures,['fig1_CCEP_' dataBase(nSubj).sub_label]))
print(figure(nSubj+20),'-vector','-depsc',fullfile(myDataPath.figures,['fig1_ERSP_' dataBase(nSubj).sub_label]))
print(figure(nSubj+30),'-vector','-depsc',fullfile(myDataPath.figures,['fig1_logSR_' dataBase(nSubj).sub_label]))

fprintf('Figures for subject %d are saved. \n',nSubj)

%housekeeping
clear nSubj

%% end of script