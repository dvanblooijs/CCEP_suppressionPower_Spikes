% Figure 1: example of ten responses to spes stimulation, an CCEP, and an
% ERSP when there is a large change in spike ratio

%% load first ccepSp03_analysis_ERs_PS_spikes.m

%% plot a figure for each subject showing the spikes before and after stimulation
% in ten epochs, in whom power suppression and a CCEP was observed

close all

cfg = dataBase(1).ERSP.cfg;

t = -cfg.epoch_prestim + 1/2048 : 1/2048: cfg.epoch_length - cfg.epoch_prestim;
stimart_spikes = t>-0.1 & t < 0.1;
stimart_CCEPs = t>=0 & t < 0.01;

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

        selResp = bothCCEP_ERSP(I(1)); % select the one stim-resp combination with the most negative spike ratio (spike decrease)

        [r,c] = ind2sub(size(IEDmat),selResp);

        selChan = dataBase(subj).IEDs.IEDch(c);
        selStim = r;

        figure(subj),
        for n = 1:10

            signal_orig = -1*squeeze(dataBase(subj).cc_epoch_sorted(selChan,n,selStim,:))+n*1400;
            signal_spikes = signal_orig;
            signal_spikes(stimart_spikes) = NaN;
            signal_nonart = signal_orig;
            signal_nonart(stimart_CCEPs) = NaN;

            hold on
            scatter(t(dataBase(subj).IEDs.totEpochspikessamp{r,c,n}),...
                signal_spikes(dataBase(subj).IEDs.totEpochspikessamp{r,c,n}),'o',...
                'MarkerEdgeColor','k','MarkerFaceColor',[220/256 220/256 255/256])
            plot(t,signal_orig,'b','Color',[220/256 220/256 255/256]) % original signal
            plot(t,signal_nonart,'b','Color',[45/256 45/256 255/256]) % signal without artefact
            plot(t,signal_spikes,'b','Color',[155/256 155/256 255/256]) % signal for spikes


        end
        title(sprintf('%s: %s-%s, %s', ...
            dataBase(subj).sub_label,dataBase(subj).cc_stimchans{selStim,:},dataBase(subj).ch{selChan}))
        xlabel('Time (s)')
        ylim([0 16000])

        % plot CCEP
        signal_avg_orig = squeeze(dataBase(subj).cc_epoch_sorted_avg(selChan,selStim,:));
        signal_avg_orig2 = mean(squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)));
        signal_sterr = std(squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)))./sqrt(10);
        lower_err = -1*signal_avg_orig2 - signal_sterr;
        upper_err = -1*signal_avg_orig2 + signal_sterr;
        signal_avg = signal_avg_orig;
        signal_avg(stimart_CCEPs) = NaN;

        figure(subj+10)
        hold on
        fill([t t(end:-1:1)],[upper_err lower_err(end:-1:1)],[150/256 150/256 255/256],'EdgeColor',[1 1 1])
        %         plot(t,-1*squeeze(dataBase(subj).cc_epoch_sorted(selChan,:,selStim,:)),'Color',[150/256 150/256 255/256])
        plot(t, -1*signal_avg_orig,'Color',[220/256 220/256 255/256],'LineWidth',2)
        plot(t, -1*signal_avg,'Color',[45/256 45/256 255/256],'LineWidth',2)

        title(sprintf('%s: Averaged CCEP: %s-%s, %s',...
            dataBase(subj).sub_label,dataBase(subj).cc_stimchans{selStim,:},dataBase(subj).ch{selChan}))
        ylim([-1500 3500])
        xlim([-0.1 0.2])
        xlabel('Time (s)')

        n = subj+20;
        ERSPall = dataBase(subj).ERSP;
        [fig1, axes1] = plot_ERSP(ERSPall,selStim,selChan,n);
        colorbar(axes1)
    end

end

%%
subj = 8;
print(figure(subj),'-painters','-depsc',fullfile(myDataPath.Figures,['fig1_IED_' dataBase(subj).sub_label]))
print(figure(subj+10),'-painters','-depsc',fullfile(myDataPath.Figures,['fig1_CCEP_' dataBase(subj).sub_label]))
print(figure(subj+20),'-painters','-depsc',fullfile(myDataPath.Figures,['fig1_ERSP_' dataBase(subj).sub_label]))

