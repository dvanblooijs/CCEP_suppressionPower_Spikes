function [fig1, axes1] = plot_ERSP(ERSPall,stimp,chan,n)

% ERSP = ERSPall.allERSP;
ERSPboot = ERSPall.allERSPboot;
times = ERSPall.times;
freqs = ERSPall.freqs;
ch = ERSPall.ch;
cc_stimchans = ERSPall.cc_stimchans;

%Plot image manually
fig1 = figure(n);
colormap(jet);
axes1 = axes('Parent',fig1);
hold(axes1,'on');

% subplot(1,2,1)
% a = ERSP{stimp,chan};
% imagesc(times,freqs,a,[-15,15])
% colormap jet
% title(sprintf('ERSP stimpair = %g, chan = %g without bootstrapping',stimp,chan))
% ax = gca;
% ax.YDir = 'normal';
% % ax.YTick = min(freqs):50:max(freqs);
% % ax.YTickLabels = max(freqs):-50:min(freqs);
% % ax.XTick = 1:20:size(times,2);
% % ax.XTickLabels = round(times(ax.XTick),1,'significant');
% xlabel('Time(ms)')
% ylabel('Frequency (Hz)')
% 

% subplot(1,2,2),
a = ERSPboot{stimp,chan};
imagesc([times(1) times(end)],[freqs(1) freqs(end)],a,[-15,15])
title(sprintf('ERSP stimpair = %s-%s, chan = %s with bootstrapping',cc_stimchans{stimp,1},cc_stimchans{stimp,2},ch{chan}))

axes1.YDir = 'normal';
pixelWidth = (times(end)-times(1))/(size(times,2)-1);
pixelHeight = (freqs(end)-freqs(1))/(size(freqs,2)-1);
axes1.XLim = [times(1)-0.5*pixelWidth times(end)+0.5*pixelWidth];
axes1.YLim = [freqs(1)-0.5*pixelHeight freqs(end)+0.5*pixelHeight];
% ax.YTick = min(freqs):50:max(freqs);
% ax.YTickLabels = max(freqs):-50:min(freqs);
% ax.XTick = 1:20:size(times,2);
% ax.XTickLabels = round(times(ax.XTick),1,'significant');
xlabel('Time(ms)')
ylabel('Frequency (Hz)')

end