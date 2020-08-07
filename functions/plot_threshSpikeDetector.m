function plot_threshSpikeDetector(dataBase,train_threshold,train_thresh,train_t_dif)

abs_thresh = NaN(4,1); abs_thresh_chan = cell(4,1); thresh_SD = NaN(4,2); thresh_t = NaN(4,2);

for subj = 1:size(dataBase,2)
    [SD,t] = find(train_threshold(subj).F(:,:,1) == max(max(train_threshold(subj).F(:,:,1))),1,'first');
    
    thresh_SD(subj,1) = train_thresh(SD);
    thresh_t(subj,1) = train_t_dif(t);
    
    abs_thresh_chan{subj} = dataBase(subj).spikes.Pharmat_norm(:,3) + thresh_SD(subj,1) * dataBase(subj).spikes.Pharmat_norm(:,2);
    
    [SD,t] = find(train_threshold(subj).F(:,:,2) == max(max(train_threshold(subj).F(:,:,2))),1,'first');
    
    thresh_SD(subj,2) = train_thresh(SD);
    thresh_t(subj,2) = train_t_dif(t);
    
    abs_thresh(subj,1) = dataBase(subj).spikes.Pharmatall_norm(3) + thresh_SD(subj,2) * dataBase(subj).spikes.Pharmatall_norm(2);
    
end

figure('Position',[680 287 1050 680]),
subplot(1,3,1)
scatter(1:4,abs_thresh(:,1),'b','filled')
hold on
for subj=1:4
scatter(subj*ones(size(abs_thresh_chan{subj})),abs_thresh_chan{subj},'r')
end
hold off
title('Absolute threshold')
% ylim([0 round(max(abs_thresh(:))+0.1*max(abs_thresh(:)))])
legend('Per patient','Per channel','Location','Best')
ylabel('Threshold value')
xlim([0, size(dataBase,2)+1])
ax = gca;
ax.XTick = 0:size(dataBase,2)+1;
ax.XTickLabel = {' ', dataBase.sub_label,' '};
ax.XTickLabelRotation = -45;

subplot(1,3,2)
scatter(1:4,thresh_SD(:,2),'b','filled')
hold on
scatter(1:4,thresh_SD(:,1),'r')
hold off
title('Scale threshold')
ylim([min(train_thresh), max(train_thresh)])
legend('Per patient','Per channel','Location','Best')
ylabel('Threshold value')
xlim([0, size(dataBase,2)+1])
ax = gca;
ax.XTick = 0:size(dataBase,2)+1;
ax.XTickLabel = {' ', dataBase.sub_label,' '};
ax.XTickLabelRotation = -45;

subplot(1,3,3)
scatter(1:4,thresh_t(:,2),'b','filled')
hold on
scatter(1:4,thresh_t(:,1),'r')
hold off
title('Time dif threshold')
ylim([0 1])
legend('Per patient','Per channel','Location','Best')
ylabel('Threshold value')
xlim([0, size(dataBase,2)+1])
ax = gca;
ax.XTick = 0:size(dataBase,2)+1;
ax.XTickLabel = {' ', dataBase.sub_label,' '};
ax.XTickLabelRotation = -45;
end