function plot_DataRerefStimfree(dataBase,subj,IED)

fs = dataBase(subj).ccep_header.Fs;

figure(subj),
subplot(2,1,1)
plot(1/fs:1/fs:size(dataBase(subj).data_avg,2)/fs,dataBase(subj).data_avg+500)
hold on
plot(1/fs:1/fs:size(dataBase(subj).data_avg,2)/fs,dataBase(subj).data(IED,:),'LineWidth',2,'Color',[0.8 0.8 0.8])
plot(1/fs:1/fs:size(dataBase(subj).data_avg,2)/fs,dataBase(subj).data_rerefnoStimArt(IED,:),'Color','k')
hold off
xlabel('Time (s)')
ylabel('uV')
legend('AVG sig','orig sig','reref w/o stimart sig')

subplot(2,1,2),
plot(1/fs:1/fs:size(dataBase(subj).data_avg,2)/fs,dataBase(subj).data_avg+500)
hold on
plot(1/fs:1/fs:size(dataBase(subj).data_avg,2)/fs,dataBase(subj).data(IED,:),'LineWidth',2,'Color',[0.8 0.8 0.8])
plot(1/fs:1/fs:size(dataBase(subj).data_avg,2)/fs,dataBase(subj).data_rerefnoStimArt(IED,:),'Color','k')
hold off
xlim([190 220])
xlabel('Time (s)')
ylabel('uV')
legend('AVG sig','orig sig','reref w/o stimart sig')
title(sprintf('Signal of channel %s of %s',dataBase(subj).ch{IED},dataBase(subj).sub_label))
end