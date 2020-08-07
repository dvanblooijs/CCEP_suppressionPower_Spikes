
function plot_HeatmapThreshSpikeDetector(dataBase,train_thresholds_all,train_t_dif,train_thresh)

for subj = 1:size(dataBase,2)
    mode = {'sens','prec','F'};
    
    h=figure(subj);
    set(h,'Position',[258 301 1472 666],'Name',dataBase(subj).sub_label)
    sgtitle(dataBase(subj).sub_label)
    for i=1:size(mode,2)
        
        subplot(2,3,i)
        imagesc(squeeze(train_thresholds_all.(mode{i})(subj,:,:,1)))
        ylabel('SD thresh')
        xlabel('t dif')
        ax = gca;
        ax.XTickLabel = train_t_dif(ax.XTick);
        ax.YTickLabel = train_thresh(ax.YTick);
        title(sprintf('%s per channel',mode{i}))
        colorbar
        
        subplot(2,3,i+3)
        imagesc(squeeze(train_thresholds_all.(mode{i})(subj,:,:,2)))
        ylabel('SD thresh')
        xlabel('t dif')
        ax = gca;
        ax.XTickLabel = train_t_dif(ax.XTick);
        ax.YTickLabel = train_thresh(ax.YTick);
        title(sprintf('%s per patient',mode{i}))
        colorbar
    end
end

%%

h=figure(8);
set(h,'Position',[161 462 1568 483])

for i=1:size(mode,2)
    
    for chanpat = 1:2
        subplot(2,3,(chanpat-1)*3+i)
        imagesc(train_thresholds_all.([mode{i},'med'])(:,:,chanpat))
        ylabel('SD threshold')
        xlabel('t dif')
        ax = gca;
        ax.XTickLabel = train_t_dif(ax.XTick);
        ax.YTickLabel = train_thresh(ax.YTick);
        if chanpat==1
            title(sprintf('%s: Median thresholds combined, per channel',mode{i}))
        elseif chanpat==2
            title(sprintf('%s: Median thresholds combined, per patient',mode{i}))
        end
        colorbar
        
    end
    
end


