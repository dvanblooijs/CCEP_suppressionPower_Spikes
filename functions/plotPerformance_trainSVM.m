% This function constructs heatmaps with parameter-ranges on x- and y-axis
% and values of performances filled in each pixel.
% NEEDED PARAMETERS INPUT:
% - train_threshold     : struct containing the performance for each combination of optimizing parameter ranges

function plotPerformance_trainSVM(train_threshold,trainPar)

modes = {'L','sens','spec','prec','F'};
cmap = colormap(parula(size(train_threshold,2)+1));
figure(),

for m = 1:size(modes,2)
    subplot(2,3,m)
    
    hold on
    for subj = 1:size(train_threshold,2)
        h(subj) = plot(trainPar.C,train_threshold(subj).(modes{m}),'Color',cmap(subj,:)); %#ok<AGROW>
        allval = mean(vertcat(train_threshold(:).(modes{m})));
        h(size(train_threshold,2)+1) = plot(trainPar.C, allval,'-.','Color','k','LineWidth',1); %#ok<AGROW>
        
        if strcmp(modes{m},'L')
            [maxval,I] = min(train_threshold(subj).(modes{m}));
            [maxallval,allI] = min(allval);
            
            ymin = floor(min(min([train_threshold(:).(modes{m})]))*10)/10;
            ymax = ceil(max(max([train_threshold(:).(modes{m})]))*10)/10;
        else
            [maxval,I] = max(train_threshold(subj).(modes{m}));
            [maxallval,allI] = max(allval);
            
            ymin = 0.5;
            ymax = 1;
        end
        
        plot(trainPar.C(I),maxval,'o','MarkerFaceColor',cmap(subj,:),'MarkerEdgeColor',cmap(subj,:),'MarkerSize',3)
        plot(trainPar.C(allI),maxallval,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
    end
    
    
    if strcmp(modes{m},'sens')
        fig_title = 'Sensitivity';
    elseif strcmp(modes{m},'spec')
        fig_title = 'Specificity';
    elseif strcmp(modes{m},'prec')
        fig_title= 'Precision';
    else
        fig_title = modes{m};
    end
    
    ax = gca;
    ax.XAxis.Scale = 'log';
    ax.XTick = 10.^(0:2);
    ax.XAxis.Label.String = 'C';
    ylim([ymin ymax])
    title(fig_title)
    lgd = legend(h,{train_threshold(:).subj, 'combined'});
    lgd.Position = [0.7 0.3 0.2 0.1];
end
end
