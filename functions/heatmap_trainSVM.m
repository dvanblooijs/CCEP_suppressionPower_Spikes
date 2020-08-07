% This function constructs heatmaps with parameter-ranges on x- and y-axis
% and values of performances filled in each pixel. 
% NEEDED PARAMETERS INPUT:
% - cfg.train           : subjects in train-dataset
% - dataBase.sub_label  : names of the subjects
% - train_threshold     : struct containing the performance for each combination of optimizing parameter ranges
% - mode                : should be L/F/sens/prec/spec --> what performance to visualize
% - par1                : what first optimized parameter to visualize
% - par2                : what second optimized parameter to visualize

function heatmap_trainSVM(cfg,dataBase,train_threshold,trainPar,mode,par1,par2)

if (strcmp(par1,'C') && strcmp(par2,'ThU')) || (strcmp(par2,'C') && strcmp(par1,'ThU'))
    val = 2;
    Xval = trainPar.C;
    Yval = trainPar.ThU;
    par1 = 'C'; par2 = 'ThU';
elseif (strcmp(par1,'ThL') && strcmp(par2,'ThU')) || strcmp(par2,'ThL') && strcmp(par1,'ThU')
    val = 3;
    Xval = trainPar.ThL;
    Yval = trainPar.ThU;
    par1 = 'ThL'; par2 = 'ThU';
elseif strcmp(par1,'ThL') && strcmp(par2,'C') || strcmp(par2,'ThL') && strcmp(par1,'C')
    val = 1;
    Xval = trainPar.C;
    Yval = trainPar.ThL;
    par1 = 'C'; par2 = 'ThL';
else 
    error('Parameters are not found')
end

if strcmp(mode,'L')
    minval = min(min(min([train_threshold(:).(mode)])));
    maxval = max(max(min([train_threshold(:).(mode)],[],val)));
    cmap = colormap(autumn);
    cmap = flipud(cmap);
else
    if val ~= 2
        minval = min(min(max([train_threshold(:).(mode)],[],val)));
        maxval = max(max(max([train_threshold(:).(mode)])));
        
    elseif val ==2
        allmax = NaN(size(cfg.train));
        for subj = cfg.train
            allmax(subj) = min(min(max([train_threshold(subj).(mode)],[],val)));
        end
        minval = min(allmax);
        maxval = max(max(max([train_threshold(:).(mode)])));
        
    end
    cmap = colormap(autumn);
end

figure,
for subj=1:size(train_threshold,2)
    if size(train_threshold,2)>1
        
        subplot(1,3,subj)
        fig_title = sprintf('%s, %s',dataBase(subj).sub_label,mode);
    elseif size(train_threshold,2) ==1
        fig_title = 'Performance for train dataset combined';
    end
    if strcmp(mode,'L')
        a = squeeze(min(train_threshold(subj).(mode),[],val)); % uses the minimal error for gridC
    else
        a = squeeze(max(train_threshold(subj).(mode),[],val)); % uses the minimal error for gridC
    end
    
    imagesc(a,[minval, maxval])
    colormap(cmap)
    ax = gca;
    ax.XTick = 1:length(Xval);
    ax.XTickLabel = Xval;
    ax.YTick = 1:length(Yval);
    ax.YTickLabel = Yval;
    title(fig_title)
    xlabel(par1)
    ylabel(par2)
colorbar
end
