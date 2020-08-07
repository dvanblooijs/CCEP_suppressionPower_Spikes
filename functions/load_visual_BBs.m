function dataBase = load_visual_BBs(dataBase,cfg)

D = dir(cfg.dir_visrate);
filenames = {D.name};

for subj = 1 : size(dataBase,2)
        
    % pre-allocation
    BS_visscores = NaN([2,size(dataBase(subj).ERSP.cc_stimchans,1),size(dataBase(subj).ERSP.ch,1)]);    
    stimorder_visscores = cell([2,size(dataBase(subj).ERSP.cc_stimchans,1),1]);  
    chan_visscores = cell([2,size(dataBase(subj).ERSP.ch,1),1]);    
    
    files = find(contains(filenames,dataBase(subj).sub_label));
    if ~isempty(files)
        for i=1:size(files,2)
            % load visual scoring from excel-file
            [~,txt,raw] = xlsread(fullfile(cfg.dir_visrate,filenames{files(i)}));
            % remove first row with stimulation pairs (see txt)
            raw(1,:) = [];
            % remove first column with numbers of electrodes
            raw(:,1) = [];
            % remove NaN-values
            raw(all(cellfun(@(x) any(isnan(x)),raw),2),:) = [];
            raw = celltomat(raw);
            
            BS_visscores(i,:,:) = raw;
            stimorder_visscores(i,:,:) = txt(2:end,1);
            chan_visscores(i,:,:) = txt(1,2:end)';
        end
        
        % extract stimulus pairs from visual ERSP scorings
        if isequal(stimorder_visscores(1,:), stimorder_visscores(2,:)) 
            
            splitcell = cellfun( @(x) strsplit(x,'-'),stimorder_visscores(1,:),'UniformOutput',false);
            
            stimp_visscores = cell(size(splitcell,2),size(splitcell{1},2));
            stimpnum_visscores = NaN(size(splitcell,2),size(splitcell{1},2));
            for n=1:size(splitcell,2)
                for m = 1:size(splitcell{n},2)
                    stimp_visscores{n,m} = splitcell{n}{m};
                    if ~isempty(find(strcmpi(splitcell{n}{m},dataBase(subj).ERSP.ch), 1))
                        stimpnum_visscores(n,m) = find(strcmpi(splitcell{n}{m},dataBase(subj).ERSP.ch));
                    else
                        fprintf('ERROR: %s: %s is not found!\n',dataBase(subj).sub_label,splitcell{n}{m})
                    end
                end
            end
            
            dataBase(subj).vis_scores.stimorder = stimorder_visscores(1,:)';
            dataBase(subj).vis_scores.stimpnum = stimpnum_visscores;
            dataBase(subj).vis_scores.stimp = stimp_visscores;
        else
            disp('ERROR: order of stimulus pairs in visual scores differs')
        end
        
        % extract channels from visual ERSP scoring
        if isequal(chan_visscores(1,:),chan_visscores(2,:))
                        
            channum_visscores = NaN(size(chan_visscores,2),1);
            for n=1:size(chan_visscores,2)
                    if ~isempty(find(strcmpi(chan_visscores{1,n},dataBase(subj).ERSP.ch), 1))
                        channum_visscores(n,1) = find(strcmpi(chan_visscores{1,n},dataBase(subj).ERSP.ch));
                    else
                        fprintf('ERROR: %s: %s is not found!\n',dataBase(subj).sub_label,chan_visscores{1,n})
                    end
            end
            
            dataBase(subj).vis_scores.chan = chan_visscores(1,:)';
            dataBase(subj).vis_scores.channum = channum_visscores;

        else
            disp('ERROR: order of channels in visual scores differs')
        end        
        
        dataBase(subj).vis_scores.BS = BS_visscores;
    else
        disp('No visual scores for this patient')
    end
end