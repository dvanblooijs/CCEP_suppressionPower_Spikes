function dataBase = visualRating_tfspes(dataBase,subj,lookPic, myDataPath,endstimp)
%%
if ~isempty(myDataPath)
    filefolder = fullfile(myDataPath.ERSPoutput,dataBase(subj).sub_label,dataBase(subj).ses_label,dataBase(subj).run_label);
    if ~exist(filefolder,'dir')
        mkdir(filefolder)
    end
    save_check = 1;
    
else
    save_check = 0;
    
end

times = dataBase(subj).ERSP.times;
freqs = dataBase(subj).ERSP.freqs;

t(1) = find(times<0,1,'first');
t(2) = find(times>0,1,'first');

ERSPdet = dataBase(subj).ERSPdet;

% ERSPdet.cc_stimchans = dataBase(subj).ERSP.cc_stimchans;
% ERSPdet.cc_stimsets = dataBase(subj).ERSP.cc_stimsets;
% ERSPdet.ch = dataBase(subj).ERSP.ch;


% if no ERSPs are checked
if sum(strcmp(fieldnames(dataBase(subj)),'ERSPdet'))==0
    ERSPdet.checked = zeros(size(dataBase(subj).ERSP.allERSPboot));
elseif sum(strcmp(fieldnames(dataBase(subj).ERSPdet),'checked')) == 0
    ERSPdet.checked = zeros(size(dataBase(subj).ERSP.allERSPboot));
end

n=numel(ERSPdet.detected(1:endstimp,:))+1;

for stimp = endstimp+1:size(dataBase(subj).ERSP.allERSPboot,1)
    for chan = 1:size(dataBase(subj).ERSP.allERSPboot,2)
        
        if ~ismember(chan,dataBase(subj).ERSP.cc_stimsets(stimp,:)) && lookPic(stimp,chan)==1 % && detected(stimp,chan) == 1
            [fig1,  axes1] = plot_ERSP(dataBase(subj).ERSP,stimp,chan);
            figchild = 0;
            alreadyClicked = 0;
            hold on
            
            C = [0.7 0.7 0.7]; % light grey
            
            if isfield(dataBase(subj).ERSPdet, 'Area')
                
                for win=1:2
                    if dataBase(subj).ERSPdet.Area{win}(stimp,chan) >0
                        pixelWidth = (times(end)-times(1))/(size(times,2)-1);
                        pixelHeigth = (freqs(end)-freqs(1))/(size(freqs,2)-1);
                        
                        % plot bounding box on largest power suppression
                        if dataBase(subj).ERSPdet.tStart{win}(stimp,chan) > 0.5
                            x1 = (times(floor(dataBase(subj).ERSPdet.tStart{win}(stimp,chan)+t(win)-1)) + times(ceil(dataBase(subj).ERSPdet.tStart{win}(stimp,chan)+t(win)-1)))/2;
                        else
                            x1 = times(t(win))-0.5*pixelWidth;
                        end
                        
                        if dataBase(subj).ERSPdet.fStart{win}(stimp,chan) > 0.5
                            y1 = (freqs(floor(dataBase(subj).ERSPdet.fStart{win}(stimp,chan))) + freqs(ceil(dataBase(subj).ERSPdet.fStart{win}(stimp,chan))))/2;
                        else
                            y1 = freqs(1)-0.5*pixelHeigth;
                        end
                        wdth = dataBase(subj).ERSPdet.tWidth{win}(stimp,chan)*pixelWidth;
                        h = dataBase(subj).ERSPdet.fWidth{win}(stimp,chan)*pixelHeigth;
                        
                        rectangle('Parent',axes1,'Position',[x1 y1 wdth h],'EdgeColor',C);
                        figchild = figchild +1;
                    end
                end
            end
            
            perc = n / size(dataBase(subj).ERSP.allERSPboot(:),1) *100;
            
            s = input(sprintf('%2.1f %% %s --- stimpair = %s-%s chan = %s --- Do you observe power suppression (y/n) and correct detection (y/yn)? [y/yn/n]: ',...
                perc,dataBase(subj).sub_label, dataBase(subj).ERSP.cc_stimchans{stimp,:},dataBase(subj).ERSP.ch{chan}),'s');
            
            if strcmp(s,'yn') || strcmp(s,'cyn')
                tStart = dataBase(subj).ERSPdet.statsPost{stimp,chan}(:,1);
                fStart = dataBase(subj).ERSPdet.statsPost{stimp,chan}(:,2);
                tWidth = dataBase(subj).ERSPdet.statsPost{stimp,chan}(:,3);
                fWidth = dataBase(subj).ERSPdet.statsPost{stimp,chan}(:,4);
                
                % plot bounding box on all power suppression
                x1 = zeros(size(tStart)); y1 = zeros(size(tStart));
                wdth = zeros(size(tStart)); h = zeros(size(tStart));
                for i = 1:size(tStart,1)
                    if tStart(i) > 0.5
                        x1(i,1) = (times(floor(tStart(i)+t(2)-1)) + times(ceil(tStart(i)+t(2)-1)))/2;
                    else
                        x1(i,1) = times(1)-0.5*pixelWidth;
                    end
                    
                    if fStart(i) > 0.5
                        y1(i,1) = (freqs(floor(fStart(i))) + freqs(ceil(fStart(i))))/2;
                    else
                        y1(i,1) = freqs(1)-0.5*pixelHeigth;
                    end
                    wdth(i,1) = tWidth(i)*pixelWidth;
                    h(i,1) = fWidth(i)*pixelHeigth;
                    
                    rectangle('Parent',axes1,'Position',[x1(i) y1(i) wdth(i) h(i)],'EdgeColor',C);
                end
                
                currkey = 0;
                while ~strcmp(currkey, 'c')
                    
                    w = waitforbuttonpress;
                    
                    if w == 0 % 0 for mouse click, 1 for button press
                        
                        % if one area is alreay selected, set this to light
                        % grey again, when another area is selected
                        if alreadyClicked == 1
                            axes1.Children(end-(figchild+idx)).EdgeColor = [0.7 0.7 0.7];
                        end
                        
                        cp = get(gca,'CurrentPoint');
                        
                        idx = find(x1 < cp(1,1) & x1+wdth > cp(1,1) &...
                            y1 < cp(1,2) & y1+h > cp(1,2));
                        
                        if all(idx == 0)
                            fprintf('No region is selected, try again! \n')
                            alreadyClicked = 0;
                        else
                            axes1.Children(end-(figchild+idx)).EdgeColor = [0.3 0.3 0.3];
                            alreadyClicked = 1;
                            
                            fprintf('You selected: tStart = %3.1fs, fStart = %3.1fHz, tWidth = %3.1fs, fWidth = %3.1fHz \n',...
                                x1(idx), y1(idx), wdth(idx), h(idx))
                            
                            fprintf('Press "c" when correct \n')
                        end
                    elseif w == 1
                        currkey = get(gcf,'CurrentCharacter');
                    end
                end
                
                ERSPdet.checked(stimp,chan) = 1;
                
                fprintf('Selection was: tStart = %3.1fsamples, fStart = %3.1fsamples, tWidth = %3.1fsamples, fWidth = %3.1fsamples \n',...
                    dataBase(subj).ERSPdet.tStart{2}(stimp,chan), ...
                    dataBase(subj).ERSPdet.fStart{2}(stimp,chan), ...
                    dataBase(subj).ERSPdet.tWidth{2}(stimp,chan),...
                    dataBase(subj).ERSPdet.fWidth{2}(stimp,chan))
                
                % update the features of this blue blobs in case the wrong
                % selection was made
                dataBase(subj).ERSPdet.tStart{2}(stimp,chan) = tStart(idx);
                dataBase(subj).ERSPdet.fStart{2}(stimp,chan) = fStart(idx);
                dataBase(subj).ERSPdet.tWidth{2}(stimp,chan) = tWidth(idx);
                dataBase(subj).ERSPdet.fWidth{2}(stimp,chan) = fWidth(idx);
                
                fprintf('Selection is now saved as: tStart = %3.1fsamples, fStart = %3.1fsamples, tWidth = %3.1fsamples, fWidth = %3.1fsamples \n \n',...
                    dataBase(subj).ERSPdet.tStart{2}(stimp,chan), ...
                    dataBase(subj).ERSPdet.fStart{2}(stimp,chan), ...
                    dataBase(subj).ERSPdet.tWidth{2}(stimp,chan),...
                    dataBase(subj).ERSPdet.fWidth{2}(stimp,chan))
                
            elseif strcmp(s,'y') || strcmp(s,'cy')
                ERSPdet.checked(stimp,chan) = 1;
                
            else
                ERSPdet.checked(stimp,chan) = 0;
                
            end
            close(fig1)
        end
        n=n+1;
        
    end
    
    % save also till which stimpair visual ERSPS were checked.
    ERSPdet.checkUntilStimp = stimp;
    
    if save_check == 1
        filename = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPdet.mat'];
        
        % save file during scoring in case of error
        save(fullfile(filefolder,filename),'-struct','ERSPdet');
    end
end

%%
if save_check == 1
    filename = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).task_label,'_',dataBase(subj).run_label,'_ERSPdet.mat'];

    % save file after scoring for completeness
    save(fullfile(filefolder,filename),'-struct','ERSPdet');
end

dataBase(subj).ERSPdet.checked = ERSPdet.checked;
dataBase(subj).ERSPdet.checkUntilStimp = ERSPdet.checkUntilStimp;

disp('Visual rating in completed')

end
