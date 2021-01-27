% made by Michelle vd Stoel
% 2018
% made BIDS compatible by Dorien van Blooijs
% september 2019

function detPar = getfeaturesTrain(subject)

% extract variables from subject
allERSP = subject.ERSP.allERSPboot;
stimpchan = subject.ERSP.cc_stimsets;
idx_ch_bad = subject.idx_ch_bad;
times = subject.ERSP.times;

% Divide hys in pre- and post stimulation
time(1,:) = times<0;
time(2,:) = times>0;

%%% variables needed to plot ERSP with bounding box
% freqs = subject.ERSP.freqs;
% pixelWidth = (times(end)-times(1))/(size(times,2)-1);
% pixelHeigth = (freqs(end)-freqs(1))/(size(freqs,2)-1);
% t(1) = find(times<0,1,'first');
% t(2) = find(times>0,1,'first');
% C = [0.9 0.9 0.9]; % light grey

% pre-allocation
Area = cell(2,1);
tStart = cell(2,1);
fStart = cell(2,1);
tWidth = cell(2,1);
fWidth = cell(2,1);

for stimpair = 1:size(allERSP,1) % for each stimulation pair
    for chan = 1:size(allERSP,2) % for each recording electrode
                
        % delineate the power suppression in bootstrapped ERSP
        ERSP = allERSP{stimpair,chan};
        ERSPfilt = imgaussfilt(ERSP,0.5,'Filterdomain','spatial');

        idx_ERSP = ERSPfilt<0;

        % %%% plot ERSP (optionally)
        % plot_ERSP(subject.ERSP,stimpair,chan)
        % hold on
                
        % pre-allocation
        hys_div = NaN(2,size(ERSPfilt,1),size(ERSPfilt,2)/2);
        for win = 1:2 % pre and post stimulation
            hys_div(win,:,:) = idx_ERSP(:,time(win,:));

            partfig = ERSPfilt(:,time(win,:));
            
            %%% Get area and other features
            stats = regionprops(logical(squeeze(hys_div(win,:,:))),'Area','BoundingBox','FilledImage','PixelIdxList');
                                  
            %%% Get area and duration in a double           
            if ismember(chan,stimpchan(stimpair,:)) % fill with NaNs when in stimulus pair
                Area{win}(stimpair,chan) = NaN;
                tStart{win}(stimpair,chan) = NaN;
                fStart{win}(stimpair,chan) = NaN;
                tWidth{win}(stimpair,chan) = NaN;
                fWidth{win}(stimpair,chan) = NaN;
                
            elseif idx_ch_bad(chan) == 1 % fill with NaNs when bad channel (noisy)
                Area{win}(stimpair,chan) = NaN;
                tStart{win}(stimpair,chan) = NaN;
                fStart{win}(stimpair,chan) = NaN;
                tWidth{win}(stimpair,chan) = NaN;
                fWidth{win}(stimpair,chan) = NaN;
                
            elseif isempty(stats)
                Area{win}(stimpair,chan) = 0;
                tStart{win}(stimpair,chan) = 0;
                fStart{win}(stimpair,chan) = 0;
                tWidth{win}(stimpair,chan) = 0;
                fWidth{win}(stimpair,chan) = 0;
                
            else
                
                % determine mean intensity of all pixels in FilledImage
                for i=1:size(stats,1)
                    stats(i).meanIntensity = mean(partfig(stats(i).PixelIdxList));
                end
                
                % find largest area with largest mean intensity
                [~,sort_idx] = sort([stats.Area].*abs([stats.meanIntensity]),'descend');
                idx = sort_idx(1);
                
                Area{win}(stimpair,chan) = stats(idx).Area;
                tStart{win}(stimpair,chan) = stats(idx).BoundingBox(1);
                fStart{win}(stimpair,chan) = stats(idx).BoundingBox(2);
                tWidth{win}(stimpair,chan) = stats(idx).BoundingBox(3);
                fWidth{win}(stimpair,chan) = stats(idx).BoundingBox(4);
            end
                        
%             %%% plot bounding box in ERSP plot (optionally)
%             if ~isnan(tStart{win}(stimpair,chan)) 
%                 %%% plot bounding box on largest power suppression
%                 if tStart{win}(stimpair,chan) > 0.5
%                     x1 = (times(floor(tStart{win}(stimpair,chan)+t(win)-1)) + times(ceil(tStart{win}(stimpair,chan)+t(win)-1)))/2;
%                 else
%                     x1 = times(1)-0.5*pixelWidth;
%                 end
%                 
%                 if fStart{win}(stimpair,chan) > 0.5
%                     y1 = (freqs(floor(fStart{win}(stimpair,chan))) + freqs(ceil(fStart{win}(stimpair,chan))))/2;
%                 else
%                     y1 = freqs(1)-0.5*pixelHeigth;
%                 end
%                 w = tWidth{win}(stimpair,chan)*pixelWidth;
%                 h = fWidth{win}(stimpair,chan)*pixelHeigth;
%                 
%                 rectangle('Position',[x1 y1 w h],'EdgeColor',C),
%                 text(x1+0.5*w,y1+0.5*h,num2str(Area{win}(stimpair,chan)))
%             end

        end
%         hold off
%         pause
    end
end

% concatenate into two vectors [channels*stimuli x 2] (pre- and post-stimulus)
detPar.Area_conc = [reshape(Area{1},numel(Area{1}),1),reshape(Area{2},numel(Area{2}),1)];
detPar.tStart_conc = [reshape(tStart{1},numel(tStart{1}),1),reshape(tStart{2},numel(tStart{2}),1)];
detPar.fStart_conc = [reshape(fStart{1},numel(fStart{1}),1),reshape(fStart{2},numel(fStart{2}),1)];
detPar.tWidth_conc = [reshape(tWidth{1},numel(tWidth{1}),1),reshape(tWidth{2},numel(tWidth{2}),1)];
detPar.fWidth_conc = [reshape(fWidth{1},numel(fWidth{1}),1),reshape(fWidth{2},numel(fWidth{2}),1)];

