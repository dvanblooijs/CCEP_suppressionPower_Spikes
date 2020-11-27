% made by Michelle vd Stoel
% 2018
% made BIDS compatible by Dorien van Blooijs
% september 2019

function [Area_conc, tStart_conc, fStart_conc, tWidth_conc, fWidth_conc, ...
    Area, tStart,  fStart, tWidth, fWidth] = getfeaturesTrain(subject)

times = subject.ERSP.times;
freqs = subject.ERSP.freqs;
pixelWidth = (times(end)-times(1))/(size(times,2)-1);
pixelHeigth = (freqs(end)-freqs(1))/(size(freqs,2)-1);

% Divide hys in pre- and post stimulation
time(1,:) = times<0;
time(2,:) = times>0;
t(1) = find(times<0,1,'first');
t(2) = find(times>0,1,'first');

allERSP = subject.ERSP.allERSPboot;
stimpchan = subject.ERSP.cc_stimsets;

% pre-allocation
Area = cell(2,1);
tStart = cell(2,1);
fStart = cell(2,1);
tWidth = cell(2,1);
fWidth = cell(2,1);

for stimpair=1:size(allERSP,1)                            % for each stimulation pair
    for chan=1:size(allERSP,2)                        % for each recording electrode
                
        % delineate the power suppression in bootstrapped ERSP
        ERSP = allERSP{stimpair,chan};
        idx_ERSP = ERSP<0;
        C = [0.9 0.9 0.9]; % light grey

%         plot_ERSP(subject.ERSP,stimpair,chan)
%         hold on
                
        % pre-allocation
        hys_div = NaN(2,size(ERSP,1),size(ERSP,2)/2);
        for win = 1:2 % pre and post stimulation
            hys_div(win,:,:) = idx_ERSP(:,time(win,:));
            
            %%Get area and other features
            stats = regionprops(logical(squeeze(hys_div(win,:,:))),'Area','BoundingBox','FilledImage');
                                  
            %%Get area and duration in a double           
            if ismember(chan,stimpchan(stimpair,:))           % don't take recording which is in stimulus pair
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
                [~,sort_idx] = sort([stats.Area],'descend');
                
                x1_lo = stats(sort_idx(1)).BoundingBox(1);
                x1_hi = stats(sort_idx(1)).BoundingBox(1) + stats(sort_idx(1)).BoundingBox(3);
                x2_lo = stats(sort_idx(2)).BoundingBox(1);
                x2_hi = stats(sort_idx(2)).BoundingBox(1) + stats(sort_idx(2)).BoundingBox(3);
                y1_lo = stats(sort_idx(1)).BoundingBox(2);
                y1_hi = stats(sort_idx(1)).BoundingBox(2) + stats(sort_idx(1)).BoundingBox(4);
                y2_lo = stats(sort_idx(2)).BoundingBox(2);
                y2_hi = stats(sort_idx(2)).BoundingBox(2) + stats(sort_idx(2)).BoundingBox(4);
               
                
                % if the largest and second largest area are in the same
                % time window, it might be that these are split due to the
                % continuous 50Hz noise
                if x2_lo < x1_hi && x2_hi > x1_lo && ...% the second largest area must start before the largest area end, and must stop after the largest area start
                        (max([y1_lo y2_lo]) - min([y1_hi y2_hi])) <50 % the difference between the higher and lower freq must be <50Hz
                
                    for i=1:2
                        y_lo = stats(sort_idx(i)).BoundingBox(2);
                        y_hi = stats(sort_idx(i)).BoundingBox(2) + stats(sort_idx(i)).BoundingBox(4);
                        x_lo = stats(sort_idx(i)).BoundingBox(1);
                        x_hi = stats(sort_idx(i)).BoundingBox(1) + stats(sort_idx(i)).BoundingBox(3);
                        
                        image = double(stats(sort_idx(i)).FilledImage);
                        image(image==1) = 2;
                        
                        hys_div(win,ceil(y_lo):floor(y_hi),ceil(x_lo):floor(x_hi)) = image;
                    end
                    
                    stats = regionprops(squeeze(hys_div(win,:,:)),'Area','BoundingBox','FilledImage');
                    
                    idx = 2;
                else
                    
                    idx = sort_idx(1);
                end
                    
                Area{win}(stimpair,chan) = stats(idx).Area;
                tStart{win}(stimpair,chan) = stats(idx).BoundingBox(1);
                fStart{win}(stimpair,chan) = stats(idx).BoundingBox(2);
                tWidth{win}(stimpair,chan) = stats(idx).BoundingBox(3);
                fWidth{win}(stimpair,chan) = stats(idx).BoundingBox(4);
            end
                        
%             if ~isnan(tStart{win}(stimpair,chan)) 
%                 % plot bounding box on largest power suppression
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
%             end
        end
%         hold off
%         pause
    end
end

Area_conc = [reshape(Area{1},numel(Area{1}),1),reshape(Area{2},numel(Area{2}),1)];
tStart_conc = [reshape(tStart{1},numel(tStart{1}),1),reshape(tStart{2},numel(tStart{2}),1)];
fStart_conc = [reshape(fStart{1},numel(fStart{1}),1),reshape(fStart{2},numel(fStart{2}),1)];
tWidth_conc = [reshape(tWidth{1},numel(tWidth{1}),1),reshape(tWidth{2},numel(tWidth{2}),1)];
fWidth_conc = [reshape(fWidth{1},numel(fWidth{1}),1),reshape(fWidth{2},numel(fWidth{2}),1)];

