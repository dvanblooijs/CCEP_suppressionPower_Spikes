% made by Michelle vd Stoel
% 2018
% made BIDS compatible by Dorien van Blooijs
% september 2019

function [Area_conc, tStart_conc, fStart_conc, tWidth_conc, fWidth_conc, ...
    Area, tStart,  fStart, tWidth, fWidth] = getfeaturesTrain(subject, ThL, ThU)

t = subject.ERSP.times;

% Divide hys in pre- and post stimulation
time(1,:) = t<0;
time(2,:) = t>0;

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
        
        % hysteresis3d assumes non-negative image. We want to detect "power suppression
        % (negative values). Therefor, we multiply the ERSP with -1
        % and remove all values <0.
        ERSP = allERSP{stimpair,chan};
        ERSP2 = -1* ERSP;
        ERSP2(ERSP2<0) = 0;
        
        % delineate blue blobs using hysteresis3d
        t1=ThL;                                          % Lower threshold
        t2=ThU;                                          % Upper threshold
        conn=8;                                             % Connectivity
        [~,hys]=hysteresis3d(ERSP2,t1,t2,conn);
        
        % pre-allocation
        hys_div = NaN(2,size(ERSP2,1),size(ERSP2,2)/2);
        for win = 1:2 % pre and post stimulation
            hys_div(win,:,:) = hys(:,time(win,:));
            
            %%Get area and other features
            stats = regionprops(logical(squeeze(hys_div(win,:,:))),'Area','BoundingBox');
                                  
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
                [val,idx] = max([stats.Area]);
                Area{win}(stimpair,chan) = val;
                tStart{win}(stimpair,chan) = stats(idx).BoundingBox(1);
                fStart{win}(stimpair,chan) = stats(idx).BoundingBox(2);
                tWidth{win}(stimpair,chan) = stats(idx).BoundingBox(3);
                fWidth{win}(stimpair,chan) = stats(idx).BoundingBox(4);
            end
        end
                
    end
end

Area_conc = [reshape(Area{1},numel(Area{1}),1),reshape(Area{2},numel(Area{2}),1)];
tStart_conc = [reshape(tStart{1},numel(tStart{1}),1),reshape(tStart{2},numel(tStart{2}),1)];
fStart_conc = [reshape(fStart{1},numel(fStart{1}),1),reshape(fStart{2},numel(fStart{2}),1)];
tWidth_conc = [reshape(tWidth{1},numel(tWidth{1}),1),reshape(tWidth{2},numel(tWidth{2}),1)];
fWidth_conc = [reshape(fWidth{1},numel(fWidth{1}),1),reshape(fWidth{2},numel(fWidth{2}),1)];

