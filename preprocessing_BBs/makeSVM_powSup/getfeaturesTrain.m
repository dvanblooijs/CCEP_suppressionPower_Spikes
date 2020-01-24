% made by Michelle vd Stoel
% 2018
% made BIDS compatible by Dorien van Blooijs
% september 2019

function [D,A,Dmat,Amat] = getfeaturesTrain(subject, ThL, ThU,trainPar)

t = subject.ERSP.times;
%%Divide hys in pre- and post stimulation
time(1,:) = t<0;
time(2,:) = t>0;

if strcmp(trainPar.boot,'yes')
    allERSP = subject.ERSP.allERSPboot;
elseif strcmp(trainPar.boot,'no')
        allERSP = subject.ERSP.allERSP;
else
    error('bootstrapping on/of is not defined')
end
stimpchan = subject.ERSP.cc_stimsets;

%         Ts=t(2)-t(1);
%         Fs=1/Ts;
Ts = median(unique(diff(t)))/1000; % /1000 to convert from ms to s
Fs = 1/Ts;
Dmat = cell(2,1);
Amat = cell(2,1);

for stimpair=1:size(allERSP,1)                            % for each stimulation pair
    for chan=1:size(allERSP,2)                        % for each recording electrode
        
        % hysteresis3d assumes non-negative image. We want to detect "blue
        % blobs" (negative values). Therefor, I multiply the ERSP with -1
        % and remove all values <0.
        ERSP = allERSP{stimpair,chan};
        ERSP2 = -1* ERSP;
        ERSP2(ERSP2<0) = 0;
        
        % delineate blue blobs using hysteresis3d
        t1=ThL;                                          % Lower threshold
        t2=ThU;                                          % Upper threshold
        conn=8;                                             % Connectivity
        [~,hys]=hysteresis3d(ERSP2,t1,t2,conn);
%         allhys{stimpair,chan}= hys;
%         alltri{stimpair,chan}= tri;
        
        hys_div = NaN(2,size(ERSP2,1),size(ERSP2,2)/2);
        for win=1:2 % pre and post stimulation
            hys_div(win,:,:) = hys(:,time(win,:));
            
            %%Get area and other features
            stats=regionprops(squeeze(hys_div(win,:,:)),'Area', 'ConvexHull', 'EquivDiameter', 'EulerNumber', 'FilledImage', 'FilledArea', 'Image');
            
            %%Get maximum duration pre and post-stimulation
            for i=1:size(stats,1)                               % for each area per picture get duration/area
                
                img = stats(i).Image;               
                stats(i).duration = max(sum(img,2))*Ts;          % get row with maximum duration
                
            end
            statsPrePost(win).stats = stats;
                      
            %%Get area and duration in a double
            l = stimpchan(stimpair,:);
            
            if chan==l(1) || chan==l(2)                   % don't take recording which is in stimulus pair
                Dmat{win}(stimpair,chan) = NaN;
                Amat{win}(stimpair,chan) = NaN;
                
            elseif isempty(stats)
                Dmat{win}(stimpair,chan) = 0;
                Amat{win}(stimpair,chan) = 0;
                
            else
                [val,idx] = max([stats.Area]);
                Dmat{win}(stimpair,chan) = stats(idx).duration;
                Amat{win}(stimpair,chan) = val;
            end
        end
        
        allstats{stimpair,chan} = statsPrePost;

%         fprintf('Stimpair %d of %d, channel %d of %d\n',stimpair, size(allERSP,1),chan,size(allERSP,2))
        
    end
end

D = [reshape(Dmat{1},numel(Dmat{1}),1),reshape(Dmat{2},numel(Dmat{2}),1)];
A = [reshape(Amat{1},numel(Amat{1}),1),reshape(Amat{2},numel(Amat{2}),1)];

