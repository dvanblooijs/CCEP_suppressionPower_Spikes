% made by Michelle vd Stoel
% 2018
% made BIDS compatible by Dorien van Blooijs
% september 2019

function [S,D,A]=getfeaturesTrain(subject, ThL, ThU)

t = subject.ERSP.times;
allERSP = subject.ERSP.allERSPboot;
BS_score = subject.BS_visscores;
stimpchan = subject.ERSP.cc_stimsets;

%         Ts=t(2)-t(1);
%         Fs=1/Ts;
Ts = median(unique(diff(t)))/1000; % /1000 to convert from ms to s
Fs = 1/Ts;

ii=1;

for n=1:size(allERSP,1)                            % for each stimulation pair
    for m=1:size(allERSP,2)                        % for each recording electrode
        
        % hysteresis3d assumes non-negative image. We want to detect "blue
        % blobs" (negative values). Therefor, I multiply the ERSP with -1
        % and remove all values <0. 
        ERSP = allERSP{n,m};
        ERSP2 = -1* ERSP;
        ERSP2(ERSP2<0) = 0;
        
        % delineate blue blobs using hysteresis3d
        t1=ThL;                                          % Lower threshold
        t2=ThU;                                          % Upper threshold
        conn=8;                                             % Connectivity
        [tri,hys]=hysteresis3d(ERSP2,t1,t2,conn);
        allhys{n,m}= hys;
        alltri{n,m}= tri;
                
        %%Get area and other features
        stats=regionprops(hys,'Area', 'ConvexHull', 'EquivDiameter', 'EulerNumber', 'FilledImage', 'FilledArea', 'Image');
        
        %%Get duration
        for i=1:size(stats,1)                               % for each area per picture get duration/area
            img = stats(i).Image;
            for j=1:size(img,1)                             % for each row in the image
                index=find(img(j,:)==1);
                if isempty(index)
                    duration(j) = 0;
                else
                    tstart=index(1);
                    tend=index(end);
                    duration(j) = (tend-tstart)/Fs;
                end
                stats(i).duration = max(duration);          % get row with maximum duration
            end
            clear duration
        end
        allstats{n,m} = stats;
        clear ERSP2 hys
        
        %Get indexes of similar scoring
        if BS_score(1,n,m)==BS_score(2,n,m)
            S(ii) = BS_score(1,n,m);
            
        else
            S(ii) = 0;
        end
        
        %%Get area and duration in a double
        l = stimpchan(n,:);
        
        if isempty(stats)
            D(ii)=0;
            A(ii)=0;
            
        elseif m==l(1) || m==l(2)                   % don't take recording which is in stimulus pair
            D(ii)=NaN;
            A(ii)=NaN;
            
        else
            [val,idx] = max([stats.Area]);
            D(ii)=stats(idx).duration;
            A(ii)=val;
        end
        
        ii=ii+1;
    end
end