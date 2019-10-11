% made by Michelle vd Stoel
% 2018
% made BIDS compatible by Dorien van Blooijs
% september 2019

function [S,D,A]=getfeaturesTrain(subject, ThL, ThU)

t = subject.ERSP.times;
%%Divide hys in pre- and post stimulation
time(1,:) = t<0;
time(2,:) = t>0;

allERSP = subject.ERSP.allERSPboot;
BS_score = subject.BS_visscores;
stimpchan = subject.ERSP.cc_stimsets;

%         Ts=t(2)-t(1);
%         Fs=1/Ts;
Ts = median(unique(diff(t)))/1000; % /1000 to convert from ms to s
Fs = 1/Ts;

ii=1;

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
        [tri,hys]=hysteresis3d(ERSP2,t1,t2,conn);
        allhys{stimpair,chan}= hys;
        alltri{stimpair,chan}= tri;
        
        
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
            
            %Get stimulus pair in visual scoring
            if isequal(subject.ERSP.ch,subject.chan_visscores)
                if isequal(subject.ERSP.cc_stimsets,subject.stimpnum_visscores)
                    stimp = stimpair;
                else
                    stimp = find(subject.ERSP.cc_stimsets(stimpair,1)==subject.stimpnum_visscores(:,1) & subject.ERSP.cc_stimsets(stimpair,2)==subject.stimpnum_visscores(:,2));
                end
            else
                stimp =find(strcmpi([subject.ERSP.cc_stimchans{stimpair,1},'-',subject.ERSP.cc_stimchans{stimpair,2}],subject.stimorder_visscores));
            end
            
            %Get indexes of similar scoring           
            if BS_score(1,stimp,chan)==BS_score(2,stimp,chan)
                S(ii,win) = BS_score(1,stimp,chan);
            else
                S(ii,win) = 0;
            end
            
            %%Get area and duration in a double
            l = stimpchan(stimpair,:);
            
            if chan==l(1) || chan==l(2)                   % don't take recording which is in stimulus pair
                D(ii,win)=NaN;
                A(ii,win)=NaN;
                
            elseif isempty(stats)
                D(ii,win)=0;
                A(ii,win)=0;
                
            else
                [val,idx] = max([stats.Area]);
                D(ii,win)=stats(idx).duration;
                A(ii,win)=val;
            end
        end
        ii=ii+1;
        
        allstats{stimpair,chan} = statsPrePost;

        fprintf('Stimpair %d of %d, channel %d of %d\n',stimpair, size(allERSP,1),chan,size(allERSP,2))
        
    end
end