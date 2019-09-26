function [D,A,I,C,F,Fstart,Fend,Tstart]=getfeatures(t,f, allERSP, ThL, ThU, stimpchan)

%F,Fstart,Fend

ii=1;

for n=1:size(allERSP,1)                            % for each stimulation pair
    for m=1:size(allERSP,2)                        % for each recording electrode
        ERSP2 = allERSP{n,m};
        t1=ThL;                                          % Lower threshold
        t2=ThU;                                          % Upper threshold, ook proberen voor 10-15
        conn=8;                                             % Connectivity
        [tri,hys]=hysteresis3d(ERSP2,t1,t2,conn);
        allhys{n,m}= hys;
        alltri{n,m}= tri;
        
        %%Get area and other features
        stats=regionprops(hys,ERSP2,'Area', 'Image','MeanIntensity');
        
        %%Get duration
        Ts=t(2)-t(1);
        Fs=1/Ts;
        
        for i=1:size(stats,1)                               % for each area per picture get duration/area
            img = stats(i).Image;
            for j=1:size(img,1)                             % for each row in the image
                index=find(img(j,:)==1);
                if isempty(index)
                    duration(j) = 0;
                    tstart(j) = 0;
                else
                    tstart(j)=index(1);
                    tend=index(end);
                    duration(j) = (tend-tstart(j))/Fs;
                end
                [stats(i).duration,K] = max(duration);          % get row with maximum duration
                stats(i).tstart = tstart(K)/Fs;
            end
            clear duration
            
            for jj=1:size(img,2)                            % for each row in the image
                index2=find(img(:,jj)==1);
                if isempty(index2)
                    freq_range(jj) = 0;
                else
                    fstart(jj)=f(index2(1));
                    fend(jj)=f(index2(end));
                    freq_range(jj) = fend(jj)-fstart(jj);   % sample frequency = 1
                end
                [stats(i).freq,Ix] = max(freq_range);        % get row with maximum duration
                stats(i).fend = fend(Ix);
                stats(i).fstart = fstart(Ix);
            end
            clear freq_range fend fstart
        end
        allstats{n,m} = stats;
        clear ERSP2 hys
        
        %%Get area and duration in a double
        l = stimpchan(n,:);
        
        if isempty(stats)
            D(ii)=0;
            A(ii)=0;
            I(ii)=0;
            C{n,m}=0;
            F(ii)=0;
            Fend(ii)=0;
            Fstart(ii)=0;
            
        elseif m==l(1) || m==l(2)                   % don't take recording which is in stimulus pair
            D(ii)=0;
            A(ii)=0;
            I(ii)=0;
            C{n,m}=0;
            F(ii)=0;
            Fend(ii)=0;
            Fstart(ii)=0;
        else
            [val,idx] = max([stats.Area]);
            D(ii)=stats(idx).duration;
            A(ii)=val;
            I(ii)=stats(idx).MeanIntensity;
            C{n,m}=A(ii);
            F(ii)=stats(idx).freq;
            Fend(ii)=stats(idx).fend;
            Fstart(ii)=stats(idx).fstart;
            Tstart(ii)=stats(idx).tstart;
        end
        
        ii=ii+1;
    end
end