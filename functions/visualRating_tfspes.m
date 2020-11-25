function dataBase = visualRating_tfspes(dataBase,subj)

times = dataBase(subj).ERSP.times;
freqs = dataBase(subj).ERSP.freqs;

t(1) = find(times<0,1,'first');
t(2) = find(times>0,1,'first');

if isfield(dataBase(subj).ERSP,'detected')
    detected = dataBase(subj).ERSP.detected;
else
    detected = ones(size(dataBase(subj).ERSP.allERSPboot));
end
    
vis_BB = zeros(size(dataBase(subj).ERSP.allERSPboot));
n=1;

for stimp = 1:size(dataBase(subj).ERSP.allERSPboot,1)
    for chan = 1:size(dataBase(subj).ERSP.allERSPboot,2)
        
        if ~ismember(chan,dataBase(subj).ERSP.cc_stimsets(stimp,:)) % && detected(stimp,chan) == 1
            plot_ERSP(dataBase(subj).ERSP,stimp,chan)
            hold on
            
            if detected(stimp,chan) == 1
               C = [0.7 0 0]; % darkred
            else
                C = [0.9 0.9 0.9]; % light grey
            end
            
            for win=1:2
                if dataBase(subj).ERSP.Area{win}(stimp,chan) >0
                    pixelWidth = (times(end)-times(1))/(size(times,2)-1);
                    pixelHeigth = (freqs(end)-freqs(1))/(size(freqs,2)-1);
                    
                    if dataBase(subj).ERSP.tStart{win}(stimp,chan) > 0.5
                        x1 = (times(floor(dataBase(subj).ERSP.tStart{win}(stimp,chan)+t(win)-1)) + times(ceil(dataBase(subj).ERSP.tStart{win}(stimp,chan)+t(win)-1)))/2;
                    else
                        x1 = times(1)-0.5*pixelWidth;
                    end
                    
                    if dataBase(subj).ERSP.fStart{win}(stimp,chan) > 0.5
                        y1 = (freqs(floor(dataBase(subj).ERSP.fStart{win}(stimp,chan))) + freqs(ceil(dataBase(subj).ERSP.fStart{win}(stimp,chan))))/2;
                    else
                        y1=freqs(1)-0.5*pixelHeigth;    
                    end
                    w = dataBase(subj).ERSP.tWidth{win}(stimp,chan)*pixelWidth;
                    h = dataBase(subj).ERSP.fWidth{win}(stimp,chan)*pixelHeigth;
                    
                    rectangle('Position',[x1 y1 w h],'EdgeColor',C),
                end
            end
            hold off,

            perc = n / size(dataBase(subj).ERSP.allERSPboot(:),1) *100;
            
            s = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this an ERSP? [y/n]: ',perc,dataBase(subj).ERSP.cc_stimchans{stimp,:},dataBase(subj).ERSP.ch{chan}),'s');
            
            if strcmp(s,'y')
                vis_BB(stimp,chan) = 1;
            end
        end
        n=n+1;
        
    end
end

dataBase(subj).ERSP.checked = vis_BB;

disp('Visual rating in completed')


% dataBase.ERSP.checked = NaN(size(chans,2),size(stimps,2));
% for stimp = stimps
%         for chan = chans
%             figure(1),
%             a=dataBase(1).ERSP.allERSPboot{stimp,chan};
%             a=flipud(a);
%             imagesc(dataBase(1).ERSP.times,dataBase(1).ERSP.freqs,a,[-15,15])
%             title(sprintf('Stimulated: %s-%s, response in: %s',dataBase(1).cc_stimchans{stimp,1},dataBase(1).cc_stimchans{stimp,2},dataBase(1).ch{chan}))
%             xlabel('Time(ms)')
%             ax = gca;
%             ax.YTick = min(dataBase(1).ERSP.freqs):50:max(dataBase(1).ERSP.freqs);
%             ax.YTickLabels = max(dataBase(1).ERSP.freqs):-50:min(dataBase(1).ERSP.freqs);
%             ylabel('Frequency (Hz)')
%             colormap jet
%             x = input('BB? [y/n] ','s');
%             if strcmp(x,'y')
%                 dataBase.ERSP.checked(chan,stimp) = 1 ;
%             else
%                 dataBase.ERSP.checked(chan,stimp) = 0 ;
%             end
%         end
% end
% 

end
