%% Peak detection according to Gaspard2013
% author: Dorien van Blooijs
% dec2018: construct script
% july 2019: adapted to BIDS readible

function [M,norm_M, Pharmat_norm, Pharmatall_norm] = findMahalanobisDist(data_stimfree,fs)

% step 1a: filterstop data between 45-55Hz
[b,a] = butter(4,[45/(fs/2) 55/(fs/2)],'stop'); %8th order butterworth filter between 45-55Hz

filt_data_transftemp = filtfilt(b,a,data_stimfree');
filt_datatemp = filt_data_transftemp';
disp('...The data is band-stop filtered between 45-55Hz...')

% step 1b: filter data between 10-70Hz
[b,a] = butter(4,[10/(fs/2) 70/(fs/2)],'bandpass'); %8th order butterworth filter between 10-70Hz

filt_data_transf = filtfilt(b,a,filt_datatemp');
filt_data = filt_data_transf';
disp('...The data is band-pass filtered between 10-70Hz...')

% step 2: Teager energy
clear E1_all
E1_all(:,1) = zeros(size(filt_data,1),1);
E1_all(:,size(filt_data,2)) = zeros(size(filt_data,1),1);

h=waitbar(0,'Calculating the Teager energy...');
for chan = 1:size(filt_data,1)
    for i=2:size(filt_data,2)-1
        E1_all(chan,i) = filt_data(chan,i)^2 - filt_data(chan,i-1)*filt_data(chan,i+1);
    end
    waitbar(chan/size(filt_data,1),h,'Calculating the Teager energy...');
end
close(h)

% step 3: upslope and downslope measure --> first derivative
filt_deriv_data = diff(filt_data,1,2);
shiftU = 10;
U_all = [zeros(size(filt_deriv_data,1),shiftU) filt_deriv_data(:,1:end-shiftU)];
shiftD = 10;
D_all = [-1*filt_deriv_data(:,shiftD:end) zeros(size(filt_deriv_data,1),shiftD-1)];
disp('...The upslope and downslope measure is calculated...')

% step 4: mahalanobis distance
% step 4A: vector of E, U, D
h=waitbar(0,'Constructing vector of Energy, up- and downslope...');
V = NaN(size(filt_data,1),3,size(U_all,2));
for chan=1:size(filt_data,1)
    U = U_all(chan,:);
    D = D_all(chan,:);
    E1 = E1_all(chan,1:end-1);
    
    % per tijdssample heb ik 3 waardes:
    V(chan,:,:) = [E1;U;D]; %V = [channels x E,U,D x time]
    waitbar(chan/size(filt_data,1),h,'Constructing vector of Energy, up- and downslope...');
end
close(h)

% step 4B: mean and covariance of variables
meanV = mean(mean(V,3),1);

covEU = cov(squeeze(V(:,1,:)), squeeze(V(:,2,:)));
covED = cov(squeeze(V(:,1,:)), squeeze(V(:,3,:)));
covUD = cov(squeeze(V(:,2,:)), squeeze(V(:,3,:)));
covV = [covEU(1,1), covEU(1,2), covED(1,2);...
    covEU(2,1), covEU(2,2), covUD(1,2);...
    covED(2,1), covUD(2,1), covUD(2,2)];
disp('...The covariance matrix is calculated...')

% step 4C: Mahalanobis distance
M = zeros(size(V,1),size(V,3));
h=waitbar(0,'Calculating the Mahalanobis distance...');
for chan = 1:size(V,1)
    for t=1:size(V,3)
        M(chan,t) = sqrt((V(chan,:,t)-meanV)*covV*(V(chan,:,t)-meanV)'); % moet 1 scalar worden voor elke [chan] x [t]
    end
    waitbar(chan/size(V,1),h,'Calculating the Mahalanobis distance...');
end
close(h)
disp('The Mahalanobis distance is calculated!') % --> vector = [channels x timesamples]

% % step 4D: delete values of M when electrode is stimulated! --> during
% artefact removal is this signal changed to 0, so no spikes can be found
% here

% for chan=1:size(V,1)
%
%     % determine in which stimulations the electrode is stimulated
%     stimchans = ceil(find([SPESconfig(:).stimulus.stimnum] == keepchan(chan))/2);
%     % for each stimulation, return M to median
%     for stimps = 1:size(stimchans,2)
%         M(chan,SPESconfig.stimulus(stimchans(stimps)).startsamp(1):SPESconfig.stimulus(stimchans(stimps)).startsamp(end)+4*fs)=median(M(chan,:));
%     end
%
% end
%
% step 5: set threshold for Mahalanobis distance to be a spike
% Mmin = 10000;
% Mmax = 100000000;
% dif = round(fs*0.045); % minimale piek = 20ms (=40samp), dus er moet minimaal 40samp tussen zitten-->
% maar dan nog veel FP detecties, dus verzetten: een piek duurt 20-70 ms =
% gemiddeld 45 ms ==
% clear Mrange

%% normalize M

norm_M = M./(sqrt(sum(M.^2,2)));

%% calculate parameter estimates

disp('Calculating parameter estimates takes some time')
Pharmat_norm = NaN(size(norm_M,1),3);
for chan=1:size(filt_data,1)
    
    Pharmat_norm(chan,:) = gevfit(norm_M(chan,:)); % returns 95% confidence intervals for parameter estimates
end

Pharmatall_norm = gevfit(norm_M(:));

disp('Calculation is done!')

%     Mthresh(j) = Pharmat(3) +SD*Pharmat(2); %2SD
%     
%     % per kanaal, locatie van M's en de waardes van M's die binnen grenzen zitten
%     Mrange(j).chan = chan;
%     Mrange(j).locs = find(M(chan,:)>Mthresh(j));
%     Mrange(j).val = M(chan,Mrange(j).locs);
%     Mrange(j).thresh = Mthresh(j);
%     
%     % per kanaal nu opdelen in groepjes
%     Mrange(j).difflocs = diff(Mrange(j).locs);
%     % dit is iedere keer de laatste sample van een groepje :)
%     Mrange(j).lastsampgroups = [Mrange(j).locs(diff(Mrange(j).locs)>dif) Mrange(j).locs(end)];
%     % en een eerste sample van een groepje
%     Mrange(j).firstsampgroups = [Mrange(j).locs(1) Mrange(j).locs(find(diff(Mrange(j).locs)>dif)+1)];
%     
%     % groups{1,:} = samplenummers met Mrange
%     %     Mrange(j).groups{1,1} = Mrange(j).locs(1):Mrange(j).lastsampgroups(1);
%     % groups{2,:} = Mvalues bij samplenummers
%     %     Mrange(j).groups{2,1} = M(chan,Mrange(j).locs:Mrange(j).lastsampgroups(1));
%     
%     for i=1:size(Mrange(j).lastsampgroups,2)
%         % samplenummers met Mrange
%         Mrange(j).groups{1,i} = Mrange(j).firstsampgroups(i):Mrange(j).lastsampgroups(i);
%         % Mvalues bij samplenummers
%         Mrange(j).groups{2,i} = M(chan,Mrange(j).firstsampgroups(i):Mrange(j).lastsampgroups(i));
%         % uV values bij samplenummers
%         Mrange(j).groups{3,i} = data_stimfree(chan,Mrange(j).firstsampgroups(i):Mrange(j).lastsampgroups(i));
%         % sample bij max Mvalue
%         Mrange(j).groups{4,i} = Mrange(j).groups{1,i}(Mrange(j).groups{2,i}==max(Mrange(j).groups{2,i}));
%         % max Mvalue
%         Mrange(j).groups{5,i} = max(Mrange(j).groups{2,i});
%         % location with max uV value
%         Mrange(j).groups{6,i} = Mrange(j).groups{1,i}(find(Mrange(j).groups{3,i} == max(Mrange(j).groups{3,i}),1,'first'));
%         % max uV value
%         Mrange(j).groups{7,i} = max(Mrange(j).groups{3,i});
%         
%         % remove all saturated results (these are mainly stimulation effect)
%         % and spikes with an amplitude < 0 uV
%         if Mrange(j).groups{7,i} < 0 || Mrange(j).groups{7,i} >2900
%             Mrange(j).groups{6,i} = [];
%         end
%     end
%     
%     
%     j=j+1;
% end

end