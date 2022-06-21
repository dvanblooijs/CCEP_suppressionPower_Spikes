function detIED = findCortSpikes(data_stimfree,fs,Mvalues,Pharmat,thresh_t_dif,thresh_SD)

dif = round(thresh_t_dif * fs);
detIED = cell(1,size(Mvalues,1));

for chan = 1:size(Mvalues,1)
    
    Mthresh = (Pharmat(chan,3) + thresh_SD*Pharmat(chan,2));
    
    % per kanaal, locatie van M's en de waardes van M's die binnen grenzen zitten
    locs = find(Mvalues(chan,:)>Mthresh);
    
    if ~isempty(locs)
        % dit is iedere keer de laatste sample van een groepje
        lastsampgroups = [locs(diff(locs)>dif) locs(end)];
        % en een eerste sample van een groepje
        firstsampgroups = [locs(1) locs(find(diff(locs)>dif)+1)];
        
        sampnums = cell(1,size(lastsampgroups,2)); uV_sampnums = cell(1,size(lastsampgroups,2));
        sampnum_minuV = NaN(1,size(lastsampgroups,2)); minuV = NaN(1,size(lastsampgroups,2));
        for i=1:size(lastsampgroups,2)
            % samplenummers met Mrange
            sampnums{1,i} = firstsampgroups(i):lastsampgroups(i);
            % uV values bij samplenummers
            uV_sampnums{1,i} = data_stimfree(chan,firstsampgroups(i):lastsampgroups(i));
            % location with min uV value (spikes are negative peaks)
            sampnum_minuV(1,i) = sampnums{1,i}(find(uV_sampnums{1,i} == min(uV_sampnums{1,i}),1,'first'));
            % min uV value
            minuV(1,i) = max(uV_sampnums{1,i});
            
            % remove all saturated results (these are mainly stimulation effect)
            % and spikes with an amplitude < 0 uV
            if minuV(1,i) < -2900
                sampnum_minuV(1,i) = NaN;
                minuV(1,i)  = NaN;
            end
        end
        
    end
    detIED{chan} = sampnum_minuV;
end


end