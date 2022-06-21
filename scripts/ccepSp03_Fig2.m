%% plotting a figure with the brain and electrodes on the brain
% based on script in
% https://github.com/MultimodalNeuroimagingLab/mnl_ccepBids/scripts/makeFig1A_plotMNI.m



%% load first ccepSp03_analysis_ERs_PS_spikes.m

%% patient settings

files = dir(myDataPath.dataPath);
idx_subj = contains({files(:).name},'sub-');
files_subj = files(idx_subj);
cfg = struct([]);

for subj = 1:size(files_subj,1)

    cfg(subj).sub_labels = files_subj(subj).name;

    files = dir(fullfile(files_subj(subj).folder,files_subj(subj).name));
    idx_ses = contains({files(:).name},'ses-');
    files_ses = files(idx_ses);

    cfg(subj).ses_label = files_ses(1).name;

    cfg(subj).task_label = 'task-SPESclin';

    files = dir(fullfile(files_ses(1).folder,files_ses(1).name,'ieeg'));
    idx_eeg = contains({files(:).name},'.eeg');
    files_eeg = files(idx_eeg);

    for run = 1:size(files_eeg,1)
        runTemp = extractBetween(files_eeg(run).name,'run-','_ieeg');
        cfg(subj).run_label{run} = ['run-', runTemp{1}];
    end
end

%% add electrodes.tsv

for subj = 1:size(cfg,2)
    
    % preprocess electrode positions
    if iscell(dataBase(subj).tb_electrodes.x)
        elecmatrix = NaN(size(dataBase(subj).tb_electrodes,1),3);
        for ll = 1:size(dataBase(subj).tb_electrodes,1)
            if ~isequal(dataBase(subj).tb_electrodes.x{ll},'n/a')
                elecmatrix(ll,:) = [str2double(dataBase(subj).tb_electrodes.x{ll}), ...
                    str2double(dataBase(subj).tb_electrodes.y{ll}),...
                    str2double(dataBase(subj).tb_electrodes.z{ll})];
            end
        end
    else
        elecmatrix = [dataBase(subj).tb_electrodes.x,...
            dataBase(subj).tb_electrodes.y,...
            dataBase(subj).tb_electrodes.z];
    end

    if any(contains(fieldnames(dataBase(subj).tb_electrodes),'ied'))
        IEDmat = strcmp(dataBase(subj).tb_electrodes.ied,'yes');
    else
        IEDmat = false(size(dataBase(subj).tb_electrodes,1),1);
    end
    
    dataBase(subj).elecmatrix = elecmatrix; %#ok<SAGROW> 
    dataBase(subj).IEDmat = IEDmat;  %#ok<SAGROW> 

    % normalize CCEPnumber and ERSPnumber according to total number of
    % electrodes
    CCEPmat = dataBase(subj).CCEPmat; % [stim x chan]
    numstim = sum(isnan(CCEPmat),1); % how often each channel is stimulated
    CCEPelec = sum(CCEPmat,1,'omitnan')./(repmat(size(CCEPmat,1),1,size(CCEPmat,2))-numstim); % normalized CCEPs
    dataBase(subj).CCEPmatnorm = CCEPelec'; %#ok<SAGROW>

    ERSPmat = dataBase(subj).ERSPmat; %[stim x chan]
    numstim = sum(isnan(ERSPmat),1); % how often each channel is stimulated
    ERSPelec = sum(ERSPmat,1,'omitnan')./(repmat(size(ERSPmat,1),1,size(ERSPmat,2))-numstim); % normalized ERSPs
    dataBase(subj).ERSPmatnorm = ERSPelec'; %#ok<SAGROW>

    % calculate response electrodes with both CCEP and ERSP (Jaccard index)
    ERSPCCEPcomb = ERSPmat + CCEPmat;
    CCEPnum = sum(CCEPmat,1,'omitnan');
    ERSPnum = sum(ERSPmat,1,'omitnan');
    both = sum(ERSPCCEPcomb==2);
    Jaccard = both./(CCEPnum+ERSPnum-both);
    Jaccard(isnan(Jaccard)) = 1; % because it is NaN when /0 --> 1, because there is total similarity when no events were detected in both situations
    dataBase(subj).Jaccard = Jaccard; %#ok<SAGROW>
end

disp('All electrodes.tsv are loaded')

%% load mni305 pial
% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

%% add all electrodes left or right hemisphere into one 
% variable: allmni_coords 

allmni_coords = [];
allmni_labels = [];
all_hemi = [];
all_IED = [];
all_ERSPelec = [];
all_CCEPelec = [];
all_percCCEPERSP = [];
all_Jaccard = [];

for subj = 1:size(dataBase,2)

    allmni_coords = [allmni_coords; dataBase(subj).elecmatrix]; %#ok<AGROW>
    all_hemi = [all_hemi; dataBase(subj).tb_electrodes.hemisphere]; %#ok<AGROW>
    all_IED = [all_IED; dataBase(subj).IEDmat];  %#ok<AGROW>
    all_ERSPelec = [all_ERSPelec; dataBase(subj).ERSPmatnorm]; %#ok<AGROW> 
    all_CCEPelec = [all_CCEPelec; dataBase(subj).CCEPmatnorm]; %#ok<AGROW> 
%     all_percCCEPERSP = [all_percCCEPERSP, dataBase(subj).percCCEPERSP]; %#ok<AGROW> 
       all_Jaccard = [all_Jaccard, dataBase(subj).Jaccard]; %#ok<AGROW>
end

%% set colorscale

cmap = parula(101); % 0-1 with steps of 0.01

%% Plot figure with left and right pial with electrodes in mni space 
% and normalized number of CCEPs

v_d = [270 0];

figure;
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
tH = ieeg_RenderGifti(gl); %#ok<NASGU>
% 
% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

hold on, 
for elec = 1:size(els,1)
    if strcmpi(all_hemi(elec),'L')
        
        cmapnum = round(all_CCEPelec(elec)*100)+1;
        plot3(els(elec,1), els(elec,2),els(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',20);

    end
end

plot3(els(ismember(all_hemi,'L') & all_IED==1,1), ...
    els(ismember(all_hemi,'L') & all_IED==1,2), ...
    els(ismember(all_hemi,'L') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

% set(tH,'FaceAlpha',.5) % make transparent
ieeg_viewLight(v_d(1),v_d(2))

% h=colorbar;
% limits = h.Limits;
% h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
% h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'}; 

% title('Normalized number of CCEPs')

figureName = sprintf('%s/fig2_brainL_normCCEP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)


%% Plot figure with right pial with electrodes in mni space
% and normalized CCEP

v_d = [96 6];

figure
gr.faces = Rmnipial_face+1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);
ieeg_RenderGifti(gr); 

% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

hold on
for elec = 1:size(els,1)
    if strcmpi(all_hemi(elec),'R')

        cmapnum = round(all_CCEPelec(elec)*100)+1;
        plot3(els(elec,1), els(elec,2),els(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',20);

    end
end

plot3(els(ismember(all_hemi,'R') & all_IED==1,1), ...
    els(ismember(all_hemi,'R') & all_IED==1,2), ...
    els(ismember(all_hemi,'R') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

% set(tH,'FaceAlpha',.5) % make transparent
ieeg_viewLight(v_d(1),v_d(2))

% h=colorbar;
% limits = h.Limits;
% h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
% h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'}; %NOG IETS AANPASSEN!!
% 
% title('Normalized number of CCEPs')

figureName = sprintf('%s/fig2_brainR_normCCEP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

%% Plot figure with left and right pial with electrodes in mni space 
% and normalized number of ERSPs

v_d = [270 0];

figure;
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
ieeg_RenderGifti(gl); 
% 
% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

hold on, 
for elec = 1:size(els,1)
    if strcmpi(all_hemi(elec),'L')
        
        cmapnum = round(all_ERSPelec(elec)*100)+1;
        plot3(els(elec,1), els(elec,2),els(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',20);

    end
end

plot3(els(ismember(all_hemi,'L') & all_IED==1,1), ...
    els(ismember(all_hemi,'L') & all_IED==1,2), ...
    els(ismember(all_hemi,'L') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

% set(tH,'FaceAlpha',.5) % make transparent
ieeg_viewLight(v_d(1),v_d(2))

% h=colorbar;
% limits = h.Limits;
% h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
% h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'}; 
% 
% title('Normalized number of ERSPs')

figureName = sprintf('%s/fig2_brainL_normERSP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)


%% Plot figure with right pial with electrodes in mni space
% and normalized ERSP

v_d = [96 6];

figure
gr.faces = Rmnipial_face+1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);
ieeg_RenderGifti(gr); 

% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

hold on
for elec = 1:size(els,1)
    if strcmpi(all_hemi(elec),'R')

        cmapnum = round(all_ERSPelec(elec)*100)+1;
        plot3(els(elec,1), els(elec,2),els(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',20);

    end
end

plot3(els(ismember(all_hemi,'R') & all_IED==1,1), ...
    els(ismember(all_hemi,'R') & all_IED==1,2), ...
    els(ismember(all_hemi,'R') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

% set(tH,'FaceAlpha',.5) % make transparent
ieeg_viewLight(v_d(1),v_d(2))

% h=colorbar;
% limits = h.Limits;
% h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
% h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'}; %NOG IETS AANPASSEN!!
% 
% title('Normalized number of ERSPs')

figureName = sprintf('%s/fig2_brainR_normERSP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

%% Plot figure with left and right pial with electrodes in mni space 
% and Jaccard index

v_d = [270 0];

F = figure;
gl.faces = Lmnipial_face+1;
gl.vertices = Lmnipial_vert;
gl = gifti(gl);
tH = ieeg_RenderGifti(gl); %#ok<NASGU>
% 
% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

hold on, 
for elec = 1:size(els,1)
    if strcmpi(all_hemi(elec),'L')
        
        cmapnum = round(all_Jaccard(elec)*100)+1;
        plot3(els(elec,1), els(elec,2),els(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',20);

    end
end

plot3(els(ismember(all_hemi,'L') & all_IED==1,1), ...
    els(ismember(all_hemi,'L') & all_IED==1,2), ...
    els(ismember(all_hemi,'L') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

% set(tH,'FaceAlpha',.5) % make transparent
ieeg_viewLight(v_d(1),v_d(2))

% h=colorbar;
% limits = h.Limits;
% h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
% h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'}; 
% 
% title('Jaccard similarity coefficient')

figureName = sprintf('%s/fig2_brainL_Jaccard',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)


%% Plot figure with right pial with electrodes in mni space
% and Jaccard index

v_d = [96 6];

figure
gr.faces = Rmnipial_face+1;
gr.vertices = Rmnipial_vert;
gr = gifti(gr);
tH = ieeg_RenderGifti(gr); 

% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      

hold on
for elec = 1:size(els,1)
    if strcmpi(all_hemi(elec),'R')

        cmapnum = round(all_Jaccard(elec)*100)+1;
        plot3(els(elec,1), els(elec,2),els(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',20);

    end
end

plot3(els(ismember(all_hemi,'R') & all_IED==1,1), ...
    els(ismember(all_hemi,'R') & all_IED==1,2), ...
    els(ismember(all_hemi,'R') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

% set(tH,'FaceAlpha',.5) % make transparent
ieeg_viewLight(v_d(1),v_d(2))

% h=colorbar;
% limits = h.Limits;
% h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
% h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'}; %NOG IETS AANPASSEN!!

% title('Jaccard similarity coefficient')

figureName = sprintf('%s/fig2_brainR_Jaccard',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

%% plot colorbar
figure,
h=colorbar;
limits = h.Limits;
h.Ticks = limits(1):(limits(2)-limits(1))/10:limits(2);
h.TickLabels = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};

figureName = sprintf('%s/fig2_colorbar',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)



