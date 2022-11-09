%% plotting a figure with the brain and electrodes on the brain
% based on script in
% https://github.com/MultimodalNeuroimagingLab/mnl_ccepBids/scripts/makeFig1A_plotMNI.m

%% first run ccepSp03_analysis_ERs_PS_spikes.m

%% add electrodes.tsv

for subj = 1:size(dataBase,2)
    
    % convert electrode positions into a matrix
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

disp('All electrodes positions are converted to a matrix.')

% housekeeping
clear both CCEPelec CCEPmat CCEPnum elecmatrix ERSPCCEPcomb ERSPelec 
clear ERSPmat ERSPnum IEDmat Jaccard ll numstim subj

%% load mni305 pial
% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

% housekeeping
clear FSsubjectsdir

%% add all electrodes left or right hemisphere into one 
% variable: allmni_coords 

all_mnicoords = [];
all_hemi = [];
all_IED = [];
all_ERSPnormelec = [];
all_CCEPnormelec = [];
% all_percCCEPERSP = [];
all_Jaccard = [];

for subj = 1:size(dataBase,2)

    all_mnicoords = [all_mnicoords; dataBase(subj).elecmatrix]; %#ok<AGROW>
    all_hemi = [all_hemi; dataBase(subj).tb_electrodes.hemisphere]; %#ok<AGROW>
    all_IED = [all_IED; dataBase(subj).IEDmat];  %#ok<AGROW>
    all_ERSPnormelec = [all_ERSPnormelec; dataBase(subj).ERSPmatnorm]; %#ok<AGROW> 
    all_CCEPnormelec = [all_CCEPnormelec; dataBase(subj).CCEPmatnorm]; %#ok<AGROW> 
%     all_percCCEPERSP = [all_percCCEPERSP, dataBase(subj).percCCEPERSP]; %#ok<AGROW> 
       all_Jaccard = [all_Jaccard, dataBase(subj).Jaccard]; %#ok<AGROW>
end

% housekeeping
clear subj

%% set characteristics in figures

cmap = parula(101); % 0-1 with steps of 0.01
mkrsz = 30;

% view left and right
v_dL = [270 0];
v_dR = [96 6];

% make sure electrodes pop out: Left and Right
a_offsetLview = .1*max(abs(all_mnicoords(:,1)))*[cosd(v_dL(1)-90)*cosd(v_dL(2)) sind(v_dL(1)-90)*cosd(v_dL(2)) sind(v_dL(2))];
all_elecLview = all_mnicoords+repmat(a_offsetLview,size(all_mnicoords,1),1);      

a_offsetRview = .1*max(abs(all_mnicoords(:,1)))*[cosd(v_dR(1)-90)*cosd(v_dR(2)) sind(v_dR(1)-90)*cosd(v_dR(2)) sind(v_dR(2))];
all_elecRview = all_mnicoords+repmat(a_offsetRview,size(all_mnicoords,1),1);      

% brain surface
gL.faces = Lmnipial_face+1;
gL.vertices = Lmnipial_vert;
gL = gifti(gL);

gR.faces = Rmnipial_face+1;
gR.vertices = Rmnipial_vert;
gR = gifti(gR);

% housekeeping
clear Rmnipial_face Rmnipial_vert Lmnipial_face Lmnipial_vert

%% Plot figure with left pial with electrodes in mni space 
% and normalized number of CCEPs

figure;
ieeg_RenderGifti(gL); 

hold on, 
for elec = 1:size(all_elecLview,1)
    if strcmpi(all_hemi(elec),'L')
        
        cmapnum = round(all_CCEPnormelec(elec)*100)+1;
        plot3(all_elecLview(elec,1), all_elecLview(elec,2),all_elecLview(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',mkrsz);

    end
end

plot3(all_elecLview(ismember(all_hemi,'L') & all_IED==1,1), ...
    all_elecLview(ismember(all_hemi,'L') & all_IED==1,2), ...
    all_elecLview(ismember(all_hemi,'L') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256], ...
    'MarkerSize',0.4*mkrsz) % color is almost white

ieeg_viewLight(v_dL(1),v_dL(2))

% title('Normalized number of CCEPs')

figureName = sprintf('%s/fig2_brainL_normCCEP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

% housekeeping
clear ans cmapnum elec figureName

%% Plot figure with right pial with electrodes in mni space
% and normalized number of CCEP

figure
ieeg_RenderGifti(gR); 

hold on
for elec = 1:size(all_elecRview,1)
    if strcmpi(all_hemi(elec),'R')

        cmapnum = round(all_CCEPnormelec(elec)*100)+1;
        plot3(all_elecRview(elec,1), all_elecRview(elec,2),all_elecRview(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',mkrsz);

    end
end

plot3(all_elecRview(ismember(all_hemi,'R') & all_IED==1,1), ...
    all_elecRview(ismember(all_hemi,'R') & all_IED==1,2), ...
    all_elecRview(ismember(all_hemi,'R') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256], ...
    'MarkerSize',0.4*mkrsz)

ieeg_viewLight(v_dR(1),v_dR(2))

% title('Normalized number of CCEPs')

figureName = sprintf('%s/fig2_brainR_normCCEP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

clear ans cmapnum elec figureName

%% Plot figure with left and right pial with electrodes in mni space 
% and normalized number of ERSPs

figure;
ieeg_RenderGifti(gL); 

hold on, 
for elec = 1:size(all_elecLview,1)
    if strcmpi(all_hemi(elec),'L')
        
        cmapnum = round(all_ERSPnormelec(elec)*100)+1;
        plot3(all_elecLview(elec,1), all_elecLview(elec,2),all_elecLview(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',mkrsz);

    end
end

plot3(all_elecLview(ismember(all_hemi,'L') & all_IED==1,1), ...
    all_elecLview(ismember(all_hemi,'L') & all_IED==1,2), ...
    all_elecLview(ismember(all_hemi,'L') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256], ...
    'MarkerSize',0.4*mkrsz)

ieeg_viewLight(v_dL(1),v_dL(2))

% title('Normalized number of ERSPs')

figureName = sprintf('%s/fig2_brainL_normERSP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

% housekeeping
clear ans cmapnum elec figureName

%% Plot figure with right pial with electrodes in mni space
% and normalized ERSP

figure
ieeg_RenderGifti(gR); 

hold on
for elec = 1:size(all_elecRview,1)
    if strcmpi(all_hemi(elec),'R')

        cmapnum = round(all_ERSPnormelec(elec)*100)+1;
        plot3(all_elecRview(elec,1), all_elecRview(elec,2),all_elecRview(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',mkrsz);

    end
end

plot3(all_elecRview(ismember(all_hemi,'R') & all_IED==1,1), ...
    all_elecRview(ismember(all_hemi,'R') & all_IED==1,2), ...
    all_elecRview(ismember(all_hemi,'R') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',0.4*mkrsz)

ieeg_viewLight(v_dR(1),v_dR(2))

% title('Normalized number of ERSPs')

figureName = sprintf('%s/fig2_brainR_normERSP',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

% housekeeping
clear ans cmapnum elec figureName

%% Plot figure with left pial with electrodes in mni space 
% and Jaccard index

figure;
ieeg_RenderGifti(gL); 

hold on, 
for elec = 1:size(all_elecLview,1)
    if strcmpi(all_hemi(elec),'L')
        
        cmapnum = round(all_Jaccard(elec)*100)+1;
        plot3(all_elecLview(elec,1), all_elecLview(elec,2),all_elecLview(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',mkrsz);

    end
end

plot3(all_elecLview(ismember(all_hemi,'L') & all_IED==1,1), ...
    all_elecLview(ismember(all_hemi,'L') & all_IED==1,2), ...
    all_elecLview(ismember(all_hemi,'L') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',0.4*mkrsz)

ieeg_viewLight(v_dL(1),v_dL(2))

% title('Jaccard similarity coefficient')

figureName = sprintf('%s/fig2_brainL_Jaccard',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

% housekeeping
clear ans cmapnum elec figureName

%% Plot figure with right pial with electrodes in mni space
% and Jaccard index

figure
ieeg_RenderGifti(gR); 

hold on
for elec = 1:size(all_elecRview,1)
    if strcmpi(all_hemi(elec),'R')

        cmapnum = round(all_Jaccard(elec)*100)+1;
        plot3(all_elecRview(elec,1), all_elecRview(elec,2),all_elecRview(elec,3),...
            '.','Color',cmap(cmapnum,:),'MarkerSize',mkrsz);

    end
end

plot3(all_elecRview(ismember(all_hemi,'R') & all_IED==1,1), ...
    all_elecRview(ismember(all_hemi,'R') & all_IED==1,2), ...
    all_elecRview(ismember(all_hemi,'R') & all_IED==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',0.4*mkrsz)

ieeg_viewLight(v_dR(1),v_dR(2))

% title('Jaccard similarity coefficient')

figureName = sprintf('%s/fig2_brainR_Jaccard',...
    myDataPath.Figures);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

% housekeeping
clear ans cmapnum elec figureName


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

% housekeeping
clear a_offsetLview a_offsetRview all_CCEPnormelec all_elecLview all_elecRview
clear all_ERSPnormelec all_hemi all_IED all_Jaccard all_mnicoords figureName gL
clear gR h limits mkrsz v_dL v_dR
