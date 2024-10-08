%% plotting a figure with the brain and electrodes on the brain
% based on script in
% https://github.com/MultimodalNeuroimagingLab/mnl_ccepBids/scripts/makeFig1A_plotMNI.m

%% first run ccepSp03_analysis_ERs_PS_spikes.m
close all
clc

%% add electrodes.tsv

for nSubj = 1:size(dataBase,2)
    
    % convert electrode positions into a matrix
    if iscell(dataBase(nSubj).tb_electrodes.x)
    
        elecmatrix = NaN(size(dataBase(nSubj).tb_electrodes,1),3);
        for ll = 1:size(dataBase(nSubj).tb_electrodes,1)
            if ~isequal(dataBase(nSubj).tb_electrodes.x{ll},'n/a')
                elecmatrix(ll,:) = [str2double(dataBase(nSubj).tb_electrodes.x{ll}), ...
                    str2double(dataBase(nSubj).tb_electrodes.y{ll}),...
                    str2double(dataBase(nSubj).tb_electrodes.z{ll})];
            end
        end
    else
        
        elecmatrix = [dataBase(nSubj).tb_electrodes.x,...
            dataBase(nSubj).tb_electrodes.y,...
            dataBase(nSubj).tb_electrodes.z];
    end
    
    dataBase(nSubj).elecmatrix = elecmatrix; %#ok<SAGROW> 
end

disp('All electrodes positions are converted to a matrix.')

% housekeeping
clear elecmatrix ll nSubj

%% load mni305 pial
% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.proj_diroutput,'freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

% housekeeping
clear FSsubjectsdir 

%% plot subject electrodes on mni brain

close all

nSubj = 8; % in Fig1A of the article number 8 is used
mkrsz = 30;

% get electrodes info
elCoords = dataBase(nSubj).elecmatrix;

% get hemisphere for each electrode
hemi = dataBase(nSubj).tb_electrodes.hemisphere;

if any(contains(fieldnames(dataBase(nSubj).tb_electrodes),'ied'))
    IEDelec = strcmp(dataBase(nSubj).tb_electrodes.ied,'yes');
else
    IEDelec = false(size(dataBase(nSubj).tb_electrodes,1),1);
end

% set the view for the correct hemisphere
if isequal(hemi{1},'L')
    g.faces = Lmnipial_face+1; % correct for zero index
    g.vertices = Lmnipial_vert;
    v_d = ([270 0]);
elseif isequal(hemi{1},'R')
    g.faces = Rmnipial_face+1; % correct for zero index
    g.vertices = Rmnipial_vert;
    v_d = ([90 0]);
end

% make the electrodes pop out of the brain cortex
a_offset = .2*max(abs(elCoords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = elCoords+repmat(a_offset,size(elCoords,1),1);      

figure
ieeg_RenderGifti(g);

hold on

% all electrodes:
plot3(els(:,1), ...
    els(:,2), ...
    els(:,3), ...
    '.','Color', 'k','MarkerSize',mkrsz)

% stimulated and response electrode:
stimelec = contains(dataBase(nSubj).tb_electrodes.name,{'C23','C24'});
respelec = contains(dataBase(nSubj).tb_electrodes.name,'C16');
plot3(els(stimelec==1,1), els(stimelec==1,2), els(stimelec==1,3), ...
    '.','Color', [184/256, 26/256, 93/256],'MarkerSize',0.8*mkrsz)
plot3(els(respelec==1,1), els(respelec==1,2), els(respelec==1,3), ...
    '.','Color', [0/256, 156/256, 180/256],'MarkerSize',0.8*mkrsz)

% electrodes with IEDs
plot3(els(IEDelec==1,1), ...
    els(IEDelec==1,2), ...
    els(IEDelec==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',0.3*mkrsz) % color is almost white, otherwise it is not saved correctly

hold off

ieeg_viewLight(v_d(1),v_d(2))

figureName = sprintf('%s/fig1a_brain_%s',...
    myDataPath.figures,dataBase(nSubj).sub_label);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)

fprintf('Figure is saved as .png in \n %s \n',figureName)

% housekeeping
clear a_offset ans elCoords els figureName g hemi IEDelec Lmnipial_face 
clear Lmnipial_vert mkrsz respelec Rmnipial_face Rmnipial_vert nSubj v_d stimelec

%% end of script