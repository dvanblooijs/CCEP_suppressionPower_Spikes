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
    
    dataBase(subj).elecmatrix = elecmatrix; %#ok<SAGROW> 
end

disp('All electrodes.tsv are loaded')

%% load mni305 pial
% Freesurfer subjects directory
FSsubjectsdir = fullfile(myDataPath.dataPath,'derivatives','freesurfer');

% load mni305 pial
[Lmnipial_vert,Lmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
[Rmnipial_vert,Rmnipial_face] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));

%% plot subject electrodes on mni brain

close all

subj = 8; % in Fig1A of the article number 8 is used

% get electrodes info
elCoords = dataBase(subj).elecmatrix;

% get hemisphere for each electrode
hemi = dataBase(subj).tb_electrodes.hemisphere;

if any(contains(fieldnames(dataBase(subj).tb_electrodes),'ied'))
    IEDelec = strcmp(dataBase(subj).tb_electrodes.ied,'yes');
else
    IEDelec = false(size(dataBase(subj).tb_electrodes,1),1);
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
    '.','Color', 'k','MarkerSize',30)

% stimulated and response electrode:
stimelec = contains(dataBase(subj).tb_electrodes.name,{'C23','C24'});
respelec = contains(dataBase(subj).tb_electrodes.name,'C16');
plot3(els(stimelec==1,1), els(stimelec==1,2), els(stimelec==1,3), ...
    '.','Color', [184/256, 26/256, 93/256],'MarkerSize',25)
plot3(els(respelec==1,1), els(respelec==1,2), els(respelec==1,3), ...
    '.','Color', [0/256, 156/256, 180/256],'MarkerSize',25)

% electrodes with IEDs
plot3(els(IEDelec==1,1), ...
    els(IEDelec==1,2), ...
    els(IEDelec==1,3), ...
    '.','Color', [240/256, 240/256, 240/256],'MarkerSize',7)

hold off

ieeg_viewLight(v_d(1),v_d(2))

figureName = sprintf('%s/fig1a_brain_%s',...
    myDataPath.Figures,dataBase(subj).sub_label);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
print('-painters','-depsc',figureName)

fprintf('Figure is saved as .png and .eps in \n %s \n',figureName)