clc; clear; close all

patchRecord_beh = crawlForFiles;
dataLoc = '~/Downloads/v4-7a/patch';

%% get/parse all data
for ii=1:length(patchRecord_beh)
    exptName = patchRecord_beh(ii).path(1+find(patchRecord_beh(ii).path=='/',1,'last'):(find(patchRecord_beh(ii).path=='.',1,'last')-1));
    disp(exptName)
    clearvars exptProps patchProps trialProps
    if ~exist([dataLoc '/' exptName '_beh.mat'],'file')
        disp('... getting remote')
        load(patchRecord_beh(ii).path)
        behav = [behav.trialData];

        if sum([behav.analyze])>0
            exptProps.name = exptName;
            exptProps.nTrial = length(behav);
            exptProps.nAnalyze = sum([behav.analyze]);
    
            % remove all the breaks for now
            behav(~[behav.analyze]) = [];
    
            exptProps.patch_n = behav(1).fixation.n;
            exptProps.patch_T_reach = behav(1).fixation.reach;
            exptProps.patch_T_stay = behav(1).fixation.stay;
            exptProps.reward_prob = behav(1).reward.prob;
            exptProps.reward_amount = behav(1).reward.time;
            exptProps.stim_x = behav(1).stim.x;
            exptProps.stim_y = behav(1).stim.y;
            exptProps.stim_s = behav(1).stim.s;
    
            for pp = 1:behav(1).fixation.n
                patchProps(pp).fix_color = behav(1).fixation.color(pp,:);
                patchProps(pp).fix_x = behav(1).fixation.x(pp);
                patchProps(pp).fix_y = behav(1).fixation.y(pp);
                patchProps(pp).stim_set = behav(1).stim.set(pp);
                patchProps(pp).stim_contrast = behav(1).stim.contrast(pp);
                patchProps(pp).sigm_beta = [behav(1).stim.on_min(pp) behav(1).stim.on_max(pp) behav(1).stim.on_mid(pp) behav(1).stim.on_slp(pp)];
            end
            
            for tt=1:length(behav)
                trialProps(tt).id = behav(tt).id;
                trialProps(tt).patch = behav(tt).patch;
                trialProps(tt).ISI = behav(tt).ISI;
                trialProps(tt).stim_num = behav(tt).stim.num;
                trialProps(tt).stim_ontime = behav(tt).stim.on;
            end
            save([dataLoc '/' exptName '_beh.mat'],'exptProps','patchProps','trialProps')
        end
    else
        load([dataLoc '/' exptName '_beh.mat'],'exptProps','patchProps','trialProps')
    end
    
    if exist('exptProps','var')
        patchRecord_beh(ii).name = exptName;
        patchRecord_beh(ii).nTrials = exptProps.nTrial;
        patchRecord_beh(ii).nAnalyze = exptProps.nAnalyze;
        patchRecord_beh(ii).nPatch = exptProps.patch_n;
        patchRecord_beh(ii).stimSets = [patchProps.stim_set];
        patchRecord_beh(ii).patchFixCount = histcounts([trialProps.patch]);
    end
end
patchRecord_beh(cellfun(@isempty,{patchRecord_beh.name})) = [];
save('data/patchRecord_beh.mat','patchRecord_beh')

%% helpers
function patchRecord = crawlForFiles
    filelist = cell(1);
    rootDir = '/Volumes/colada/Ram/data/neural/v4-7a';
    dirs = dir([rootDir '/25*']);
    for ii=1:length(dirs)
        files_grid = dir([rootDir '/' dirs(ii).name '/*zippy*patch*mat*']);
        filelist = [filelist cellfun(@(x) strrep([rootDir '/' dirs(ii).name '/' x],'.nev',''),{files_grid.name},'UniformOutput',false)];
        % if length(files_nev) ~= length(files_grid)
        %     for jj=1:length(files_nev)
        %         files_nev.name(jj)
        %     end
        % end
    end
    filelist(1) = [];
    patchRecord = cell2struct(filelist,'path',1);
end
