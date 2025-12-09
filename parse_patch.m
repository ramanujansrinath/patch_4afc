clc; clear; close all

patchRecord = crawlForFiles;
dataLoc = '~/Downloads/v4-7a/patch';

%% get/parse all data
for ii=1:length(patchRecord)
    exptName = patchRecord(ii).path(1+find(patchRecord(ii).path=='/',1,'last'):end);
    disp(exptName)
    if ~exist([dataLoc '/' exptName '_dense.mat'],'file')
        disp('... getting remote')
        load([patchRecord(ii).path '_trialinfo'])
        load(patchRecord(ii).path)
        
        clearvars params resp resp_base
        for tt=1:length(trials)
            if trials{tt}.good
                params(tt).set = trials{tt}.set;
                params(tt).num = trials{tt}.num;
                params(tt).x = trials{tt}.x;
                params(tt).y = trials{tt}.y;
                params(tt).s = trials{tt}.s;
                params(tt).good = trials{tt}.good;
                params(tt).spikes = trials{tt}.spikes(:,[1 3]);
                params(tt).spikes(:,2) = params(tt).spikes(:,2)-trials{tt}.stimstart;
    
                sp = params(tt).spikes(:,1);
                t_sp = params(tt).spikes(:,2);
    
                respDur = trials{tt}.stimend-trials{tt}.stimstart;
                respLims = [0.05 respDur];
                baseLims = [-0.1 0];
                
                sp_idx = t_sp > baseLims(1) & t_sp < respLims(2); % find all spikes in the window
                sp_ch = sp(sp_idx); % the channels the spikes were on
                [sp_ch,idx] = sort(sp_ch); % sort by channel
                sp_ts = t_sp(sp_idx); % get spike timestamps (in s)
                sp_ts = sp_ts(idx); % sort in the same way as channels
                
                params(tt).resp      = histcounts(sp_ch(sp_ts > respLims(1)),0.5:1:128.5)/respDur;
                params(tt).resp_base = histcounts(sp_ch(sp_ts < baseLims(2)),0.5:1:128.5)/0.1;
                
                resp(tt,:) = params(tt).resp;
                resp_base(tt,:) = params(tt).resp_base;
            end
        end

        save([patchRecord(ii).path '_dense'],'params','resp','resp_base')
        save([dataLoc '/' exptName '_dense.mat'],'params','resp','resp_base')
    else
        load([dataLoc '/' exptName '_dense.mat'],'params','resp','resp_base')
    end

    patchRecord(ii).name = exptName;
    patchRecord(ii).nTrials = length(params);
    patchRecord(ii).nCorrect = sum([params.good]);
    patchRecord(ii).set = unique([params.set]);
    patchRecord(ii).num = unique([params.num]);
    patchRecord(ii).x = unique([params.x]);
    patchRecord(ii).y = unique([params.y]);
    patchRecord(ii).s = unique([params.s]);
    patchRecord(ii).rfPos = [];
    
    % if actual rf map expt, save the RFs
    if length(unique([params.x]))>1
        if ~exist([dataLoc '/' exptName '_rf.mat'],'file')
            disp('... saving rf')
            rfv4 = rfmap_save(params,resp,eid);
            save([dataLoc '/' exptName '_rf.mat'],'rfv4')
        else
            load([dataLoc '/' exptName '_rf.mat'],'rfv4')
        end
        % rfmap_plot(rfv4)
        patchRecord(ii).rfPos = rfv4.pos;
    end
end
% save('data/mapRecord.mat','mapRecord')
% mapRecord = mapRecord([mapRecord.nTrials]>100);

%% helpers
function patchRecord = crawlForFiles
    filelist = cell(1);
    rootDir = '/Volumes/colada/Ram/data/neural/v4-7a';
    dirs = dir([rootDir '/25*']);
    for ii=1:length(dirs)
        files_grid = dir([rootDir '/' dirs(ii).name '/*patch*nev*']);
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
