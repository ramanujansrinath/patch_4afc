clc; clear; close all

mapRecord = crawlForFiles;
dataLoc = '~/Downloads/v4-7a';

%% get/parse all data
for ii=1:length(mapRecord)
    exptName = mapRecord(ii).path(1+find(mapRecord(ii).path=='/',1,'last'):end);
    disp(exptName)
    if ~exist([dataLoc '/' exptName '_dense.mat'],'file')
        disp('... getting remote')
        load([mapRecord(ii).path '_trialinfo'])
        load(mapRecord(ii).path)
        
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

        save([mapRecord(ii).path '_dense'],'params','resp','resp_base')
        save([dataLoc '/' exptName '_dense.mat'],'params','resp','resp_base')
    else
        load([dataLoc '/' exptName '_dense.mat'],'params','resp','resp_base')
    end

    mapRecord(ii).name = exptName;
    mapRecord(ii).nTrials = length(params);
    mapRecord(ii).nCorrect = sum([params.good]);
    mapRecord(ii).set = unique([params.set]);
    mapRecord(ii).num = unique([params.num]);
    mapRecord(ii).x = unique([params.x]);
    mapRecord(ii).y = unique([params.y]);
    mapRecord(ii).s = unique([params.s]);
    mapRecord(ii).rfPos = [];
    
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
        mapRecord(ii).rfPos = rfv4.pos;
    end
end
% save('data/mapRecord.mat','mapRecord')
% mapRecord = mapRecord([mapRecord.nTrials]>100);

%% grid decoding
% eid = [3*ones(1,32) ones(1,64) 3*ones(1,32)];
% v4sites = resp(:,eid==1);
% v4sites_base = resp_base(:,eid==1);
% badTrials = cellfun(@isempty,{params.set});
% v4sites(badTrials,:) = [];
% v4sites_base(badTrials,:) = [];
% params_edit = params;
% params_edit(badTrials) = [];
% sets = unique([params_edit.set]);
% 
% % [params,idx] = sortStruct(params,'num');
% % v4sites = v4sites(idx,:);
% % params = params(idx);
% v4sites_red = pca(v4sites','NumComponents',10);
% 
% stimVals = [mat2vec(repmat(1:5,5,1)) mat2vec(repmat(1:5,1,5))];
% colIds = stimVals([params_edit.num],1);
% shpIds = stimVals([params_edit.num],2);
% 
% [dec_col,r_col] = decodeFn(colIds,v4sites_red);
% [dec_shp,r_shp] = decodeFn(shpIds,v4sites_red);
% 
% figure;
% subplot(221); notBoxPlot(dec_col,colIds); fixPlot(gca,[0 6],[0 6],'col','dec col',1:5,1:5)
% subplot(222); notBoxPlot(dec_shp,shpIds); fixPlot(gca,[0 6],[0 6],'shp','dec shp',1:5,1:5)
% 
% subplot(223); scatter3(v4sites_red(:,1),v4sites_red(:,2),v4sites_red(:,3),20,colIds,'filled')
% subplot(224); scatter3(v4sites_red(:,1),v4sites_red(:,2),v4sites_red(:,3),20,shpIds,'filled')


%% decoding for non-map files
% v4sites = resp(:,eid==1);
% v4sites_base = resp_base(:,eid==1);
% badTrials = cellfun(@isempty,{params.set});
% v4sites(badTrials,:) = [];
% v4sites_base(badTrials,:) = [];
% params(badTrials) = [];
% sets = unique([params.set]);
% 
% [params,idx] = sortStruct(params,'num');
% v4sites = v4sites(idx,:);
% 
% figure; hold on
% for ii=1:length(sets)
%     idx = [params.set]==sets(ii);
%     [dec_num,rs] = decodeFn([params(idx).num]',v4sites(idx,:));
%     plot([params(idx).num],dec_num,'.');
% end
% fixPlot(gca,[0 max([params.num])+1],[0 max([params.num])+1],'num','dec',1:max([params.num]),1:max([params.num]))

%% 3d score
% figure('color','w','pos',[440,231,609,566])
% [id,~,grp] = unique([[params.set]' [params.num]'],'rows');
% mResp = groupsummary(v4sites,grp,'mean');
% r3d = mean(mResp(id(:,1)==1,:),1);
% r2d = mean(mResp(id(:,1)==2,:),1);
% r2d_rand = mean(mResp(id(:,1)==3,:),1);
% scv4 = (r3d-r2d)./max([r3d;r2d]);
% scv4_vrand = (r3d-r2d_rand)./max([r3d;r2d_rand]);
% subplot(221); h = cdfplot(scv4); h.LineWidth = 2; hold on;
% h = cdfplot(scv4_vrand); h.LineWidth = 2; hold on;
% line([0 0],[0 1],'color','k','linestyle','--','linewidth',2)
% fixPlot(gca,[-0.7 0.7],[0 1],'3d score','prob',-0.5:0.25:0.5,0:0.25:1,'3d score')
% 
% subplot(222); 
% histogram(scv4,linspace(-1,1,35),'DisplayStyle','stairs','linewidth',2); hold on;
% histogram(scv4_vrand,linspace(-1,1,35),'DisplayStyle','stairs','linewidth',2);
% % histogram(scv1,linspace(-1,1,35),'DisplayStyle','stairs','linewidth',2);
% % histogram(scv1_vrand,linspace(-1,1,35),'DisplayStyle','stairs','linewidth',2);
% line([0 0],[0 20],'color','k','linestyle','--','linewidth',2)
% fixPlot(gca,[-0.7 0.7],[0 20],'3d score','count',-0.5:0.25:0.5,0:10:100,'3d score',{'match','rand'})
% legend('Location','northwest')
% 
% % score and resp scatter between 2d match and rand
% mResp = groupsummary(v4sites,grp,'mean');
% r3d = mean(mResp(id(:,1)==1,:),1);
% r2d = mean(mResp(id(:,1)==2,:),1);
% r2d_rand = mean(mResp(id(:,1)==3,:),1);
% 
% subplot(223);
% plot(r3d,r2d,'.','markersize',15); hold on
% plot(r3d,r2d_rand,'.','markersize',15)
% line([-20 250],[-20 250],'color','k','linewidth',2,'linestyle','--')
% fixPlot(gca,[-10 250],[-10 250],'3d resp','2d resp',0:50:150,0:50:150,'v4 mua responses',{'match','rand'})
% 
% 
% subplot(224);
% plot(scv4,scv4_vrand,'k.','markersize',15); hold on
% fixPlot(gca,[-0.6 0.7],[-0.6 0.7],'score with match','score with rand',-0.5:0.25:0.5,-0.5:0.25:0.5,'v4 scores')
% line([-20 150],[-20 150],'color','k','linewidth',2,'linestyle','--')

%% helpers
function [yh,corr_y] = decodeFn(y,X)
    yh = nan(size(y));
    for ii=1:length(y)
        xval = X(ii,:);
        ytrain = y; ytrain(ii) = [];
        xtrain = X; xtrain(ii,:) = [];
        
        yh(ii) = xval*regress(ytrain,xtrain);
    end
    corr_y = corr(yh,y);
end

function rfv4 = rfmap_save(params,resp,eid)    
    badTrials = cellfun(@isempty,{params.set});
    resp(badTrials,:) = [];
    params(badTrials) = [];

    resp(~[params.good],:) = [];
    params(~[params.good]) = [];
    v4sites = resp(:,eid == 1);

    % STA v4
    [pos,~,grp] = unique([[params.x]' [params.y]' [params.s]'],'rows');
    mResp = groupsummary(v4sites,grp,'mean');
    rfv4 = getRF_sta(mResp,pos,zeros(1,64));
    xx = median(rfv4.staFit_beta(:,2));
    yy = median(rfv4.staFit_beta(:,4));
    rr(1) = median(abs(rfv4.staFit_beta(:,3)));
    rr(2) = median(abs(rfv4.staFit_beta(:,5)));
    rr = min(rr)/2;
    rfv4.pos = [xx yy rr];
    rfv4.posDeg = rad2deg(atan(([xx yy rr]*(609.6/sqrt(1920^2 + 1080^2)))/540));
end

function rfmap_plot(rfv4)
    pos = [280 -130 280];
    % pos = [-274 -198 340];
    stimPos = [pos(1)-pos(3)/2 pos(2)-pos(3)/2 pos(3) pos(3)];

    figure('color','w','pos',[147,86,754,711]); 
    ha = tight_subplot(8,8,0.02,0.02,0.02);
    for ii=1:64
        imagesc(rfv4.x(:),rfv4.y(:),rfv4.sta(:,:,ii),'parent',ha(ii));
        hold(ha(ii),'on');
        axis(ha(ii),'image','off');
        set(ha(ii),'ydir','normal');
        plot(0,0,'w.','parent',ha(ii),'markersize',10);
        rectangle('Position',stimPos,'LineWidth',2,'EdgeColor','w','parent',ha(ii))
%         drawellipse(ha(ii),'center',rfv4.staFit_beta(ii,[2 4]),'SemiAxes',abs(rfv4.staFit_beta(ii,[3 5])/2),'Color','k','facealpha',0,'InteractionsAllowed','none','selected',false,'linewidth',1);
    end
    
    figure('color','w','pos',[909,350,254,445]);
    ha = tight_subplot(2,1); axes(ha(1))
    imagesc(rfv4.x(:),rfv4.y(:),mean(rfv4.sta,3));
    axis image; set(gca,'ydir','normal'); hold on;
    plot(0,0,'k.','parent',gca,'markersize',10);
    rectangle('Position',stimPos,'LineWidth',2,'EdgeColor','w')
    axis([-200 960 -540 200])
    title(['V4: pix ' num2str(round(rfv4.pos))])
    
    axes(ha(2)); 
 
    % clf; hold on;
    for ii=1:64
        plot(gca,rfv4.staFit_beta(ii,2),rfv4.staFit_beta(ii,4),'.','markersize',15,'color',[0.2 0.9 0.3]*0.8);
        drawellipse(gca,'center',rfv4.staFit_beta(ii,[2 4]),'SemiAxes',abs(rfv4.staFit_beta(ii,[3 5])/2),'Color',[0.2 0.9 0.3],'facealpha',0,'InteractionsAllowed','none','selected',false,'linewidth',1,'edgealpha',0.5);
    end
    plot(0,0,'r.','markersize',20);
    rectangle('Position',stimPos,'LineWidth',2,'EdgeColor','k')
    axis([-200 960 -540 200])
   
    rectangle('Position',[rfv4.pos(1)-rfv4.pos(3) rfv4.pos(2)-rfv4.pos(3) rfv4.pos(3)*2 rfv4.pos(3)*2],'LineWidth',2,'EdgeColor','b')
    
    % get 5 degree ticks
    ticks = 5*(200./rad2deg(atan((200*(609.6/sqrt(1920^2 + 1080^2)))/540)));
    fixPlot(gca,[-200 800],[-500 200],'','',-4*ticks:ticks:ticks*4,-4*ticks:ticks:ticks*4)
    axis normal
    axis equal
end


function mapRecord = crawlForFiles
    filelist = cell(1);
    rootDir = '/Volumes/colada/Ram/data/neural/v4-7a';
    dirs = dir([rootDir '/2510*']);
    for ii=1:length(dirs)
        files_grid = dir([rootDir '/' dirs(ii).name '/*map*nev*']);
        filelist = [filelist cellfun(@(x) strrep([rootDir '/' dirs(ii).name '/' x],'.nev',''),{files_grid.name},'UniformOutput',false)];
        % if length(files_nev) ~= length(files_grid)
        %     for jj=1:length(files_nev)
        %         files_nev.name(jj)
        %     end
        % end
    end
    filelist(1) = [];
    mapRecord = cell2struct(filelist,'path',1);
end

function exptRecord = crawlForFiles_local
    files_grid = dir('*grid*dense*');
    files_grid = cellfun(@(x) [pwd '/' strrep(x,'_dense.mat','')],{files_grid.name},'UniformOutput',false);
    exptRecord = cell2struct(files_grid,'path',1);
end