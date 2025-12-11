clc; clear; close all
addpath('dep')
load('data/patchRecord_beh.mat','patchRecord_beh')
dataLoc = '~/Downloads/v4-7a/patch';

patchRecord_beh(cellfun(@isempty,{patchRecord_beh.name})) = [];
patchRecord_beh([patchRecord_beh.nAnalyze]<400) = [];

for rr=1:length(patchRecord_beh)
    load([dataLoc '/' patchRecord_beh(rr).name '_beh.mat'],'exptProps','patchProps','trialProps')

    exptTime = sum([trialProps.ISI]+[trialProps.stim_ontime]);
    for pp=1:exptProps.patch_n
        idx = [trialProps.patch]==pp;
        patchProps(pp).idx = pp;
        patchProps(pp).totalfixtime = (sum([trialProps(idx).ISI] + [trialProps(idx).stim_ontime]))./exptTime;

        patchProps(pp).sigm_beta_min = patchProps(pp).sigm_beta(1);
        patchProps(pp).sigm_beta_max = patchProps(pp).sigm_beta(2);
        patchProps(pp).sigm_beta_mid = patchProps(pp).sigm_beta(3);
        patchProps(pp).sigm_beta_slp = patchProps(pp).sigm_beta(4);

        patchProps(pp).col_r = patchProps(pp).fix_color(1);
        patchProps(pp).col_g = patchProps(pp).fix_color(2);
        patchProps(pp).col_b = patchProps(pp).fix_color(3);
    end
    if rr==1
        allpatchprops = patchProps;
    else
        allpatchprops = [allpatchprops patchProps];
    end
end
propdwelltime = [allpatchprops.totalfixtime]';

allpatchprops = rmfield(allpatchprops,'fix_color');
allpatchprops = rmfield(allpatchprops,'sigm_beta');
allpatchprops = rmfield(allpatchprops,'stim_contrast');
allpatchprops = rmfield(allpatchprops,'totalfixtime');

%%
mdl = fitlm(struct2table(allpatchprops),propdwelltime,'linear','Intercept',false);

figure('Position',[168,73,876,591],'color','w')
rs = plotRx(subplot(231),[allpatchprops.idx]',propdwelltime,1);
fixPlot(gca,[0 5],[0 0.6],'patch number','dwell proportion',1:4,0:0.2:0.6,['r2 = ' num2str(round(rs,3))])
legend off

rs = plotRx(subplot(234),[allpatchprops.stim_set]',propdwelltime,0);
fixPlot(gca,[0 21],[0 0.6],'stimulus set','dwell proportion',[1 5:5:20],0:0.2:0.6,['r2 = ' num2str(round(rs,3))])
legend off


rs = plotRx(subplot(232),[allpatchprops.sigm_beta_min]',propdwelltime,0);
fixPlot(gca,[150 450],[0 0.6],'sigmoid bias','dwell proportion',150:50:400,0:0.2:0.6,['r2 = ' num2str(round(rs,3))])
legend off

rs = plotRx(subplot(233),[allpatchprops.sigm_beta_max]',propdwelltime,0);
fixPlot(gca,[700 1400],[0 0.6],'sigmoid maximum','dwell proportion',700:350:1400,0:0.2:0.6,['r2 = ' num2str(round(rs,3))])
legend off

rs = plotRx(subplot(235),[allpatchprops.sigm_beta_mid]',propdwelltime,1);
fixPlot(gca,[2 11],[0 0.6],'sigmoid x-mid','dwell proportion',2:2:10,0:0.2:0.6,['r2 = ' num2str(round(rs,3))])
legend off

idx = [allpatchprops.sigm_beta_slp]>2;
allpatchprops(idx) = [];
propdwelltime(idx) = [];

rs = plotRx(subplot(236),[allpatchprops.sigm_beta_slp]',propdwelltime,0.5);
fixPlot(gca,[0.1 1.2],[0 0.6],'sigmoid slope','dwell proportion',0.5:0.5:2,0:0.2:0.6,['r2 = ' num2str(round(rs,3))])
legend off

%%
function rs = plotRx(h,x,y,jit)
    axes(h); hold(h,'on');
    mdl = fitlm(x,y);
    plot(h,mdl)
    hc = get(h,'children');
    hc(1).Color = 'k'; hc(1).LineWidth = 2;
    hc(2).Color = 'k'; hc(2).LineWidth = 2;
    
    set(hc(3),'Marker','o','MarkerFaceColor','w','MarkerSize',6,'Color','k')

    if jit > 0
        x = x + jit*2*(rand(length(x),1)-0.5)/5;
        hc(3).XData = x;
    end
    
    rs = mdl.Rsquared.Adjusted;
end