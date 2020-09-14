% [quantif, data] = BestResp_BeadsSelec_v5_draft(data, 0, [0.2 0.2], 1, [0.01 0.05], 1, 1, [])

function [quantif, data] = BestResp_BeadsSelec_v5(data, selection, RespLag, comp_mthd, alpha, IndivROI, saveFig, SessionData)
% % This code runs through the visual stim until it finds one significant
% % vis stim are order based on their median response

% selection: 0, no selection; 1, bead-based
% RespLag: time window (s) to use to determine the best vis stim (respective to stim onset and offset, respectively) 
% comp_mthd: 1, difference bbtwn dFF; 2, tbd, for future usage?
% alpha = 2-element vector for 1) significant visual response and 2) effect of  opto

global mainDir
global saveName

%% Preambule
fprintf('\n')
if saveFig
    mainDir =  data.Info.mainDir;
    if exist(mainDir,'dir')
        disp(mainDir)
    else
        mainDir = uigetdir('X:\camille.mazo\2P_processed','Remap main dir');
        fprintf('Main directory has changed to %s \n',mainDir)
        data.Info.mainDir = mainDir;
    end
    saveDir = [data.Info.mainDir, filesep 'Figures\PopPlots'];
    saveIndivROIsDir = [data.Info.mainDir, filesep 'Figures\IndivROIs'];
    if ~ exist(saveDir)
        mkdir(saveDir)
    end
    if ~ exist(saveIndivROIsDir)
        mkdir(saveIndivROIsDir)
    end
        saveName = data.Info.saveName;
end


if selection == 1
    if ~isfield(data,'beadID') % Only for verisions of "IDfying_redCells.m" anterior to June 15, 2020
        for i = 1:size(data.beads_pos,1)
            ROInum = data.beads_pos(i,1);
            data.beads_pos(i,3) = sum(data.s2p.iscell(1:ROInum));
        end
        selec = zeros(1,size(data.data_sorted,1));
        selec(data.beads_pos(:,3)) = 1;
    else
        selec =  data.beadID.bead_pos;
    end
elseif selection == 0
    selec = ones(1,size(data.data_sorted,1));
end
% selec = logical(selec);

timeStimInd = data.Info.timeStimInd;
timeVect = data.Info.timeVect;
xvect = data.Info.xvect;
xstim_opto = data.Info.xstim_opto;
tAna = timeVect>data.Info.xstim_vis(1)+RespLag(1) & timeVect<data.Info.xstim_vis(2)+RespLag(2);
tBase = timeVect(1) & timeVect<data.Info.xstim_vis(1);
nROIs = size(data.data_sorted,1); 
nBeadsPos = sum(selec);
nBeadsNeg = sum(~selec);

color_bp = [0.8500 0.3250 0.0980];
color_bn = [0      0.4470 0.7410];

% oldfolder = cd([mainDir, filesep 'analysis']);
% dt = datestr(now,'yyyymmdd');
% diary off; diary(['BestResp_BeadsSelec_' dt]);
% cd(oldfolder)

disp(saveName)
disp(['Baseline, frames: ', num2str(find(tBase,1)), ' - ', num2str(find(tBase,1,'last'))])
disp(['Response Window, frames: ', num2str(find(tAna,1)), ' - ', num2str(find(tAna,1,'last'))])

%% 1. Best significant visual response

%% 1.1 ------------ find best significant vis resp ------------------------
meanResp = squeeze(nanmean(median(data.data_sorted(:,tAna,2:9,:,:,:),6,'omitnan'),2));
if ndims(meanResp) == 4 
   meanResp_vis = meanResp(:,:,:,1);
elseif  ndims(meanResp) == 3 % if only one SF/TF couple
   meanResp_vis = meanResp(:,:,1);
end

% h_vis = NaN(size(data.data_sorted,1),1);
% if ndims(meanResp) == 4
%     h_vis2 = zeros(size(data.data_sorted,1),2);
% else
%     h_vis2 = zeros(size(data.data_sorted,1),1);
% end

n_run = NaN(size(data.data_sorted,1),1);
h_vis = false(1,size(data.data_sorted,1));
p_vis = false(1,size(data.data_sorted,1));
% selVisStim = NaN(size(data.data_sorted,1),size(data.data_sorted,2),2,size(data.data_sorted,6));
for i = 1:size(data.data_sorted,1)
    if ndims(meanResp) == 4 % if 2 SF/TF couple
        test2 = cat(2,meanResp_vis(i,:,1),meanResp_vis(i,:,2));
    else
        test2 = meanResp_vis(i,:);
    end
    [~,I] = sort(test2,'descend');
    h1 = 0; ii = 1;
    while h1 == 0 && ii < 4 %size(I,2)+1
        if I(ii) <= size(meanResp_vis,2)
%             disp(['ROI', num2str(i),'; SF = 1'])
            vis_resp = squeeze(mean(data.data_sorted(i,tAna,I(ii)+1,1,1,:),2));
            vis_base = squeeze(mean(data.data_sorted(i,tBase,I(ii)+1,1,1,:),2));
            [h1, p1] = ttest(vis_resp, vis_base,'Alpha',alpha(1));
            sf = 1;
        else
%             disp(['ROI', num2str(i),'SF 2'])
            vis_resp = squeeze(mean(data.data_sorted(i,tAna,I(ii)+1-size(meanResp_vis,2),2,1,:),2));
            vis_base = squeeze(mean(data.data_sorted(i,tBase,I(ii)+1-size(meanResp_vis,2),2,1,:),2));
            [h1, p1] = ttest(vis_resp, vis_base,'Alpha',alpha(1));
            sf = 2;
        end
        nrun = ii;
        ii = ii+1;
    end
    n_run(i) = nrun;
    h_vis(i) = h1; p_vis(i) = p1;
%     h_vis2(i,sf) = h1;
    if h1
%         selVisStim = true
        selVisStim(i,:,:,:) = squeeze(data.data_sorted(i,:,I(ii-1)-(sf-1)*size(meanResp_vis,2)+1,sf,:,:));
    else
%                 selVisStim = false    
        selVisStim(i,:,:,:) = NaN(1,size(data.data_sorted,2),2,size(data.data_sorted,6));
    end
end
% disp(['Fraction of cell passing: ', num2str(sum(h_vis)/size(h_vis,2),2)])
figure;histogram(n_run(h_vis));
xl = xlim;yl = ylim;
ylabel('# of cells'); xlabel('# of run');
xticks([1 2 3]);
title('Number of run to pass')
text(xl(1)+0.1,yl(2),['Fraction of cell passing: ', num2str(sum(h_vis)/size(h_vis,2),2)],...
    'HorizontalAlignment','left','VerticalAlignment','top');
if saveFig
screen2png([saveDir filesep saveName,'_ImpactOfLaserBestVis_Scatter']);
end

nVisResp(1) = sum(h_vis(selec));
nVisResp(2) = sum(h_vis(~selec));

%% 1.2 ------------ scatter plot ------------------------------------------
h_vis = logical(h_vis);
figure; hold on 
meanResp_selec = squeeze(nanmean(median(selVisStim(:,tAna,:,:),4,'omitnan'),2));
[~,p_all] = ttest(meanResp_selec(h_vis,1),meanResp_selec(h_vis,2));
[~,p_bp] = ttest(meanResp_selec(selec,1),meanResp_selec(selec,2));
[~,p_bn] = ttest(meanResp_selec(~selec,1),meanResp_selec(~selec,2));
s1 = scatter(meanResp_selec(~selec,1),meanResp_selec(~selec,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = 0.2;
scatter(meanResp_selec(h_vis&selec,1),meanResp_selec(h_vis&selec,2),...
   'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
xl = xlim; yl = ylim;
axismax = max(xl(2),yl(2));
axismin = min(xl(1),yl(1));
xlim([axismin axismax]); ylim([axismin axismax])
plot(xlim,ylim,'k:');
xl = xlim; yl = ylim;
title([{'Effect Of Laser on Best Significant Vis Resp'},{['Significant Resp: ttest w/ alpha = ',num2str(alpha(1))]}])
text(xl(1)-0.1*xl(1),yl(2)-0.05*yl(2),['red: bead pos, n: ', num2str(nVisResp(1)),'/', num2str(nBeadsPos),' = ',num2str(nVisResp(1)/nBeadsPos,2),'; ttest, p = ',num2str(p_bp,2)],...
    'HorizontalAlignment','Left');
text(xl(1)-0.1*xl(1),yl(2)-0.1*yl(2),['blue: bead neg, n: ', num2str(nVisResp(2)),'/', num2str(nBeadsNeg),' = ',num2str(nVisResp(2)/nBeadsNeg,2),'; ttest, p = ',num2str(p_bn,2)],...
    'HorizontalAlignment','Left');
text(xl(1)-0.1*xl(1),yl(2)-0.15*yl(2),['(population t-test, p = ',num2str(p_all,2), ', n = ', num2str(sum(h_vis)),')'],...
    'HorizontalAlignment','Left');
xlabel({'dF/F','Best Vis Stim'}); ylabel({'Same Vis Stim+Photostim','dF/F'})
if saveFig
    saveas(gcf,[saveDir filesep saveName,'_ImpactOfLaserBestVis_Scatter'],'fig');
    screen2png([saveDir filesep saveName,'_ImpactOfLaserBestVis_Scatter']);
end

%% 1.3 ------------ bar graph ---------------------------------------------
%% 1.3.1 Response categorization
AllTrialsResp_selec = squeeze(nanmean(selVisStim(:,tAna,:,:),2));
h_opto = NaN(1,size(data.data_sorted,1)); p_opto = NaN(1,size(data.data_sorted,1));
pos = zeros(1,size(data.data_sorted,1)); neg = zeros(1,size(data.data_sorted,1)); nochange = zeros(1,size(data.data_sorted,1));
for i = 1 : size(selVisStim,1)
    if h_vis(i)
        [h2, p2] = ttest(AllTrialsResp_selec(i,1,:),AllTrialsResp_selec(i,2,:),'Alpha',alpha(2));
        h_opto(i) = h2; p_opto(i) = p2;
        if h2 && meanResp_selec(i,2) > meanResp_selec(i,1)
            pos(i)= true;
        elseif  h2 && meanResp_selec(i,2) < meanResp_selec(i,1)
            neg(i) = true;
        elseif h2 == 0
            nochange(i) = true;
        else
            debug(i) = true;
            disp('bug. this cannot be')
        end
    end
end
final_beads_vis = [sum(pos(selec),'omitnan') sum(neg(selec),'omitnan') sum(nochange(selec),'omitnan');...
    sum(pos(~selec),'omitnan') sum(neg(~selec),'omitnan') sum(nochange(~selec),'omitnan')];

x1 = [repmat('a',nVisResp(1),1); repmat('b',nVisResp(2),1)];
x2 = [repmat(1,final_beads_vis(1,1),1) ; repmat(2,final_beads_vis(1,2),1); repmat(3,final_beads_vis(1,3),1);...
    repmat(1,final_beads_vis(2,1),1) ; repmat(2,final_beads_vis(2,2),1); repmat(3,final_beads_vis(2,3),1)];
[~,chi2stat,pval] = crosstab(x1,x2);
nVisResp = repmat(nVisResp',1,3);
final_beads_vis = final_beads_vis./ nVisResp;

figure;
subplot(1,3,1); hold on
b1 = bar(final_beads_vis,'stacked');
legend('excited','inhibited','no change')
yVal = [0 0];
for i = 1:3
    yVal = yVal+double(b1(i).YData);
text(double(b1(i).XData), yVal,...
   {num2str(final_beads_vis(1,i),2), num2str(final_beads_vis(2,i),2)},...
    'HorizontalAlignment','center','VerticalAlignment','top')
end
xticks([1 2]); xticklabels({'Beads+','Beads-'})
text(1, 1, ['n = ',num2str(nVisResp(1,1))],...
    'HorizontalAlignment','Center','VerticalAlignment','bottom')
text(2, 1, ['n = ',num2str(nVisResp(2,2))],...
    'HorizontalAlignment','Center','VerticalAlignment','bottom')
title({'Response Categorization',['chi2 stat = ',num2str(chi2stat,2),'; p = ',num2str(pval,2)]})
ylabel('Fraction')
ylim([0 1.05]);

%% 1.3.2 Response Magnitude, average
subplot(1,3,2); hold on
ImpactOfLaser_Vis = meanResp_selec(:,2) - meanResp_selec(:,1);
[~,p] = ttest2(ImpactOfLaser_Vis(selec), ImpactOfLaser_Vis(~selec));

delta_dFF(1,1) = nanmean(ImpactOfLaser_Vis(selec));
delta_dFF(1,2) = nanmean(ImpactOfLaser_Vis(pos&selec));
delta_dFF(1,3) = nanmean(ImpactOfLaser_Vis(~pos&~nochange&selec));
delta_dFF(1,4) = nanmean(ImpactOfLaser_Vis(nochange&selec));
delta_dFF(2,1) = nanmean(ImpactOfLaser_Vis(~selec));
delta_dFF(2,2) = nanmean(ImpactOfLaser_Vis(pos&~selec));
delta_dFF(2,3) = nanmean(ImpactOfLaser_Vis(~pos&~nochange&~selec));
delta_dFF(2,4) = nanmean(ImpactOfLaser_Vis(nochange&~selec));

n(1,1) = sum(~isnan(ImpactOfLaser_Vis(selec)));
n(1,2) = sum(~isnan(ImpactOfLaser_Vis(pos&selec)));
n(1,3) = sum(~isnan(ImpactOfLaser_Vis(~pos&~nochange&selec)));
n(1,4) = sum(~isnan(ImpactOfLaser_Vis(nochange&selec)));
n(2,1) = sum(~isnan(ImpactOfLaser_Vis(~selec)));
n(2,2) = sum(~isnan(ImpactOfLaser_Vis(pos&~selec)));
n(2,3) = sum(~isnan(ImpactOfLaser_Vis(~pos&~nochange&~selec)));
n(2,4) = sum(~isnan(ImpactOfLaser_Vis(nochange&~selec)));

delta_dFF_sem(1,1) = std(ImpactOfLaser_Vis(selec),'omitnan')./sqrt(n(1,1));
delta_dFF_sem(1,2) = std(ImpactOfLaser_Vis(pos&selec),'omitnan')./sqrt(n(1,2));
delta_dFF_sem(1,3) = std(ImpactOfLaser_Vis(~pos&~nochange&selec),'omitnan')./sqrt(n(1,3));
delta_dFF_sem(1,4) = std(ImpactOfLaser_Vis(nochange&selec),'omitnan')./sqrt(n(1,4));
delta_dFF_sem(2,1) = std(ImpactOfLaser_Vis(~selec),'omitnan')./sqrt(n(2,1));
delta_dFF_sem(2,2) = std(ImpactOfLaser_Vis(pos&~selec),'omitnan')./sqrt(n(2,2));
delta_dFF_sem(2,3) = std(ImpactOfLaser_Vis(~pos&~nochange&~selec),'omitnan')./sqrt(n(2,3));
delta_dFF_sem(2,4) = std(ImpactOfLaser_Vis(nochange&~selec),'omitnan')./sqrt(n(2,4));

% b1 = bar([nanmean(ImpactOfLaser_Vis(selec));nanmean(ImpactOfLaser_Vis(~selec))]);
% b1.FaceColor = 'flat';
% b1.CData(1,:) = color_bp;
% b1.CData(2,:) = color_bn;

b = bar(delta_dFF');
xticks([1 2 3 4]); xticklabels({'All','Pos','Neg', 'No Change'})
drawnow
for i = 1:size(delta_dFF,2)
    plot([i+b(1).XOffset i+b(1).XOffset],[delta_dFF(1,i)-delta_dFF_sem(1,i) delta_dFF(1,i)+delta_dFF_sem(1,i)],'k-');
    plot([i-b(1).XOffset i-b(1).XOffset],[delta_dFF(2,i)-delta_dFF_sem(2,i) delta_dFF(2,i)+delta_dFF_sem(2,i)],'k-');
end
b(1).FaceColor = color_bp;
b(2).FaceColor = color_bn;
ylabel('delta dFF')
legend(b,{'Beads+','Beads-'},'Location','southwest')
% yl= min(ylim);
% yl = [yl yl yl];
for i = 1:length(b)
    xtips1 = double(b(i).XData);
    ytips1 = double(b(i).YData);
    labels1 = string([n(i,1) n(i,2) n(i,3)  n(i,4)]);
    text(xtips1+b(i).XOffset,ytips1,labels1,...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'Color',[0 0 0])
end
title({'Response Magnitude, Averaged',['ttest, p = ',num2str(p,2)]})

%% 1.3.3 Response Magnitude, all cells
subplot(1,3,3); hold on
x = ones(1,size(ImpactOfLaser_Vis(selec&h_vis),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec&h_vis),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 2*ones(1,size(ImpactOfLaser_Vis(selec&h_vis&pos),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec&h_vis&pos),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 3*ones(1,size(ImpactOfLaser_Vis(selec&h_vis&neg),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec&h_vis&neg),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 4*ones(1,size(ImpactOfLaser_Vis(selec&h_vis&nochange),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec&h_vis&nochange),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);

x = ones(1,size(ImpactOfLaser_Vis(~selec&h_vis),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec&h_vis),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
x = 2*ones(1,size(ImpactOfLaser_Vis(~selec&h_vis&pos),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec&h_vis&pos),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
x = 3*ones(1,size(ImpactOfLaser_Vis(~selec&h_vis&neg),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec&h_vis&neg),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
x = 4*ones(1,size(ImpactOfLaser_Vis(~selec&h_vis&nochange),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec&h_vis&nochange),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);

% % plot(1+b(1).XOffset,ImpactOfLaser_Vis(selec&h_vis),'r.')
% try
%     plot(2+b(1).XOffset,ImpactOfLaser_Vis(selec&h_vis&pos),'r.')
% end
% try
%     plot(3+b(1).XOffset,ImpactOfLaser_Vis(selec&h_vis&neg),'r.')
% end
% try
%     plot(4+b(1).XOffset,ImpactOfLaser_Vis(selec&h_vis&nochange),'r.')
% end
% plot(1-b(1).XOffset,ImpactOfLaser_Vis(~selec&h_vis),'b.');
% plot(2-b(1).XOffset,ImpactOfLaser_Vis(~selec&h_vis&pos),'b.')
% plot(3-b(1).XOffset,ImpactOfLaser_Vis(~selec&h_vis&neg),'b.')
% plot(4-b(1).XOffset,ImpactOfLaser_Vis(~selec&h_vis&nochange),'b.')

for i = 1:size(delta_dFF,2)
    plot([i+b(1).XOffset-0.1 i+b(1).XOffset+0.1], [delta_dFF(1,i) delta_dFF(1,i)],'k','LineWidth',1);
    plot([i+b(2).XOffset-0.1 i+b(2).XOffset+0.1], [delta_dFF(2,i) delta_dFF(2,i)],'k','LineWidth',1);
end
xticks([1 2 3 4]);xticklabels({'All','Pos','Neg','No change'})
xlim([.5 4.5])
plot(xlim,[0  0],'k:');
title('Response Magnitude, All Data')

sgtitle([{'Impact of Laser on Best Vis Resp'},{'Significant Cells Only'}])
set(gcf,'Units','Normalized','Position',[0.01 0.45 0.95 0.45])
if saveFig
    saveas(gcf,[saveDir filesep saveName,'_ImpactOfLaserBestVis_Diff_Bar'],'fig');
    screen2png([saveDir filesep saveName,'_ImpactOfLaserBestVis_Diff_Bar']);
end

%% 1.4 ------------ save --------------------------------------------------

quantif.Info = data.Info;
quantif.Info.RespLag = RespLag;
quantif.Info.alpha = alpha;

if selection == 1
quantif.bead_selection = selec;
end
quantif.bestVisResp.h_vis = h_vis;
quantif.bestVisResp.nVisResp = nVisResp;
quantif.bestVisResp.selVisStim = selVisStim;
quantif.bestVisResp.meanResp_selec = meanResp_selec;
quantif.bestVisResp.h_opto = h_opto;
quantif.bestVisResp.pos = pos;
quantif.bestVisResp.neg = neg;
quantif.bestVisResp.nochange = nochange;
quantif.bestVisResp.ImpactOfLaser = ImpactOfLaser_Vis;
quantif.bestVisResp.deltadFF = delta_dFF;

data.bestVisResp = quantif.bestVisResp;
data.bestVisResp.Info.RespLag = RespLag;
data.bestVisResp.Info.alpha = alpha;

%% 2. Spontaneous activity

%% 2.1 ------------ scatter plot ------------------------------------------
meanResp_spont = squeeze(nanmean(nanmean(data.data_sorted(:,tAna,1,1,:,:),6),2));
figure; hold on
s1 = scatter(meanResp_spont(~selec,1),meanResp_spont(~selec,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = 0.2;
scatter(meanResp_spont(h_vis&selec,1),meanResp_spont(h_vis&selec,2),...
   'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');

[~,p_all] = ttest(meanResp_spont(:,1),meanResp_spont(:,2));
[~,p_bp] = ttest(meanResp_spont(selec,1),meanResp_spont(selec,2));
[~,p_bn] = ttest(meanResp_spont(~selec,1),meanResp_spont(~selec,2));
xl = xlim; yl = ylim;
axismax = max(xl(2),yl(2));
axismin = min(xl(1),yl(1));
xlim([axismin axismax]); ylim([axismin axismax])
plot(xlim,ylim,'k:');
clear xlim ylim
xlim(xl); ylim(yl)
xlabel({'dF/F','Spont Activity'}); ylabel({'Photostim','dF/F'})
text(xl(1)-0.1*xl(1),yl(2)-0.1*yl(2),['bead pos; t-test, p = ',num2str(p_bp,2), ', n = ', num2str(sum(selec))]); 
text(xl(1)-0.1*xl(1),yl(2)-0.2*yl(2),['bead neg; t-test, p = ',num2str(p_bn,2), ', n = ', num2str(sum(~selec))]); 
text(xl(1)-0.1*xl(1),yl(2)-0.3*yl(2),['population t-test, p = ',num2str(p_all,2), ', n = ', num2str(size(data.data_sorted,1))]); 
title('Effect Of Laser on Spont Activity, All Cells')
if saveFig
    saveas(gcf,[saveDir filesep saveName,'_ImpactOfLaserSpont_Scatter'],'fig');
    screen2png([saveDir filesep saveName,'_ImpactOfLaserSpont_Scatter']);
end

%% 2.2 ------------ bar graph ---------------------------------------------
ImpactOfLaser_Spont = meanResp_spont(:,2) - meanResp_spont(:,1);
delta_dFF_spont(1) = mean(ImpactOfLaser_Spont(selec));
delta_dFF_spont(2) = mean(ImpactOfLaser_Spont(~selec));
n(1,1) = size(ImpactOfLaser_Spont(selec),1);
n(2,1) = size(ImpactOfLaser_Spont(~selec),1);
delta_dFF_spont_sem(1) = std(ImpactOfLaser_Spont(selec))./sqrt(n(1,1));
delta_dFF_spont_sem(2) = std(ImpactOfLaser_Spont(~selec))./sqrt(n(2,1));
[~,p] = ttest2(ImpactOfLaser_Spont(selec),ImpactOfLaser_Spont(~selec));

figure;
subplot(1,3,2); hold on
b2 = bar(delta_dFF_spont);
for i = 1:2
    plot([i i],[delta_dFF_spont(i)-delta_dFF_spont_sem(i) delta_dFF_spont(i)+delta_dFF_spont_sem(i)],'k-');
end
xticks([1 2]); xticklabels({'Beads+','Beads-'})
b2.FaceColor = 'flat';
b2.CData(1,:) = color_bp;
b2.CData(2,:) = color_bn;
text([1 2],double(b2.YData),string([n(1,1) n(2,1)]),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Color',[0 0 0])
ylabel('delta dFF')
title({'Magnitude of the effect, averaged',['ttest, p = ',num2str(p,2)]})

subplot(1,3,3); hold on
x = ones(1,size(ImpactOfLaser_Spont(selec),1));
scatter(x,ImpactOfLaser_Spont(selec),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 2*ones(1,size(ImpactOfLaser_Spont(~selec),1));
scatter(x,ImpactOfLaser_Spont(~selec),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
% for i = 1:size(ImpactOfLaser_Spont,1)
% if any(i == data.beads_pos(:,3))
%     plot(1,ImpactOfLaser_Spont(i),'Color',color_bp,'Marker','.')
% else
%     plot(2,ImpactOfLaser_Spont(i),'Color',color_bn,'Marker','.')
% end
% end
plot([0.9 1.1], [mean(ImpactOfLaser_Spont(selec)) mean(ImpactOfLaser_Spont(selec))],'k', 'LineWidth',1)
plot([1.9 2.1], [mean(ImpactOfLaser_Spont(~selec)) mean(ImpactOfLaser_Spont(~selec))],'k', 'LineWidth',1)
xticks([1 2]); xticklabels({'Beads+','Beads-'})
xlim([0.5 2.5])
title('Magnitude of the effect, all data') 

Resp_Laser = squeeze(nanmean(data.data_sorted(:,tAna,1,1,2,:),2));
Base_Laser = squeeze(nanmean(data.data_sorted(:,tBase,1,1,2,:),2));
SpontR_Laser = squeeze(nanmean(data.data_sorted(:,tAna,1,1,1,:),2));
SpontB_Laser = squeeze(nanmean(data.data_sorted(:,tBase,1,1,1,:),2));

exc_spont = false(2,size(data.data_sorted,1));
inhib_spont = false(2,size(data.data_sorted,1));
nothing_spont = false(2,size(data.data_sorted,1));
for i = 1:size(data.data_sorted,1)
    [h1,s1] = ttest(Resp_Laser(i,:),Base_Laser(i,:),'Alpha',alpha(2));
    h_all(1,i) = h1;     p_all(1,i) = s1;
    [h2,p2] = ttest(SpontR_Laser(i,:),SpontB_Laser(i,:),'Alpha',alpha(2));
    h_all(2,i) = h2;     p_all(2,i) = p2;
    if h1==1 && mean(Resp_Laser(i,:)) > mean(Base_Laser(i,:))
        exc_spont(1,i) = true;
    elseif h1==1 && mean(Resp_Laser(i,:)) < mean(Base_Laser(i,:))
        inhib_spont(1,i) = true;
    elseif h1 == 0
        nothing_spont(1,i) = true;
    end
        if h2==1 && mean(SpontR_Laser(i,:)) > mean(SpontB_Laser(i,:))
        exc_spont(2,i) = true;
    elseif h2==1 && mean(SpontR_Laser(i,:)) < mean(SpontB_Laser(i,:))
        inhib_spont(2,i) = true;
    elseif h2 == 0
        nothing_spont(2,i) = true;
    end
end
final_vector = [sum(exc_spont(1,:))/nROIs sum(inhib_spont(1,:))/nROIs sum(nothing_spont(1,:))/nROIs;...
    sum(exc_spont(2,:))/nROIs sum(inhib_spont(2,:))/nROIs sum(nothing_spont(2,:))/nROIs];

clear x1 x2 chi2stat pval;
final_beads = [sum(exc_spont(1,selec)) sum(inhib_spont(1,selec)) sum(nothing_spont(1,selec));...
    sum(exc_spont(1,~selec)) sum(inhib_spont(1,~selec)) sum(nothing_spont(1,~selec))];

x1 = [repmat('a',n(1,1),1); repmat('b',n(2,1),1)];
x2 = [repmat(1,final_beads(1,1),1) ; repmat(2,final_beads(1,2),1); repmat(3,final_beads(1,3),1);...
    repmat(1,final_beads(2,1),1) ; repmat(2,final_beads(2,2),1); repmat(3,final_beads(2,3),1)];
[~,chi2stat,pval] = crosstab(x1,x2);

nBeadsPos = repmat(nBeadsPos,1,3);
nBeadsNeg = repmat(nBeadsNeg,1,3);
final_beads_ratio = final_beads./[nBeadsPos;nBeadsNeg];

subplot(1,3,1)
b2 = bar(final_beads_ratio,'stacked');
legend('excited','inhibited','no change')
xticklabels({'Beads+','Beads-'})
ylabel('Fraction')
yVal  = [0 0];
for i = 1:3
    yVal = yVal+double(b2(i).YData);
    text(double(b2(i).XData), yVal,...
   {num2str(b2(i).YData(1),2), num2str(b2(i).YData(2),2)},...
    'HorizontalAlignment','center','VerticalAlignment','top')
end
text(1,1,['n = ',num2str(nBeadsPos(1))],...
    'HorizontalAlignment','Center','VerticalAlignment','bottom');
text(2,1,['n = ',num2str(nBeadsNeg(2))],...
    'HorizontalAlignment','Center','VerticalAlignment','bottom');
ylim([0  1.05]);
title({'Response Categorization',['chi2 stat = ', num2str(chi2stat,2), ' ; p = ', num2str(pval,2)]})

sgtitle('Impact of Laser on Spontaneous Activity')
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.95 0.4])

if saveFig
saveas(gcf,[saveDir filesep saveName,'_ImpactOfLaserSpont_StackedBar'],'fig');
screen2png([saveDir filesep saveName,'_ImpactOfLaserSpont_StackedBar']);
end

figure; hold on
bar(final_vector,'stacked')
legend('excited','inhibited','no change')
xticks([1 2]); xticklabels({'Laser','NoLaser'})
title('Determining False Positive rate')
screen2png([saveDir filesep saveName,'_FalsePosRate_StackedBar']);

%% 2.3 ------------ save---------------------------------------------------

quantif.spont.meanResp = meanResp_spont;
quantif.spont.ImpactOfLaser = ImpactOfLaser_Spont;
quantif.spont.delta_dFF = delta_dFF_spont;
quantif.spont.h_opto = h_all(1,:);
quantif.spont.final_beads = final_beads;

data.spont = quantif.spont;

save([mainDir, filesep 'analysis\quantif_', saveName],'quantif');
save([mainDir, filesep 'analysis\', saveName],'data','-append');
disp('data mean resp saved');

%% Individual ROIs
if IndivROI
    data = PlotIndivROIs(data, [], RespLag, selection, alpha(1),1);
end
clear global mainDir
clear global saveName
% dairy off
end