% load('E:\TempData\CMloop19\20200306_pos5bis\analysis\quantif_CMloop19_pos5bis_run3.mat')
% load('E:\TempData\CMloop19\20200306_pos5\analysis\quantif_CMloop19_pos5_run2.mat')

% pop_quantif({q_18_1,q_18_2,q_18_3_90um,q_18_3_150um,q_18_3_1,q_18_3_2,q_19_1,q_19_2,q_19_3,q_19_4,q_19_5_1,q_19_5_2,q_19_6_1,q_19_6_2,q_22_1_1,q_22_1_2,q_22_1_3,q_22_1_4,q_22_2},1)
function pop_quantif(data_list,saveFig)
%% Preanbule

magn_th = 0.15; %0.15
% ############################################

tv = datestr(now, 'yyyy_mm_dd');
fprintf('\n')

if saveFig
    mainDir = uigetdir('E:\TempData\CMloop','Main directory');
    oldfolder = cd([mainDir]);
    diary off; diary(['PopQuantif_' tv]);
    cd(oldfolder)
    
    fprintf('Main directory is %s \n',mainDir)

    saveDir = [mainDir, '\Figures\PopPlots\ImpactOfLaser'];
    if ~ exist(saveDir)
        mkdir(saveDir)
    end
end

% check for duplicates
for i = 1:length(data_list)
   sessionName{i} = data_list{i}.Info.saveName;
end
duplicates = ismember(sessionName, sessionName{1});
duplicates = duplicates(2:end);
if sum(duplicates)
    disp('two times the same session')
else
disp(['number of sessions to analyze: ', num2str(length(data_list))])
end
disp(tv)
for i = 1:length(data_list)
   disp(sessionName{i})
end
color_bp = [0.8500 0.3250 0.0980];
color_bn = [0      0.4470 0.7410];

%  some definiitions
meanResp_vis = [];
h_vis = [];

meanResp_spont = [];
selec = [];
final_beads = zeros(2,3);
ImpactOfLaser_Vis = []; ImpactOfLaser_Spont = [];
pos = []; neg =[]; nochange = [];
for i = 1:length(data_list)
    try
        selec = [selec ; data_list{i}.bead_selection']; % legacy version (before TG L5 data - June 2020)
    catch
        selec = [selec ; data_list{i}.beadID.bead_pos'];
    end

    meanResp_vis = [meanResp_vis ; data_list{i}.bestVisResp.meanResp_selec];
    h_vis = [h_vis data_list{i}.bestVisResp.h_vis];
    pos = [pos data_list{i}.bestVisResp.pos];
    neg = [neg data_list{i}.bestVisResp.neg];
    nochange = [nochange data_list{i}.bestVisResp.nochange];
    ImpactOfLaser_Vis = [ImpactOfLaser_Vis ; data_list{i}.bestVisResp.ImpactOfLaser];
    
    meanResp_spont = [meanResp_spont; data_list{i}.spont.meanResp];
    final_beads = final_beads + data_list{i}.spont.final_beads;
    ImpactOfLaser_Spont = [ImpactOfLaser_Spont ; data_list{i}.spont.ImpactOfLaser];
end
h_vis = logical(h_vis); selec = logical(selec);
respMagn = false(size(meanResp_vis,1),1);
respMagn(meanResp_vis(:,1) > magn_th) = true; % respMagn is de facto selection only h_vis cells

nVisResp(1) = sum(h_vis(selec&respMagn));
nVisResp(2) = sum(h_vis(~selec&respMagn));
nBeadsPos = sum(selec);
nBeadsNeg = sum(~selec);

%% Best Significant Visual Response
try
    alpha = data_list{1}.Info.alpha;
catch
   alpha = data_list{1}.bestVisResp.Info.alpha;
end

mean_bp(1) = mean(meanResp_vis(selec&respMagn,1));
mean_bp(2) = mean(meanResp_vis(selec&respMagn,2));
mean_bn(1) = mean(meanResp_vis(~selec&respMagn,1));
mean_bn(2) = mean(meanResp_vis(~selec&respMagn,2));
sd_bp(1) = std(meanResp_vis(selec&respMagn,1));
sd_bp(2) = std(meanResp_vis(selec&respMagn,2));
sd_bn(1) = std(meanResp_vis(~selec&respMagn,1));
sd_bn(2) = std(meanResp_vis(~selec&respMagn,2));
[~,p_VisResp] = ttest2(meanResp_vis(selec&respMagn,1),meanResp_vis(~selec&respMagn,1));

% ---------------------- scatter plot -------------------------------------
figure;
subplot(1,2,1); hold on
% meanResp_selec = squeeze(nanmean(median(selVisStim(:,tAna,:,:),4,'omitnan'),2));
[~,p_all] = ttest(meanResp_vis(h_vis'&respMagn,1),meanResp_vis(h_vis'&respMagn,2));
[~,p_bp] = ttest(meanResp_vis(selec&respMagn,1),meanResp_vis(selec&respMagn,2));
[~,p_bn] = ttest(meanResp_vis(~selec&respMagn,1),meanResp_vis(~selec&respMagn,2));
s1 = scatter(meanResp_vis(~selec&respMagn,1),meanResp_vis(~selec&respMagn,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = 0.2;
s2 = scatter(meanResp_vis(selec&respMagn,1),meanResp_vis(selec&respMagn,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2.MarkerFaceAlpha = 0.5;
plot([mean_bn(1)-sd_bn(1)/2 mean_bn(1)+sd_bn(1)/2],[mean_bn(2) mean_bn(2)],'b')
plot([mean_bn(1) mean_bn(1)],[mean_bn(2)-sd_bn(2)/2 mean_bn(2)+sd_bn(2)/2],'b')
plot([mean_bp(1)-sd_bp(1)/2 mean_bp(1)+sd_bp(1)/2],[mean_bp(2) mean_bp(2)],'r')
plot([mean_bp(1) mean_bp(1)],[mean_bp(2)-sd_bp(2)/2 mean_bp(2)+sd_bp(2)/2],'r')
xl = xlim; yl = ylim;
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
plot([-axismax axismax],[-axismax  axismax ],'k:');
clear xlim ylim
xlim(xl); ylim(yl)
plot([magn_th magn_th],ylim,'k:')
xlabel({'dF/F','Best Vis Stim'}); ylabel({'Same Vis Stim+Photostim','dF/F'})
text(magn_th, yl(1),num2str(magn_th),...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
title([{'Effect Of Laser on Best Significant Vis Resp'},{['Significant Resp: ttest w/ alpha = ',num2str(alpha(1))]}])
subplot(1,2,2); hold on
xl = xlim; yl = ylim;
text(xl(1),yl(2),...
    [{'bead-'},...
    {['vis resp: n = ' num2str(nVisResp(2)),'/', num2str(nBeadsNeg),' ; ',num2str(nVisResp(2)/nBeadsNeg,2)]},...
    {'mean +- sd:'},...
    {['laser OFF, : ' num2str(mean_bn(1),3) '+- ' num2str(sd_bn(1),3)]},...
    {['laser ON, : ' num2str(mean_bn(2),3) '+- ' num2str(sd_bn(2),3)]},...
    {['ttest, p = ',num2str(p_bn,3)]}],...
    'HorizontalAlignment','Left','VerticalAlignment','top')
 text(xl(2),yl(2)/2,...
    [{'bead+'},...
    {['vis resp: n = ' num2str(nVisResp(1)),'/', num2str(nBeadsPos),' ; ',num2str(nVisResp(1)/nBeadsPos,2)]},...
    {'mean +- sd:'},...
    {['laser OFF, : ' num2str(mean_bp(1),3) '+- ' num2str(sd_bp(1),3)]},...
    {['laser ON, : ' num2str(mean_bp(2),3) '+- ' num2str(sd_bp(2),3)]},...
    {['ttest, p = ',num2str(p_bp,3)]}],...
    'HorizontalAlignment','right','VerticalAlignment','top')   
text(xl(1),yl(1),...
    [{['population t-test, p = ',num2str(p_all,2), ', n = ', num2str(sum(h_vis'&respMagn))]},...
    {['Vis. Resp. Magnitude, bead+ vs bead-; ttest: p = ', num2str(p_VisResp,3)]}],...
    'HorizontalAlignment','left','VerticalAlignment','bottom');
axis off
set(gcf,'Units','Normalized','Position',[0.4 0.4 0.4 0.4])
if saveFig
    saveas(gcf,[saveDir filesep tv '_ImpactOfLaserBestVis_Scatter'],'fig');
    screen2png([saveDir filesep tv '_ImpactOfLaserBestVis_Scatter']);
end
% -------------------------------------------------------------------------

% ---------------------- bar graph ----------------------------------------
% Response categorization
final_beads_vis = [sum(pos(selec&respMagn),'omitnan') sum(neg(selec&respMagn),'omitnan') sum(nochange(selec&respMagn),'omitnan');...
    sum(pos(~selec&respMagn),'omitnan') sum(neg(~selec&respMagn),'omitnan') sum(nochange(~selec&respMagn),'omitnan')];

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
title({'Response Categorization (alpha = 0.05)',['chi2 stat = ',num2str(chi2stat,2),'; p = ',num2str(pval,2)]})
ylabel('Fraction')
ylim([0 1.05]);

% Response Magnitude, average
subplot(1,3,2); hold on
[~,p_all] = ttest2(ImpactOfLaser_Vis(selec'&respMagn'), ImpactOfLaser_Vis(~selec'&respMagn'));
[~,p_inhib]= ttest2(ImpactOfLaser_Vis(~pos&~nochange&selec'&respMagn'), ImpactOfLaser_Vis(~pos&~nochange&~selec'&respMagn'));

delta_dFF(1,1) = nanmean(ImpactOfLaser_Vis(selec'&respMagn'));
delta_dFF(1,2) = nanmean(ImpactOfLaser_Vis(pos&selec'&respMagn'));
delta_dFF(1,3) = nanmean(ImpactOfLaser_Vis(~pos&~nochange&selec'&respMagn'));
delta_dFF(1,4) = nanmean(ImpactOfLaser_Vis(nochange&selec'&respMagn'));
delta_dFF(2,1) = nanmean(ImpactOfLaser_Vis(~selec'&respMagn'));
delta_dFF(2,2) = nanmean(ImpactOfLaser_Vis(pos&~selec'&respMagn'));
delta_dFF(2,3) = nanmean(ImpactOfLaser_Vis(~pos&~nochange&~selec'&respMagn'));
delta_dFF(2,4) = nanmean(ImpactOfLaser_Vis(nochange&~selec'&respMagn'));

n(1,1) = sum(~isnan(ImpactOfLaser_Vis(selec'&respMagn')));
n(1,2) = sum(~isnan(ImpactOfLaser_Vis(pos&selec'&respMagn')));
n(1,3) = sum(~isnan(ImpactOfLaser_Vis(~pos&~nochange&selec'&respMagn')));
n(1,4) = sum(~isnan(ImpactOfLaser_Vis(nochange&selec'&respMagn')));
n(2,1) = sum(~isnan(ImpactOfLaser_Vis(~selec'&respMagn')));
n(2,2) = sum(~isnan(ImpactOfLaser_Vis(pos&~selec'&respMagn')));
n(2,3) = sum(~isnan(ImpactOfLaser_Vis(~pos&~nochange&~selec'&respMagn')));
n(2,4) = sum(~isnan(ImpactOfLaser_Vis(nochange&~selec'&respMagn')));

delta_dFF_sem(1,1) = std(ImpactOfLaser_Vis(selec'&respMagn'),'omitnan')./sqrt(n(1,1));
delta_dFF_sem(1,2) = std(ImpactOfLaser_Vis(pos&selec'&respMagn'),'omitnan')./sqrt(n(1,2));
delta_dFF_sem(1,3) = std(ImpactOfLaser_Vis(~pos&~nochange&selec'&respMagn'),'omitnan')./sqrt(n(1,3));
delta_dFF_sem(1,4) = std(ImpactOfLaser_Vis(nochange&selec'&respMagn'),'omitnan')./sqrt(n(1,4));
delta_dFF_sem(2,1) = std(ImpactOfLaser_Vis(~selec'&respMagn'),'omitnan')./sqrt(n(2,1));
delta_dFF_sem(2,2) = std(ImpactOfLaser_Vis(pos&~selec'&respMagn'),'omitnan')./sqrt(n(2,2));
delta_dFF_sem(2,3) = std(ImpactOfLaser_Vis(~pos&~nochange&~selec'&respMagn'),'omitnan')./sqrt(n(2,3));
delta_dFF_sem(2,4) = std(ImpactOfLaser_Vis(nochange&~selec'&respMagn'),'omitnan')./sqrt(n(2,4));


b = bar(delta_dFF');
xticks([1 2 3 4]); xticklabels({'All','Excited','Inhibited', 'No Change'})
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
yl = ylim;
text(1,yl(2),['ttest, p = ',num2str(p_all,2)],...
    'HorizontalAlignment','center','VerticalAlignment','top')
text(3,yl(2),['ttest, p = ',num2str(p_inhib,2)],...
    'HorizontalAlignment','center','VerticalAlignment','top')
title({'Response Magnitude, Averaged'})

% Response Magnitude, all cells
subplot(1,3,3); hold on
x = ones(1,size(ImpactOfLaser_Vis(selec'&h_vis&respMagn'),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec'&h_vis&respMagn'),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 2*ones(1,size(ImpactOfLaser_Vis(selec'&h_vis&pos&respMagn'),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec'&h_vis&pos&respMagn'),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 3*ones(1,size(ImpactOfLaser_Vis(selec'&h_vis&neg&respMagn'),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec'&h_vis&neg&respMagn'),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
x = 4*ones(1,size(ImpactOfLaser_Vis(selec'&h_vis&nochange&respMagn'),1))+b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(selec'&h_vis&nochange&respMagn'),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);

x = ones(1,size(ImpactOfLaser_Vis(~selec'&h_vis&respMagn'),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec'&h_vis&respMagn'),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
x = 2*ones(1,size(ImpactOfLaser_Vis(~selec'&h_vis&pos&respMagn'),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec'&h_vis&pos&respMagn'),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
x = 3*ones(1,size(ImpactOfLaser_Vis(~selec'&h_vis&neg&respMagn'),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec'&h_vis&neg&respMagn'),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
x = 4*ones(1,size(ImpactOfLaser_Vis(~selec'&h_vis&nochange&respMagn'),1))-b(1).XOffset;
scatter(x,ImpactOfLaser_Vis(~selec'&h_vis&nochange&respMagn'),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);


for i = 1:size(delta_dFF,2)
    plot([i+b(1).XOffset-0.1 i+b(1).XOffset+0.1], [delta_dFF(1,i) delta_dFF(1,i)],'k','LineWidth',1);
    plot([i+b(2).XOffset-0.1 i+b(2).XOffset+0.1], [delta_dFF(2,i) delta_dFF(2,i)],'k','LineWidth',1);
end
xticks([1 2 3 4]);xticklabels({'All','Excited','Inhibited','No change'})
xlim([.5 4.5])
plot([.5 4.5],[0 0],':','Color',[0.5 0.5 0.5])
title('Response Magnitude, All Data')

sgtitle([{'Impact of Laser on Best Vis Resp'},{'Only Significant Cells + Magn th'}])
set(gcf,'Units','Normalized','Position',[0.01 0.45 0.95 0.45])
if saveFig
    saveas(gcf,[saveDir filesep tv '_ImpactOfLaserBestVis_Diff_Bar'],'fig');
    screen2png([saveDir filesep tv '_ImpactOfLaserBestVis_Diff_Bar']);
end

% -------------------------------------------------------------------------

%% Spontaneous activity
[~,p_all] = ttest(meanResp_spont(:,1),meanResp_spont(:,2));
[~,p_bp] = ttest(meanResp_spont(selec,1),meanResp_spont(selec,2));
[~,p_bn] = ttest(meanResp_spont(~selec,1),meanResp_spont(~selec,2));

mean_bp(1) = mean(meanResp_spont(selec,1));
mean_bp(2) = mean(meanResp_spont(selec,2));
mean_bn(1) = mean(meanResp_spont(~selec,1));
mean_bn(2) = mean(meanResp_spont(~selec,2));
sd_bp(1) = std(meanResp_spont(selec,1));
sd_bp(2) = std(meanResp_spont(selec,2));
sd_bn(1) = std(meanResp_spont(~selec,1));
sd_bn(2) = std(meanResp_spont(~selec,2));

% ---------------------- scatter plot -------------------------------------
figure;
subplot(3,3,1);hold on
s1 = scatter(squeeze(meanResp_spont(~selec,1)),squeeze(meanResp_spont(~selec,2)),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = .2;
plot([mean_bn(1)+sd_bn(1)/2 mean_bn(1)-sd_bn(1)/2],[mean_bn(2) mean_bn(2)],'Color',[1 1 1])
plot([mean_bn(1) mean_bn(1)],[mean_bn(2)+sd_bn(2)/2 mean_bn(2)-sd_bn(2)/2],'Color',[1 1 1])
xl = xlim; yl = ylim;
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
axismin = min(min(xl), min(yl));
xlim([axismin axismax]) ; ylim([axismin axismax])
plot([-axismax axismax],[-axismax  axismax ],'k:');
axis equal
title('bead -')

subplot(3,3,3);hold on
s2 = scatter(squeeze(meanResp_spont(selec,1)),squeeze(meanResp_spont(selec,2)),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2.MarkerFaceAlpha = .2;
plot([mean_bp(1)+sd_bp(1)/2 mean_bp(1)-sd_bp(1)/2],[mean_bp(2) mean_bp(2)],'Color',[1 1 1])
plot([mean_bp(1) mean_bp(1)],[mean_bp(2)+sd_bp(2)/2 mean_bp(2)-sd_bp(2)/2],'Color',[1 1 1])
xl = xlim; yl = ylim;
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
axismin = min(min(xl), min(yl));
xlim([axismin axismax]) ; ylim([axismin axismax])
plot([-axismax axismax],[-axismax  axismax ],'k:');
axis equal
title('bead +')
yyaxis right
yticks([])
ylabel({'Plotted:','All Cells'})

subplot(3,3,2);hold on
s1 = scatter(squeeze(meanResp_spont(~selec,1)),squeeze(meanResp_spont(~selec,2)),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = .2;
s2 = scatter(squeeze(meanResp_spont(selec,1)),squeeze(meanResp_spont(selec,2)),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2.MarkerFaceAlpha = .2;
plot([mean_bn(1)+sd_bn(1)/2 mean_bn(1)-sd_bn(1)/2],[mean_bn(2) mean_bn(2)],'b')
plot([mean_bn(1) mean_bn(1)],[mean_bn(2)+sd_bn(2)/2 mean_bn(2)-sd_bn(2)/2],'b')
plot([mean_bp(1)+sd_bp(1)/2 mean_bp(1)-sd_bp(1)/2],[mean_bp(2) mean_bp(2)],'r')
plot([mean_bp(1) mean_bp(1)],[mean_bp(2)+sd_bp(2)/2 mean_bp(2)-sd_bp(2)/2],'r')
xl = xlim; yl = ylim;
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
axismin = min(min(xl), min(yl));
xlim([axismin axismax]) ; ylim([axismin axismax])
plot([-axismax axismax],[-axismax  axismax ],'k:');
axis equal
title('bead - vs bead+')

subplot(3,3,4);hold on
s1 = scatter(squeeze(meanResp_spont(~selec,1)),squeeze(meanResp_spont(~selec,2)),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = .2;
plot([mean_bn(1)+sd_bn(1)/2 mean_bn(1)-sd_bn(1)/2],[mean_bn(2) mean_bn(2)],'Color',[0 0 0])
plot([mean_bn(1) mean_bn(1)],[mean_bn(2)+sd_bn(2)/2 mean_bn(2)-sd_bn(2)/2],'Color',[0 0 0])
xl = mean_bn(1)+3*sd_bn(1); yl = mean_bn(1)+3*sd_bn(1);
xl = [-xl xl]; yl=[-yl yl];
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
axismin = min(min(xl), min(yl));
xlim([axismin axismax]) ; ylim([axismin axismax])
plot([-axismax axismax],[-axismax  axismax ],'k:');
axis equal
xlabel({'Laser OFF','dF/F'}); ylabel({'Laser ON','dF/F'})

subplot(3,3,6);hold on
s2 = scatter(squeeze(meanResp_spont(selec,1)),squeeze(meanResp_spont(selec,2)),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2.MarkerFaceAlpha = .2;
plot([mean_bp(1)+sd_bp(1)/2 mean_bp(1)-sd_bp(1)/2],[mean_bp(2) mean_bp(2)],'Color',[0 0 0])
plot([mean_bp(1) mean_bp(1)],[mean_bp(2)+sd_bp(2)/2 mean_bp(2)-sd_bp(2)/2],'Color',[0 0 0])
xl = mean_bp(1)+3*sd_bp(1); yl = mean_bp(1)+3*sd_bp(1);
xl = [-xl xl]; yl=[-yl yl];
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
axismin = min(min(xl), min(yl));
xlim([axismin axismax]) ; ylim([axismin axismax])
plot([-axismax axismax],[-axismax  axismax ],'k:');
axis equal
yyaxis right
yticks([])
ylabel({'Plotted:','Mean +- 3*sd of laser off'});

subplot(3,3,5);hold on
s1 = scatter(squeeze(meanResp_spont(~selec,1)),squeeze(meanResp_spont(~selec,2)),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = .2;
s2 = scatter(squeeze(meanResp_spont(selec,1)),squeeze(meanResp_spont(selec,2)),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2.MarkerFaceAlpha = .2;
plot([mean_bn(1)+sd_bn(1)/2 mean_bn(1)-sd_bn(1)/2],[mean_bn(2) mean_bn(2)],'b')
plot([mean_bn(1) mean_bn(1)],[mean_bn(2)+sd_bn(2)/2 mean_bn(2)-sd_bn(2)/2],'b')
plot([mean_bp(1)+sd_bp(1)/2 mean_bp(1)-sd_bp(1)/2],[mean_bp(2) mean_bp(2)],'r')
plot([mean_bp(1) mean_bp(1)],[mean_bp(2)+sd_bp(2)/2 mean_bp(2)-sd_bp(2)/2],'r')
xl = mean_bn(1)+3*sd_bn(1); yl = mean_bn(1)+3*sd_bn(1);
xl = [-xl xl]; yl=[-yl yl];
xmax = max(abs(xl(1)),xl(2)); ymax = max(abs(yl(1)),yl(2));
axismax = max(xmax,ymax);
axismin = min(min(xl), min(yl));
xlim([axismin axismax]) ; ylim([axismin axismax])
plot([-axismax axismax],[-axismax  axismax ],'k:');
axis equal

subplot(3,3,7)
xl = xlim; yl = ylim;
text(xl(1),yl(1),...
    [{'mean +- sd'},...
    {['laser OFF: ',num2str(mean_bn(1),3), ' +- ', num2str(sd_bn(1),3)]},...
    {['laser ON: ',num2str(mean_bn(2),3), ' +- ', num2str(sd_bn(2),3)]},...
    {['t-test, p = ',num2str(p_bn,2), ', n = ', num2str(sum(~selec))]}]);
axis off
subplot(3,3,9)
xl = xlim; yl = ylim;
text(xl(1),yl(1),...
    [{'mean +- sd'},...
    {['laser OFF: ',num2str(mean_bp(1),3), ' +- ', num2str(sd_bp(1),3)]},...
    {['laser ON: ',num2str(mean_bp(2),3), ' +- ', num2str(sd_bp(2),3)]},...
    {['t-test, p = ',num2str(p_bp,2), ', n = ', num2str(sum(selec))]}]);
axis off
subplot(3,3,8)
xl = xlim; yl = ylim;
text(xl(1),yl(1),...
    {['pop. t-test, p = ',num2str(p_all,2), ', n = ', num2str(size(meanResp_spont,1))]});
axis off
set(gcf,'Units','Normalized','Position',[0.05 0.05 0.4 0.5])
sgtitle('Effect Of Laser on Spont Activity, All Cells')
if saveFig
    saveas(gcf,[saveDir filesep tv '_ImpactOfLaserSpont_Scatter'],'fig');
    screen2png([saveDir filesep tv '_ImpactOfLaserSpont_Scatter']);
end
% -------------------------------------------------------------------------

% ---------------------- histogram ----------------------------------------
deltaLaser{1} = meanResp_spont(~selec,1)-meanResp_spont(~selec,2);
deltaLaser{2} = meanResp_spont(selec,1)-meanResp_spont(selec,2);

figure;
subplot(3,2,1);
h1 = histogram(deltaLaser{1});
h1.EdgeColor = 'none';
ylabel('Cell #')
title('bead - ')

subplot(3,2,5);
h2 = histogram(deltaLaser{2});
h2.BinWidth = h1.BinWidth;
h2.FaceColor = color_bp;
h2.EdgeColor = 'none';
xlabel({'Net effect of laser','All data'})
xl = xlim; yl = ylim;
text(xl(1),yl(2),{'excitation'},...
    'HorizontalAlignment','Left','VerticalAlignment','Top')
text(xl(2),yl(2),{'inhibition'},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
ylabel('Cell #')
title('bead + ')

subplot(3,2,3); hold on
h3 = histogram(deltaLaser{1},'Normalization','probability');
h4 = histogram(deltaLaser{2},'Normalization','probability');
h3.BinWidth = h1.BinWidth;
h4.BinWidth = h1.BinWidth;
h3.EdgeColor = 'none';
h4.EdgeColor = 'none';
ylabel('Probability')
title('bead+ vs bead-, normalized')

subplot(3,2,2); hold on
h1 = histogram(deltaLaser{1});
h1.EdgeColor = 'none';
ylabel('Cell #')
yl = ylim;
plot([0 0],yl,'k:')
xlim([mean_bn(1)-3*sd_bn(1) mean_bn(1)+3*sd_bn(1)])
title('bead - ')

subplot(3,2,6);  hold on
h2 = histogram(deltaLaser{2});
h2.BinWidth = h1.BinWidth;
h2.FaceColor = color_bp;
h2.EdgeColor = 'none';
yl = ylim;
plot([0 0],yl,'k:')
xlabel({'Net effect of laser','Zoom in: mean +- 3*sd'})
xl = xlim; yl = ylim;
text(xl(1),yl(2),{'excitation'},...
    'HorizontalAlignment','Left','VerticalAlignment','Top')
text(xl(2),yl(2),{'inhibition'},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
ylabel('Cell #')
xlim([mean_bp(1)-3*sd_bp(1) mean_bp(1)+3*sd_bp(1)])
title('bead + ')
xl = xlim;
text(double(xl(1)),yl(2),{'excitation'},...
    'HorizontalAlignment','Left','VerticalAlignment','Top')
text(double(xl(2)),yl(2),{'inhibition'},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
sgtitle('Effect Of Laser on Spont Activity, All Cells')

subplot(3,2,4); hold on
h3 = histogram(deltaLaser{1},'Normalization','probability');
h4 = histogram(deltaLaser{2},'Normalization','probability');
h3.BinWidth = h1.BinWidth;
h4.BinWidth = h1.BinWidth;
h3.EdgeColor = 'none';
h4.EdgeColor = 'none';
ylabel('Probability')
title('bead+ vs bead-, normalized')
xlim([mean_bn(1)-3*sd_bn(1) mean_bn(1)+3*sd_bn(1)])
yl = ylim;
plot([0 0],yl,'k:')
if saveFig
    saveas(gcf,[saveDir filesep tv '_ImpactOfLaserSpont_Histo'],'fig');
    screen2png([saveDir filesep tv '_ImpactOfLaserSpont_Histo']);
end
% -------------------------------------------------------------------------

% ---------------------- bar graph ----------------------------------------
% Calculations
delta_dFF_spont(1) = mean(ImpactOfLaser_Spont(selec));
delta_dFF_spont(2) = mean(ImpactOfLaser_Spont(~selec));
n(1,1) = size(ImpactOfLaser_Spont(selec),1);
n(2,1) = size(ImpactOfLaser_Spont(~selec),1);
delta_dFF_spont_sem(1) = std(ImpactOfLaser_Spont(selec))./sqrt(n(1,1));
delta_dFF_spont_sem(2) = std(ImpactOfLaser_Spont(~selec))./sqrt(n(2,1));

% Stacked
x1 = [repmat('a',n(1,1),1); repmat('b',n(2,1),1)];
x2 = [repmat(1,final_beads(1,1),1) ; repmat(2,final_beads(1,2),1); repmat(3,final_beads(1,3),1);...
    repmat(1,final_beads(2,1),1) ; repmat(2,final_beads(2,2),1); repmat(3,final_beads(2,3),1)];
[~,chi2stat,pval] = crosstab(x1,x2);

final_beads_ratio = final_beads./[sum(selec);sum(~selec)];
figure;
subplot(1,3,1); hold on
b1 = bar(final_beads_ratio,'stacked');
legend('excited','inhibited','no change')
xticks([1 2])
xticklabels({'Beads+','Beads-'})
ylabel('Fraction')
yVal  = [0 0];
for i = 1:3
    yVal = yVal+double(b1(i).YData);
    text(double(b1(i).XData), yVal,...
        {num2str(b1(i).YData(1),2), num2str(b1(i).YData(2),2)},...
        'HorizontalAlignment','center','VerticalAlignment','top','Color',[1 1 1])
end
text(1,1,['n = ',num2str(sum(selec))],...
    'HorizontalAlignment','Center','VerticalAlignment','bottom');
text(2,1,['n = ',num2str(sum(~selec))],...
    'HorizontalAlignment','Center','VerticalAlignment','bottom');
ylim([0  1.05]);
title({'Response Categorization (alpha = 0.05)',['chi2 stat = ', num2str(chi2stat,2), ' ; p = ', num2str(pval,2)]})


% Box whisker plot
% boxplot_data = [ImpactOfLaser_Spont(selec) ; ImpactOfLaser_Spont(~selec)];
% boxplot_grouping = [ones(size(ImpactOfLaser_Spont(selec),1),1); 2*ones(size(ImpactOfLaser_Spont(~selec),1),1)];
% subplot(1,3,2); hold on
% boxplot(boxplot_data, boxplot_grouping,...
%     'Notch','on','Labels',{'bead+','bead-'});
% xl = xlim;
% plot([xl(1) xl(2)],[0 0],'k--');

% Bar graph
[~,p] = ttest2(ImpactOfLaser_Spont(selec),ImpactOfLaser_Spont(~selec));

subplot(1,3,2); hold on
b2 = bar(delta_dFF_spont);
for i = 1:2
    plot([i i],[delta_dFF_spont(i)-delta_dFF_spont_sem(i) delta_dFF_spont(i)+delta_dFF_spont_sem(i)],'k-');
end
xticks([1 2]); xticklabels({'Beads+','Beads-'}); xtickangle(45)
b2.FaceColor = 'flat';
b2.CData(1,:) = color_bp;
b2.CData(2,:) = color_bn;
text([1 2],double(b2.YData),string([n(1,1) n(2,1)]),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Color',[0 0 0])
ylabel('delta dFF')
title({'Magnitude of the effect, averaged',['ttest, p = ',num2str(p,2)]})

subplot(1,3,3); hold on
scatter(ones(size(ImpactOfLaser_Spont(selec),1),1),ImpactOfLaser_Spont(selec),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
% s1.MarkerFaceAlpha = .2;
scatter(2*ones(size(ImpactOfLaser_Spont(~selec),1),1),ImpactOfLaser_Spont(~selec),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
xlim([0 3])
plot(xlim,[0 0],'k:')
xticks([1 2]); xticklabels({'Beads+','Beads-'}); xtickangle(45)
title('Magnitude of the effect, all cells')

sgtitle('Impact of Laser on Spontaneous Activity')
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.95 0.4])
if saveFig
    saveas(gcf,[saveDir filesep tv '_ImpactOfLaserSpont_Diff_Bar'],'fig');
    screen2png([saveDir filesep tv '_ImpactOfLaserSpont_Diff_Bar']);
end
% -------------------------------------------------------------------------

%% z-score (see Isacson's preprint on contralateral auditory axons)
for i = 1:length(data_list)
    timeStimInd = data_list{1, i}.Info.timeStimInd;
    timeStimInd = timeStimInd(1:size(data_list{i}.bestVisResp.selVisStim,2));
    timeBaseInd = data_list{1, 2}.Info.timeBaseInd;
    timeBaseInd = timeBaseInd(1:size(data_list{i}.bestVisResp.selVisStim,2));
    vis{i} = squeeze(mean(data_list{i}.bestVisResp.selVisStim(:,timeStimInd,:,1:20),2));
    spont{i} = squeeze(mean(data_list{i}.bestVisResp.selVisStim(:,timeBaseInd,:,1:20),2));
%     spont{i} = data_list{i}.spont.allTrials(h_vis,timeStimInd,:,1:20);
end
vis = cat(1,vis{:});
spont = cat(1,spont{:});

vis_mean_bp = mean(vis(selec&respMagn,:,:),3);
vis_mean_bn = mean(vis(~selec&respMagn,:,:),3);
vis_sd_bp = std(vis(selec&respMagn,:,:),[],3);
vis_sd_bn = std(vis(~selec&respMagn,:,:),[],3);
spont_mean_bp = mean(spont(selec&respMagn,:,:),3);
spont_mean_bn = mean(spont(~selec&respMagn,:,:),3);
spont_sd_bp = std(spont(selec&respMagn,:,:),[],3);
spont_sd_bn = std(spont(~selec&respMagn,:,:),[],3);

z_score_bp = (vis_mean_bp - spont_mean_bp)./sqrt((vis_sd_bp.^2)./20+(spont_sd_bp.^2)./20);
z_score_bn = (vis_mean_bn - spont_mean_bn)./sqrt((vis_sd_bn.^2)./20+(spont_sd_bn.^2)./20);
zsc_mean_bp = mean(z_score_bp);
zsc_mean_bn = mean(z_score_bn);
[~,p1]=ttest(z_score_bp(:,1),z_score_bp(:,2));
[~,p2]=ttest(z_score_bn(:,1),z_score_bn(:,2));

figure;
subplot(3,3,[1,2,4,5]); hold on
s1 = scatter(z_score_bn(:,1),z_score_bn(:,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = 0.5;
s2 = scatter(z_score_bp(:,1),z_score_bp(:,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2.MarkerFaceAlpha = 0.5;
xlabel('Laser OFF'); ylabel('Laser ON');
xl = xlim; yl = ylim;
axismax = max(xl(2),yl(2));
axismin = min(xl(1),yl(1));
plot([axismin axismax],[axismin axismax],'k--')
plot(zsc_mean_bn(1),zsc_mean_bn(2),'b+')
plot(zsc_mean_bp(1),zsc_mean_bp(2),'r+')
title({'z-score'})
subplot(3,3,7);
xl = xlim; yl = ylim;
text(xl(1),yl(2)/2,[{'bead-'},...
    {['ttest, p = ',num2str(p2,3)]},...
    {'mean +- sd'},...
    {['Laser OFF: ',num2str(zsc_mean_bn(1),3)]},...
    {['Laser ON: ',num2str(zsc_mean_bn(2),3)]}]);
axis off
subplot(3,3,8);
xl = xlim; yl = ylim;
text(xl(1),yl(2)/2,[{'bead+'},...
    {['ttest, p = ',num2str(p1,3)]},...
    {'mean +- sd'},...
    {['Laser OFF: ',num2str(zsc_mean_bp(1),3)]},...
    {['Laser ON: ',num2str(zsc_mean_bp(2),3)]}]);
axis off

subplot(3,3,3); hold on
zsc_diff_bn = z_score_bn(:,1) - z_score_bn(:,2);
h1 = histogram(zsc_diff_bn);
h1.EdgeColor = 'none'; h1.FaceColor = color_bn;
yl = ylim;
plot([0 0],yl,'k:')
ylabel('Cell #')
title('bead-')

subplot(3,3,9); hold on
zsc_diff_bp = z_score_bp(:,1) - z_score_bp(:,2);
h2 = histogram(zsc_diff_bp);
h2.EdgeColor = 'none'; h2.FaceColor = color_bp;
yl = ylim;
plot([0 0],yl,'k:')
ylabel('Cell #')
title('bead+')

subplot(3,3,6); hold on
[~,p] = ttest2(zsc_diff_bn,zsc_diff_bp);
h1 = histogram(zsc_diff_bn,'Normalization','probability');
h2 = histogram(zsc_diff_bp,'Normalization','probability');
h2.BinWidth = h1.BinWidth;
h1.EdgeColor = 'none'; h1.FaceColor = color_bn;
h2.EdgeColor = 'none'; h2.FaceColor = color_bp;
ylabel('Probability')
title({['ttest, p = ',num2str(p,3)]})
yl = ylim;
plot([0 0],yl,'k:')

sgtitle('Discriminability. LaserOFF vs LaserON')
set(gcf,'Position',[680   429   576   549])
if saveFig
    saveas(gcf,[saveDir filesep tv '_ImpactOfLaser_zScore'],'fig');
    screen2png([saveDir filesep tv '_ImpactOfLaser_zScore']);
end

%% Adaptation
ImpactOfLaser = squeeze(vis(:,2,:) - vis(:,1,:));
ImpactOfLaser_mean_bp = mean(ImpactOfLaser(respMagn&selec,:));
ImpactOfLaser_mean_bn = mean(ImpactOfLaser(respMagn&~selec,:));
ImpactOfLaser_sd_bp = std(ImpactOfLaser(respMagn&selec,:));
ImpactOfLaser_sd_bn = std(ImpactOfLaser(respMagn&~selec,:));
ImpactOfLaser_sem_bp = ImpactOfLaser_sd_bp./sqrt(size(ImpactOfLaser(respMagn&selec,:),1));
ImpactOfLaser_sem_bn = ImpactOfLaser_sd_bn./sqrt(size(ImpactOfLaser(respMagn&~selec,:),1));

figure;
subplot(1,2,1);hold on
plot(ImpactOfLaser_mean_bp,'Color',color_bp)
plot(ImpactOfLaser_mean_bn,'Color',color_bn)
% fill([0:20],[ImpactOfLaser_mean_bn-ImpactOfLaser_sem_bn flip(ImpactOfLaser_mean_bn+ImpactOfLaser_sem_bn)],'Color',color_bn)
plot([ImpactOfLaser_mean_bn-ImpactOfLaser_sem_bn],'b:')
plot([ImpactOfLaser_mean_bn+ImpactOfLaser_sem_bn],'b:')
plot([ImpactOfLaser_mean_bp-ImpactOfLaser_sem_bp],'r:')
plot([ImpactOfLaser_mean_bp+ImpactOfLaser_sem_bp],'r:')

subplot(1,2,2); hold on
% binned_bp(1) = mean(ImpactOfLaser_mean_bp(1:10));
% binned_bp(2) = mean(ImpactOfLaser_mean_bp(11:20));
% binned_bn(1) = mean(ImpactOfLaser_mean_bn(1:10));
% binned_bn(2) = mean(ImpactOfLaser_mean_bn(11:20));
% bar([binned_bn; binned_bp]')
binned_bp_test(:,1) = mean(ImpactOfLaser(respMagn&selec,1:10),2);
binned_bp_test(:,2) = mean(ImpactOfLaser(respMagn&selec,11:20),2);
binned_bp_test_mean = mean(binned_bp_test);
binned_bp_test_sem = std(binned_bp_test)./sqrt(size(binned_bp_test,1));
binned_bn_test(:,1) = mean(ImpactOfLaser(respMagn&~selec,1:10),2);
binned_bn_test(:,2) = mean(ImpactOfLaser(respMagn&~selec,11:20),2);
binned_bn_test_mean = mean(binned_bn_test);
binned_bn_test_sem = std(binned_bn_test)./sqrt(size(binned_bn_test,1));

b1 = bar([binned_bn_test_mean ; binned_bp_test_mean]');
drawnow
offset = b1.XOffset;
plot([1+offset 1+offset],...
    [binned_bn_test_mean(1)-binned_bn_test_sem(1) binned_bn_test_mean(1)+binned_bn_test_sem(1)],...
    'k-')
plot([2+offset 2+offset],...
    [binned_bn_test_mean(2)-binned_bn_test_sem(2) binned_bn_test_mean(2)+binned_bn_test_sem(2)],...
    'k-')
plot([1-offset 1-offset],...
    [binned_bp_test_mean(1)-binned_bp_test_sem(1) binned_bp_test_mean(1)+binned_bp_test_sem(1)],...
    'k-')
plot([2-offset 2-offset],...
    [binned_bp_test_mean(2)-binned_bp_test_sem(2) binned_bp_test_mean(2)+binned_bp_test_sem(2)],...
    'k-')
%%
diary off;
end