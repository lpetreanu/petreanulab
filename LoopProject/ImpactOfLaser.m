% Identify the beads+ and beads- cells first
% Select ROIs based on "iscell" and plot the impact of laser based on their
% projection ID (i.e, beads+ or -)

%%%%%%%%%%%%%%%%%%% -- Preambule -- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except data
mainDir = data.Info.mainDir;
% mainDir = uigetdir('E:\TempData\CMloop6\20200121_day8\CMloop6_day8_L2', 'main directory');
response_lag = 1; % (s). For selecting the cells using a ttest. time after vis stim OFFSET to use
plot_indivROI = 1;
saveIndivROI_Flag = 1;
saveFigFlag = 1;
saveVarFlag = 1;
a = 0.01; % for t-test (only for the one detecting vis respon atm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
if saveIndivROI_Flag
    saveDir = [mainDir, '\Figures\IndivROI'];
    disp(mainDir)
    if ~ exist(saveDir)
        mkdir(saveDir)
    end
end
%%%%%%%%%%%%%%%%%%%% ------------------------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% basic calculations
% beads = [data.cellsID(~isnan(data.cellsID(:,2)),1) data.cellsID(~isnan(data.cellsID(:,2)),2)];
% beads_pos2 = [data.cellsID(~isnan(data.cellsID(:,2)),1)];
% isnan(data.cellsID(:,2)) = 0
% beads_temp = data.cellsID(:,2);
for i = 1:size(data.cellsID(:,2),1)
    if isnan(data.cellsID(i,2))
        beads_pos(i) = false;
        beads_neg(i) = false;
    else
        beads_pos(i) = logical(data.cellsID(i,2));
        beads_neg(i) = logical(~data.cellsID(i,2));
    end
end

cells = data.s2p.iscell;

Vis_Laser_mean = mean(data.visNopto.Vis_Laser(cells,:,:,:),4);
Vis_Laser_sem = std(data.visNopto.Vis_Laser(cells,:,:,:),[],4)/sqrt(size(data.visNopto.Vis_Laser(cells,:,:,:),4));

Vis_NoLaser_mean = mean(data.visNopto.Vis_NoLaser(cells,:,:,:),4);
Vis_NoLaser_sem = std(data.visNopto.Vis_NoLaser(cells,:,:,:),[],4)/sqrt(size(data.visNopto.Vis_NoLaser(cells,:,:,:),4));

Vis_Laser_meanResp = squeeze(mean(Vis_Laser_mean(:,data.Info.timeStimInd,:),2));
Vis_Laser_meanResp_sem = squeeze(mean(Vis_Laser_sem(:,data.Info.timeStimInd,:),2));

Vis_NoLaser_meanResp = squeeze(mean(Vis_NoLaser_mean(:,data.Info.timeStimInd,:),2));
Vis_NoLaser_meanResp_sem = squeeze(mean(Vis_NoLaser_sem(:,data.Info.timeStimInd,:),2));

Vis_meanResp_diff = Vis_Laser_meanResp - Vis_NoLaser_meanResp;
Vis_meanResp_ratio = Vis_Laser_meanResp ./ Vis_NoLaser_meanResp;

%% Population average for all the cells
Pop_avgDiff = mean(Vis_meanResp_diff,1);
Pop_semDiff = std(Vis_meanResp_diff,[],1)./sqrt(size(Vis_meanResp_diff,1));

Pop_avgDiv = mean(Vis_meanResp_ratio,1);
Pop_semDiv = std(Vis_meanResp_ratio,[],1)./sqrt(size(Vis_meanResp_ratio,1));

%% Population averages for beads+ vs beads-
% % average across cells
% Beads_pos_Pop_avgDiff = mean(Vis_meanResp_diff(beads_pos,:),1);
% Beads_pos_Pop_semDiff = std(Vis_meanResp_diff(beads_pos,:),[],1)./sqrt(size(Vis_meanResp_diff(beads_pos,:),1));
% Beads_neg_Pop_avgDiff = mean(Vis_meanResp_diff(beads_neg,:),1);
% Beads_neg_Pop_semDiff = std(Vis_meanResp_diff(beads_neg,:),[],1)./sqrt(size(Vis_meanResp_diff(beads_neg,:),1));
% 
% figure; hold on; 
% bar([Beads_neg_Pop_avgDiff; Beads_pos_Pop_avgDiff]')
% plot([1 1],[Beads_neg_Pop_avgDiff-Beads_neg_Pop_semDiff Beads_neg_Pop_avgDiff+Beads_neg_Pop_semDiff],'k');
% plot([2 2],[Beads_pos_Pop_avgDiff-Beads_neg_Pop_semDiff Beads_pos_Pop_avgDiff+Beads_neg_Pop_semDiff],'k')
% 
% % average across orientation
% BeadsPos_PopGrandAvg_Diff = mean(Beads_pos_Pop_avgDiff);
% BeadsNeg_PopGrandAvg_Diff = mean(Beads_neg_Pop_avgDiff);
% figure; hold on
% bar([BeadsNeg_PopGrandAvg_Diff BeadsPos_PopGrandAvg_Diff]);

% average difference per cell, then average across cells
BeadsPos_Avg_Diff = mean(Vis_meanResp_diff(beads_pos,:),2);
BeadsNeg_Avg_Diff = mean(Vis_meanResp_diff(beads_neg,:),2);
[h, p] = ttest2(BeadsPos_Avg_Diff, BeadsNeg_Avg_Diff)

BeadsPos_PopAvg_Diff = mean(BeadsPos_Avg_Diff);
BeadsPos_PopSem_Diff = std(BeadsPos_Avg_Diff)/sqrt(size(BeadsPos_Avg_Diff,1));
BeadsNeg_PopAvg_Diff = mean(BeadsNeg_Avg_Diff);
BeadsNeg_PopSem_Diff = std(BeadsNeg_Avg_Diff)/sqrt(size(BeadsNeg_Avg_Diff,1));

figure; 
subplot(1,2,1); hold on; 
b1 = bar(1,BeadsPos_PopAvg_Diff,'r');
b2 = bar(2,BeadsNeg_PopAvg_Diff,'b');
plot([0.9 0.9],[BeadsPos_PopAvg_Diff-BeadsPos_PopSem_Diff BeadsPos_PopAvg_Diff+BeadsPos_PopSem_Diff],'k');
plot([1.9 1.9],[BeadsNeg_PopAvg_Diff-BeadsNeg_PopSem_Diff BeadsNeg_PopAvg_Diff+BeadsNeg_PopSem_Diff],'k');
xticks([1 2])
xticklabels({'beads+','beads-'})
% legend([b1 b2],{'beads+','beads-'})
plot(1.1,BeadsPos_Avg_Diff,'ko');
plot(2.1,BeadsNeg_Avg_Diff,'ko');
yl = ylim;
text(1,yl(2),num2str(size(BeadsPos_Avg_Diff,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top')
text(2,yl(2),num2str(size(BeadsNeg_Avg_Diff,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top')
title(num2str(p,3))
sgtitle('Effect of opto on the averaged visual stim')

subplot(1,2,2); hold on; 
BeadsPos_Avg(:,1) = mean(Vis_NoLaser_meanResp(beads_pos,:),2);
BeadsPos_Avg(:,2) = mean(Vis_Laser_meanResp(beads_pos,:),2);
BeadsNeg_Avg(:,1) = mean(Vis_NoLaser_meanResp(beads_neg,:),2);
BeadsNeg_Avg(:,2) = mean(Vis_Laser_meanResp(beads_neg,:),2);
plot(1,BeadsPos_Avg(:,1),'ro')
plot(2,BeadsPos_Avg(:,2),'ro')
plot(3,BeadsNeg_Avg(:,1),'bo')
plot(4,BeadsNeg_Avg(:,2),'bo')
xlim([0 5])
xticks([1 2 3 4])
xticklabels({'No Opto','Opto','No Opto','Opto'})
xtickangle(45)
for i = 1:size(BeadsPos_Avg,1)
    plot([1 2],[BeadsPos_Avg(i,1), BeadsPos_Avg(i,2)],'k-');
    
end
for i = 1:size(BeadsNeg_Avg,1)
    plot([3 4],[BeadsNeg_Avg(i,1), BeadsNeg_Avg(i,2)],'k-');
end

figure; hold on
plot(BeadsPos_Avg(:,1),BeadsPos_Avg(:,2),'ro')
plot(BeadsNeg_Avg(:,1),BeadsNeg_Avg(:,2),'bo')
yl = ylim;
xl = xlim;
plot([xl(1) -xl(1)],[yl(1) -yl(1)],'k-')
plot([xl(1) xl(2)],[0 0],'k:')
plot([0 0],[yl(1) yl(2)],'k:')
xlabel('Grating, dF/F')
ylabel('Grating+Opt, dF/F')
title('Impact of laser on each cell s average response to grating')

figure;
subplot(1,2,1); hold on
plot(Vis_NoLaser_meanResp(beads_pos,:),Vis_Laser_meanResp(beads_pos,:),'ro')
plot(Vis_NoLaser_meanResp(beads_neg,:),Vis_Laser_meanResp(beads_neg,:),'bo')
yl = ylim;
xl = xlim;
ax_max = max(max(abs(yl(1)),abs(yl(2))),max(abs(xl(1)),abs(xl(2))));
plot([-ax_max ax_max],[-ax_max ax_max],'k-')
axis equal
plot([-ax_max ax_max],[0 0],'k:')
plot([0 0],[-ax_max ax_max],'k:')
xlabel('Grating, dF/F')
ylabel('Grating+Opt, dF/F')
subplot(1,2,2); hold on
plot(Vis_NoLaser_meanResp(beads_pos,:),Vis_meanResp_diff(beads_pos,:),'ro')
plot(Vis_NoLaser_meanResp(beads_neg,:),Vis_meanResp_diff(beads_neg,:),'bo')
yl = ylim;
xl = xlim;
ax_max = max(max(abs(yl(1)),abs(yl(2))),max(abs(xl(1)),abs(xl(2))));
plot([ax_max -ax_max],[-ax_max ax_max],'k--')
axis equal
plot([-ax_max ax_max],[0 0],'k:')
plot([0 0],[-ax_max ax_max],'k:')
xlabel('Grating, dF/F')
ylabel('Grating+Opt - Grating, dF/F')
sgtitle('Impact of laser on each cell-grating pair')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effect on best stim only
[bestResp, index] = sort(Vis_NoLaser_meanResp,2, 'descend');
for i = 1:size(Vis_NoLaser_meanResp,1)
    Best_NoLaser(i,1) = Vis_NoLaser_meanResp(i,index(i,1));
    Best_Laser(i,1) = Vis_Laser_meanResp(i,index(i,1));
end
BeadsPos_Best_NoLaser = Best_NoLaser(beads_pos);
BeadsPos_Best_Laser = Best_Laser(beads_pos);
BeadsNeg_Best_NoLaser = Best_NoLaser(beads_neg);
BeadsNeg_Best_Laser = Best_Laser(beads_neg);

BeadsPos_BestStim_Diff = BeadsPos_Best_Laser - BeadsPos_Best_NoLaser;
BeadsNeg_BestStim_Diff = BeadsNeg_Best_Laser - BeadsNeg_Best_NoLaser;
% BeadsPos_BestStim_Diff = bestResp(beads_pos,1);
% BeadsNeg_BestStim_Diff = bestResp(beads_neg,1);
[h, p] = ttest2(BeadsPos_BestStim_Diff, BeadsNeg_BestStim_Diff)

BeadsPos_BestStim_PopAvg_Diff = mean(BeadsPos_BestStim_Diff);
BeadsPos_BestStim_PopSem_Diff = std(BeadsPos_BestStim_Diff)/sqrt(size(BeadsPos_BestStim_Diff,1));
BeadsNeg_BestStim_PopAvg_Diff = mean(BeadsNeg_BestStim_Diff);
BeadsNeg_BestStim_PopSem_Diff = std(BeadsNeg_BestStim_Diff)/sqrt(size(BeadsNeg_BestStim_Diff,1));

figure; 
subplot(1,2,1); hold on
bar(1,BeadsPos_BestStim_PopAvg_Diff,'r');
bar(2,BeadsNeg_BestStim_PopAvg_Diff,'b');
plot([1 1],[BeadsPos_BestStim_PopAvg_Diff-BeadsPos_BestStim_PopSem_Diff BeadsPos_BestStim_PopAvg_Diff+BeadsPos_BestStim_PopSem_Diff],'k')
plot([2 2],[BeadsNeg_BestStim_PopAvg_Diff-BeadsNeg_BestStim_PopSem_Diff BeadsNeg_BestStim_PopAvg_Diff+BeadsNeg_BestStim_PopSem_Diff],'k')
xticks([1 2])
xticklabels({'beads+','beads-'})
% legend([b1 b2],{'beads+','beads-'})
plot(1.1,BeadsPos_BestStim_Diff,'ko');
plot(2.1,BeadsNeg_BestStim_Diff,'ko');
yl = ylim;
text(1,yl(2),num2str(size(BeadsPos_Avg_Diff,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top')
text(2,yl(2),num2str(size(BeadsNeg_Avg_Diff,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top')
title(['Difference, p = ',num2str(p,3)])
ylabel('Grating+Opto - Opto, dF/F')
subplot(1,2,2); hold on
plot(1,BeadsPos_Best_NoLaser,'ro')
plot(2,BeadsPos_Best_Laser,'ro')
plot(3,BeadsNeg_Best_NoLaser,'bo')
plot(4,BeadsNeg_Best_Laser,'bo')
xlim([0 5])
xticks([1 2 3 4])
xticklabels({'Grating','Grating+Opto','Grating','Grating+Opto'})
xtickangle(45)
for i = 1:size(BeadsPos_Best_NoLaser,1)
    plot([1 2],[BeadsPos_Best_NoLaser(i), BeadsPos_Best_Laser(i)],'k-');
end
for i = 1:size(BeadsNeg_Best_NoLaser,1)
    plot([3 4],[BeadsNeg_Best_NoLaser(i), BeadsNeg_Best_Laser(i)],'k-');
end
ylabel('dF/F')
title('Raw responses')
sgtitle('Effect of opto on best visual stim')

figure;
subplot(1,2,1); hold on
plot(BeadsPos_Best_NoLaser,BeadsPos_Best_Laser,'ro')
plot(BeadsNeg_Best_NoLaser,BeadsNeg_Best_Laser,'bo')
yl = ylim;
xl = xlim;
ax_max = max(max(abs(yl(1)),abs(yl(2))),max(abs(xl(1)),abs(xl(2))));
plot([-ax_max ax_max],[-ax_max ax_max],'k-')
axis equal
plot([-ax_max ax_max],[0 0],'k:')
plot([0 0],[-ax_max ax_max],'k:')
% plot([xl(1) xl(2)],[0 0],'k:')
% plot([0 0],[yl(1) yl(2)],'k:')
xlabel('Grating, dF/F')
ylabel('Grating+Opt, dF/F')
subplot(1,2,2); hold on
plot(BeadsPos_Best_NoLaser,BeadsPos_BestStim_Diff,'ro')
plot(BeadsNeg_Best_NoLaser,BeadsNeg_BestStim_Diff,'bo')
yl = ylim;
xl = xlim;
ax_max = max(max(abs(yl(1)),abs(yl(2))),max(abs(xl(1)),abs(xl(2))));
plot([-ax_max ax_max],[ax_max -ax_max],'k-')
axis equal
plot([-ax_max ax_max],[0 0],'k:')
plot([0 0],[-ax_max ax_max],'k:')
xlabel('Grating, dF/F')
ylabel('Grating+Opt - Grating, dF/F')
sgtitle('Impact of laser on each cell s best response to grating')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effect on worst stim only (proxy for spont. activity)
BeadsPos_WorstStim_Diff = bestResp(beads_pos,8);
BeadsNeg_WorstStim_Diff = bestResp(beads_neg,8);

BeadsPos_WorstStim_PopAvg_Diff = mean(BeadsPos_WorstStim_Diff);
BeadsPos_WorstStim_PopSem_Diff = std(BeadsPos_WorstStim_Diff)/sqrt(size(BeadsPos_WorstStim_Diff,1));
BeadsNeg_WorstStim_PopAvg_Diff = mean(BeadsNeg_WorstStim_Diff);
BeadsNeg_WorstStim_PopSem_Diff = std(BeadsNeg_WorstStim_Diff)/sqrt(size(BeadsNeg_WorstStim_Diff,1));

figure; hold on
bar(1,BeadsPos_WorstStim_PopAvg_Diff,'r');
bar(2,BeadsNeg_WorstStim_PopAvg_Diff,'b');
plot([1 1],[BeadsPos_WorstStim_PopAvg_Diff-BeadsPos_WorstStim_PopSem_Diff BeadsPos_WorstStim_PopAvg_Diff+BeadsPos_WorstStim_PopSem_Diff],'k')
plot([2 2],[BeadsNeg_WorstStim_PopAvg_Diff-BeadsNeg_WorstStim_PopSem_Diff BeadsNeg_WorstStim_PopAvg_Diff+BeadsNeg_WorstStim_PopSem_Diff],'k')
xticks([1 2])
xticklabels({'beads+','beads-'})
% legend([b1 b2],{'beads+','beads-'})
yl = ylim;
text(1,yl(2),num2str(size(BeadsPos_Avg_Diff,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top')
text(2,yl(2),num2str(size(BeadsNeg_Avg_Diff,1)),...
    'HorizontalAlignment','center','VerticalAlignment','top')
title('Effect of opto on worst visual stim (proxy for spont. activity')

[h, p] = ttest2(BeadsPos_WorstStim_Diff, BeadsNeg_WorstStim_Diff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% timeVect = data.Info.timeVect;
% xvect = data.Info.xvect;
% 
% directions = data.Info.directions; % directions = [0 45 90 135 180 225 270 315];
% 
% % Sort the data by prefered visual stimulus w/out laser
% [bestResp, index] = sort(Vis_NoLaser_meanResp,2, 'descend');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;bar([mean(Vis_Laser_meanResp,1);mean(Vis_NoLaser_meanResp,1)]')
% cells = find(iscell(:,1));
% cells(:,2) = NaN(size(cells,1),1);
