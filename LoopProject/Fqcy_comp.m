% For all ROIs
% Plot a bunch of graph, effect of opto on each direction stimulus
% Select the cells based on their visual responses
% Sort by best stimulus
% Effect of opto on spontaneous activity
% etc...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRepPerStimType = 10;
nStinCond = 3;
selection = true; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainDir = data.Info.mainDir;
% mainDir = uigetdir('E:\TempData\CMloop6\20200121_day8\CMloop6_day8_L2', 'main directory');
saveFigFlag = 1;
saveVarFlag = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
if saveFigFlag
    saveDir = [mainDir, '\Figures\Fqcy_comp'];
    disp(mainDir)
    if ~ exist(saveDir)
        mkdir(saveDir)
    end
end

timeVect = data.Info.timeVect;
xvect = data.Info.xvect;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directions = data.Info.directions; % directions = [0 45 90 135 180 225 270 315];
if selection
    ROIs = logical(iscell(:,1));
else
    ROIs = logical(ones(size(iscell,1),1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Opto on Vis Stim
clear Vis_Laser_Trialmean Vis_NoLaser_Trialmean Vis_Laser_meanResp Vis_NoLaser_meanResp
counter = 1;
for i = 1:nStinCond
    Vis_Laser_Trialmean(:,:,:,i) = mean(data.visNopto.Vis_Laser(ROIs,:,:,counter:i*nRepPerStimType),4);
    Vis_NoLaser_Trialmean(:,:,:,i) = mean(data.visNopto.Vis_NoLaser(ROIs,:,:,counter:i*nRepPerStimType),4);
    counter = counter + nRepPerStimType;
end
Vis_Laser_meanResp = squeeze(mean(Vis_Laser_Trialmean(:,data.Info.timeStimInd,:,:),2));
Vis_NoLaser_meanResp = squeeze(mean(Vis_NoLaser_Trialmean(:,data.Info.timeStimInd,:,:),2));

Vis_meanResp_diff = Vis_Laser_meanResp - Vis_NoLaser_meanResp;

Pop_avgDiff = squeeze(mean(Vis_meanResp_diff,1));
Pop_semDiff = squeeze(std(Vis_meanResp_diff,[],1)./sqrt(size(Vis_meanResp_diff,1)));

figure;
cmean = mean2(squeeze(Vis_meanResp_diff(:,:,1)));
cmax = cmean+5*std2(squeeze(Vis_meanResp_diff(:,:,1)));
for i = 1:nStinCond
subplot(1,nStinCond,i); 
imagesc(squeeze(Vis_meanResp_diff(:,:,i)),[- cmax cmax]);
colormap('diverging_bwr_40_95_c42_n256');
end
if saveFigFlag
    saveas(gcf,[saveDir '\AllROIs_Vis_HM'],'fig');
    screen2png([mainDir '\AllROIs_Vis_HM']);
end

figure;
subplot(1,2,1); hold on;
bar(directions,Pop_avgDiff);
yl = ylim;
subplot(1,2,2); hold on;
bar(mean(Pop_avgDiff));
% plot([1 1],[mean(Pop_avg40Diff)-std(Pop_avg40Diff)./sqrt(size(Pop_avg40Diff,2)) mean(Pop_avg40Diff)+std(Pop_avg40Diff)./sqrt(size(Pop_avg40Diff,2))],'k')
% plot([2 2],[mean(Pop_avg80Diff)-std(Pop_avg80Diff)./sqrt(size(Pop_avg80Diff,2)) mean(Pop_avg80Diff)+std(Pop_avg80Diff)./sqrt(size(Pop_avg80Diff,2))],'k')
ylim(yl);
if saveFigFlag
    saveas(gcf,[saveDir '\AllROIs_Vis_Histo'],'fig');
    screen2png([mainDir '\AllROIs_Vis_Histo']);
end

%% Opto on spont
clear NoVis_Laser_Trialmean NoVis_NoLaser_Trialmean NoVis_Laser_meanResp NoVis_NoLaser_meanResp
counter = 1;
for i = 1:nStinCond
    NoVis_Laser_Trialmean(:,:,i) = mean(data.OptoOnly.Laser_NoVis(ROIs,:,counter:i*nRepPerStimType),3);
    NoVis_NoLaser_Trialmean(:,:,i) = mean(data.OptoOnly.NoLaser_NoVis(ROIs,:,counter:i*nRepPerStimType),3);
    counter = counter + nRepPerStimType;
end

NoVis_Laser_meanResp = squeeze(mean(NoVis_Laser_Trialmean(:,data.Info.timeStimInd,:),2));
NoVis_NoLaser_meanResp = squeeze(mean(NoVis_NoLaser_Trialmean(:,data.Info.timeStimInd,:),2));

NoVis_Laser_meanResp_diff = NoVis_Laser_meanResp - NoVis_NoLaser_meanResp;

Pop_avg_NoVis_Diff = mean(NoVis_Laser_meanResp_diff);
Pop_sem_NoVis_Diff = std(NoVis_Laser_meanResp_diff)./sqrt(size(NoVis_Laser_meanResp_diff,2));

% % to edit
% figure;
% cmean = max(mean2(NoVis_meanResp40_diff),mean2(NoVis_Laser_meanResp80_diff));
% cmax = cmean+5*std2(NoVis_meanResp40_diff);
% subplot(1,2,1); imagesc(NoVis_meanResp40_diff,[- cmax cmax]);
% colormap('diverging_bwr_40_95_c42_n256'); title('40 W')
% subplot(1,2,2); imagesc(NoVis_Laser_meanResp80_diff,[- cmax cmax]);
% colormap('diverging_bwr_40_95_c42_n256'); title('80 W')
%     saveas(gcf,[saveDir '\AllROIs_NoVis_HM'],'fig');
%     screen2png([mainDir '\AllROIs_NoVis_HM']);
%

figure; hold on;
bar(Pop_avg_NoVis_Diff);
% % to edit
% plot([1 1],[Pop_avg_NoVis_40Diff-Pop_sem_NoVis_40Diff Pop_avg_NoVis_40Diff+Pop_sem_NoVis_40Diff],'k')
% plot([2 2],[Pop_avg_NoVis_80Diff-Pop_sem_NoVis_80Diff Pop_avg_NoVis_80Diff+Pop_sem_NoVis_80Diff],'k')
%     saveas(gcf,[saveDir '\AllROIs_NoVis_Histo'],'fig');
%     screen2png([mainDir '\AllROIs_NoVis_Histo']);
%     
%     figure; hold on
%     plot(NoVis_NoLaser80_meanResp,NoVis_Laser80_meanResp,'bo')
%     plot(NoVis_NoLaser40_meanResp,NoVis_Laser40_meanResp,'ro')
%     axis equal
%     plot([-0.1 0.1],[-0.1 0.1])
%     
%     
%         figure; hold on
%     plot(Vis_NoLaser80_meanResp,Vis_Laser80_meanResp,'bo')
%     plot(Vis_NoLaser40_meanResp,Vis_Laser40_meanResp,'ro')
%     axis equal
%     plot([-0.1 2],[-0.1 2])
%     
%     figure; hold on
%     for i = 1:size(Vis_meanResp80_diff,1)
%        plot([1  2],[mean(Vis_meanResp80_diff(i,:),2) mean(Vis_meanResp40_diff(i,:),2)],'k-')
%     end
%     xlim([0 3])
    