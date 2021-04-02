%%%%   For ALL ROIs (not selecting based on iscell) %%%%%%%%%%%%%%%%%%%%%%%
%  Plot a bunch of graph
%      - effect of opto on each direction stimulus
%      - Select the cells based on their visual responses
%      - Sort by best stimulus
%      - Effect of opto on spontaneous activity
%      - etc...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except data
mainDir = data.Info.mainDir;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vis_Laser_mean = mean(data.visNopto.Vis_Laser,4);
Vis_Laser_sem = std(data.visNopto.Vis_Laser,[],4)/sqrt(size(data.visNopto.Vis_Laser,4));

Vis_NoLaser_mean = mean(data.visNopto.Vis_NoLaser,4);
Vis_NoLaser_sem = std(data.visNopto.Vis_NoLaser,[],4)/sqrt(size(data.visNopto.Vis_NoLaser,4));

Vis_Laser_meanResp = squeeze(mean(Vis_Laser_mean(:,data.Info.timeStimInd,:),2));
Vis_Laser_meanResp_sem = squeeze(mean(Vis_Laser_sem(:,data.Info.timeStimInd,:),2));

Vis_NoLaser_meanResp = squeeze(mean(Vis_NoLaser_mean(:,data.Info.timeStimInd,:),2));
Vis_NoLaser_meanResp_sem = squeeze(mean(Vis_NoLaser_sem(:,data.Info.timeStimInd,:),2));

Vis_meanResp_diff = Vis_Laser_meanResp - Vis_NoLaser_meanResp;
Vis_meanResp_ratio = Vis_Laser_meanResp ./ Vis_NoLaser_meanResp;

timeVect = data.Info.timeVect;
xvect = data.Info.xvect;

directions = data.Info.directions; % directions = [0 45 90 135 180 225 270 315];

% Sort the data by prefered visual stimulus w/out laser
[bestResp, index] = sort(Vis_NoLaser_meanResp,2, 'descend');

%%
% % plot the ROIs with best visual resp
% for i  = 1:10
%     figure;
%     for ii = 1:4
%         subplot(4,1,ii); hold on
%         plot(Vis_Laser_mean(data.visNopto.Idx(i),:,ii)','r')
%         plot(Vis_Laser_mean(data.visNopto.Idx(i),:,ii)'+Vis_Laser_sem(data.visNopto.Idx(i),:,ii)','r:')
%         plot(Vis_Laser_mean(data.visNopto.Idx(i),:,ii)'-Vis_Laser_sem(data.visNopto.Idx(i),:,ii)','r:')
%
%         plot(Vis_NoLaser_mean(data.visNopto.Idx(i),:,ii)','k')
%         plot(Vis_NoLaser_mean(data.visNopto.Idx(i),:,ii)'+Vis_NoLaser_sem(data.visNopto.Idx(i),:,ii)','k:')
%         plot(Vis_NoLaser_mean(data.visNopto.Idx(i),:,ii)'-Vis_NoLaser_sem(data.visNopto.Idx(i),:,ii)','k:')
%     end
%     saveas(gcf,[saveDir,'\VisStim_ResponsiveMost_ROI',num2str(data.visNopto.Idx(i)) ,'.tif'])
% end

%%
% % plot best ROIs for effect of photostim on visual stim
% % control > photostim
% for k = 1:data.Info.nVisStimTypes
%     for i  = 1:5
%         figure;
%         for ii = 1:data.Info.nVisStimTypes
%             subplot(data.Info.nVisStimTypes,1,ii); hold on
%             plot(Vis_Laser_mean(data.visNopto.Idx_diff(i,k),:,ii)','r')
%             plot(Vis_Laser_mean(data.visNopto.Idx_diff(i,k),:,ii)'+Vis_Laser_sem(data.visNopto.Idx_diff(i,k),:,ii)','r:')
%             plot(Vis_Laser_mean(data.visNopto.Idx_diff(i,k),:,ii)'-Vis_Laser_sem(data.visNopto.Idx_diff(i,k),:,ii)','r:')
%
%             plot(Vis_NoLaser_mean(data.visNopto.Idx_diff(i,k),:,ii)','k')
%             plot(Vis_NoLaser_mean(data.visNopto.Idx_diff(i,k),:,ii)'+Vis_NoLaser_sem(data.visNopto.Idx_diff(i,k),:,ii)','k:')
%             plot(Vis_NoLaser_mean(data.visNopto.Idx_diff(i,k),:,ii)'-Vis_NoLaser_sem(data.visNopto.Idx_diff(i,k),:,ii)','k:')
%             if k == ii
%                 box on
%             end
%         end
%         sgtitle(['ROI#',num2str(data.visNopto.Idx_diff(i,k))])
%         saveas(gcf,[saveDir,'\VisStim_ResponsiveMost_ROI',num2str(data.visNopto.Idx_diff(i,k)) ,'.tif'])
%     end
% end

% % plot best ROIs for effect of photostim on visual stim, sort according
% % control < photostim
%
% for k = 1:data.Info.nVisStimTypes
%     for i  = 1:3
%         figure;
%         for ii = 1:data.Info.nVisStimTypes
%             subplot(data.Info.nVisStimTypes,1,ii); hold on
%             plot(Vis_Laser_mean(data.visNopto.Idx_diff(end-i+1,k),:,ii)','r')
%             plot(Vis_Laser_mean(data.visNopto.Idx_diff(end-i+1,k),:,ii)'+Vis_Laser_sem(data.visNopto.Idx_diff(end-i+1,k),:,ii)','r:')
%             plot(Vis_Laser_mean(data.visNopto.Idx_diff(end-i+1,k),:,ii)'-Vis_Laser_sem(data.visNopto.Idx_diff(end-i+1,k),:,ii)','r:')
%
%             plot(Vis_NoLaser_mean(data.visNopto.Idx_diff(end-i+1,k),:,ii)','k')
%             plot(Vis_NoLaser_mean(data.visNopto.Idx_diff(end-i+1,k),:,ii)'+Vis_NoLaser_sem(data.visNopto.Idx_diff(end-i+1,k),:,ii)','k:')
%             plot(Vis_NoLaser_mean(data.visNopto.Idx_diff(end-i+1,k),:,ii)'-Vis_NoLaser_sem(data.visNopto.Idx_diff(end-i+1,k),:,ii)','k:')
%             if k == ii
%                 box on
%             end
%         end
%         sgtitle(['ROI#',num2str(data.visNopto.Idx(end-i+1,k))])
%         %         saveas(gcf,[saveDir,'\VisStim_ResponsiveMost_ROI',num2str(data.visNopto.Idx(end-i+1,k)) ,'.tif'])
%
%         figure; hold on
%         h1 = plot(directions,Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),:),'ro-');
%         h2 = plot(directions,Vis_NoLaser_meanResp(data.visNopto.Idx_diff(end-i+1,k),:),'ko-');
%         xlabel('Directions')
%         xticks(directions)
%         ylabel('dF/F')
%         plot(directions,Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),:)+Vis_Laser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),:),'r:')
%         plot(directions,Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),:)-Vis_Laser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),:),'r:')
%         plot(directions,Vis_NoLaser_meanResp(data.visNopto.Idx_diff(end-i+1,k),:)+Vis_NoLaser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),:),'k:')
%         plot(directions,Vis_NoLaser_meanResp(data.visNopto.Idx_diff(end-i+1,k),:)-Vis_NoLaser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),:),'k:')
%         plot(directions(k),Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),k), 'ro',  'MarkerFaceColor',[1,0,0])
%
%         legend([h1 h2],{'Photostim','Control'})
%         sgtitle(['ROI#',num2str(data.visNopto.Idx(end-i+1,k))])
%     end
% end

%%
% % plot the ROIs with biggest effect of laser, trials sorted by response
% % magnitude in baseline


% for k =1:data.Info.nVisStimTypes
%     for i = 1:3
%         figure; hold on
%         h1 = plot(Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:)),'ro-');
%         h2 = plot(Vis_NoLaser_meanResp(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:)),'ko-');
%
%         plot(Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:))+Vis_Laser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,1),:)),'r:')
%         plot(Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:))-Vis_Laser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,1),:)),'r:')
%         plot(Vis_NoLaser_meanResp(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:))+Vis_NoLaser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,1),:)),'k:')
%         plot(Vis_NoLaser_meanResp(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:))-Vis_NoLaser_meanResp_sem(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,1),:)),'k:')
%
%         h3 = plot(Vis_meanResp_diff(data.visNopto.Idx_diff(end-i+1,k),index(data.visNopto.Idx_diff(end-i+1,k),:)),'bo-');
%
% %         % pry need to use "find" here, to plot a filled red dot where the
% %         % biggest impact of light is
% %         % something like find(directions == index
% %         plot(directions(index(data.visNopto.Idx_diff(end-i+1,k),k),Vis_Laser_meanResp(data.visNopto.Idx_diff(end-i+1,k),k), 'ro',  'MarkerFaceColor',[1,0,0])
%
%         legend([h1 h2 h3],{'Photostim','Control','Diff'})
%         sgtitle(['ROI#',num2str(data.visNopto.Idx(end-i+1,k))])
%         xticklabels([directions(index(data.visNopto.Idx_diff(end-i+1,k),:))])
%         xlabel('Directions, sorted')
%         ylabel('dF/F')
%         % h1.XticksLabel = (directions(index(data.visNopto.Idx_diff(end-i+1,1),:)))
%     end
% end
%
% figure;
% imagesc(Vis_meanResp_diff(index));
% plot(mean(Vis_meanResp_diff(index),1))

%% Select the ROIs with at least 1 visual stim where photostim had an effect

clear row col h_all p_all
time_ttest = timeVect>data.Info.xstim_vis(1) & timeVect<data.Info.xstim_vis(2)+response_lag;
VISnOPTO = squeeze(mean(data.visNopto.Vis_Laser(:,time_ttest,:,:),2));
VISonly = squeeze(mean(data.visNopto.Vis_NoLaser(:,time_ttest,:,:),2));
for i = 1:data.Info.nROIs
    for ii = 1:data.Info.nVisStimTypes
        [h,p] = ttest(squeeze(VISnOPTO(i,ii,:)),squeeze(VISonly(i,ii,:)),'Alpha',a);
        h_all(i,ii) = h;
        p_all(i,ii) = p;
    end
end
[row, col] =  find(h_all);
n_respROI = size(unique(row),1)/size(h_all,1);
Photostim_respROIs = unique(row);
disp([num2str(n_respROI*100,4),'% ROI whose visual response is modulate by opto, at least for one condition; alpha = ' num2str(a)])

% Vis_meanResp_diff_all = reshape(Vis_meanResp_diff,size(Vis_meanResp_diff,1)*data.Info.nVisStimTypes,1);
% figure; ecdf(Vis_meanResp_diff_all);

% plot them
if plot_indivROI
    for  i = 1:data.Info.nROIs
        if any(i == unique(row))
            ymax = max(max(max(Vis_Laser_mean(i,:,:)+Vis_Laser_sem(i,:,:))),max(max(Vis_NoLaser_mean(i,:,:)+Vis_NoLaser_sem(i,:,:))));
            ymin = min(min(min(Vis_Laser_mean(i,:,:)-Vis_Laser_sem(i,:,:))),min(min(Vis_NoLaser_mean(i,:,:)-Vis_NoLaser_sem(i,:,:))));
            figure;
            for ii = 1:data.Info.nVisStimTypes
                subplot(2,data.Info.nVisStimTypes/2,ii); hold on
                h1 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                    [Vis_Laser_mean(i,xvect,ii)+Vis_Laser_sem(i,xvect,ii) flip(Vis_Laser_mean(i,xvect,ii)-Vis_Laser_sem(i,xvect,ii))],[1 0 0],'LineStyle','none');
                set(h1,'facealpha',.5)
                h2 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                    [Vis_NoLaser_mean(i,xvect,ii)+Vis_NoLaser_sem(i,xvect,ii) flip(Vis_NoLaser_mean(i,xvect,ii)-Vis_NoLaser_sem(i,xvect,ii))],[0 0 0],'LineStyle','none');
                set(h2,'facealpha',.5)
                plot(timeVect(xvect),Vis_Laser_mean(i,xvect,ii),'r')
                plot(timeVect(xvect),Vis_NoLaser_mean(i,xvect,ii),'k')
                ylim([ymin-0.05 ymax+0.05])
                
                plot([data.Info.xstim_opto],[ymax ymax],'r','LineWidth',5)
                plot([data.Info.xstim_vis],[ymax+0.05 ymax+0.05],'k','LineWidth',5)
                %             plot([data.Info.xstim_opto(1) data.Info.xstim_opto(1)],[ymin ymax],'k:')
                %             plot([data.Info.xstim_opto(2) data.Info.xstim_opto(2)],[ymin ymax],'k:')
                
                if ii == data.Info.nVisStimTypes/2+1
                    xlabel('Time from photostim onset (s)');
                    ylabel('dF/F')
                end
                blah = find(i == row(:));
                if any(col(blah) == ii)
                    title(num2str(p_all(i,ii)));
                    box on
                end
                
            end
            sgtitle(['ROI#', num2str(i)]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.04, 0.80, 0.5]);
            
            if saveFigFlag
                saveas(gcf,[saveDir,'\',num2str(i),'_Traces.tif'])
                pause(0.2); close
            end
            
            figure;
            subplot(2,1,1); hold on
            plot(directions,Vis_Laser_meanResp(i,:),'ro-')
            plot(directions,Vis_NoLaser_meanResp(i,:),'ko-')
            plot(directions, Vis_meanResp_diff(i,:),'bo-')
            plot(directions,Vis_Laser_meanResp(i,:)+Vis_Laser_meanResp_sem(i,:),'r:')
            plot(directions,Vis_Laser_meanResp(i,:)-Vis_Laser_meanResp_sem(i,:),'r:')
            plot(directions,Vis_NoLaser_meanResp(i,:)+Vis_NoLaser_meanResp_sem(i,:),'k:')
            plot(directions,Vis_NoLaser_meanResp(i,:)-Vis_NoLaser_meanResp_sem(i,:),'k:')
            %         xlabel('Directions')
            xticks(directions)
            ylabel('dF/F')
            title('Tuning, absolute directions')
            
            subplot(2,1,2); hold on
            plot(Vis_Laser_meanResp(i,index(i,:)),'ro-')
            plot(Vis_NoLaser_meanResp(i,index(i,:)),'ko-')
            plot(Vis_meanResp_diff(i,index(i,:)),'bo-')
            plot(Vis_Laser_meanResp(i,index(i,:))+Vis_Laser_meanResp_sem(i,index(i,:)),'r:')
            plot(Vis_Laser_meanResp(i,index(i,:))-Vis_Laser_meanResp_sem(i,index(i,:)),'r:')
            plot(Vis_NoLaser_meanResp(i,index(i,:))+Vis_NoLaser_meanResp_sem(i,index(i,:)),'k:')
            plot(Vis_NoLaser_meanResp(i,index(i,:))-Vis_NoLaser_meanResp_sem(i,index(i,:)),'k:')
            xlabel('Direction #')
            xticklabels(directions(index(i,:)))
            ylabel('dF/F')
            title('Tuning, absolute directions')
            title('Tuning, sorted')
            
            sgtitle(['ROI#', num2str(i)]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.85, 0.04, 0.1, 0.5]);
            %         text(-6.5,0.55,[ 'exp vs light ctrl ' num2str(p,2)])
            if saveIndivROI_Flag
                saveas(gcf,[saveDir,'\',num2str(i),'_Directions_Tuning.tif'])
                pause(0.2); close
            end
        end
        
    end
end
%% Select the ROIs with at least 1 visual stim eliciting significant response

clear row col
% time_ttest = timeVect>data.Info.xstim_vis(1) & timeVect<data.Info.xstim_vis(2)+response_lag;
numberOfDataPoints = size(find(time_ttest),2);
VISonly = squeeze(mean(data.visNopto.Vis_NoLaser(:,time_ttest,:,:),2));
VISonly_base = squeeze(mean(data.visNopto.Vis_NoLaser(:,data.Info.timeBaseInd(1:numberOfDataPoints),:,:),2));

for i = 1:data.Info.nROIs
    for ii = 1:data.Info.nVisStimTypes
        [h,p] = ttest(squeeze(VISonly(i,ii,:)),squeeze(VISonly_base(i,ii,:)),'Alpha',a);
        h_all_Vis(i,ii) = h;
        p_all_Vis(i,ii) = p;
    end
end
[row, col] =  find(h_all_Vis);
respROIs = unique(row);
n_respROI = size(unique(row),1)/data.Info.nROIs;
disp([num2str(n_respROI*100,4),'% ROI with visual response, at least for one condition'])

for i = 1:data.Info.nROIs
    k(i) = size(find(h_all_Vis(i,:)),2);
end
disp(['ROIs significantly responded in average to ' num2str(mean(k),2),' stimuli ; aplha =' num2str(a)])

%% plot the summary of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1) Real directions (not sorted by elicited activity)
Pop_avgDiff = mean(Vis_meanResp_diff,1);
Pop_semDiff = std(Vis_meanResp_diff,[],1)./sqrt(size(Vis_meanResp_diff,1));

Pop_avgDiv = mean(Vis_meanResp_ratio,1);
Pop_semDiv = std(Vis_meanResp_ratio,[],1)./sqrt(size(Vis_meanResp_ratio,1));

% 1.1) All ROIs, all directions
AllROIs_Diff_DirNotSorted = figure;
subplot(1,4,1:3); hold on;
b1 = bar(directions,Pop_avgDiff);
for i = 1:size(directions,2)
    plot([directions(i) directions(i)],[Pop_avgDiff(i)-Pop_semDiff(i) Pop_avgDiff(i)+Pop_semDiff(i)], 'k' )
end
title('All ROIs, modulation by photostim')
ylabel({'Opto+Vis - Vis Stim','dF/F'})
xticks(directions)
xlabel('Grating Direction')
sem_max = max(max(Pop_semDiff),abs(min(Pop_semDiff)));
ymax = max(b1.YData); ymin = min(b1.YData);
ylimit = max(ymax+sem_max, abs(ymin)+sem_max);
ylim([-ylimit-1/10*ylimit ylimit+1/10*ylimit])
yl = ylim;
subplot(1,4,4); hold on
b2 = bar(mean(Pop_avgDiff));
plot([1 1],[mean(Pop_avgDiff)-std(Pop_avgDiff)./sqrt(size(Pop_avgDiff,2)) mean(Pop_avgDiff)+std(Pop_avgDiff)./sqrt(size(Pop_avgDiff,2))],'k');
b2(1).BaseValue = 0;
ylim(yl)
title('Average across directions')
clear sem_max  ymax ymin ylimit

% % Does not make sense to do that for all ROIs, as a lot of visResp are gonna
% % be 0
% AllROIs_Div_DirNotSorted = figure;
% subplot(1,4,1:3);hold on;
% b1 = bar(directions,Pop_avgDiv);
% for i = 1:size(directions,2)
%     plot([directions(i) directions(i)],[Pop_avgDiv(i)-Pop_semDiv(i) Pop_avgDiv(i)+Pop_semDiv(i)], 'k' )
% end
% b1(1).BaseValue = 1;
% title('All ROIs, modulation by photostim')
% ylabel({'Opto+Vis / Vis Stim','dF/F'})
% xticks(directions)
% xlabel('Grating Direction')
% sem_max = max(max(Pop_semDiv),abs(min(Pop_semDiv)));
% ymax = max(b1.YData); ymin = min(b1.YData);
% ylimit = max(ymax+sem_max, abs(ymin)+sem_max);
% ylim([-ylimit-1/10*ylimit ylimit+1/10*ylimit])
% yl = ylim;
% subplot(1,4,4); hold on
% b2 = bar(mean(Pop_avgDiv));
% plot([1 1],[mean(Pop_avgDiv)-std(Pop_avgDiv)./sqrt(size(Pop_avgDiv,2)) mean(Pop_avgDiv)+std(Pop_avgDiv)./sqrt(size(Pop_avgDiv,2))],'k');
% b2(1).BaseValue = 1;
% ylim(yl)
% title('Average across directions')
% clear sem_max  ymax ymin ylimit
% % axis off

% figure; hold on;
% bar(directions,Pop_avgDiv)
% for i = 1:size(directions,2)
%     plot([directions(i) directions(i)],[Pop_avgDiv(i)-Pop_semDiv(i) Pop_avgDiv(i)+Pop_semDiv(i)], 'k' )
% end
% title('Population average modulation by light')
% ylabel('Opto/Vis Stim')
% xticks(directions)
% xlabel('Grating Direction')


% 1.2) Responsive ROIs only, all of their directions

% h_all = logical(h_all);
% ResponsiveROI_diff = Vis_meanResp_diff(h_all);
% figure;ecdf(ResponsiveROI_diff);
% hold on
% plot([0 0],[0 1],'k:');
% plot([-1.5 0],[0.5 0.5],'k:');
% title({'Selection on OPTO effect','"Opto+visual" - "visual"'})
% 
% figure;ecdf(ResponsiveROI_div);
% ResponsiveROI_div = Vis_meanResp_ratio(h_all);
% title({'Selection on OPTO effect','"Opto+visual" / "visual"'})

h_all_Vis = logical(h_all_Vis);
ResponsiveROI_diff = Vis_meanResp_diff(h_all_Vis);
figure;ecdf(ResponsiveROI_diff);
hold on
plot([0 0],[0 1],'k:');
plot([-1.5 0],[0.5 0.5],'k:');
title({'Selection on VIS responses','"Opto+visual" - "visual"'})

ResponsiveROI_div = Vis_meanResp_ratio(h_all_Vis);
figure;ecdf(ResponsiveROI_div);
title({'Selection on VIS responses','"Opto+visual" / "visual"'})

% 1.3) Only {ROI,direction} eliciting significant response...
% % Attention! here i am selecting by the effect of opto, not visual stim
ResponsiveCombo = cell(1,size(directions,2));
ResponsiveCombo_div = cell(1,size(directions,2));
counter = 1;
% for i = 1:data.Info.nROIs
%     %     if any(i == unique(row))
%     if any(h_all(i,:))
%         idx = find(h_all(i,:));
%         for ii = 1:length(idx)
%             ResponsiveCombo{idx(ii)}(counter) = Vis_meanResp_diff(i,idx(ii));
%             ResponsiveCombo_div{idx(ii)}(counter) = Vis_meanResp_ratio(i,idx(ii));
%         end
%         counter=counter+1;
%         clear idx
%     end
% end

for i = 1:data.Info.nROIs
    %     if any(i == unique(row))
    if any(h_all_Vis(i,:))
        idx = find(h_all_Vis(i,:));
        for ii = 1:length(idx)
            ResponsiveCombo{idx(ii)}(counter) = Vis_meanResp_diff(i,idx(ii));
            ResponsiveCombo_div{idx(ii)}(counter) = Vis_meanResp_ratio(i,idx(ii));
        end
        counter=counter+1;
        clear idx
    end
end

for i = 1:size(directions,2)
    idx = find(ResponsiveCombo{i});
    PopAvg(i) = mean(ResponsiveCombo{i}(idx));
    PopSem(i) = std(ResponsiveCombo{i}(idx))./sqrt(size(idx,2));
    clear idx
    
    idx_pos = find(ResponsiveCombo{i}>0);
    PosResp_avg(i) = mean(ResponsiveCombo{i}(idx_pos));
    PosResp_sem(i)= std(ResponsiveCombo{i}(idx_pos))./sqrt(size(idx_pos,2));
    nPosModROIs(1,i) = size(idx_pos,2);
    
    idx_neg = find(ResponsiveCombo{i}<0);
    NegResp_avg(i) = mean(ResponsiveCombo{i}(idx_neg));
    NegResp_sem(i) = std(ResponsiveCombo{i}(idx_neg))./sqrt(size(idx_neg,2));
    nNegModROIs(1,i) = size(idx_neg,2);

    clear idx_pos ; clear idx_neg
     
    idx = find(ResponsiveCombo_div{i});
    PopAvg_div(i) = mean(ResponsiveCombo_div{i}(idx));
    PopSem_div(i) = std(ResponsiveCombo_div{i}(idx))./sqrt(size(idx,2));
    
    idx_pos = find(ResponsiveCombo_div{i}>1);
    PosResp_avg_div(i) = mean(ResponsiveCombo_div{i}(idx_pos));
    PosResp_sem_div(i)= std(ResponsiveCombo_div{i}(idx_pos))./sqrt(size(idx_pos,2));
    nPosModROIs(2,i) = size(idx_pos,2);

    idx_neg = find(ResponsiveCombo_div{i}<1 & ResponsiveCombo_div{i}~=0);
    NegResp_avg_div(i) = mean(ResponsiveCombo_div{i}(idx_neg));
    NegResp_sem_div(i) = std(ResponsiveCombo_div{i}(idx_neg))./sqrt(size(idx_neg,2));
    nNegModROIs(2,i) = size(idx_neg,2);

    clear idx ; clear idx_pos ; clear idx_neg
end
Resp_sorted = cat(1,PosResp_avg, abs(NegResp_avg));
Resp_sorted_div = cat(1,PosResp_avg_div, abs(NegResp_avg_div));

% ...plot the subtraction
ResponsiveROIs_Diff_allDir_NegPosSep_DirNotSorted = figure;
subplot(1,4,1:3); hold on;
bar(directions,PosResp_avg)
bar(directions,NegResp_avg)
for i = 1:size(directions,2)
    plot([directions(i) directions(i)],[PosResp_avg(i)-PosResp_sem(i) PosResp_avg(i)+PosResp_sem(i)], 'k' )
    plot([directions(i) directions(i)],[NegResp_avg(i)-NegResp_sem(i) NegResp_avg(i)+NegResp_sem(i)], 'k' )
end
title('Responsive population average modulation by light')
ylabel('Opto-Vis Stim')
xticks(directions)
xlabel('Grating Direction')
legend({'Positive','Negative'})
yl = ylim;
ylimit = max(abs(yl(1)), yl(2));
ylim([-ylimit ylimit]);
for i = 1:data.Info.nVisStimTypes
    text(directions(i),yl(2),num2str(nPosModROIs(1,i)),'HorizontalAlignment','center')
    text(directions(i),yl(1),num2str(nNegModROIs(1,i)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom')
end
subplot(1,4,4); hold on
bar(mean(PosResp_avg));
bar(mean(NegResp_avg));
plot([1 1],[mean(PosResp_avg)-std(PosResp_avg)./sqrt(size(PosResp_avg,1)) mean(PosResp_avg)+std(PosResp_avg)./sqrt(size(PosResp_avg,1))],'k')
plot([1 1],[mean(NegResp_avg)-std(NegResp_avg)./sqrt(size(NegResp_avg,1)) mean(NegResp_avg)+std(NegResp_avg)./sqrt(size(NegResp_avg,1))],'k')
ylim([-ylimit ylimit]);
title('Average across Dir')
text(1,yl(2),num2str(sum(nPosModROIs(1,:))),...
    'HorizontalAlignment','center','Color','b')
text(1,yl(1),num2str(sum(nNegModROIs(1,:))),...
    'HorizontalAlignment','center','Color','r','VerticalAlignment','bottom')

ResponsiveROIs_Diff_allDir_DirNotSorted = figure;
subplot(1,4,1:3); hold on;
bar(directions,PopAvg)
for i = 1:size(directions,2)
    plot([directions(i) directions(i)],[PopAvg(i)-PopSem(i) PopAvg(i)+PopSem(i)], 'k' )
end
title('Responsive population average modulation by light')
ylabel('Opto-Vis Stim')
xticks(directions)
xlabel('Grating Direction')
yl = ylim;
ylimit = max(abs(yl(1)), yl(2));
ylim([-ylimit ylimit]);
subplot(1,4,4); hold on
bar(mean(PopAvg));
plot([1 1],[mean(PopAvg)-std(PopAvg)./sqrt(size(PopAvg,1)) mean(PopAvg)+std(PopAvg)./sqrt(size(PopAvg,1))],'k')
ylim([-ylimit ylimit]);
title('Average across Dir')

ResponsiveROIs_Diff_allDir_NegPosSepAbs_DirNotSorted = figure;
subplot(1,4,1:3); hold on;
b1 = bar(directions,Resp_sorted');
for i = 1:size(directions,2)
    plot([directions(i)-7 directions(i)-7],[PosResp_avg(i)-PosResp_sem(i) PosResp_avg(i)+PosResp_sem(i)], 'k' )
    plot([directions(i)+7 directions(i)+7],[abs(NegResp_avg(i))-NegResp_sem(i) abs(NegResp_avg(i))+NegResp_sem(i)], 'k' )
end
title('Responsive population average modulation by light')
ylabel('Opto-Vis Stim')
xticks(directions)
xlabel('Grating Direction')
legend({'Positive','Negative'})
yl = ylim;
ylimit = max(abs(yl(1)), yl(2));
ylim([-ylimit ylimit]);
subplot(1,4,4); hold on
b2 = bar(mean(Resp_sorted,2));
b2.FaceColor = 'flat';
b2.CData(2,:) = [0.8500    0.3250    0.0980];
plot([1 1],[mean(PosResp_avg)-std(PosResp_avg)./sqrt(size(PosResp_avg,1)) mean(PosResp_avg)+std(PosResp_avg)./sqrt(size(PosResp_avg,1))],'k')
plot([2 2],[mean(abs(NegResp_avg))-std(NegResp_avg)./sqrt(size(NegResp_avg,1)) mean(abs(NegResp_avg))+std(NegResp_avg)./sqrt(size(NegResp_avg,1))],'k')
ylim([-ylimit ylimit]);
title('Average across Dir')

% ...plot the division
ResponsiveROIs_Div_allDir_NegPosSep_DirNotSorted = figure;
subplot(1,4,1:3); hold on;
b1 = bar(directions,PosResp_avg_div);
b2 = bar(directions,NegResp_avg_div);
for i = 1:size(directions,2)
    plot([directions(i) directions(i)],[PosResp_avg_div(i)-PosResp_sem_div(i) PosResp_avg_div(i)+PosResp_sem_div(i)], 'k' )
    plot([directions(i) directions(i)],[NegResp_avg_div(i)-NegResp_sem_div(i) NegResp_avg_div(i)+NegResp_sem_div(i)], 'k' )
end
title('Responsive population average modulation by light')
ylabel('"Opto+Vis"/"Vis" Stim')
xticks(directions)
xlabel('Grating Direction')
b1.BaseValue = 1;
b2.BaseValue = 1;
lgd = legend({'Increase','Decrease'});
title(lgd,'Effect on visual response');
yl = ylim;
ylimit = max(1-yl(1), yl(2)-1);
ylim([1-ylimit ylimit+1]);
for i = 1:data.Info.nVisStimTypes
    text(directions(i),yl(2),num2str(nPosModROIs(2,i)),...
        'HorizontalAlignment','center','Color','b')
    text(directions(i),yl(1),num2str(nNegModROIs(2,i)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','Color','r')
end
subplot(1,4,4); hold on
b3 = bar(mean(PosResp_avg_div));
b4 = bar(mean(NegResp_avg_div));
plot([1 1],[mean(PosResp_avg_div)-std(PosResp_avg_div)./sqrt(size(PosResp_avg_div,1)) mean(PosResp_avg_div)+std(PosResp_avg_div)./sqrt(size(PosResp_avg_div,1))],'k')
plot([1 1],[mean(NegResp_avg_div)-std(NegResp_avg_div)./sqrt(size(NegResp_avg_div,1)) mean(NegResp_avg_div)+std(NegResp_avg_div)./sqrt(size(NegResp_avg_div,1))],'k')
b3.BaseValue = 1;
b4.BaseValue = 1;
ylim([1-ylimit ylimit+1]);
title('Average across Dir')
text(1,yl(2),num2str(sum(nPosModROIs(2,:))),...
    'HorizontalAlignment','center','Color','b')
text(1,yl(1),num2str(sum(nNegModROIs(2,:))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Color','r')

ResponsiveROIs_Div_allDir_DirNotSorted = figure;
subplot(1,4,1:3); hold on;hold on;
b1 = bar(directions,PopAvg_div);
for i = 1:size(directions,2)
    plot([directions(i) directions(i)],[PopAvg_div(i)-PopSem_div(i) PopAvg_div(i)+PopSem_div(i)], 'k' )
end
title('Responsive population average modulation by light')
ylabel('"Opto+Vis"/"Vis" Stim')
xticks(directions)
xlabel('Grating Direction')
b1.BaseValue = 1;
yl = ylim;
ylimit = max(1-yl(1), yl(2)-1);
ylim([1-ylimit ylimit+1]);
subplot(1,4,4); hold on
b2 = bar(mean(PopAvg_div));
plot([1 1],[mean(PopAvg_div)-std(PopAvg_div)./sqrt(size(PopAvg_div,1)) mean(PopAvg_div)+std(PopAvg)./sqrt(size(PopAvg_div,1))],'k')
b2.BaseValue = 1;
ylim([1-ylimit ylimit+1]);
title('Average across Dir')
clear b1 b2

% ResponsiveROIs_Div_allDir_NegPosSepAbs_DirNotSorted = figure;
% subplot(1,4,1:3); hold on;hold on;
% b1 = bar(directions,Resp_sorted_div');
% for i = 1:size(directions,2)
%     plot([directions(i)-7 directions(i)-7],[PosResp_avg_div(i)-PosResp_sem_div(i) PosResp_avg_div(i)+PosResp_sem_div(i)], 'k' )
%     plot([directions(i)+7 directions(i)+7],[2-NegResp_avg_div(i)-NegResp_sem_div(i) 2-NegResp_avg_div(i)+NegResp_sem_div(i)], 'k' )
% end
% title('Responsive population average modulation by light')
% ylabel('"Opto+Vis"/"Vis" Stim')
% xticks(directions)
% xlabel('Grating Direction')
% legend({'Positive','Negative'})
% b1(1).BaseValue = 1;
% yl = ylim;
% ylimit = max(1-yl(1), yl(2)-1);
% ylim([1-ylimit ylimit+1]);
% subplot(1,4,4); hold on
% b2 = bar(mean(Resp_sorted_div,2));
% b2.FaceColor = 'flat';
% b2.CData(2,:) = [0.8500    0.3250    0.0980];
% plot([1 1],[mean(PosResp_avg_div)-std(PosResp_avg_div)./sqrt(size(PosResp_avg_div,1)) mean(PosResp_avg_div)+std(PosResp_avg_div)./sqrt(size(PosResp_avg_div,1))],'k')
% plot([2 2],[2-mean(NegResp_avg_div)-std(NegResp_avg_div)./sqrt(size(NegResp_avg_div,1)) 2-mean(abs(NegResp_avg_div))+std(NegResp_avg_div)./sqrt(size(NegResp_avg_div,1))],'k')
% b2(1).BaseValue = 1;
% ylim([1-ylimit ylimit+1]);
% title('Average across Dir')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2) Sorted by best directions
% 2.1) All ROIs
for i = 1:data.Info.nROIs
    Diff_best(i,:) = Vis_meanResp_diff(i,index(i,:));
    %     Diff_best_sem(i,:) = std(Vis_meanResp_diff(i,index(i,:));
    Ratio_best(i,:) = Vis_meanResp_ratio(i,index(i,:));
end

AllROI_Diff_DirSorted = figure('Units', 'Normalized', 'OuterPosition',[0 0.5 0.33 0.5]);
hold on
plot(mean(Diff_best,1),'r');
plot(mean(Diff_best,1)+(std(Diff_best,[],1)./sqrt(size(Diff_best,1))),'r:')
plot(mean(Diff_best,1)-(std(Diff_best,[],1)./sqrt(size(Diff_best,1))),'r:')
plot([1 size(directions,2)], [0 0],'k:')
xlabel('Preferred direction #')
ylabel({'Opto-Ctrl','dF/F'})
title('All ROIs, sorted by preferred direction')

% % should not plot the division when not selecting the ROIs
% AllROI_Div_DirSorted = figure('Units', 'Normalized', 'OuterPosition',[0.33 0.5 0.33 0.5]);
% hold on
% plot(mean(Ratio_best,1),'r');
% plot(mean(Ratio_best,1)+(std(Ratio_best,[],1)./sqrt(size(Ratio_best,1))),'r:')
% plot(mean(Ratio_best,1)-(std(Ratio_best,[],1)./sqrt(size(Ratio_best,1))),'r:')
% plot([1 size(directions,2)], [1 1],'k:')
% xlabel('Best direction #')
% ylabel({'Opto/Ctrl','dF/F'})
% title('All ROIs, sorted by best direction')

% 2.2) Responsive ROIs only
ResponsiveROI_diff_best = NaN(data.Info.nROIs,data.Info.nVisStimTypes);
ResponsiveROI_div_best = NaN(data.Info.nROIs,data.Info.nVisStimTypes);
for i = 1:data.Info.nROIs
    h_all_Vis2(i,:) =  h_all_Vis(i,index(i,:));
    for ii = 1:data.Info.nVisStimTypes
        if h_all_Vis2(i,ii)
            ResponsiveROI_diff_best(i,ii) = Vis_meanResp_diff(i,ii);
            ResponsiveROI_div_best(i,ii) = Vis_meanResp_ratio(i,ii);
        end
    end
end

for ii = 1:data.Info.nVisStimTypes
    nRespROIsPerVisStim(ii) = size(find(~isnan(ResponsiveROI_diff_best(:,ii))),1);
end

Pop_ResponsiveROI_diff_best = nanmean(ResponsiveROI_diff_best,1);
Pop_ResponsiveROI_div_best = nanmean(ResponsiveROI_div_best,1);

RespROI_Diff_BoxPlot = figure('Units', 'Normalized', 'OuterPosition',[0 0.5 0.33 0.5]);
boxplot(ResponsiveROI_diff_best,'Notch','on','OutlierSize',1);
ylim([-0.5 0.5])
hold on; plot([1 8],[0 0],'k:')
xlabel('Preferrred Visual Stim #')
ylabel('Opto+Vis - Vis')
title('Impact of opto, sign vis responses')
for i = 1:data.Info.nVisStimTypes
    text(i,0.4,num2str(nRespROIsPerVisStim(i)),'HorizontalAlignment','center')
end

RespROI_Div_BoxPlot = figure('Units', 'Normalized', 'OuterPosition',[0.33 0.5 0.33 0.5]);
b = boxplot(ResponsiveROI_div_best,'Notch','on','OutlierSize',1);
ylim([-6 6])
hold on; plot([1 8],[1 1],'k:')
xlabel('Preferrred Visual Stim #')
ylabel('Opto+Vis / Vis')
title('Impact of opto, sign vis responses')
for i = 1:data.Info.nVisStimTypes
    text(i,4,num2str(nRespROIsPerVisStim(i)),'HorizontalAlignment','center')
end


RespROI_Div_HM = figure('Position',[100 100 500 900]);
climit = max(max(max(ResponsiveROI_div_best)),min(min(ResponsiveROI_div_best)));
% imagesc(ResponsiveROI_div_best,[-climit climit]);
imagesc(ResponsiveROI_div_best,'AlphaData', ~isnan(ResponsiveROI_div_best))
caxis([-climit climit])
colormap('diverging_bwr_40_95_c42_n256')
c = colorbar;
c.Label.String = 'dFF';
xlabel('Preferrred Visual Stim #')
xticks(1:data.Info.nVisStimTypes)
ylabel('ROI #')
title('Response ratio, Opto+Vis / Vis')
set(gca,'Color',[.5 .5 .5])

RespROI_Diff_HM = figure('Position',[600 100 500 900]);
climit = max(max(max(ResponsiveROI_diff_best)),min(min(ResponsiveROI_diff_best)));
imagesc(ResponsiveROI_diff_best,'AlphaData', ~isnan(ResponsiveROI_diff_best))
caxis([-climit climit])
colormap('diverging_bwr_40_95_c42_n256')
c = colorbar;
c.Label.String = 'dFF';
% c.Label.VerticalAlignment = 'middle';
% c.Position = [0.92 0.11 0.03 0.75]; %[x y width height]
title('Response diff, Opto+Vis - Vis')
xticks(1:data.Info.nVisStimTypes)
xlabel('Preferred Stimulus #')
ylabel('ROI #')
set(gca,'Color',[.5 .5 .5])
% [r,c] = find(isnan(ResponsiveROI_diff_best)) ;
% hold on ; 
% plot(c, r, '+','MarkerSize',2) ;

%% save figures
if saveFigFlag
    if ~ exist([mainDir '\Figures\PopPlots'])
        mkdir([mainDir '\Figures\PopPlots'])
    end
    figure(AllROIs_Diff_DirNotSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\AllROI_Diff_DirNotSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\AllROI_Diff_DirNotSorted']);
%     figure(AllROIs_Div_DirNotSorted)
%     saveas(gcf,[mainDir '\Figures\PopPlots\AllROI_Div_DirNotSorted'],'fig');
%     screen2png([mainDir '\Figures\PopPlots\AllROI_Div_DirNotSorted']);
    
    figure(ResponsiveROIs_Diff_allDir_NegPosSep_DirNotSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\ResponsiveROIs_Diff_NegPosSep_DirNotSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\ResponsiveROIs_Diff_NegPosSep_DirNotSorted']);
    figure(ResponsiveROIs_Diff_allDir_DirNotSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\ResponsiveROIs_Diff_DirNotSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\ResponsiveROIs_Diff_DirNotSorted']);
    figure(ResponsiveROIs_Diff_allDir_NegPosSepAbs_DirNotSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\ResponsiveROIs_Diff_NegPosSepAbs_DirNotSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\ResponsiveROIs_Diff_NegPosSepAbs_DirNotSorted']);
    
    figure(ResponsiveROIs_Div_allDir_NegPosSep_DirNotSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\ResponsiveROIs_Div_NegPosSep_DirNotSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\ResponsiveROIs_Div_NegPosSep_DirNotSorted']);
    figure(ResponsiveROIs_Div_allDir_DirNotSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\ResponsiveROIs_Div_DirNotSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\ResponsiveROIs_Div_DirNotSorted']);
%     figure(ResponsiveROIs_Div_allDir_NegPosSepAbs_DirNotSorted)
%     saveas(gcf,[mainDir '\Figures\PopPlots\ResponsiveROIs_Div_NegPosSepAbs_DirNotSorted'],'fig');
%     screen2png([mainDir '\Figures\PopPlots\ResponsiveROIs_Div_NegPosSepAbs_DirNotSorted']);
%     
    figure(AllROI_Diff_DirSorted)
    saveas(gcf,[mainDir '\Figures\PopPlots\AllROI_Diff_DirSorted'],'fig');
    screen2png([mainDir '\Figures\PopPlots\AllROI_Diff_DirSorted']);
%     figure(AllROI_Div_DirSorted)
%     saveas(gcf,[mainDir '\Figures\PopPlots\AllROI_Div_DirSorted'],'fig');
%     screen2png([mainDir '\Figures\PopPlots\AllROI_Div_DirSorted']);
    
    figure(RespROI_Diff_BoxPlot)
    saveas(gcf,[mainDir '\Figures\PopPlots\RespROI_Diff_DirSorted_BoxPlot'],'fig');
    screen2png([mainDir '\Figures\PopPlots\RespROI_Diff_DirSorted_BoxPlot']);
    figure(RespROI_Div_BoxPlot)
    saveas(gcf,[mainDir '\Figures\PopPlots\RespROI_Div_DirSorted_BoxPlot'],'fig');
    screen2png([mainDir '\Figures\PopPlots\RespROI_Div_DirSorted_BoxPlot']);
    
    figure(RespROI_Diff_HM)
    saveas(gcf,[mainDir '\Figures\PopPlots\RespROI_Diff_DirSorted_HM'],'fig');
    screen2png([mainDir '\Figures\PopPlots\RespROI_Diff_DirSorted_HM']);
    figure(RespROI_Div_HM)
    saveas(gcf,[mainDir '\Figures\PopPlots\RespROI_Div_DirSorted_HM'],'fig');
    screen2png([mainDir '\Figures\PopPlots\RespROI_Div_DirSorted_HM']);
end

if saveVarFlag
    data.VisResponsiveROIs.List = respROIs;
    data.VisResponsiveROIs.h = h_all_Vis;
    data.VisResponsiveROIs.p = p_all_Vis;
    
    data.PhotostimResponsiveROIs.List = Photostim_respROIs;
    data.PhotostimResponsiveROIs.h = h_all;
    data.PhotostimResponsiveROIs.p = p_all;
    try
        save([mainDir, '\analysis\', data.Info.saveName],'data')
    catch
        save([mainDir, '\analysis\data_CMloop6_day8_L2_popDone'],'data')
    end
end

%% Effect of opto stim on baseline activity
if data.Info.NoVisStim_exist
    clear row col h_OptoOnly_all p_OptoOnly_all
    
    NoVis_Laser_mean = mean(data.OptoOnly.Laser_NoVis,3);
    NoVis_Laser_sem = std(data.OptoOnly.Laser_NoVis,[],3)/sqrt(size(data.OptoOnly.Laser_NoVis,3));
    
    NoVis_NoLaser_mean = mean(data.OptoOnly.NoLaser_NoVis,3);
    NoVis_NoLaser_sem = std(data.OptoOnly.NoLaser_NoVis,[],3)/sqrt(size(data.OptoOnly.NoLaser_NoVis,3));
    
    OPTO = squeeze(mean(data.OptoOnly.Laser_NoVis(:,time_ttest,:),2));
    SPONTANEOUS = squeeze(mean(data.OptoOnly.NoLaser_NoVis(:,time_ttest,:),2));
    for i = 1:size(OPTO,1)
        [h,p] = ttest(squeeze(OPTO(i,:)),squeeze(SPONTANEOUS(i,:)),a);
        h_OptoOnly_all(i) = h;
        p_OptoOnly_all(i) = p;
    end
    [row] =  find(h_OptoOnly_all);
    n_respROI = size(row,2)/size(h_OptoOnly_all,2);
    disp([num2str(n_respROI*100,4),'% ROI respond to optogenetic'])
    
    if plot_indivROI      
        for  i = 1:size(row,2)
            figure; hold on
            h1 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                [NoVis_Laser_mean(row(i),xvect)+NoVis_Laser_sem(row(i),xvect) flip(NoVis_Laser_mean(row(i),xvect)-NoVis_Laser_sem(row(i),xvect))],[1 0 0],'LineStyle','none');
            set(h1,'facealpha',.5)
            h2 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                [NoVis_NoLaser_mean(row(i),xvect)+NoVis_NoLaser_sem(row(i),xvect) flip(NoVis_NoLaser_mean(row(i),xvect)-NoVis_NoLaser_sem(row(i),xvect))],[0 0 0],'LineStyle','none');
            set(h2,'facealpha',.5)
            plot(timeVect(xvect),NoVis_Laser_mean(row(i),xvect),'r')
            plot(timeVect(xvect),NoVis_NoLaser_mean(row(i),xvect),'k')
            yl = ylim;
            xlabel('Time for phosotim onset (s)')
            ylabel('dF/F')
            plot([0 data.Info.visStimLength],[yl(2) yl(2)],'Color',[1 0 0],'LineWidth',8)
            title({['Effect of Opto Only'], ['ROI # ', num2str(row(i)) ,', p=' ,num2str(p_OptoOnly_all(row(i)),2)]});
            if saveIndivROI_Flag
                saveas(gcf,[saveDir,'\StimOnBase_ResponsiveROI',num2str(row(i)) ,'.tif'])
                pause(0.2)
                close
            end
        end
    end

NoVis_meanResp_diff = NoVis_Laser_mean - NoVis_NoLaser_mean;
NoVis_meanResp_ratio = NoVis_Laser_mean ./ NoVis_NoLaser_mean;
NoVis_meanResp_index = (NoVis_Laser_mean - NoVis_NoLaser_mean)./(NoVis_Laser_mean + NoVis_NoLaser_mean);

counter = 1;
% % could use a logical here
for i = 1:data.Info.nROIs
    %     if any(i == unique(row))
    if h_OptoOnly_all(i) == 1
        Responsive_OptoOnSpont(counter) = mean(NoVis_meanResp_diff(i,time_ttest),2);
        if Responsive_OptoOnSpont(counter)>0
            ResponsivePos_OptoOnSpont(counter) = true;
        end
        counter=counter+1;
    end
end
OptoOnSpont = figure; hold on
bar(1,mean(Responsive_OptoOnSpont));
plot([1 1],[mean(Responsive_OptoOnSpont)-std(Responsive_OptoOnSpont)./sqrt(size(Responsive_OptoOnSpont,2)) mean(Responsive_OptoOnSpont)+std(Responsive_OptoOnSpont)./sqrt(size(Responsive_OptoOnSpont,2))],'k')
yl=ylim;
text(1,yl(2),num2str(size(Responsive_OptoOnSpont,2)),'HorizontalAlignment','center')
ylabel('dF/F, Opto - Spont')
xticks([1 2])
xticklabels({'Avg','Segregated'})
if exist('ResponsivePos_OptoOnSpont','var')
    bar(2,mean(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont)));
    plot([2 2],[mean(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont))-std(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont))./sqrt(size(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont),2)) mean(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont))+std(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont))./sqrt(size(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont),2))],'k')
    bar(2,mean(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont)));
    plot([2 2],[mean(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont))-std(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont))./sqrt(size(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont),2)) mean(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont))+std(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont))./sqrt(size(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont),2))],'k')
text(2,yl(2),num2str(size(Responsive_OptoOnSpont(ResponsivePos_OptoOnSpont),2)),'HorizontalAlignment','center')
text(2,yl(1),num2str(size(Responsive_OptoOnSpont(~ResponsivePos_OptoOnSpont),2)),...
    'HorizontalAlignment','center','VerticalAlignment','bottom')
end
% legend({'Avg','Positive','Negative'})
title('Opto Stim on Spont activity, resp. ROIs only')

if saveFigFlag
    figure(OptoOnSpont)
    saveas(gcf,[mainDir '\Figures\PopPlots\OptoOnSpont'],'fig');
    screen2png([mainDir '\Figures\PopPlots\OptoOnSpont']);
end
end