function data = PlotIndivROIs(data, SessionData, RespLag, selection, alpha,IndivROIFigFlag)
% Plot responses for all cells.
% Plot direction tuning if there is a significant visual reponse

%% Preamble
fprintf('\n')
tv = datestr(now, 'yyyy_mm_dd');

% RespLag = [0.2 0.2];
% selection = 1;
% alpha = 0.01;
% saveFigFlag = 0;

try
    mainDir = data.Info.mainDir;
catch
    mainDir = uigetdir('E:\TempData','Map Main Dir');
end
saveIndivROIsDir = [mainDir filesep 'Figures' filesep 'IndivROIs'];
try
    saveName = data.Info.saveName;
catch
    answer = inputdlg('savename');
    saveName = answer{1,1};
end
% saveName = [saveName '_' tv];

oldfolder = cd([mainDir, '\analysis']);
diaryName = ['AnalysisLog_tuning_' tv];
diary off; diary(diaryName)
cd(oldfolder)
disp([tv,newline,...
    'Response Window, Stim ON/OFF Offsets: ',num2str(RespLag(1)),'/',num2str(RespLag(2)),newline,...
    'selection: ' num2str(selection) newline,...
    'threshold for significant visual resp: ' num2str(alpha)])

timeStimInd = data.Info.timeStimInd;
timeVect = data.Info.timeVect;
xvect = data.Info.xvect;
xstim_opto = data.Info.xstim_opto;
xstim_vis = data.Info.xstim_vis;
tAna = timeVect>data.Info.xstim_vis(1)+RespLag(1) & timeVect<data.Info.xstim_vis(2)+RespLag(2);
tBase = timeVect(1) & timeVect<data.Info.xstim_vis(1);

color_bp = [0.8500 0.3250 0.0980];
color_bn = [0      0.4470 0.7410];


plotDirOrder=[5 2 3 6 9 8 7 4 1]; % in the middle of the plot is gonna be the no-vis trial
try
    sfValues = data.Info.sfValues;
    tfValues = data.Info.tfValues;
catch
    sfValues = SessionData.TrialSettings.GUI1_LED_CMa.spatFreq;
    tfValues = SessionData.TrialSettings.GUI1_LED_CMa.tempFreq;
end

dirValues = data.Info.directions;
nSF = length(sfValues);
nTF = length(tfValues);
nDir = length(dirValues);
captionText =[];
for j =  1:nSF
    captionText=cat(2,captionText,{['sf=', num2str(sfValues(j)) ,'cpd, tf=', num2str(tfValues(j)), 'Hz' ]});
end

if selection == 1
    if ~isfield(data,'beadID') % Only for verisions of "IDfying_redCells.m" anterior to June 15, 2020
        for i = 1:size(data.beads_pos,1)
            ROInum = data.beads_pos(i,1);
            data.beads_pos(i,3) = sum(data.s2p.iscell(1:ROInum));
        end
        selec = zeros(1,size(data.data_sorted,1));
        %    for i = 1:length(data.beads_pos)
        selec(data.beads_pos(:,3)) = 1;
        %    end
    else
        selec =  data.beadID.bead_pos;
    end
elseif selection == 0
    selec = ones(1,size(data.data_sorted,1));
end
selec = logical(selec);

%% Caluclate the impact of laser on spont or visual-evoked activity
meanResp = nanmean(data.data_sorted,6);
stdResp = std(data.data_sorted,[],6,'omitnan')./sqrt(size(data.data_sorted,6));
if nSF == 2
    RespAvg = squeeze(nanmean(meanResp(:,tAna,:,:,:),2));
    StdAvg = squeeze(nanmean(stdResp(:,tAna,:,:,:),2));
elseif nSF == 1
    RespAvg = nanmean(meanResp(:,tAna,:,:,:),2);
    StdAvg = nanmean(stdResp(:,tAna,:,:,:),2);
    RespAvg = squeeze(repmat(RespAvg,1,1,1,2,1));
    StdAvg = squeeze(repmat(StdAvg,1,1,1,2,1));
    RespAvg(:,:,2,:) = NaN(1,1);
    StdAvg(:,:,2,:) = NaN(1,1);
end

%% Plot data for individual ROI
dir_tuning = NaN(size(meanResp,1),nDir+1,nSF,2);
% if isfield(data.bestVisResp,'h_vis')
%     h_vis = data.bestVisResp.h_vis;
% else
h_vis = zeros(size(meanResp,1),nSF);
% end

MF = zeros(size(dir_tuning,1),nSF,2);
% Error = zeros(size(y,1),nSF,2); diff = zeros(length(directions),2);
gof = zeros(size(dir_tuning,1),nSF,2);
sigma = zeros(size(dir_tuning,1),nSF,2);

for i = 1:size(meanResp,1)
    if selec(i)
        color = color_bp;
    else
        color = color_bn;
    end
    allResp = cat(2,RespAvg(i,2:9,1,1),RespAvg(i,2:9,2,1));
    
    %     if ~find(h_vis)
    [~,I_all] = sort(allResp,'descend','MissingPlacement','last');
    
    h1 = 0; j = 1;
    while h1 == 0 && j < 4 %size(I,2)+1
        if I_all(j)+1 <= size(RespAvg,2)
            %             disp(['ROI', num2str(i),'; SF = 1'])
            vis_resp = squeeze(mean(data.data_sorted(i,tAna,I_all(j)+1,1,1,:),2));
            vis_base = squeeze(mean(data.data_sorted(i,tBase,I_all(j)+1,1,1,:),2));
            [h1, p1] = ttest(vis_resp, vis_base,'Alpha',alpha);
            sf = 1;
        else
            %             disp(['ROI', num2str(i),'SF 2'])
            vis_resp = squeeze(mean(data.data_sorted(i,tAna,I_all(j)+1-(size(RespAvg,2)-1),2,1,:),2));
            vis_base = squeeze(mean(data.data_sorted(i,tBase,I_all(j)+1-(size(RespAvg,2)-1),2,1,:),2));
            [h1, p1] = ttest(vis_resp, vis_base,'Alpha',alpha);
            sf = 2;
        end
        j = j+1;
    end
    h_vis(i,sf) = h1;
    %     end
    if h1
        n_run(i) = j-1;
        p_vis(i) = p1;
        if sf==1
            BestRespIdx = I_all(j-1)+1;
        elseif sf==2
            BestRespIdx = I_all(j-1)+1-(size(RespAvg,2)-1);
        end
    else
        BestRespIdx = 0;
    end
    
    for iii = 1:nSF
        if IndivROIFigFlag
            maxYval = max(max(max(meanResp(i,:,:,iii,:)+stdResp(i,:,:,iii,:),[],2)));
            minYval = min(min(min(meanResp(i,:,:,iii,:)-stdResp(i,:,:,iii,:),[],2)));
            figure;
            for ii = 1:nDir+1 % actual # of vis stim directions + no vis stim trial
                subtightplot(3,3,plotDirOrder(ii),[],0.1,0.1); hold on
                g1 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                    [meanResp(i,xvect,ii,iii,1)+stdResp(i,xvect,ii,iii,1)...
                    flip(meanResp(i,xvect,ii,iii,1)-stdResp(i,xvect,ii,iii,1))],...
                    [0 0 0],'LineStyle','none');
                set(g1,'facealpha',.5)
                g2 = fill([timeVect(xvect) flip(timeVect(xvect))],...
                    [meanResp(i,xvect,ii,iii,2)+stdResp(i,xvect,ii,iii,2)...
                    flip(meanResp(i,xvect,ii,iii,2)-stdResp(i,xvect,ii,iii,2))],...
                    color,'LineStyle','none');
                set(g2,'facealpha',.5)
                g3 = plot(timeVect(xvect),meanResp(i,xvect,ii,iii,1),'k');
                g4 = plot(timeVect(xvect),meanResp(i,xvect,ii,iii,2),'Color',color);
                ylim([minYval maxYval]);
                xlim([timeVect(find(xvect,1,'first')) timeVect(find(xvect,1,'last'))])
                xl = xlim; yl = ylim;
                plot([xl(1) xl(2)],[0 0],'--k')
                plot([xstim_vis(1) xstim_vis(2)],[yl(2)-2*0.1*yl(2) yl(2)-2*0.1*yl(2)],'r','LineWidth',4)
                if ii == 1
                    leg = legend([g3,g4],{'No Laser','Laser'},'Location','Northeast');
                    legend('boxoff')
                    leg.ItemTokenSize = [10,1];
                else
                    plot([xstim_opto(1) xstim_opto(2)],[yl(2)-0.1*yl(2) yl(2)-0.1*yl(2)],'k','LineWidth',4)
                end
                if plotDirOrder(ii) == 7
                    xlabel('Time (s)')
                    ylabel('dF/F')
                else
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])
                end
                if ii == BestRespIdx && iii == sf
                    box on
                    text(double(xl(2)),double(yl(2)),['p = ',num2str(p1,2)],...
                        'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8)
                end
            end
            sgtitle({['ROI ' num2str(i) ' | Mean traces'],captionText{iii}})
            screen2png([saveIndivROIsDir '\ROI' num2str(i) '_' num2str(2*iii-1)]);
            pause(0.1)
            close
        end
        
        
        [~, I]= sort(RespAvg(i,2:9,iii,1),'descend');
        Dir_sorted = [mod(I(1)-4,8) mod(I(1)-3,8) mod(I(1)-2,8) mod(I(1)-1,8) I(1) mod(I(1)+1,8) mod(I(1)+2,8) mod(I(1)+3,8) mod(I(1)+4,8)];
        k = find(Dir_sorted == 0);
        Dir_sorted(k) = 8;
        Dir_sorted = Dir_sorted+1;
        
        if IndivROIFigFlag && h_vis(i,iii)
            if iii == sf
                figure;
                subplot(2,2,1); hold on
                h1 = fill([1:8 flip(1:8)],[RespAvg(i,2:9,iii,1)+StdAvg(i,2:9,iii,1)...
                    flip(RespAvg(i,2:9,iii,1)-StdAvg(i,2:9,iii,1))],...
                    [0 0 0],'LineStyle','none');
                set(h1,'facealpha',.5)
                h2 = fill([1:8 flip(1:8)],[RespAvg(i,2:9,iii,2)+StdAvg(i,2:9,iii,2)...
                    flip(RespAvg(i,2:9,iii,2)-StdAvg(i,2:9,iii,2))],...
                    color,'LineStyle','none');
                set(h2,'facealpha',.5)
                p1 = plot(RespAvg(i,2:9,iii,1),'ko-');
                p2 = plot(RespAvg(i,2:9,iii,2),'ro-','Color',color);
                xticks(1:8); xticklabels(dirValues); xtickangle(45)
                xlim([0 9])
                ylabel('dFF')
                legend([p1, p2],{'Laser OFF','Laser ON'})
                legend box off
                title('Real data')
                subplot(2,2,3); hold on
                h1 = fill([1:9 flip(1:9)],[RespAvg(i,Dir_sorted,iii,1)+StdAvg(i,Dir_sorted,iii,1)...
                    flip(RespAvg(i,Dir_sorted,iii,1)-StdAvg(i,Dir_sorted,iii,1))],...
                    [0 0 0],'LineStyle','none');
                set(h1,'facealpha',.5)
                h2 = fill([1:9 flip(1:9)],[RespAvg(i,Dir_sorted,iii,2)+StdAvg(i,Dir_sorted,iii,2)...
                    flip(RespAvg(i,Dir_sorted,iii,2)-StdAvg(i,Dir_sorted,iii,2))],...
                    color,'LineStyle','none');
                set(h2,'facealpha',.5)
                plot(RespAvg(i,Dir_sorted,iii,1),'ko-')
                plot(RespAvg(i,Dir_sorted,iii,2),'o-','Color',color)
                xticks(1:9); xticklabels({dirValues(Dir_sorted-1),'180'});xtickangle(45)
                xlabel('Direction of Gratings');
                ylabel('dFF')
                title('Dir. shifted data')
                
                %             subplot(3,1,3); hold on
                %             scatter(RespAvg(i,2:9,iii,1),ImpactOfLaser(i,2:9,iii,1),...
                %                 'MarkerFaceColor',color,'MarkerEdgeColor','none');
                %             LaserImpactmdl= fitlm(RespAvg(i,2:9,iii,1)',ImpactOfLaser(i,2:9,iii,1)','linear');
                %             [LaserImpactpred, LaserImpactci] = predict(LaserImpactmdl,RespAvg(i,2:9,iii,1)');
                %             plot(RespAvg(i,2:9,iii,1)',LaserImpactpred,'k-');
                %             plot(RespAvg(i,2:9,iii,1)',LaserImpactci, 'color',[0.5 0.5 0.5],'LineStyle',':');
                %             xlabel('Magnitude of Visual Response (dFF)')
                %             ylabel({'Impact of laser' ,'(delta dFF)'})
                %             title({['slope: ' num2str(table2array(LaserImpactmdl.Coefficients(2,1)),3) ' | adjusted r2 = ' num2str(LaserImpactmdl.Rsquared.Adjusted,3)]})
                %
                %             sgtitle({['ROI ' num2str(i) ' | Tuning and laser magn = f(magn vis resp)'],captionText{iii}})
                %             screen2png([saveIndivROIsDir '\ROI' num2str(i) '_' num2str(iii*2)]);
                %             pause(0.1)
                %             close
            end
        end
        dir_tuning(i,:,iii,:) = RespAvg(i,Dir_sorted,iii,:);
    end
    
    % -----------------------------
    if nSF == 2
y = squeeze(circshift(dir_tuning(i,1:8,:,:),-2,2));
    else
    y(1:8,1,:)=circshift(dir_tuning(i,1:8,:,:),-2,2);
    end
    directions = -90:45:225;
    
    options = fitoptions('gauss2');
    options.Lower = [0 0 0 0 180 0];
    options.Upper = [Inf 0 Inf Inf 180 Inf];
    
    for  iii = 1:nSF
        if h_vis(i,iii)
            y1 = squeeze(y(:,iii,1));
            [f1, gof1] = fit(directions.',y1,'gauss2',options);
            y2 = squeeze(y(:,iii,2));
            [f2, gof2] = fit(directions.',y2,'gauss2',options);
            
            a1(i,iii,1) = f1.a1;
            a2(i,iii,1) = f1.a2;
            a1(i,iii,2) = f2.a1;
            a2(i,iii,2) = f2.a2;
            gof(i,iii,1) = gof1.rsquare;
            gof(i,iii,2) = gof2.rsquare;
            sigma(i,iii,1) = f1.c1;
            sigma(i,iii,2) = f2.c1;
            
            if IndivROIFigFlag
                subplot(2,2,2); hold on
                plot(f1,'k--',directions,y1,'k')
                xlabel(''); ylabel('')
                xticks(-90:45:225); xticklabels('')
                title(['Fit Laser OFF | adj. R2 = ',num2str(gof1.adjrsquare,4)])
                subplot(2,2,4); hold on
                plot(f2,'r--',directions,y2,'r')
                xlabel('angle from best direction'); ylabel('')
                xticks(-90:45:225); xtickangle(45)
                title(['Fit Laser ON | adj. R2 = ',num2str(gof2.adjrsquare,4)])
                sgtitle({['ROI ' num2str(i) ' | ',captionText{iii}]})
                %                 screen2png([saveIndivROIsDir '\ROI' num2str(i) '_' num2str(iii)]);
                %                 pause(0.1)
                %                 close
                %                 sgtitle({['ROI ' num2str(i) ' | Tuning and laser magn = f(magn vis resp)'],captionText{iii}})
                screen2png([saveIndivROIsDir '\ROI' num2str(i) '_' num2str(iii*2)]);
                pause(0.1)
                close
            end
            
            MF(i,iii,1) = (f2.a1-f1.a1)/(f1.a1+f2.a1);
            MF(i,iii,2) = (f2.a2-f1.a2)/(f1.a2+f2.a2);
            
        end
    end
end

data.fit.MF = MF;
% data.fit.Error = Error;
data.fit.gof = gof;
data.fit.sigma = sigma;h_vis = logical(h_vis);
data.fit.h_vis = h_vis;
data.fit.alpha = alpha;

data.dir_tuning = dir_tuning;

% a = exist([mainDir, '\analysis\', saveName,'.mat'],'file');
% if a  == 2
%     saveName = [saveName '_' tv];
% end
save([mainDir, '\analysis\', saveName],'data','-append')
% save([mainDir, filesep 'analysis' filesep 'dir_tuning_' saveName],'tuning')
disp(['data tuning and fit saved with name: ' saveName]);

diary
diary off

clear global mainDir
clear global saveName
end

