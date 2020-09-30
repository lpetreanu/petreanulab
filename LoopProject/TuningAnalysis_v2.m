% load('E:\TempData\CMloop\tuning_CMloopAll_30April2020.mat')

function TuningData = TuningAnalysis_v2(data,do_MF,looped_cells,saveFig,varargin)
%% dir_tuning does not contain the no-visual stim stim type

%% Preambule - set parameters
DirSelection = 1;
% 0: no selection based on cell direction selectivity
% 1: selection based on the error of the 2 gaussian fit
% 2: selection based on ANOVA (to be implemented)
error_th = 0.2; % for the error of the fit

Fig_IndivROI = 0; % if "PlotIndivROIs.m" has been run, these plots have been generated already

% - fetch the data
dir_tuning = data.dir_tuning;
try
    h_vis = cat(1,data.fit(:).h_vis);
    %     h_vis = data.fit.h_vis;    
catch
%     h_vis = data.fit.h_vis;    
    h_vis = data.h_vis;
end

try
    bead = data.s2p.bead_pos;
catch % if "BestResp_BeadsSelec_v5.m" ran before 15 June 2020 version
    bead = false(size(dir_tuning,1),1);
    if isfield(data, 'beads_pos')
        beads_pos = data.beads_pos;
        bead(beads_pos(:,3)) = true;
    else
        bead = false(size(dir_tuning,1),1);
    end
end
% dirValues = data.Info.directions;
sfValues = data.Info.sfValues;
tfValues = data.Info.tfValues;
nSF = length(sfValues);
if nSF == 2
    hVis_cat = cat(1,h_vis(:,1),h_vis(:,2));
    bead_cat = [bead ; bead];
else
    hVis_cat = h_vis;
    bead_cat = bead;
end


if Fig_IndivROI
    captionText =[];
    for j =  1:nSF
        captionText=cat(2,captionText,{['sf=', num2str(sfValues(j)) ,'cpd, tf=', num2str(tfValues(j)), 'Hz' ]});
    end
    mainDir = uigetdir('E:\TempData','Choose Main dir');
    saveIndivROIsDir = [mainDir filesep 'Figures' filesep 'IndivROIs'];
end

if looped_cells
color_bp = [0.8500 0.3250 0.0980];
color_bn = [0      0.4470 0.7410];
else
 color_bp = [0.4940 0.1840 0.5560]; 
 color_bn = color_bp;
end
%%

%% 1) Modulation Factor
if ~isfield(data,'fit') || do_MF
    y = circshift(dir_tuning(:,1:8,:,:),-2,2);
    directions = -90:45:225;
    
    options = fitoptions('gauss2');
    options.Lower = [0 0 0 0 180 0];
    options.Upper = [Inf 0 Inf Inf 180 Inf];
    
    MF = zeros(size(y,1),nSF,2); Error = zeros(size(y,1),nSF,2); diff = zeros(length(directions),2);
    gof = zeros(size(y,1),nSF,2);
    sigma = zeros(size(y,1),nSF,2);
    f = waitbar(0,'modeling of individual ROI in progress');
    for i = 1:size(y,1)
        waitbar(i/size(y,1),f,[num2str(i/size(y,1)*100,2),'%'])
        for  ii = 1:nSF
            if h_vis(i,ii)
                y1 = squeeze(y(i,:,ii,1));
                [f1, gof1] = fit(directions.',y1.','gauss2',options);
                y2 = squeeze(y(i,:,ii,2));
                [f2, gof2] = fit(directions.',y2.','gauss2',options);
                
                gof(i,ii,1) = gof1.rsquare;
                gof(i,ii,2) = gof2.rsquare;
                sigma(i,ii,1) = f1.c1;
                sigma(i,ii,2) = f2.c1;
                
                if Fig_IndivROI
                    figure
                    subplot(1,2,1); hold on
                    plot(f1,'k--',directions,y1,'k')
                    xlabel('directions')
                    xticks(-90:45:225); xtickangle(45)
                    ylabel('dF/F')
                    title(['Laser OFF | adj. R2 = ',num2str(gof1.adjrsquare,4)])
                    subplot(1,2,2); hold on
                    plot(f2,'r--',directions,y2,'r')
                    xlabel('directions')
                    xticks(-90:45:225); xtickangle(45)
                    title(['Laser ON | adj. R2 = ',num2str(gof2.adjrsquare,4)])
                    sgtitle({['ROI ' num2str(i) ' | ',captionText{ii}]})
                    screen2png([saveIndivROIsDir '\ROI' num2str(i) '_' num2str(ii)]);
                    pause(0.1)
                    close
                end
                
                MF(i,ii,1) = (f2.a1-f1.a1)/(f1.a1+f2.a1);
                MF(i,ii,2) = (f2.a2-f1.a2)/(f1.a2+f2.a2);
                
                for iii = 1:length(directions)
                    diff(iii,1) = (y1(iii) - f1(directions(iii)))^2;
                    diff(iii,2) = (y2(iii) - f2(directions(iii)))^2;
                end
                
                Error(i,ii,1) = sum(diff(:,1))./sum((y1-mean(y1)).^2);
                Error(i,ii,2) = sum(diff(:,2))./sum((y2-mean(y2)).^2);
            end
        end
    end
    close(f)
    
    % bla = cat(1,Error(:,1,1),Error(:,2,1));
    % bla2 = cat(1,Error(:,1,2),Error(:,2,2));
    % gofLOFF_cat = cat(1,gof(:,1,1),gof(:,2,1));
    % gofLON_cat = cat(1,gof(:,1,2),gof(:,2,2));
    % figure;
    % hold on
    % plot(gofLOFF_cat(hVis_cat),1-bla(hVis_cat),'ko')
    % plot(gofLON_cat(hVis_cat),1-bla2(hVis_cat),'ro')
    data.fit.MF = MF;
    data.fit.Error = Error;
    data.fit.gof = gof;
    data.fit.sigma = sigma;
else
    MF = cat(1,data.fit(:).MF);
%     Error = data.fit.Error;
    gof = cat(1,data.fit(:).gof);
    sigma = cat(1,data.fit(:).sigma);
end
% fitselec = Error(:,:,1) < error_th & Error(:,:,2) < error_th;
fitselec = gof(:,:,1) > 1-error_th & gof(:,:,2) > 1-error_th;

%%  This was plotting the modulation factor (MF) without segregating bead+
% % and bead-

% fit_selected = mean(MF(h_vis(:,1)&k,:));
% fit_selected_sem = std(MF(h_vis(:,1)&k,:))./sqrt(sum(h_vis(:,1)&k));
% all = mean(MF(h_vis(:,1),:));
% all_sem = std(MF(h_vis(:,1),:))./sqrt(sum(h_vis(:,1)));

% figure;
% % Bead+ and bead- not segregated
% subplot(1,4,1); hold on
% b = bar([fit_selected' all']);
% ylabel('Opto Modulation Factor')
% xticks([1 2]); xticklabels({'Pref.','Antipref.'})
% xlabel('direction')
% for i = 1:2
% plot([i+b(1).XOffset i+b(1).XOffset],...
%     [fit_selected(i)+fit_selected_sem(i) fit_selected(i)-fit_selected_sem(i)],'k-')
% plot([i+b(2).XOffset i+b(2).XOffset],...
%     [all(i)+all_sem(i) all(i)-all_sem(i)],'k-')
% end
% scatter(ones(1,sum(h_vis(:,1)&k))+2*b(1).XOffset,MF(h_vis(:,1)&k,1),...
%     'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% scatter(2*ones(1,sum(h_vis(:,1)&k))+2*b(1).XOffset,MF(h_vis(:,1)&k,2),...
%     'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% scatter(ones(1,size(y(h_vis(:,1)),1))+2*b(2).XOffset,MF(h_vis(:,1),1),...
%     'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% scatter(2*ones(1,size(y(h_vis(:,1)),1))+2*b(2).XOffset,MF(h_vis(:,1),2),...
%     'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% legend(b,{['selected data, n = ',num2str(size(k,1))],['all data, n = ',num2str(size(y,1))]},...
%     'Location','southwest')

% subplot(1,4,2); hold on
% xl = xlim; yl = ylim;
% text({'mean +- sem',[num2str(fit_selected(1),3) ' +- ' num2str(fit_selected_sem(1),3)],...
%     [num2str(all(1),3) ' +- ' num2str(all_sem(1),3)]})
%%  segregate bead+ and bead-

if nSF == 2
    MF_cat = squeeze(cat(1,MF(:,1,:),MF(:,2,:)));
    k_cat = cat(1,fitselec(:,1),fitselec(:,2));
else
    MF_cat = squeeze(MF);
    k_cat = fitselec;
end

fit_sel_Bpos = mean(MF_cat(hVis_cat&k_cat&bead_cat,:),1);
fit_sel_Bpos_sem = std(MF_cat(hVis_cat&k_cat&bead_cat,:),[],1)./sqrt(sum(hVis_cat&k_cat&bead_cat));
fit_sel_Bneg = mean(MF_cat(hVis_cat&k_cat&~bead_cat,:),1);
fit_sel_Bneg_sem = std(MF_cat(hVis_cat&k_cat&~bead_cat,:),[],1)./sqrt(sum(hVis_cat&k_cat&~bead_cat));
all_Bpos = mean(MF_cat(hVis_cat&bead_cat,:),1);
all_Bneg = mean(MF_cat(hVis_cat&~bead_cat,:),1);
all_Bpos_sem = std(MF_cat(hVis_cat&bead_cat,:),[],1)./sqrt(sum(hVis_cat&bead_cat));
all_Bneg_sem = std(MF_cat(hVis_cat&~bead_cat,:),[],1)./sqrt(sum(hVis_cat&~bead_cat));

MF_Figure = figure;
if looped_cells
    if isfield(data, 'beads_pos') && sum(hVis_cat&k_cat&bead_cat,1)>=4
        if ~adtest(MF_cat(hVis_cat&k_cat&bead_cat,1)) && ~adtest(MF_cat(hVis_cat&k_cat&~bead_cat,1))
            [~,p]= ttest2(MF_cat(hVis_cat&k_cat&bead_cat,:),MF_cat(hVis_cat&k_cat&~bead_cat,:));
            testName = 'ttest';
        else
            [p,~]= ranksum(MF_cat(hVis_cat&k_cat&bead_cat,1),MF_cat(hVis_cat&k_cat&~bead_cat,1));
            testName = 'Wilcoxon';
        end
    else
        testName = ''; p=[];
    end
else
        testName = 'No stat'; p=[];
end

subplot(1,4,1); hold on
b2 = bar([fit_sel_Bpos' fit_sel_Bneg'],'FaceColor','flat');
drawnow
for i = 1:2
    plot([i+b2(1).XOffset i+b2(1).XOffset],...
        [fit_sel_Bpos(i)+fit_sel_Bpos_sem(i) fit_sel_Bpos(i)-fit_sel_Bpos_sem(i)],'k-')
    plot([i+b2(2).XOffset i+b2(2).XOffset],...
        [fit_sel_Bneg(i)+fit_sel_Bneg_sem(i) fit_sel_Bneg(i)-fit_sel_Bneg_sem(i)],'k-')
end
b2(1).CData = color_bp;
b2(2).CData = color_bn;
scatter(ones(1,sum(hVis_cat&k_cat&bead_cat))+2*b2(1).XOffset,MF_cat(hVis_cat&k_cat&bead_cat,1),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(2*ones(1,sum(hVis_cat&k_cat&bead_cat))+2*b2(1).XOffset,MF_cat(hVis_cat&k_cat&bead_cat,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(ones(1,sum(hVis_cat&k_cat&~bead_cat))+2*b2(2).XOffset,MF_cat(hVis_cat&k_cat&~bead_cat,1),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
scatter(2*ones(1,sum(hVis_cat&k_cat&~bead_cat))+2*b2(2).XOffset,MF_cat(hVis_cat&k_cat&~bead_cat,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
lg1 = legend(b2,{['bead+, n = ',num2str(sum(hVis_cat&k_cat&bead_cat))],...
    ['bead-, n = ',num2str(sum(hVis_cat&k_cat&~bead_cat))]},...
    'Location','none');
xticks([1 2]); xticklabels({'Pref.','Antipref.'})
xlabel('direction')
title({['Selected data, (fit error < ' num2str(error_th) ')'],[testName ', p = ' num2str(p,3)]})
ylabel('Opto Modulation Factor')

subplot(1,4,3); hold on
if looped_cells
    if isfield(data, 'beads_pos') && sum(hVis_cat&k_cat&bead_cat,1)>=4
        if ~adtest(MF_cat(hVis_cat&bead_cat,1)) && ~adtest(MF_cat(hVis_cat&~bead_cat,1))
            [~,p]= ttest2(MF_cat(hVis_cat&bead_cat,:),MF_cat(hVis_cat&~bead_cat,:));
        else
            p = ranksum(MF_cat(hVis_cat&bead_cat,1),MF_cat(hVis_cat&~bead_cat,1));
        end
    end
else
    p = 'NaN';
end

b3 = bar([all_Bpos' all_Bneg'],'FaceColor','flat');
drawnow
for i = 1:2
    plot([i+b3(1).XOffset i+b3(1).XOffset],...
        [all_Bpos(i)+all_Bpos_sem(i) all_Bpos(i)-all_Bpos_sem(i)],'k-')
    plot([i+b3(2).XOffset i+b3(2).XOffset],...
        [all_Bneg(i)+all_Bneg_sem(i) all_Bneg(i)-all_Bneg_sem(i)],'k-')
end
b3(1).CData = color_bp; b3(2).CData = color_bn;
scatter(ones(1,sum(hVis_cat&bead_cat))+2*b3(1).XOffset,MF_cat(hVis_cat&bead_cat,1),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(2*ones(1,sum(hVis_cat&bead_cat))+2*b3(1).XOffset,MF_cat(hVis_cat&bead_cat,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(ones(1,sum(hVis_cat&~bead_cat))+2*b3(2).XOffset,MF_cat(hVis_cat&~bead_cat,1),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
scatter(2*ones(1,sum(hVis_cat&~bead_cat))+2*b3(2).XOffset,MF_cat(hVis_cat&~bead_cat,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
lg2 = legend(b3,{['bead+, n = ',num2str(sum(hVis_cat&bead_cat))],...
    ['bead-, n = ',num2str(sum(hVis_cat&~bead_cat))]},...
    'Location','none');
xticks([1 2]); xticklabels({'Pref.','Antipref.'})
xlabel('direction')
title({'All data',['ttest, p = ',num2str(p,3)]})

subplot(1,4,2); hold on
axis off
lg1.Position = [0.3 0.8 0.10 0.1];
lg1.ItemTokenSize = [10,4];
lg2.Position = [0.4 0.2 0.10 0.1];
lg2.ItemTokenSize = [10,4];
%
% tuning.modulation_factor.MF = MF;
% tuning.modulation_factor.error_selec = k;
subplot(1,4,4); hold on
if adtest(MF_cat(hVis_cat&k_cat,1)) 
    p=signrank(MF_cat(hVis_cat&k_cat,1));
    testName = 'Wilcoxon';
else
    [~,p]=ttest(MF_cat(hVis_cat&k_cat,1));
    testName = 'ttest';
end

violin({MF_cat(hVis_cat&k_cat,1),MF_cat(hVis_cat&k_cat,2)},'facecolor',[0.5 0.5 0.5],'edgecolor','none','FaceAlpha',1);
plot(xlim,[0 0],'k--')
xlabel('direction','FontSize',10)
xticks([1 2])
xticklabels([{'Pref.'},{'Antipref.'}])
title({'All fit selected cells',[testName ', p = ' num2str(p,3)]},'FontSize',10);
set(gcf,'Units','Normalized','Position',[0.05 0.05 0.7 0.5])

%%

%% Starting the proper tuning analysis now. Selection
if DirSelection == 1 % select cells based on the double gaussian
    dirSel = k_cat;
elseif DirSelection == 0 % don't select cells based on their tuning properties
    dirSel = true(size(hVis_cat,1),1);
elseif DirSelection == 2
    %     sel = anovaSel_cat;
end
%%

%% 2) Direction tuning
% normalize amplitude to pref vis stim in Laser OFF control condition
max1 = dir_tuning(:,5,:,1);
max1 = repmat(max1,1,9,1,2);
norm_ctrl = dir_tuning./max1;
normCtrl_merged = squeeze(cat(1,norm_ctrl(:,:,1,:),norm_ctrl(:,:,2,:)));

%% 2.1 segregate SF/TF, but not bead+ vs bead- (check SF/TF stim properties)
% test(:,1,:) = squeeze(nanmean(norm_ctrl(h_vis(:,1),:,1,:),1));
% test(:,2,:) = squeeze(nanmean(norm_ctrl(h_vis(:,2),:,2,:),1));
%
% sem(:,1,:) = std(norm_ctrl(h_vis(:,1),:,1,:),[],1,'omitnan')./sqrt(sum(h_vis(:,1)));
% sem(:,2,:) = std(norm_ctrl(h_vis(:,2),:,2,:),[],1,'omitnan')./sqrt(sum(h_vis(:,2)));
%
% merged_mean = squeeze(nanmean(normCtrl_merged(hVis_cat,:,:)));
% merged_sem = squeeze(std(normCtrl_merged(hVis_cat,:,:)))./sqrt(sum(hVis_cat));
%
% DirTunAll_Fig = figure;
% subplot(1,4,1); hold on
% p1 = plot(test(:,1,1),'r','LineWidth',1);
% p2 = plot(test(:,1,2),'ro--','LineWidth',1);
% p3 = plot(test(:,2,1),'b','LineWidth',1);
% p4 = plot(test(:,2,2),'bo--','LineWidth',1);
% p5 = plot(merged_mean(:,1),'k','LineWidth',2);
% p6 = plot(merged_mean(:,2),'k--','LineWidth',2);
%
% for i = 1:size(test,1)
%     plot([i i], [test(i,1,1)+sem(i,1,1) test(i,1,1)-sem(i,1,1)],'r-')
%     plot([i i], [test(i,1,2)+sem(i,1,2) test(i,1,2)-sem(i,1,2)],'r-')
%     plot([i i], [test(i,2,1)+sem(i,2,1) test(i,2,1)-sem(i,2,1)],'b')
%     plot([i i], [test(i,2,2)+sem(i,2,2) test(i,2,2)-sem(i,2,2)],'b-')
%     plot([i i], [merged_mean(i,1)+merged_sem(i,1) merged_mean(i,1)-merged_sem(i,1)],'k-')
%     plot([i i], [merged_mean(i,2)+merged_sem(i,2) merged_mean(i,2)-merged_sem(i,2)],'k-')
% end
% xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
% title('Averaged data')
% xl=xlim;yl=ylim;
% text(p1.XData(1),yl(2),['all data, n = ',num2str(sum(hVis_cat))],...
%     'HorizontalAlignment','Left','VerticalAlignment','top')
% lgd = legend([p1,p2,p3,p4,p5,p6],{'SF/TF 1, no laser','SF/TF 1, laser','SF/TF 2, no laser','SF/TF 2, laser',...
%     'All SFs/TFs, no laser','All SFs/TFs, laser'},'location','none');
% subplot(1,4,2); axis off
% lgd.Position =  [0.3 0.5 0.15 0.15];
%
%
% subplot(1,4,3); hold on
% plot(test(:,1,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
% plot(test(:,1,2),'r','LineWidth',2)
% xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
% title('SF/TF 1, all data')
%
%
% subplot(1,4,4); hold on
% plot(test(:,2,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
% plot(test(:,2,2),'b','LineWidth',2)
% xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
% title('SF/TF 2, all data')
%
% for i = 1:sum(h_vis(:,1))
%     if h_vis(i,1)
%         subplot(1,4,3)
%         lh = plot(norm_ctrl(i,:,1,2));
%         lh.Color=[1,0,0,0.2];
%     end
%     if h_vis(i,2)
%         subplot(1,4,4)
%         lh = plot(norm_ctrl(i,:,2,2));
%         lh.Color=[0,0,1,0.2];
%     end
% end
%
% subplot(1,4,3)
% xl=xlim;yl=ylim;
% text(lh.XData(1),yl(2),['n = ',num2str(sum(h_vis(:,1)))],...
%     'HorizontalAlignment','Left','VerticalAlignment','top')
%
% subplot(1,4,4)
% xl=xlim;yl=ylim;
% text(lh.XData(1),yl(2),['n = ',num2str(sum(h_vis(:,2)))],...
%     'HorizontalAlignment','Left','VerticalAlignment','top')
%
% set(gcf,'Units','Normalized','Position',[0.05 0.4 0.7 0.5])
% sgtitle('Impact of Opto stim of FB on tuning curve')
%%

%% 2.2 Segregate bead+ and bead neg
% some calculations
tuning_bp = squeeze(mean(normCtrl_merged(hVis_cat&bead_cat&dirSel,:,:)));
tuning_bn = squeeze(mean(normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,:)));
tuning_bp_sem = squeeze(std(normCtrl_merged(hVis_cat&bead_cat&dirSel,:,:))./sqrt(sum(hVis_cat&bead_cat&dirSel)));
tuning_bn_sem = squeeze(std(normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,:))./sqrt(sum(hVis_cat&~bead_cat&dirSel)));

deltaFR = squeeze(dir_tuning(:,5,:,2) - dir_tuning(:,5,:,1));
ratioFR = squeeze(dir_tuning(:,5,:,2) ./ dir_tuning(:,5,:,1));
deltaFR = cat(1,deltaFR(:,1),deltaFR(:,2));
ratioFR = cat(1,ratioFR(:,1),ratioFR(:,2));
%  ------------- polar representation of tuning ---------------------------
%%  ------------ pop descriptive directionality stats ---------------------
% max2 = dir_tuning(:,5,:,2);
% max2 = repmat(max2,1,9,1,2);
% norm_opto = dir_tuning./max2;
% normOpto_cat = squeeze(cat(1,norm_opto(:,:,1,:),norm_opto(:,:,2,:)));

PolarFig = figure;
pos1 = [0.05 0.2 0.4 0.75];
subplot('Position',pos1);
p1 = polarplot(deg2rad(0:45:360),tuning_bn(:,1),'k');
hold on
p2 =polarplot(deg2rad(0:45:360),tuning_bn(:,2),'Color',color_bn);
% polarplot(deg2rad(0:45:360),squeeze(mean(normOpto_cat(hVis_cat&~bead_cat,:,2))),'Color',color_bn);
% % apparently need to re-align to pref direction BEFORE normalization ----
pos2 = [0.56 0.2 0.4 0.75];
subplot('Position',pos2);
p3 = polarplot(deg2rad(0:45:360),tuning_bp(:,1),'k');
hold on
p4 = polarplot(deg2rad(0:45:360),tuning_bp(:,2),'Color',color_bp);

alpha = deg2rad(45:45:360);
for i = 1:2
    w1 = tuning_bn(2:9,i);
    r1(i) = sum(w1.*exp(1i*alpha'),1);
    mu1(i) = angle(r1(i));
    r1(i) = abs(r1(i))./sum(w1,1);
    subplot('Position',pos1);
    if i ==1
        polarplot([mu1(i) mu1(i)],[0 r1(i)],'k');
    elseif i ==2
        polarplot([mu1(i) mu1(i)],[0 r1(i)],'Color',color_bn);
    end
    mu1(i) = rad2deg(mu1(i));
    
    w2 = tuning_bp(2:9,i);
    r2(i) = sum(w2.*exp(1i*alpha'),1);
    mu2(i) = angle(r2(i));
    r2(i) = abs(r2(i))./sum(w2,1);
    subplot('Position',pos2);
    if i == 1
        polarplot([mu2(i) mu2(i)],[0 r2(i)],'k');
    elseif i ==2
        polarplot([mu2(i) mu2(i)],[0 r2(i)],'Color',color_bp);
    end
    mu2(i) = rad2deg(mu2(i));
end
subplot('Position',pos1);
lg1 = legend([p1,p2],{'Laser OFF','Laser ON'});
title('beads-')
subplot('Position',pos2);
lg2 = legend([p3,p4],{'Laser OFF','Laser ON'});
title('beads+')
pos3 = [0.05 0.05 0.9 0.1];
subplot('Position',pos3); axis off
% subplot(2,2,[3,4]); axis off
lg1.Position = [0.2 0.05 0.15 0.05];
lg1.ItemTokenSize = [20,2];
lg2.Position = [0.65 0.05 0.15 0.05];
lg2.ItemTokenSize = [20,2];
xl = xlim; yl = ylim;
text(xl(1),yl(2),{[sprintf('vector length = %.3g, angle = %.4g \n',r1(1),mu1(1)),...
    sprintf('vector length = %.3g, angle = %.4g',r1(2),mu1(2))]},...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(xl(2),yl(2),{[sprintf('vector length = %.3g, angle = %.4g \n',r2(1),mu2(1)),...
    sprintf('vector length = %.3g, angle = %.4g',r2(2),mu2(2))]},...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
sgtitle('Polar tuning curve, Laser OFF norm.')


[~,I] = max(dir_tuning(:,:,1,2),[],2);
Iopto(:,1) = I;
[~,I] = max(dir_tuning(:,:,2,2),[],2);
Iopto(:,2) = I;
for i = 1:size(dir_tuning,1)
    for ii = 1:nSF
        shift = 5-Iopto(i,ii);
        temp(i,:,ii) = circshift(dir_tuning(i,:,ii,2),shift,2);
    end
end
max2 = temp(:,5,:);
max2 = repmat(max2,1,9,1);
norm_opto = temp./max2;
normOpto_cat = squeeze(cat(1,norm_opto(:,:,1,:),norm_opto(:,:,2,:)));
normOptoCat_bn = squeeze(mean(normOpto_cat(hVis_cat&~bead_cat&dirSel,:)));
normOptoCat_bp = squeeze(mean(normOpto_cat(hVis_cat&bead_cat&dirSel,:)));
clear temp

PolarFig2 = figure;
pos1 = [0.05 0.2 0.4 0.7]; % [left bottom width height]
subplot('Position',pos1);
p1 = polarplot(deg2rad(0:45:360),tuning_bn(:,1),'k');
hold on
p2 = polarplot(deg2rad(0:45:360),normOptoCat_bn,'Color',color_bn);

pos2 = [0.56 0.2 0.4 0.7];
subplot('Position',pos2);
p3 = polarplot(deg2rad(0:45:360),tuning_bp(:,1),'k');
hold on
p4 = polarplot(deg2rad(0:45:360),normOptoCat_bp,'Color',color_bp);

w1 = normOptoCat_bn(2:9);
r1(3) = sum(w1'.*exp(1i*alpha'),1);
mu1(3) = angle(r1(3));
r1(3) = abs(r1(3))./sum(w1,2);
subplot('Position',pos1);
polarplot([deg2rad(mu1(1)) deg2rad(mu1(1))],[0 r1(1)],'k');
polarplot([mu1(3) mu1(3)],[0 r1(3)],'Color',color_bn);
mu1(3) = rad2deg(mu1(3));

w2 = normOptoCat_bp(2:9);
r2(3) = sum(w2'.*exp(1i*alpha'),1);
mu2(3) = angle(r2(3));
r2(3) = abs(r2(3))./sum(w2,2);
subplot('Position',pos2);
polarplot([deg2rad(mu2(1)) deg2rad(mu2(1))],[0 r2(1)],'k');
polarplot([mu2(3) mu2(3)],[0 r2(3)],'Color',color_bp);
mu2(3) = rad2deg(mu2(3));
subplot('Position',pos1);
lg1 = legend([p1,p2],{'Laser OFF','Laser ON'});
title('beads-')
subplot('Position',pos2);
lg2 = legend([p3,p4],{'Laser OFF','Laser ON'});
title('beads+')

pos3 = [0.05 0.05 0.9 0.1];
subplot('Position',pos3); axis off
% subplot(2,2,[3,4]); axis off
lg1.Position = [0.2 0.05 0.15 0.05];
lg1.ItemTokenSize = [20,2];
lg2.Position = [0.65 0.05 0.15 0.05];
lg2.ItemTokenSize = [20,2];
xl = xlim; yl = ylim;
text(xl(1),yl(2),{[sprintf('vector length = %.3g, angle = %.4g \n',r1(1),mu1(1)),...
    sprintf('vector length = %.3g, angle = %.4g',r1(3),mu1(3))]},...
    'HorizontalAlignment','left','VerticalAlignment','bottom')
text(xl(2),yl(2),{[sprintf('vector length = %.3g, angle = %.4g \n',r2(1),mu2(1)),...
    sprintf('vector length = %.3g, angle = %.4g',r2(3),mu2(3))]},...
    'HorizontalAlignment','right','VerticalAlignment','bottom')
sgtitle('Polar tuning curve, Opto shifted and norm')
%%
%% 2.2 Suite
DirTun_Bead_Fig = figure;
subplot(1,4,1); hold on
p7 = plot(tuning_bp(:,1),'-','Color',color_bp);
p8 = plot(tuning_bp(:,2),'--','Color',color_bp);
p9 = plot(tuning_bn(:,1),'-','Color',color_bn);
p10 = plot(tuning_bn(:,2),'--','Color',color_bn);
for i = 1:9
    plot([i i],[tuning_bp(i,1)+tuning_bp_sem(i,1) tuning_bp(i,1)-tuning_bp_sem(i,1)],'-','Color',color_bp)
    plot([i i],[tuning_bp(i,2)+tuning_bp_sem(i,2) tuning_bp(i,2)-tuning_bp_sem(i,2)],'-','Color',color_bp)
    plot([i i],[tuning_bn(i,1)+tuning_bn_sem(i,1) tuning_bn(i,1)-tuning_bn_sem(i,1)],'-','Color',color_bn)
    plot([i i],[tuning_bn(i,2)+tuning_bn_sem(i,2) tuning_bn(i,2)-tuning_bn_sem(i,2)],'-','Color',color_bn)
end
xlim([0.5 9.5]); xticks(1:2:9);
xticklabels({'-180','-90','Pref.','+90','+180'});
xlabel('Directions'); ylabel('dFF, norm')
lg3 = legend([p7, p8, p9, p10],{['bead+ laser OFF, n = ' num2str(sum(hVis_cat&bead_cat&dirSel))],...
    'bead+ laser ON',...
    ['bead- laser OFF, n = ' num2str(sum(hVis_cat&~bead_cat&dirSel))],...
    'bead- laser ON'});
title('Session Average')

subplot(1,4,3); hold on
if isfield(data, 'beads_pos')
    bla = normCtrl_merged(hVis_cat&bead_cat&dirSel,:,:);
for i = 1:sum(hVis_cat&bead_cat&dirSel)
    lh1 = plot(squeeze(bla(i,:,2)),'-','Color',color_bp);
    lh1.Color=[color_bp 0.1];
end
p11 = plot(tuning_bp(:,1),'k-','LineWidth',2);
p12 = plot(tuning_bp(:,2),'k--','LineWidth',2);
xlim([0.5 9.5]); xticks(1:2:9);
xticklabels({'-180','-90','Pref.','+90','+180'});
xlabel('Directions'); ylabel('dFF, norm')
title('All bead+ cells')
xl = xlim; yl = ylim;
% text(p7.XData(1),yl(2),['n = ',num2str(sum(hVis_cat&bead_cat))],...
%     'HorizontalAlignment','Left','VerticalAlignment','top')
clear bla
lg4 = legend([p11,p12,lh1],{'laser OFF, average','laser ON, average','laser ON, all data'});
end

if looped_cells
    subplot(1,4,4); hold on
    bla = normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,:);
    for i = 1:sum(hVis_cat&~bead_cat&dirSel)
        lh2 = plot(squeeze(bla(i,:,2)),'-','Color',color_bn);
        lh2.Color=[color_bn 0.05];
    end
    plot(tuning_bn(:,1),'k-','LineWidth',2);
    plot(tuning_bn(:,2),'k--','LineWidth',2);
    xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'});
    xlabel('Directions'); ylabel('dFF, norm')
    title('All bead- cells')
    ylim([yl(1) yl(2)]);
    lg5 = legend(lh2,'laser ON, all data');
end

subplot(1,4,2); axis off
lg3.Position = [0.3 0.7 0.13 0.1]; lg3.ItemTokenSize = [20,4];
lg4.Position = [0.3 0.4 0.13 0.1]; lg4.ItemTokenSize = [20,3];
lg5.Position = [0.3 0.1 0.13 0.1]; lg5.ItemTokenSize = [20,2];

set(gcf,'Units','Normalized','Position',[0.05 0.4 0.7 0.3])
%%
%% 2.3.1 quantification of vect. length for each cell ---------------
tuning_bp_all = normCtrl_merged(hVis_cat&bead_cat&dirSel,1:8,:);
for i = 1:size(tuning_bp_all,1)
    for ii = 1:2
        w = squeeze(tuning_bp_all(i,:,ii));
        w( w < 0) = 0;
        r = sum(w.*exp(1i*alpha),2);
        mu_bp(i,ii) = angle(r);
        r_bp(i,ii) = abs(r)./sum(w,2);
    end
end
rBp_mean = nanmean(r_bp);
rBp_sem = std(r_bp,[],1,'omitnan')./sqrt(size(r_bp,1));

if looped_cells
    tuning_bn_all = normCtrl_merged(hVis_cat&~bead_cat&dirSel,1:8,:);
    for i = 1:size(tuning_bn_all,1)
        for ii = 1:2
            w = squeeze(tuning_bn_all(i,:,ii));
            w( w < 0) = 0;
            r = sum(w.*exp(1i*alpha),2);
            mu_bn(i,ii) = angle(r);
            r_bn(i,ii) = abs(r)./sum(w,2);
        end
    end
else
    r_bn = zeros(size(r_bp));
    r_bn = zeros(size(r_bp));
end
rBn_mean = nanmean(r_bn);
rBn_sem = std(r_bn,[],1,'omitnan')./sqrt(size(r_bn,1));
    
VectorLength_Fig = figure;
subplot(1,3,1); hold on
b = bar([rBn_mean' rBp_mean']);
drawnow
for i = 1:2
    plot([i+b(1).XOffset i+b(1).XOffset],[rBn_mean(i)+rBn_sem(i) rBn_mean(i)-rBn_sem(i)],'k')
    plot([i+b(2).XOffset i+b(2).XOffset],[rBp_mean(i)+rBp_sem(i) rBp_mean(i)-rBp_sem(i)],'k')
end
legend(b,{'bead-','bead+'},'Location','southwest')
xticks([1 2]); xticklabels({'Laser OFF','Laser ON'})
ylabel('Vector length')
title('Forcing neg values to 0');

subplot(1,3,2); hold on
delta_rBn = r_bn(:,2) - r_bn(:,1);
delta_rBp = r_bp(:,2) - r_bp(:,1);
mean_deltaBn = nanmean(delta_rBn);
sem_deltaBn = std(delta_rBn,[],1,'omitnan')./sqrt(size(delta_rBn,1));
mean_deltaBp = nanmean(delta_rBp);
sem_deltaBp = std(delta_rBp,[],1,'omitnan')./sqrt(size(delta_rBp,1));
b2 = bar([mean_deltaBn mean_deltaBp],'FaceColor', 'flat');
drawnow
plot([1 1],[mean_deltaBn+sem_deltaBn mean_deltaBn-sem_deltaBn],'k')
plot([2 2],[mean_deltaBp+sem_deltaBp mean_deltaBp-sem_deltaBp],'k')
b2.CData(1,:) = color_bn;
b2.CData(2,:) = color_bp;
xticks([1 2]); xticklabels({'bead-','bead+'}); xtickangle(45)
ylabel('delta(Vector length)')
scatter(0.8*ones(size(delta_rBn,1),1),delta_rBn,...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
scatter(2.2*ones(size(delta_rBp,1),1),delta_rBp,...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
if isfield(data, 'beads_pos')
if adtest(delta_rBn) && adtest(delta_rBp)
   p = ranksum(delta_rBn,delta_rBp);
    testName = 'wilcoxon';
else
[~,p]=ttest2(delta_rBn,delta_rBp);
testName = 'ttest';
end
end
yl = ylim;
title([testName,', p = ' num2str(p)])
set(gcf,'Units','Normalized','Position',[0.4 0.4 0.4 0.3])

subplot(1,3,3); hold on
if looped_cells
r_all = cat(1,r_bn,r_bp); delta_rall = r_all(:,2)-r_all(:,1);
else
    delta_rall = r_bp(:,2)-r_bp(:,1);
end
violin(delta_rall,'xlabel',{'All Cells'},'facecolor',[0.5 .5 .5],'edgecolor','none');
plot(xlim,[0 0],'k--');
if adtest(delta_rall)
    p = signrank(delta_rall);
    testName = 'wilcoxon';
else
    [~,p] = ttest(delta_rall);
    testName = 'ttest';
end
ylim(yl)
title([testName,', p = ',num2str(p,3)])
% title('Forcing neg values to 0');
%%

%% 2.3.2 Vector length as a function of dFF change (laser ON vs laser OFF)
VectLengthVsAct_Fig = figure;
if looped_cells
subplot(2,3,1); hold on
scatter(deltaFR(hVis_cat&~bead_cat&dirSel),delta_rBn,...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
Bn_mdl= fitlm(deltaFR(hVis_cat&~bead_cat&dirSel),delta_rBn,'linear');
[Bn_pred, Bn_ci] = predict(Bn_mdl,deltaFR(hVis_cat&~bead_cat&dirSel));
plot(deltaFR(hVis_cat&~bead_cat&dirSel),Bn_pred,'-','color',color_bn);
plot(deltaFR(hVis_cat&~bead_cat&dirSel),Bn_ci, 'color',color_bn,'LineStyle',':');
ylabel('delta Vector Length')
title('Bead-')
xlabel('delta dFF, pref vis stim')
end

subplot(2,3,3); hold on
scatter(deltaFR(hVis_cat&bead_cat&dirSel),delta_rBp,...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
Bp_mdl= fitlm(deltaFR(hVis_cat&bead_cat&dirSel),delta_rBp,'linear');
[Bp_pred, Bp_ci] = predict(Bp_mdl,deltaFR(hVis_cat&bead_cat&dirSel));
plot(deltaFR(hVis_cat&bead_cat&dirSel),Bp_pred,'-','color',color_bp);
plot(deltaFR(hVis_cat&bead_cat&dirSel),Bp_ci, 'color',color_bp,'LineStyle',':');
title('Bead+')
xlabel('delta dFF, pref vis stim')

subplot(2,3,2); axis off
xl = xlim; yl =ylim;
if looped_cells
text(xl(1),yl(2),...
    sprintf('bead- \n slope = %.3g \n adj. r2 = %.3g',table2array(Bn_mdl.Coefficients(2,1)),Bn_mdl.Rsquared.Adjusted),...
    'HorizontalAlignment','left','VerticalAlignment','top');
end
text(xl(2),yl(1),...
    sprintf('bead+ \n slope = %.3g \n adj. r2 = %.3g',table2array(Bp_mdl.Coefficients(2,1)),Bp_mdl.Rsquared.Adjusted),...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
clear Bp_mdl Bp_pred Bp_ci
clear Bn_mdl Bn_pred Bn_ci

if looped_cells
subplot(2,3,4); hold on
scatter(ratioFR(hVis_cat&~bead_cat&dirSel),delta_rBn,...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
Bn_mdl= fitlm(ratioFR(hVis_cat&~bead_cat&dirSel),delta_rBn,'linear');
[Bn_pred, Bn_ci] = predict(Bn_mdl,deltaFR(hVis_cat&~bead_cat&dirSel));
plot(ratioFR(hVis_cat&~bead_cat&dirSel),Bn_pred,'-','color',color_bn);
plot(ratioFR(hVis_cat&~bead_cat&dirSel),Bn_ci, 'color',color_bn,'LineStyle',':');
ylabel('delta Vector Length')
xlabel('ratio dFF, pref vis stim')
end

subplot(2,3,6); hold on
scatter(ratioFR(hVis_cat&bead_cat&dirSel),delta_rBp,...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
Bp_mdl= fitlm(ratioFR(hVis_cat&bead_cat&dirSel),delta_rBp,'linear');
[Bp_pred, Bp_ci] = predict(Bp_mdl,ratioFR(hVis_cat&bead_cat&dirSel));
plot(ratioFR(hVis_cat&bead_cat&dirSel),Bp_pred,'-','color',color_bp);
plot(ratioFR(hVis_cat&bead_cat&dirSel),Bp_ci, 'color',color_bp,'LineStyle',':');
xlabel('ratio dFF, pref vis stim')

subplot(2,3,5); axis off
xl = xlim; yl =ylim;
if looped_cells
text(xl(1),yl(2),...
    sprintf('bead- \n slope = %.3g \n adj. r2 = %.3g',table2array(Bn_mdl.Coefficients(2,1)),Bn_mdl.Rsquared.Adjusted),...
    'HorizontalAlignment','left','VerticalAlignment','top');
end
text(xl(2),yl(1),...
    sprintf('bead+ \n slope = %.3g \n adj. r2 = %.3g',table2array(Bp_mdl.Coefficients(2,1)),Bp_mdl.Rsquared.Adjusted),...
    'HorizontalAlignment','right','VerticalAlignment','bottom');
clear Bp_mdl Bp_pred Bp_ci
clear Bn_mdl Bn_pred Bn_ci
%%

%% 2.3.3 Tuning width (sigma of the gaussian)
sigma_cat = squeeze(cat(1,sigma(:,1,:),sigma(:,2,:)));
sigma_bn = sigma_cat(hVis_cat&~bead_cat&dirSel,:);
sigma_bp = sigma_cat(hVis_cat&bead_cat&dirSel,:);
SigmaChange_bn = (sigma_bn(:,2) - sigma_bn(:,1))./(sigma_bn(:,2) + sigma_bn(:,1));
meanSCh_bn = mean(SigmaChange_bn);
semSCh_bn = std(SigmaChange_bn)./sqrt(sum(hVis_cat&~bead_cat&dirSel));
SigmaChange_bp = (sigma_bp(:,2) - sigma_bp(:,1))./(sigma_bp(:,2) + sigma_bp(:,1));
meanSCh_bp = mean(SigmaChange_bp);
semSCh_bp = std(SigmaChange_bp)./sqrt(sum(hVis_cat&bead_cat&dirSel));

Sigma_Fig = figure;
subplot(1,3,1); hold on
s1_bn = scatter(sigma_bn(:,1),sigma_bn(:,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
s2 = scatter(sigma_bp(:,1),sigma_bp(:,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
xlabel('Tuning width, Laser OFF')
ylabel('Tuning width, Laser ON')
xl = xlim; yl = ylim;
axis equal
ax_min = min(xl(1),yl(1)) ; ax_max = max(xl(2),yl(2));
xlim([ax_min ax_max]); ylim([ax_min ax_max])
plot([ax_min ax_max],[ax_min ax_max],'k--', 'LineWidth',0.5)
legend([s1_bn,s2],{'bead-','bead+'})
subplot(1,3,2); hold on
b = bar([meanSCh_bn; meanSCh_bp],'FaceColor','flat');
b.CData(1,:) = color_bn; b.CData(2,:) = color_bp;
xticks([1 2]); xticklabels({'bead-','bead+'}); xtickangle(45)
plot([1 1],[meanSCh_bn+semSCh_bn meanSCh_bn-semSCh_bn],'k')
plot([2 2],[meanSCh_bp+semSCh_bp meanSCh_bp-semSCh_bp],'k')
ylabel('Change in Tuning Width')

if isfield(data, 'beads_pos') && sum(SigmaChange_bp)>=4
if ~adtest(SigmaChange_bn) && ~adtest(SigmaChange_bp)
    [~,p]= ttest2(SigmaChange_bn,SigmaChange_bp);
    testName = 'ttest';
else
    [p,~]= ranksum(SigmaChange_bn,SigmaChange_bp);
    testName = 'Wilcoxon';
end
else
    testName = ''; p=[];
end
title([testName, ' p = ',num2str(p,3)])

subplot(1,3,3); hold on
SigmaChange_all = [SigmaChange_bn; SigmaChange_bp];
violin(SigmaChange_all,'xlabel',{'all cells'},'facecolor',[.5 .5 .5],'edgecolor', 'none');
plot(xlim,[0 0],'k--')
if isfield(data, 'beads_pos')
if adtest(SigmaChange_all)
   p = signrank(SigmaChange_all);
        testName = 'wilcoxon';
else
     [~,p] = ttest(SigmaChange_all);
    testName = 'ttest';
end
end
title([testName,', p = ',num2str(p,3)])

set(gcf,'Units','Normalized','Position',[0.1 0.4 0.7 0.3])

%% 2.4 Modeled inhibition on tuning curve

%% 2.4.1 Model Individual Cell
indivCell_mdl = cell(1,size(normCtrl_merged,1));
c_bp = 0; c_bn = 0;
for i = 1:size(normCtrl_merged,1)
    if hVis_cat(i) && dirSel(i)
        indivCell_mdl{i}= fitlm(normCtrl_merged(i,:,1),normCtrl_merged(i,:,2),'linear');
        %         [Bn_pred, Bn_ci] = predict(indivCell_mdl,normCtrl_merged(i,:,1));
        if bead_cat(i) == 1
            c_bp = c_bp+1;
            slope_bp(c_bp) = table2array(indivCell_mdl{i}.Coefficients(2,1));
        else
            c_bn = c_bn+1;
            slope_bn(c_bn) = table2array(indivCell_mdl{i}.Coefficients(2,1));
        end
    end
end

if looped_cells == 0
    slope_bn = slope_bp;
end

data_violin{1} = slope_bn;
data_violin{2} = slope_bp;
Slope_AllCells = figure;
subplot(1,2,1); hold on
% b = bar(1,mean(slope_bn)) ;bar(2, mean(slope_bp))
v = violin(data_violin,'xlabel',{'bead-','bead+'},'facecolor',color_bn,'edgecolor','none','FaceAlpha',1);
v(2).FaceColor = color_bp;
scatter(0.6*ones(size(slope_bn,2),1),slope_bn,...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
scatter(1.6*ones(size(slope_bp,2),1),slope_bp,...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
ylabel('slope')
plot(xlim,[1 1],'k--')
if isfield(data, 'beads_pos') && sum(slope_bp)>=4
if ~adtest(slope_bn) && ~adtest(slope_bp)
    [~,p]= ttest2(slope_bn,slope_bp);
    testName = 'ttest';
else
    [p,~]= ranksum(slope_bn,slope_bp);
    testName = 'Wilcoxon';
end
else
    testName='';
end
yl = ylim;
text(1,yl(2),['n = ',num2str(size(slope_bn,2))],...
    'HorizontalAlignment','center','VerticalAlignment','top')
text(2,yl(2),['n = ',num2str(size(slope_bp,2))],...
    'HorizontalAlignment','center','VerticalAlignment','top')
title(sprintf('Linear Fit for all cells \n%s p = %.3g',testName,p))

subplot(1,2,2); hold on
slope_all = [slope_bn slope_bp];
violin(slope_all','xlabel',{'all cells'},'facecolor',[.5 .5 .5],'edgecolor','none','FaceAlpha',1);
plot(xlim,[1 1],'k--')
if isfield(data, 'beads_pos')
if adtest(slope_all)
    p = signtest(slope_all,1);
    testName = 'wilcoxon';
else
      p = signtest(slope_all,1);
    testName = 'ttest';  
end
end
title([testName,', p = ' num2str(p,3)])

%%

%% 2.4.2 Model median response
tuning_bp = squeeze(median(normCtrl_merged(hVis_cat&bead_cat&dirSel,:,:)));
tuning_bn = squeeze(median(normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,:)));
[~,I_bn]=sort(tuning_bn(1:8,1));
[~,I_bp]=sort(tuning_bp(1:8,1));

ModelInhib_Figure1 = figure;
subplot(1,3,1); hold on
scatter(tuning_bn(I_bn,1),tuning_bn(I_bn,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',1)
xlabel('laser OFF'); ylabel('laser ON')
Bn_mdl= fitlm(tuning_bn(1:8,1),tuning_bn(1:8,2),'linear');
[Bn_pred, Bn_ci] = predict(Bn_mdl,tuning_bn(1:8,1));
plot(tuning_bn(1:8,1),Bn_pred,'-','color',color_bn);
plot(tuning_bn(1:8,1),Bn_ci, 'color',color_bn,'LineStyle',':');
axis equal
xl = xlim; yl = ylim;
plot([max(xl(1),yl(1)) min(xl(2),yl(2))],[max(xl(1),yl(1)) min(xl(2),yl(2))],'k-')
xl = xlim; yl = ylim;
text(xl(2),yl(2),{['slope: ' num2str(table2array(Bn_mdl.Coefficients(2,1)),3)],...
    ['adj. r2 = ' num2str(Bn_mdl.Rsquared.Adjusted,3)]},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
title(['bead-, n = ',num2str(sum(hVis_cat&~bead_cat&dirSel))])

subplot(1,3,3); hold on
scatter(tuning_bp(I_bp,1),tuning_bp(I_bp,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',1)
Bp_mdl= fitlm(tuning_bp(1:8,1),tuning_bp(1:8,2),'linear');
[Bp_pred, Bp_ci] = predict(Bp_mdl,tuning_bp(1:8,1));
plot(tuning_bp(1:8,1),Bp_pred,'-','color',color_bp);
plot(tuning_bp(1:8,1),Bp_ci, 'color',color_bp,'LineStyle',':');
axis equal
xl = xlim; yl = ylim;
plot([max(xl(1),yl(1)) min(xl(2),yl(2))],[max(xl(1),yl(1)) min(xl(2),yl(2))],'k-')
xl = xlim; yl = ylim;
text(xl(2),yl(2),{['slope: ' num2str(table2array(Bp_mdl.Coefficients(2,1)),3)],...
    ['adj. r2 = ' num2str(Bp_mdl.Rsquared.Adjusted,3)]},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
xlabel('laser OFF');
title(['bead+, n = ',num2str(sum(hVis_cat&bead_cat&dirSel))])

subplot(1,3,2); hold on
scatter(tuning_bn(I_bn,1),tuning_bn(I_bn,2),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',1)
scatter(tuning_bp(I_bp,1),tuning_bp(I_bp,2),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',1)
plot(tuning_bn(1:8,1),Bn_pred,'-','color',color_bn);
plot(tuning_bn(1:8,1),Bn_ci, 'color',color_bn,'LineStyle',':');
plot(tuning_bp(1:8,1),Bp_pred,'-','color',color_bp);
plot(tuning_bp(1:8,1),Bp_ci, 'color',color_bp,'LineStyle',':');
axis equal
xl = xlim; yl = ylim;
plot([max(xl(1),yl(1)) min(xl(2),yl(2))],[max(xl(1),yl(1)) min(xl(2),yl(2))],'k-')
xlabel('laser OFF');

%--------------------------- do stats ---------------------------
% i) z score
slope_bn = table2array(Bn_mdl.Coefficients(2,1));
slope_bp = table2array(Bp_mdl.Coefficients(2,1));
SEslope_bn = table2array(Bn_mdl.Coefficients(2,2));
SEslope_bp = table2array(Bp_mdl.Coefficients(2,2));
z_fischer = (slope_bn-slope_bp)/sqrt(SEslope_bn.^2+SEslope_bp.^2);
p = normcdf([-z_fischer z_fischer]);
pval = 1-(p(2)-p(1));

% ii) ANCOVA
% test1 = [tuning_bn;tuning_bp];
% group = [ones(9,1);2*ones(9,1)];
% [h,atab,ctab,stats] = aoctool(test1(:,1),test1(:,2),group);
% ancova_p = table2array(atab(4,6));
% title({'all',sprintf('z = %.3g; p = %.3g',z_fischer,pval),...
%     sprintf('ANCOVA p = %.3g',ancova_p)})
% sgtitle('median responses of the session to the 8 directions')

% iii) re-shuffle test
n_reshuffle = 300;
tuning_all = normCtrl_merged(hVis_cat&dirSel,:,:);
slope_diff = NaN(1,n_reshuffle);
SlopeShuffle_Fig = figure;
subplot(1,2,1); hold on
for i = 1:n_reshuffle
    reSamp_bp = datasample(tuning_all,sum(hVis_cat&bead_cat&dirSel),1,'Replace',false);
    reSamp_bn = datasample(tuning_all,sum(hVis_cat&~bead_cat&dirSel),1,'Replace',false);
    median_bn = squeeze(median(reSamp_bn));
    median_bp = squeeze(median(reSamp_bp));
    % fit_bn = regress(median_bn(:,2),median_bn(:,1));
    % fit_bp = regress(median_bp(:,2),median_bp(:,1));
    fit_bn = regress(median_bn(:,2),[ones(size(median_bn(:,1))) median_bn(:,1)]);
    fit_bp = regress(median_bp(:,2),[ones(size(median_bp(:,1))) median_bp(:,1)]);
    
    % plot(median_bn(:,1),median_bn(:,2),'bo')
    p1 = plot(median_bn(:,1),fit_bn(1)+fit_bn(2)*median_bn(:,1),'b-');
    p1.Color = [color_bn 0.1];
    % plot(median_bp(:,1),median_bp(:,2),'ro')
    p2 = plot(median_bp(:,1),fit_bp(1)+fit_bp(2)*median_bp(:,1),'r-');
    p2.Color = [color_bp 0.1];
    slope_diff(i) = fit_bp(2)-fit_bn(2);
end
xlabel('norm. Response Magn, Laser ON')
ylabel('norm. Response Magn, Laser OFF')
p4 = plot(tuning_bp(1:8,1),Bp_pred,'r-');
p3 = plot(tuning_bn(1:8,1),Bn_pred,'b-');
legend([p1,p2,p3,p4],{'bead-, shuffle','bead+, shuffle','bead-, data','bead+, data'},...
    'Location','Northwest')
legend('boxoff')
title(['Linear Fit | ',num2str(n_reshuffle),' resampling'])
[fi,xi]=ksdensity(slope_diff);
subplot(1,2,2); hold on
plot(xi,fi,'k');
xlabel('delta slope, {\it bead^+ - bead^-}')
ylabel('probability density estimate')
yyaxis right
histogram(slope_diff)
ylabel('# of observations')
plot([slope_bp-slope_bn slope_bp-slope_bn],ylim,'r')
legend({'p. density','data'},'Location','Southeast')
legend('boxoff')
interval = xi>slope_bp-slope_bn;
twosided_prob =(1-trapz(fi(interval))/trapz(fi))*2;
yl = ylim;
text(slope_bp-slope_bn,yl(2),sprintf('%.3g %% values \n below or above\n obs. value',twosided_prob*100),...
    'HorizontalAlignment','right','VerticalAlignment','top','Fontsize',8)
title('Difference in slopes')
sgtitle('Difference in slope, shuffle')

% -------- pour boubou ---------------------
figure;
subplot(1,2,1); hold on
plot(log10(tuning_bn(:,1)),'k')
plot(log10(tuning_bn(:,2)),'o','Color',color_bn)
pred_bn = tuning_bn(:,1)*table2array(Bn_mdl.Coefficients(2,1))+table2array(Bn_mdl.Coefficients(1,1));
plot(log10(pred_bn),'-','Color',color_bn)
ylabel('log response')
xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
title('Control')
subplot(1,2,2); hold on
plot(log10(tuning_bp(:,1)),'ko-')
plot(log10(tuning_bp(:,2)),'o','Color',color_bp)
pred_bp = tuning_bp(:,1)*table2array(Bp_mdl.Coefficients(2,1))+table2array(Bp_mdl.Coefficients(1,1));
pred_bp(pred_bp<0) = 0.001;
plot(log10(pred_bp),'-','Color',color_bp)
xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
title('Traitement')

% figure; hold on
% plot(log(tuning_bn(:,2))-log(tuning_bn(:,1)))
% plot(log(tuning_bp(:,2))-log(tuning_bp(:,1)),'Color',color_bp)
% ylabel('delta(log response)')
% xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
% legend({'Control','Traitement'})
% --------------------------- ---------- ---------------------------

% % Test Gildas
% for i = 1:8
% ecart(i) = (tuning_bp(i,1)-tuning_bp(i,2))-(tuning_bn(i,1)-tuning_bn(i,2));
% end
% ecart_sum = sum(ecart(:));
% var_x = std(0:45:315);
% var_sum = std(tuning_bp(:,1)).^2 + std(tuning_bp(:,2)).^2 + std(tuning_bn(:,1)).^2 + std(tuning_bn(:,2)).^2 ;
% z_boubou = (ecart_sum*sqrt(8)*var_x)/sqrt(var_sum);
% 
% p = normcdf([-z_boubou z_boubou]);
% pval = 1-(p(2)-p(1));


% fit_sel_Bpos_all = MF_cat(hVis_cat&k_cat&bead_cat&~DS,:);
% fit_sel_Bneg_all = MF_cat(hVis_cat&k_cat&~bead_cat&~DS,:);
% two_inhib = fit_sel_Bpos_all(:,1)<0 & fit_sel_Bpos_all(:,2)<0;
% two_inhib_neg = fit_sel_Bneg_all(:,1)<0 & fit_sel_Bneg_all(:,2)<0;
% % two_inhib = true(size(fit_sel_Bpos_all(:,1)));
% % two_inhib_neg = true(size(fit_sel_Bneg_all(:,1)));
% index_bp = (fit_sel_Bpos_all(two_inhib,1)-fit_sel_Bpos_all(two_inhib,2))./(fit_sel_Bpos_all(two_inhib,1)+fit_sel_Bpos_all(two_inhib,2));
% index_bn = (fit_sel_Bneg_all(two_inhib_neg,1)-fit_sel_Bneg_all(two_inhib_neg,2))./(fit_sel_Bneg_all(two_inhib_neg,2)+fit_sel_Bneg_all(two_inhib_neg,1));
% % index_bp = fit_sel_Bpos_all(two_inhib,2)./fit_sel_Bpos_all(two_inhib,1);
% % index_bn = fit_sel_Bneg_all(two_inhib_neg,2)./fit_sel_Bneg_all(two_inhib_neg,1);
% violin_data{1} = index_bn;
% violin_data{2} = index_bp;
% figure; hold on
% v = violin(violin_data);
% v(1).FaceColor = color_bn; v(2).FaceColor = color_bp;
% % b=bar([mean(index_bn) mean(index_bp)],'FaceColor','flat');
% % b.CData(2,:) = color_bp;
% % sem_bn = std(index_bn)./sqrt(length(index_bn));
% % sem_bp = std(index_bp)./sqrt(length(index_bp));
% % plot([1 1],[mean(index_bn)-sem_bn mean(index_bn)+sem_bn],'k-')
% % plot([2 2],[mean(index_bp)-sem_bp mean(index_bp)+sem_bp],'k-')
% xticks([1 2]); xticklabels({'bead-','bead+'})
% ylabel('Ratio of the MF for the 2 Gaussians')
% p = ranksum(index_bn,index_bp);
%     title({'OS cells & laser-induced inhibition on BOTH gaussian',['Wilcoxon, p = ',num2str(p,3)]})
% % yyaxis right
% scatter(0.8*ones(size(index_bn)),index_bn,...
%     'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
% scatter(2.2*ones(size(index_bp)),index_bp,...
%     'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
% plot(xlim,[0 0],'k:')
% text(1,1,sprintf('n = %g',length(index_bn)),...
%     'HorizontalAlignment','center','VerticalAlignment','top')
% text(2,1,sprintf('n = %g',length(index_bp)),...
%         'HorizontalAlignment','center','VerticalAlignment','top')
% % figure; hold on
% % scatter(fit_sel_Bpos_all(:,1),fit_sel_Bpos_all(:,2),...
% %     'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
% % scatter(fit_sel_Bneg_all(:,1),fit_sel_Bneg_all(:,2),...
% %     'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
% % plot(xlim,[0 0],'k-')
% % plot([0 0 ],ylim,'k-')
% 
% % sampsizepwr('z',[mean(index_bp) std(index_bp)],mean(index_bn))
% --------------------------- ---------- ---------------------------

% ---------------------- fit data from all cells --------------------
% figure; hold on
% s1_bn=size(normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,1),1);
% plot(normCtrl_merged(hVis_cat&~bead_cat&dirSel,1:8,1),normCtrl_merged(hVis_cat&~bead_cat&dirSel,1:8,2),'o','Color',color_bn);
% test = reshape(normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,1),s1_bn*9,1);
% test2 = reshape(normCtrl_merged(hVis_cat&~bead_cat&dirSel,:,2),s1_bn*9,1);
% Bn_mdl_indiv= fitlm(test,test2,'linear');
% [Bn_pred_indiv, Bn_ci_indiv] = predict(Bn_mdl_indiv,test);
% plot(test,Bn_pred_indiv,'-','color',color_bn);
% plot(test,Bn_ci_indiv, 'color',color_bn,'LineStyle',':');
% 
% s1_bp=size(normCtrl_merged(hVis_cat&bead_cat&dirSel,:,1),1);
% plot(normCtrl_merged(hVis_cat&bead_cat&dirSel,1:8,1),normCtrl_merged(hVis_cat&bead_cat&dirSel,1:8,2),'o','Color',color_bp);
% test = reshape(normCtrl_merged(hVis_cat&bead_cat&dirSel,:,1),s1_bp*9,1);
% test2 = reshape(normCtrl_merged(hVis_cat&bead_cat&dirSel,:,2),s1_bp*9,1);
% Bp_mdl_indiv= fitlm(test,test2,'linear');
% [Bp_pred_indiv, Bp_ci_indiv] = predict(Bp_mdl_indiv,test);
% plot(test,Bp_pred,_indiv'-','color',color_bp);
% plot(test,Bp_ci_indiv, 'color',color_bp,'LineStyle',':');
% --------------------------- ---------- ---------------------------

ModelInhib_Figure2 = figure;
subplot(2,3,1); hold on
p1 = plot(tuning_bn(:,1),'ko-','MarkerFaceColor','k');
p2 = plot(tuning_bn(:,2),'o','MarkerFaceColor',color_bn,'MarkerEdgeColor',color_bn);
pred_bn = tuning_bn(:,1)*table2array(Bn_mdl.Coefficients(2,1))+table2array(Bn_mdl.Coefficients(1,1));
p3 = plot(pred_bn,'-','Color',color_bn);
for i = 1:9
    plot([i i],[tuning_bn(i,1)+tuning_bn_sem(i,1) tuning_bn(i,1)-tuning_bn_sem(i,1)],'k-')
    plot([i i],[tuning_bn(i,2)+tuning_bn_sem(i,2) tuning_bn(i,2)-tuning_bn_sem(i,2)],'-','Color',color_bn)
end
xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
ylabel('Resp (norm.)')
lg1 = legend([p1,p2,p3],{'laser OFF','Laser ON, data','Laser ON, divisive model'});
title('Bead-')

subplot(2,3,3); hold on
p4 = plot(tuning_bp(:,1),'ko-','MarkerFaceColor','k');
p5 = plot(tuning_bp(:,2),'o','MarkerFaceColor',color_bp,'MarkerEdgeColor',color_bp);
pred_bp = tuning_bp(:,1)*table2array(Bp_mdl.Coefficients(2,1))+table2array(Bp_mdl.Coefficients(1,1));
p6 = plot(pred_bp,'-','Color',color_bp);
for i = 1:9
    plot([i i],[tuning_bp(i,1)+tuning_bp_sem(i,1) tuning_bp(i,1)-tuning_bp_sem(i,1)],'k-')
    plot([i i],[tuning_bp(i,2)+tuning_bp_sem(i,2) tuning_bp(i,2)-tuning_bp_sem(i,2)],'-','Color',color_bp)
end
xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'})
lg2 = legend([p4,p5,p6],{'laser OFF','Laser ON, data','Laser ON, divisive model'});
title('Bead+')

subplot(2,3,2); axis off
lg1.Position = [0.35 0.8 0.1 0.1];
lg2.Position = [0.55 0.6 0.1 0.1];
lg1.ItemTokenSize = [20,4];
lg2.ItemTokenSize = [20,4];

subplot(2,3,4); hold on
ImpactOfLaser_sub = tuning_bn(:,2) - tuning_bn(:,1);
ImpactOfLaser_div = tuning_bn(:,2)./ tuning_bn(:,1);
s1_bn = scatter(tuning_bn(:,1),ImpactOfLaser_sub,'MarkerFaceColor',color_bn,'MarkerEdgeColor','none');
s2 = scatter(tuning_bn(:,1),ImpactOfLaser_div,'MarkerFaceColor','none','MarkerEdgeColor',color_bn);

BnSub_mdl= fitlm(tuning_bn(:,1),ImpactOfLaser_sub,'linear');
[BnSub_pred, BnSub_ci] = predict(BnSub_mdl,tuning_bn(:,1));
plot(tuning_bn(:,1),BnSub_pred,'-','color',color_bn);
plot(tuning_bn(:,1),BnSub_ci, 'color',color_bn,'LineStyle',':');

BnDiv_mdl= fitlm(tuning_bn(:,1),ImpactOfLaser_div,'linear');
[BnDiv_pred, BnDiv_ci] = predict(BnDiv_mdl,tuning_bn(:,1));
plot(tuning_bn(:,1),BnDiv_pred,'-','color',color_bn);
plot(tuning_bn(:,1),BnDiv_ci, 'color',color_bn,'LineStyle',':');

lg1 = legend({'Subtractive','Divisive'},'Location','southwest');
xlabel('Response to Vis Stim (norm.)')
ylabel({'Impact of Laser'})

subplot(2,3,6); hold on
ImpactOfLaser_sub = tuning_bp(:,2) - tuning_bp(:,1);
ImpactOfLaser_div = tuning_bp(:,2)./ tuning_bp(:,1);
s1_bn = scatter(tuning_bp(:,1),ImpactOfLaser_sub,'MarkerFaceColor',color_bp,'MarkerEdgeColor','none');
s2 = scatter(tuning_bp(:,1),ImpactOfLaser_div,'MarkerFaceColor','none','MarkerEdgeColor',color_bp);

BpSub_mdl= fitlm(tuning_bp(:,1),ImpactOfLaser_sub,'linear');
[BpSub_pred, BpSub_ci] = predict(BpSub_mdl,tuning_bp(:,1));
plot(tuning_bp(:,1),BpSub_pred,'-','color',color_bp);
plot(tuning_bp(:,1),BpSub_ci, 'color',color_bp,'LineStyle',':');

BpDiv_mdl= fitlm(tuning_bp(:,1),ImpactOfLaser_div,'linear');
[BpDiv_pred, BpDiv_ci] = predict(BpDiv_mdl,tuning_bp(:,1));
plot(tuning_bp(:,1),BpDiv_pred,'-','color',color_bp);
plot(tuning_bp(:,1),BpDiv_ci, 'color',color_bp,'LineStyle',':');

lg2 = legend({'Subtractive','Divisive'},'Location','southwest');
xlabel('Response to Vis Stim (norm.)')
title('Impact of Laser')

subplot(2,3,5); axis off
xl = xlim; yl = ylim;
text(xl(1),yl(2),{'bead-',['subtractive, slope = ', num2str(table2array(BnSub_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BnSub_mdl.Rsquared.Adjusted,3)],...
    ['divisive, slope = ', num2str(table2array(BnDiv_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BnDiv_mdl.Rsquared.Adjusted,3)]},...
    'HorizontalAlignment','Left','VerticalAlignment','Top')
text(xl(2),yl(1),{'bead+',['subtractive, slope = ', num2str(table2array(BpSub_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BpSub_mdl.Rsquared.Adjusted,3)],...
    ['divisive, slope = ', num2str(table2array(BpDiv_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BpDiv_mdl.Rsquared.Adjusted,3)]},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
lg1.Position = [0.35 0.22 0.1 0.1];
lg2.Position = [0.55 0.12 0.1 0.1];
lg1.ItemTokenSize = [20,4];
lg2.ItemTokenSize = [20,4];

sgtitle('Modeled inhibition')
set(gcf,'Units','Normalized','Position',[0.05 0.4 0.7 0.5])
%%

%% 2.5 DSI
%  2.5.1 --------- calculate shift in pref direction ----------------------
% force negative values for antiprefered to 0
k1 = dir_tuning(:,1,1,1) < 0;
dir_tuning(k1,1,1,1) = 0;
k2 = dir_tuning(:,1,2,1) < 0;
dir_tuning(k2,1,2,1) = 0;

DSI(:,:,1) = (dir_tuning(:,5,:,1) - dir_tuning(:,1,:,1)) ./(dir_tuning(:,5,:,1) + dir_tuning(:,1,:,1));
[~,I] = max(dir_tuning(:,:,1,2),[],2);
Iopto(:,1) = I;
[~,I] = max(dir_tuning(:,:,2,2),[],2);
Iopto(:,2) = I;
PrefShift = [5,5]-Iopto;
PrefShift_cat = cat(1,PrefShift(:,1),PrefShift(:,2));

%% ----------------- delta(pref direction) for all cells ------------------
% regardless if they are OS or DS

% DSI_Figure1 = figure;
% subplot(3,2,1); hold on
% histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel),'FaceColor',color_bn);
% ylabel('# of cells');
% xlim([-3.5 4.5]); xticks([]);
% subplot(3,2,2); hold on
% histogram(abs(PrefShift_cat(hVis_cat&~bead_cat&dirSel)),'FaceColor',color_bn);
% xticks([]);
% subplot(3,2,3); hold on
% h1 = histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel),'FaceColor',color_bn,'Normalization','probability');
% h2 = histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel),'FaceColor',color_bp,'Normalization','probability');
% yl = ylim;
% text(1,yl(2),{['bead-: ' num2str(max(h1.Values),2)],...
%     ['bead+: ' num2str(max(h2.Values),2)]},...
%     'HorizontalAlignment','left','VerticalAlignment','top')
% ylabel('Probability');
% xlim([-3.5 4.5]); xticks([]);
% subplot(3,2,4); hold on
% histogram(abs(PrefShift_cat(hVis_cat&~bead_cat&dirSel)),'FaceColor',color_bn,'Normalization','probability');
% histogram(abs(PrefShift_cat(hVis_cat&bead_cat&dirSel)),'FaceColor',color_bp,'Normalization','probability');
% xticks([]);
% subplot(3,2,5); hold on
% histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel),'FaceColor',color_bp);
% ylabel('# of cells');
% xlim([-3.5 4.5]);
% xticks(-3:4);xticklabels({'-135','-90','-45','0','45','90','135','180'});
% xtickangle(45)
% xlabel('Shift in Pref Direction (deg)')
% subplot(3,2,6); hold on
% histogram(abs(PrefShift_cat(hVis_cat&bead_cat&dirSel)),'FaceColor',color_bp);
% xticks(0:4);xticklabels({'0','45','90','135','180'});
% xtickangle(45)
% xlabel('Abs Shift in Pref Direction (deg)')
% sgtitle('Change in pref direction, cells selected on control stim')
% %  chi2 test
% nBnCells = size(PrefShift_cat(hVis_cat&~bead_cat&dirSel),1);
% nBpCells = size(PrefShift_cat(hVis_cat&bead_cat&dirSel),1);
% dataID = [repmat('a',nBnCells,1);repmat('b',nBpCells,1)];
% nPerCat_bn = zeros(1,7); nPerCat_bp = zeros(1,7);
% for i = -3:4
%     nPerCat_bn(i+4) = length(find(PrefShift_cat(hVis_cat&~bead_cat&dirSel)==i));
%     nPerCat_bp(i+4) = length(find(PrefShift_cat(hVis_cat&bead_cat&dirSel)==i));
% end
%
% categories_bn = []; categories_bp = [];
% for i = 1:8
%     categories_bn = [categories_bn; repmat(i,nPerCat_bn(i),1)];
%     categories_bp = [categories_bp; repmat(i,nPerCat_bp(i),1)];
% end
% categories = [categories_bn;categories_bp];
% [~,chi2stat,pval] = crosstab(dataID,categories);
% subplot(3,2,3);
% title(sprintf('chi2 test, chi2stat = %.3g; p = %.3g',chi2stat ,pval))
%% ------------------------------------------------------------------------
%% 2.5.2 ------- calculate DSI using new pref direction for LASER ON ------
Iantipref = Iopto-4;
Iantipref = mod(Iantipref,8);
Iantipref(Iantipref == 0)= 8;

for i = 1:size(dir_tuning,1)
    for ii = 1:nSF
        if dir_tuning(i,Iantipref(i,ii),ii,2) <0
            dir_tuning(i,Iantipref(i,ii),ii,2) = 0;
        end
        DSI_temp(i,ii) = (dir_tuning(i,Iopto(i,ii),ii,2) - dir_tuning(i,Iantipref(i,ii),ii,2))./(dir_tuning(i,Iopto(i,ii),ii,2) + dir_tuning(i,Iantipref(i,ii),ii,2));
    end
end
DSI(:,:,2) = DSI_temp; clear DSI_temp;
%     DSI(:,:,2) = (y(:,3,:,2) - y(:,7,:,2)) ./(y(:,3,:,2) + y(:,7,:,2));

DSI_cat = squeeze(cat(1,DSI(:,1,:),DSI(:,2,:)));
DSI_bn = mean(DSI_cat(hVis_cat&~bead_cat&dirSel,:));
DSI_bn_sem = std(DSI_cat(hVis_cat&~bead_cat&dirSel,:))./sqrt(sum(hVis_cat&~bead_cat&dirSel));
DSI_bp = nanmean(DSI_cat(hVis_cat&bead_cat&dirSel,:));
DSI_bp_sem = std(DSI_cat(hVis_cat&bead_cat&dirSel,:))./sqrt(sum(hVis_cat&bead_cat&dirSel));

DSI_Figure2 = figure;
subplot(2,4,1); hold on
b1 = bar([DSI_bn' DSI_bp'],'FaceColor','flat');
b1(1).CData = color_bn;
b1(2).CData = color_bp;
drawnow
for i = 1:2
    plot([i+b1(1).XOffset i+b1(1).XOffset],[DSI_bn(i)+DSI_bn_sem(i) DSI_bn(i)-DSI_bn_sem(i)],'k')
    plot([i+b1(2).XOffset i+b1(2).XOffset],[DSI_bp(i)+DSI_bp_sem(i) DSI_bp(i)-DSI_bp_sem(i)],'k')
end
xticks([1 2]); xticklabels({'laserOFF','laserON'})
ylabel('DSI')

% subplot(2,3,2); hold on
for i = 1:2
    % plot([i-0.2 i-0.1],[DSI_bn(i) DSI_bn(i)],'k', 'LineWidth',2)
    % plot([i+0.1 i+0.2],[DSI_bp(i) DSI_bp(i)],'k', 'LineWidth',2)
    scatter(i*ones(sum(hVis_cat&~bead_cat&dirSel),1)+2*b1(1).XOffset,DSI_cat(hVis_cat&~bead_cat&dirSel,i),...
        'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.05)
    scatter(i*ones(sum(hVis_cat&bead_cat&dirSel),1)+2*b1(2).XOffset,DSI_cat(hVis_cat&bead_cat&dirSel,i),...
        'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
end
xticks([1 2]); xticklabels({'laserOFF','laserON'})
xlim([0 3])

deltaDSI = DSI_cat(:,2) - DSI_cat(:,1);
dDSI_bn = mean(deltaDSI(hVis_cat&~bead_cat&dirSel));
dDSI_bp = mean(deltaDSI(hVis_cat&bead_cat&dirSel));
dDSI_bn_sem = std(deltaDSI(hVis_cat&~bead_cat&dirSel))/sqrt(sum(hVis_cat&~bead_cat&dirSel));
dDSI_bp_sem = std(deltaDSI(hVis_cat&bead_cat&dirSel))/sqrt(sum(hVis_cat&bead_cat&dirSel));

if looped_cells
if isfield(data, 'beads_pos') && sum(hVis_cat&bead_cat&dirSel)>=4
    if ~adtest(deltaDSI(hVis_cat&~bead_cat&dirSel)) && ~adtest(deltaDSI(hVis_cat&bead_cat&dirSel))
        [~,p] = ttest2(deltaDSI(hVis_cat&~bead_cat&dirSel),deltaDSI(hVis_cat&bead_cat&dirSel));
        testName = 'ttest';
    else
        p = ranksum(deltaDSI(hVis_cat&~bead_cat&dirSel),deltaDSI(hVis_cat&bead_cat&dirSel));
        testName = 'wilcoxon';
    end
else
    testName = '';
end
else
    testName = ''; p= [];
end

subplot(2,4,3); hold on
b3 = bar([dDSI_bn dDSI_bp]','FaceColor','flat');
b3.CData(1,:) = color_bn;
b3.CData(2,:) = color_bp;
plot([1 1],[dDSI_bn+dDSI_bn_sem dDSI_bn-dDSI_bn_sem],'k-')
plot([2 2],[dDSI_bp+dDSI_bp_sem dDSI_bp-dDSI_bp_sem],'k-')
ylabel('delta DSI')
xticks([1 2]); xticklabels({'Bead-','Bead+'})
xtickangle(45)
title([testName, ', p = ' num2str(p)]);
% ------------------------- ------------ ----------------------------------

% deltaFR = squeeze(test_merged(:,3,2) - test_merged(:,3,1));
% ratioFR = squeeze(test_merged(:,3,2) ./ test_merged(:,3,1));

subplot(2,4,5); hold on
s1_bn = scatter(deltaFR(hVis_cat&~bead_cat&dirSel),deltaDSI(hVis_cat&~bead_cat&dirSel),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
s2 = scatter(deltaFR(hVis_cat&bead_cat&dirSel),deltaDSI(hVis_cat&bead_cat&dirSel),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
xlabel('delta dFF of pref vis stim')
ylabel('delta DSI')

if looped_cells
BnDSIsub_mdl= fitlm(deltaFR(hVis_cat&~bead_cat&dirSel),deltaDSI(hVis_cat&~bead_cat&dirSel),'linear');
[BnDSIsub_pred, BnDSIsub_ci] = predict(BnDSIsub_mdl,deltaFR(hVis_cat&~bead_cat&dirSel));
plot(deltaFR(hVis_cat&~bead_cat&dirSel),BnDSIsub_pred,'-','color',color_bn);
plot(deltaFR(hVis_cat&~bead_cat&dirSel),BnDSIsub_ci, 'color',color_bn,'LineStyle',':');
end

BpDSIsub_mdl= fitlm(deltaFR(hVis_cat&bead_cat&dirSel),deltaDSI(hVis_cat&bead_cat&dirSel),'linear');
[BpDSIsub_pred, BpDSIsub_ci] = predict(BpDSIsub_mdl,deltaFR(hVis_cat&bead_cat&dirSel));
plot(deltaFR(hVis_cat&bead_cat&dirSel),BpDSIsub_pred,'-','color',color_bp);
plot(deltaFR(hVis_cat&bead_cat&dirSel),BpDSIsub_ci, 'color',color_bp,'LineStyle',':');

subplot(2,4,7); hold on
s1_bn = scatter(ratioFR(hVis_cat&~bead_cat&dirSel),deltaDSI(hVis_cat&~bead_cat&dirSel),...
    'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
s2 = scatter(ratioFR(hVis_cat&bead_cat&dirSel),deltaDSI(hVis_cat&bead_cat&dirSel),...
    'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
xlabel('ratio dFF of pref vis stim')

if looped_cells
BnDSIdiv_mdl= fitlm(ratioFR(hVis_cat&~bead_cat&dirSel),deltaDSI(hVis_cat&~bead_cat&dirSel),'linear');
[BnDSIdiv_pred, BnDSIdiv_ci] = predict(BnDSIdiv_mdl,ratioFR(hVis_cat&~bead_cat&dirSel));
plot(ratioFR(hVis_cat&~bead_cat&dirSel),BnDSIdiv_pred,'-','color',color_bn);
plot(ratioFR(hVis_cat&~bead_cat&dirSel),BnDSIdiv_ci, 'color',color_bn,'LineStyle',':');
end

BpDSIdiv_mdl= fitlm(ratioFR(hVis_cat&bead_cat&dirSel),deltaDSI(hVis_cat&bead_cat&dirSel),'linear');
[BpDSIdiv_pred, BpDSIdiv_ci] = predict(BpDSIdiv_mdl,ratioFR(hVis_cat&bead_cat&dirSel));
plot(ratioFR(hVis_cat&bead_cat&dirSel),BpDSIdiv_pred,'-','color',color_bp);
plot(ratioFR(hVis_cat&bead_cat&dirSel),BpDSIdiv_ci, 'color',color_bp,'LineStyle',':');

subplot(2,4,6); axis off
xl = xlim; yl = ylim;
if looped_cells
text(xl(1),yl(2),{['bead- delta FR, slope = ', num2str(table2array(BnDSIsub_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BnDSIsub_mdl.Rsquared.Adjusted,3)],...
    ['bead+ delta FR, slope = ', num2str(table2array(BpDSIsub_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BpDSIsub_mdl.Rsquared.Adjusted,3)]},...
    'HorizontalAlignment','Left','VerticalAlignment','Top')
text(xl(2),yl(1),{['bead- ratio FR, slope = ', num2str(table2array(BnDSIdiv_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BnDSIdiv_mdl.Rsquared.Adjusted,3)],...
    ['bead+ ratio FR, slope = ', num2str(table2array(BpDSIdiv_mdl.Coefficients(2,1)),3),'; adj. r2 = ' num2str(BpDSIdiv_mdl.Rsquared.Adjusted,3)]},...
    'HorizontalAlignment','Right','VerticalAlignment','Top')
set(gcf,'Units','Normalized','Position',[0.05 0.4 0.7 0.5])
end

subplot(2,4,4); hold on
dDSI_all = deltaDSI(hVis_cat&dirSel);
violin(dDSI_all,'xlabel',{'all cells'},'facecolor',[.5 .5 .5],'edgecolor','none');
plot(xlim,[0 0],'k--')
if isfield(data, 'beads_pos') 
if adtest(dDSI_all)
    p = signtest(dDSI_all);
    testName= 'wilcoxon';
else
       p = ttest(dDSI_all);
    testName= 'ttest';    
end
end
title([testName,', p = ',num2str(p,3)]);
sgtitle('DSI');
%%

%% 2.5.3 ------ Separate cells in OS vs DS, pref dir shift ----------------
% dirSel is already selecting for cells with tuning WITH AND WITHOUT laser
DS = false(size(PrefShift_cat));
DS(DSI_cat(:,1)>0.35) = true;
% figure; hold on
% histogram(DSI_cat(hVis_cat&~bead_cat&dirSel&DS))
% histogram(DSI_cat(hVis_cat&~bead_cat&dirSel&~DS))

DSI_Figure3 = figure;
subplot(3,2,1); hold on
histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel&DS),'FaceColor',color_bn);
ylabel('# of cells');
xlim([-3.5 4.5]); xticks([]);
xl = xlim; yl = ylim;
title('DS cells')
text(xl(2),yl(2),['bead-, n = ', num2str(sum(hVis_cat&~bead_cat&dirSel&DS))],...
'HorizontalAlignment','right','VerticalAlignment','top','Fontsize',8)
subplot(3,2,2); hold on
histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel&~DS),'FaceColor',color_bn);
xlim([-3.5 4.5]); xticks([]);
xl = xlim; yl = ylim;
title('OS cells')
text(xl(2),yl(2),['bead-, n = ', num2str(sum(hVis_cat&~bead_cat&dirSel&~DS))],...
'HorizontalAlignment','right','VerticalAlignment','top','Fontsize',8)
subplot(3,2,3); hold on
h1 = histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel&DS),'FaceColor',color_bn,'Normalization','probability');
h2 = histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel&DS),'FaceColor',color_bp,'Normalization','probability');
yl = ylim;
text(1,yl(2),{['bead-: ' num2str(max(h1.Values),2)],...
    ['bead+: ' num2str(max(h2.Values),2)]},...
    'HorizontalAlignment','left','VerticalAlignment','top')
ylabel('Probability');
xlim([-3.5 4.5]); xticks([]);
subplot(3,2,4); hold on
h1 = histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel&~DS),'FaceColor',color_bn,'Normalization','probability');
h2 = histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel&~DS),'FaceColor',color_bp,'Normalization','probability');
yl = ylim;
text(1,yl(2),{['bead-: ' num2str(max(h1.Values),2)],...
    ['bead+: ' num2str(max(h2.Values),2)]},...
    'HorizontalAlignment','left','VerticalAlignment','top')
xlim([-3.5 4.5]); xticks([])
subplot(3,2,5); hold on
histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel&DS),'FaceColor',color_bp);
ylabel('# of cells');
xlim([-3.5 4.5]);
xticks(-3:4);xticklabels({'-135','-90','-45','0','45','90','135','180'});
xtickangle(45)
xlabel('Shift in Pref Direction (deg)')
xl = xlim; yl = ylim;
text(xl(2),yl(2),['bead+, n = ', num2str(sum(hVis_cat&bead_cat&dirSel&DS))],...
'HorizontalAlignment','right','VerticalAlignment','top','Fontsize',8)
subplot(3,2,6); hold on
histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel&~DS),'FaceColor',color_bp);
xlim([-3.5 4.5]);
xticks(-3:4);xticklabels({'-135','-90','-45','0','45','90','135','180'});
xtickangle(45)
xlabel('Shift in Pref Direction (deg)')
xl = xlim; yl = ylim;
text(xl(2),yl(2),['bead+, n = ', num2str(sum(hVis_cat&bead_cat&dirSel&~DS))],...
'HorizontalAlignment','right','VerticalAlignment','top','Fontsize',8)

% ----------------- plot all cells on the same graph
% subplot(3,3,3); hold on
% histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel),'FaceColor',color_bn);
% xlim([-3.5 4.5]); xticks([]);
% title('All tuned cells')
% subplot(3,3,6); hold on
% h1 = histogram(PrefShift_cat(hVis_cat&~bead_cat&dirSel),'FaceColor',color_bn,'Normalization','probability');
% h2 = histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel),'FaceColor',color_bp,'Normalization','probability');
% yl = ylim;
% text(1,yl(2),{['bead-: ' num2str(max(h1.Values),2)],...
%     ['bead+: ' num2str(max(h2.Values),2)]},...
%     'HorizontalAlignment','left','VerticalAlignment','top')
% xlim([-3.5 4.5]); xticks([]);
% subplot(3,3,9); hold on
% histogram(PrefShift_cat(hVis_cat&bead_cat&dirSel),'FaceColor',color_bp);
% xlim([-3.5 4.5]);
% xticks(-3:4);xticklabels({'-135','-90','-45','0','45','90','135','180'});
% xtickangle(45)
% xlabel('Shift in Pref Direction (deg)')
sgtitle('Change in pref direction')

% % chi2 test for DS cells
nDSnBn = size(PrefShift_cat(hVis_cat&~bead_cat&dirSel&DS),1);
nDSnBp = size(PrefShift_cat(hVis_cat&bead_cat&dirSel&DS),1);
dataID = [repmat('a',nDSnBn,1);repmat('b',nDSnBp,1)];
nPerCat_bn = zeros(1,7); nPerCat_bp = zeros(1,7);
for i = -3:4
    nPerCat_bn(i+4) = length(find(PrefShift_cat(hVis_cat&~bead_cat&dirSel&DS)==i));
    nPerCat_bp(i+4) = length(find(PrefShift_cat(hVis_cat&bead_cat&dirSel&DS)==i));
end
categories_bn = []; categories_bp = [];
for i = 1:8
    categories_bn = [categories_bn; repmat(i,nPerCat_bn(i),1)];
    categories_bp = [categories_bp; repmat(i,nPerCat_bp(i),1)];
end
categories = [categories_bn;categories_bp];
[~,chi2stat,pval] = crosstab(dataID,categories);
subplot(3,2,3);
title(sprintf('chi2stat = %.3g; p = %.3g',chi2stat ,pval))

% % chi2 test for oS cells
nDSnBn = size(PrefShift_cat(hVis_cat&~bead_cat&dirSel&~DS),1);
nDSnBp = size(PrefShift_cat(hVis_cat&bead_cat&dirSel&~DS),1);
dataID = [repmat('a',nDSnBn,1);repmat('b',nDSnBp,1)];
nPerCat_bn = zeros(1,7); nPerCat_bp = zeros(1,7);
for i = -3:4
    nPerCat_bn(i+4) = length(find(PrefShift_cat(hVis_cat&~bead_cat&dirSel&~DS)==i));
    nPerCat_bp(i+4) = length(find(PrefShift_cat(hVis_cat&bead_cat&dirSel&~DS)==i));
end
categories_bn = []; categories_bp = [];
for i = 1:8
    categories_bn = [categories_bn; repmat(i,nPerCat_bn(i),1)];
    categories_bp = [categories_bp; repmat(i,nPerCat_bp(i),1)];
end
categories = [categories_bn;categories_bp];
[~,chi2stat,pval] = crosstab(dataID,categories);
subplot(3,2,4);
title(sprintf('chi2stat = %.3g; p = %.3g',chi2stat ,pval))

%  ---------- how to compare shift to no shift? ----------
% DS_all = PrefShift_cat(hVis_cat&dirSel&DS);
% DS_all(DS_all(:,1) ~= 0) = 1;
% OS_all = PrefShift_cat(hVis_cat&dirSel&~DS);
% p = signtest(DS_all)
% ------------------------------------------------------------

%% 2.5.4 ------ Separate cells in OS vs DS, DSI ----------------
DSI_Figure4 = figure;
for i = 1:2
    data_violin{i} = DSI_cat(hVis_cat&~bead_cat&dirSel&DS,i);
    data_violin{i+2} = DSI_cat(hVis_cat&bead_cat&dirSel&DS,i);
    data_violin{i+4} = DSI_cat(hVis_cat&~bead_cat&dirSel&~DS,i);
    data_violin{i+6} = DSI_cat(hVis_cat&bead_cat&dirSel&~DS,i);
end
for i = 8:-1:1 % if you increment i, it's gonna delete the empty cell and change the dimension of your array and bug when reaching the last cell(s)
    if isempty(data_violin{1,i})
        
        data_violin(i)=[];
    end
end

pos1 = [0.1 0.55 0.9 0.4];
subplot('Position',pos1); hold on %[left bottom width height]
v = violin(data_violin,'facecolor',color_bn,'edgecolor','none','FaceAlpha',1);
% violin(data_violin,'xlabel',{'DS bead-','DS bead+','OS bead-','OS bead+'},...
%     'facecolor',color_bn,'edgecolor','none','FaceAlpha',1);
for i = 2:2:size(data_violin,2)
    v(i).FaceColor = color_bp;
end
for i = 1:2:size(data_violin,2)
    v(i).FaceAlpha = 0.2;
end
xl = xlim;
plot(xlim,[0.35 0.35],'k:')
xl = xlim;
x_size = xl(1)+xl(2);
plot([x_size/2 x_size/2],ylim,'k--')
text(xl(1),0.35,'{\it DS cells threshold}','FontSize',8,'Color',[.5 .5 .5],...
'HorizontalAlignment','left','VerticalAlignment','top')
ylabel('DSI')
legend([v(1),v(2),v(3),v(4)],{'bead- laser OFF','bead- laser ON','bead+ laser OFF','bead+ laser ON'},...
    'NumColumns',2,'FontSize',6)
legend boxoff
xticks([])
title ('DS cells  | OS cells')
clear data_violin
% subplot(2,3,2); hold on
% for i = 1:2
%     % plot([i-0.2 i-0.1],[DSI_bn(i) DSI_bn(i)],'k', 'LineWidth',2)
%     % plot([i+0.1 i+0.2],[DSI_bp(i) DSI_bp(i)],'k', 'LineWidth',2)
%     scatter(i*ones(sum(hVis_cat&~bead_cat&dirSel),1)+2*b1(1).XOffset,DSI_cat(hVis_cat&~bead_cat&dirSel,i),...
%         'MarkerFaceColor',color_bn,'MarkerEdgeColor','none','MarkerFaceAlpha',0.05)
%     scatter(i*ones(sum(hVis_cat&bead_cat&dirSel),1)+2*b1(2).XOffset,DSI_cat(hVis_cat&bead_cat&dirSel,i),...
%         'MarkerFaceColor',color_bp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% end

pos2 = [0.1 0.15 0.9 0.4];
subplot('Position',pos2); hold on
deltaDSI = DSI_cat(:,2) - DSI_cat(:,1);
data_violin{1} = deltaDSI(hVis_cat&~bead_cat&dirSel&DS,:);
data_violin{2} = deltaDSI(hVis_cat&bead_cat&dirSel&DS,:);
data_violin{3} = deltaDSI(hVis_cat&~bead_cat&dirSel&~DS,:);
data_violin{4} = deltaDSI(hVis_cat&bead_cat&dirSel&~DS,:);
for i = 4:-1:1 % if you increment i, it's gonna delete the empty cell and change the dimension of your array and bug when reaching the last cell(s)
    if isempty(data_violin{1,i})
        data_violin(i)=[];
    end
end
v = violin(data_violin,...
    'facecolor',color_bn,'edgecolor','none','FaceAlpha',1);
for i = 2:2:size(data_violin,2)
    v(i).FaceColor = color_bp;
end
plot(xlim,[0 0],'k:')
xl = xlim;
x_size = xl(1)+xl(2);
plot([x_size/2 x_size/2],ylim,'k--')
legend([v(1),v(2)],{'bead-','bead+'},...
    'NumColumns',2,'FontSize',6, 'Location', 'northeast')
legend box off
xticks([])
ylabel('delta DSI')
clear data_violin

if looped_cells
    if isfield(data, 'beads_pos') && sum(hVis_cat&bead_cat&dirSel&DS)>=4
        if ~adtest(deltaDSI(hVis_cat&~bead_cat&dirSel&DS)) && ~adtest(deltaDSI(hVis_cat&bead_cat&dirSel&DS))
            [~,p1] = ttest2(deltaDSI(hVis_cat&~bead_cat&dirSel&DS),deltaDSI(hVis_cat&bead_cat&dirSel&DS));
            testName1 = 'ttest';
        else
            p1 = ranksum(deltaDSI(hVis_cat&~bead_cat&dirSel&DS),deltaDSI(hVis_cat&bead_cat&dirSel&DS));
            testName1 = 'wilcoxon';
        end
    else
        testName1 = '';
    end
else
    testName2 = ''; p1 = [];
end
xtickangle(45)
title([testName, ', p = ' num2str(p,3)]);

if looped_cells
    if isfield(data, 'beads_pos') && sum(hVis_cat&bead_cat&dirSel&DS)>4
        if ~adtest(deltaDSI(hVis_cat&~bead_cat&dirSel&~DS)) && ~adtest(deltaDSI(hVis_cat&bead_cat&dirSel&~DS))
            [~,p2] = ttest2(deltaDSI(hVis_cat&~bead_cat&dirSel&~DS),deltaDSI(hVis_cat&bead_cat&dirSel&~DS));
            testName2 = 'ttest';
        else
            p2 = ranksum(deltaDSI(hVis_cat&~bead_cat&dirSel&~DS),deltaDSI(hVis_cat&bead_cat&dirSel&~DS));
            testName2 = 'wilcoxon';
        end
        title([testName2, ', p = ' num2str(p1,3), ' | ', testName2, ', p = ' num2str(p2,3)],...
            'FontSize',10);
    else
        testName2 = ''; p1 = [];
    end
else
    testName2 = ''; p1 = [];
end

pos3 = [0.1 0.05 0.9 0.1];
subplot('Position',pos3); axis off
if isfield(data, 'beads_pos')
    if ~adtest(deltaDSI(hVis_cat&dirSel&DS))
        [~,p1] = ttest(deltaDSI(hVis_cat&dirSel&DS));
        testName1 = 'ttest';
    else
        p1 = signtest(deltaDSI(hVis_cat&dirSel&DS));
        testName1 = 'wilcoxon';
    end
else
    testName1 = '';
end
if isfield(data, 'beads_pos') && sum(hVis_cat&dirSel&~DS)>=4
    if ~adtest(deltaDSI(hVis_cat&dirSel&~DS))
        [~,p2] = ttest(deltaDSI(hVis_cat&dirSel&~DS));
        testName2 = 'ttest';
    else
        p2 = signtest(deltaDSI(hVis_cat&dirSel&~DS));
        testName2 = 'wilcoxon';
    end
    text(0.5,0.5,{'-----------------------','All cells',['DS cells, ' testName1, ' p  =  ',num2str(p1,3) ' | OS cells, ' testName2 ' p = ' num2str(p2,3)]},...
        'HorizontalAlignment','center','Fontsize',10)
else
    text(0.5,0.5,{'-----------------------','All cells',['DS cells, ' testName1, ' p  =  ' ' | OS cells, ' testName2 ' p = ' ]},...
        'HorizontalAlignment','center','Fontsize',10)
end


%%
%% save
TuningData = data;

if saveFig
    if isfield(data.Info,'saveName') && isfield(data.Info,'mainDir')
        saveName = data.Info.saveName;
        mainDir = data.Info.mainDir;
    else
        if isempty(varargin) || nargin == 4
            mainDir = uigetdir('E:\TempData','Choose Main dir');
            answer = inputdlg('savename');
            saveName = answer{1,1};
        else
            mainDir = varargin{1};
            saveName = varargin{2};
        end
    end
    tv = datestr(now, 'yyyy_mm_dd');
    saveName = [saveName '_' 'errorTh0' num2str(error_th*10) '_' tv ];
    
    saveDir = [mainDir filesep 'Figures' filesep 'PopPlots' filesep 'DirectionTuning'];
    if DirSelection == 0
        saveDir = [saveDir filesep 'AllVisRespCells'];
    elseif DirSelection == 1
        saveDir = [saveDir filesep 'GaussianFitSelected'];
    elseif DirSelection == 2
        saveDir = [saveDir filesep 'AnovaSelected'];
    end
    if ~exist(saveDir)
        mkdir(saveDir)
    end
    
    if exist('MF_Figure')
        figure(MF_Figure)
        saveas(gcf,[saveDir filesep saveName,'_ModulationFactor'],'fig');
        screen2png([saveDir filesep saveName,'_ModulationFactor']);
    end
    
    if exist('DirTunAll_Fig')
        figure(DirTunAll_Fig)
        screen2png([saveDir filesep saveName,'_TuningCurve_SF1vsSF2']);
    end
    
    figure(PolarFig)
    saveas(gcf,[saveDir filesep saveName,'_PolarPlots'],'fig');
    screen2png([saveDir filesep saveName,'_PolarPlots']);
    
    figure(PolarFig2)
    saveas(gcf,[saveDir filesep saveName,'_PolarPlots_Norm'],'fig');
    screen2png([saveDir filesep saveName,'_PolarPlots_Norm']);
    
    figure(VectorLength_Fig)
    saveas(gcf,[saveDir filesep saveName,'_VectorLength'],'fig');
    screen2png([saveDir filesep saveName,'_VectorLength']);
    
    figure(VectLengthVsAct_Fig)
    saveas(gcf,[saveDir filesep saveName,'_VectorLengthVsAct'],'fig');
    screen2png([saveDir filesep saveName,'_VectorLengthVsAct']);
    
    figure(DirTun_Bead_Fig)
    saveas(gcf,[saveDir filesep saveName,'_TuningCurve_BpVsBn'],'fig');
    screen2png([saveDir filesep saveName,'_TuningCurve_BpVsBn']);
    
    figure(ModelInhib_Figure1)
    saveas(gcf,[saveDir filesep saveName,'_ModeledInhib_1'],'fig');
    screen2png([saveDir filesep saveName,'_ModeledInhib_1']);
    
    figure(ModelInhib_Figure2)
    saveas(gcf,[saveDir filesep saveName,'_ModeledInhib_2'],'fig');
    screen2png([saveDir filesep saveName,'_ModeledInhib_2']);
    
    if exist('DSI_Figure1')
        figure(DSI_Figure1)
        saveas(gcf,[saveDir filesep saveName,'_DSI_1'],'fig');
        screen2png([saveDir filesep saveName,'_DSI_1']);
    end
    
    figure(DSI_Figure2)
    saveas(gcf,[saveDir filesep saveName,'_DSI_2'],'fig');
    screen2png([saveDir filesep saveName,'_DSI_2']);
    
    figure(DSI_Figure3)
    saveas(gcf,[saveDir filesep saveName,'_DSI_3'],'fig');
    screen2png([saveDir filesep saveName,'_DSI_3']);
    
    figure(Sigma_Fig)
    saveas(gcf,[saveDir filesep saveName,'_Sigma'],'fig');
    screen2png([saveDir filesep saveName,'_Sigma']);
    
    figure(Slope_AllCells)
    saveas(gcf,[saveDir filesep saveName,'_Slope_AllCells'],'fig');
    screen2png([saveDir filesep saveName,'_Slope_AllCells']);
    
    figure(SlopeShuffle_Fig)
    saveas(gcf,[saveDir filesep saveName,'_Slope_Shuffle'],'fig');
    screen2png([saveDir filesep saveName,'_Slope_Shuffle']);
    
    figure(DSI_Figure4)
    saveas(gcf,[saveDir filesep saveName,'_DSI_Figure4'],'fig');
    screen2png([saveDir filesep saveName,'_DSI_Figure4']);
end

end
