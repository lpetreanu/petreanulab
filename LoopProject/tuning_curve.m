function tuning = tuning_curve(dir_tuning,h_vis,beads_pos,saveFig,varargin)
color_bp = [0.8500 0.3250 0.0980];
color_bn = [0      0.4470 0.7410];

nSF = 2;

% y = squeeze(dir_tuning(:,:,1,:));
y = circshift(dir_tuning(:,1:8,:,:),-2,2);
y(:,9,:,:) = y(:,1,:,:);
directions = -90:45:270;

options = fitoptions('gauss2');
options.Lower = [0 0 0 0 180 0];
options.Upper = [Inf 0 Inf Inf 180 Inf];

for i = 1:size(y,1)
    for  ii = 1:2 %nSF
        y1 = squeeze(y(i,:,ii,1));
        % figure; hold on
        f1 = fit(directions.',y1.','gauss2',options);
        % plot(f1,'r--',directions,y1,'r')
        
        y2 = squeeze(y(i,:,ii,2));
        f2 = fit(directions.',y2.','gauss2',options);
        % plot(f2,'b--',directions,y2,'b')
        
        MF(i,ii,1) = (f1.a1-f2.a1)/(f1.a1+f2.a1);
        MF(i,ii,2) = (f1.a2-f2.a2)/(f1.a2+f2.a2);
        
        for iii = 1:length(directions)
            diff(iii,1) = (y1(iii) - f1(directions(iii)))^2;
            diff(iii,2) = (y2(iii) - f2(directions(iii)))^2;
        end
        
        Error(i,ii,1) = sum(diff(:,1))./sum(y1.^2);
        Error(i,ii,2) = sum(diff(:,2))./sum(y2.^2);
    end
end

k = Error(:,:,1) < 0.1 & Error(:,:,2) < 0.1;

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

bead = false(size(y,1),1);
bead(beads_pos(:,3)) = true;

if nSF == 2
    MF_cat = squeeze(cat(1,MF(:,1,:),MF(:,2,:)));
    hVis_cat = cat(1,h_vis(:,1),h_vis(:,2));
    k_cat = cat(1,k(:,1),k(:,2));
    bead_cat = [bead ; bead];
else
    MC_cat = squeeze(MF);
    hVis_cat = h_vis;
    k_cat = k;
end

fit_sel_Bpos = mean(MF_cat(hVis_cat&k_cat&bead_cat,:),1);
fit_sel_Bpos_sem = std(MF_cat(hVis_cat&k_cat&bead_cat,:),[],1)./sqrt(sum(hVis_cat&k_cat&bead_cat));
fit_sel_Bneg = mean(MF_cat(hVis_cat&k_cat&~bead_cat,:),1);
fit_sel_Bneg_sem = std(MF_cat(hVis_cat&k_cat&~bead_cat,:),[],1)./sqrt(sum(hVis_cat&k_cat&~bead_cat));
all_Bpos = mean(MF_cat(hVis_cat&bead_cat,:),1);
all_Bneg = mean(MF_cat(hVis_cat&~bead_cat,:),1);
all_Bpos_sem = std(MF_cat(hVis_cat&bead_cat,:),[],1)./sqrt(sum(hVis_cat&bead_cat));
all_Bneg_sem = std(MF_cat(hVis_cat&~bead_cat,:),[],1)./sqrt(sum(hVis_cat&~bead_cat));

figure
subplot(1,3,1); hold on
b2 = bar([fit_sel_Bpos' fit_sel_Bneg'],'FaceColor','flat');
drawnow
for i = 1:2
    plot([i+b2(1).XOffset i+b2(1).XOffset],...
        [fit_sel_Bpos(i)+fit_sel_Bpos_sem(i) fit_sel_Bpos(i)-fit_sel_Bpos_sem(i)],'-','Color',color_bp)
    plot([i+b2(2).XOffset i+b2(2).XOffset],...
        [fit_sel_Bneg(i)+fit_sel_Bneg_sem(i) fit_sel_Bneg(i)-fit_sel_Bneg_sem(i)],'-','Color',color_bn)
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
title({'Selected data','(fit error < 0.1)'})
ylabel('Opto Modulation Factor')

subplot(1,3,3); hold on
b3 = bar([all_Bpos' all_Bneg'],'FaceColor','flat');
drawnow
for i = 1:2
    plot([i+b3(1).XOffset i+b3(1).XOffset],...
        [all_Bpos(i)+all_Bpos_sem(i) all_Bpos(i)-all_Bpos_sem(i)],'-','Color',color_bp)
    plot([i+b3(2).XOffset i+b3(2).XOffset],...
        [all_Bneg(i)+all_Bneg_sem(i) all_Bneg(i)-all_Bneg_sem(i)],'-','Color',color_bn)
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
title('All data')

subplot(1,3,2); hold on
axis off
lg1.Position = [0.35 0.8 0.15 0.1];
lg1.ItemTokenSize = [10,4];
lg2.Position = [0.49 0.2 0.15 0.1];
lg2.ItemTokenSize = [10,4];

tuning.modulation_factor.MF = MF;
tuning.modulation_factor.error_selec = k;

set(gcf,'Units','Normalized','Position',[0.05 0.05 0.7 0.5])
if saveFig
    if isempty(varargin) || nargin == 4
        mainDir = uigetdir('E:\TempData','Choose Main dir');
        answer = inputdlg('savename');
        saveName = answer{1,1};
    else
        mainDir = varargin{1};
        saveName = varargin{2};
    end
    saveDir = [mainDir filesep 'Figures' filesep 'PopPlots'];
    saveas(gcf,[saveDir filesep saveName,'_ModulationFactor'],'fig');
    screen2png([saveDir filesep saveName,'_ModulationFactor']);
end

%% normalize amplitude and plot all
% -------- segregate SF/TF, all cells

max1 = dir_tuning(:,5,:,1);
max1 = repmat(max1,1,9,1,2);
% counter = zeros(173,2);
% for i = 1:size(dir_tuning,1)
%     for ii = 1:2
%         if max1(i,:,ii,:) < 0
%             disp('yo')
%             max1(i,:,ii,:) = NaN(1,9,1,2);
%             counter(i,ii) = 1;
%         end
%     end
% end
test2 = dir_tuning./max1;

test(:,1,:) = squeeze(nanmean(test2(h_vis(:,1),:,1,:),1));
test(:,2,:) = squeeze(nanmean(test2(h_vis(:,2),:,2,:),1));

sem(:,1,:) = std(test2(h_vis(:,1),:,1,:),[],1,'omitnan')./sqrt(sum(h_vis(:,1)));
sem(:,2,:) = std(test2(h_vis(:,2),:,2,:),[],1,'omitnan')./sqrt(sum(h_vis(:,2)));

test_merged = squeeze(cat(1,test2(:,:,1,:),test2(:,:,2,:)));
merged_mean = squeeze(nanmean(test_merged(hVis_cat,:,:)));
merged_sem = squeeze(std(test_merged(hVis_cat,:,:)))./sqrt(sum(hVis_cat));


figure;
subplot(1,4,1); hold on
p1 = plot(test(:,1,1),'r','LineWidth',1);
p2 = plot(test(:,1,2),'ro--','LineWidth',1);
p3 = plot(test(:,2,1),'b','LineWidth',1);
p4 = plot(test(:,2,2),'bo--','LineWidth',1);
p5 = plot(merged_mean(:,1),'k','LineWidth',2);
p6 = plot(merged_mean(:,2),'k--','LineWidth',2);

for i = 1:size(test,1)
    plot([i i], [test(i,1,1)+sem(i,1,1) test(i,1,1)-sem(i,1,1)],'r-')
    plot([i i], [test(i,1,2)+sem(i,1,2) test(i,1,2)-sem(i,1,2)],'r-')
    plot([i i], [test(i,2,1)+sem(i,2,1) test(i,2,1)-sem(i,2,1)],'b')
    plot([i i], [test(i,2,2)+sem(i,2,2) test(i,2,2)-sem(i,2,2)],'b-')
    plot([i i], [merged_mean(i,1)+merged_sem(i,1) merged_mean(i,1)-merged_sem(i,1)],'k-')
    plot([i i], [merged_mean(i,2)+merged_sem(i,2) merged_mean(i,2)-merged_sem(i,2)],'k-')
end
xticks([1:2:9]); xticklabels({'-180','-90','Pref.','+90','+180'})
title('Averaged data')
xl=xlim;yl=ylim;
text(p1.XData(1),yl(2),['all data, n = ',num2str(sum(hVis_cat))],...
    'HorizontalAlignment','Left','VerticalAlignment','top')
lgd = legend([p1,p2,p3,p4,p5,p6],{'SF/TF 1, no laser','SF/TF 1, laser','SF/TF 2, no laser','SF/TF 2, laser',...
    'All SFs/TFs, no laser','All SFs/TFs, laser'},'location','none');
subplot(1,4,2); axis off
lgd.Position =  [0.3 0.5 0.15 0.15];


subplot(1,4,3); hold on
plot(test(:,1,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
plot(test(:,1,2),'r','LineWidth',2)
xticks([1:2:9]); xticklabels({'-180','-90','Pref.','+90','+180'})
title('SF/TF 1, all data')


subplot(1,4,4); hold on
plot(test(:,2,1),'Color',[0.5 0.5 0.5],'LineWidth',2)
plot(test(:,2,2),'b','LineWidth',2)
xticks([1:2:9]); xticklabels({'-180','-90','Pref.','+90','+180'})
title('SF/TF 2, all data')

for i = 1:sum(h_vis(:,1))
    if h_vis(i,1)
        subplot(1,4,3)
        lh = plot(test2(i,:,1,2));
        lh.Color=[1,0,0,0.2];
    end
    if h_vis(i,2)
        subplot(1,4,4)
        lh = plot(test2(i,:,2,2));
        lh.Color=[0,0,1,0.2];
    end
end

subplot(1,4,3)
xl=xlim;yl=ylim;
text(lh.XData(1),yl(2),['n = ',num2str(sum(h_vis(:,1)))],...
    'HorizontalAlignment','Left','VerticalAlignment','top')

subplot(1,4,4)
xl=xlim;yl=ylim;
text(lh.XData(1),yl(2),['n = ',num2str(sum(h_vis(:,2)))],...
    'HorizontalAlignment','Left','VerticalAlignment','top')

set(gcf,'Units','Normalized','Position',[0.05 0.4 0.7 0.5])
sgtitle('Impact of Opto stim of FB on tuning curve')
if saveFig
    %     if isempty(varargin) || nargin == 4
    %         mainDir = uigetdir('E:\TempData','Choose Main dir');
    %         answer = inputdlg('savename');
    %         saveName = answer{1,1};
    %     else
    %        mainDir = varargin{1};
    %        saveName = varargin{2};
    %     end
    %     saveDir = [mainDir filesep 'Figures' filesep 'PopPlots'];
%     saveas(gcf,[saveDir filesep saveName,'_tuning_curve'],'fig');
    screen2png([saveDir filesep saveName,'_tuning_curve']);
end


% -------- segregate bead+ and bead neg
tuning_bp = squeeze(mean(test_merged(hVis_cat&bead_cat,:,:)));
tuning_bn = squeeze(mean(test_merged(hVis_cat&~bead_cat,:,:)));
tuning_bp_sem = squeeze(std(test_merged(hVis_cat&bead_cat,:,:))./sqrt(sum(hVis_cat&bead_cat)));
tuning_bn_sem = squeeze(std(test_merged(hVis_cat&~bead_cat,:,:))./sqrt(sum(hVis_cat&~bead_cat)));

figure
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
lg3 = legend([p7, p8, p9, p10],{['bead+ laser OFF, n = ' num2str(sum(hVis_cat&bead_cat))],...
    'bead+ laser ON',...
    ['bead- laser OFF, n = ' num2str(sum(hVis_cat&~bead_cat))],...
    'bead- laser ON'})
title('Session Average')

subplot(1,4,3); hold on
p7 = plot(tuning_bp(:,1),'-','Color',color_bp,'LineWidth',2);
p8 = plot(tuning_bp(:,2),'--','Color',color_bp,'LineWidth',2);
bla = test_merged(hVis_cat&bead_cat,:,:);
for i = 1:sum(hVis_cat&bead_cat)
    lh = plot(squeeze(bla(i,:,2)),'-','Color',color_bp);
    lh.Color=[color_bp 0.2];
end
xlim([0.5 9.5]); xticks(1:2:9);
xticklabels({'-180','-90','Pref.','+90','+180'});
xlabel('Directions'); ylabel('dFF, norm')
title('All bead+ cells')
% xl = xlim; yl = ylim;
% text(p7.XData(1),yl(2),['n = ',num2str(sum(hVis_cat&bead_cat))],...
%     'HorizontalAlignment','Left','VerticalAlignment','top')
clear bla

subplot(1,4,4); hold on
plot(tuning_bn(:,1),'-','Color',color_bn,'LineWidth',2);
plot(tuning_bn(:,2),'--','Color',color_bn,'LineWidth',2);
bla = test_merged(hVis_cat&~bead_cat,:,:);
for i = 1:sum(hVis_cat&~bead_cat)
    lh = plot(squeeze(bla(i,:,2)),'-','Color',color_bn);
    lh.Color=[color_bn 0.2];
end
xticks(1:2:9); xticklabels({'-180','-90','Pref.','+90','+180'});
xlabel('Directions'); ylabel('dFF, norm')
title('All bead- cells')

subplot(1,4,2); axis off
lg3.Position = [0.3 0.7 0.2 0.2];
lg3.ItemTokenSize = [20,4];
set(gcf,'Units','Normalized','Position',[0.05 0.4 0.7 0.5])
if saveFig
    saveas(gcf,[saveDir filesep saveName,'_tuning_curve_BeadSeg'],'fig');
    screen2png([saveDir filesep saveName,'_tuning_curve_BeadSeg']);
end
end