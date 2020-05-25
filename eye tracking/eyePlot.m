clearvars
%% set group to analyze
fr = 15; %framerate of eye camera
co = [0 0.35 0.57;
    0.12 0.63 0.67]; % color scheme
%set some things
prot = 'natImgs';
exp = 'HVA';
pthresh = 2; %pixel cutoff for considering eye at "center"
plt = 0; %whether to plot
%% session files where at least one day has tracking
mice = {'GF338','GF338','GF339','GF341','GF343','GF343','GF343','GF344','GF344','GF346','GF346'};
cno_day = {'20190618','20190702','20190802','20190828','20190903','20190912',...
    '20190917','20190903','20190912','20190927','20191004'};
sal_day = {'20190617','20190701','20190801','20190827','20190902','20190911',...
    '20190916','20190902','20190911','20190926','20191003'};
hva_ip = [0,0,1,1,0,1,1,0,1,1,1];
hva_v1 = [0,0,0,0,1,0,0,1,0,0,0];
dds = logical([0,0,1,1,1,1,1,1,1,1,1]); %whether animal expressed DREADDS (1)
ind = find(hva_ip);

switch exp
    case 'CNO ctl'
        ind = find(~dds);
    case 'HVA'
        ind = find(hva_ip);
    case 'HVA-proj'
        ind = find(hva_v1);
    case 'all'
        ind = 1:length(hva_ip);
end

%%
sal_pup = struct(); cno_pup = struct();
for m = 1:length(mice)
    mouse = mice{m};
    
    % SALINE DAY
    dt = circshift(sal_day{m},4);
    
    % load eye tracking data
    cd(sprintf('E:\\quietFB\\eyeVids\\%s\\%s\\natImgs',mouse,dt))
    
    try
        load('pupil_data.mat')
        
        % pull out trials when tracking worked
        eye = checkTrials(mouse,circshift(dt,4),eye);
        
        % combine variables across trials
        x = vertcat(eye.px{:});
        y = vertcat(eye.py{:});
        d = vertcat(eye.pd{:});
        
        %load running (produced by trialSepByRun_xSess.m)
        load(['E:\quietFB\run\' mouse '_' sal_day{m} '_natImg_sal.mat']) %'natImg' 'grating' 'movie' 'dark'
        [tRun, xRun] = prepRun(speedTr);
        
        %plot traces
        if plt
            figure('Name',[mouse dt],'Position',[20 500 1800 400]);
            plotTraces(x,y,d,mouse,dt,fr,'k',tRun,xRun)
        end
        
        % add to struct
        sal_pup(m).d = d;
        sal_pup(m).x = x;
        sal_pup(m).y = y;
        sal_pup(m).tRun = tRun;
        sal_pup(m).xRun = xRun;
        
        % save trials when pupil is @center
        tr = centerPupil(x,pthresh,fr);
        save(sprintf('E:\\quietFB\\eye\\%s_%s_centerTrials.mat',mouse,dt),'tr')
        fprintf('%d trials at center\n',sum(tr))
    catch
        fprintf('Tracking doesn''t exist for %s %s\n',mouse,dt)
    end
    
    % CNO DAY
    dt = circshift(cno_day{m},4);
    
    % load eye tracking data
    cd(sprintf('E:\\quietFB\\eyeVids\\%s\\%s\\natImgs',mouse,dt))
    
    try
        load('pupil_data.mat')
        
        % pull out trials when tracking worked
        eye = checkTrials(mouse,circshift(dt,4),eye);
        
        % combine variables across trials
        x = vertcat(eye.px{:});
        y = vertcat(eye.py{:});
        d = vertcat(eye.pd{:});
        
        %load running (produced by trialSepByRun_xSess.m)
        load(['E:\quietFB\run\' mouse '_' cno_day{m} '_natImg_cno.mat']) %'natImg' 'grating' 'movie' 'dark'
        [tRun, xRun] = prepRun(speedTr);
        
        %plot traces
        if plt
        figure('Name',[mouse dt],'Position',[20 20 1800 400]);
        plotTraces(x,y,d,mouse,dt,fr,'m',tRun,xRun)
        end
        
        % add to struct
        cno_pup(m).d = d;
        cno_pup(m).x = x;
        cno_pup(m).y = y;
        cno_pup(m).tRun = tRun;
        cno_pup(m).xRun = xRun;
        
        % save trials when pupil is @center
        tr = centerPupil(x,pthresh,fr);
        save(sprintf('E:\\quietFB\\eye\\%s_%s_centerTrials.mat',mouse,dt),'tr')
        fprintf('%d trials at center\n',sum(tr))
    catch
        fprintf('Tracking doesn''t exist for %s %s\n',mouse,dt)
    end
end

%% pupil diameter by condition
s = vertcat(sal_pup(ind).d);
c = vertcat(cno_pup(ind).d);
s = s(~isnan(s));
c = c(~isnan(c));
[~,p] = ttest2(s,c);
figure,
h=histogram(s,50,'Normalization','probability','FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
hold on, histogram(c,h.BinEdges,'Normalization','probability','FaceColor',[1 0 1],'EdgeColor',[1 0 1])
legend('saline','cno')
xlabel('Pupil diameter [pixels]'); ylabel('Probability')
title(sprintf('%s effect on pupil diameter\nt-test p=%1.2f',exp,p))

%% saccades by condition
% sx = vertcat(sal_pup(ind).x);
% sy = vertcat(sal_pup(ind).y);
% cx = vertcat(cno_pup(ind).x);
% cy = vertcat(cno_pup(ind).y);

for j = ind
    fprintf('Looking at %s %s saline session...\n',mice{j},sal_day{j})
    x = sal_pup(j).x;
    s_on = saccadeDetector(x,fr,plt,mice{j},sal_day{j});
    sal_pup(j).s_on = s_on;
    
    fprintf('Looking at %s %s CNO session...\n',mice{j},cno_day{j})
    x = cno_pup(j).x;
    s_on = saccadeDetector(x,fr,plt,mice{j},cno_day{j});
    cno_pup(j).s_on = s_on;
     waitforbuttonpress()
end
%% plot saccade likelihood during running
f1 =  fullfig;
f2 = figure('Position',[600 70 560 900]);
count=1;
for j = ind
    %****SALINE*****
    %interpolate run
    xq = (1:length(sal_pup(j).x))/fr;
    vq = interp1(sal_pup(j).tRun,sal_pup(j).xRun,xq);
    
    % find percentage of saccades during running
    rs = sum(vq(sal_pup(j).s_on)>1); %# of saccades during running
    ts = length(sal_pup(j).s_on);
    pRun = sum(vq>1)/length(vq);
    pTracked = 100*sum(~isnan(sal_pup(j).x))/length(sal_pup(j).x);
    
    %save saccade metrics
    sal_pup(j).sacRunRatio = rs/ts/pRun;
    sal_pup(j).pTracked = pTracked;
    sal_pup(j).nrs = rs;
    sal_pup(j).ns = ts;
    
    %plot saccade and running likelihood
    figure(f2)
    subplot(length(ind),2,count)
    bar(100*[rs/ts pRun])
    set(gca,'XTickLabels',{'runSac', 'run'})
    title(sprintf('%s saline\nNsac=%i eye tracking@%2.0f%%\nRunSac:Run Ratio=%2.1f',mice{j},ts,pTracked,rs/ts/pRun))
    
    % plot saccades and running throughout session
    figure(f1)
    subplot(length(ind),2,count),plot(sal_pup(j).tRun,sal_pup(j).xRun,'k'); hold on
    stem(sal_pup(j).s_on/fr,ones(size(sal_pup(j).s_on))*10,'g')
    title(sprintf('%s saline',mice{j}))
    count = count + 1;
    
    %****CNO*****
    %interpolate run
    xq = (1:length(cno_pup(j).x))/fr;
    vq = interp1(cno_pup(j).tRun,cno_pup(j).xRun,xq);
    
    % find percentage of saccades during running
    rs = sum(vq(cno_pup(j).s_on)>1); %# of saccades during running
    ts = length(cno_pup(j).s_on);
    pRun = sum(vq>1)/length(vq);
    pTracked = 100*sum(~isnan(cno_pup(j).x))/length(cno_pup(j).x);
    
    %save saccade metrics
    cno_pup(j).sacRunRatio = rs/ts/pRun;
    cno_pup(j).pTracked = pTracked;
    cno_pup(j).nrs = rs;
    cno_pup(j).ns = ts;
    
    %plot saccade and running likelihood
    figure(f2)
    subplot(length(ind),2,count)
    bar(100*[rs/ts pRun])
    set(gca,'XTickLabels',{'runSac', 'run'})
    title(sprintf('%s cno\nNsac=%i eye tracking@%2.0f%%\nRunSac:Run Ratio=%2.1f',mice{j},ts,pTracked,rs/ts/pRun))
    
    % plot saccades and running throughout session
    figure(f1)
    subplot(length(ind),2,count),plot(cno_pup(j).tRun,cno_pup(j).xRun,'k'); hold on
    stem(cno_pup(j).s_on/fr,ones(size(cno_pup(j).s_on))*10,'g')
    title(sprintf('%s cno',mice{j}))
    if count == 2
        legend('speed','saccade')
    end
    count = count + 1;
end

% compare ratio of running saccades to running; number of saccades
[~,p] = ttest(vertcat(sal_pup(:).sacRunRatio),vertcat(cno_pup(:).sacRunRatio));
figure; subplot(1,2,1)
scatter(vertcat(sal_pup(:).sacRunRatio),vertcat(cno_pup(:).sacRunRatio),30,...
    'k','filled','MarkerFaceAlpha',0.5); hold on
plot([0 15],[0 15],'k'); axis square
xlabel('saline');ylabel('cno')
title(sprintf('Ratio of saccades during running\nto total running time\np = %1.2f',p))

[~,p] = ttest(vertcat(sal_pup(:).ns),vertcat(cno_pup(:).ns));
subplot(1,2,2)
scatter(vertcat(sal_pup(:).ns),vertcat(cno_pup(:).ns),30,...
    'k','filled','MarkerFaceAlpha',0.5); hold on
plot([0 150],[0 150],'k'); axis square
xlabel('saline');ylabel('cno')
title(sprintf('Number of saccades\np = %1.2f',p))
sal_runsac=100*nanmean(vertcat(sal_pup(:).nrs)./vertcat(sal_pup(:).ns));
cno_runsac=100*nanmean(vertcat(cno_pup(:).nrs)./vertcat(cno_pup(:).ns));

fprintf('HVA intact: %2.1f%% of saccades during running\n',sal_runsac)
fprintf('HVA silenced: %2.1f%% of saccades during running\n',cno_runsac)

%%
mouse = 'GF341'; dt = '20190828';
x = cno_pup(4).x; d = cno_pup(4).d; y = cno_pup(4).y;
x(isnan(x)) = []; y(isnan(y)) = []; d(isnan(d)) = [];
figure('Name',[mouse dt],'Position',[20 20 1800 400]);
plotTraces(lowpass(x,5,fr),y,d,mouse,dt,fr,'k',tRun,xRun)
%%
tt = [fr*333:fr*390,fr*403:fr*484];
figure; scatter(cno_pup(4).x(tt)-nanmean(cno_pup(4).x),cno_pup(4).d(tt),25, tt,'filled','MarkerFaceAlpha',0.5)
xlabel('x pos [px]'); ylabel('pup diameter [px]')
figure,histogram(x-mean(x)); xlabel('x pos [px]')
% s_on = saccadeDetector(x,fr,1);
%% OTHER PLOTS/PLOTTING FUNCTIONS
% figure('Position',[20 550 1800 500]);
% subplot(3,1,1);
% imagesc(reshape(x-nanmean(x),15,1500)); colorbar
% set(gca,'YTick',3:3:12,'YTickLabel',(3:3:12)/15)
% title('x')
% subplot(3,1,2);
% imagesc(reshape(y-nanmean(y),15,1500)); colorbar
% set(gca,'YTick',3:3:12,'YTickLabel',(3:3:12)/15)
% title('y'); ylabel('Time (s)')
% subplot(3,1,3);
% imagesc(reshape(d,15,1500)); colorbar
% set(gca,'YTick',3:3:12,'YTickLabel',(3:3:12)/15)
% title('diameter [pixels]')
% xlabel('trial #')
% colormap viridis
j = 3;
fprintf('Looking at %s %s saline session...\n',mice{j},sal_day{j})
x = sal_pup(j).x;
s_on = saccadeDetector(x,fr,plt,mice{j},sal_day{j});
sal_pup(j).s_on = s_on;
function s_on = saccadeDetector(x,fr,toPlot,mouse,dt)
%set thresholds (in multiples of s.d.)
vsd = 3;%velocity
dsd =  4;%distance
%smooth x by taking moving median over 2s window
ssx = movmedian(x,5,'omitnan');

% estimate noise
dm = x-ssx;
mdm = nanmean(dm);
sdm = nanstd(dm);

dm(dm>(mdm+2*sdm) | dm<(mdm-2*sdm)) = nan;
noise_wind = 9;
noise = circshift(movstd(dm,noise_wind,'omitnan'),floor(noise_wind/2));

% velocity: px/frame
vnoise_wind = 9;
dssx = movmean([0 diff(ssx)],7);
dssx2 = dssx;
dssx2(dssx2>2*nanstd(dssx2) | dssx2<-2*nanstd(dssx2)) = nan;
vnoise = circshift(movstd(dssx2,vnoise_wind,'omitnan'),floor(vnoise_wind/2));

% putative saccades (velocity threshold)
sacs = [0; abs(dssx)>(vsd*vnoise)];%(2.5*nanstd(dssx));

% saccade onsets (not precise)
s_on = find(diff(sacs)==1);
s_off = find(diff(sacs)==-1);
if length(s_off)<length(s_on)
    s_off = [s_off;s_on(end)];
end
% add distance threshold (>1 sd away)
bads = [];
for s = 1:length(s_on)
%     span = s_on(s):s_off(s);
     if abs(ssx(s_off(s))-ssx(s_on(s))) < dsd*noise(s_on(s)) || isnan(ssx(s_on(s)-1)) || isnan(ssx(s_off(s)))
         bads = [bads s];
     end% abs(diff(x))>noise(1:end-1)
end
s_on(bads) = [];
fprintf('Found %i saccades; %i kicked off by dist thresh\n',length(s_on),length(bads))

% plot
if toPlot
    t = (1:length(x))/fr;
    
    figure('Position',[10 500 1850 400])
    
    ax1 = subplot(2,1,1);
    plot(t,x-nanmean(x)); hold on
    plot(t,ssx-nanmean(ssx),'g:','LineWidth',1.5); 
    plot(t,sacs)
    plot(t,noise*dsd,'r')
    stem(s_on/fr,ones(size(s_on)),'k')
    legend('smooth x','raw x','crossed vel thresh','noise thresh','sac onsets')
    title(sprintf('%s %s velocity thresh= %1.1f; dist thresh = %1.1f',mouse,dt,vsd,dsd))
    
    ax2 = subplot(2,1,2);
    plot(t(1:end-1),dssx,'k'); hold on
    plot(t(1:end-1),vnoise*vsd,'m',....
        t(1:end-1),vnoise *-vsd,'m')
    legend('velocity','vel thresh')
    linkaxes([ax1 ax2],'x')
end
end
function [tRun, xRun] = prepRun(speedTr)
% create timestamps
ts = 0.005 * [0:size(speedTr,1)-1];
ts = repmat(ts(:),1,size(speedTr,2));

% change timestamps so that they are relative to start of session
cts = ts + repmat(0:size(speedTr,2)-1,size(speedTr,1),1);
cts = cts(:);

%interpolate
xq = 0:.005:ceil(max(cts));
vq = interp1(cts,speedTr(:),xq);

% subsample to 20Hz
tRun = xq(1:10:end);
xRun = vq(1:10:end);
end
function plotTraces(x,y,d,mouse,dt,fr,c,tRun,xRun)
ax1 = subplot(4,1,1);
plot((1:length(x))/fr,x-nanmean(x),c)
title([mouse ' ' dt]); ylabel('x pos [pixels]')
ax2 = subplot(4,1,2);
plot((1:length(y))/fr,y-nanmean(y),c)
ylabel('y pos [pixels]')
ax3 = subplot(4,1,3);
plot((1:length(d))/fr,d,c)
ylabel('diameter [pixels]')
xlabel('time (s)')
ax4 = subplot(4,1,4);
plot(tRun,xRun,'Color',[.5 .5 .5])
ylabel('speed [cm/s]')
xlabel('time (s)')
linkaxes([ax1 ax2 ax3 ax4],'x')
end
function eye = checkTrials(mouse,dt,eye)
nt = length(eye.px);
if strcmp(mouse,'GF338') && strcmp(dt,'20190701')
    good = 1:70;
elseif strcmp(mouse,'GF338') && strcmp(dt,'20190703')
    good = 1:113;
elseif strcmp(mouse,'GF339') && strcmp(dt,'20190802')
    good = 1:50;
elseif strcmp(mouse,'GF341') && strcmp(dt,'20190828')
    good = [1:287,291:nt];
elseif strcmp(mouse,'GF343') && strcmp(dt,'20190903')
    good = [4:13,26:34,42,46:51,55-86];
elseif strcmp(mouse,'GF343') && strcmp(dt,'20190911')
    good = [1:205,209:224,226:243,245:259,261:343,345:nt];
elseif strcmp(mouse,'GF343') && strcmp(dt,'20190916')
    good = [1:107,111:189,205:226,231:300,302:nt];
elseif strcmp(mouse,'GF344') && strcmp(dt,'20190902')
    good = [1:20,22:nt];
elseif strcmp(mouse,'GF344') && strcmp(dt,'20190903')
    good = [4:32,34:76];
elseif strcmp(mouse,'GF346') && strcmp(dt,'20190926')
    good = [1:43,45:47,49,53:57,61:nt];
elseif strcmp(mouse,'GF346') && strcmp(dt,'20190927')
    good = [1:14,16:132,150:159,174:341];
else
    good = 1:nt;
end

for t = 1:nt
    if ~ismember(t,good)
        eye.px{t} = nan(60,1);
        eye.py{t} = nan(60,1);
        eye.pd{t} = nan(60,1);
    end
end

end
function tr = centerPupil(x,pthresh,fr)
%return trials when pupil is at center throughout visual stimulus
%presentation
% set # of trials
nt = 1500;

% center x
x = x-nanmean(x);

% cuts/extends to appropriate length
if length(x) <= nt*fr
    xx = nan(nt*fr,1);
    xx(1:length(x)) = x;
else
    xx = x(1:nt*fr);
end

%reshapes, removes values not during stimulus presentation
xxx = reshape(xx,fr,nt);
xxx([1:floor(0.25*fr),ceil(0.75*fr):fr],:) = nan;

%find trials where pupil x position is less than pthresh away from mean
tr = sum(abs(xxx)>pthresh)==0;
end