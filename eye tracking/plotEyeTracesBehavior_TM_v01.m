function [eyeTracesSorted,eyeInfo]=plotEyeTracesBehavior_TM_v01(fTraces,behavior,results,param)

%% Start new presentation
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
warning('off','all')
mouse=behavior.header.mouse;
date=behavior.header.date;
colorMapTemporal = [12 178 232;0 255 250;12 232 165]/255;
colorMapNasal = [232 184 0;255 150 0;255 75 0]/255;
color_coher = {colorMapNasal,colorMapTemporal};

%% HEADER
% Tiago G. Marques
%fTraces: deltaF/F traces
%ROInumber: ROI number to analyze
%behavior: behavior file produced by loadBehavior_TM_v08
%Tmax: maximum trial time
%Tana: [Tmin Tmax] - to analyze traces

multFreq=10;
tWindowConcat = 100;
cleanEyeTraces=1;
eyeCalib=1/94;
d_eye=1.2;
fo_baseline=5;
sac_thr=2;

%% Inputs
% Times
tTail = -2;
Tbase=[-1,0];
Ttrial = [-1.5 3.5];
Tstim = [0 1.5];

limStim = [-0.5 2];
xTickStim = [0 1 2];

TstimPlot = [0 1.5];

rtBins = [0.5 1 1.5 2];

beThrsh=[2 1 1 2];
beNSmooth = [2 3 3 2];
beDuration = [2 1 1 1];
beBoutsName = {'speedBouts','leftLickBouts','rightLickBouts'};
nROCbins = 500;

%%
%%% DO NOT EDIT
%%

%%
%%% PREPARATIONS
%%

%% Auxiliary preparations

%% Gets parameteres from input structures
% Number of trials
numberOfTraces = behavior.header.n_trials;
numberOfEyeTraces=length(results.byTrial.nframes);
if numberOfEyeTraces<numberOfTraces
    numberOfTraces=numberOfEyeTraces;
end

% dt for traces (should be the same)
dt = behavior.header.dt;

% Parameter values
paramValues = behavior.values.(param);
nValues = length(paramValues);
% Trial times
stimOnset = behavior.times.stim_on(1:numberOfTraces);

%% Trials indices
% All trials
selectedTrials = behavior.trials.selected;
discontinuities = find(diff(selectedTrials)==1)+1;
discontinuities = [1; discontinuities];

for j=1:2 % For each side
    trialsSide{j}=behavior.trials.side==2-j;
    trialsASide{j}=circshift(trialsSide{j},1);
    trialsASide{j}(discontinuities)=0;
end

for i=1:2 % for each outcome
    trialsOut{i} = behavior.parameters.outcome==i;
end

for k =1:nValues % for each Coherence
    trialsCoher{k} = behavior.parameters.(param)==paramValues(k);
end

for p=1:2 % for each choice
    trialsChoice{p}=behavior.trials.choice==2-p;
    trialsAChoice{p}=circshift(trialsChoice{p},1);
    trialsAChoice{p}(discontinuities)=0;
end


% Removes timeouts from selected trials
%selectedTrials((behavior.parameters.outcome==3)|changeOfMind==1)=false;
selectedTrialsNum = find(selectedTrials);
nSelectedTrials=sum(selectedTrials);

% Master trials sorted
[trialsSorted, ~] = getTrialsMaster({'Side','Out','Choice','Coher'},...
    {2,2,2,nValues},{trialsSide,trialsOut,trialsChoice,trialsCoher},selectedTrials,0);

%% Trial values
trialsMasterSide=behavior.trials.side*(-2)+1;
trialsMasterOut=behavior.parameters.outcome;
trialsMasterSideSel=trialsMasterSide(selectedTrials);
trialsMasterOutSel=trialsMasterOut(selectedTrials);

%% Formats Traces
% Ca traces padded
if not(isempty(fTraces))
    [~,~,tracesLength,~,~,~,~,~] =...
        formatCaTraces(fTraces);
    tracesLength(end)=3;
    [eyeTracesPadded,~,eyeTimes] =...
        formatEyeTraces(results,multFreq,cleanEyeTraces,dt,tracesLength);
else
    [eyeTracesPadded,~,eyeTimes,tracesLength] =...
        formatEyeTraces(results,multFreq,cleanEyeTraces,dt);
end

[beTracesPadded,beTracesConcat,times,~,~] =...
    formatBeTraces_v4(behavior,0,tracesLength,beNSmooth);

trialDuration=tracesLength*dt;


% Adds tail of last trial for plotting traces
for b=1:3
    [beTracesTail{b},timeVectFull]=addTailTraces_v3(beTracesPadded{b},tTail,dt,tracesLength);
end

% eyeTracesTail{2}=addTailTracesBaseline_v3(eyeTracesPadded{2},tTail,dt,tracesLength);
% eyeTracesTail{3}=addTailTracesBaseline_v3(eyeTracesPadded{3},tTail,dt,tracesLength);

% Converts px to º in eye variables
for e=1:3
    eyeTracesTail{e}=addTailTraces_v3(eyeTracesPadded{e},tTail,dt,tracesLength);
    eyeTracesTail{e}=asin(eyeTracesTail{e}*eyeCalib)*180/pi;
end

% Windows
trialWindow = [roundn(stimOnset,2) roundn(stimOnset,2)]+repmat(Ttrial,numberOfTraces,1);
% Number of frames
trialTimeLength = floor((Ttrial(2)-Ttrial(1))/dt);

% Indices
indexTracesTrial = false(length(timeVectFull),numberOfTraces);

for i=1:numberOfTraces
    indexTracesTrial(find(timeVectFull>trialWindow(i,1),trialTimeLength,'first'),i) = true;
end
% Time values aligned with stim onset
trialTimeValues = timeVectFull(find(timeVectFull>trialWindow(1,1),trialTimeLength,'first'))-mode(roundn(stimOnset,2));
indexTracesBase=(trialTimeValues>Tbase(1))&(trialTimeValues<Tbase(2));


% Eye traces
for b=1:3
    beTrialTraces{b} = getSelTraces(beTracesTail{b},indexTracesTrial);
    beTracesSorted{b} = getTracesMaster(beTrialTraces{b},1,{'Side','Out','Choice','Coher'},...
        {2,2,2,nValues},{trialsSide,trialsOut,trialsChoice,trialsCoher},selectedTrials,0);
end


% Eye traces
for e=1:3
    eyeTrialTraces{e} = getSelTraces(eyeTracesTail{e},indexTracesTrial);
end
%% Eye events (saccades)
eyeTimesTail=addTailEvents_v1(eyeTimes,trialWindow,dt,tracesLength,mode(roundn(stimOnset,2)));
eyeTimesTailSel=eyeTimesTail(selectedTrials);

eyeEventNames=fieldnames(eyeTimesTailSel{1});
for j=1:length(eyeEventNames)
    eyeEventArray.(eyeEventNames{j})=[];
end
for i=1:length(eyeTimesTailSel);
    for j=1:length(eyeEventNames)
        eyeEventArray.(eyeEventNames{j})=[eyeEventArray.(eyeEventNames{j}); eyeTimesTailSel{i}.(eyeEventNames{j})];
    end
end
eyeTrialTracesCorr{1}=eyeTrialTraces{1};

for i=1:size(eyeTrialTraces{2},2)
    xTrace=eyeTrialTraces{2}(:,i);
    
    dXTrace=diff(xTrace);
    
    x_sac(:,i)=abs(dXTrace)>sac_thr;
    ind_x_sac=find(x_sac(:,i));
    
    for j=1:length(ind_x_sac)
        sac_jump=xTrace(ind_x_sac(j)+1)-xTrace(ind_x_sac(j));
        xTrace(ind_x_sac(j)+1:end)=xTrace(ind_x_sac(j)+1:end)-sac_jump;
    end
    eyeTrialTracesCorr{2}(:,i)=xTrace;
end
eyeTrialTracesCorr{3}=eyeTrialTraces{3};
eyeTrialTracesCorr{4}=[x_sac; zeros(1,size(x_sac,2))]*100;
for i=1:size(eyeTrialTraces{2},2)
    eyeTrialTracesCorr{2}(:,i)=eyeTrialTracesCorr{2}(:,i)-mean(eyeTrialTracesCorr{2}(indexTracesBase,i));
    eyeTrialTracesCorr{3}(:,i)=eyeTrialTracesCorr{3}(:,i)-mean(eyeTrialTracesCorr{3}(indexTracesBase,i));

    eyeTrialTraces{2}(:,i)=eyeTrialTraces{2}(:,i)-mean(eyeTrialTraces{2}(indexTracesBase,i));
    eyeTrialTraces{3}(:,i)=eyeTrialTraces{3}(:,i)-mean(eyeTrialTraces{3}(indexTracesBase,i));
end


% eyeTrialTracesCorr{1}=eyeTrialTraces{1};
% for i=1:length(eyeTimesTail)
%     if not(isempty(eyeTimesTail{i}.x_saccades))
%         tempTrace=eyeTrialTraces{2}(:,i);
%         for j=1:length(eyeTimesTail{i}.x_saccades)
%             ind= find(abs(eyeTimesTail{i}.x_saccades(j)-trialTimeValues)<0.001);
%             if ind<length(eyeTimesTail{i}.x_saccades)
%                 sac_jump=tempTrace(ind+1)-tempTrace(ind);
%                 tempTrace(ind+1:end)=tempTrace(ind+1:end)-sac_jump;
%             end
%         end
%         eyeTrialTracesCorr{2}(:,i)=tempTrace;
%     else
%         eyeTrialTracesCorr{2}(:,i)=eyeTrialTraces{2}(:,i);
%     end
% end
% eyeTrialTracesCorr{3}=eyeTrialTraces{3};

for e=1:4
    eyeTracesSorted{e} = getTracesMaster(eyeTrialTracesCorr{e},1,{'Side','Out','Choice','Coher'},...
        {2,2,2,nValues},{trialsSide,trialsOut,trialsChoice,trialsCoher},selectedTrials,0);
end

%%
%%% FIGURES PREPARATIONS COMMON
%%

%% Builds colormaps according to number of parameter values
orange_color = [255 140 0]/255;
purple_color = [75 0 130]/255;
color_choice = {'red','green'};
color_outcome = {'blue',orange_color};
color_stim = {'cyan','magenta'};

%% Legends
% Sides
sideLegend = {'270º', '90º'};
% Outcomes
outLegend = {'Correct', 'Wrong'};
% Choice
choiceLegend = {'Left','Right'};
% Parameter values
for j=1:nValues
    paramLegend{j} = num2str(paramValues(j));
end

% Behavioral parameters
beLegend = {'Speed','LL Rate','RL Rate'};
eyeLegend = {'Radius','xx','yy'};

%% X-axis limits
limTrial = [trialTimeValues(1)  trialTimeValues(end)];
paramRange = (paramValues(end)-paramValues(1));
limParam = [paramValues(1)-paramRange/10 paramValues(end)+paramRange/10];

%% Limits for traces for defining yLim for behavior plots
maxBeTraces{1} = max(cellfun(@max,beTracesSorted{1}.Choice.mean(:)));
maxBeTraces{2} = max(max(cellfun(@max, beTracesSorted{3}.Choice.mean(:))),...
    max(cellfun(@max,beTracesSorted{2}.Choice.mean(:))));
maxBeTraces{3}=maxBeTraces{2};

minBeTraces{1} = min(cellfun(@min, beTracesSorted{1}.Choice.mean(:)));
minBeTraces{2} = min(min(cellfun(@min, beTracesSorted{3}.Choice.mean(:))),...
    min(cellfun(@min,beTracesSorted{2}.Choice.mean(:))));
minBeTraces{3}=minBeTraces{2};
for b=1:3
    limBeTraces{b} = [min(-maxBeTraces{b}*0.1,minBeTraces{b}*1.1) maxBeTraces{b}*1.1];
    if isnan(limBeTraces{b})
        limBeTraces{b} =[0 1];
    end
end

for e=1:3
%     maxEyeTraces{e} = max(cellfun(@nanmax,eyeTracesSorted{e}.All.mean(:)));
%     minEyeTraces{e} = min(cellfun(@nanmin, eyeTracesSorted{e}.All.mean(:)));
%     rangeEyeTraces{e}=max(cellfun(@nanmax,eyeTracesSorted{e}.All.mean(:)))-...
%         min(cellfun(@nanmin,eyeTracesSorted{e}.All.mean(:)));
%     meanEyeTraces{e}=(max(cellfun(@nanmax,eyeTracesSorted{e}.All.mean(:)))+...
%         min(cellfun(@nanmin,eyeTracesSorted{e}.All.mean(:))))/2;
    maxEyeTraces{e} = nanmax(eyeTracesSorted{e}.All.mean);
    minEyeTraces{e} = nanmin(eyeTracesSorted{e}.All.mean);
    rangeEyeTraces{e}=nanmax(eyeTracesSorted{e}.All.mean)-...
        nanmin(eyeTracesSorted{e}.All.mean);
    meanEyeTraces{e}=(nanmax(eyeTracesSorted{e}.All.mean)+...
        nanmin(eyeTracesSorted{e}.All.mean))/2;
end
rangeEyeTraces{2}=max([rangeEyeTraces{2:3}]);
rangeEyeTraces{3}=rangeEyeTraces{2};
for e=1:3
    limEyeTraces{e} = [-rangeEyeTraces{e}*0.55 +rangeEyeTraces{e}*0.55]+meanEyeTraces{e};
end

btBin=0.5;
sac_bins=(Ttrial(1):btBin:Ttrial(2))';

x_sac_hist=histc(eyeEventArray.x_saccades,sac_bins)/btBin/sum(selectedTrials);

max_x_sac_hist=max(x_sac_hist);

figure;
subtightplot(4,1,1,[],0.1,0.1)
ylim(limBeTraces{1})
plotStimTime_v3(gca,TstimPlot)
hold on
plotMeanFandSE(beTracesSorted{1}.Choice.mean{1},beTracesSorted{e}.Choice.errs{1},trialTimeValues,color_choice{1},'-',1);
plotMeanFandSE(beTracesSorted{1}.Choice.mean{2},beTracesSorted{e}.Choice.errs{2},trialTimeValues,color_choice{2},'-',1);
ylabel(beLegend{1},'FontSize',12)
hold off
box off
set(gca,'xtick',[])
set(gca,'XLim',limTrial);
ylim(limBeTraces{1})


for e=1:3
    subtightplot(4,1,e+1,[],0.1,0.1)
    ylim(limEyeTraces{e})
    plotStimTime_v3(gca,TstimPlot)
    hold on
    plotMeanFandSE(eyeTracesSorted{e}.Choice.mean{1},eyeTracesSorted{e}.Choice.errs{1},trialTimeValues,color_choice{1},'-',1);
    plotMeanFandSE(eyeTracesSorted{e}.Choice.mean{2},eyeTracesSorted{e}.Choice.errs{2},trialTimeValues,color_choice{2},'-',1);
    ylabel(eyeLegend{e},'FontSize',12)
    hold off
    box off
    if e<3
        set(gca,'xtick',[])
    else
        xlabel('Time (s)','FontSize',12)
    end
    set(gca,'XLim',limTrial);
    ylim(limEyeTraces{e})
end

figure;
subtightplot(4,1,1,[],0.1,0.1)
ylim(limBeTraces{1})
plotStimTime_v3(gca,TstimPlot)
hold on
plotMeanFandSE(beTracesSorted{1}.SideCoher.mean{1,3},beTracesSorted{e}.SideCoher.errs{1,3},trialTimeValues,color_stim{1},'-',1);
plotMeanFandSE(beTracesSorted{1}.SideCoher.mean{2,3},beTracesSorted{e}.SideCoher.errs{2,3},trialTimeValues,color_stim{2},'-',1);
ylabel(beLegend{1},'FontSize',12)
hold off
box off
set(gca,'xtick',[])
set(gca,'XLim',limTrial);
ylim(limBeTraces{1})

for e=1:3
    subtightplot(4,1,e+1,[],0.1,0.1)
    ylim(limEyeTraces{e})
    plotStimTime_v3(gca,TstimPlot)
    hold on
    plotMeanFandSE(eyeTracesSorted{e}.SideCoher.mean{1,3},eyeTracesSorted{e}.SideCoher.errs{1,3},trialTimeValues,color_stim{1},'-',1);
    plotMeanFandSE(eyeTracesSorted{e}.SideCoher.mean{2,3},eyeTracesSorted{e}.SideCoher.errs{2,3},trialTimeValues,color_stim{2},'-',1);
    ylabel(eyeLegend{e},'FontSize',12)
    hold off
    box off
    if e<3
        set(gca,'xtick',[])
    else
        xlabel('Time (s)','FontSize',12)
    end
    set(gca,'XLim',limTrial);
    ylim(limEyeTraces{e})
end

figure;
for e=1:3
    subtightplot(3,1,e,[],0.1,0.1)
    ylim(limEyeTraces{e})
    plotStimTime_v3(gca,TstimPlot)
    hold on
    plotMeanFandSE(eyeTracesSorted{e}.Out.mean{1},eyeTracesSorted{e}.Out.errs{1},trialTimeValues,color_outcome{1},'-',1);
    plotMeanFandSE(eyeTracesSorted{e}.Out.mean{2},eyeTracesSorted{e}.Out.errs{2},trialTimeValues,color_outcome{2},'-',1);
    ylabel(eyeLegend{e},'FontSize',12)
    hold off
    box off
    if e<3
        set(gca,'xtick',[])
    else
        xlabel('Time (s)','FontSize',12)
    end
    set(gca,'XLim',limTrial);
    ylim(limEyeTraces{e})
end

% figure
% ylim([0 max_x_sac_hist*1.2])
% xlim(Ttrial)
% plotStimTime_v3(gca,[0 1.5]);
% hold on
% plot(sac_bins(1:(end-1))+dt/2,x_sac_hist(1:(end-1)));
% legend boxoff
% ylabel('Saccades (Hz)')
% xlabel('Time (s)')

figure
ylim([0 length(eyeTimesTail)])
plotStimTime_v3(gca,[0 1.5])
hold on
for i=1:length(eyeTimesTail)
    for j=1:length(eyeTimesTail{i}.x_saccades)
        plot([eyeTimesTail{i}.x_saccades(j) eyeTimesTail{i}.x_saccades(j)], [i-1 i],'-g','LineWidth',2);
    end
    for j=1:length(eyeTimesTail{i}.y_saccades)
        plot([eyeTimesTail{i}.y_saccades(j) eyeTimesTail{i}.y_saccades(j)], [i-1 i],'-r','LineWidth',2);
    end
end
ylabel('# Trials')
xlabel('Time (s)')
legend boxoff
title('Saccade raster plot')

figure;
ylim([-10 10])
plotStimTime_v3(gca,[0 1.5])
hold on
for j=1:2
    plotMeanSeTrace(eyeTracesSorted{2}.SideCoher.mean{j,3},eyeTracesSorted{2}.SideCoher.errs{j,3},trialTimeValues,color_coher{j}(end,:),1,'-');
    plot([0 1.5],[0 ((1-j)*2+1)*25*1.5],'--','color',color_coher{j}(end,:),'LineWidth',2);
end
ylabel('x position (º)','FontSize',14)
xlabel('Time (s)','FontSize',14)
hold on
plot(Ttrial,[0 0],'--k','LineWidth',1);
xlabel('Time (s)','FontSize',14)
set(gca,'XLim',Ttrial);
box off
hold off

figure;
ylim([-10 10])
plotStimTime_v3(gca,[0 1.5])
hold on
for j=1:2
    plotMeanSeTrace(eyeTracesSorted{3}.SideCoher.mean{j,3},eyeTracesSorted{3}.SideCoher.errs{j,3},trialTimeValues,color_coher{j}(end,:),1,'-');
    plot([0 1.5],[0 ((1-j)*2+1)*25*1.5],'--','color',color_coher{j}(end,:),'LineWidth',2);
end
ylabel('y position (º)','FontSize',14)
xlabel('Time (s)','FontSize',14)
hold on
plot(Ttrial,[0 0],'--k','LineWidth',1);
xlabel('Time (s)','FontSize',14)
set(gca,'XLim',Ttrial);
box off
hold off

figure;
ylim([-0.2 25])
plotStimTime_v3(gca,[0 1.5])
hold on
for j=1:2
    plotMeanSeTrace(eyeTracesSorted{4}.SideCoher.mean{j,3},eyeTracesSorted{4}.SideCoher.errs{j,3},trialTimeValues,color_coher{j}(end,:),1,'-');
end
ylabel('x saccades (%)','FontSize',14)
xlabel('Time (s)','FontSize',14)
hold on
plot(Ttrial,[0 0],'--k','LineWidth',1);
xlabel('Time (s)','FontSize',14)
set(gca,'XLim',Ttrial);
box off
hold off

figure;
ylim([10 30])
plotStimTime_v3(gca,[0 1.5])
hold on
for j=1:2
    plotMeanSeTrace(eyeTracesSorted{1}.SideCoher.mean{j,3},eyeTracesSorted{1}.SideCoher.errs{j,3},trialTimeValues,color_coher{j}(end,:),1,'-');
end
ylabel('Pupil size (%)','FontSize',14)
xlabel('Time (s)','FontSize',14)
hold on
plot(Ttrial,[0 0],'--k','LineWidth',1);
xlabel('Time (s)','FontSize',14)
set(gca,'XLim',Ttrial);
box off
hold off

figure;
for j=1:2
    subtightplot(2,1,j,[],0.1,0.1)
    ylim([-5 5])
    plotStimTime_v3(gca,[0 1.5])
    hold on
    plot(trialTimeValues,eyeTrialTraces{2}(:,trialsSorted.SideCoher{j,end}),'-','color',color_coher{j}(end,:));
    plot(trialTimeValues,eyeTracesSorted{2}.SideCoher.mean{j,3},'-','color','k','LineWidth',2);
    %plot(trialTimeValues,eyeTrialTraces{2}(:,102),'-','color','y');
    plot([0 1.5],[0 ((1-j)*2+1)*25*1.5],'--','color',color_coher{j}(end,:),'LineWidth',2);
    plot(Ttrial,[0 0],'--k','LineWidth',1);
    set(gca,'XLim',Ttrial);
    box off
    hold off
    ylabel('x position (º)','FontSize',14)
    if j==1
        set(gca,'xtick',[])
    else
        xlabel('Time (s)','FontSize',14)
    end
end

figure;
for j=1:2
    subtightplot(2,1,j,[],0.1,0.1)
    ylim([-10 10])
    plotStimTime_v3(gca,[0 1.5])
    hold on
    plot(trialTimeValues,eyeTrialTracesCorr{2}(:,trialsSorted.SideCoher{j,end}),'-','color',color_coher{j}(end,:));
    plot([0 1.5],[0 ((1-j)*2+1)*25*1.5],'--','color',color_coher{j}(end,:),'LineWidth',2);
    plot(Ttrial,[0 0],'--k','LineWidth',1);
    set(gca,'XLim',Ttrial);
    box off
    hold off
    ylabel('x position (º)','FontSize',14)
    if j==1
        set(gca,'xtick',[])
    else
        xlabel('Time (s)','FontSize',14)
    end
end

eyeInfo.trialTimeValues=trialTimeValues;
% figure;
% ylim([-25 25])
% plotStimTime_v3(gca,[0 1.5])
% hold on
% plot(trialTimeValues,eyeTrialTracesCorr{2}(:,35),'-','color',color_coher{2}(end,:));
% plot(trialTimeValues,eyeTrialTracesCorr{2}(:,85),'-','color',color_coher{1}(end,:));
% plot([0 1.5],[0 25*1.5],'-','color',color_coher{1}(end,:));
% plot([0 1.5],[0 -25*1.5],'-','color',color_coher{2}(end,:));
% ylabel('x position (º)','FontSize',14)
% xlabel('Time (s)','FontSize',14)
% hold on
% plot(Ttrial,[0 0],'--k','LineWidth',1);
% xlabel('Time (s)','FontSize',14)
% set(gca,'XLim',Ttrial);
% box off
% hold off


