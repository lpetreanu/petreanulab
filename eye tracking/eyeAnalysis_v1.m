%whatDoesTheEyeSay?

load('eye_fixedTrialNums.mat')

nt = size(eye.pd,1);
tstimstart = round(.2*45);
tstimend = round(1.7*45);
for i= 1:nt
    if ~isempty(eye.pd{i})
        nfr(i) = size(eye.pd{i},1);
        avgDiamWT(i) = nanmean(eye.pd{i}); %average diameter in whole trial
        avgDiamST(i) = nanmean(eye.pd{i}(tstimstart:tstimend)); %average diameter during stimpresentation
        traces(i,:) = eye.pd{i}(1:230);
    end
end

load('behavior.mat')
choice = behavior.trials.choice; %1: left; 0: right
correct = behavior.trials.correct;
wrong = behavior.trials.wrong;

%% plot some things
figure,
title('Traces of pupil diameter, all trials and mean')
plot((1:230)/45,traces,'Color',[.9 .9 .9])
hold on;
plot((1:230)/45,nanmean(traces),'k','LineWidth',2)

k  = 111;
figure,
subplot(2,1,1)
imagesc(mROI(k).tracesOLclean)
hold on
subplot(2,1,2)
imagesc(traces)

