% correct eye tracking movies for skipped triggers
% load F

frMult = 5; %frame multuplier constant for eye videos. 3 up to MF183; 5 after
sessDir = {'F:\MF131\2P\070415','F:\MF132\2P\060415','F:\MF132\2P\020515',...
    'f:\MF132\2P\150415','F:\MF122\2P\180215','F:\MF122\2P\100315',...
    'F:\MF183\150116','F:\MF183\180116','F:\MF183\220116',...
    'F:\MF183\230116','F:\MF183\250116','F:\MF183\260116',...
    'F:\MF209\2P\050516','F:\MF209\2P\060516','F:\MF209\2P\070516',...
    'F:\MF209\2P\080516','F:\MF209\2P\100516','F:\MF209\2P\110516',...
    'F:\MF209\2P\140516','F:\MF209\2P\190516','F:\MF209\2P\230516',};


 cd(sessDir{21})
load('FAcrossTrials_20a.mat')
% F = FAcrossTrials_Reg_1MF209_050516_6a_;

nfr = [];
for i = 1:length(F.data)
    nfr(i) = size(F.data{i},2);
end

cd F:\MF209\2P\EyeVideos\230516
base = 'MF209_230516_20a';
fn = ls([base '*.avi']); 
eyenfr = [];
for i = 1:size(fn,1)
    mov = VideoReader(sprintf([base '%d.avi'],i-1));
    eyenfr(i) = mov.Duration*mov.FrameRate;
end

%% check frame lengths 
eyenfrtmp = eyenfr;
for i = 1:length(nfr)
    if trialsSkipped
        
        % somehow find closest trial...
        
    end
    frdif = abs(nfr(i)-eyenfrtmp(i)/frMult);
    if nfr(i) ~=0
        if i < length(nfr)
            if frdif > frMult % check if difference between # of frames is different between scanimage and eye tracking file
                disp(['Eye tracking file ' num2str(i) ' is too long'])
                frdif2 = abs((nfr(i) + nfr(i+1))-eyenfrtmp(i)/frMult); %check if eye movie nfr is explained by 2 trials
                if frdif2 > frMult
                    disp(['Eye tracking file ' num2str(i) ' is longer than 2 trials'])
                    frdif3 = ((nfr(i) + nfr(i+1) + nfr(i+2))-eyenfrtmp(i)/frMult); %check if eye movie nfr is explained by 2 trials
                    if frdif3 < -1*frMult
                        disp(['Eye tracking file ' num2str(i) ' is longer than 3 trials'])
                    else
                        eyenfrtmp(i+3:end+2) = eyenfrtmp(i+1:end);
                        eyenfrtmp(i+1:i+2) = nan;
                    end
                else
                    eyenfrtmp(i+2:end+1) = eyenfrtmp(i+1:end);
                    eyenfrtmp(i+1) = nan;
                end
            end
        end
        trialsSkipped = 0;
    else
        trialsSkipped = 1;
    end
end




%% change file names of eye tracking files to match trial #

goodTrials = find(~isnan(eyenfrtmp));
goodTrials = 186:288;
fixCsv = 1;
for i = 1:length(goodTrials)
    newFn = sprintf([base '_FIXED_%03d.avi'],goodTrials(i));
    movefile(sprintf([base '%d.avi'],i-1),newFn);
    if fixCsv
        movefile(sprintf([base '%d.csv'],i-1),sprintf([base '_FIXED_%03d.csv'],goodTrials(i)));
    end
end