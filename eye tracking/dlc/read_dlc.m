% analyze results of Deep Lab Cut for eye tracking
% DLC run in rdk_2afc_eye_tracking.ipynb
% Marina March 2021

% TO DO: set threshold for points using likelihood ?0.8?
%%
dirs = {'E:\\Imaging2AFCPop\\eye tracking\\MF131_070415\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\MF122_180215\\',...%'E:\\Imaging2AFCPop\\eye tracking\\MF183\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\050516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\060616\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\070516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\080516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\100516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\110516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\190516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\230516\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\MF183_150116\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\MF183_180116\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\MF183_220116\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\MF183_230116\\',...
       'E:\\Imaging2AFCPop\\eye tracking\\MF183_260116\\'};
base = {'ec_MF131_070415_eye_FIXED_',... %03d
       'ec_MF122_180215_eye_4b_FIXED_',...%'E:\\Imaging2AFCPop\\eye tracking\\MF183\\',...
       'ec_MF209_050516_6a__FIXED_',...
       'MF209_060516_7a_FIXED_',...
       'MF209_070516_8b_FIXED_',...
       'MF209_080516_9a_FIXED_',...
       'MF209_100516_11a_FIXED_',...
       'MF209_110516_12a_FIXED_',...
       'MF209_190516_18a_FIXED_',...
       'MF209_230516_20aFIXED_',...
       'MF183_150116_Loc9a_eye_',...
       'MF183_180116_Loc11a_eye_',...
       'MF183_220116_Loc15a_eye_',...
       'MF183_230116_Loc16a_eye_',...
       'MF183_260116_Loc18a_eye_'};
%% post-processing on DLC results: load files, fit ellipse, save
for d = 1:length(dirs)% loop through sessions
    cd(dirs{d})
    files = dir([base{d} '*DLC_resnet_50_EyeTrackingMar1shuffle1_200000.csv']);
    strfmt = [base{d} '%03dDLC_resnet_50_EyeTrackingMar1shuffle1_200000.csv'];

    nt = length(files); % number of trials
    fprintf('Analyzing data from folder:\n%s\n',dirs{d})
    
    % get trial numbers from file (fixed first with fixEyeVids)
    trNum = nan(1,nt);
    for t = 1:nt % loop through trials
        trNum(t) = sscanf(files(t).name,strfmt);
    end
    lastTrial = max(trNum);
    firstTrial = min(trNum);
    
    area = {lastTrial,1}; px_fit = {lastTrial,1}; py_fit = {lastTrial,1};
    cx_fit = {lastTrial,1}; cy_fit = {lastTrial,1};
    for t = firstTrial:lastTrial % loop through trials
        fn = sprintf('%s%03dDLC_resnet_50_EyeTrackingMar1shuffle1_200000.csv',base{d},t);
        
        if mod(t,50)==0
            fprintf('fitting trial %i of %i\n',t,lastTrial)
        end
        if exist(fn,'file') %checks whether file for this trial exists
            M = readtable(fn,'headerlines',1);
            
            % extract x, y, likelihood
            %pupil
            px = cellfun(@str2double,table2array(M(2:end,2:3:23)));
            py = cellfun(@str2double,table2array(M(2:end,3:3:24)));
            pl = cellfun(@str2double,table2array(M(2:end,4:3:25)));
            
            %corneal reflection
            cx = cellfun(@str2double,table2array(M(2:end,32:3:35)));
            cy = cellfun(@str2double,table2array(M(2:end,33:3:36)));
            cl = cellfun(@str2double,table2array(M(2:end,34:3:37)));
            
            % fit ellipse
            nfr = size(px,1);
            area{t} = zeros(nfr,1);
            for f = 1:nfr % loop through frames
                e = fit_ellipse(px(f,:),py(f,:),'n');
                if isempty(e.a)
                    area{t}(f) = nan;
                    px_fit{t}(f) = nan;
                    py_fit{t}(f) = nan;
                else
                    area{t}(f) = get_area(e);
                    px_fit{t}(f) = e.X0_in;
                    py_fit{t}(f) = e.Y0_in;
                end
            end
            
            cx_fit{t} = mean(cx,2);
            cy_fit{t} = mean(cy,2);
        else
            fprintf('No eye tracking for trial %i\n',t)
        end
    end
    %save results
    save('eye_dlc_results.mat','area','px_fit','py_fit','cx_fit','cy_fit','files')
end
%% plot on image
% fn = sprintf('ec_MF209_050516_6a_%iDLC_resnet_50_EyeTrackingMar1shuffle1_200000.csv',1);
% M = readtable(fn,'headerlines',1);
% px = cellfun(@str2double,table2array(M(2:end,2:3:23)));
% py = cellfun(@str2double,table2array(M(2:end,3:3:24)));
% f = 1; %frame
% e = fit_ellipse(px(1,:),py(f,:),'n');
% vr = VideoReader('E:\Imaging2AFCPop\eye tracking\050516\ec_MF209_050516_6a_0.avi');
% eyeFr = read(vr,f);
% figure
% imagesc(eyeFr)
% axis equal; axis off
% colormap gray
% hold on
% scatter(e.X0_in,e.Y0_in,40,'rs','filled')
% plot(e.rotated_ellipse(1,:),e.rotated_ellipse(2,:))
 
function a = get_area(e)
a = pi*e.a*e.b;
end