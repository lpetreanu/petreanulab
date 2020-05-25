%process batch of eye files, save in structure called eye.

%% set analysis 
%displayFlag (only for pupil)
% 0 : don't display or save resulted fit as movie (default)
% 1 : display and save movie with fitted pupil
% 2 : save, but don't display movie with fitted pupil
%     movie saved as 'fit_ + filename'
displayFlag = 2;
toAnalyze = 'pupil';%'pupil' %whether to track pupil or corneal reflection

warning('off','MATLAB:hardcopy:DeprecatedHardcopyFunction')
warning('off','MATLAB:print:DeprecateZbuffer')

%% read in files to analyze
[exfn,path]=uigetfile('.avi','Select one of the files to analyze');
cd(path)
 d = dir('*.avi'); %get list of all .avi files in directory
%  allfiles = sortrows({d.name}');
allfiles = {d.name};
 
[s,v] = listdlg('PromptString','Select files to analyze:', 'OKString', 'OK',...
    'SelectionMode','multiple', 'Name', 'Select a File','ListString',allfiles);
names = allfiles(s);
nfiles = length(names);

%% set parameters
[paramsfn]=uigetfile('.mat','Select file with algorithm parameters [Cancel if none]');
if ~(paramsfn)
    params.pMax = 40; %max radius of pupil(in pixels) (length of starburst rays)
    params.eyeInit = [75 57]; %set initial position of pupil at beginning of session by eye
    params.enhance = 1; %set to 1 to boost contrast; otherwise, leave 0
    params.startFrame = 1;
    params.sthresh = .02;  %threshold for sobel filter; higher means more restrictive (more restrictive means it finds less edges)
    params.pthresh = 15;  %percent threshold for darkest pixels (to count as pupil)
    params.cfthresh = 100; %threshold for points to be counted as edge of corneal reflection
    params.mvthresh = 10;   %max movement in pixels allowed from frame to frame
    params.deltapdthresh = 5; % max change in pupil diameter allowed frame to frame
    params.filtSize = 2; % size of circular averaging filter
    params.s = 'c';
    save([toAnalyze '_params.mat'],'params')
else
    load([path filesep paramsfn])
end

%% run the eyetracker

%initialize some variables
switch toAnalyze
    case 'pupil'
        eye.px = cell(nfiles,1);
        eye.py = cell(nfiles,1);
        eye.pd = cell(nfiles,1);
    case 'CR'
        cr.cx = cell(nfiles,1);
        cr.cy = cell(nfiles,1);
end

%loop through files
for i = 1:nfiles
    disp(names{i})
    % trialNum assumes file name has been fixed to match trial number (padded with 0, first trial is 1)
    trialNum = str2double(names{i}(end-6:end-4));
   
    switch toAnalyze
        case 'pupil'
            if i > 1 && ~isnan(eye.px{trialNum-1}(end)) % use eye position from end of previous trial
                lastfitind = find(~isnan(eye.px{trialNum-1}),1,'last'); %find most recent fit
                params.eyeInit = round([eye.px{trialNum-1}(lastfitind) eye.py{trialNum-1}(lastfitind)]);
            end
            [ pd,px,py ] = analyzeThatPupil( params,[path filesep names{i}],displayFlag);
        case 'CR'
            [ cx,cy ] = analyzeThatCR( params,[path filesep names{i}],0);
    end
    
    %save output to structure
    switch toAnalyze
        case 'pupil'
            eye.px{trialNum} = px(:,2); eye.py{trialNum} = py(:,2); eye.pd{trialNum} = pd;
        case 'CR'
            cr.cx{trialNum} = cx(:,2); cr.cy{trialNum} = cy(:,2);
    end
end
%% save stuff
switch toAnalyze
    case 'pupil'
        save('pupil_data.mat','eye')
    case 'CR'
        save('cr_data.mat','cr')
end