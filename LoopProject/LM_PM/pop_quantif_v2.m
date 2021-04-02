% run "gui.m" first to generate ops_analysis
function pop_quantif_v2(ops_analysis)
global MainDir


%% Preambule
color_bp = [0.8500 0.3250 0.0980];
color_bn = [0      0.4470 0.7410];

magn_th      = ops_analysis.resp_th;
saveFig      = ops_analysis.saveFlag;
alphaVal     = ops_analysis.alpha;
tResp        = ops_analysis.tResp;
r2_th        = ops_analysis.rsq_th;
rebuild_data = ops_analysis.redo_data;
redo_selec   = ops_analysis.selection;

tv = datestr(now, 'yyyy_mm_dd');
fprintf('\n')

% if ~exist(MainDir)
MainDir = uigetdir('E:\Loop\L2','Master directory');
% else
%     MainDir = uigetdir(MainDir,'Main position directory');
% end
saveDir_plots = [MainDir filesep 'plot'];
if ~exist(saveDir_plots,'dir')
    mkdir(saveDir_plots)
end
oldfolder = cd([MainDir]);
diary off; diary(['PopQuantif_v2_' tv]);
cd(oldfolder)

%% Load data
if rebuild_data
    
    LMdir = dir([MainDir filesep 'LM']);
    str = sortrows({LMdir.name}');
    [s,~] = listdlg('PromptString','Select animals to analyze:', 'OKString', 'OK',...
        'SelectionMode','multiple',...
        'ListString', str, 'Name', 'Select a File');
    LManimals = str(s);
    
    PMdir = dir([MainDir filesep 'PM']);
    str = sortrows({PMdir.name}');
    [s,~] = listdlg('PromptString','Select animals to analyze:', 'OKString', 'OK',...
        'SelectionMode','multiple',...
        'ListString', str, 'Name', 'Select a File');
    PManimals = str(s);
    
    n_animals_LM = length(LManimals);
    n_animals_PM = length(PManimals);
    
    clear s v str
    
    disp('FORCED to rebuild data structure')
    % - LM
    fprintf(2,'LM animals\n')
    for i = 1:n_animals_LM
        disp(LManimals{i})
        pos_directories = dir([LMdir(1).folder filesep LManimals{i}]);
        match = {}; n_pos_LM(i) = 0;
        for ii = 1:length(pos_directories)
            if ~isempty(regexp(pos_directories(ii).name,'pos','ONCE'))
                match{end+1} = pos_directories(ii).name;
                n_pos_LM(i) = n_pos_LM(i) +1;
                index(i,n_pos_LM(i)) = ii;
            end
        end
        for ii = 1:n_pos_LM(i)
            path = [pos_directories(index(i,ii)).folder filesep pos_directories(index(i,ii)).name];
            disp(pos_directories(index(i,ii)).name)
            file = dir(fullfile([path filesep], '*.mat'));
            temp=load([file.folder filesep file.name]);
            data_LM{i,ii} = temp.data.data_sorted;
            info_LM{i,ii} = temp.data.Info;
            info_LM{i,ii}.beads_pos = temp.data.beads_pos;
            info_LM{i,ii}.pos = file;
            
        end
        %     animal_dir{i} = pos_directories;
    end
    
    % - PM
    fprintf(2,'PM animals\n')
    for i = 1:n_animals_PM
        disp(PManimals{i})
        pos_directories = dir([PMdir(1).folder filesep PManimals{i}]);
        match = {}; n_pos_PM(i) = 0;
        for ii = 1:length(pos_directories)
            if ~isempty(regexp(pos_directories(ii).name,'pos','ONCE'))
                match{end+1} = pos_directories(ii).name;
                n_pos_PM(i) = n_pos_PM(i) +1;
                index(i,n_pos_PM(i)) = ii;
            end
        end
        for ii = 1:n_pos_PM(i)
            path = [pos_directories(index(i,ii)).folder filesep pos_directories(index(i,ii)).name];
            disp(pos_directories(index(i,ii)).name)
            file = dir(fullfile([path filesep], '*.mat'));
            temp=load([file.folder filesep file.name]);
            data_PM{i,ii} = temp.data.data_sorted;
            info_PM{i,ii} = temp.data.Info;
            info_PM{i,ii}.beads_pos = temp.data.beads_pos;
            info_PM{i,ii}.pos = file;
        end
        %     animal_dir{i} = pos_directories;
    end
    clear temp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save data structures
    disp('saving data structures ...')
    save([LMdir(1).folder filesep 'dataLM_' tv], 'data_LM');
    save([PMdir(1).folder filesep 'dataPM_' tv], 'data_PM');
    save([LMdir(1).folder filesep 'infoLM_' tv], 'info_LM');
    save([PMdir(1).folder filesep 'infoPM_' tv], 'info_PM');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    disp('load data structures')
    LMdir = dir([MainDir filesep 'LM']);
    PMdir = dir([MainDir filesep 'PM']);
    LM_fileName = uigetfile([LMdir(1).folder],'Select LM data and info','MultiSelect','on');
    PM_fileName = uigetfile([PMdir(1).folder],'Select PM data and info','MultiSelect','on');
    
    LM_data = LM_fileName{1};
    LM_info = LM_fileName{2};
    temp = load([LMdir filesep LM_data]);
    data_LM = temp.data_LM;
    temp = load([LMdir filesep LM_info]);
    info_LM = temp.info_LM;
    
    PM_data = PM_fileName{1};
    PM_info = PM_fileName{2};
    temp = load([PMdir filesep PM_data]);
    data_PM = temp.data_PM;
    temp = load([PMdir filesep PM_info]);
    info_PM = temp.info_PM;
    clear temp
end


%% Select responsive ROIs
if redo_selec
    fprintf('FORCED to redo data selection with alpha = %g and amp_th = %g, window(s): %g - %g\n',...
        alphaVal,magn_th,tResp(1),tResp(2))
    
    n_pos{1} = sum(~cellfun('isempty',data_LM'));
    n_pos{2} = sum(~cellfun('isempty',data_PM'));
    options = fitoptions('gauss2');
    options.Lower = [0 0 0 0 180 0];
    options.Upper = [Inf 0 Inf Inf 180 Inf];
    directions = -90:45:225;
    
    for j = 1:2 % for LM and PM data
        disp(['dataset' num2str(j)])
        
        medianResp = cell(size(n_pos{j},2),max(n_pos{j}));
        h_all = cell(size(n_pos{j},2),max(n_pos{j}));
        tuned = cell(size(n_pos{j},2),max(n_pos{j}));
        bestSF = cell(size(n_pos{j},2),max(n_pos{j}));
        beads_pos = cell(size(n_pos{j},2),max(n_pos{j}));
        for i = 1:size(n_pos{j},2) % for all animals
            disp(['animal' num2str(i)])
            for ii = 1:n_pos{j}(i) % for all positions
                if j == 1
                    data = data_LM;
                    info = info_LM;
                else
                    data = data_PM;
                    info = info_PM;
                end
                nROIs = size(data{i,ii},1);
                h_all{i,ii} = false(2,nROIs);
                tuned{i,ii} = false(2,nROIs);
                bestSF{i,ii} = false(2,nROIs);
                
                beads_pos{i,ii} = info{i,ii}.beads_pos;
                
                nSF = size(info{i,ii}.sfValues,2);
                nDir = size(info{i,ii}.directions,2);
                t_resp = false(1,size(data{i,ii},2));
                t_base = false(1,size(data{i,ii},2));
                t_resp(info{i,ii}.timeVect(info{i,ii}.xvect)>tResp(1) & info{i,ii}.timeVect((info{i,ii}.xvect))<info{i,ii}.xstim_vis(2)+tResp(2))=true;
                t_base(info{i,ii}.timeVect(info{i,ii}.xvect)<info{i,ii}.xstim_vis(1)) = true;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Vis_Resp = squeeze(mean(data{i,ii}(:,t_resp,2:9,:,1,:),2)); % non-opto trials and vis stim, to check
                Vis_Base = squeeze(mean(data{i,ii}(:,t_base,2:9,:,1,:),2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if nSF == 1
                    medianVis_Resp = median(Vis_Resp,3);
                    for iii = 1:nROIs
                        [m,idx]=max(medianVis_Resp(iii,:));
                        if m > magn_th
                            h=ttest(Vis_Resp(iii,idx,:),Vis_Base(iii,idx,:),alphaVal);
                        end
                        if h
                            Dir_sorted = [mod(idx-2,8) mod(idx-1,8) idx mod(idx+1,8) mod(idx+2,8) mod(idx+3,8) mod(idx+4,8) mod(idx+5,8)];
                            Dir_sorted(Dir_sorted == 0) = 8;
                            
                            
                            y1 = double(squeeze(medianVis_Resp(iii,Dir_sorted)));
                            y1 = y1-min(y1);
                            [~, gof1] = fit(directions.',y1','gauss2',options);
                            if gof1.rsquare > r2_th
                                tuned{i,ii}(1,iii) = true;
                            end
                        end
                        h_all{i,ii}(1,iii) = h;
                        %                     idx_LM{i,ii}(1,iii) = idx;
                    end
                    medianResp{i,ii} = mean(mean(data{i,ii}(:,t_resp,[1 Dir_sorted+1],:,:,:),2),6);
                    medianResp{i,ii}(:,:,:,2,:) = NaN(nROIs,1,9,1,2);
                    medianResp{i,ii} = squeeze(medianResp{i,ii});
                else
                    medianVis_Resp = median(Vis_Resp,4);
                    for iii = 1:nROIs
                        for iiii = 1:nSF
                            [m,idx]=max(medianVis_Resp(iii,iiii));
                            if m > magn_th
                                h=ttest(Vis_Resp(iii,idx,:),Vis_Base(iii,idx,:),alphaVal);
                            end
                            h_all{i,ii}(iiii,iii) = h;
                            if h
                                Dir_sorted = [mod(idx-2,8) mod(idx-1,8) idx mod(idx+1,8) mod(idx+2,8) mod(idx+3,8) mod(idx+4,8) mod(idx+5,8)];
                                Dir_sorted(Dir_sorted == 0) = 8;
                                
                                
                                y1 = double(squeeze(medianVis_Resp(iii,Dir_sorted)));
                                y1 = y1-min(y1);
                                [~, gof1] = fit(directions.',y1','gauss2',options);
                                if gof1.rsquare > r2_th
                                    tuned{i,ii}(iiii,iii) = true;
                                end
                            end
                            medianResp{i,ii} = squeeze(mean(mean(data{i,ii}(:,t_resp,[1 Dir_sorted+1],iiii,:,:),2),6));
                        end
                        Best_SFTF= reshape(medianVis_Resp,nROIs,nDir*nSF);
                        [~,idx_SF] = max(Best_SFTF);
                        if idx_SF < 9
                            bestSF{i,ii}(1,iii) = true;
                        else
                            bestSF{i,ii}(2,iii) = true;
                        end
                    end
                end
            end
            
            
            % if saveFig
            % %     saveas(gcf,[saveDir filesep tv '_ImpactOfLaserSpont_Diff_Bar'],'fig');
            %     screen2png([saveDir filesep tv '_Adaptation']);
            % end
            %%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save data structures
        selection.medianResp = medianResp;
        selection.bestSF = bestSF;
        selection.tuned = tuned;
        selection.h_all = h_all;
        selection.bead_pos = beads_pos;
        
        disp('saving ROI selection data ...')
        if j == 1
            save([LMdir(1).folder filesep 'selectionLM_' tv], 'selection');
        elseif j == 2
            save([PMdir(1).folder filesep 'selectionPM_' tv], 'selection');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
diary off;
end