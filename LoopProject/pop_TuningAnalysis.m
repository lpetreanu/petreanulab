% pop_TuningAnalysis({data_18_p1,data_18_p2,data_18_p3_1,data_18_p3_2,data_18_p3_3,data_18_p3_4,data_19_p1,data_19_p2,data_19_p3,data_19_p4,data_19_p5_1,data_19_p5_2,data_19_p6_1,data_19_p6_2,data_22_p1_1,data_22_p1_2,data_22_p1_3,data_22_p1_4,data_22_p2},1,1)
function TuningData = pop_TuningAnalysis(data_list,do_MF,looped_cells,saveFig)
% looped_cells = 0 if bead+ cells have not been systematically identified in the dataset 

for i = 1:length(data_list)
    nSF(i) = length(data_list{i}.Info.sfValues);
    nDir(i) = length(data_list{i}.Info.directions);
    nROIs(i) = size(data_list{i}.dir_tuning,1);
    if looped_cells
    nBp(i) = size(data_list{i}.beads_pos,1);
    end
end
disp(['Data from ' num2str(i) ' sessions'])
[max_nSF,I]=max(nSF);
k1 = find(nSF~=max_nSF);
% k2 = find(nDir~=max(nDir));
if ~isempty(k1)
    disp('All sessions do not have the same parameter space')
    for ii = 1:length(k1)
        data_list{k1(ii)}.dir_tuning = repmat(data_list{k1(ii)}.dir_tuning,1,1,max(nSF),1);
        data_list{k1(ii)}.h_vis(:,2) = false(size(data_list{k1(ii)}.h_vis,1),1);
%         data_list{k1(ii)}.h_vis = repmat(data_list{k1(ii)}.h_vis,1,max(nSF));
    end
else
    disp('All sessions have the same parameter space')
end

if looped_cells == 1
    try % legacy version (before TG L5 data - June 2020)
        for i = 1:length(data_list)
            bp_temp{i} = false(nROIs(i),1);
            bp_temp{i}(data_list{i}.beads_pos(:,3)) = true;
        end
    catch
        for i = 1:length(data_list)
            bp_temp{i} = data_list{i}.beadID.bead_pos';
        end
    end
else
    for i = 1:length(data_list)
        bp_temp{i} = true(nROIs(i),1);
        if isfield(data_list{i},'beads_pos')
        data_list{i} = rmfield(data_list{i},{'beads_pos','beadID'});
        end
    end
    
end
bp_temp = cat(1,bp_temp{:});
beads_pos = find(bp_temp);
beads_pos = repmat(beads_pos,1,3);

%         for i = 1:length(data_list)
% fit_temp{i} = data_list{i}
%             
%         end
        
temp = cat(1,data_list{:});
data_all.dir_tuning = cat(1,temp(:).dir_tuning);
data_all.fit = cat(1,temp(:).fit);
try % legacy version (before TG L5 data - June 2020)
    data_all.h_vis = cat(1,temp(:).h_vis);
catch
%     h_vis = cell(length(temp),1);
    for session = 1:length(data_list)
        h_vis{session} = temp(session).fit.h_vis;
    end
    data_all.h_vis = cat(1,h_vis{:});
end
data_all.beads_pos = beads_pos;
data_all.Info.sfValues = data_list{I}.Info.sfValues;
data_all.Info.tfValues = data_list{I}.Info.tfValues;

% data_all.MoreInfo = cat(1,temp(:).Info);

TuningData = TuningAnalysis_v2(data_all,do_MF,looped_cells,saveFig);
end