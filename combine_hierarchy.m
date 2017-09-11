function [hier_json_file,hier_xls_file] = combine_hierarchy(stats_json_file,study_info)
% TO DO
%
%

%%%Load json
[json_path,json_fn,~]=fileparts(stats_json_file);
js = loadjson(stats_json_file);
js_data = js{1};
js_data = [js_data{:}];

%%%%List regions
info = loadjson(study_info);
if isfield(info,'hier_dir')
    list_reg_dir = info.hier_dir;
    list_reg = dir([list_reg_dir '\ABAHier2015_*']);
else
    error('combine_hierarchy:MissingField',...
        'The json file %s is missing the hier_dir field',...
        study_info);
end
fprintf(1,'Creating %d hierarchies\n',length(list_reg));
for iL = 1:length(list_reg)
    idx = xlsread(fullfile(list_reg_dir,list_reg(iL).name));
    iC=0;
    idx_sub = [];
    for iI = 1 : length(idx)
        if ~isempty(find([js_data(:).reg_idx]==idx(iI),1))
            iC=iC+1;
            idx_sub(iC) = idx(iI);
            structure(iC)=js_data(([js_data(:).reg_idx]==idx(iI)));
        end
    end
    if ~exist('structure','var')
        warning('The region indexed %d could not be found in the file %s. Skipping...\n', idx(iI),stats_json_file)
    else
        hier_region(iL).reg_name = list_reg(iL).name(strfind(list_reg(iL).name,'_')+1:end-5);
        hier_region(iL).reg_idx_full  = sprintfc('%d',idx);
        hier_region(iL).reg_idx  = sprintfc('%d',idx_sub);
        hier_region(iL).reg_pxl = sum([structure(:).reg_pxl]);
        hier_region(iL).reg_area = sum([structure(:).reg_area]);
        hier_region(iL).reg_area_unit = structure(end).reg_area_unit;
        hier_region(iL).obj_cnt = sum([structure(:).obj_cnt]);
        hier_region(iL).obj_reg_cnt_ratio = hier_region(iL).obj_cnt ./ hier_region(iL).reg_area;
        hier_region(iL).obj_pxl = sum([structure(:).obj_pxl]);
        hier_region(iL).obj_area = sum([structure(:).obj_area]);
        hier_region(iL).obj_area_unit = structure(end).obj_area_unit;
        hier_region(iL).obj_reg_area_ratio = hier_region(iL).obj_area ./ hier_region(iL).reg_area;
        %
        clear structure
    end
end

%% Write results
hier_json_file = fullfile(json_path,sprintf('%s_hier.json',json_fn));
savejson('',hier_region,hier_json_file);
plaques = struct2table(hier_region);
hier_xls_file = fullfile(json_path,sprintf('%s_hier.xlsx',json_fn));
writetable(plaques,hier_xls_file);
% writetable(regions,fullfile(output_dir,'crossseeds_m287_regions_area_122.xlsx'));