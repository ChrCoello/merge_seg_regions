function [stats_json_filename,stats_xls_filename] = combine_obj_reg(varargin)
%COMBINE_OBJ_SEG Combines objects from same regions
% Takes as input JSON file saved by quantifyFullDataset
%  - List individual regions in the JSON file and calculate total area over
%  all the sections of each region
%  - List individual regions where there are objects and in each of these
%  regions count how many objects there are (and other characteristics)
% 
%
% CC

% Parse inputs
if nargin==0
    [file_name,path_name] = uigetfile('*.json','Please provide with the JSON file containing the result of the quantify_dataset.m');
    json_file = fullfile(path_name,file_name);
elseif nargin==1
    json_file = varargin{1};
else
    error('Too many inputs');
end

[json_path,json_filename,~] = fileparts(json_file);
json_data = loadjson(json_file);
% all base region
obj_cell = json_data.objects{1};
reg_cell = json_data.regions{1};
%
obj_struct = [obj_cell{:}];
reg_struct = [reg_cell{:}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at the regions
reg_name_lst = {reg_struct(:).name};
[reg_name_lst_unq,idx_reg_name_lst_unq,~] = unique(reg_name_lst,'sorted');
%
reg_area_lst = [reg_struct(:).area];
reg_area_units_lst = {reg_struct(:).area_units};
%
reg_pixel_lst = [reg_struct(:).pixel];
%
reg_idx_lst = [reg_struct(:).idx];
reg_idx_lst_unq = reg_idx_lst(idx_reg_name_lst_unq);
%
reg_clr_lst = vertcat(reg_struct(:).rgb);
reg_clr_lst_unq = reg_clr_lst(idx_reg_name_lst_unq,:);
%
% Check units
if length(unique(reg_area_units_lst))==1
    reg_area_unit = unique(reg_area_units_lst);
    reg_area_unit = reg_area_unit{:};
else
    error('Discrepency in area units for regions. You should check that');
end
% Init
n_unq_reg = length(reg_idx_lst_unq);
reg_empty = struct('name','',...
    'idx',0,...
    'pxl',0,...
    'area',0,...
    'area_unit',reg_area_unit);
reg_info = repmat(reg_empty,n_unq_reg,1);
% Loop on regions
for iR = 1 : n_unq_reg
    reg_info(iR).name = reg_name_lst_unq{iR};
    reg_info(iR).idx  = reg_idx_lst_unq(iR);
    reg_info(iR).rgb  = reg_clr_lst_unq(iR,:);
    reg_info(iR).pxl  = sum(reg_pixel_lst(reg_idx_lst==reg_idx_lst_unq(iR)));
    reg_info(iR).area = sum(reg_area_lst(reg_idx_lst==reg_idx_lst_unq(iR)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at the objects 
obj_name_lst = {obj_struct(:).region_name};
[obj_name_lst_unq,idx_obj_name_lst_unq,~] = unique(obj_name_lst,'sorted');
%
obj_area_lst       = [obj_struct(:).object_area];
obj_area_units_lst = {obj_struct(:).object_area_units};
%
obj_pixel_lst = [obj_struct(:).object_pixel];
%
obj_idx_lst     = [obj_struct(:).region_idx];
obj_idx_lst_unq = obj_idx_lst(idx_obj_name_lst_unq);
%
obj_coord_lst = vertcat(obj_struct(:).object_centroid_atlas);
obj_coord_lst_unq = obj_coord_lst(idx_obj_name_lst_unq,:);
%
% obj_mean_r = [obj_struct(:).object_meanR];
% obj_mean_g = [obj_struct(:).object_meanG];
% obj_mean_b = [obj_struct(:).object_meanB];
%
% Check units
if length(unique(obj_area_units_lst))==1
    obj_area_unit = unique(obj_area_units_lst);
    obj_area_unit = obj_area_unit{:};
else
    error('Discrepency in area units for regions. You should check that');
end
% Init
n_unq_obj = length(obj_idx_lst_unq);
obj_empty = struct('name','',...
    'idx',0,...
    'cnt',0,...
    'pxl',0,...
    'area',0,...
    'area_unit',obj_area_unit);
obj_info = repmat(obj_empty,n_unq_obj,1);
%
for iR = 1 : n_unq_obj
    %
    obj_info(iR).name  = obj_name_lst_unq{iR};
    obj_info(iR).idx   = obj_idx_lst_unq(iR);
    obj_info(iR).cnt   = length(find(obj_idx_lst==obj_idx_lst_unq(iR)));
    obj_info(iR).pxl   = sum(obj_pixel_lst(obj_idx_lst==obj_idx_lst_unq(iR)));
    obj_info(iR).area  = sum(obj_area_lst(obj_idx_lst==obj_idx_lst_unq(iR)));
    obj_info(iR).coord = obj_coord_lst(obj_idx_lst==obj_idx_lst_unq(iR),:);
    %
%     obj_info(iR).mean_red   = mean(obj_mean_r(obj_idx_lst==obj_idx_lst_unq(iR)));
%     obj_info(iR).mean_green = mean(obj_mean_g(obj_idx_lst==obj_idx_lst_unq(iR)));
%     obj_info(iR).mean_blue  = mean(obj_mean_b(obj_idx_lst==obj_idx_lst_unq(iR)));
end
% combine
n_reg = length(reg_info);
stats_empty = struct('reg_name','',...
    'reg_idx',0,...
    'reg_pxl',0,...
    'reg_area',0,...
    'reg_area_unit',reg_area_unit,...
    'obj_cnt',0,...
    'obj_reg_cnt_ratio',0,...
    'obj_pxl',0,...
    'obj_coord',[],...
    'obj_area',0,...
    'obj_area_unit',obj_area_unit,...
    'obj_reg_area_ratio',0);
%
stats = repmat(stats_empty,n_reg,1);
%
for iR = 1 : n_reg
    idx_obj_in_reg = find([obj_info(:).idx]==reg_info(iR).idx);
    stats(iR).reg_name = reg_name_lst_unq{iR};
    stats(iR).reg_idx  = reg_info(iR).idx;
    stats(iR).reg_pxl  = reg_info(iR).pxl;
    stats(iR).reg_area = reg_info(iR).area;
    stats(iR).reg_rgb  = reg_info(iR).rgb; 
    stats(iR).reg_area_unit = reg_info(iR).area_unit;
    if ~(idx_obj_in_reg==0)
        %
        stats(iR).obj_cnt = obj_info(idx_obj_in_reg).cnt;
        stats(iR).obj_reg_cnt_ratio = stats(iR).obj_cnt./stats(iR).reg_area;
        %
        stats(iR).obj_pxl = obj_info(idx_obj_in_reg).pxl;
        %
        stats(iR).obj_coord = obj_info(idx_obj_in_reg).coord;
        %
        stats(iR).obj_area = obj_info(idx_obj_in_reg).area;
        stats(iR).obj_area_unit = obj_info(idx_obj_in_reg).area_unit;
        stats(iR).obj_reg_area_ratio = stats(iR).obj_area./stats(iR).reg_area;
    end
    %
end
% save data json
stats_json_filename = fullfile(json_path,sprintf('%s_stats.json',json_filename));
savejson('',stats,stats_json_filename);
% save data excel
stats_xls_filename = fullfile(json_path,sprintf('%s_stats.xlsx',json_filename));
writetable(struct2table(stats),stats_xls_filename);
return

