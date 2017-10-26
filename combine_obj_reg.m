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

% Check JSON
checkJson();

% Load object list
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
fprintf(1,'Checking units\n');
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
fprintf(1,'Looping on the regions\n');
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
fprintf(1,'Looping on the objects\n');
unq_obj = 0;
for iR = 1 : n_unq_obj
    %
    obj_info(iR).name  = obj_name_lst_unq{iR};
    obj_info(iR).idx   = obj_idx_lst_unq(iR);
    obj_info(iR).cnt   = length(find(obj_idx_lst==obj_idx_lst_unq(iR)));
    obj_info(iR).pxl   = sum(obj_pixel_lst(obj_idx_lst==obj_idx_lst_unq(iR)));
    obj_info(iR).area  = sum(obj_area_lst(obj_idx_lst==obj_idx_lst_unq(iR)));
    obj_info(iR).coord = obj_coord_lst(obj_idx_lst==obj_idx_lst_unq(iR),:);
    unq_obj = unq_obj + length(obj_info(iR).coord);
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
region_struct.input_fn = json_file;
region_struct.date_analysis = datestr(now);
region_struct.stats = repmat(stats_empty,n_reg,1);
region_struct.study_info = json_data.study_info;
%
fprintf(1,'Combining %d objects in %d regions\n',unq_obj,n_reg);
for iR = 1 : n_reg
    idx_obj_in_reg = find([obj_info(:).idx]==reg_info(iR).idx);
    region_struct.stats(iR).reg_name = reg_name_lst_unq{iR};
    region_struct.stats(iR).reg_idx  = reg_info(iR).idx;
    region_struct.stats(iR).reg_pxl  = reg_info(iR).pxl;
    region_struct.stats(iR).reg_area = reg_info(iR).area;
    region_struct.stats(iR).reg_rgb  = reg_info(iR).rgb; 
    region_struct.stats(iR).reg_area_unit = reg_info(iR).area_unit;
    if ~(idx_obj_in_reg==0)
        %
        region_struct.stats(iR).obj_cnt = obj_info(idx_obj_in_reg).cnt;
        region_struct.stats(iR).obj_reg_cnt_ratio = region_struct.stats(iR).obj_cnt./region_struct.stats(iR).reg_area;
        %
        region_struct.stats(iR).obj_pxl = obj_info(idx_obj_in_reg).pxl;
        %
        region_struct.stats(iR).obj_coord = obj_info(idx_obj_in_reg).coord;
        %
        region_struct.stats(iR).obj_area = obj_info(idx_obj_in_reg).area;
        region_struct.stats(iR).obj_area_unit = obj_info(idx_obj_in_reg).area_unit;
        region_struct.stats(iR).obj_reg_area_ratio = region_struct.stats(iR).obj_area./region_struct.stats(iR).reg_area;
    end
    %
end
% save data json
stats_json_filename = fullfile(json_data.study_info.output_dir,...
    sprintf('%s_objects_per_region.json',json_data.study_info.study_name));
savejson('',region_struct,stats_json_filename);
% save data excel
output_dir_xls = fullfile(json_data.study_info.output_dir,'excel');
if ~exist(output_dir_xls,'dir')
    mkdir(output_dir_xls);
end
stats_xls_filename = fullfile(output_dir_xls,...
    sprintf('%s_objects_per_region.xlsx',json_data.study_info.study_name));
writetable(struct2table(region_struct.stats),stats_xls_filename);
return

function checkJson()
if ~(exist('loadjson','file')==2)
    try
        %install local
        if exist('jsonlab-1.5','dir')
            addpath(genpath(fullfile(pwd,'jsonlab-1.5')));
        else  %install on the nesys server
            if exist('Z:\NESYS_Tools\Matlab\jsonlab-1.5','dir')
                addpath(genpath('Z:\NESYS_Tools\Matlab\jsonlab-1.5'));
            elseif exist('Y:\NESYS_Tools\Matlab\jsonlab-1.5','dir')
                addpath(genpath('Y:\NESYS_Tools\Matlab\jsonlab-1.5'));
            end
        end
        if ~(exist('loadjson','file')==2)
            error('quantify_dataset:add_jsonlab',['Program tried to add the ',...
                'JSONlab package without success. Please follow instruction',...
                'on the README document to install the package']);
        end
    catch
        error('quantify_dataset:add_jsonlab',['Program tried to add the ',...
            'JSONlab package without success. Please follow instruction',...
            'on the README document to install the package']);
    end
else
    fprintf(1,'\nJSONlab toolbox is detected.\n');
end
return

