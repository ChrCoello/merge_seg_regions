function [hier_json_file,hier_xls_file] = combine_hierarchy(stats_json_file,varargin)
% TO DO
%
%

%%%Check JSONlab
checkJson();

%%%Load json
[~,json_fn,~]=fileparts(stats_json_file);
js = loadjson(stats_json_file);
js_data = js.stats{1};
js_data = [js_data{:}];

%%%%List regions
info = js.study_info;
if isfield(info,'hier_dir')
    list_reg_dir = info.hier_dir;
    list_reg = dir([list_reg_dir '\ABAHier2015_*']);
elseif nargin>1
    info = loadjson(varargin{1});
    list_reg_dir = info.hier_dir;
    list_reg = dir([list_reg_dir '\ABAHier2015_*']);
else
    error('combine_hierarchy:MissingField',...
        'The hier_dir field is missing from the study info information');
end
fprintf(1,'Read %d hierarchies\n',length(list_reg));
for iL = 1:length(list_reg)
    hier_fn = fullfile(list_reg_dir,list_reg(iL).name);
    idx = xlsread(hier_fn);
    iC = 0;
    idx_sub = [];
    for iI = 1 : length(idx)
        if ~isempty(find([js_data(:).reg_idx]==idx(iI),1))
            iC=iC+1;
            idx_sub(iC) = idx(iI);
            structure(iC)=js_data(([js_data(:).reg_idx]==idx(iI)));
        end
    end
    if ~exist('structure','var')
        warning('No regions from hierarchy file %s could be found in the file %s. Skipping...\n',...
            hier_fn,stats_json_file)
    else
       
        hier_nm = list_reg(iL).name(strfind(list_reg(iL).name,'_')+1:end-5);
        fprintf(1,' - Creating hierarchy %d / %d : %s\n',iL,length(list_reg),hier_nm);
        %
        hier_region(iL).reg_name = hier_nm;
        hier_region(iL).reg_idx_full    = sprintfc('%d',idx);
        hier_region(iL).reg_idx         = sprintfc('%d',idx_sub);
        hier_region(iL).reg_pxl         = sum([structure(:).reg_pxl]);
        hier_region(iL).reg_area        = sum([structure(:).reg_area]);
        hier_region(iL).reg_area_unit   = structure(end).reg_area_unit;
        hier_region(iL).obj_cnt         = sum([structure(:).obj_cnt]);
        hier_region(iL).obj_reg_cnt_ratio = hier_region(iL).obj_cnt ./ hier_region(iL).reg_area;
        hier_region(iL).obj_pxl         = sum([structure(:).obj_pxl]);
        hier_region(iL).obj_area        = sum([structure(:).obj_area]);
        hier_region(iL).obj_area_unit   = structure(end).obj_area_unit;
        hier_region(iL).obj_reg_area_ratio = hier_region(iL).obj_area ./ hier_region(iL).reg_area;
        %
        clear structure
    end
end

%%% Write results
%as json
hier_json_file = fullfile(info.output_dir,...
    sprintf('%s_objects_per_hierarchy.json',info.study_name));
savejson('',hier_region,hier_json_file);
%as excel
output_dir_xls = fullfile(info.output_dir,'excel');
if ~exist(output_dir_xls,'dir')
    mkdir(output_dir_xls);
end
plaques = struct2table(hier_region);
hier_xls_file = fullfile(output_dir_xls,sprintf('%s_objects_per_hierarchy.xlsx',json_fn));
writetable(plaques,hier_xls_file);

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