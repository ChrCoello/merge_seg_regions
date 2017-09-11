function [output_json,output_xls_obj,output_xls_reg] = quantify_dataset(study_info_json)
%
%
% TO DO
%
% Check that the JSONlab is in the path and do something about it
if ~(exist('loadjson','file')==2)
    try 
        %install local
        if exist('jsonlab-1.5','dir')
            addpath(genpath(fullpath(pwd,'jsonlab-1.5')));
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
    fprintf(1,'\nJSONlab toolbox is detected.');
end
study_info = loadjson(study_info_json);
% Add parsing inputs (one day)
study_name     = study_info.study_name;
atlas_dir      = study_info.atlas_dir;
atlas_lbl_file = study_info.atlas_lbl_file;
seg_dir        = study_info.seg_dir;
obj_lbl        = study_info.obj_lbl;
slice_dir      = study_info.slice_dir;
output_dir     = study_info.output_dir;
%

%% Get info on system
tic;
sysinfo.username = getenv('username');
% Computer type
sysinfo.platformUsed  = computer;
%OS used
sysinfo.osType        = system_dependent('getos');
% MATLAB
sysinfo.matlabVersion = version;
%%
GlobalStats.date_analysis = datestr(now);
GlobalStats.system_info = sysinfo;
GlobalStats.plaques = [];
GlobalStats.regions = [];
%
%%
seg_dir_ctn = dir([seg_dir '*.png']);
%
n_slice = length(seg_dir_ctn);
fprintf(1,'\nTotal number of section detected: %d',n_slice);
%%
atlas_dir_ctn = dir(atlas_dir);
atlas_dir_ctn = atlas_dir_ctn(~[atlas_dir_ctn(:).isdir]);
%%
slice_dir_ctn = dir([slice_dir '*.png']);
%
%% Loop on all the images
for iS = n_slice:n_slice %1 : n_slice
    %
    fprintf(1,'\n -- Analyzing slice #%d / %d',iS,n_slice);
    %
    seg_name_orig = seg_dir_ctn(iS).name;
    % Remove Object Prediction from the file name if its finds it
    if strfind(seg_name_orig,'Object Prediction')
        seg_name = seg_name_orig(1:strfind(seg_name_orig,'Object Prediction')-2);
    end
    seg_name_file = fullfile(seg_dir,seg_name_orig);
    %
    atlas_slices = atlas_dir_ctn(~cellfun('isempty',strfind({atlas_dir_ctn(:).name},seg_name)));
    %
    atlas_name = atlas_slices(~cellfun('isempty',strfind({atlas_slices(:).name},'.bin'))).name;
    atlas_name_file = fullfile(atlas_dir,atlas_name);
    %
    slice_name = slice_dir_ctn(~cellfun('isempty',strfind({slice_dir_ctn(:).name},seg_name))).name;
    slice_name_file = fullfile(slice_dir,slice_name);
    slice_json      = fullfile(slice_dir,[slice_name(1:end-4) '.json']);
    slice_txt       = fullfile(slice_dir,[slice_name(1:end-4) '.txt']);
    
    [obj_stats,reg_stats] = quantify_single_section(atlas_name_file,seg_name_file,...
        slice_name_file,atlas_lbl_file,slice_json,slice_txt,output_dir,obj_lbl);
    % Concatenate the individual objects with the ones from preivous
    % sections
    GlobalStats.plaques = vertcat(GlobalStats.plaques,obj_stats);
    GlobalStats.regions = vertcat(GlobalStats.regions,reg_stats);
    fprintf('\n -- Analyzing slice #%d / %d -- done \n',iS,n_slice);
end
%% Write output
%as json
output_json = fullfile(output_dir,[study_name '_obj_reg_data.json']);
savejson('',GlobalStats,output_json);
%as excel
plaques = struct2table(GlobalStats.plaques);
regions = struct2table(GlobalStats.regions);
%
output_xls_obj = fullfile(output_dir,[study_name '_obj.xlsx']);
output_xls_reg = fullfile(output_dir,[study_name '_reg.xlsx']);
%
writetable(plaques,output_xls_obj);
writetable(regions,output_xls_reg);
%
t_p = toc;
fprintf(1,'Analysis completed in %.0f seconds\n',t_p);
return