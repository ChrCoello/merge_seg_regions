function [output_json,study_info_json] = quantify_dataset(varargin)
%
%
% TO DO
% Input : study_info_json
%
% Parse inputs
if nargin==0
    [file_name,path_name] = uigetfile('*.json','Please provide with the JSON file containing the study information (i.e study_info.json)');
    study_info_json = fullfile(path_name,file_name);
elseif nargin==1
    study_info_json = varargin{1};
else
    error('Too many inputs');
end

% % Check if RemoveSheet123 is there and add it if possible. If not, then
% % don't use the function in the code
% if ~(exist('RemoveSheet123','file')==2)
%     if exist('RemoveSheet123','dir')
%         addpath(fullfile(pwd,'RemoveSheet123'));
%     end
% end
% Check that the JSONlab is in the path and do something about it
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
    fprintf(1,'\nJSONlab toolbox is detected.');
end
study_info = loadjson(study_info_json);
% Add parsing inputs (one day)
study_name     = study_info.study_name;
atlas_dir      = study_info.atlas_dir;
atlas_lbl_file = study_info.atlas_lbl_file;
atlas_xml_file = study_info.atlas_xml_file;
seg_dir        = study_info.seg_dir;
obj_lbl        = study_info.obj_lbl;
slice_dir      = study_info.slice_dir;
output_dir     = study_info.output_dir;
original_dir   = study_info.original_dir;
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
GlobalStats.objects = [];
GlobalStats.regions = [];
%
%%
seg_dir_ctn = dir([seg_dir '*.png']);
%
n_slice = length(seg_dir_ctn);
GlobalStats.n_slices = n_slice;
GlobalStats.n_objects = NaN;
GlobalStats.n_regions = NaN;
fprintf(1,'\nTotal number of section detected: %d',n_slice);
%%
atlas_dir_ctn = dir(atlas_dir);
atlas_dir_ctn = atlas_dir_ctn(~[atlas_dir_ctn(:).isdir]);
%%
slice_dir_ctn = dir([slice_dir '*.png']);
%
output_xls_obj_ind = fullfile(output_dir,[study_name '_obj_ind.xlsx']);
output_xls_reg_ind = fullfile(output_dir,[study_name '_reg_ind.xlsx']);

%% Coordinates
atlas_json_fn    = xmlcoord2jsonmat(atlas_xml_file);
atlas_coord_json = loadjson(atlas_json_fn);
sections_coord = [atlas_coord_json.slice{:}];

%% Original dir
ori_dir_ctn = dir([slice_dir '*.tif']);

%% Loop on all the images
n_objects = 0;
n_regions = 0;
for iS = n_slice:n_slice
    %
    fprintf(1,'\n -- Analyzing slice #%d / %d',iS,n_slice);
    %
    seg_name_orig = seg_dir_ctn(iS).name;
    % Remove Object Prediction from the file name if its finds it
    if strfind(seg_name_orig,'Object Prediction')
        seg_name = seg_name_orig(1:strfind(seg_name_orig,'Object Prediction')-2);
    end
    seg_name_file = fullfile(seg_dir,seg_name_orig);
    % Atlas
    atlas_slices = atlas_dir_ctn(~cellfun('isempty',strfind({atlas_dir_ctn(:).name},seg_name)));
    %
    atlas_name = atlas_slices(~cellfun('isempty',strfind({atlas_slices(:).name},'.bin'))).name;
    atlas_name_file = fullfile(atlas_dir,atlas_name);
    % Section
    slice_name = slice_dir_ctn(~cellfun('isempty',strfind({slice_dir_ctn(:).name},seg_name))).name;
    slice_name_file = fullfile(slice_dir,slice_name);
    %
    slice_txt       = fullfile(slice_dir,[slice_name(1:end-4) '.txt']);
    if ~exist(slice_txt,'file')
        % orignal tif
        try
        ori_name = ori_dir_ctn(~cellfun('isempty',strfind({ori_dir_ctn(:).name},seg_name))).name;
        ori_metadata = imfinfo(fullfile(original_dir,ori_name));
        metadata.width  = ori_metadata.Width;
        metadata.height = ori_metadata.Height;
        xresolution = unit_convert(ori_metadata.XResolution,ori_metadata.ResolutionUnit,'um');
        yresolution = unit_convert(ori_metadata.YResolution,ori_metadata.ResolutionUnit,'um');
        %         metadata.resolution_unit = 'pixel/um';
        metadata.x_pixel_size = 1/xresolution;
        metadata.y_pixel_size = 1/yresolution;
        metadata.pixel_size_unit = 'um/pixel';
        catch
            error('quantify_dataset:fetching_metadata',...
                'Unable to fetch the appropriate metadata from original tiff file');
        end
    else
        txt_section = load_txt(slice_txt);
        if txt_section{2,3}~=txt_section{3,3}
            error('non isotropic pixel resolution');
        else
            metadata.width  = str2double(txt_section{4,4});
            metadata.height = str2double(txt_section{5,4});
            metadata.x_pixel_size = txt_section{2,3};
            metadata.y_pixel_size = txt_section{3,3};
            metadata.pixel_size_unit = 'um';
        end
    end
    % Coord
    metadata.pixel_to_atlas_mat = sections_coord(~cellfun('isempty',...
        strfind({sections_coord(:).filename},seg_name))).transf_mat;

    
    [obj_stats,reg_stats] = quantify_single_section(atlas_name_file,seg_name_file,...
        slice_name_file,atlas_lbl_file,output_dir,obj_lbl,metadata);
    % Concatenate the individual objects with the ones from preivous
    % sections
    n_objects =+ length(obj_stats);
    n_regions =+ length(reg_stats);
    %
    GlobalStats.objects = vertcat(GlobalStats.objects,obj_stats);
    GlobalStats.regions = vertcat(GlobalStats.regions,reg_stats);
    %as excel
    objects = struct2table(obj_stats);
    regions = struct2table(reg_stats);
    %
    warning('off','MATLAB:xlswrite:AddSheet');
    writetable(objects,output_xls_obj_ind,'Sheet',slice_name(end-7:end-4));
    writetable(regions,output_xls_reg_ind,'Sheet',slice_name(end-7:end-4));
    %
    fprintf('\n -- Analyzing slice #%d / %d -- done',iS,n_slice);
end
% if exist('RemoveSheet123','file')==2
%     RemoveSheet123(output_xls_obj_ind);
%     RemoveSheet123(output_xls_obj_ind);
% end
%% Write output
GlobalStats.n_objects = n_objects;
GlobalStats.n_regions = n_regions;
%as json
output_json = fullfile(output_dir,[study_name '_obj_reg_data.json']);
savejson('',GlobalStats,output_json);
%as excel
objects = struct2table(GlobalStats.objects);
regions = struct2table(GlobalStats.regions);
%
output_xls_obj = fullfile(output_dir,[study_name '_obj.xlsx']);
output_xls_reg = fullfile(output_dir,[study_name '_reg.xlsx']);
%
writetable(objects,output_xls_obj);
writetable(regions,output_xls_reg);
%
t_p = toc;
fprintf(1,'\n\nAnalysis completed in %.0f seconds\n',t_p);
return

function struct_txt = load_txt(filename, startRow, endRow)
%IMPORTFILE Import text from the automatically generated text file that
%contains the resolution of the original file
%
% Example:
%   tg2576m2871D1s054 = importfile('tg2576_m287_1D1_s054.txt', 2, 7);
%
%    See also TEXTSCAN.
% Auto-generated by MATLAB on 2017/08/15 11:39:55

% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 2;
    endRow = 7;
end

% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Converts text in the input cell array to numbers. Replaced non-numeric
% text with NaN.
rawData = dataArray{3};
for row=1:size(rawData, 1)
    % Create a regular expression to detect and remove non-numeric prefixes and
    % suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData{row}, regexstr, 'names');
        numbers = result.numbers;

        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if any(numbers==',')
            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(numbers, thousandsRegExp, 'once'))
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % Convert numeric text to numbers.
        if ~invalidThousandsSeparator
            numbers = textscan(strrep(numbers, ',', ''), '%f');
            numericData(row, 3) = numbers{1};
            raw{row, 3} = numbers{1};
        end
    catch me
    end
end

%
struct_txt = raw;

return

function [outData] = unit_convert(inData,inUnits,outUnits)
% 
conversionFactors = {
    'length',...
    {'m','cm','Centimeter','mm','um','nm','km','in','ft','yd','mile','AU'},...
    [1 1e-2 1e-2 1e-3 1e-6 1e-9 1e3 2.54e-2 12*2.54e-2 36*2.54e-2 1760*36*2.54e-2 1.49598e11]};
    
inUnits  = strtrim(inUnits);
outUnits = strtrim(outUnits);
%
if ~ismember(inUnits,conversionFactors{2}) || ~ismember(outUnits,conversionFactors{2})
    error('UnitConvert:IncompatibleUnits','%s%s%s%s%s','Specified input (',inUnits,') and output (',outUnits,') units are incompatible.')
end
%
factorInToStandard  = conversionFactors{3}(strcmp(inUnits,conversionFactors{2}));
factorStandardToOut = 1/conversionFactors{3}(strcmp(outUnits,conversionFactors{2}));
factorInToOut = factorInToStandard*factorStandardToOut;
%
outData       = factorInToOut*inData;
