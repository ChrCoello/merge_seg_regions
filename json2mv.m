function json2mv(stats_json_file,clr_fix)
% TO DO
%
%
if nargin<2
    clr_fix = false;
end
%%%check Jsonlab
checkJson();

%%%load json
[json_path,json_fn,~]=fileparts(stats_json_file);
js = loadjson(stats_json_file);
js_data = js.stats{1};
js_data = [js_data{:}];

%%% output folder
output_dir = fullfile(json_path,'mv');
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

fprintf(1,'Creating %d Meshview compatibles text files : one file per region\n',length(js_data));
fid0 = fopen(fullfile(js.study_info.output_dir,[js.study_info.study_name '_meshview.txt']),'w+');
output_dir_txt = fullfile(js.study_info.output_dir,'meshview_ind_reg');
if ~exist(output_dir_txt,'dir')
    mkdir(output_dir_txt);
end
for iL = 1:length(js_data)
    cur_reg = js_data(iL);
    if ~isempty(cur_reg.obj_coord)
        % remove weird chararcters
        cur_reg.reg_name(regexp(cur_reg.reg_name,'\W'))='_';
        % Create the txt file
        fid1 = fopen(fullfile(output_dir_txt,[cur_reg.reg_name '.txt']),'w+');
        % RGBA 1 0 0 1 # RGBA
        cur_clr = cur_reg.reg_rgb/255;
        if clr_fix
            cur_clr = [1 0 0];
        end
        fprintf(fid1,'RGBA %0.3f %0.3f %0.3f 1 # RGBA\n',cur_clr);
        fprintf(fid1,'%d,%d,%d\n',round(cur_reg.obj_coord'));
        fclose(fid1);
        %
        fprintf(fid0,'RGBA %0.3f %0.3f %0.3f 1 # RGBA\n',cur_clr);
        fprintf(fid0,'%d,%d,%d\n',round(cur_reg.obj_coord'));
        %
    end
end
fclose(fid0);
% having fun
% web('http://www.nesys.uio.no/MeshGen/MeshView.html?bitlas=ABAv3.bitlas', '-browser')
%

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
