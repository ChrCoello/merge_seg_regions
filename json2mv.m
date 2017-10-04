function json2mv(stats_json_file, clr_fix)
% TO DO
%
%
if nargin<2
    clr_fix = false;
end
%%%Load json
[json_path,json_fn,~]=fileparts(stats_json_file);
js = loadjson(stats_json_file);
js_data = js{1};
js_data = [js_data{:}];

%%% output folder
output_dir = fullfile(json_path,'mv');
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

fprintf(1,'Creating %d Meshview compatibles text files : one file per region\n',length(js_data));
fid0 = fopen(fullfile(output_dir,['00_' json_fn '.txt']),'w+');
for iL = 1:length(js_data)
    cur_reg = js_data(iL);
    if ~isempty(cur_reg.obj_coord)
        % remove weird chararcters
        cur_reg.reg_name(regexp(cur_reg.reg_name,'\W'))='_';
        % Create the txt file
        fid1 = fopen(fullfile(output_dir,[cur_reg.reg_name '.txt']),'w+');
        % RGBA 1 0 0 1 # RGBA
        cur_clr = cur_reg.reg_clr/255;
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
web('http://www.nesys.uio.no/MeshGen/MeshView.html?bitlas=ABAv3.bitlas', '-browser')
%
