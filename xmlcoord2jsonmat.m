function json_fn  = xmlcoord2jsonmat(xml_filename)
% TO DO COMMENT
% Change & to ,
fid = fopen(xml_filename,'rt') ;
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;
% replace string S1 with string S2
Y = strrep(X, '&', ',') ;
[pp,tt,~]=fileparts(xml_filename);
tmp_xml = fullfile(pp,[tt '_withcomma.xml']);
fid2 = fopen(tmp_xml,'wt') ;
fwrite(fid2,Y) ;
fclose (fid2) ;

% xml2struct : yay we love it
S=xml2struct(tmp_xml);
delete(tmp_xml);

%
json_struct.(S.Name)=struct(S.Attributes(1).Name,S.Attributes(1).Value,...
    S.Attributes(2).Name,S.Attributes(2).Value,...
    S.Attributes(3).Name,S.Attributes(3).Value);

% S children
C=S.Children;
data=C(~cellfun('isempty',strfind({C(:).Name},'slice')));

for iR = 1:length(data)
    json_struct.slice(iR).nr  = data(iR).Attributes(~cellfun('isempty',strfind({data(iR).Attributes(:).Name},'nr'))).Value;
    json_struct.slice(iR).filename = data(iR).Attributes(~cellfun('isempty',strfind({data(iR).Attributes(:).Name},'filename'))).Value;
    json_struct.slice(iR).height = data(iR).Attributes(~cellfun('isempty',strfind({data(iR).Attributes(:).Name},'height'))).Value;
    json_struct.slice(iR).width  = data(iR).Attributes(~cellfun('isempty',strfind({data(iR).Attributes(:).Name},'width'))).Value;
    llt=data(iR).Attributes(~cellfun('isempty',strfind({data(iR).Attributes(:).Name},'anchoring'))).Value;
    equal_sym=strfind(llt,'=');
    comma_sym=[0 strfind(llt,',') length(llt)];
    for iT=1:length(equal_sym)
        json_struct.slice(iR).(llt(comma_sym(iT)+1:equal_sym(iT)-1))=str2double(llt((equal_sym(iT)+1):(comma_sym(iT+1)-1)));
    end
    json_struct.slice(iR).transf_mat = [json_struct.slice(iR).ux json_struct.slice(iR).ox json_struct.slice(iR).vx;...
        json_struct.slice(iR).uy json_struct.slice(iR).oy json_struct.slice(iR).vy;...
        json_struct.slice(iR).uz json_struct.slice(iR).oz json_struct.slice(iR).vz];
end
json_fn = fullfile(pp,[tt '.json']);
savejson('',json_struct,json_fn);
return