% sgao 07.04.16
% write riverine nutrient input into NetCDF
function riv_nut_ncwrite_model(coordi_f,coordj_f,scenario,data_ocegrid,resol_tag,vstr)

% define dimensions 
IDM=coordi_f;
JDM=coordj_f;
%
%
FILENAME=['river_nutrients_GNEWS2000' scenario '_' char(resol_tag) '_' vstr '.nc'];
%
VARNAME{1}='Qact';
VARLONGNAME{1}='Annual actual (anthropogenically disturbed) discharge at the basin mouth';
VARUNITS{1}='kg m-2 s-1';
%
VARNAME{2}='Qnat';
VARLONGNAME{2}='Annual natural discharge at the basin mouth';
VARUNITS{2}='kg m-2 s-1';
%
VARNAME{3}='DIN';
VARLONGNAME{3}='dissolved inorganic nitrogen';
VARUNITS{3}='kmol m-2 yr-1';
%
VARNAME{4}='DIP';
VARLONGNAME{4}='dissolved inorganic phosphorous';
VARUNITS{4}='kmol m-2 yr-1';
%
VARNAME{5}='DSi';
VARLONGNAME{5}='dissolved silicon';
VARUNITS{5}='kmol m-2 yr-1';
%
VARNAME{6}='DIC';
VARLONGNAME{6}='dissolved inorganic carbon';
VARUNITS{6}='kmol m-2 yr-1';
%
VARNAME{7}='DET';
VARLONGNAME{7}='detritus';
VARUNITS{7}='kmol P m-2 yr-1';
%
VARNAME{8}='DOC';
VARLONGNAME{8}='dissolved organic carbon';
VARUNITS{8}='kmol m-2 yr-1';
%
VARNAME{9}='Fe';
VARLONGNAME{9}='dissolved iron';
VARUNITS{9}='kmol m-2 yr-1';

num_var = length(VARNAME);

% delete output file if exists
delete(FILENAME);

% write data to file 
for n = 1:num_var
    data = reshape(data_ocegrid(:,n),IDM,[]);
    nccreate(FILENAME,VARNAME{n},'Format','classic','Dimensions',{'x' IDM 'y' JDM});
    ncwrite(FILENAME,VARNAME{n},data);
    ncwriteatt(FILENAME,VARNAME{n},'long_name',VARLONGNAME{n});
    ncwriteatt(FILENAME,VARNAME{n},'units',VARUNITS{n});
end

end
