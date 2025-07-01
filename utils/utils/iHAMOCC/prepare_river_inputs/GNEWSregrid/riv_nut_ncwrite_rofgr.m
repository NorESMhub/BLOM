% J. Schwinger 20.09.2024: write riverine nutrient input on runoff-grid into NetCDF
function riv_nut_ncwrite_rofgr(idm,jdm,lon,lat,data_rofgr,scenario,vstr)

global use_octave
if ~isempty(use_octave) && use_octave
    import_netcdf
end

%
filename=['river_nutrients_GNEWS2000' scenario '_r05_' vstr '.nc'];
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
delete(filename);

% Define netCDF file
%cmode  = netcdf.getConstant('NETCDF4');
cmode  = netcdf.getConstant('NC_64BIT_DATA');
cmode2 = netcdf.getConstant('CLASSIC_MODEL');
cmode3 = netcdf.getConstant('CLOBBER');
cmode  = bitor(bitor(cmode,cmode2),cmode3);
ncid   = netcdf.create(filename,cmode);

did_lon  = netcdf.defDim(ncid,'lon',idm);
did_lat  = netcdf.defDim(ncid,'lat',jdm);
vid_lon  = netcdf.defVar(ncid,'lon', 'NC_DOUBLE', did_lon);
vid_lat  = netcdf.defVar(ncid,'lat', 'NC_DOUBLE', did_lat);

netcdf.putAtt(ncid,vid_lon,'long_name','longitude');
netcdf.putAtt(ncid,vid_lat,'long_name','latitude');
netcdf.putAtt(ncid,vid_lon,'units','degrees_east');
netcdf.putAtt(ncid,vid_lat,'units','degrees_north');

vid = zeros(num_var,1);
for n=1:num_var
    vid(n)=netcdf.defVar(ncid,VARNAME{n},'NC_DOUBLE',[did_lon did_lat]);
    netcdf.putAtt(ncid,vid(n),'long_name',VARLONGNAME{n});
    netcdf.putAtt(ncid,vid(n),'units',VARUNITS{n});
end

% Global attributes
vid_g = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,vid_g,'Source','GNEWS2 Dataset version 1.0. 2014-11-12');
netcdf.putAtt(ncid,vid_g,'Creation_date',datestr(date,29));
netcdf.putAtt(ncid,vid_g,'Created_by','Jörg Schwinger, jorg.schwinger0@gmail.com');
netcdf.putAtt(ncid,vid_g,'Conventions','CF-1.8');

netcdf.endDef(ncid)

for n = 1:num_var
    netcdf.putVar(ncid,vid_lon,lon)
    netcdf.putVar(ncid,vid_lat,lat)
    data = reshape(data_rofgr(:,n),idm,[]);
    netcdf.putVar(ncid,vid(n),data);
end


netcdf.close(ncid)



end
