% Jörg Schwinger 20.09.2024:  interpolate nutrient input first on runoff grid
function [idm,jdm,lon,lat,gndata_rofgr] = regrid_riv_rofgr(gnlon,gnlat,gndata,scenarios,resol_tag)

idm=720;
jdm=360;

% read map file
switch char(resol_tag)
	case 'gx1v6'
	resol_postfix={'_e1000r300_090226'};
	case 'tnx1v2'
	resol_postfix={'_e1000r300_140828'};
	case 'tnx1v3'
	resol_postfix={'_e1000r300_170621'};
	case 'tnx1v4'
	resol_postfix={'_e1000r300_170609'};
	case 'tnx1v1'
	resol_postfix={'_e1000r300_120120'};
	case 'tnx0.25v1'
	resol_postfix={'_e300r100_130930'};
	case 'tnx0.25v3'
	resol_postfix={'_e300r100_170629'};
	case 'tnx0.25v4'
	resol_postfix={'_e300r100_170629'};
	case 'tnx2v1'
	resol_postfix={'_e1000r300_130206'};
end 
MAPFILE=['map_r05_to_' char(resol_tag)  char(resol_postfix) '.nc']; 
xc_a=ncread(MAPFILE,'xc_a');
yc_a=ncread(MAPFILE,'yc_a');
area_a=ncread(MAPFILE,'area_a')*6.37122e6^2; % convert units from rad^2 to m^2 
mask_a=ncread(MAPFILE,'mask_a'); 

% extract 1d longitude and latitude
lon = xc_a(1:idm);
lat = yc_a(1:idm:idm*jdm);


% distribute nutrient flux on 0.5 x 0.5 degree runoff grid (find nearest grid cell)
nnut=size(gndata,2); 
gndata_rofgr=zeros(length(area_a),nnut);

% mask_a seems to be unused
%oceanpoints=find(mask_a==1);
oceanpoints=find(mask_a==1|mask_a==0);

for n=1:length(gnlon) 

  if mod(n,1000)==0
    disp(int2str(n))
  end
  
  % compute distance between river mouth and all runoff-grid cells 
  dist=earthdist(gnlon(n),gnlat(n),xc_a(oceanpoints),yc_a(oceanpoints)); 
   
  % compute minimum distance and corresponding indices 
  [mindist,ij]=min(dist);
  ij=oceanpoints(ij);  

  % add nutrient flux to grid cell and divide by grid cell area (flux must be per m2)
  gndata_rofgr(ij,:)=gndata_rofgr(ij,:)+gndata(n,:)/area_a(ij);

end

% compute global totals
tot_gn=sum(gndata)
tot_runoffgrid=area_a'*gndata_rofgr


end
