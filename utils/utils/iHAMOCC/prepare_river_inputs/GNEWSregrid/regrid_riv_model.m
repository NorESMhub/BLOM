% Shuang Gao 16.08.17, Jörg Schwinger 20.09.2024: interpolate nutrient input first to runoff grid 
% and then onto ocean model grid (spread out as runoff)
function [idm,jdm,gndata_oceangrid] = regrid_riv_model(gnlon,gnlat,gndata,scenarios,resol_tag)


% read NorESM coordinates and land/ocean mask 
GRIDFILE=['grid_' char(resol_tag) '.nc']
lon=ncread(GRIDFILE,'plon');
lat=ncread(GRIDFILE,'plat');
mask=ncread(GRIDFILE,'pmask');
nreg=ncreadatt(GRIDFILE,'/','nreg'); % 1=dipolar 2=tripolar

% clip pole grid cell row (use only for tnx-grids!)
if nreg==2
  lon=lon(:,1:end-1);
  lat=lat(:,1:end-1);
  mask=mask(:,1:end-1);
end

% set i and j dimensions of model grid
[idm,jdm]=size(mask);

% read map file
switch char(resol_tag)
	case 'gx1v6'
	resol_postfix={'_e1000r300_090226'};
	case 'tnx2v1'
	resol_postfix={'_e1000r300_130206'};
	case 'tnx1v2'
	resol_postfix={'_e1000r300_140828'};
	case 'tnx1v4'
	resol_postfix={'_e1000r300_170609'};
	case 'tnx0.5v1'
	resol_postfix={'_e300r100_20240702'};
	case 'tnx0.25v4'
	resol_postfix={'_e300r100_170629'};
end 
MAPFILE=['map_r05_to_' char(resol_tag)  char(resol_postfix) '.nc']; 
xc_a=ncread(MAPFILE,'xc_a');
yc_a=ncread(MAPFILE,'yc_a');
area_a=ncread(MAPFILE,'area_a')*6.37122e6^2; % convert units from rad^2 to m^2 
mask_a=ncread(MAPFILE,'mask_a'); 


% distribute nutrient flux on 0.5degX0.5deg runoff grid (find nearest grid cell)
nnut=size(gndata,2); 
gndata_runoffgrid=zeros(length(area_a),nnut);
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
  gndata_runoffgrid(ij,:)=gndata_runoffgrid(ij,:)+gndata(n,:)/area_a(ij);
end

% plot first nutrient on runoff grid 
% figure(1)
% clf
% rfield=reshape(gndata_runoffgrid(:,1),720,360);
% rlon=reshape(xc_a,720,360);
% rlat=reshape(yc_a,720,360);
% imagesc(rlon(:,1),rlat(1,:),rfield');
% colorbar
% caxis([0 1]*1e-7);
% axis xy 

% interpolate from runoff grid to ocean grid 
area_b=ncread(MAPFILE,'area_b')*6.37122e6^2; % area of ocean grid cells
S=ncread(MAPFILE,'S'); % interpolation weights 
col=ncread(MAPFILE,'col'); % interpolation index
row=ncread(MAPFILE,'row'); % interpolation index
gndata_oceangrid=zeros(idm*jdm,nnut);
for n=1:length(S) 
  if mod(n,100000)==0 
    disp(int2str(n))
  end 
  gndata_oceangrid(row(n),:)=gndata_oceangrid(row(n),:)+S(n)*gndata_runoffgrid(col(n),:);
end 

%plot first nutrient on ocean grid 
% figure(2)
% clf
% ofield=reshape(gndata_oceangrid(:,1),idm,jdm);
% imagesc(ofield');
% colorbar
% caxis([0 1]*1e-7);
% axis xy

% compute global totals
tot_gn=sum(gndata)
tot_runoffgrid=area_a'*gndata_runoffgrid
tot_oceangrid=area_b'*gndata_oceangrid

% add extra line of cells for tripolar grids 
if nreg==2
  gndata_oceangrid=reshape(gndata_oceangrid,idm,jdm,[]);
  gndata_oceangrid=cat(2,gndata_oceangrid,flipdim(gndata_oceangrid(:,end,:),1));
  jdm=jdm+1;
  gndata_oceangrid=reshape(gndata_oceangrid,idm*jdm,[]);
end 


end
