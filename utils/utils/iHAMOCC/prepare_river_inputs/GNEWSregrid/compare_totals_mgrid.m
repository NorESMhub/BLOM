% J. Schwinger 20.09.2024 
% Test that files created by the riv_nut_mkinpfile-script have the 
% same amount of nutrient loads


clear all


global use_octave

use_octave=true;

if use_octave
    pkg load netcdf
	import_netcdf
end

%file1='river_nutrients_GNEWS2000_tnx1v4_20170820.nc'
file1='river_nutrients_GNEWS2000c00_tnx1v4_20250220.nc'
%file2='river_nutrients_GNEWS2000c70_tnx1v4_20240921.nc'
%file2='river_nutrients_GNEWS2000c70_tnx0.5v1_20250220.nc'
file2='river_nutrients_GNEWS2000c00_tnx0.5v1_20250220.nc'


grid1='grid_tnx1v4_20170622.nc'
%grid2='grid_tnx1v4_20170622.nc'
grid2='grid_tnx0.5v1_20240702.nc'
vars = {'Qact','Qnat','DET','DIN','DIP','DSi','DIC','DOC','Fe'}
nvars= length(vars)

area1=ncread(grid1,'parea');
area2=ncread(grid2,'parea');
mask1=ncread(grid1,'pmask');
mask2=ncread(grid2,'pmask');

for iv=1:nvars

    var1=ncread(file1,vars{iv});
    var2=ncread(file2,vars{iv});

    disp([vars{iv} ': '] )
    tot1=sum(var1(:).*double(mask1(:)).*area1(:)) 
    tot2=sum(var2(:).*double(mask2(:)).*area2(:))
    rdiff=(tot1-tot2)/tot1
    
end
