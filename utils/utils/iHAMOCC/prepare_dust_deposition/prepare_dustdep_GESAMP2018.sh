#! /usr/bin/bash
#---------------------------------------------------------------------------------------------------
# This script generates dust/iron deposition input for BLOM/HAMOCC using data from the 
# GESAMP atmospheric iron deposition model intercomparison study (Myriokefalitakis et al. 2018, 
# doi:10.5194/bg-15-6659-2018). The data files have been downloaded from 
#
#      https://ecpl.chemistry.uoc.gr/GESAMP/ 
#
# This script generates files in two formats
#
#  1) Data on the original 1x1 degree grid to be used with the dust-stream functionally
#     (online time and spatial interpolation)
#
#  2) Data pre-interpolated to the model grid by using bi-linear interpolation.
#
# In both cases units are converted to kg Fe/m2/s and a new time axis in units of 
# "days since 1850-01-01" is inserted.
#----------------------------------------------------------------------------------------------------

#gridspec=tnx2v1
#gridspec=tnx1v4
#gridspec=tnx0.5v1
#gridspec=tnx0.25v4
gridspec=tnx0.125v4

datetag=20250516

# The script expects a link to the grid file (without date-tag in the link name) in the current directory
gridfile=grid_${gridspec}.nc
tinfile=../GESAMP_dust/gesamp.ensemble.tfe_dep.monthly.present_r360x180.nc
linfile=../GESAMP_dust/gesamp.ensemble.lfe_dep.monthly.present_r360x180.nc
outfile_ogrid=dustdep_GESAMP2018_r360x180_${datetag}.nc
outfile_mgrid=dustdep_GESAMP2018_${gridspec}_${datetag}.nc


rm $outfile_ogrid tmp1.nc tmp2.nc

#-----------------------------------------------
# Create file on the original 1x1 degree grid
#-----------------------------------------------
# Remove degenrate "lev" dimension and copy TFe_dep and LFe_dep into one file
ncap2 -A -s 'TFe[time,lat,lon]=double(TFe_dep(:,0,:,:))'  $tinfile tmp1.nc
ncap2 -A -s 'LFe[time,lat,lon]=double(LFe_dep(:,0,:,:))'  $linfile tmp2.nc
ncks  -A -C -x -v lev,TFe_dep tmp1.nc $outfile_ogrid
ncks  -A -C -x -v lev,LFe_dep tmp2.nc $outfile_ogrid


# Convert units to kg Fe/m2/s and specify a reasonable time axis
ncap2 -A -s 'TFe=TFe/1000.0/1000.0*12.0/365.0/86400.0' $outfile_ogrid $outfile_ogrid
ncap2 -A -s 'LFe=LFe/1000.0/1000.0*12.0/365.0/86400.0' $outfile_ogrid $outfile_ogrid
ncap2 -A -s 'time={15.5,45,74.5,105,135.5,166,196.5,227.5,258,288.5,319,349.5};' $outfile_ogrid $outfile_ogrid
ncap2 -A -s 'defdim("bnds",2);time_bnds[time,bnds]={ 0.,31.,31.,59.,59.,90.,90.,120.,120.,151.,151.,181.,181.,212.,212.,243.,243.,273.,273.,304.,304.,334.,334.,365.};'  $outfile_ogrid $outfile_ogrid

# Adjust attributes
ncatted -a units,TFe,o,c,'kg/m2/s'                         $outfile_ogrid
ncatted -a units,LFe,o,c,'kg/m2/s'                         $outfile_ogrid
ncatted -a long_name,TFe,o,c,'total iron deposition'       $outfile_ogrid
ncatted -a long_name,LFe,o,c,'soluble iron deposition'     $outfile_ogrid
ncatted -a max,TFe,d,,                                     $outfile_ogrid
ncatted -a min,TFe,d,,                                     $outfile_ogrid
ncatted -a max,LFe,d,,                                     $outfile_ogrid
ncatted -a min,LFe,d,,                                     $outfile_ogrid
ncatted -a units,time,o,c,'days since 1850-01-01 00:00:00' $outfile_ogrid
ncatted -a calendar,time,o,c,'noleap'                      $outfile_ogrid
ncatted -a bounds,time,o,c,'time_bnds'                     $outfile_ogrid


rm $outfile_mgrid tmp1.nc tmp2.nc

#------------------------------------------------
# Create file on model grid
#------------------------------------------------
# Prepare file on model grid using plon, plat from the grid-file
# and the modified time axis from the file on 'ogrid'
ncks -A -v plon,plat      -o $outfile_mgrid $gridfile
ncks -A -v time,time_bnds -o $outfile_mgrid $outfile_ogrid
ncatted -a corners,plon,d,,                 $outfile_mgrid
ncatted -a corners,plat,d,,                 $outfile_mgrid
ncap2 -A -s 'TFe[time,y,x]=0.0'             $outfile_mgrid
ncatted -a coordinates,TFe,c,c,'plon plat'  $outfile_mgrid
ncap2 -A -s 'LFe[time,y,x]=0.0'             $outfile_mgrid
ncatted -a coordinates,LFe,c,c,'plon plat'  $outfile_mgrid

# Bilinear interpolation to model grid
cdo remapbil,$outfile_mgrid  $outfile_ogrid $outfile_mgrid


