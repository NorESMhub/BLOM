#!/bin/bash

module load CDO/2.0.6-gompi-2022a
module load NCO/5.1.9-iomkl-2022a
version=$(date '+%Y%m%d')

GRID="tnx0.5v1"
DUST_IN="dustdep_mhw2006_T42.nc"

basepath="/cluster/shared/noresm/inputdata/ocn/blom/grid/"

echo "$GRID"
if [ ${GRID} == 'tnx2v1' ]; then
  GRIDFILE="grid_${GRID}_20130206.nc"
elif [ ${GRID} == 'tnx1v4' ]; then
  GRIDFILE="grid_${GRID}_20170622.nc"
elif [ ${GRID} == 'tnx0.5v1' ]; then
  GRIDFILE="grid_tnx0.5v1_20240702.nc"
elif [ ${GRID} == 'tnx0.25v4' ]; then
  GRIDFILE="grid_tnx0.25v4_20170622.nc"
elif [ ${GRID} == 'tnx0.125v4' ]; then
  GRIDFILE="grid_tnx0.125v4_20221013.nc"
fi

GRIDFILE=${basepath}/${GRIDFILE}

#-----------------------------------------------------------------generate gridfile:
# get the respective variables:
ncks -A -v plon,plat,parea,pclon,pclat $GRIDFILE gridtmp.nc
# rename plon and plat
ncrename -v plon,lon -v plat,lat gridtmp.nc 
# add lon, lat variables tas coordinates to parea
ncatted -a coordinates,parea,c,c,"lon lat" gridtmp.nc
# rename attribute name corners to bounds (cf-conventions)
ncrename -a lon@corners,bounds -a lat@corners,bounds gridtmp.nc
    
# re-order coordinates to enable cdo to read/interpret them 
ncpdq -O --rdr=y,x,nv gridtmp.nc gridfile_${GRID}.nc  
    
# check the structure of the gridfile
cdo verifygrid gridfile_${GRID}.nc


cdo gencon,gridfile_${GRID}.nc ${DUST_IN} weights_${GRID}.nc

# unit conversion factor kg/m2/s to kg/m2/month: 30*86400 = 2592000
cdo -setattribute,DUST@units="kg/m2/month" \
    -setattribute,DUST@long_name="dust_deposition" \
    -setattribute,DUST@coordinates="plat plon"  \
    -setname,DUST -mulc,2592000 -remap,gridfile_${GRID}.nc,weights_${GRID}.nc $DUST_IN dustdep_mhw2006_${GRID}_${version}.nc

ncrename -h -O -v lat,plat dustdep_mhw2006_${GRID}_${version}.nc
ncrename -h -O -v lon,plon dustdep_mhw2006_${GRID}_${version}.nc
ncrename -h -O -v lon_bnds,plon_bnds dustdep_mhw2006_${GRID}_${version}.nc
ncrename -h -O -v lat_bnds,plat_bnds dustdep_mhw2006_${GRID}_${version}.nc
ncatted -a bounds,plat,o,c,"plat_bnds" dustdep_mhw2006_${GRID}_${version}.nc
ncatted -a bounds,plon,o,c,"plon_bnds" dustdep_mhw2006_${GRID}_${version}.nc

rm gridtmp.nc
rm weights_${GRID}.nc
rm gridfile_${GRID}.nc 
