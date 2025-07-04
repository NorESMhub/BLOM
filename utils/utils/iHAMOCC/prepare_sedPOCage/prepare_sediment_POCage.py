import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def auto_detect_grid(ds):
  '''
  Auto-detect horizontal grid information
  We assume that x,y are always present in the (BLOM) restart file as lon,lat dimensions
  ''' 
  sx = ds.sizes['x']
  sy = ds.sizes['y']
  
  if sx == 180 and sy == 193:
    grid = 'tnx2'
  elif sx == 360 and sy == 385:
    grid = 'tnx1'
  elif sx == 720 and sy == 641:
    grid = 'tnx0.5'
  elif sx == 1440 and sy == 1153:
    grid = 'tnx0.25'
  elif sx == 2880 and sy == 2165:
    grid = 'tnx0.125_20200722'
  elif sx == 2880 and sy == 2161:
    grid = 'tnx0.125_20221013'
  else:
    # Currently unclear, which grid version refers to 0.125...  
    print('Grid undefined with lon, lat dimensions: '+str(sx)+' x ' +str(sy))
    print('Please update the data - exit now.')
    exit()
  
  print('Grid name successfully recognized and set to: '+grid)
  return grid

def assign_sed_dz(data):
    '''
    Double-checked sediment vertical extension of grid boxes
    '''
    ks     = 12
    dzs    = np.zeros(ks+1) 
        
    # Grid box extensions
    dzs[0] = 0.001
    dzs[1] = 0.003
    dzs[2] = 0.005
    dzs[3] = 0.007
    dzs[4] = 0.009
    dzs[5] = 0.011
    dzs[6] = 0.013
    dzs[7] = 0.015
    dzs[8] = 0.017
    dzs[9] = 0.019
    dzs[10] = 0.021
    dzs[11] = 0.023
    dzs[12] = 0.025
        
    # depth at pressure points:
    seddw    = np.zeros(ks)
    seddw[0] = 0.5*dzs[0]
    for k in range(1, ks):
        seddw[k] = 0.5 * (dzs[k] + dzs[k+1])
 
    sedz = np.cumsum(seddw)
    xrdzs = xr.DataArray(data=dzs[0:-1],dims=['ks'],
                         name='seddz',
                         attrs={'long_name':'Sediment grid box vertical extension','units':'m'})
    xrzs  = xr.DataArray(data=sedz,dims=['ks'],
                         name='sedz',
                         attrs={'long_name':'Sediment depth at pressure points','units':'m'})
        
    data = xr.merge([data,xrdzs.to_dataset(),xrzs.to_dataset()])
    return data  


#===================================================================================================
# File that holds climatlogical sedimentation fluxes (carflx_bot,calflx_bot,bsiflx_bot,dustflx_bot)
# The units in that file require to be:
#   - carflx_bot: mol C m-2 s-1 
#   - calflx_bot: mol Ca m-2 s-1
#   - bsiflx_bot: mol Si m-2 s-1
#   - dustflx_bot: g m-2 s-1

clim_file = 'flx_bot.nc'

# in case of testing this out for single column mode, 
# gridbox (1,1) will be set with values from xloc,yloc (usually masked, since Antarctic)
# NOTE: xloc,yloc currently randomly chosen in tnx2 grid (somewhere Pacific)
PREP_SINGLE_COLUMN = 0
xloc = 157
yloc = 90


#---------------------------------------------------------------------------------------------------
# Open the climatological dataset
ds   = xr.open_dataset(clim_file)
grid = auto_detect_grid(ds)

#det_molP2kgPOC = 3166./122./1000. # This would be the stoichiometry version
det_molC2kgPOC = 30./1000.  # kg/mol iHAMOCC: orgwei/(1000 mol/kmol)
rho_det        = 1000.      # kg/m3  iHAMOCC: orgdens
calc_mol2kg    = 100./1000. # kg/mol iHAMOCC: calcwei/(1000 mol/kmol)
rho_calc       = 2600.      # kg/m3  iHAMOCC: calcdens
opal_mol2kg    = 60./1000.  # kg/mol iHAMOCC: opalwei/(1000 mol/kmol)
rho_opal       = 2200.      # kg/m3  iHAMOCC: opaldens
dust_g2kg      = 1./1000.   # kg/g   
rho_dust       = 2600.      # kg/m3  iHAMOCC: claydens

# Calculating the total volume fluxes to sediment [m3/s]
if 'dustflx_bot' in list(ds.keys()):
  ds = ds.assign(FV_tot = ds.carflx_bot * det_molC2kgPOC / rho_det  
                        + ds.calflx_bot * calc_mol2kg    / rho_calc
                        + ds.bsiflx_bot * opal_mol2kg    / rho_opal
                        + ds.dustflx_bot* dust_g2kg      / rho_dust)
else:
  print('\n\n   !!! WARNING: dustflx_bot not found - overestimating minimum age\n')    
  ds = ds.assign(FV_tot = ds.carflx_bot * det_molC2kgPOC / rho_det  
                        + ds.calflx_bot * calc_mol2kg    / rho_calc
                        + ds.bsiflx_bot * opal_mol2kg    / rho_opal)

year = 86400.*360. # one year in seconds
ds   = assign_sed_dz(ds)

# We here calculate the minimum age of the individual sediment layers
# Minimum, since we do not account for dissolution and remineralization here
# In case that dustflx_bot is not available in the climatology file, 
# this increases the layer age
ds   = ds.assign(sedPOCage = ds.sedz/(ds.FV_tot.squeeze()*year))
ds   = ds.assign(prorca_mavg = ds.carflx_bot.squeeze()/1000./122.)
ds['sedPOCage'].attrs['long_name']   = 'Minimum_POC_age'
ds['sedPOCage'].attrs['units']       = 'yr'
ds['prorca_mavg'].attrs['long_name'] = 'Climatological_mean_POP_sedimentation_flux'
ds['prorca_mavg'].attrs['units']     = 'kmol P m-2 s-1'

# Saving the init file for sediment quality
if PREP_SINGLE_COLUMN:
  dss = ds.copy()
  dss['sedPOCage'].loc[{'x':0,'y':0}] = dss['sedPOCage'].loc[{'x':xloc,'y':yloc}]
  dss['prorca_mavg'].loc[{'x':0,'y':0}] = dss['prorca_mavg'].loc[{'x':xloc,'y':yloc}]

  fname = 'sedment_quality_init_'+grid+'_single_column.nc'
  ds[['sedPOCage','prorca_mavg']].to_netcdf(fname)
else:
  fname = 'sedment_quality_init_'+grid+'.nc'
  ds[['sedPOCage','prorca_mavg']].to_netcdf(fname)

# Plot age of last layer to see regions of highest sediment age
fig = plt.figure()
plt.pcolormesh(np.log10(ds.sedPOCage.loc[{'ks':ds.sedPOCage.sizes['ks']-1}].squeeze()) ,rasterized=True)
#plt.pcolormesh(ds.age.squeeze() ,rasterized=True)
cb=plt.colorbar()
cb.ax.set_ylabel(r'log(age (yr)) ')
plt.title(r'Max. sediment age based on iHAMOCC fluxes')
figname = 'sediment-residence-time.pdf'
fig.savefig(figname)

# Output some info to screen:
print('\n=======================================================================================\n')
print('The following data was saved to the output file:')
print(ds[['sedPOCage','prorca_mavg']])
print('\n=======================================================================================\n')
print('    Figure with minimum sediment residence time saved:   '+figname)
print('    The following input file for iHAMOCC is generated:   '+fname+'\n\n')


