import xarray as xr
import numpy as np
import scipy.interpolate as sc_interp

import cdo 
import os as os
import yaml as yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
import argparse as argparse
import socket

# own imports
from blom_eos import rho as blom_rho
from blom_eos import p_alpha as blom_p_alpha


cdo = cdo.Cdo()

# --------------------------------------------------------------------------------------------------
def get_socket():
  '''
  Get the socket name and return short form 
  '''
  socket_name  = socket.gethostname()
  # names:
  #   - betzy: login-1.betzy.sigma2.no
  #   - local: uib-d93nzh3
  #   - nird:  

  if 'betzy' in socket_name:
    return 'betzy'
   
  elif 'nird' in socket_name:
    return 'nird'

  elif 'ipcc' in socket_name:
    return 'nird'

  else:
    print('machine_settings.py: Socket name '+socket_name+' not found in list - exit now')
    exit()

# --------------------------------------------------------------------------------------------------
def check_file_existence(fil):
  if os.path.exists(fil):
    pass
  else:
    print('File: \n   '+fil+'\ndoes not exist - exit now.')
    exit()

# --------------------------------------------------------------------------------------------------
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
  else:
    print('Grid undefined with lon, lat dimensions: '+str(sx)+' x ' +str(sy))
    print('Please update the data - exit now.')
    exit()
  
  print('Grid name successfully recognized and set to: '+grid)
  return grid

# --------------------------------------------------------------------------------------------------
def auto_detect_units(ds,v):
  '''
  Auto-detect the units used in BLOM
  (check global mean of sigma)
  '''
  m = ds[v].mean().squeeze().values
  if m < 1:
    units = 'CGS'
  else: 
    units = 'MKS'
  print('Mean sigma: '+str(m)+'   - settings BLOM units to: '+units)
  return units 

# --------------------------------------------------------------------------------------------------
def grid_specs(grid,comp):
  # Naming conventions used in restart files
  machine = get_socket()
  match machine:
    case 'betzy':
      gridbasepath = '/cluster/shared/noresm/inputdata/ocn/blom/grid/'
    case 'nird':
      gridbasepath = '/projects/NS2980K/data/inputdata/ocn/blom/grid/'

  # NOTE: we currently make no use of the here given vdim names - up to change/to be improved
  grid_des = {'tnx2':{'gridfile':os.path.join(gridbasepath,'grid_tnx2v1_20130206.nc'),
                      'grid':   {'lat':'plat','lon':'plon'},
                      'blom':   {'lat':'lat',
                                 'lon':'lon',
                                 'latdim':'y',
                                 'londim':'x',
                                 'vdim': ['kk2'],
                                }, 
                       'hamocc':{'lat':'lat',
                                 'lon':'lon',
                                 'latdim':'lat',
                                 'londim':'lon',
                                 'vdim':['depth2']
                                },
                      },
               'tnx1':{'gridfile':os.path.join(gridbasepath,'grid_tnx1v4_20170622.nc'),
                      'grid':   {'lat':'plat','lon':'plon'},
                      'blom':   {'lat':'lat',
                                 'lon':'lon',
                                 'latdim':'y',
                                 'londim':'x',
                                 'vdim': ['kk2'],
                                }, 
                       'hamocc':{'lat':'lat',
                                 'lon':'lon',
                                 'latdim':'lat',
                                 'londim':'lon',
                                 'vdim':['depth2']
                                },
                      },
               'tnx0.5':{'gridfile':os.path.join(gridbasepath,'grid_tnx0.5v1_20240702.nc'),
                      'grid':   {'lat':'plat','lon':'plon'},
                      'blom':   {'lat':'lat',
                                 'lon':'lon',
                                 'latdim':'y',
                                 'londim':'x',
                                 'vdim': ['kk2'],
                                }, 
                       'hamocc':{'lat':'lat',
                                 'lon':'lon',
                                 'latdim':'lat',
                                 'londim':'lon',
                                 'vdim':['depth2']
                                },
                      },

               }

  return grid_des[grid]['gridfile'],grid_des[grid]['grid']['lon'],grid_des[grid]['grid']['lat'],   \
         grid_des[grid][comp]['lon'],grid_des[grid][comp]['lat'],                                  \
         grid_des[grid][comp]['londim'],grid_des[grid][comp]['latdim']

# --------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parsing the arguments (if provided)
    '''
    # enabling to show a quick usage information when called via: python myfun.py --help
    parser=argparse.ArgumentParser(description='iHAMOCC restart file remapping tool ...',          \
                                   epilog='... and some usage information.')
    parser.add_argument('yml',                                                                     \
                        nargs='?',                                                                 \
                        default='no_ymlfile_provided',                                             \
                        help='  - please provide a yml file' )
    args=parser.parse_args()
    if args.yml == 'no_ymlfile_provided':
      print('No yml-file provided - please provide one - see yml_template.yml - exit now')
      exit() 
    return args

# --------------------------------------------------------------------------------------------------
def read_yml(fil):
    '''
    Read a yml file
    '''
    if os.path.exists(fil):  
      yml = yaml.load(open(fil), Loader=Loader)
      #for key, value in yaml.load(open('test.txt'))['instance'].iteritems():
      #print key, value
      return yml
    else:
      print('yml_conf_file:\n    '+fil+'\ndoes not exist - please provide a valid file name in the config file\nexit now')
      exit()

# --------------------------------------------------------------------------------------------------
def _ufunc_1dvertical_interpolation_na(ds,depth_orig,depth_interp):
  # https://github.com/pydata/xarray/issues/3931
 
  # since there might be layers with NaN values and interp1d doesn't like those, 
  # we only chose finite values to interpolate/extrapolate the values to new pressure points 
  valid = np.isfinite(ds) 
  if valid.sum() > 1:
    f=sc_interp.interp1d(depth_orig[valid],ds[valid],fill_value='extrapolate')
    interp = f(depth_interp)
  elif valid.sum() == 1: 
    interp = np.ones(len(depth_interp))*ds[valid]
  else:
    interp = np.zeros(len(depth_interp))*np.nan
  return interp

# --------------------------------------------------------------------------------------------------
def xr_vertical_nan_fill(ds,v,mask,levels='dp',depth='depth',maskname='pmask'):
  '''
  UN-USED - now using  _ufunc_1dvertical_interpolation_na instead
  '''
  print('Vertical NaN filling via linear interpolation for ' + v + ' ... (could be speeded up via apply_ufunc usage)')
  da = ds[v].copy()
  co = ds[depth].values.astype('int')
  for i in range(0,ds.sizes['x']):
    for j in range(0,ds.sizes['y']):
       tmp = ds[v].isel(x=i,y=j)
       if(mask[maskname].isel(x=i,y=j) == 1): # operate only on ocean points
         tmp[depth] = ds[levels].isel(x=i,y=j)+np.linspace(0.0,0.001,ds.sizes[depth])
         tmp        = tmp.interpolate_na(dim=depth, method="linear", fill_value="extrapolate",use_coordinate=depth,limit=ds.sizes[depth]-2,keep_attrs=True)
         tmp[depth] = co
       da.loc[{'x':i,'y':j,depth:co}]=tmp  
  da = da.to_dataset()
  print('.... done')
  return da

# --------------------------------------------------------------------------------------------------
def get_P_mks2cgs(units):
  # pressure coefficient converting CGS to MKS
  match units:
    case 'MKS':
      P_mks2cgs = 1.
    case 'CGS':
      P_mks2cgs = 1e1
    case _:
      print('Units: '+units+' not known to get_P_mks2cgs - exit now')
      exit()
  return P_mks2cgs

# --------------------------------------------------------------------------------------------------
def drop_lonlat(new_restart,bgc_lon=None,bgc_lat=None):
  # drop the lon,lat variables to better conform to original iHAMOCC restart file
  try: 
    new_restart = new_restart.drop_vars([bgc_lon,bgc_lat])
  except:
    pass
  return new_restart

# --------------------------------------------------------------------------------------------------
def drop_depth(new_restart,depth='depth'):
  try: 
     new_restart = new_restart.drop_vars([depth])
  except:
     # In case of isopycnic interpolation, no depth variable is added - so just pass
     print(new_restart)
  return new_restart

# --------------------------------------------------------------------------------------------------
def store_var(new_restart,dstmp,var,v,vdim):
  if v in list(dstmp.keys()):
    # if variable in dstmp, the variable is now ready to be put into the final restart file
    var = xr.concat([dstmp[v],var[v]],vdim)
    new_restart = xr.merge([new_restart,var])
  else:
    # first processed time step
    dstmp = xr.merge([dstmp,var])
  return new_restart,dstmp

# --------------------------------------------------------------------------------------------------
def clean_tmpfiles(tmpfiles):
  # Clean up of temporary cdo files
  cdo.cleanTempDir()
  
  # Cleaning temporary files in the folder
  for f in tmpfiles:
    try: 
      print('Deleting temporary file: '+f)
      os.remove(f) 
    except:
      pass

# --------------------------------------------------------------------------------------------------
def rename_restart_dimensions(ds,bgc_londim=None,phy_londim=None,bgc_latdim=None,phy_latdim=None):
  new_restart = xr.Dataset()
        
  list_keys = list(ds.keys())
  for v in list_keys:
    # rename the dimensions (questinable, if needed in this way) 
    # - cannot be done at the dataset level(?!) - tried, but failed
    var = ds[v].rename({phy_londim:bgc_londim,phy_latdim:bgc_latdim}).to_dataset(name=v)
          
    # merge the datasets
    new_restart = xr.merge([new_restart,var])  

  return new_restart

# --------------------------------------------------------------------------------------------------
def reset_depth_dim(ds,coords,zdim='depth2',old='lev'):
  dsnew = xr.Dataset()
  list_keys = list(ds.keys())
  for v in list_keys:
     if zdim == old:
       var = ds[v].to_dataset()
     else: 
       # First rename the dimension name 
       var = ds[v].to_dataset().rename_dims({old:zdim})
     
       # Reset dimension
       var[zdim] = coords#.astype('int')
       

     # Now, the old, cdo-added coordinates can be dropped
     try: 
       var = var.drop_vars([old])
     except:
       pass

     dsnew = xr.merge([dsnew,var])
  return dsnew


# --------------------------------------------------------------------------------------------------
def level_interp(ds,mask,maskfile,mvar,bgc_londim=None,phy_londim=None,bgc_latdim=None,phy_latdim=None):
  # This interpolation can be on isopycnals or simply regular grid like in the sediment
  # Problem is that the format of restart file is not cdo-conform 
  # - hence, we try to circumnavigate this by 
  # new_restart = cdo.remapbil('remap_tnx1.nc',input=ds_tnx2_hamocc.rename_dims({'lon':'x','lat','y'})
  new_restart = xr.Dataset()
        
  list_keys = list(ds.keys())
  for v in list_keys:
    # Retrieve individual variable
    var = ds[v]
      
    # rename dimensions for remapping
    var = var.to_dataset().rename_dims({bgc_londim:phy_londim,bgc_latdim:phy_latdim})
    # add time
    #var = var.expand_dims(dim='time',axis=0)

    # ----- perform the remapping
    # First set land (=0) to missing values, fill with nearest neighbour and then remap bilinear
    var_nn       = cdo.setmisstonn(' -setctomiss,0 ',input=var,returnXArray=v).to_dataset()
    var_remapped = cdo.remapbil(maskfile,input=var_nn,returnXArray=v)
        
    # mask the land
    var_remapped = (var_remapped*mask[mvar]).to_dataset(name=v)

    # make sure that units and long name is provided - caution: currently newly added attribute: _FillValue
    var_remapped[v].attrs['units']         = ds[v].attrs['units']
    var_remapped[v].attrs['long_name']     = ds[v].attrs['long_name']
    var_remapped[v].attrs['missing_value'] = 99999. # somehow, missing value is lost in ds... - reset it manually here
     
    # merge the datasets
    new_restart = xr.merge([new_restart,var_remapped])  
  return new_restart 


# --------------------------------------------------------------------------------------------------
def adjust_wc_inventory(c,f,rho_c,rho_f,Vc,Vf,vdim_phy='kk2',vdim_bgc='depth2',ccoords=np.linspace(0,53,53),fcoords=np.linspace(0,53,53),bgc_londim='lon',bgc_latdim='lat',phy_londim='x',phy_latdim='y'):
  '''
  Adjusting the water column inventory to the one found in the coarse model version
 
  '''

  # Reset bgc dim names to physics dim names
  c   = rename_restart_dimensions(c,bgc_londim=phy_londim,phy_londim=bgc_londim,bgc_latdim=phy_latdim,phy_latdim=bgc_latdim)
  c   = reset_depth_dim(c,ccoords,zdim=vdim_phy,old=vdim_bgc)
  
  f   = rename_restart_dimensions(f,bgc_londim=phy_londim,phy_londim=bgc_londim,bgc_latdim=phy_latdim,phy_latdim=bgc_latdim)
  f   = reset_depth_dim(f,fcoords,zdim=vdim_phy,old=vdim_bgc)
  
  # Make sure the vertical dimensions correspond to each other 
  Vc      = reset_depth_dim(Vc.to_dataset(name='vol'),ccoords,zdim=vdim_phy,old='kk2')
  Vf      = reset_depth_dim(Vf.to_dataset(name='vol'),ccoords,zdim=vdim_phy,old='kk2')
  rho_c   = reset_depth_dim(rho_c.to_dataset(name='rho'),ccoords,zdim=vdim_phy,old='kk2')
  rho_f   = reset_depth_dim(rho_f.to_dataset(name='rho'),ccoords,zdim=vdim_phy,old='kk2')
 
  I_c  = (c*rho_c['rho']*Vc['vol']).sum() # Inventory in the original file
  I_f  = (f*rho_f['rho']*Vf['vol']).sum() # Inventory in the interpolated restart file
  fct  = I_c/I_f                          # Adjustment factor

  adjusted =  f * fct       # carry out the inventory adjustment

  # Provide some info to screen
  v = str(list(adjusted.keys())[0])
  fct_val = fct.to_array().squeeze().values
  print('Adjusted: '+v+' by factor: '+ str(fct_val))

  new_inv = (adjusted*rho_f['rho']*Vf['vol']).sum().to_array().squeeze().values
  old_inv_coarse = I_c.to_array().squeeze().values
  print('Check:    '+str(new_inv)+'   ==   '+str(old_inv_coarse)                                   \
        + '  mismatch: '+ str(new_inv - old_inv_coarse)                                            \
        + '  => '+  str(np.abs(100.*(new_inv - old_inv_coarse)/old_inv_coarse))+'%') 

  # Reverse dims resetting - change dim names back to original bgc dim names
  adjusted  = rename_restart_dimensions(adjusted,bgc_londim=bgc_londim,phy_londim=phy_londim,bgc_latdim=bgc_latdim,phy_latdim=phy_latdim)
  adjusted  = reset_depth_dim(adjusted,fcoords,zdim=vdim_bgc,old=vdim_phy)
  
  return adjusted,fct_val


# --------------------------------------------------------------------------------------------------
def bgc_rho(ds,dp='dp',T='temp',S='saln',units='CGS'):
  '''
  Calculate the used bgc rho in iHAMOCC - following mo_intfcblom
  '''
  # Calculating the pressure following blom2hamocc 

  match units: 
    case 'MKS':
     rho0      = 1e3 # Reference value for density [kg/m3]
     #P_mks2cgs = 1. # Pressure coefficient converting CGS to MKS
    case 'CGS':
     rho0      = 1. # Reference value for density [g/cm3]
     #P_mks2cgs = 1.e1 # Pressure coefficient converting CGS to MKS
    case _:
      print('Unknown units in bgc_rho:  '+units+'    -exit now')
      exit()

  P_mks2cgs = get_P_mks2cgs(units)

  # Calculate the pressure at individual layers (top = 0, first = dp,...)
  rp  = ds[dp].cumsum(dim='kk2') - ds[dp] 
 
  # Get ldp - differentiate between case 1 and 2+3 here
  ldp = xr.where(ds[dp] == 0.,1,ds[dp])

  # Calculate density according to BLOM
  rho = blom_rho(rp,ds[T],ds[S],units=units)
  pa  = ldp/rho # case 1+2 in iHAMOCC
  
  # Take care about the third case in blom2hamocc
  p2  = rp + ds[dp] 
  pa2 = blom_p_alpha(rp,p2,ds[T],ds[S])

  # Merge all the three cases together
  pa = xr.where(ds[dp] < 1e-3*P_mks2cgs, pa, pa2)

  bgcrho = ldp/(pa*rho0)
  
  return bgcrho, pa

# --------------------------------------------------------------------------------------------------
def gridvolume(grid,dz,area='parea'):
  vol = dz*grid[area]
  return vol  

# --------------------------------------------------------------------------------------------------
def dp2dz(pa,units='CGS'):
  '''
   Calculate the vertical extension of the sigma pressure layer in terms of units of meters
  '''
  match units:
    case 'MKS':
      rho0 = 1e3
      onem = 9806. # 1 m in units of pressure [kg m-1 s-2].
    case 'CGS':
      rho0 = 1.
      onem = 98060. # 1 m in units of pressure [g cm-1 s-2].
    case _:
      print('Unknown units in dp2dz:  '+units+'    -exit now')
      exit()

  dz = rho0*pa/onem

  return dz


# --------------------------------------------------------------------------------------------------
def dp2p(ds,v='dp',vdim='kk2',sstart=0,send=53):
  '''
  Sum up the delta pressure values to get the absolute pressure levels
  '''
  ds = ds[v][{vdim:slice(sstart,send)}].to_dataset()
  p  = ds[v].cumsum(dim=vdim)
  with xr.set_options(keep_attrs=True): 
    # since the above may set land to masked values, we set the pressure on land to zero for now 
    p  =xr.where(np.isfinite(p),p,0.)
  return p

# --------------------------------------------------------------------------------------------------
def dp2mid_layer_pressure(ds,v='dp',vdim='kk2',sstart=0,send=53,P_cgs2mks=1.):
  '''
  Calculate the absolute pressure at the center of the grid cells
  ''' 
  ds     = ds[v][{vdim:slice(sstart,send)}].to_dataset()

  # calculate mid-layer pressure and convert units to MKS
  mid_p  = (ds[v].cumsum(dim=vdim) - 0.5*ds[v]) * P_cgs2mks
  with xr.set_options(keep_attrs=True): 
    # since the above may set land to masked values, we set the pressure on land to zero for now 
    mid_p  = xr.where(np.isfinite(mid_p),mid_p,0.)
  return mid_p.squeeze().to_dataset()

# --------------------------------------------------------------------------------------------------
def add_grid_info(ds,grid,mlon='lon',mlon_dim='x',mlat='lat',mlat_dim='y',glon='plon',glat='plat'):
  '''
  add grid information for remapping, etc.
  '''
  try:
    ds  = ds.assign_coords({mlat:([mlat_dim,mlon_dim], grid[glat].values.squeeze()),
                            mlon:([mlat_dim,mlon_dim], grid[glon].values.squeeze())})
    # setting the attributes is key for later remapping!
    ds[mlon].attrs['units']='degrees_east'
    ds[mlat].attrs['units']='degrees_north'
  except:
    print('\n!!! Something went wrong when trying to attach lon/lat information to dataset !!!')
    print('!!! maybe the gridfile and the loaded dataset are not compatible ')
    print('!!! - e.g. different lon/lat dims x,y')
    print('!!! --- grid dims information')
    print(grid.dims)
    print('!!! --- dataset dims information')
    print(ds.dims)
    print('\n   - exit now')
    exit()
  return ds

