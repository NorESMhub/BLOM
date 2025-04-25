import xarray as xr
import numpy as np
import cdo 
import os
from validate import Validator

cdo = cdo.Cdo()

from iHAMOCC_restart_remapping_functions import dp2p,dp2mid_layer_pressure, add_grid_info,         \
                  grid_specs, get_P_mks2cgs, gridvolume, bgc_rho, dp2dz, level_interp,             \
                  adjust_wc_inventory, rename_restart_dimensions, reset_depth_dim, drop_lonlat,    \
                  drop_depth,store_var,clean_tmpfiles, parse_args,read_yml,auto_detect_units,      \
                  auto_detect_grid,check_file_existence,xr_vertical_nan_fill,                      \
                  _ufunc_1dvertical_interpolation_na

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
def main():

  # ===== START OF DEFINITIONS

  inputs = parse_args()

  yml = read_yml(inputs.yml)
  for item in yml.keys():
    print('Reading yml-settings: '+item+'    '+str(yml[item]))
  
  vtor = Validator()

  # Getting restart file information from the yml file
  coarse_blom   = vtor.check('string',yml['coarse_blom'])       # Coarse resolution BLOM restart file
  coarse_hamocc = vtor.check('string',yml['coarse_hamocc'])     # Coarse resolution iHAMOCC restart file
  fine_blom     = vtor.check('string',yml['fine_blom'])         # BLOM restart file to which the iHAMOCC coarse file shall be remapped

  # InterpCase provides different ways to interpolate the data from tnx2 to tnx1
  #  - isopycnic_interp:   just do an interpolation along the isopycnic layers
  #  - sigma_plvl_interp:  interpolate pressure and data to fine grid and then vertically interpolate 
  #                        to different pressure levels of the fine sigma coordinates
  #  - levitus_plvl_inter: interpolate to levitus pressure/depth levels, interpolate horizontally,
  #                        and performa vertical interpolation to isopycnal pressure levels
  #  - skip_interp:        enables to pass the previously written new restart file directly 
  #                        on to the AdjustInventory routine, if needed 
  #                        Requires:
  #                            - previously written new restart file
  #                            - coarse HAMOCC file information
  #                            - coarse BLOM file information
  #                            - fine BLOM file information
  
  InterpCase = vtor.check('string',yml['InterpCase'])

  # Perform an inventory adjustment (e.g. thickness of isopycnals can vary, leading to different inventories
  AdjustInventory = vtor.check('boolean',yml['AdjustInventory'])
  

  # The input restart file that will be adjusted, in case that no interpolation is carried out
  # else, the newly generated restart file is adjusted 
  restart_file_to_adjust  = vtor.check('string',yml['restart_file_to_adjust'])
      

  # The time step, when the diagnostic iHAMOCC variables were written 
  # As a simple test, remap to original BLOM grid and check the differences
  # They look minor and the shock due to regridding is likely more than the shock due 
  # wrong time step for the diagnostics (while they seem to be used unchanged after restart
  # which makes it likely, that always the same time step is at the end of the year reached)    
  diag_tstep  = vtor.check('integer',yml['diag_tstep'])
  
  #-----------------------------------------------------
  # NOTE: everything below here is not rigorously used - up to improvement 
  # - some basic setting, could enter the grid definition at some point
  coarse_blom_pvar   =  'dp'    
  fine_blom_pvar     =  'dp'     

  # Naming convention used in the grid file
  grid_lat           = 'plat'
  grid_lon           = 'plon'
  maskname           = 'pmask'
 
  # NOTE: The following could be used to allow further generalization in the vertical coordinates/dims
  #       Currently, a number of vertical dim names are hard-coded. The following are currently not used: 
  cblom_wcz  =  'kk'
  cblom_wcz2 = 'kk2'
  fblom_wcz  =  'kk'
  fblom_wcz2 = 'kk2'
 
  chamocc_wcz   = 'depth'
  chamocc_wcz2  = 'depth2'
  chamocc_sedz  = 'nks'
  chamocc_sedz2 = 'nks2' 
  chamocc_bt    = 'tlvl2'
  fhamocc_wcz   = 'depth'
  fhamocc_wcz2  = 'depth2'
  fhamocc_sedz  = 'nks'
  fhamocc_sedz2 = 'nks2' 
  fhamocc_bt    = 'tlvl2'

  # ===== END OF DEFINITIONS - BELOW START THE MAIN PROGRAM

  #=================================================================================================
  #=================================================================================================
  #=================================================================================================
  # Check if the provided files exist
  check_file_existence(coarse_blom)
  check_file_existence(coarse_hamocc)
  check_file_existence(fine_blom)
  
  # Open BLOM/iHAMOCC restart files
  ds_coarse_blom   = xr.open_dataset(coarse_blom)
  ds_coarse_hamocc = xr.open_dataset(coarse_hamocc) #[['satoxy','hi']]
  ds_fine_blom     = xr.open_dataset(fine_blom)

  # Simple auto-detection of units used in BLOM 
  coarse_blom_units = auto_detect_units(ds_coarse_blom,'sigma')
  fine_blom_units   = auto_detect_units(ds_fine_blom,'sigma')

  # Auto-detect the grid name (tnx2, tnx1,...) - potentially to be updated
  coarse_grid_name = auto_detect_grid(ds_coarse_blom)
  fine_grid_name   = auto_detect_grid(ds_fine_blom)
  
  adjustable_restart_file = None # when no restart file interpolation is carried out, this will be set
                                 # to restart_file_to_adjust else the current new provided restart 
                                 # file is adjusted 
  provide_new_restart     = True # default

  # Get some grid specs
  coarse_gridfile, coarse_grid_lon,coarse_grid_lat,                                                \
  coarse_blom_lon,coarse_blom_lat, coarse_blom_londim,coarse_blom_latdim = grid_specs(coarse_grid_name,'blom')
  coarse_gridfile, coarse_grid_lon,coarse_grid_lat,                                                \
  coarse_hamocc_lon,coarse_hamocc_lat, coarse_hamocc_londim,coarse_hamocc_latdim = grid_specs(coarse_grid_name,'hamocc')
  
  fine_gridfile, fine_grid_lon, fine_grid_lat,                                                     \
  fine_blom_lon, fine_blom_lat, fine_blom_londim, fine_blom_latdim = grid_specs(fine_grid_name,'blom')
  fine_gridfile, fine_grid_lon, fine_grid_lat,                                                     \
  fine_hamocc_lon, fine_hamocc_lat, fine_hamocc_londim, fine_hamocc_latdim = grid_specs(fine_grid_name,'hamocc')

  # Conversion factor. NOTE: for vertical interpolation, we define the reference pressures in Pa (kg/(m2 s))
  coarse_P_cgs2mks = 1./get_P_mks2cgs(coarse_blom_units)
  fine_P_cgs2mks   = 1./get_P_mks2cgs(fine_blom_units)

  # Collect info on temporary files to delete them at the end
  tmpfiles = [] 

  # Define a mask file name that will be written
  out_maskfile_fine  = 'maskfile_'+fine_grid_name+'.nc'
  tmpfiles.append(out_maskfile_fine)

  # Define new restart file name:
  out_new_restart_file = 'new_restart-HAMOCC_'+fine_grid_name+'_'+InterpCase+'.nc'
   
  # Check gridfile existence and open gridfiles 
  check_file_existence(coarse_gridfile)
  check_file_existence(fine_gridfile)
  grid_coarse      = xr.open_dataset(coarse_gridfile)
  grid_fine        = xr.open_dataset(fine_gridfile)

  # Get a mask field for masking
  mask_coarse    = grid_coarse[maskname].to_dataset()
  mask_fine      = grid_fine[maskname].to_dataset()

  # Print basic info on dataset
  print('=============================================================================== '+coarse_grid_name+' BLOM')
  print(ds_coarse_blom)
  print('============================================================================= '+coarse_grid_name+' HAMOCC')
  print(ds_coarse_hamocc)
  print('=============================================================================== '+fine_grid_name+' BLOM')
  print(ds_fine_blom)


  # Get number of vertical layer information from dimensions 
  N_coarse_wclayers  = ds_coarse_blom.sizes['kk'] # Number of sigma layers in the water column
  N_fine_wclayers    = ds_fine_blom.sizes['kk']
  try: 
    N_coarse_sedlayers = ds_coarse_hamocc.sizes['nks'] # Number of sediment layers - we currently do not support different sediment layer numbers
  except:
    pass  
  # =============================== add grid information (lon, lat info)
  ds_coarse_blom   = add_grid_info(ds_coarse_blom,grid_coarse,mlon=coarse_blom_lon,mlon_dim=coarse_blom_londim,\
                                 mlat=coarse_blom_lat,mlat_dim=coarse_blom_latdim,glon=coarse_grid_lon,glat=coarse_grid_lat)
  ds_coarse_hamocc = add_grid_info(ds_coarse_hamocc,grid_coarse,mlon=coarse_hamocc_lon,mlon_dim=coarse_hamocc_londim,\
                                 mlat=coarse_hamocc_lat,mlat_dim=coarse_hamocc_latdim,glon=coarse_grid_lon,glat=coarse_grid_lat)
  ds_fine_blom   = add_grid_info(ds_fine_blom,grid_fine,mlon=fine_blom_lon,mlon_dim=fine_blom_londim,\
                                 mlat=fine_blom_lat,mlat_dim=fine_blom_latdim,glon=fine_grid_lon,glat=fine_grid_lat)
  mask_fine      = add_grid_info(mask_fine,grid_fine,mlon=fine_blom_lon,mlon_dim=fine_blom_londim,\
                                 mlat=fine_blom_lat,mlat_dim=fine_blom_latdim,glon=fine_grid_lon,glat=fine_grid_lat)

  # Write out maskfile for fine grid - for regridding purposes 
  mask_fine.to_netcdf(out_maskfile_fine)

  # ===================================================== Start with the remapping and interpolation
  match InterpCase:
    # ================================================================================================
    # Three ways to check out:
    #    1. interpolate data to tnx1 grid - assuming that concentrations are mainly driven by along isopycnic
    #       flows/mixing, less by gravitational sinking and cross-isopycnic mixing, 
    #       assuming same number and sigma of vertical layers, potential issues with empty layers
    #    2. interpolate tnx2 plevels to tnx1 grid, interpolate data to tnx1 grid,  then remap 3D to tnx1 levels
    #       this remapping assumes that ocean depth (and thus vertical mixing+gravitational sinking)
    #       plays similarly a big role as along isopycnic transport
    #    3. interpolate tnx2 data first to regular p-levels, then interpolate data, and remap to 3D tnx1 levels
    #       similar to 2., but slightly more in line with the internal remapping of BLOM/HAMOCC for output
    #       while a higher numerical diffusion is allowed through linear interpolation 

    #-----------------------------------------------------------------------------------------------
    case 'isopycnic_interp': # Refers to case 1. in the description above 
      print('================================================ Perform InterpCase: isopycnic_interp')
      new_restart = level_interp(ds_coarse_hamocc,mask_fine,out_maskfile_fine,maskname,        \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
      new_restart = rename_restart_dimensions(new_restart,                                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
      print('\n    Could be done/improved via:')
      print('\n      - Vertical interpolation where NaN values instead of misstonn filling')
      
    #-----------------------------------------------------------------------------------------------
    case 'sigma_plvl_interp': # Refers to case 2. in the description above 
      print('=============================================== Perform InterpCase: sigma_plvl_interp')
      #=================================================================================================
      # Explanations for cdo interp3D:
      # ===============================
      # cdo intlevel3d,icoordinate infile1 infile2 outfile
      # icoordinate contains a single 3d variable, which represents the target 3d vertical coordinate
      # infile1     contains the source data, which the vertical coordinate from icoordinate belongs to
      # infile2     only contains the source 3d height levels
      # intlevels3d requires same lon/lat grid size and montonically increasing depth axis
      #=================================================================================================
      
      # Init empty (temporary) datasets           
      new_restart = xr.Dataset()
      depth2_tmp  = xr.Dataset()  
      nks2_tmp    = xr.Dataset()  
      tlvl2_tmp   = xr.Dataset() 
      for tstep in range(0,2): 
        # Since we deal with two time steps, where pressure levels are close to each other 
        # and we interpolate in the vertical, we need to split up the data into those time steps
        # which are concatenated along the vertical dimension in BLOM, while not always in iHAMOCC,
        # where sometimes two time steps exist
        
        # define the slicing in the vertical dimension:
        coarse_sstart = tstep     * N_coarse_wclayers
        coarse_send   = (tstep+1) * N_coarse_wclayers
        fine_sstart   = tstep     * N_fine_wclayers
        fine_send     = (tstep+1) * N_fine_wclayers
        try: 
          coarse_sstart_sed = tstep     * N_coarse_sedlayers
          coarse_send_sed   = (tstep+1) * N_coarse_sedlayers
        except:
          pass
        # =================== Generate absolute pressure level values 
        # Calculate the absolute pressure at the centers of the pressure layers, 
        # since concencentrations are centered at pressure points
        # - !!!! single depth grid here (at one time point) !!!!!
        plevels_coarse = dp2mid_layer_pressure(ds_coarse_blom,v=coarse_blom_pvar,vdim='kk2',       \
                                 sstart=coarse_sstart,send=coarse_send,P_cgs2mks=coarse_P_cgs2mks )
        plevels_fine   = dp2mid_layer_pressure(ds_fine_blom,v=fine_blom_pvar,vdim='kk2',           \
                                 sstart=fine_sstart,send=fine_send,P_cgs2mks=fine_P_cgs2mks)

        # Save the pressure levels for interpolation purposes
        plevels_coarse_file = 'plevels_'+coarse_grid_name+'_coarse_tstep-'+str(tstep)+'.nc'
        plevels_fine_file   = 'plevels_'+fine_grid_name+'_fine_tstep-'+str(tstep)+'.nc'
        plevels_coarse.squeeze().to_netcdf(plevels_coarse_file)
        plevels_fine.squeeze().to_netcdf(plevels_fine_file)
        tmpfiles.append(plevels_coarse_file)
        tmpfiles.append(plevels_fine_file)
        
        # Interpolate the coarse resolution pressures to finer grid 
        pinterp_file = 'plevels_'+coarse_grid_name+'-to-'+fine_grid_name+'_tstep-'+str(tstep)+'_interp.nc'
        pinterp      = cdo.remapbil(out_maskfile_fine,input=plevels_coarse_file,returnXArray='dp').to_dataset() 
        pinterp.squeeze().to_netcdf(pinterp_file)
        tmpfiles.append(pinterp_file)

        # intlevel3D requires monotonically increasing pressure/depth levels, 
        # which is not the case on land - hence, we add here some epsilon value to both, 
        # coarse and fine resolution pressure levels
        plevels_coarse_adjusted_file = 'plevels_'+coarse_grid_name+'_tstep-'+str(tstep)+'_coarse_peps-adjusted.nc'
        pinterp_adjusted_file      = 'plevels_'+coarse_grid_name+'-to-'+fine_grid_name             \
                                      +'_tstep-'+str(tstep)+'_interp_peps-adjusted.nc'
        plevels_fine_adjusted_file = 'plevels_'+fine_grid_name+'_tstep-'+str(tstep)+'_fine_peps-adjusted.nc'
        coarse_epspress = np.linspace(0.001,0.1,N_coarse_wclayers) 
        interp_peps     = xr.Dataset(data_vars=dict(p_eps = (['kk2'],coarse_epspress,{'units':''})))
        fine_epspress   = np.linspace(0.001,0.1,N_fine_wclayers) 
        fine_peps       = xr.Dataset(data_vars=dict(p_eps = (['kk2'],fine_epspress,{'units':''})))
        
        # On land, the pressure levels need to be monotontically rising,
        # so we add the peps only on land points      
        coarse_peps = interp_peps*xr.where(mask_coarse[maskname]==0,1,0)
        interp_peps = interp_peps*xr.where(mask_fine[maskname]==0,1,0)
        fine_peps   = fine_peps*xr.where(mask_fine[maskname]==0,1,0)
        
        ( plevels_coarse  + coarse_peps['p_eps'] ).to_netcdf(plevels_coarse_adjusted_file)
        ( pinterp      + interp_peps['p_eps'] ).to_netcdf(pinterp_adjusted_file)
        ( plevels_fine + fine_peps['p_eps']   ).to_netcdf(plevels_fine_adjusted_file)
        tmpfiles.append(pinterp_adjusted_file)
        tmpfiles.append(plevels_coarse_adjusted_file)
        tmpfiles.append(plevels_fine_adjusted_file)

        # ========== Start variable interpolation - each variable individually
        # Separate the water column and the sediment/burial via dimension name
        # For sediment/burial variables, just interpolate
        # For water column, we do the 3D vertical remapping
        # Header of iHAMOCC restart file:
        # netcdf orig_tnx2_hamocc {
        # dimensions:
	#    lon = 180 ;
	#    lat = 193 ;
	#    depth = 53 ;    # water column one time step
	#    depth2 = 106 ;  # water column two time steps
	#    nks = 12 ;      # sediment one time step
	#    nks2 = 24 ;     # sediment two time steps
	#    tlvl2 = 2 ;     # burial, two time steps
        # NOTE: we still need to figure out, which tstep to use for variables with depth=53

        lkeys = list(ds_coarse_hamocc.keys())
        for v in lkeys:
          # Remapping the individual variables for individual time steps
          print('Remapping t='+str(tstep)+' ' + v +'   '+str(ds_coarse_hamocc[v].dims)+ '   '      \
                                                        +str(ds_coarse_hamocc[v].attrs))
          if 'depth2' in ds_coarse_hamocc[v].dims:
            # BGC water column data with two time levels
            # Interpolate BGC data to finer grid
            lev_interp = level_interp(ds_coarse_hamocc[v].isel(depth2=slice(coarse_sstart,coarse_send)).to_dataset(),\
                                        mask_fine,out_maskfile_fine,maskname,\
                                        bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim, \
                                        bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
            
            # Write interpolated variable to tmp-file for subsequent veryical interpolation
            tmp = 'tmp-'+v+'_tstep-'+str(tstep)+'.nc'
            lev_interp.to_netcdf(tmp)
            tmpfiles.append(tmp)

            # perform a vertical interpolation, followed by negative and NaN-nearest neighbor filling
            # to avoid e.g. issues at the coast
            #var = cdo.setmisstonn(' -setrtomiss,-9999999.,0. -intlevel3d,'+plevels_fine_adjusted_file  \
            #                      + ' selvar,'+v,input=cdo.copy(lev_interp, returnCdf=True)+' '   \
            #                      + pinterp_adjusted_file,returnXArray=v)
            var = cdo.setmisstonn(' -setrtomiss,-99999999999.,0. -intlevel3d,'+plevels_fine_adjusted_file,\
                            input=' -selvar,'+v +' '+tmp+' '+pinterp_adjusted_file,returnXArray=v)
            var = (var*xr.where(mask_fine[maskname]==0,np.nan,1)).to_dataset(name=v) # masking data with land
            var = rename_restart_dimensions(var,                                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
            var = reset_depth_dim(var,np.linspace(fine_sstart,fine_send-1,N_fine_wclayers),zdim='depth2',old='lev') 
            # make sure that units and long name are provided - caution: currently newly added attribute: _FillValue
            var[v].attrs['units']         = ds_coarse_hamocc[v].attrs['units']
            var[v].attrs['long_name']     = ds_coarse_hamocc[v].attrs['long_name']
            var[v].attrs['missing_value'] = 99999. # somehow, missing value is lost in ds - reset it manually here
            new_restart,depth2_tmp = store_var(new_restart,depth2_tmp,var,v,'depth2')                 

          elif 'depth' in ds_coarse_hamocc[v].dims:
            print('!!!! depth: Requires time step check !!!!!')
            if diag_tstep == tstep:
              # Prepare tmp variable for vertical NaN-interpolation 
              # - fill range -inf to 1e-8 (currently smallest value for diagnostic vars)
              if v == 'hi':
                minval=0.
              else:
                minval=1e-8
              dstmp = xr.where(ds_coarse_hamocc[v]>minval,ds_coarse_hamocc[v],np.NaN).to_dataset()
              tmp_p = plevels_coarse  + coarse_peps['p_eps']
              tmp_p = rename_restart_dimensions(tmp_p,                                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
              tmp_p = reset_depth_dim(tmp_p,np.linspace(0,N_coarse_wclayers-1,N_coarse_wclayers),zdim='depth',old='kk2') 
              dstmp = xr.merge([dstmp,tmp_p])
              dstmp = rename_restart_dimensions(dstmp,                                         \
                                     bgc_londim=coarse_blom_londim,phy_londim=coarse_hamocc_londim,\
                                     bgc_latdim=coarse_blom_latdim,phy_latdim=coarse_hamocc_latdim)
              dstmp=dstmp.rename({'y':'lat','x':'lon'})
              dstmp=dstmp.rename_dims({'lat':'y','lon':'x'})
              
              # Perform the vertical interpolation - fill range -inf to 1e-8 (currently smallest value for diagnostic vars)
              # dstmp = xr_vertical_nan_fill(dstmp,v,mask_coarse,levels='dp')
              # For each input field to applied ufunc, one needs to define the core dims
              dstmp = xr.apply_ufunc(_ufunc_1dvertical_interpolation_na,dstmp[v],dstmp['dp'],dstmp['dp'],\
                                     input_core_dims=[['depth'],['depth'],['depth']],                    \
                                     output_core_dims=[['depth']],vectorize=True)
              dstmp=dstmp.to_dataset(name=v)
              dstmp = dstmp.transpose('depth','y','x')
              dstmp = dstmp.rename({'y':'lat','x':'lon'})
              dstmp = drop_depth(dstmp,depth='depth')
              dstmp.lon.attrs['units'] = ds_coarse_hamocc.lon.attrs['units']
              dstmp.lat.attrs['units'] = ds_coarse_hamocc.lat.attrs['units']
              dstmp[v].attrs['units']  = ds_coarse_hamocc[v].attrs['units']
              dstmp[v].attrs['long_name']  = ds_coarse_hamocc[v].attrs['long_name']
              # Interpolate BGC data to finer grid
              #lev_interp = level_interp(ds_coarse_hamocc[v].to_dataset(),\
              lev_interp = level_interp(dstmp,\
                                        mask_fine,out_maskfile_fine,maskname,\
                                        bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim, \
                                        bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
              tmp = 'tmp-'+v+'_tstep-'+str(tstep)+'.nc'
              lev_interp.to_netcdf(tmp)
              tmpfiles.append(tmp)
              
              # perform a vertical interpolation, followed by negative and NaN-neartest neighbor filling
              # to avoid e.g. issues at the coast
              var = cdo.setmisstonn(' -setrtomiss,-99999999999.,0. -intlevel3d,'+plevels_fine_adjusted_file,\
                              input=' -selvar,'+v +' '+tmp+' '+pinterp_adjusted_file,returnXArray=v)
              var = (var*xr.where(mask_fine[maskname]==0,np.nan,1)).to_dataset(name=v) # masking data with land
              var = rename_restart_dimensions(var,                                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
              var = reset_depth_dim(var,np.linspace(0,N_fine_wclayers-1,N_fine_wclayers),zdim='depth',old='lev') 
              # make sure that units and long name are provided - caution: currently newly added attribute: _FillValue
              var[v].attrs['units']         = ds_coarse_hamocc[v].attrs['units']
              var[v].attrs['long_name']     = ds_coarse_hamocc[v].attrs['long_name']
              var[v].attrs['missing_value'] = 99999. # somehow, missing value is lost in ds - reset it manually here
              new_restart = xr.merge([new_restart,var])

          elif 'nks2' in ds_coarse_hamocc[v].dims:
            lev_interp = level_interp(ds_coarse_hamocc[v].isel(nks2=slice(coarse_sstart_sed,coarse_send_sed)).to_dataset(),\
                                       mask_fine,out_maskfile_fine,maskname,\
                                       bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim, \
                                       bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
            var = rename_restart_dimensions(lev_interp,                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
            new_restart,nks2_tmp = store_var(new_restart,nks2_tmp,var,v,'nks2')

          elif 'nks' in ds_coarse_hamocc[v].dims:
            print('!!!! nks: Requires time step check !!!!!')
            if diag_tstep == tstep: 
               lev_interp = level_interp(ds_coarse_hamocc[v].to_dataset(),\
                                         mask_fine,out_maskfile_fine,maskname,\
                                         bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim, \
                                         bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
               var = rename_restart_dimensions(lev_interp,                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
               new_restart = xr.merge([new_restart,var])

          elif 'tlvl2' in ds_coarse_hamocc[v].dims:
            lev_interp = level_interp(ds_coarse_hamocc[v].isel(tlvl2=tstep).to_dataset(),\
                                       mask_fine,out_maskfile_fine,maskname,\
                                       bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim, \
                                       bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
            var = rename_restart_dimensions(lev_interp,                         \
                                     bgc_londim=coarse_hamocc_londim,phy_londim=coarse_blom_londim,\
                                     bgc_latdim=coarse_hamocc_latdim,phy_latdim=coarse_blom_latdim)
            new_restart,tlvl2_tmp = store_var(new_restart,tlvl2_tmp,var,v,'tlvl2')
            
          else: 
            print('===== Unknown dimensions for variable: '+v)
            print('Some info on the variable:')
            print(ds_coarse_hamocc[v])
            print('===== exit now')
            exit()
        # end v in lkeys
      #end tstep      
    # end case sigma_plvl_interp

    #-----------------------------------------------------------------------------------------------
    case 'levitus_plvl_interp': # Refers to case 3. in the description above 
      print('============================================= Perform InterpCase: levitus_plvl_interp')
      print('Not yet implemented - please try another InterpCase - exit now')
      exit()
 
    #-----------------------------------------------------------------------------------------------
    case _: # default case
      print('InterpCase: '+InterpCase+' undefined - pass on to adjusting the restart file inventory, if requested')
      provide_new_restart     = False 
      adjustable_restart_file = restart_file_to_adjust 
      

  # end match case ---------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------
  if provide_new_restart:
    # remove lon, lat, depths variables from new restart file ds
    new_restart = drop_lonlat(new_restart,bgc_lon=coarse_hamocc_londim,bgc_lat=coarse_hamocc_latdim)
    new_restart = drop_depth(new_restart,depth='depth2')
    new_restart = drop_depth(new_restart,depth='depth')
  
    # write global attributes:
    new_restart.attrs = ds_coarse_hamocc.attrs
    new_restart.attrs['history'] = 'Original HAMOCC restart file: '+coarse_hamocc+' remapped to BLOM '\
                                   + fine_blom                                                     \
                                   + ' pressure levels using 3Dinterp_restart_remapping.py; '      \
                                   + 'coarse pressure level file used: '+coarse_blom+'; '          \
                                   + 'coarse gridfile used: '+coarse_gridfile+'; '                 \
                                   + 'fine gridfile used:'+fine_gridfile+'; '                      \
                                   + 'Original file history: '+ds_coarse_hamocc.attrs['history']
 
    # ----------------------------------------- Eventually write the restart file
    new_restart.to_netcdf(out_new_restart_file,format="NETCDF3_64BIT")
  
    # Clean up temporary files
    clean_tmpfiles(tmpfiles)


    if coarse_grid_name == fine_grid_name and  coarse_blom  == fine_blom:
      # Perform a difference between the original and the new file
      # Some minor differences are expected due to monotinicity issue
      print('\n====================================================================================')
      print('\nIt looks as if the old and new grid as well as the pressure information files are equal')
      print('This allows the check of the goodness of interpolation, remapping, etc.                ')
      print(' - performing this now')
      to_check  = new_restart - ds_coarse_hamocc
      checkfile = 'check_'+ coarse_grid_name+'-remapped-to-'+fine_grid_name+'-delta.nc'
      to_check.to_netcdf(checkfile)
      print('\n    Check the file:    '+checkfile+'    for differences')
      print('Expect some minor differences due to a small offset introduced in absolute pressure')
      print('to enforce monotonicity of the vertical dimension (also on land)') 
 

    print('\n====================================================================================')
    print('\n   The following data set has been written to the new restart file:               \n')
    print(new_restart)
    print('\n------------------------------------------------------------------------------------')
    print('\n     Note: currently, the diagnostic variables time step is undefined                 \
           \n           - this will affect the first time step of iHAMOCC,                         \
           \n           while the effects should be minor compared to uncertainties due to         \
           \n           remapping and interpolations                                             ')
    print('\n     New restart file: '+out_new_restart_file)
    print('\n====================================================================================')
    print('\n   Still to do:')
    print('\n   - sort out time levels                                                             \
           \n   - need for vertical interpolation due to empty layers instead of nn-filling?')

  # Close ds
  ds_coarse_blom.close()
  ds_coarse_hamocc.close()
  ds_fine_blom.close()
  grid_coarse.close()
  grid_fine.close()
  

  if AdjustInventory:
    print('\n\n----- Carrying out inventory adjustment ----')
    # Load the new restart file, BLOM old & fine, gridfiles
    adjustable_restart_file = adjustable_restart_file                                              \
                              if adjustable_restart_file is not None else out_new_restart_file
    new_restart      = xr.open_dataset(adjustable_restart_file)
    ds_coarse_blom   = xr.open_dataset(coarse_blom)
    ds_coarse_hamocc = xr.open_dataset(coarse_hamocc)
    grid_coarse      = xr.open_dataset(coarse_gridfile)
    ds_fine_blom     = xr.open_dataset(fine_blom)
    grid_fine        = xr.open_dataset(fine_gridfile)

    # Initialize empty (temporary) datasets
    adjusted_restart = xr.Dataset()
    tmp              = xr.Dataset()
    for tstep in range(0,2): 
        # define the slicing in the vertical dimension:
        coarse_sstart = tstep     * N_coarse_wclayers
        coarse_send   = (tstep+1) * N_coarse_wclayers
        fine_sstart   = tstep     * N_fine_wclayers
        fine_send     = (tstep+1) * N_fine_wclayers
        try:  
          coarse_sstart_sed = tstep     * N_coarse_sedlayers
          coarse_send_sed   = (tstep+1) * N_coarse_sedlayers
        except: # No sediment in restart file
          pass

        # Convert dp into meter and calculate grid cell volume
        bgcrho_coarse, pa_coarse = bgc_rho(ds_coarse_blom[['dp','temp','saln']].isel(kk2=slice(coarse_sstart,coarse_send)),units=coarse_blom_units)
        dzcoarse                 = dp2dz(pa_coarse,units=coarse_blom_units)
        Vcoarse                  = gridvolume(grid_coarse,dzcoarse)

        bgcrho_fine, pa_fine     = bgc_rho(ds_fine_blom[['dp','temp','saln']].isel(kk2=slice(fine_sstart,fine_send)),units=fine_blom_units)
        dzfine                   = dp2dz(pa_fine,units=fine_blom_units)
        Vfine                    = gridvolume(grid_fine,dzfine)

        # Start adjusting the water column variables
        lkeys = list(new_restart.keys())
        for v in lkeys:
          if 'depth2' in new_restart[v].dims:
            # BGC water column data with two time levels
            var_adjusted,fct = adjust_wc_inventory(ds_coarse_hamocc[v].isel(depth2=slice(coarse_sstart,coarse_send)).to_dataset(),\
                                            new_restart[v].isel(depth2=slice(fine_sstart,fine_send)).to_dataset(),     \
                                            bgcrho_coarse,bgcrho_fine,Vcoarse,Vfine,
                                            ccoords=np.linspace(coarse_sstart,coarse_send-1,N_coarse_wclayers),\
                                            fcoords=np.linspace(fine_sstart,fine_send-1,N_fine_wclayers),vdim_phy='kk2',vdim_bgc='depth2',\
                                            bgc_londim=coarse_hamocc_londim,bgc_latdim=coarse_hamocc_latdim,\
                                            phy_londim=coarse_blom_londim,phy_latdim=coarse_blom_latdim)
            var_adjusted[v].attrs['units']         = ds_coarse_hamocc[v].attrs['units']
            var_adjusted[v].attrs['long_name']     = ds_coarse_hamocc[v].attrs['long_name']
            var_adjusted[v].attrs['missing_value'] = 99999. # somehow, missing value is lost in ds - reset it manually here
 
            adjusted_restart,tmp = store_var(adjusted_restart,tmp,var_adjusted,v,'depth2')
          elif 'depth' in new_restart[v].dims:
            print('!!!! depth: Requires time step check !!!!!')
            if diag_tstep == tstep: 
              var_adjusted,fct = adjust_wc_inventory(ds_coarse_hamocc[v].to_dataset(),\
                                            new_restart[v].to_dataset(),     \
                                            bgcrho_coarse,bgcrho_fine,Vcoarse,Vfine, \
                                            ccoords=np.linspace(0,N_coarse_wclayers-1,N_coarse_wclayers),\
                                            fcoords=np.linspace(0,N_fine_wclayers-1,N_fine_wclayers),vdim_phy='kk',vdim_bgc='depth',\
                                            bgc_londim=coarse_hamocc_londim,bgc_latdim=coarse_hamocc_latdim,\
                                            phy_londim=coarse_blom_londim,phy_latdim=coarse_blom_latdim)
              var_adjusted[v].attrs['units']         = ds_coarse_hamocc[v].attrs['units']
              var_adjusted[v].attrs['long_name']     = ds_coarse_hamocc[v].attrs['long_name']
              var_adjusted[v].attrs['missing_value'] = 99999. # somehow, missing value is lost in ds - reset it manually here
 
              adjusted_restart = xr.merge([adjusted_restart,var_adjusted])

    print('For now, do not reset inventories for sediment - just take it as it is')
    for v in lkeys:
      if v not in list(adjusted_restart.keys()):
        print('We are not adjusting sediment variable: '+v)
        adjusted_restart = xr.merge([adjusted_restart,new_restart[v].to_dataset()])
        adjusted_restart[v].attrs['missing_value'] = 99999.
        adjusted_restart[v].encoding = {} # For some reasons needed - clean-up of encoding to not run
                                          # into conflicting _FillValue versus missing_value issue
                                          # when eventually writing to netcdf       

    # Drop some coordinate variables, if present
    adjusted_restart = drop_depth(adjusted_restart,depth='depth')
    adjusted_restart = drop_depth(adjusted_restart,depth='depth2')
   
    # Add global attributes to adjusted restart file 
    adjusted_restart.attrs = new_restart.attrs
    adjusted_restart.attrs['history'] = 'Adjusted inventory of '+adjustable_restart_file +'; '     \
                                   +' using HAMOCC coarse restart file: '+coarse_hamocc+'; '       \
                                   + 'coarse pressure level file used: '+coarse_blom+'; '          \
                                   + ' pressure levels using 3Dinterp_restart_remapping.py; '      \
                                   + 'coarse gridfile used: '+coarse_gridfile+'; '                 \
                                   + 'fine gridfile used:'+fine_gridfile+'; '                      \
                                   + 'fine pressure level restart file: '+fine_blom+'; '           \
                                   + 'Original file history: '+new_restart.attrs['history']



    # Write out the adjusted restart file
    fil,ext=os.path.splitext(adjustable_restart_file)
    out_file_adjusted = os.path.join(fil+'_adjusted'+ext)
    adjusted_restart.to_netcdf(out_file_adjusted,format="NETCDF3_64BIT")    
    print('\n====================================================================================')
    print('\n   The following data set has been written to the new adjusted restart file:      \n')
    print(adjusted_restart)
    print('\n------------------------------------------------------------------------------------')
    print('\n   New adjusted restart file written: '+out_file_adjusted+'\n\n')
    print('To be done:                                                                             \
           \n     - re-check, document and potentially generalize further                          \
          ')
    

# END OF MAIN
#===================================================================================================
if __name__=='__main__':
    main()
