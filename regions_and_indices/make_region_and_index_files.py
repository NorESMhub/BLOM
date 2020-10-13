import numpy as np
import xarray as xr

xdim=208
ydim=512
zdim=53
#
mask=np.ones((ydim,xdim))
mask[ydim//2:,:] = 2
mask[0,:]=-1
mask[-1,:]=-1
#
region = xr.DataArray(mask.astype(int),dims=('y','x'),name='region')
region_names = xr.DataArray(['south','north'],dims=('regions'),name='region_names')
out = xr.merge([region,region_names])
#
out.to_netcdf('ocean_regions.nc',encoding={'region_names':{'char_dim_name':'strlen','_FillValue':None},'region':{'_FillValue': -1}},format='NETCDF3_CLASSIC')
#
# write user defined sections for transport calculations
f=open('section_index.dat','w')

f.write('Name: middle_section\n')
for j in range(1,ydim-1):
    f.write(str(xdim//2+1).rjust(3)+' '+str(j+1).rjust(3)+' 1  0\n')

f.close()

# write sections for meridional heat transports
f1=open('mertra_index.dat','w')
for j in range(1,ydim-1):
    f1.write('Section '+(str(j+1)+'.000').rjust(7)+' '+str(xdim)+'\n')
    for i in range(xdim):
        f1.write(str(i+1).rjust(3)+' '+str(j+1).rjust(3)+' 0  1\n')

f1.close()

# create a dummy inicon.nc - only the sigma coordinate is important
#
sigma=np.ones(zdim)*34.59 #
for k in range(3,zdim+1): 
    sigma[k-1]=sigma[k-2] + 0.05*k*(1.+np.cos(np.pi*(k-1.)/zdim))/zdim

#sigma2 = np.array([27.22 , 27.72 , 28.202, 28.681, 29.158, 29.632, 30.102, 30.567, 31.026,
#       31.477, 31.92 , 32.352, 32.772, 33.176, 33.564, 33.932, 34.279, 34.602,
#       34.9  , 35.172, 35.417, 35.637, 35.832, 36.003, 36.153, 36.284, 36.398,
#       36.497, 36.584, 36.66 , 36.728, 36.789, 36.843, 36.893, 36.939, 36.982,
#       37.022, 37.06 , 37.096, 37.131, 37.166, 37.199, 37.231, 37.264, 37.295,
#       37.327, 37.358, 37.388, 37.419, 37.45 , 37.48 , 37.58 , 37.8])
#
sigma3d = np.swapaxes(np.tile(sigma,(ydim,xdim,1)).T,1,2)
#zdim = len(sigma)
dum = np.zeros((zdim,ydim,xdim))
dz_dum = xr.DataArray(dum,dims={'z':zdim,'y':ydim,'x':xdim},name='dz')
T_dum  = xr.DataArray(dum,dims={'z':zdim,'y':ydim,'x':xdim},name='temp')
S_dum  = xr.DataArray(dum,dims={'z':zdim,'y':ydim,'x':xdim},name='saln')
sigma_dum = xr.DataArray(sigma3d,dims={'z':zdim,'y':ydim,'x':xdim},name='sigma')
#
ini_out = xr.merge([sigma_dum,T_dum,S_dum,dz_dum])
ini_out.to_netcdf('inicon.nc',format='NETCDF3_CLASSIC',encoding={'dz':{'_FillValue':None},'temp':{'_FillValue':None},'saln':{'_FillValue':None},'sigma':{'_FillValue':None}})
