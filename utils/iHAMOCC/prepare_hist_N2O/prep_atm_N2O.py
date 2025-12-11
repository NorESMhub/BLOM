import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
"""
Simple script to read in the N2O concentration file, 
calculate a mean (of latitudinally non-varying) latitudinal data
and print the values to a text-csv file - and append already the working precision 
(apart from the last value, where this needs to be added manually) 
"""
n2ovar='N2O_LBC'

fil = '/cluster/shared/noresm/inputdata/atm/waccm/lb/LBC_1750-2015_CMIP6_GlobAnnAvg_c180926.nc'

ds = xr.open_dataset(fil)

molpmol2ppb = 1e12

# Print the var for checking values
print(ds[n2ovar].isel(time=1).values)
print(ds[n2ovar].mean(dim='lat').values*molpmol2ppb)

# Length of the data set and start and end time
print(len(ds.time))
print(ds.time.min().values)
print(ds.time.max().values)

# Save the values in one row to csv file, while appending _rp for Fortran working precision
np.savetxt('historical_N2O.csv', (ds[n2ovar].mean(dim='lat').values*molpmol2ppb).reshape(1, -1)   , delimiter='_rp, ',fmt='%6.5f')

# Quick plot to investigate the data
if 1:
  plt.figure()
  (molpmol2ppb*ds[n2ovar]).plot.pcolormesh(rasterized=True,add_colorbar=True)
  plt.show()
