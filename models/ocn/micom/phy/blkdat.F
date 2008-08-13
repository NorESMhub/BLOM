      block data blkdat
c
c --- NERSC version
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
c
c --- layer reference potential densities (sigma units):
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- s i g m a _ 2
      data sigmar/24.000,25.000,26.000,27.000,27.980,
     .            28.880,29.680,30.400,31.070,31.720,
     .            32.370,33.020,33.670,34.320,34.940,
     .            35.500,35.954,36.285,36.530,36.712,
     .            36.846,36.946,37.020,37.074,37.115,
     .            37.145,37.175,37.205,37.235,37.265,
     .            37.295,37.325,37.355,37.385,37.415,
     .            37.440,37.460,37.480,37.580,37.800/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- 'thkdff' = diffusion velocity (cm/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (cm/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (cm/s) for temp/salin. mixing
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
      data thkdff/1.0/,veldff/1.0/,temdff/0.35/,viscos/0.5/
c
c --- 'diapyc' = diapycnal diffusivity times buoyancy freq. (cm**2/s**2)
c --- 'vertmx' = diffusion velocity (cm/s) for mom.mixing across mix.layr.base
c --- 'thkmin' = minimum mixed-layer thickness (m)
c --- 'acurcy' = permissible roundoff error in column integral calc.
      data diapyc/1.8e-3/,vertmx/2.e-3/
      data thkmin/20./,acurcy/1.e-10/
c
c --- slip=+1  for free-slip boundary cond., slip=-1  for non-slip cond.
c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
c --- 'thkbot' = thickness of bottom boundary layer (m)
c --- 'sigjmp' = minimum density jump at mixed-layer bottom (theta units)
      data slip/-1./,cbar/10./,thkbot/10./,sigjmp/2.e-5/
c
c --- weights for time smoothing
ccc      data wuv1,wuv2/.5,.25/
      data wuv1,wuv2/.75,.125/
      data wts1,wts2/.875,.0625/
ccc      data wts1,wts2/.9375,.03125/
ccc      data wts1,wts2/.96875,.015625/
      data wbaro/.125/
ccc      data wbaro/.0625/
      data wpgf/.25/
c
c --- layer thicknesses in units of pressure (g/cm/sec^2):
      data tenm,onem,tencm,onecm,onemm/980600.,98060.,9806.,980.6,98.06/
      data radian/57.295779513/,pi/3.1415926536/
c
c --- 'g'      = gravitational acceleration
c --- 'spcifh' = specific heat of sea water (j/g/deg)
c --- 'rhoa_r' = reference air density (g/cm**3)
c --- 'cd_r'   = reference transfer coefficient of momentum
c --- 'ch_r'   = reference transfer coefficient of sensible heat
c --- 'ce_r'   = reference transfer coefficient of tent heat
c --- 'wg2_r'  = reference gustiness squared (cm**2/s**2)
c --- 'alpha0' = reference value of specific volume (cm**3/g)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c --- 'raddep' = maximum depth of light penetration (m)
c --- 'redfac' = red fraction of light aborbed in mixed layer (jerlov 1)
c --- 'betabl' = blue light extinction coefficient (m) (jerlov 1)
      data g/980.6/,spcifh/3.99/,rhoa_r/1.225e-3/,cd_r/1.2e-3/,
     .     ch_r/1.1e-3/,ce_r/1.3e-3/,alpha0/1./,epsil/1.e-11/,
     .     wg2_r/1.e3/,raddep/500./,redfac/.58/,betabl/23./
c
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
      data itest,jtest/192,192/
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- thermo      use thermodynamic forcing functions
c
      data thermo/.true./,csdiag/.false./,cnsvdi/.false./
c
c --- use 'huge' to initialize array portions that the code should never access
      data huge/1.e33/
c
c --- i/o file names
c
      data path /'./ '/
      data path1/'./ '/
      data path2/'./ '/
c
      end
c
c
c> Revision history:
c>
c> Mar. 1999 - made eqn-of-state coefficients conform to Brydon & Sun
c> Oct. 2000 - changed units of -thkbot- from pressure units to meters
