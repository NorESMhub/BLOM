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
c --- 'veldff' = diffusion velocity (cm/s) for momentum dissipation
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
      data veldff/2.0/,viscos/0.5/
c
c --- 'thkmin' = minimum mixed-layer thickness (m)
c --- 'acurcy' = permissible roundoff error in column integral calc.
      data thkmin/10./,acurcy/1.e-10/
c
c --- slip=+1  for free-slip boundary cond., slip=-1  for non-slip cond.
c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
c --- 'cb'     = coefficient of quadratic bottom friction
c --- 'thkbot' = thickness of bottom boundary layer (m)
      data slip/-1./,cbar/5./,cb/.002/,thkbot/10./
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
c --- 'rearth' = radius of the earth (cm)
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
      data g/980.6/,rearth/6.37122e8/,spcifh/3.99/,rhoa_r/1.225e-3/,
     .     cd_r/1.2e-3/,ch_r/1.1e-3/,ce_r/1.3e-3/,alpha0/1./,
     .     epsil/1.e-11/,wg2_r/1.e3/,raddep/500./,redfac/.58/,
     .     betabl/23./
c
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
      data itest,jtest/192,192/
c
c ---  s w i t c h e s    (if set to .true., then...)
c
      data csdiag/.false./,cnsvdi/.true./
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
