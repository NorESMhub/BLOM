      block data frcdat
c
      use mod_xc
c
      implicit none
c
#include "common_forc.h"
c
      data
     . albw_d   /.065/    ! --                 albedo over water; diff light
     .,rhowat   /1024./   ! kg / m^3           app. density of surf. sea water
     .,sref     /34.65/   ! g / kg             global ref. surface salinity
     .,nrfets   /7/       ! days               e-folding time scale for
                          !                    detrainment of runoff reservoar
c
      end
