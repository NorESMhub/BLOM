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
     .,t0deg    /273.15/  ! K                  zero deg celcius in K
     .,nrfets   /7/       ! days               e-folding time scale for
                          !                    detrainment of runoff reservoar
c
      end
