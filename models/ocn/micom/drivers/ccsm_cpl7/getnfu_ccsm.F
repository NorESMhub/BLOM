      subroutine getnfu_ccsm
c
c --- ------------------------------------------------------------------
c --- Get file unit from CCSM system
c --- ------------------------------------------------------------------
c
      use mod_xc
      use shr_file_mod, only: shr_file_getUnit
c
      implicit none
c
      nfu=shr_file_getUnit()
c
      return
      end
