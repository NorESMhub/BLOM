      subroutine setlogunit

      ! Uses modules

      use mod_xc
      use shr_file_mod, only: shr_file_getUnit, shr_file_setIO

      implicit none

      logical exists

      ! ----------------------------------------------------------------
      ! Redirect standard out to a log file if requested
      ! ----------------------------------------------------------------
      if (mnproc.eq.1) then
         inquire(file = 'ocn_modelio.nml', exist = exists)
         if (exists) then
            lp = shr_file_getUnit()
            call shr_file_setIO('ocn_modelio.nml', lp)
         endif
      endif

      end subroutine setlogunit
