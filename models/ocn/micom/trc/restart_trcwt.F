      subroutine restart_trcwt(rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Write tracer state to restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ocn*(*)
c
#ifdef HAMOCC
      call restart_hamoccwt(rstfnm_ocn)
#endif
c
      return
      end
