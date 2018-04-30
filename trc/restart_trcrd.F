      subroutine restart_trcrd(rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Read tracer state from restart files
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ocn*(*)
c
      call restart_ocntrcrd(rstfnm_ocn)
c
#ifdef HAMOCC
      call restart_hamoccrd(rstfnm_ocn)
#endif
c
      return
      end
