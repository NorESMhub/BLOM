      subroutine restart_trcrd(nday1,rstfnm_ocn)
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
      integer nday1
c
      call restart_ocntrcrd(rstfnm_ocn)
c
#ifdef HAMOCC
      call restart_hamoccrd(nday1,rstfnm_ocn)
#endif
c
      return
      end
