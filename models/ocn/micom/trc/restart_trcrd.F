      subroutine restart_trcrd(nday1,rstfnm_ext)
c
c --- ------------------------------------------------------------------
c --- Read tracer initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ext*(*)
      integer nday1
c
#ifdef HAMOCC
      call restart_hamoccrd(nday1,rstfnm_ext)
#endif
c
      return
      end
