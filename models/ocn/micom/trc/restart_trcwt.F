      subroutine restart_trcwt(rstfnm_ext)
c
c --- ------------------------------------------------------------------
c --- Write tracer state to restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ext*(*)
c
#ifdef HAMOCC
      call restart_hamoccwt(rstfnm_ext)
#endif
c
      return
      end
