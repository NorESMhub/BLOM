      subroutine initrc
c
c --- ------------------------------------------------------------------
c --- initialization of tracers
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#ifdef HAMOCC
      call hamocc_init
#endif
#ifdef IDLAGE
      call idlage_init
#endif
c
      return
      end
