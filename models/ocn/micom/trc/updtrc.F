      subroutine updtrc(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- update tracers due to non-passive processes
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
c
      integer m,n,mm,nn,k1m,k1n
c
#ifdef HAMOCC
      call hamocc_step(m,n,mm,nn,k1m,k1n)
#endif
#ifdef IDLAGE
      call idlage_step
#endif
c
      return
      end
