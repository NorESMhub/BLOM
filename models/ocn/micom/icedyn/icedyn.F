      subroutine icedyn(m,n,mm,nn,k1m,k1n)
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_ice.h"
#include "param_atm.h"
#include "common_atm.h"
#include "common_clndr.h"
c
c --- ------------------------------------------------------------------
c --- Call sea ice dynamic module the first time step of a new day
c --- ------------------------------------------------------------------
c
      if (mod(nstep,nstep_in_day).eq.1) then
        if (mnproc.eq.1) then
          write (lp,*) 'seaice...'
        endif
        call mod2ice(m,n,k1m,k1n,ztx,mty,ub,vb,u,v,sealv,scuxi,scvyi,
     .               hicem,ficem,onecm)
        if (mnproc.eq.1) then
          call seaice
        endif
        call ice2mod(ustari,tauxice,tauyice,uicem,vicem,hicem,ficem,itp)
        call xctilr(uicem, 1,1, 2,2, halo_uv)
        call xctilr(vicem, 1,1, 2,2, halo_vv)
        if (csdiag) then
          call chksummsk(ustari,ip,1,'ustari')
          call chksummsk(tauxice,iu,1,'tauxice')
          call chksummsk(tauyice,iv,1,'tauyice')
          call chksummsk(uicem,iu,1,'uicem')
          call chksummsk(vicem,iv,1,'vicem')
        endif
      endif
c
c --- ------------------------------------------------------------------
c --- Do horizontal advection and diffusion of sea ice scalar fields
c --- ------------------------------------------------------------------
c
      call iadvdf
c
      return
      end
