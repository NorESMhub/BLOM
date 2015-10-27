      subroutine icaux_ben02
c
c --- ------------------------------------------------------------------
c --- Define initial conditions specifically when using the Bentsen and
c --- Drange (2002) forcing and thermodynamic sea ice
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02, only: tsi_tda,tml_tda,sml_tda,alb_tda,fice_tda,ntda
      use mod_thdysi, only: tsrfm,ticem,albs_f,fice_max
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_seaice.h"
c
      integer i,j,l
c
c --- ------------------------------------------------------------------
c --- Define initial conditions for thermodynamic sea ice model and
c --- accumulation variables
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (ricclm(i,j,12).lt..05) then
            hicem(i,j)=0.
            ficem(i,j)=0.
            hsnwm(i,j)=0.
            tsrfm(i,j)=temp(i,j,1)+t0deg
            tsi_tda(i,j)=temp(i,j,1)+t0deg
          else
            hicem(i,j)=2.
            ficem(i,j)=min(fice_max,ricclm(i,j,12)) ! assuming start in January
            hsnwm(i,j)=0.1
            tsrfm(i,j)=sstclm(i,j,1)
            tsi_tda(i,j)=(sstclm(i,j,1)
     .                   -(1.-ficem(i,j))*(temp(i,j,1)+t0deg))
     .                   /ficem(i,j)
            iagem(i,j)=0.
          endif
          ticem(i,j)=tsrfm(i,j)
          tml_tda(i,j)=temp(i,j,1)+t0deg
          sml_tda(i,j)=saln(i,j,1)
          alb_tda(i,j)=albs_f*ficem(i,j)+albw_d*(1.-ficem(i,j))
          fice_tda(i,j)=ficem(i,j)
        enddo
        enddo
      enddo
c
      ntda=1
c
      return
      end
