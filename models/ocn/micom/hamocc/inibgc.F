      subroutine inibgc
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_bgc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_geo.h"
c
      integer i,j,k,l,nt
c  set bgc timestep to nphys*dt   
      bgcdt=nphys*baclin

      do j=1,jj
      do i=1,ii
c initialise 2D arrays for HAMOCC interface
      bgc_dx(i,j)=0.0
      bgc_dy(i,j)=0.0
      bgc_slp(i,j)=0.0
      omask (i,j)=0.0
      pglon (i,j)=plon(i,j)
      pglat (i,j)=plat(i,j)
c
c initialise 3D arrays for HAMOCC interface
      do k=1,kk
c
      do nt=itrbgc,itrbgc+ntrbgc-1
         trc(i,j,k,nt)=0.0
      enddo

      bgc_dp  (i,j,k)=0.0
      bgc_dpio(i,j,k)=0.0
      bgc_rho (i,j,k)=0.0
      bgc_t   (i,j,k)=0.0
      bgc_s   (i,j,k)=0.0
      bgc_pu  (i,j,k)=0.0
      bgc_pw  (i,j,k)=0.0
      enddo
      bgc_pu  (i,j,kk+1)=0.0
      bgc_pw  (i,j,kk+1)=0.0
      enddo
      enddo
c
      do 10 j=1,jj
      do 10 l=1,isp(j)
      do 10 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
      omask(i,j)=1.
 10   continue

      return
      end
