      subroutine inibgc
c
c --- version 2.8 (NERSC version)
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_bgc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_forc.h"
#include "common_clim.h"
#include "common_ice.h"
c
      integer i,j,k,l,nt
cka  set bgc timestep to nphys*dt   
      bgcdt=nphys*baclin

      do j=1,jj
      do i=1,ii
cka initialise 2D arrays for HAMOCC interface
      bgc_dx  (i,j)=0.0
      bgc_dy  (i,j)=0.0
      omask(i,j)=0.
c        
      do k=1,kk
c
      do nt=itrbgc,itrbgc+ntrbgc-1
         trc(i,j,k,nt)=0.
      enddo
c
cka initialise 3D arrays for HAMOCC interface
      bgc_dp  (i,j,k)=0.
      bgc_dpio(i,j,k)=0.
      bgc_t   (i,j,k)=0.
      bgc_s   (i,j,k)=0.
      bgc_pu  (i,j,k)=0.
      bgc_pw  (i,j,k)=0.
      enddo
      bgc_pu  (i,j,kk+1)=0.
      bgc_pw  (i,j,kk+1)=0.
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
c
c
c> Revision history:
c>
c> May  1997 - eliminated "+ thbase" from definition of thmix (loop 10)
c> Mar. 1998 - added -th3d-
c> Aug. 2000 - added diagnostic count of static instabilities
