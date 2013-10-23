      subroutine micom2hamocc(nnb,nb)
!***********************************************************************
!
!**** *SUBROUTINE micom2hammoc* - Interface between MICOM and HAMOCC.
!
!     K. Assmann        *GFI, UiB        initial version
!     J. Schwinger      *GFI, UiB        2013-04-22
!      - 
!
!  
!     Purpose
!     -------
!     - 
!
!**********************************************************************
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_clim.h"
#include "common_ice.h"
#include "common_bgc.h"
c
      integer       :: nnb,nb

      integer       :: i,j,k,l,kn
      real          :: p1,p2,ldp,th,s,pa
      real          :: rho,p_alpha
      real          :: rp(idm,jdm,kdm+1)
c
c$OMP PARALLEL DO
c
c --- calculate pressure at interfaces (necesarry since p has  
c --- not been calculated at restart)
      rp(:,:,:) = 0.0
      do j=1,jj
      do k=1,kk
      kn=k+nnb
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        rp(i,j,k+1)=rp(i,j,k)+dp(i,j,kn)
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c
c
c --- ------------------------------------------------------------------
c --- 2D fields
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do j=1,jj
      do i=1,ii
c
c --- - dimension of grid box in meters
        bgc_dx(i,j)=scpx(i,j)/1.e2
        bgc_dy(i,j)=scpy(i,j)/1.e2
c
      enddo
      enddo
c$OMP END PARALLEL DO
c
c
c$OMP PARALLEL DO
      do j=1,jj
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- - short wave radiation
        bgc_swr (i,j)=swa(i,j)
c --- - sea ice fraction
        bgc_fice(i,j)=ficem(i,j)
c --- - absolute wind speed
        bgc_awnd(i,j)=abswnd(i,j)
c --- - atmospheric co2 concentration
        bgc_atmco2(i,j)=atmco2(i,j)
c --- - sea level pressure
        bgc_slp (i,j)=slp(i,j)
c --- - set depth of surface interface to zero
        bgc_pw(i,j,1)=0.
c --- - depth of layer kk+1 centre
        bgc_pu(i,j,kk+1)=9000.
c
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- 3D fields
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do j=1,jj
      do k=1,kk
      kn=k+nnb
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- - integrated specific volume
        th  = temp(i,j,kn)
        s   = saln(i,j,kn)
        p1  = rp(i,j,k)
        if(dp(i,j,kn) == 0.0) then
          ldp = 1.0
          pa  = ldp/rho(p1,th,s)
        else if(dp(i,j,kn) <  1.0e-2) then
          ldp = dp(i,j,kn)
          pa  = ldp/rho(p1,th,s)
        else
          ldp = dp(i,j,kn)
          p2  = p1+ldp
          pa  = p_alpha(p1,p2,th,s)
        end if
c
c --- - temp and saln 
        bgc_t(i,j,k)=th
        bgc_s(i,j,k)=s
c
c --- - density in g/cm^3
        bgc_rho(i,j,k)=ldp/pa
c
c --- - layer thickness in meters 
        bgc_dp(i,j,k)=0.0
        if(dp(i,j,kn).ne.0.0) bgc_dp(i,j,k)=pa/onem
c     
c --- - inverse layer thickness in meters 
        bgc_dpio(i,j,k)=-99.999999999 
        if(dp(i,j,kn).ne.0.0) bgc_dpio(i,j,k)=onem/pa
c
c --- - depth of layer interfaces
        bgc_pw(i,j,k+1)=bgc_pw(i,j,k)+bgc_dp(i,j,k)
c
c --- - depth of layer centres
        bgc_pu(i,j,k)=bgc_pw(i,j,k)+0.5*bgc_dp(i,j,k)
c
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c
      ldtday  =mod(nstep,nstep_in_day)/nphys
      ldtmonth=((nday-1)*nstep_in_day/nphys)+ldtday
      kpndtrun=(nstep2-nstep1)/nphys
c
      return
      end
