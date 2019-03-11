      subroutine micom2hamocc(m,n,mm,nn)
c***********************************************************************
c
c**** *SUBROUTINE micom2hammoc* - Interface between MICOM and HAMOCC.
c
c     K. Assmann        *GFI, UiB        initial version
c     J. Schwinger      *GFI, UiB        2013-04-22
c      - 
c
c     Modified
c     --------
c     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
c     - removed inverse of layer thickness
c     - added sediment bypass preprocessor option
c  
c     Purpose
c     -------
c     - 
c
c**********************************************************************
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra,natm
      use mo_carbch,     only: ocetra,atm,atm2
      use mo_sedmnt,     only: sedlay,powtra,sedhpl,burial
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_seaice.h"
#include "common_bgc.h"
#ifndef sedbypass
#include "common_bgcs.h"
#endif
c
      integer       :: m,n,mm,nn

      integer       :: i,j,k,l,nns,kn
      real          :: p1,p2,ldp,th,s,pa
      real          :: rho,p_alpha
      real          :: rp(idm,jdm,kdm+1)
c
      nns=(n-1)*ks
c
c$OMP PARALLEL DO PRIVATE(kn)
c
c --- calculate pressure at interfaces (necesarry since p has  
c --- not been calculated at restart)
      rp(:,:,:) = 0.0
      do j=1,jj
      do k=1,kk
      kn=k+nn
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
c$OMP PARALLEL DO PRIVATE(kn,th,s,p1,p2,ldp,pa)
      do j=1,jj
      do k=1,kk
      kn=k+nn
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
c --- ------------------------------------------------------------------
c --- - return if restart (HAMOCC fields are not allocated yet)
c --- ------------------------------------------------------------------
      if( .not. allocated(ocetra) ) return
c
c --- ------------------------------------------------------------------
c --- pass tracer fields from ocean model; convert mol/kg -> kmol/m^3
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
      do k=1,kk
      kn=k+nn
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        ocetra(i,j,k,:)=trc(i,j,kn,itrbgc:itrbgc+ntrbgc-1)*bgc_rho(i,j,k)
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- pass sediments fields (a two time-level copy of sediment fields
c --- is kept outside HAMOCC)
c --- ------------------------------------------------------------------
c
#ifndef sedbypass
      nns=(n-1)*ks
c
c$OMP PARALLEL DO PRIVATE(kn)
c
      do j=1,jj
      do k=1,ks
      kn=k+nns
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        sedlay(i,j,k,:) = sedlay2(i,j,kn,:)
        powtra(i,j,k,:) = powtra2(i,j,kn,:)
        burial(i,j,:)   = burial2(i,j,n,:)
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
#endif
c
c --- ------------------------------------------------------------------
c --- pass atmosphere fields if required (a two time-level copy of
c --- atmosphere fields is kept outside HAMOCC)
c --- ------------------------------------------------------------------
c
#if defined(BOXATM) || defined(DIFFAT)
c$OMP PARALLEL DO
      do j=1,jj
      do i=1,ii
        atm(i,j,:) = atm2(i,j,n,:)
      enddo
      enddo
c$OMP END PARALLEL DO
#endif
c
c --- ------------------------------------------------------------------
c --- set control parameters
c --- ------------------------------------------------------------------
c
      ldtday   = mod(nstep,nstep_in_day)
      ldtmonth = (nday-1)*nstep_in_day + ldtday
      kpndtrun = nstep2-nstep1
c
      return
      end
