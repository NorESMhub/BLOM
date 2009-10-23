      subroutine micom2hamocc(nnb,nb)
c
c --- version 2.8 (NERSC version)
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
      integer i,j,k,l,kn
      integer nnb,nb
c
      real intpol
      external intpol
c
c$OMP PARALLEL DO
      do j=1,jj
      do i=1,ii
c
cka dimension of grid box in meters
        bgc_dx(i,j)=scpx(i,j)/1.e2
        bgc_dy(i,j)=scpy(i,j)/1.e2
c
      enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
cka sea surface height for tracer inventory
c        bgc_slvl(i,j)=sealv(i,j)/100.   
cka short wave radiation
        bgc_swr (i,j)=swa(i,j)
cka sea ice fraction
        bgc_fice(i,j)=ficem(i,j)
cka absolute wind speed
        bgc_awnd(i,j)=abswnd(i,j)
cka atmospheric co2 concentration
        bgc_atmco2(i,j)=atmco2(i,j)
cka sea level pressure
c        bgc_slp (i,j)=slp(i,j)
cka snow thickness
c        bgc_hsnw(i,j)=hsnwm(i,j)
cka ice thickness
c        bgc_hice(i,j)=hicem(i,j)
cka set surface interface to zero
        bgc_pw(i,j,1)=0.
cka depth of layer kk+1 centre
        bgc_pu(i,j,kk+1)=9000.
c
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
      do k=1,kk
      kn=k+nnb
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
cka layer thickness in meters 
        bgc_dp(i,j,k)=dp(i,j,kn)/(onem*1.027)
        bgc_dpio(i,j,k)=-99.999999999 
        if(bgc_dp(i,j,k).ne.0.0) bgc_dpio(i,j,k)=1./bgc_dp(i,j,k)
c
cka temp and saln 
        bgc_t(i,j,k)=temp(i,j,kn)
        bgc_s(i,j,k)=saln(i,j,kn)
c
cka depth of layer interfaces
        bgc_pw(i,j,k+1)=bgc_pw(i,j,k)+(dp(i,j,kn)/(onem*1.027))
c
cka depth of layer centres
        bgc_pu(i,j,k)=bgc_pw(i,j,k)+0.5*(dp(i,j,kn)/(onem*1.027))
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
c

