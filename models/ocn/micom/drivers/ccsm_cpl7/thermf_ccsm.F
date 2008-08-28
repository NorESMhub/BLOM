      subroutine thermf_ccsm(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. To be used when coupled to CCSM
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_geo.h"
#include "common_ice.h"
c
      integer i,j,l
      real fwflx,totsfl,totwfl,sflxc
c
      real swtfrz
      external swtfrz
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- ------------------------------------------------------------------
c --- --- Fresh water and salt fluxes
c --- ------------------------------------------------------------------
c
c --- --- Copy runoff flux to rnfflx for consistency with uncoupled
c --- --- configuration
          rnfflx(i,j)=rnf(i,j)
c
c --- --- Fresh water flux to the ocean [m/s]
          fwflx=(eva(i,j)+lip(i,j)+sop(i,j)+rnfflx(i,j))*1.e-3
c
c --- --- Virtual salt flux from ocean [kg m-2 s-1] (positive upwards)
          salflx(i,j)=saln(i,j,k1n)*fwflx
c
c --- --- Accumulate area weighted virtual salt flux and fresh water
c --- --- flux
          util1(i,j)=salflx(i,j)*scp2(i,j)
          util2(i,j)=fwflx*scp2(i,j)
c
c --- ------------------------------------------------------------------
c --- --- Heat fluxes
c --- ------------------------------------------------------------------
c
c --- --- Freezing/melting potential [W m-2]. A positive flux means the
c --- --- ocean surface has a temperature below freezing temperature and
c --- --- must be heated
          frzmlt(i,j)=(swtfrz(saln(i,j,k1n))-temp(i,j,k1n))
     .                *spcifh*dp(i,j,k1n)/(g*delt1)*1.e4
c
c --- --- Total heat flux from the ocean [W cm-2] (positive upwards).
c --- --- Apply the freezing/melting potential if it cause a heating of
c --- --- the ocean, store the energy used as frazil ice production
          surflx(i,j)=-(swa(i,j)+nsf(i,j)+max(0.,frzmlt(i,j)))*1.e-4
c
c --- --- short wave heat flux from the ocean [W cm-2] (positive upwards)
          sswflx(i,j)=-swa(i,j)*1.e-4
c
c --- ------------------------------------------------------------------
c --- --- Relaxation fluxes
c --- ------------------------------------------------------------------
c
c --- --- set to zero!
          surrlx(i,j)=0.
          salrlx(i,j)=0.
c
c --- -------------------------------------------------------------------
c --- --- Friction velocity (cm/s)
c --- -------------------------------------------------------------------
c
          ustar(i,j)=ustarw(i,j)*1.e2
c
        enddo
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- Compute correction to the virtual salt flux so it is globally
c --- consistent with a salt flux based on some reference salinity and
c --- add salt flux recieved from coupler. Also convert salt flux to
c --- MICOM units [10e-3 g cm-2 s-1]
c --- ------------------------------------------------------------------
c
      call xcsum(totsfl,util1,ip)
      call xcsum(totwfl,util2,ip)
c
c --- Correction for salt flux
      sflxc=(srefg*totwfl-totsfl)/area
      if (mnproc.eq.1) then
        write (lp,*) 'thermf: totsfl/area,sflxc',totsfl/area,sflxc
      endif
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          salflx(i,j)=(salflx(i,j)+sflxc-sfl(i,j))*1.e2/alpha0
        enddo
        enddo
      enddo
c
      return
      end
