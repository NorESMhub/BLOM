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
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
      real, dimension(ntr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ttrsf,ttrav
      real tottrsf,tottrav,trflxc
#endif
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
c --- --- Freezing/melting potential [J m-2]. A positive flux means the ocean
c --- --- surface has a temperature below freezing temperature and must
c --- --- be heated. Note the freezing potential is multiplied by 1/2
c --- --- due to the leap-frog time stepping. The melting potential uses
c --- --- time averaged quantities since it is not accumulated.
          frzpot(i,j)=max(0.,swtfrz(saln(i,j,k1n))-temp(i,j,k1n))
     .                *spcifh*dp(i,j,k1n)/(2.*g)*1.e4
          mltpot(i,j)=min(0.,swtfrz(.5*(saln(i,j,k1m)+saln(i,j,k1n)))
     .                      -.5*(temp(i,j,k1m)+temp(i,j,k1n)))
     .                *spcifh*.5*(dp(i,j,k1m)+dp(i,j,k1n))/g*1.e4
c
c --- --- Total heat flux from the ocean [W cm-2] (positive upwards).
c --- --- Apply the freezing potential
          surflx(i,j)=-(swa(i,j)+nsf(i,j)+frzpot(i,j)/baclin)*1.e-4
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
#ifdef TRC
c --- ------------------------------------------------------------------
c --- --- Virtual tracer fluxes
c --- ------------------------------------------------------------------
c
          do nt=1,ntrc
            trflx(nt,i,j)=trc(i,j,1,nt)*fwflx
            ttrsf(nt,i,j)=trflx(nt,i,j)*scp2(i,j)
            ttrav(nt,i,j)=trc(i,j,1,nt)*scp2(i,j)
          enddo
#endif
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
#ifdef TRC
c --- ------------------------------------------------------------------
c --- Compute correction to the virtual tracer flux
c --- ------------------------------------------------------------------
c
      do nt=1,ntrc
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=ttrsf(nt,i,j)
            util2(i,j)=ttrav(nt,i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call xcsum(tottrsf,util1,ip)
        call xcsum(tottrav,util2,ip)
c
        tottrav=tottrav/area
c
        trflxc=(tottrav*totwfl-tottrsf)/area
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            trflx(nt,i,j)=(trflx(nt,i,j)+trflxc)*1.e2/alpha0
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c 
      enddo
#endif
c
      return
      end
