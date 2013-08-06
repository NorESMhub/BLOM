      subroutine difest(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate layer interface, isopycnal, and diapycnal diffusivities
c --- ------------------------------------------------------------------
#define DIAG
#undef DIAG
c
      use mod_xc
#ifdef DIAG
      use mod_nctools
#endif
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_ice.h"
#include "common_forc.h"
#include "common_geo.h"
c
c --- parameters:
c ---   iidtyp - type of interface and isopycnal diffusivities. If
c ---            iidtyp=1 the diffusivities are diffusive velocities
c ---            multiplied by the local horizontal grid scale, if
c ---            iidtyp=2 the diffusivities are parameterized according
c ---            to Eden and Greatbatch (2008).
c ---   bdmtyp - type of background diapycnal mixing. If bdmtyp=1
c ---            the background diffusivity is a constant divided by the
c ---            Brunt-Vaisala frequency, if bdmtyp=2 the background
c ---            diffusivity is constant.
c ---   bdmldp - If bdmldp=1, make the background mixing latitude
c ---            dependent according to Gregg et al. (2003).
c ---   tdmflg - If tdmflg=1, apply tidally driven diapycnal mixing.
c ---   iwdflg - If iwdflg=1, reduce background diapycnal diffusivity
c ---            due to internal wave damping under sea-ice.
c ---   dpbmin - smallest layer thickness allowed in evaluating
c ---            local gradient richardson number [g/cm/s**2].
c ---   drhomn - minimum density difference in evaluations the
c ---            Brunt-Vaisala frequency and the local gradient
c ---            Richardson number [g/cm*3].
c ---   thkdff - diffusive velocity for thickness diffusion [cm/s].
c ---   temdff - diffusive velocity for tracer isopycnal diffusion
c ---            [cm/s].
c ---   ri0    - critical gradient richardson number for shear driven
c ---            vertical mixing [].
c ---   nus0   - maximum shear driven diapycnal diffusivity
c ---            [cm**2/s].
c ---   nug0   - maximum gravity current diapycnal diffusivity
c ---            [cm**2/s].
c ---   drho0  - critical local interface density difference [g/cm**3]
c ---   nuls0  - maximum diapycnal diffusivity applied when local
c ---            stability is weak [cm**2/s].
c ---   bdmc1  - background diapycnal diffusivity times buoyancy
c ---            frequency [cm**2/s**2].
c ---   bdmc2  - background diapycnal diffusivity [cm**2/s].
c ---   iwdfac - internal wave dissipation factor under sea ice [].
c ---   dmxeff - diapycnal mixing efficiency [].
c ---   tdmq   - tidal dissipation efficiency [].
c ---   tdmls0 - tidal driven mixing length scale below critical
c ---            latitude [g/cm/s**2].
c ---   tdmls1 - tidal driven mixing length scale above critical
c ---            latitude [g/cm/s**2].
c ---   tdclat - critical latitude for tide M2 propagation [].
c ---   tddlat - latitudinal transition zone for different tidal driven
c ---            mixing length scales near the critical latitude.
c ---   cori30 - coriolis parameter at 30N [1/s].
c ---   bvf0   - reference stratification in the parameterization of
c ---            latitude dependent background diapycnal mixing [1/s].
c ---   nubmin - minimum background diapycnal diffusivity [cm**2/s].
c ---   dpgc   - thickness of region near the bottom where the maximum
c ---            diffusivity is increased due to gravity current mixing
c ---            processes [g/cm/s**2].
c ---   dpdiav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate lateral diffusivities in the
c ---            non-isopycnic layers [g/cm/s**2].
c ---   dpddav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate diapycnal diffusivities in the
c ---            non-isopycnic layers [g/cm/s**2].
c ---   dpnbav - thickness of region near the bottom used to estimate
c ---            bottom Brunt-Vaisala frequency [g/cm/s**2].
      integer iidtyp,bdmtyp,bdmldp,tdmflg,iwdflg
      real dptmin,dpbmin,drhomn,thkdff,temdff,
     .     ri0,nus0,nug0,drho0,nuls0,bdmc1,bdmc2,iwdfac,dmxeff,
     .     tdmq,tdmls0,tdmls1,tdclat,tddlat,cori30,bvf0,nubmin,dpgc,
     .     dpdiav,dpddav,dpnbav,ustmin,kappa,bfeps,sleps,zetas,as,cs
      parameter (iidtyp=2,bdmtyp=2,bdmldp=1,tdmflg=1,iwdflg=1,
     .           dptmin=98060.,dpbmin=980.6,drhomn=6.e-6,thkdff=.5,
     .           temdff=.35,ri0=.7,nus0=50.,nug0=2500.,
     .           drho0=6.e-6,nuls0=500.,bdmc1=5.e-4,bdmc2=.1,iwdfac=.06,
     .           dmxeff=.2,tdmq=1./3.,tdmls0=500.*98060.,
     .           tdmls1=100.*98060.,tdclat=74.5,tddlat=3.,
     .           cori30=7.2722e-5,bvf0=5.24e-3,nubmin=.001,
     .           dpgc=300.*98060.,dpdiav=100.*98060.,dpddav=10.*98060.,
     .           dpnbav=250.*98060.,ustmin=.1,kappa=.4,bfeps=1.e-12,
     .           sleps=.1,zetas=-1.,as=-28.86,cs=98.96)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  dv2
      real, dimension(1-nbdy:idm+nbdy,kdm) ::
     .  du2,drho,bvfsq,bvf,rig
      real, dimension(1-nbdy:idm+nbdy) ::
     .  tup,pup,sup,cr,bcrrd,afeql,bvfbot,dps,dfints,dfddsu,dfddsl
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  mdv2
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  kfil,kmax
      integer, dimension(1-nbdy:idm+nbdy,kdm) ::
     .  mdu2
      integer, dimension(1-nbdy:idm+nbdy) ::
     .  kfpl,klpl
      integer i,j,k,l,kn
      real q,plo,tlo,slo,egr,rhisc,els,nus,nub,nut,nuls,tdmf,ust,mols,h,
     .     sg,zeta,phis,ws
c
      real rho,sig
      external rho,sig
c
#ifdef TKE
#  include "param_tke.h"
#  include "common_tke.h"
      real gls_c3,tke_prod,tke_buoy,tke_epsilon,gls_prod,gls_buoy,
     .     gls_diss,Ls_unlmt,Ls_lmt,tke_Q,gls_Q,Gm,Gh,Sm,Sh,cff,ql
#endif
#ifdef DIAG
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  tmp3d_1,tmp3d_2,tmp3d_3,tmp3d_4,tmp3d_5,tmp3d_6,
     .  tmp3d_7,tmp3d_8,tmp3d_9,tmp3d_10
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  tmp2d_1,tmp2d_2,tmp2d_3,tmp2d_4,tmp2d_5
c
      tmp3d_1=0.
      tmp3d_2=0.
      tmp3d_3=0.
      tmp3d_4=0.
      tmp3d_5=0.
      tmp3d_6=0.
      tmp3d_7=0.
      tmp3d_8=0.
      tmp3d_9=0.
      tmp3d_10=0.
      tmp2d_1=0.
      tmp2d_2=0.
      tmp2d_3=0.
      tmp2d_4=0.
      tmp2d_5=0.
#endif
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=-2,jj+3
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(-2,ifp(j,l)),min(ii+3,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(u, 1,2*kk, 2,2, halo_uv)
      call xctilr(v, 1,2*kk, 2,2, halo_vv)
c
c --- Locate the range of layers to be considered in the computation of
c --- diffusivities.
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          kmax(i,j)=1
          do k=3,kk
            if (dp(i,j,k).gt.dpbmin) kmax(i,j)=k
          enddo
          if (kfpla(i,j,n).ge.kmax(i,j)) then
            kfil(i,j)=kfpla(i,j,n)+1
          else
            if (sigma(i,j,kfpla(i,j,n)+nn).lt.
     .          .5*(sigmar(i,j,kfpla(i,j,n)  )
     .             +sigmar(i,j,kfpla(i,j,n)+1))) then
              kfil(i,j)=kfpla(i,j,n)+1
            else
              kfil(i,j)=kfpla(i,j,n)+2
            endif
          endif
        enddo
        enddo
      enddo
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kfil(i,j)
#ifdef DIAG
          tmp2d_1(i,j)=kfpla(i,j,n)
          tmp2d_2(i,j)=kfil(i,j)
#endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kfil(i,j)=nint(util1(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute squared vertical velocity gradients of v-component
c$OMP PARALLEL DO PRIVATE(kfpl,klpl,kn,q,tup)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          kfpl(i)=kk+1
          klpl(i)=1
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (dpv(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,4,-1
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k.ge.max(kfil(i,j-1),kfil(i,j)).and.
     .          dpv(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
        do k=1,kk
          kn=k+nn
          do i=1,ii
            dv2(i,j,k)=0.
            mdv2(i,j,k)=0
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k.ge.kfpl(i).and.k.le.klpl(i).and.
     .          klpl(i)-kfpl(i).ge.1) then
              if     (k.eq.kfpl(i)) then
                q=v(i,j,kn+1)-v(i,j,kn)
                q=q*q
                dv2(i,j,k)=q
                tup(i)=q
              elseif (k.lt.klpl(i)) then
                q=v(i,j,kn+1)-v(i,j,kn)
                q=q*q
                dv2(i,j,k)=.5*(tup(i)+q)
                tup(i)=q
              else
                dv2(i,j,k)=tup(i)
              endif
              mdv2(i,j,k)=1
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ kfpl,klpl,kn,du2,mdu2,q,tup,pup,sup,cr,plo,tlo,slo,bcrrd,afeql,
c$OMP+ bvfbot,dps,drho,bvfsq,bvf,dfints,dfddsu,dfddsl,rig,egr,rhisc,els,
c$OMP+ nus,nub,nut,nuls,tdmf,ust,mols,h,sg,zeta,phis,ws)
      do j=1,jj
c
c ----- Compute squared vertical velocity gradients of u-component
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          kfpl(i)=kk+1
          klpl(i)=1
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (dpu(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,4,-1
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k.ge.min(kfil(i-1,j),kfil(i,j)).and.
     .          dpu(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
        do k=1,kk
          kn=k+nn
          do i=1,ii+1
            du2(i,k)=0.
            mdu2(i,k)=0
          enddo
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k.ge.kfpl(i).and.k.le.klpl(i).and.
     .          klpl(i)-kfpl(i).ge.1) then
              if     (k.eq.kfpl(i)) then
                q=u(i,j,kn+1)-u(i,j,kn)
                q=q*q
                du2(i,k)=q
                tup(i)=q
              elseif (k.lt.klpl(i)) then
                q=u(i,j,kn+1)-u(i,j,kn)
                q=q*q
                du2(i,k)=.5*(tup(i)+q)
                tup(i)=q
              else
                du2(i,k)=tup(i)
              endif
              mdu2(i,k)=1
            endif
          enddo
          enddo
        enddo
c
c ----- Compute the first baroclinic rossby radius of deformation using
c ----- the WKB approximation by Chelton at al. (1998).
c ----- !!! Could include top layer in computation !!!
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pup(i)=.5*(3.*p(i,j,3)-p(i,j,min(kk,kfpla(i,j,n))+1))
          kn=2+nn
          tup(i)=temp(i,j,kn)
          sup(i)=saln(i,j,kn)
          cr(i)=0.
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfpla(i,j,n)) then
              if (p(i,j,kk+1)-p(i,j,k+1).lt.epsil) then
                plo=p(i,j,kk+1)
              else
                plo=.5*(p(i,j,k)+p(i,j,k+1))
              endif
              tlo=temp(i,j,kn)
              slo=saln(i,j,kn)
              cr(i)=cr(i)
     .             +sqrt(max(0.,(rho(p(i,j,k),tlo,slo)
     .                          -rho(p(i,j,k),tup(i),sup(i)))
     .                          *(plo-pup(i))))
              pup(i)=plo
              tup(i)=tlo
              sup(i)=slo
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cr(i)=alpha0*cr(i)/pi
          bcrrd(i)=min(cr(i)/max(1.e-24,abs(coriop(i,j))),
     .                 sqrt(cr(i)/(2.*max(1.e-24,betafp(i,j)))))
#ifdef DIAG
          tmp2d_3(i,j)=bcrrd(i)
#endif
          afeql(i)=max(abs(coriop(i,j)),sqrt(2.*betafp(i,j)*cr(i)))
#ifdef DIAG
          tmp2d_4(i,j)=afeql(i)
#endif
        enddo
        enddo
c
c ----- Compute local gradient richardson number and Brunt-Vaisala
c ----- frequency.
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          bvfbot(i)=0.
          dps(i)=0.
        enddo
        enddo
        do k=4,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .          kmax(i,j)-kfil(i,j).ge.1) then
              if     (k.eq.kfil(i,j)) then
                q=max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1))
     .                  -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drho(i,k)=q
                tup(i)=q
              elseif (k.lt.kmax(i,j)) then
                q=max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1))
     .                  -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drho(i,k)=2.*tup(i)*q/max(1.e-14,tup(i)+q)
                tup(i)=q
              else
                drho(i,k)=tup(i)
              endif
c
c --- ------- Brunt-Vaisala frequency squared
              bvfsq(i,k)=g*g*max(drhomn,drho(i,k))
     .                   /max(epsil,dp(i,j,kn))
#ifdef DIAG
              tmp3d_1(i,j,k)=bvfsq(i,k)
#endif
c
c --- ------- Brunt-Vaisala frequency
              bvf(i,k)=sqrt(bvfsq(i,k))
c
              q=(mdu2(i,k)*du2(i,k)+mdu2(i+1,k)*du2(i+1,k))
     .          /max(1,mdu2(i,k)+mdu2(i+1,k))
     .         +(mdv2(i,j,k)*dv2(i,j,k)+mdv2(i,j+1,k)*dv2(i,j+1,k))
     .          /max(1,mdv2(i,j,k)+mdv2(i,j+1,k))
c
c --- ------- Local gradient Richardson number
              rig(i,k)=alpha0*alpha0*max(drhomn,drho(i,k))*dp(i,j,kn)
     .                 /max(1.e-9,q)
#ifdef DIAG
              tmp3d_2(i,j,k)=rig(i,k)
#endif
#ifdef TKE
              if (dp(i,j,kn).gt.dpbmin) then
                Buoy(i,j,k)=-difdia(i,j,k)*bvfsq(i,k)
                h=max(onem,dp(i,j,kn))*alpha0/g
                Shear2(i,j,k)=max(1.e-9,q)/(h*h)
                Prod(i,j,k)=difdia(i,j,k)*Pr_t*Shear2(i,j,k)
              else
                Buoy(i,j,k)=0.
                Shear2(i,j,k)=1.e-9
                Prod(i,j,k)=0.
              endif
#endif
c
c --- ------- Accumulate Brunt-Vaisala frequency in a region near the
c --- ------- bottom
              q=max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpnbav,p(i,j,k)))
              if (q.gt.0.) then
                bvfbot(i)=bvfbot(i)+bvf(i,k)*q
                dps(i)=dps(i)+q
              endif
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i).gt.0.) then
            bvfbot(i)=bvfbot(i)/dps(i)
#ifdef DIAG
            tmp2d_5(i,j)=bvfbot(i)
#endif
          endif
        enddo
        enddo
c
c --- - Compute diffusivity weigth to reduce eddy diffusivity when the
c --- - Rossby radius is resolved by the grid.
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          q=bcrrd(i)/sqrt(scpx(i,j)*scpx(i,j)+scpy(i,j)*scpy(i,j))
          difwgt(i,j)=1./(1.+q**4)
        enddo
        enddo
c
c --- ------------------------------------------------------------------
c --- - Compute layer interface and isopycnal diffusivities
c --- ------------------------------------------------------------------
c
        if (iidtyp.eq.1) then
c
c --- --- Type 1: Diffusivities are diffusive velocities multiplied by
c --- --- the local horizontal grid scale.
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=sqrt(scp2(i,j))
            difint(i,j,1)=thkdff*q
            difiso(i,j,1)=temdff*q
          enddo
          enddo
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              difint(i,j,k)=difint(i,j,1)
              difiso(i,j,k)=difiso(i,j,1)
            enddo
            enddo
          enddo
c
        else
c
c --- --- Type 2: Diffusivities are parameterized according to Eden and
c --- --- Greatbatch (2008).
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            difint(i,j,1)=egmndf
            dfints(i)=0.
            dps(i)=0.
          enddo
          enddo
          do k=2,kk
            kn=k+nn
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
c --- --------- Eady growth rate.
                egr=afeql(i)/sqrt(rig(i,k)+eggam)
c
c --- --------- Rhines scale.
                rhisc=egr/max(1.e-24,betafp(i,j))
c
c --- --------- Eddy length scale.
                els=min(bcrrd(i),rhisc)
c
c --- --------- Layer interface diffusivity.
                difint(i,j,k)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,egc*egr*els*els*difwgt(i,j)))
#ifdef DIAG
                tmp3d_3(i,j,k)=egr
                tmp3d_4(i,j,k)=rhisc
                tmp3d_5(i,j,k)=els
#endif
c
c --- --------- Accumulate diffusivities in a region below the first
c --- --------- physical layer
                q=max(0.,
     .                min(p(i,j,kfil(i,j))+dpdiav,p(i,j,k+1))-p(i,j,k))
                dps(i)=dps(i)+q
                dfints(i)=dfints(i)+difint(i,j,k)*q
c
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (dps(i).gt.0.) then
              dfints(i)=dfints(i)/dps(i)
            else
              dfints(i)=egmndf
            endif
          enddo
          enddo
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.lt.kfil(i,j)) then
                difint(i,j,k)=dfints(i)
              endif
              difiso(i,j,k)=difint(i,j,k)*egidfq
            enddo
            enddo
          enddo
c
        endif
c
c --- ------------------------------------------------------------------
c --- - Compute diapycnal diffusivity.
c --- ------------------------------------------------------------------
c
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          difdia(i,j,1)=nubmin
          dfddsu(i)=0.
          dfddsl(i)=0.
          dps(i)=0.
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .          kmax(i,j)-kfil(i,j).ge.1) then
c
c --- ------- Background diapycnal mixing.
              if     (bdmtyp.eq.1) then
c
c --- --------- Type 1: Background diffusivity is a constant divided by
c --- --------- Brunt-Vaisala frequency.
                nub=bdmc1/bvf(i,k)
              elseif (bdmtyp.eq.2) then
c
c --- --------- Type 2: Background diffusivity is a constant
                nub=bdmc2
              else
                nub=0.
              endif
              if (iwdflg.eq.1) then
                nub=nub*(1.+(iwdfac-1.)*ficem(i,j))
              endif
c
c --- ------- Latitude dependency of background diapycnal mixing
              if (bdmldp.eq.1) then
                q=max(1.e-9,abs(coriop(i,j)))
                nub=nub*q/cori30*log(2.*bvf0/q)/log(2.*bvf0/cori30)
              endif
c
              nub=max(nubmin,nub)
c
#ifndef TKE
c --- ------- Shear driven diapycnal mixing.
              if (rig(i,k).lt.ri0) then
c
c --- --------- Maximum diffusivity is increased near the bottom to
c --- --------- provide additional mixing of gravity currents.
                q=(p(i,j,kk+1)-p(i,j,k)+.5*dp(i,j,kn))
     .            /min(dpgc,.5*p(i,j,kk+1))
                q=max(0.,1.-q*q)
                q=q*q*q
                nus=q*nug0+(1.-q)*nus0
c
c --- --------- Parameterization of diffusivity as a function of local
c --- --------- gradient richardson number.
                q=rig(i,k)/ri0
                q=max(0.,1.-q*q)
                nus=nus*q*q*q
              else
                nus=0.
              endif
#else
              if (bvfsq(i,k).gt.0.) then  ! stable stratification
                gls_c3=gls_c3minus
              else                        ! unstable stratification
                gls_c3=gls_c3plus
              endif
#  ifndef GLS
              gls_psi(i,j,kn)=max((gls_c1*Prod(i,j,k)
     .                            +gls_c3*Buoy(i,j,k))/gls_c2,
     .                            gls_psi_min)
#  endif              
              tke_epsilon=cmu_fac2*tke(i,j,kn)**(1.5+gls_m/gls_n)
     .                    *gls_psi(i,j,kn)**(-1./gls_n)
              tke_prod=Prod(i,j,k)
              tke_buoy=Buoy(i,j,k)
              tke_Q=tke_epsilon/tke(i,j,kn)
#  ifdef GLS
              gls_prod=(gls_psi(i,j,kn)/tke(i,j,kn))*gls_c1*Prod(i,j,k)
              gls_buoy=(gls_psi(i,j,kn)/tke(i,j,kn))*gls_c3*Buoy(i,j,k)
              gls_diss=(gls_psi(i,j,kn)/tke(i,j,kn))*gls_c2*tke_epsilon
              gls_Q=gls_diss/gls_psi(i,j,kn)
#  endif
#  ifdef GLS
              if (gls_prod+gls_buoy.ge.0.) then
                gls_psi(i,j,kn)=
     .            (gls_psi(i,j,kn)+delt1*(gls_prod+gls_buoy))
     .            /(1.+delt1*gls_Q)
              else
                gls_psi(i,j,kn)=
     .            (gls_psi(i,j,kn)+delt1*gls_prod)
     .            /(1.+delt1*(gls_Q-(gls_buoy/gls_psi(i,j,kn))))
              endif
              gls_psi(i,j,kn)=max(gls_psi(i,j,kn),gls_psi_min)
#  endif
              q=.56**(.5*gls_n)*gls_cmu0**gls_p
     .          *tke(i,j,kn)**(gls_m+.5*gls_n)
     .          *bvf(i,k)**(-gls_n)
              if (gls_n.lt.0.) then
                gls_psi(i,j,kn)=max(gls_psi(i,j,kn),q)
              else
                gls_psi(i,j,kn)=min(gls_psi(i,j,kn),q)
              endif
c
              tke_epsilon=cmu_fac2*tke(i,j,kn)**(1.5+gls_m/gls_n)
     .                    *gls_psi(i,j,kn)**(-1./gls_n)
              tke_Q=tke_epsilon/tke(i,j,kn)
c
              if (tke_prod+tke_buoy.ge.0.) then
                tke(i,j,kn)=(tke(i,j,kn)+delt1*(tke_prod+tke_buoy))
     .                      /(1.+delt1*tke_Q)
              else
                tke(i,j,kn)=(tke(i,j,kn)+delt1*tke_prod)
     .                      /(1.+delt1*(tke_Q-(tke_buoy/tke(i,j,kn))))
                tke(i,j,kn)=max(tke(i,j,kn),tke_min)
              endif
c
c --- ------- Bottom Boundary Conditions
              if (k.eq.kmax(i,j)) then
                ust=max(ustarb(i,j),ustmin)
                tke(i,j,kn)=max(tke_min,(ust/gls_cmu0)**2)
#  ifdef GLS
                gls_psi(i,j,kn)=max(gls_psi_min,
     .                              (gls_cmu0**(gls_p-2.*gls_m))
     .                              *(ust**(2.*gls_m))
     .                              *(kappa*1.e2)**gls_n)
#  endif
              endif
c
              Ls_unlmt=max(Ls_unlmt_min,
     .                     cmu_fac1*gls_psi(i,j,kn)**(gls_exp1)
     .                     *tke(i,j,kn)**(-tke_exp1))

              if (bvfsq(i,k).gt.0.) then  ! stable stratification
                Ls_lmt=min(Ls_unlmt,
     .                     sqrt(.56*tke(i,j,kn)/max(bvfsq(i,k),1.e-10)))
c               Ls_lmt=Ls_unlmt
              else                        ! unstable stratification
                Ls_lmt=Ls_unlmt
              endif
c
c --- ------- Compute nondimensional stability functions for tracers
c --- ------- (Sh) and momentum (Sm). Canuto-A
              Gh=min(gls_Gh0,-bvfsq(i,k)*Ls_lmt*Ls_lmt/(2.*tke(i,j,kn)))
              Gh=min(Gh,(Gh-(Gh-gls_Ghcri)**2)
     .                  /(Gh+gls_Gh0-2.*gls_Ghcri))
              Gh=max(Gh,gls_Ghmin)
c
c --- ------- Compute shear number.
              Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh**2))
     .           /(gls_b2-gls_b4*gls_fac6*Gh)
              Gm=min(Gm,Shear2(i,j,k)*Ls_lmt*Ls_lmt/(2.*tke(i,j,kn)))
c
c --- ------- Compute stability functions
              cff=gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm
     .           +gls_b3*gls_fac6**2*Gh**2-gls_b4*gls_fac6**2*Gh*Gm
     .           +gls_b5*gls_fac6**2*Gm*Gm
              Sm=(gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/cff
              Sh=(gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/cff
              Sm=max(Sm,0.)
              Sh=max(Sh,0.)
c
c --- ------- Relate Canuto stability to MICOM notation
              Sm=Sm*cmu_fac3/gls_cmu0**3
              Sh=Sh*cmu_fac3/gls_cmu0**3
c
              ql=sqrt2*.5*(Ls_lmt+L_scale(i,j,k))*sqrt(tke(i,j,kn))
c
              nus=Sh*ql
c             nus=Sh*(tke(i,j,k)*tke(i,j,k))/gls_psi(i,j,k)
              L_scale(i,j,k)=max(Ls_lmt,Ls_unlmt_min)
#  ifdef GLS
c
c --- ------- Recompute gls based on limited length scale
              gls_psi(i,j,kn)=max(gls_cmu0**gls_p*tke(i,j,kn)**gls_m
     .                            *L_scale(i,j,k)**gls_n,gls_psi_min)
#  endif
#endif
c
c --- ------- Tidally driven diapycnal mixing
              if (tdmflg.eq.1) then
                q=.5*(tanh(4.*(abs(plat(i,j))-tdclat)/tddlat-2.)+1.)
                q=(1.-q)*tdmls0+q*tdmls1
                if (dp(i,j,kn).lt.epsil) then
                  tdmf=exp(p(i,j,k)/q)/(q*(exp(p(i,j,kk+1)/q)-1.))
                else
                  tdmf=(exp(p(i,j,k+1)/q)-exp(p(i,j,k)/q))
     .                 /(dp(i,j,kn)*(exp(p(i,j,kk+1)/q)-1.))
                endif
                nut=g*tdmq*dmxeff*twedon(i,j)*bvfbot(i)*tdmf/bvfsq(i,k)
              else
                nut=0.
              endif
c
c --- ------- Diapycnal mixing when local stability is weak
              if (drho(i,k).lt.drho0) then
                q=drho(i,k)/drho0
                q=max(0.,1.-q*q)
                nuls=nuls0*q*q*q
              else
               nuls=0.
              endif
c
c --- ------- Total diapycnal diffusivity.
              difdia(i,j,k)=nub+nus+nut+nuls
#ifdef DIAG
              tmp3d_7(i,j,k)=nub
              tmp3d_8(i,j,k)=nus
              tmp3d_9(i,j,k)=nut
              tmp3d_10(i,j,k)=nuls
#endif
c
c --- ------- Accumulate diffusivities in a region below the first
c --- ------- physical layer
              q=max(0.,min(p(i,j,kfil(i,j))+dpddav,p(i,j,k+1))-p(i,j,k))
              dps(i)=dps(i)+q
              dfddsu(i)=dfddsu(i)+nub*q
              dfddsl(i)=dfddsl(i)+difdia(i,j,k)*q
c
            else
              difdia(i,j,k)=difdia(i,j,k-1)
#ifdef TKE
c             tke(i,j,kn)=tke_min
c             L_scale(i,j,k)=Ls_unlmt_min
              tke(i,j,kn)=tke(i,j,kn-1)
              L_scale(i,j,k)=L_scale(i,j,k-1)
#  ifdef GLS
c             gls_psi(i,j,kn)=gls_psi_min
              gls_psi(i,j,kn)=gls_psi(i,j,kn-1)
#  endif
#endif
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i).gt.0.) then
            dfddsu(i)=dfddsu(i)/dps(i)
            dfddsl(i)=dfddsl(i)/dps(i)
          else
            dfddsu(i)=nubmin
            dfddsl(i)=nubmin
          endif
        enddo
        enddo
        do k=2,kk-1
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.lt.kfil(i,j)) then
              if (k.gt.2.and.kfil(i,j).le.kk.and.
     .            p(i,j,min(kk,kfil(i,j)))-p(i,j,3).gt.epsil) then
                q=.5*(p(i,j,k+1)+p(i,j,k))
                difdia(i,j,k)=((q-p(i,j,3))*dfddsl(i)
     .                        +(p(i,j,kfil(i,j))-q)*dfddsu(i))
     .                        /(p(i,j,kfil(i,j))-p(i,j,3))
              else
                difdia(i,j,k)=dfddsu(i)
              endif
            endif
          enddo
          enddo
        enddo
c
c --- - Diffusivity at the lower interface of the top layer
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Lower bounded friction velocity
          ust=max(ustmin,ustar(i,j))
c
c --- --- Monin-Obukhov length scale
          mols=ust**3
     .         /(kappa*sign(max(abs(buoyfl(i,j)),bfeps),-buoyfl(i,j)))
c
c --- --- Mixed layer thickness
          h=(p(i,j,3)-p(i,j,1))/onecm
c
c --- --- Dimensionless vertical coordinate in the boundary layer
          sg=(p(i,j,2)-p(i,j,1))/(p(i,j,3)-p(i,j,1))
c
c --- --- Velocity scale
          if (mols.lt.0.) then
            zeta=min(sleps,sg)*h/mols
            if (zeta.gt.zetas) then
              phis=(1.-16.*zeta)**(-1./2.)
            else
              phis=(as-cs*zeta)**(-1./3.)
            endif
          else
            zeta=sg*h/mols
            phis=1.+5.*zeta
          endif
          ws=kappa*ust/phis
c
          difdia(i,j,1)=h*ws*sg*(1.-sg)**2
        enddo
        enddo
c
      enddo
c$OMP END PARALLEL DO
c
#ifdef DIAG
      call ncfopn('difest.nc','w','c',1)
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncdims('z',kdm)
      call ncwrtr('coriop','x y',coriop,ip,1,1.,0.,8)
      call ncwrtr('betafp','x y',betafp,ip,1,1.,0.,8)
      call ncwrtr('kfpla','x y',tmp2d_1,ip,1,1.,0.,8)
      call ncwrtr('kfil','x y',tmp2d_2,ip,1,1.,0.,8)
      call ncwrtr('bcrrd','x y',tmp2d_3,ip,1,1.,0.,8)
      call ncwrtr('afeql','x y',tmp2d_4,ip,1,1.,0.,8)
      call ncwrtr('bvfbot','x y',tmp2d_5,ip,1,1.,0.,8)
      call ncwrtr('dp','x y z',dp(1-nbdy,1-nbdy,k1n),ip,1,1.,0.,8)
      call ncwrtr('temp','x y z',temp(1-nbdy,1-nbdy,k1n),ip,1,1.,0.,8)
      call ncwrtr('saln','x y z',saln(1-nbdy,1-nbdy,k1n),ip,1,1.,0.,8)
      call ncwrtr('u','x y z',u(1-nbdy,1-nbdy,k1n),iu,1,1.,0.,8)
      call ncwrtr('v','x y z',v(1-nbdy,1-nbdy,k1n),iv,1,1.,0.,8)
      call ncwrtr('bvfsq','x y z',tmp3d_1,ip,1,1.,0.,8)
      call ncwrtr('rig','x y z',tmp3d_2,ip,1,1.,0.,8)
      call ncwrtr('egr','x y z',tmp3d_3,ip,1,1.,0.,8)
      call ncwrtr('rhisc','x y z',tmp3d_4,ip,1,1.,0.,8)
      call ncwrtr('els','x y z',tmp3d_5,ip,1,1.,0.,8)
      call ncwrtr('difint','x y z',difint,ip,1,1.,0.,8)
      call ncwrtr('difdia','x y z',difdia,ip,1,1.,0.,8)
      call ncwrtr('nub','x y z',tmp3d_7,ip,1,1.,0.,8)
      call ncwrtr('nus','x y z',tmp3d_8,ip,1,1.,0.,8)
      call ncwrtr('nut','x y z',tmp3d_9,ip,1,1.,0.,8)
      call ncwrtr('nuls','x y z',tmp3d_10,ip,1,1.,0.,8)
#  ifdef TKE
      call ncwrtr('tke','x y z',tke(1-nbdy,1-nbdy,k1n),ip,1,1.,0.,8)
#    ifdef GLS
      call ncwrtr('gls_psi','x y z',gls_psi(1-nbdy,1-nbdy,k1n),
     .            ip,1,1.,0.,8)
#    endif
#  endif
      call ncfcls
      call xcstop('(difest)')
             stop '(difest)'
#endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest:'
        endif
        call chksummsk(difint,ip,kk,'difint')
        call chksummsk(difiso,ip,kk,'difiso')
        call chksummsk(difdia,ip,kk,'difdia')
#ifdef TKE
        call chksummsk(tke,ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(gls_psi,ip,2*kk,'gls_psi')
#  endif
#endif
      endif
c
      return
      end
