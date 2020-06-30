! ------------------------------------------------------------------------------
! Copyright (C) 2009-2020 Mats Bentsen, Mehmet Ilicak
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

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
      use mod_dia, only : iotype
#endif
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_seaice.h"
#include "common_forc.h"
#include "common_geo.h"
c
c --- parameters:
c ---   iidtyp - type of interface and isopycnal diffusivities. If
c ---            iidtyp=1 the diffusivities are diffusive velocities
c ---            multiplied by the local horizontal grid scale, if
c ---            iidtyp=2 the diffusivities are parameterized according
c ---            to Eden and Greatbatch (2008).
c ---   bdmldp - If bdmldp=1, make the background mixing latitude
c ---            dependent according to Gregg et al. (2003).
c ---   tdmflg - If tdmflg=1, apply tidally driven diapycnal mixing.
c ---   iwdflg - If iwdflg=1, reduce background diapycnal diffusivity
c ---            due to internal wave damping under sea-ice.
c ---   ipfac  - inertial period factor for the velocity averaging time
c ---            scale [].
c ---   dpbmin - smallest layer thickness allowed in evaluating
c ---            local gradient richardson number [g/cm/s**2].
c ---   drhomn - minimum density difference in evaluations the
c ---            Brunt-Vaisala frequency and the local gradient
c ---            Richardson number [g/cm*3].
c ---   thkdff - diffusive velocity for thickness diffusion [cm/s].
c ---   temdff - diffusive velocity for tracer isopycnal diffusion
c ---            [cm/s].
c ---   nu0    - diapycnal diffusivity when range of isopycnic physical
c ---            layers is restricted [cm**2/s].
c ---   nus0   - maximum shear driven diapycnal diffusivity
c ---            [cm**2/s].
c ---   nug0   - maximum gravity current diapycnal diffusivity
c ---            [cm**2/s].
c ---   drho0  - critical local interface density difference [g/cm**3]
c ---   nuls0  - maximum diapycnal diffusivity applied when local
c ---            stability is weak [cm**2/s].
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
c ---   tkepls - length scale of surface TKE penetration beneath the
c ---            mixed layer [g/cm/s**2]
c ---   niwls  - near-inertial waves driven mixing length scale
c ---            beneath the mixed layer [g/cm/s**2].
c ---   cori10 - coriolis parameter at 10N [1/s].
c ---   cori30 - coriolis parameter at 30N [1/s].
c ---   bvf0   - reference stratification in the parameterization of
c ---            latitude dependent background diapycnal mixing [1/s].
c ---   nubmin - minimum background diapycnal diffusivity [cm**2/s].
c ---   dpgc   - thickness of region near the bottom where the maximum
c ---            diffusivity is increased due to gravity current mixing
c ---            processes [g/cm/s**2].
c ---   dpgrav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate upper ocean Eady growth rate
c ---            [g/cm/s**2].
c ---   dpdiav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate lateral diffusivities in the
c ---            non-isopycnic layers [g/cm/s**2].
c ---   dpddav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate diapycnal diffusivities in the
c ---            non-isopycnic layers [g/cm/s**2].
c ---   dpnbav - thickness of region near the bottom used to estimate
c ---            bottom Brunt-Vaisala frequency [g/cm/s**2].
      integer iidtyp,bdmldp,tdmflg,iwdflg
      real ipfac,dptmin,dpbmin,drhomn,thkdff,temdff,nu0,nus0,nug0,drho0,
     .     nuls0,iwdfac,dmxeff,tdmq,tdmls0,tdmls1,tdclat,tddlat,tkepls,
     .     niwls,cori10,cori30,bvf0,nubmin,dpgc,dpgrav,dpdiav,dpddav,
     .     dpnbav,ustmin,kappa,bfeps,sleps,zetas,as,cs
      parameter (iidtyp=2,bdmldp=1,tdmflg=1,iwdflg=1,ipfac=2.,
     .           dptmin=98060.,dpbmin=980.6,drhomn=6.e-6,
     .           thkdff=.5,temdff=.35,nu0=.1,nus0=50.,nug0=2500.,
     .           drho0=6.e-6,nuls0=500.,iwdfac=.06,dmxeff=.2,tdmq=1./3.,
     .           tdmls0=500.*98060.,tdmls1=100.*98060.,tdclat=74.5,
     .           tddlat=3.,tkepls=20.*98060.,niwls=300.*98060.,
     .           cori10=2.5256e-5,cori30=7.2722e-5,bvf0=5.24e-3,
     .           nubmin=.01,dpgc=300.*98060.,dpgrav=100.*98060.,
     .           dpdiav=100.*98060.,dpddav=10.*98060.,
     .           dpnbav=250.*98060.,ustmin=.1,kappa=.4,bfeps=1.e-12,
     .           sleps=.1,zetas=-1.,as=-28.86,cs=98.96)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  dv2
      real, dimension(1-nbdy:idm+nbdy,kdm) ::
     .  du2,drho,bvfsq,bvf,rig,egr
      real, dimension(1-nbdy:idm+nbdy) ::
     .  tup,pup,sup,cr,bcrrd,afeql,bvfbot,dps,egrs,egrup,dfints,urmse,
     .  cpse,dfddsu,dfddsl
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  mskv
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  kfil,kmax
      integer, dimension(1-nbdy:idm+nbdy,kdm) ::
     .  msku
      integer, dimension(1-nbdy:idm+nbdy) ::
     .  kfpl,klpl
      integer i,j,k,l,kn,mmm
      real ubt,uml1t,uml2t,q,uml1a,uml2a,vbt,vml1t,vml2t,vml1a,vml2a,
     .     plo,tlo,slo,rhisc,els,egrlo,umnsc,esfac,nus,nub,nut,nuls,
     .     vsf,nusm,ust,mols,h,sg,zeta,phis,ws
c
      real rho,sig
      external rho,sig
c
#if defined(TRC) && defined(TKE)
#  include "param_trc.h"
#  include "common_trc.h"
#  include "param_tke.h"
#  include "common_tke.h"
      real gls_c3,tke_prod,tke_buoy,tke_epsilon,gls_prod,gls_buoy,
     .     gls_diss,Ls_unlmt,Ls_lmt,tke_Q,gls_Q,Gm,Gh,Sm,Sh,cff,ql
#endif
#ifdef DIAG
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  tmp3d_1,tmp3d_2,tmp3d_3,tmp3d_4,tmp3d_5,tmp3d_6,
     .  tmp3d_7,tmp3d_8,tmp3d_9,tmp3d_10,tmp3d_11,tmp3d_12
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  tmp2d_1,tmp2d_2,tmp2d_3,tmp2d_4,tmp2d_5,tmp2d_6,tmp2d_7
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
      tmp3d_11=0.
      tmp3d_12=0.
      tmp2d_1=0.
      tmp2d_2=0.
      tmp2d_3=0.
      tmp2d_4=0.
      tmp2d_5=0.
      tmp2d_6=0.
      tmp2d_7=0.
#endif
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
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
      call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
c
c --- ------------------------------------------------------------------
c --- Estimate energy input by near-inertial motions by considering the
c --- mixed layer kinetic energy change over a model time step of
c --- velocity deviations from time averaged velocities. The velocity
c --- averaging time scale is the inertial period times -ipfac- with
c --- the inertial period limited to the one at 10N/S.
c --- ------------------------------------------------------------------
c
      mmm=(m-1)*2
c
c$OMP PARALLEL DO PRIVATE(l,i,ubt,uml1t,uml2t,q,uml1a,uml2a)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
c
          ubt=ubflxs_p(i,j,m)*dlt/(delt1*scuy(i,j)*pbu(i,j,m))
          uml1t=u(i,j,1+mm)+ubt
          uml2t=u(i,j,2+mm)+ubt
c
          q=delt1*max(cori10,abs(.5*(coriop(i-1,j)+coriop(i,j))))
     .      /(ipfac*2.*pi)
          umlres(i,j,1)=umlres(i,j,1)+uml1t
          uml1a=umlres(i,j,1)*q
          umlres(i,j,1)=umlres(i,j,1)*(1.-q)
          umlres(i,j,2)=umlres(i,j,2)+uml2t
          uml2a=umlres(i,j,2)*q
          umlres(i,j,2)=umlres(i,j,2)*(1.-q)
c
          util1(i,j)=((uml1t-uml1a)**2-(uml(i,j,1+mmm)-uml1a)**2)
     .               *dpu(i,j,1+mm)
     .              +((uml2t-uml2a)**2-(uml(i,j,2+mmm)-uml2a)**2)
     .               *dpu(i,j,2+mm)

          uml(i,j,1+mmm)=uml1t
          uml(i,j,2+mmm)=uml2t
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,vbt,vml1t,vml2t,q,vml1a,vml2a)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c
          vbt=vbflxs_p(i,j,m)*dlt/(delt1*scvx(i,j)*pbv(i,j,m))
          vml1t=v(i,j,1+mm)+vbt
          vml2t=v(i,j,2+mm)+vbt
c
          q=delt1*max(cori10,abs(.5*(coriop(i,j-1)+coriop(i,j))))
     .      /(ipfac*2.*pi)
          vmlres(i,j,1)=vmlres(i,j,1)+vml1t
          vml1a=vmlres(i,j,1)*q
          vmlres(i,j,1)=vmlres(i,j,1)*(1.-q)
          vmlres(i,j,2)=vmlres(i,j,2)+vml2t
          vml2a=vmlres(i,j,2)*q
          vmlres(i,j,2)=vmlres(i,j,2)*(1.-q)
c
          util2(i,j)=((vml1t-vml1a)**2-(vml(i,j,1+mmm)-vml1a)**2)
     .               *dpv(i,j,1+mm)
     .              +((vml2t-vml2a)**2-(vml(i,j,2+mmm)-vml2a)**2)
     .               *dpv(i,j,2+mm)

          vml(i,j,1+mmm)=vml1t
          vml(i,j,2+mmm)=vml2t
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          idkedt(i,j)=abs((util1(i,j)*iu(i,j)+util1(i+1,j)*iu(i+1,j))
     .                    /max(1,iu(i,j)+iu(i+1,j))
     .                   +(util2(i,j)*iv(i,j)+util2(i,j+1)*iv(i,j+1))
     .                    /max(1,iv(i,j)+iv(i,j+1)))
     .                *alpha0/(2.*g*delt1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Estimate friction velocity cubed.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          ustar3(i,j)=ustar(i,j)**3
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Locate the range of layers to be considered in the computation of
c --- diffusivities.
      do j=0,jj+1
        do i=0,ii+1
          kmax(i,j)=0
        enddo
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmax(i,j)=1
          do k=3,kk
            kn=k+nn
            if (dp(i,j,kn).gt.dpbmin) kmax(i,j)=k
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
c$OMP PARALLEL DO PRIVATE(l,i)
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
c$OMP PARALLEL DO PRIVATE(l,i)
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
c$OMP PARALLEL DO PRIVATE(l,i,kfpl,klpl,k,kn,q,tup)
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
            mskv(i,j,k)=0
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
              mskv(i,j,k)=1
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,kfpl,klpl,k,kn,du2,msku,q,tup,pup,sup,cr,plo,tlo,slo,bcrrd,
c$OMP+ afeql,bvfbot,dps,drho,bvfsq,bvf,rig,egrs,egr,egrup,egrlo,dfints,
c$OMP+ rhisc,els,urmse,cpse,umnsc,esfac,dfddsu,dfddsl,nub,nus,ust,vsf,
c$OMP+ nut,nuls,nusm,mols,h,sg,zeta,phis,ws
#if defined(TRC) && defined(TKE)
c$OMP+ ,gls_c3,tke_epsilon,tke_prod,tke_buoy,tke_Q,Ls_unlmt,Ls_lmt,Gh,
c$OMP+ Gm,cff,Sm,Sh,ql
# ifdef GLS
c$OMP+ ,gls_prod,gls_buoy,gls_diss,gls_Q
# endif
#endif
c$OMP+ )
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
            msku(i,k)=0
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
              msku(i,k)=1
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
          bcrrd(i)=
     .      sqrt(cr(i)*cr(i)
     .           /max(coriop(i,j)*coriop(i,j)+2.*betafp(i,j)*cr(i),
     .                1.e-24))
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
              q=(msku(i,k)*du2(i,k)+msku(i+1,k)*du2(i+1,k))
     .          /max(1,msku(i,k)+msku(i+1,k))
     .         +(mskv(i,j,k)*dv2(i,j,k)+mskv(i,j+1,k)*dv2(i,j+1,k))
     .          /max(1,mskv(i,j,k)+mskv(i,j+1,k))
c
c --- ------- Local gradient Richardson number
              rig(i,k)=alpha0*alpha0*max(drhomn,drho(i,k))*dp(i,j,kn)
     .                 /max(1.e-9,q)
#ifdef DIAG
              tmp3d_2(i,j,k)=rig(i,k)
#endif
#if defined(TRC) && defined(TKE)
              if (dp(i,j,kn).gt.dpbmin) then
                Buoy(i,j,k)=-difdia(i,j,k)*bvfsq(i,k)
                h=max(onem,dp(i,j,kn))*alpha0/g
c               h=max(onem*1e-8,dp(i,j,kn))*alpha0/g
c               h=max(onemm,dp(i,j,kn))*alpha0/g
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
        if     (edwmth.eq.'smooth') then
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j)
     .                         +scpy(i,j)*scpy(i,j)))
            difwgt(i,j)=1./(1.+.25*q**4)
          enddo
          enddo
        elseif (edwmth.eq.'step') then
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j)
     .                         +scpy(i,j)*scpy(i,j)))
            if (q.le.2.) then
              difwgt(i,j)=1.
            else
              difwgt(i,j)=0.
            endif
          enddo
          enddo
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' edwmth=',trim(edwmth),
     .                        ' is unsupported!'
          endif
          call xcstop('(difest)')
                 stop '(difest)'
        endif
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
c
c --- --- Eady growth rate.
          if (edsprs) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              egrs(i)=0.
              dps(i)=0.
            enddo
            enddo
          endif
          if     (edritp.eq.'shear') then
            do k=2,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .              kmax(i,j)-kfil(i,j).ge.1) then
                  egr(i,k)=afeql(i)/sqrt(rig(i,k)+eggam)
#ifdef DIAG
                  tmp3d_3(i,j,k)=egr(i,k)
#endif
                  if (edsprs) then
                    q=max(0.,min(p(i,j,kfil(i,j))+dpgrav,
     .                           p(i,j,k+1))-p(i,j,k))
                    dps(i)=dps(i)+q
                    egrs(i)=egrs(i)+egr(i,k)*q
                  endif
                endif
              enddo
              enddo
            enddo
          elseif (edritp.eq.'large scale') then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (kmax(i,j)-kfil(i,j).ge.1) then
                k=kfil(i,j)
                if     (kmax(i-1,j).ge.k.and.kmax(i+1,j).ge.k) then
                  q=.25*(nnslpx(i,j,k)+nnslpx(i+1,j,k))**2
                elseif (kmax(i-1,j).ge.k) then
                  q=nnslpx(i,j,k)**2
                elseif (kmax(i+1,j).ge.k) then
                  q=nnslpx(i+1,j,k)**2
                else
                  q=0.
                endif
                if     (kmax(i,j-1).ge.k.and.kmax(i,j+1).ge.k) then
                  q=q+.25*(nnslpy(i,j,k)+nnslpy(i,j+1,k))**2
                elseif (kmax(i,j-1).ge.k) then
                  q=q+nnslpy(i,j,k)**2
                elseif (kmax(i,j+1).ge.k) then
                  q=q+nnslpy(i,j+1,k)**2
                endif
                egrup(i)=sqrt(q)
              endif
            enddo
            enddo
            do k=2,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (kmax(i,j)-kfil(i,j).ge.1) then
                  if     (k.ge.kfil(i,j).and.k.lt.kmax(i,j)) then
                    if     (kmax(i-1,j).gt.k.and.kmax(i+1,j).gt.k) then
                      q=.25*(nnslpx(i,j,k+1)+nnslpx(i+1,j,k+1))**2
                    elseif (kmax(i-1,j).gt.k) then
                      q=nnslpx(i,j,k+1)**2
                    elseif (kmax(i+1,j).gt.k) then
                      q=nnslpx(i+1,j,k+1)**2
                    else
                      q=0.
                    endif
                    if     (kmax(i,j-1).gt.k.and.kmax(i,j+1).gt.k) then
                      q=q+.25*(nnslpy(i,j,k+1)+nnslpy(i,j+1,k+1))**2
                    elseif (kmax(i,j-1).gt.k) then
                      q=q+nnslpy(i,j,k+1)**2
                    elseif (kmax(i,j+1).gt.k) then
                      q=q+nnslpy(i,j+1,k+1)**2
                    endif
                    egrlo=sqrt(q)
                    egr(i,k)=.5*(egrup(i)+egrlo)
                    egrup(i)=egrlo
#ifdef DIAG
                    tmp3d_3(i,j,k)=egr(i,k)
#endif
                    if (edsprs) then
                      q=max(0.,min(p(i,j,kfil(i,j))+dpgrav,
     .                             p(i,j,k+1))-p(i,j,k))
                      dps(i)=dps(i)+q
                      egrs(i)=egrs(i)+egr(i,k)*q
                    endif
                  elseif (k.eq.kmax(i,j)) then
                    egr(i,k)=egr(i,k-1)
#ifdef DIAG
                    tmp3d_3(i,j,k)=egr(i,k)
#endif
                  endif
                endif
              enddo
              enddo
            enddo
          else
            if (mnproc.eq.1) then
              write (lp,'(3a)') ' edritp=',trim(edritp),
     .                          ' is unsupported!'
            endif
            call xcstop('(difest)')
                   stop '(difest)'
          endif
          if (edsprs) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (dps(i).gt.0.) then
                egrs(i)=egrs(i)/dps(i)
              else
                egrs(i)=0.
              endif
            enddo
            enddo
          endif
c
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            difint(i,j,1)=egmndf
            dfints(i)=0.
            dps(i)=0.
          enddo
          enddo
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
c --- --------- Rhines scale.
                rhisc=egr(i,k)/max(1.e-24,betafp(i,j))
c
c --- --------- Eddy length scale.
                els=max(eglsmn,min(bcrrd(i),rhisc))
c
c --- --------- Temporary layer interface diffusivity.
                difint(i,j,k)=egc*egr(i,k)*els*els
#ifdef DIAG
                tmp3d_4(i,j,k)=rhisc
                tmp3d_5(i,j,k)=els
#endif
c
c --- --------- Accumulate diffusivities in a region below the first
c --- --------- physical layer.
                q=max(0.,min(p(i,j,kfil(i,j))+dpdiav,
     .                       p(i,j,k+1))-p(i,j,k))
                dps(i)=dps(i)+q
                dfints(i)=dfints(i)+difint(i,j,k)*q
c
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
c
c --- --- Apply eddy diffusivity limiting, suppression when the Rossby
c --- --- radius is resolved by the grid, and suppression away from
c --- --- steering levels if requested.
c
c --- --- Eddy diffusivity modification of surface non-isopycnic
c --- --- layers.
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
            if (edsprs) then
c
c --- ------- RMS eddy velocity estimated from K = Gamma*u_rms*L, where
c --- ------- a mixing efficiency of Gamma = 0.35 is used (Klocker and
c --- ------- Abernathey, 2014).
              rhisc=egrs(i)/max(1.e-24,betafp(i,j))
              els=max(eglsmn,min(bcrrd(i),rhisc))
              urmse(i)=2.86*egc*egrs(i)*els
#ifdef DIAG
              tmp2d_6(i,j)=urmse(i)
#endif
c
c --- ------- Zonal eddy phase speed minus zonal barotropic velocity
c --- ------- with a lower bound of -20 cm s-1.
              cpse(i)=max(-20.,-betafp(i,j)*bcrrd(i)**2)
#ifdef DIAG
              tmp2d_7(i,j)=cpse(i)
#endif
c
            endif
c
            if (dps(i).gt.0.) then
c
              if (edsprs) then
c
c --- --------- Zonal mixed layer velocity minus eddy phase speed. Note
c --- --------- that only the baroclinic component is used since the
c --- --------- barotropic velocity is subtracted from the estimate of
c --- --------- eddy phase speed.
                if     (ip(i-1,j)+ip(i+1,j).eq.2) then
                  q=.5*((u(i  ,j,1+nn)*dpu(i  ,j,1+nn)
     .                  +u(i  ,j,2+nn)*dpu(i  ,j,2+nn))
     .                  /(dpu(i  ,j,1+nn)+dpu(i  ,j,2+nn))
     .                 +(u(i+1,j,1+nn)*dpu(i+1,j,1+nn)
     .                  +u(i+1,j,2+nn)*dpu(i+1,j,2+nn))
     .                  /(dpu(i+1,j,1+nn)+dpu(i+1,j,2+nn)))
                elseif (ip(i-1,j).eq.1) then
                  q=(u(i  ,j,1+nn)*dpu(i  ,j,1+nn)
     .              +u(i  ,j,2+nn)*dpu(i  ,j,2+nn))
     .              /(dpu(i  ,j,1+nn)+dpu(i  ,j,2+nn))
                elseif (ip(i+1,j).eq.1) then
                  q=(u(i+1,j,1+nn)*dpu(i+1,j,1+nn)
     .              +u(i+1,j,2+nn)*dpu(i+1,j,2+nn))
     .              /(dpu(i+1,j,1+nn)+dpu(i+1,j,2+nn))
                else
                  q=0.
                endif
                umnsc=q*cosang(i,j)
                if     (ip(i,j-1)+ip(i,j+1).eq.2) then
                  q=.5*((v(i,j  ,1+nn)*dpv(i,j  ,1+nn)
     .                  +v(i,j  ,2+nn)*dpv(i,j  ,2+nn))
     .                  /(dpv(i,j  ,1+nn)+dpv(i,j  ,2+nn))
     .                 +(v(i,j+1,1+nn)*dpv(i,j+1,1+nn)
     .                  +v(i,j+1,2+nn)*dpv(i,j+1,2+nn))
     .                  /(dpv(i,j+1,1+nn)+dpv(i,j+1,2+nn)))
                elseif (ip(i,j-1).eq.1) then
                  q=(v(i,j  ,1+nn)*dpv(i,j  ,1+nn)
     .              +v(i,j  ,2+nn)*dpv(i,j  ,2+nn))
     .              /(dpv(i,j  ,1+nn)+dpv(i,j  ,2+nn))
                elseif (ip(i,j+1).eq.1) then
                  q=(v(i,j+1,1+nn)*dpv(i,j+1,1+nn)
     .              +v(i,j+1,2+nn)*dpv(i,j+1,2+nn))
     .              /(dpv(i,j+1,1+nn)+dpv(i,j+1,2+nn))
                else
                  q=0.
                endif
                umnsc=umnsc-q*sinang(i,j)-cpse(i)
#ifdef DIAG
                tmp3d_6(i,j,1)=umnsc
#endif
c
c --- --------- Eddy mixing suppresion factor where lower bounds of
c --- --------- zonal velocity minus eddy phase speed and absolute value
c --- --------- of RMS eddy velocity is set to -20 cm s-1 and 5 cm s-1,
c --- --------- respectively.
                esfac=1./(1.+4.*(umnsc/max(5.,abs(urmse(i))))**2)
#ifdef DIAG
                tmp3d_7(i,j,1)=esfac
#endif
c
              else
                esfac=1.
              endif
c
              dfints(i)=dfints(i)/dps(i)
              dfints(i)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,dfints(i)*difwgt(i,j)*esfac))
            else
              dfints(i)=egmndf
            endif
          enddo
          enddo
c
c --- --- Eddy diffusivity modification of isopycnic layers.
          do k=2,kk
            kn=k+nn
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
                if (edsprs) then
c
c --- ----------- Zonal velocity minus eddy phase speed.
                  umnsc=
     .              (msku(i,k)  *u(i,j,kn)+msku(i+1,k)  *u(i+1,j,kn))
     .              /max(1,msku(i,k)  +msku(i+1,k)  )*cosang(i,j)
     .             -(mskv(i,j,k)*v(i,j,kn)+mskv(i,j+1,k)*v(i,j+1,kn))
     .              /max(1,mskv(i,j,k)+mskv(i,j+1,k))*sinang(i,j)
     .             -cpse(i)
#ifdef DIAG
                  tmp3d_6(i,j,k)=umnsc
#endif
c
c --- ----------- Eddy mixing suppresion factor.
                  esfac=1./(1.+4.*(umnsc/max(5.,abs(urmse(i))))**2)
#ifdef DIAG
                  tmp3d_7(i,j,k)=esfac
#endif
c
                else
                  esfac=1.
                endif
c
                difint(i,j,k)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,difint(i,j,k)*difwgt(i,j)*esfac))
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
c
c --- --- Set isopycnal tracer diffusivity proportional to the layer
c --- --- interface diffusivity by the factor EGIDFQ.
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
          difdia(i,j,1)=nu0
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
#if !defined(TRC) || !defined(TKE)
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
              trc(i,j,kn,itrgls)=max((gls_c1*Prod(i,j,k)
     .                               +gls_c3*Buoy(i,j,k))/gls_c2,
     .                               gls_psi_min)
#  endif              
              tke_epsilon=cmu_fac2*trc(i,j,kn,itrtke)**(1.5+gls_m/gls_n)
     .                    *trc(i,j,kn,itrgls)**(-1./gls_n)
              tke_prod=Prod(i,j,k)
              tke_buoy=Buoy(i,j,k)
              tke_Q=tke_epsilon/trc(i,j,kn,itrtke)
#  ifdef GLS
              gls_prod=(trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke))
     .                 *gls_c1*Prod(i,j,k)
              gls_buoy=(trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke))
     .                 *gls_c3*Buoy(i,j,k)
              gls_diss=(trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke))
     .                 *gls_c2*tke_epsilon
              gls_Q=gls_diss/trc(i,j,kn,itrgls)
#  endif
#  ifdef GLS
              if (gls_prod+gls_buoy.ge.0.) then
                trc(i,j,kn,itrgls)=
     .            (trc(i,j,kn,itrgls)+delt1*(gls_prod+gls_buoy))
     .            /(1.+delt1*gls_Q)
              else
                trc(i,j,kn,itrgls)=
     .            (trc(i,j,kn,itrgls)+delt1*gls_prod)
     .            /(1.+delt1*(gls_Q-(gls_buoy/trc(i,j,kn,itrgls))))
              endif
              trc(i,j,kn,itrgls)=max(trc(i,j,kn,itrgls),gls_psi_min)
              q=.56**(.5*gls_n)*gls_cmu0**gls_p
     .          *trc(i,j,kn,itrtke)**(gls_m+.5*gls_n)
     .          *bvf(i,k)**(-gls_n)
              if (gls_n.lt.0.) then
                trc(i,j,kn,itrgls)=max(trc(i,j,kn,itrgls),q)
              else
                trc(i,j,kn,itrgls)=min(trc(i,j,kn,itrgls),q)
              endif
#  endif
c
              tke_epsilon=cmu_fac2*trc(i,j,kn,itrtke)**(1.5+gls_m/gls_n)
     .                    *trc(i,j,kn,itrgls)**(-1./gls_n)
              tke_Q=tke_epsilon/trc(i,j,kn,itrtke)
c
              if (tke_prod+tke_buoy.ge.0.) then
                trc(i,j,kn,itrtke)=
     .            (trc(i,j,kn,itrtke)+delt1*(tke_prod+tke_buoy))
     .            /(1.+delt1*tke_Q)
              else
                trc(i,j,kn,itrtke)=
     .            (trc(i,j,kn,itrtke)+delt1*tke_prod)
     .            /(1.+delt1*(tke_Q-(tke_buoy/trc(i,j,kn,itrtke))))
                trc(i,j,kn,itrtke)=max(trc(i,j,kn,itrtke),tke_min)
              endif
c
c --- ------- Penetration of surface TKE below mixed layer.
              if (tkepf.gt.0.) then
                if (dp(i,j,kn).lt.epsil) then
                  q=exp(-p(i,j,k)/tkepls)
                else
                  q=tkepls*(exp(-p(i,j,k  )/tkepls)
     .                     -exp(-p(i,j,k+1)/tkepls))/dp(i,j,kn)
                endif
                trc(i,j,kn,itrtke)=trc(i,j,kn,itrtke)
     .                            +67.83*tkepf*q*ustar(i,j)**2
              endif
c
c --- ------- Set TKE and GLS to prescribed minimum values in surface
c --- ------- mixed layers and thin layers
              if (dp(i,j,kn).lt.epsil) then
                trc(i,j,kn,itrtke)=tke_min
                trc(i,j,kn,itrgls)=gls_psi_min
              endif
              trc(i,j,1+nn,itrtke)=tke_min
              trc(i,j,2+nn,itrtke)=tke_min
              trc(i,j,1+nn,itrgls)=gls_psi_min
              trc(i,j,2+nn,itrgls)=gls_psi_min
c
c --- ------- Bottom Boundary Conditions
              if (k.eq.kmax(i,j)) then
                ust=max(ustarb(i,j),ustmin)
                trc(i,j,kn,itrtke)=max(tke_min,(ust/gls_cmu0)**2)
#  ifdef GLS
                trc(i,j,kn,itrgls)=max(gls_psi_min,
     .                                 (gls_cmu0**(gls_p-2.*gls_m))
     .                                 *(ust**(2.*gls_m))
     .                                 *(kappa*1.e2)**gls_n)
#  endif
              endif
c
              Ls_unlmt=max(Ls_unlmt_min,
     .                     cmu_fac1*trc(i,j,kn,itrgls)**(gls_exp1)
     .                     *trc(i,j,kn,itrtke)**(-tke_exp1))

              if (bvfsq(i,k).gt.0.) then  ! stable stratification
c               Ls_lmt=min(Ls_unlmt,
c    .                     sqrt(.56*trc(i,j,kn,itrtke)
c    .                          /max(bvfsq(i,k),1.e-10)))
 
                Ls_lmt=min(Ls_unlmt,trc(i,j,kn,itrtke)**(-gls_m/gls_n)
     .                    *trc(i,j,kn,itrgls)**gls_n)
c               Ls_lmt=Ls_unlmt
              else                        ! unstable stratification
                Ls_lmt=Ls_unlmt
              endif
c
c --- ------- Compute nondimensional stability functions for tracers
c --- ------- (Sh) and momentum (Sm). Canuto-A
              Gh=min(gls_Gh0,-bvfsq(i,k)*Ls_lmt*Ls_lmt
     .                        /(2.*trc(i,j,kn,itrtke)))
              Gh=min(Gh,(Gh-(Gh-gls_Ghcri)**2)
     .                  /(Gh+gls_Gh0-2.*gls_Ghcri))
              Gh=max(Gh,gls_Ghmin)
              Gh=min(Gh,gls_Gh0)
c
c --- ------- Compute shear number.
              Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh**2))
     .           /(gls_b2-gls_b4*gls_fac6*Gh)
              Gm=min(Gm,Shear2(i,j,k)*Ls_lmt*Ls_lmt
     .                  /(2.*trc(i,j,kn,itrtke)))
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
c --- ------- Relate Canuto stability to BLOM notation
              Sm=Sm*cmu_fac3/gls_cmu0**3
              Sh=Sh*cmu_fac3/gls_cmu0**3
c
              ql=sqrt2*(Ls_lmt)
     .           *sqrt(trc(i,j,kn,itrtke))
c             ql=sqrt2*.5*(Ls_lmt+L_scale(i,j,k))
c    .           *sqrt(trc(i,j,kn,itrtke))
c
c             nus=Sh*ql
c             nus=min(0.1*ql,4.05*nug0)
              nus=min(Sh*ql,4.05*nug0)
c             nus=Sh*(trc(i,j,k,itrtke)*trc(i,j,k,itrtke))
c    .            /trc(i,j,k,itrgls)
              L_scale(i,j,k)=max(Ls_lmt,Ls_unlmt_min)
#  ifdef GLS
c
c --- ------- Recompute gls based on limited length scale
              trc(i,j,kn,itrgls)=
     .          max(gls_cmu0**gls_p*trc(i,j,kn,itrtke)**gls_m
     .              *L_scale(i,j,k)**gls_n,gls_psi_min)
#  endif
#endif
c
c --- ------- Tidally driven diapycnal mixing
              if (tdmflg.eq.1) then
                q=.5*(tanh(4.*(abs(plat(i,j))-tdclat)/tddlat-2.)+1.)
                q=(1.-q)*tdmls0+q*tdmls1
                if (dp(i,j,kn).lt.epsil) then
                  vsf=exp(p(i,j,k)/q)/(q*(exp(p(i,j,kk+1)/q)-1.))
                else
                  vsf=(exp(p(i,j,k+1)/q)-exp(p(i,j,k)/q))
     .                /(dp(i,j,kn)*(exp(p(i,j,kk+1)/q)-1.))
                endif
                nut=g*tdmq*dmxeff*twedon(i,j)*bvfbot(i)*vsf/bvfsq(i,k)
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
              tmp3d_8(i,j,k)=nub
              tmp3d_9(i,j,k)=nus
              tmp3d_10(i,j,k)=nut
              tmp3d_11(i,j,k)=nuls
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
#if defined(TRC) && defined(TKE)
c             trc(i,j,kn,itrtke)=tke_min
c             L_scale(i,j,k)=Ls_unlmt_min
              trc(i,j,kn,itrtke)=trc(i,j,kn-1,itrtke)
              L_scale(i,j,k)=L_scale(i,j,k-1)
#  ifdef GLS
c             trc(i,j,kn,itrgls)=gls_psi_min
              trc(i,j,kn,itrgls)=trc(i,j,kn-1,itrgls)
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
            dfddsu(i)=nu0
            dfddsl(i)=nu0
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
c --- - Diapycnal diffusivity beneath mixed layer by dissipation of
c --- - energy originating from near-inertial waves.
        do k=2,kk-1
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.le.kmax(i,j).and.kmax(i,j)-kfil(i,j).ge.1) then
              q=niwls
              if (k.eq.2.or.dp(i,j,kn).lt.epsil) then
                vsf=exp((p(i,j,3)-p(i,j,k+1))/q)
     .              /(q*(1.-exp((p(i,j,3)-p(i,j,kk+1))/q)))
              else
                vsf=(exp((p(i,j,3)-p(i,j,k  ))/q)
     .              -exp((p(i,j,3)-p(i,j,k+1))/q))
     .              /(dp(i,j,kn)*(1.-exp((p(i,j,3)-p(i,j,kk+1))/q)))
              endif
              nusm=g*niwgf*(1.-niwbf)*niwlf*dmxeff*idkedt(i,j)*vsf
     .             /(alpha0*bvfsq(i,max(k,kfil(i,j))))
              difdia(i,j,k)=difdia(i,j,k)+nusm
#ifdef DIAG
              tmp3d_12(i,j,k)=nusm
#endif
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
      call ncfopn('difest.nc','w','c',1,iotype)
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncdims('z',kdm)
      call ncdefvar('coriop','x y',ndouble,8)
      call ncdefvar('betafp','x y',ndouble,8)
      call ncdefvar('kfpla','x y',ndouble,8)
      call ncdefvar('kfil','x y',ndouble,8)
      call ncdefvar('bcrrd','x y',ndouble,8)
      call ncdefvar('afeql','x y',ndouble,8)
      call ncdefvar('bvfbot','x y',ndouble,8)
      if (edsprs) then
        call ncdefvar('urmse','x y',ndouble,8)
        call ncdefvar('cpse','x y',ndouble,8)
      endif
      call ncdefvar('dp','x y z',ndouble,8)
      call ncdefvar('temp','x y z',ndouble,8)
      call ncdefvar('saln','x y z',ndouble,8)
      call ncdefvar('u','x y z',ndouble,8)
      call ncdefvar('v','x y z',ndouble,8)
      call ncdefvar('bvfsq','x y z',ndouble,8)
      call ncdefvar('rig','x y z',ndouble,8)
      call ncdefvar('egr','x y z',ndouble,8)
      call ncdefvar('rhisc','x y z',ndouble,8)
      call ncdefvar('els','x y z',ndouble,8)
      if (edsprs) then
        call ncdefvar('umnsc','x y z',ndouble,8)
        call ncdefvar('esfac','x y z',ndouble,8)
      endif
      call ncdefvar('difint','x y z',ndouble,8)
      call ncdefvar('difdia','x y z',ndouble,8)
      call ncdefvar('nub','x y z',ndouble,8)
      call ncdefvar('nus','x y z',ndouble,8)
      call ncdefvar('nut','x y z',ndouble,8)
      call ncdefvar('nuls','x y z',ndouble,8)
      call ncdefvar('nusm','x y z',ndouble,8)
#if defined(TRC) && defined(TKE)
      call ncdefvar('tke','x y z',ndouble,8)
#    ifdef GLS
      call ncdefvar('gls_psi','x y z',ndouble,8)
#    endif
#  endif
      call ncedef  
      call ncwrtr('coriop','x y',coriop,ip,1,1.,0.,8)
      call ncwrtr('betafp','x y',betafp,ip,1,1.,0.,8)
      call ncwrtr('kfpla','x y',tmp2d_1,ip,1,1.,0.,8)
      call ncwrtr('kfil','x y',tmp2d_2,ip,1,1.,0.,8)
      call ncwrtr('bcrrd','x y',tmp2d_3,ip,1,1.,0.,8)
      call ncwrtr('afeql','x y',tmp2d_4,ip,1,1.,0.,8)
      call ncwrtr('bvfbot','x y',tmp2d_5,ip,1,1.,0.,8)
      if (edsprs) then
        call ncwrtr('urmse','x y',tmp2d_6,ip,1,1.,0.,8)
        call ncwrtr('cpse','x y',tmp2d_7,ip,1,1.,0.,8)
      endif
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
      if (edsprs) then
        call ncwrtr('umnsc','x y z',tmp3d_6,ip,1,1.,0.,8)
        call ncwrtr('esfac','x y z',tmp3d_7,ip,1,1.,0.,8)
      endif
      call ncwrtr('difint','x y z',difint,ip,1,1.,0.,8)
      call ncwrtr('difdia','x y z',difdia,ip,1,1.,0.,8)
      call ncwrtr('nub','x y z',tmp3d_8,ip,1,1.,0.,8)
      call ncwrtr('nus','x y z',tmp3d_9,ip,1,1.,0.,8)
      call ncwrtr('nut','x y z',tmp3d_10,ip,1,1.,0.,8)
      call ncwrtr('nuls','x y z',tmp3d_11,ip,1,1.,0.,8)
      call ncwrtr('nusm','x y z',tmp3d_12,ip,1,1.,0.,8)
#if defined(TRC) && defined(TKE)
      call ncwrtr('tke','x y z',trc(1-nbdy,1-nbdy,k1n,itrtke),
     .            ip,1,1.,0.,8)
#    ifdef GLS
      call ncwrtr('gls_psi','x y z',trc(1-nbdy,1-nbdy,k1n,itrgls),
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
        call chksummsk(idkedt,ip,1,'idkedt')
        call chksummsk(ustar3,ip,1,'ustar3')
        call chksummsk(difint,ip,kk,'difint')
        call chksummsk(difiso,ip,kk,'difiso')
        call chksummsk(difdia,ip,kk,'difdia')
#if defined(TRC) && defined(TKE)
        call chksummsk(trc(1-nbdy,1-nbdy,1,itrtke),ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(trc(1-nbdy,1-nbdy,1,itrgls),ip,2*kk,'gls_psi')
#  endif
#endif
      endif
c
      return
      end
