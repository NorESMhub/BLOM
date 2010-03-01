      subroutine difest(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate layer interface, isopycnal, and diapycnal diffusivities
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
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
c ---   dpbmin - smallest layer thickness allowed in evaluating
c ---            local gradient richardson number [g/cm/s**2].
c ---   bvfeps - small value of Brunt-Vaisala frequency [1/s].
c ---   thkdff - diffusive velocity for thickness diffusion [cm/s].
c ---   temdff - diffusive velocity for tracer isopycnal diffusion
c ---            [cm/s].
c ---   egc    - the parameter c in the Eden and Greatbatch (2008)
c ---            parameterization [].
c ---   eggam  - the parameter gamma in the Eden and Greatbatch (2008)
c ---            parameterization [].
c ---   egmndf - minimum diffusivity in the Eden and Greatbatch (2008)
c ---            parameterization [cm**2/s].
c ---   egmxdf - maximum diffusivity in the Eden and Greatbatch (2008)
c ---            parameterization [cm**2/s].
c ---   egidfq - factor relating the isopycnal diffusivity to the layer
c ---            interface diffusivity in the Eden and Greatbatch (2008)
c ---            parameterization. egidfq=difint/difiso [].
c ---   ri0    - critical gradient richardson number for shear driven
c ---            vertical mixing [].
c ---   nus0   - maximum shear driven diapycnal diffusivity
c ---            [cm**2/s].
c ---   nug0   - maximum gravity current diapycnal diffusivity
c ---            [cm**2/s].
c ---   bdmc1  - background diapycnal diffusivity times buoyancy
c ---            frequency [cm**2/s**2].
c ---   bdmc2  - background diapycnal diffusivity [cm**2/s].
c ---   tdmq   - tidal dissipation efficiency [].
c ---   tdmgam - tidal driven mixing efficiency [].
c ---   tdmls  - tidal driven mixing length scale [g/cm/s**2].
c ---   cori30 - coriolis parameter at 30N [1/s].
c ---   bvf0   - reference stratification in the parameterization of
c ---            latitude dependent background diapycnal mixing [1/s].
c ---   nubmin - minimum diapycnal diffusivity [cm**2/s].
c ---   nubmax - maximum background diapycnal diffusivity [cm**2/s].
c ---   nusmax - maximum diapycnal diffusivity in the first physical
c ---            layer and the mixed layer [cm**2/s].
c ---   dpgc   - thickness of region near the bottom where the maximum
c ---            diffusivity is increased due to gravity current mixing
c ---            processes [g/cm/s**2].
c ---   dpdfav - thickness of region below the first physical layer used
c ---            to estimate diffusivities in first physical layer and
c ---            the mixed layer [g/cm/s**2].
c ---   dpnbav - thickness of region near the bottom used to estimate
c ---            bottom Brunt-Vaisala frequency [g/cm/s**2].
      integer iidtyp,bdmtyp,bdmldp,tdmflg
      real dptmin,dpbmin,bvfeps,thkdff,temdff,egc,eggam,egmndf,egmxdf,
     .     egidfq,ri0,nus0,nug0,bdmc1,bdmc2,tdmq,tdmgam,tdmls,cori30,
     .     bvf0,nubmin,nubmax,nusmax,dpgc,dpdfav,dpnbav
      parameter (iidtyp=2,bdmtyp=1,bdmldp=1,tdmflg=1,dptmin=98060.,
     .           dpbmin=980.6,bvfeps=1.e-9,thkdff=.5,temdff=.35,egc=1.,
     .           eggam=200.,egmndf=100.e4,egmxdf=1500.e4,egidfq=1.,
     .           ri0=.7,nus0=50.,nug0=2500.,bdmc1=1.2e-3,bdmc2=.1,
     .           tdmq=1./3.,tdmgam=.2,tdmls=500.*98060.,
     .           cori30=7.2722e-5,bvf0=5.24e-3,nubmin=.001,nubmax=1000.,
     .           nusmax=10.,dpgc=300.*98060.,dpdfav=100.*98060.,
     .           dpnbav=250.*98060.)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  dv2
      real, dimension(1-nbdy:idm+nbdy,kdm) ::
     .  du2,bvfsq,bvf,rig
      real, dimension(1-nbdy:idm+nbdy) ::
     .  tup,pup,sup,cr,bcrrd,afeql,bvfbot,dps,dfints,dfisos,dfdias
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  mdv2
      integer, dimension(1-nbdy:idm+nbdy,kdm) ::
     .  mdu2
      integer, dimension(1-nbdy:idm+nbdy) ::
     .  kfpl,klpl
      integer i,j,k,l,kn
      real q,drho,plo,tlo,slo,egr,rhisc,els,nus,nub,nut,tdmf
c
      real rho,sig
      external rho,sig
c
cdiag real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
cdiag.  tmp3d_1,tmp3d_2,tmp3d_3,tmp3d_4,tmp3d_5
cdiag real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
cdiag.  tmp2d_1,tmp2d_2,tmp2d_3,tmp2d_4,tmp2d_5
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
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kfpla(i,j,n)
cdiag     tmp2d_1(i,j)=kfpla(i,j,n)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kfpla(i,j,n)=nint(util1(i,j))
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
        do k=2,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (dpv(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,3,-1
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k.gt.min(kfpla(i,j-1,n),kfpla(i,j,n))+1.and.
     .          dpv(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
        do k=3,kk
          do i=1,ii
            dv2(i,j,k)=0.
            mdv2(i,j,k)=0
          enddo
        enddo
        do k=3,kk
          kn=k+nn
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
c$OMP+ bvfbot,dps,drho,bvfsq,bvf,dfints,dfisos,dfdias,rig,egr,
c$OMP+ rhisc,els,nus,nub,nut,tdmf)
      do j=1,jj
c
c ----- Compute squared vertical velocity gradients of u-component
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          kfpl(i)=kk+1
          klpl(i)=1
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (dpu(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,3,-1
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k.gt.min(kfpla(i-1,j,n),kfpla(i,j,n))+1.and.
     .          dpu(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
        do k=3,kk
          do i=1,ii+1
            du2(i,k)=0.
            mdu2(i,k)=0
          enddo
        enddo
        do k=3,kk
          kn=k+nn
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
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pup(i)=.5*(3.*p(i,j,2)-p(i,j,min(kk,kfpla(i,j,n))+1))
          tup(i)=temp(i,j,k1n)
          sup(i)=saln(i,j,k1n)
          cr(i)=0.
        enddo
        enddo
        do k=2,kk
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
     .             +sqrt(max(0.,rho(p(i,j,k),tlo,slo)
     .                         -rho(p(i,j,k),tup(i),sup(i)))
     .                   *(plo-pup(i)))
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
cdiag     tmp2d_2(i,j)=bcrrd(i)
          afeql(i)=max(abs(coriop(i,j)),sqrt(2.*betafp(i,j)*cr(i)))
cdiag     tmp2d_3(i,j)=afeql(i)
        enddo
        enddo
c
c ----- Locate range of physical interior layers
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          kfpl(i)=kk+1
          klpl(i)=1
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (dp(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,3,-1
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.gt.kfpla(i,j,n)+1.and.dp(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
c
c ----- Compute local gradient richardson number and Brunt-Vaisala
c ----- frequency squared.
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          bvfbot(i)=0.
          dps(i)=0.
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfpl(i).and.k.le.klpl(i).and.
     .          klpl(i)-kfpl(i).ge.1) then
              if     (k.eq.kfpl(i)) then
                q=max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1))
     .                  -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drho=q
                tup(i)=q
              elseif (k.lt.klpl(i)) then
                q=max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1))
     .                  -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drho=2.*tup(i)*q/max(1.e-14,tup(i)+q)
                tup(i)=q
              else
                drho=tup(i)
              endif
c
c --- ------- Brunt-Vaisala frequency squared
              bvfsq(i,k)=max(bvfeps*bvfeps,
     .                       g*g*drho/max(epsil,dp(i,j,kn)))
cdiag         tmp3d_1(i,j,k)=bvfsq(i,k)
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
              rig(i,k)=alpha0*alpha0*drho*dp(i,j,kn)/max(1.e-9,q)
cdiag         tmp3d_2(i,j,k)=rig(i,k)
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
c
c --- - Compute latitude dependent background diffusivity and
c --- - diffusivity weigth to reduce eddy diffusivity when the Rossby
c --- - radius is resolved by the grid
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          q=bcrrd(i)/sqrt(scpx(i,j)*scpx(i,j)+scpy(i,j)*scpy(i,j))
          difwgt(i,j)=1./(1.+q**4)
        enddo
        enddo
c
c --- - Compute layer interface, isopycnal, and diapycnal diffusivities
c       endif
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (iidtyp.eq.1) then
            q=sqrt(scp2(i,j))
            difint(i,j,1)=thkdff*q
            difiso(i,j,1)=temdff*q
          else
            difint(i,j,1)=egmndf
            difiso(i,j,1)=difint(i,j,1)*egidfq
          endif
          difdia(i,j,1)=nubmin
          dfints(i)=0.
          dfisos(i)=0.
          dfdias(i)=0.
          if (dps(i).gt.0.) then
            bvfbot(i)=bvfbot(i)/dps(i)
          endif
          dps(i)=0.
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if(k.ge.kfpl(i).and.k.le.klpl(i).and.
     .         klpl(i)-kfpl(i).ge.1) then
c
c --- ------------------------------------------------------------------
c --- ------- Layer interface and isopycnal diffusivities
c --- ------------------------------------------------------------------
c
              if (iidtyp.eq.1) then
c
c --- --------- Type 1: Diffusivities are diffusive velocities
c --- --------- multiplied by the local horizontal grid scale.
                difint(i,j,k)=difint(i,j,1)
                difiso(i,j,k)=difiso(i,j,1)
c
              else
c
c --- --------- Type 2: Diffusivities are parameterized according
c --- --------- to Eden and Greatbatch (2008).
c
c --- --------- Eady growth rate.
                egr=afeql(i)/sqrt(rig(i,k)+eggam)
cdiag           tmp3d_3(i,j,k)=egr
c
c --- --------- Rhines scale.
                rhisc=egr/max(1.e-24,betafp(i,j))
cdiag           tmp3d_4(i,j,k)=rhisc
c
c --- --------- Eddy length scale.
                els=min(bcrrd(i),rhisc)
cdiag           tmp3d_5(i,j,k)=els
c
c --- --------- Layer interface and isopycnic diffusivity.
                difint(i,j,k)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,egc*egr*els*els*difwgt(i,j)))
                difiso(i,j,k)=difint(i,j,k)*egidfq
c
              endif
c
c --- ------------------------------------------------------------------
c --- ------- Shear driven diapycnal mixing.
c --- ------------------------------------------------------------------
c
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
c
c --- ------------------------------------------------------------------
c --- ------- Background diapycnal mixing.
c --- ------------------------------------------------------------------
c
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
c
c --- ------------------------------------------------------------------
c --- ------- Latitude dependency of background diapycnal mixing
c --- ------------------------------------------------------------------
c
              if (bdmldp.eq.1) then
                q=max(1.e-9,abs(coriop(i,j)))
                nub=nub*q/cori30*log(2.*bvf(i,k)/q)/log(2.*bvf0/cori30)
              endif
c
c --- ------------------------------------------------------------------
c --- ------- Tidally driven diapycnal mixing
c --- ------------------------------------------------------------------
c
              if (tdmflg.eq.1) then
                if (dp(i,j,kn).lt.epsil) then
                  tdmf=exp(p(i,j,k)/tdmls)
     .                 /(tdmls*(exp(p(i,j,kk+1)/tdmls)-1.))
                else
                  tdmf=(exp(p(i,j,k+1)/tdmls)-exp(p(i,j,k)/tdmls))
     .                 /(dp(i,j,kn)*(exp(p(i,j,kk+1)/tdmls)-1.))
                endif
                nut=g*tdmq*tdmgam*twedon(i,j)*bvfbot(i)*tdmf/bvfsq(i,k)
              else
                nut=0.
              endif
c
c --- ------- Total diapycnal diffusivity.
              difdia(i,j,k)=nus+max(nubmin,min(nubmax,nub+nut))
c
c --- ------- Accumulate diffusivities in a region below the first
c --- ------- physical layer
              q=max(0.,min(p(i,j,kfpl(i))+dpdfav,p(i,j,k+1))-p(i,j,k))
              dps(i)=dps(i)+q
              dfints(i)=dfints(i)+difint(i,j,k)*q
              dfisos(i)=dfisos(i)+difiso(i,j,k)*q
              dfdias(i)=dfdias(i)+max(nubmin,min(nubmax,nub))*q
c
            else
              difint(i,j,k)=difint(i,j,k-1)
              difiso(i,j,k)=difiso(i,j,k-1)
              difdia(i,j,k)=difdia(i,j,k-1)
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i).gt.0.) then
            q=1./dps(i)
            dfints(i)=dfints(i)*q
            dfisos(i)=dfisos(i)*q
            dfdias(i)=min(nusmax,dfdias(i)*q)
          else
            if (iidtyp.eq.1) then
              q=sqrt(scp2(i,j))
              dfints(i)=thkdff*q
              dfisos(i)=temdff*q
            else
              dfints(i)=egmndf
              dfisos(i)=dfints(i)*egidfq
            endif
            dfdias(i)=nubmin
          endif
        enddo
        enddo
        do k=kk-1,1,-1
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.lt.kfpl(i)) then
              difint(i,j,k)=dfints(i)
              difiso(i,j,k)=dfisos(i)
              difdia(i,j,k)=dfdias(i)
            endif
          enddo
          enddo
        enddo
c
      enddo
c$OMP END PARALLEL DO
c
cdiag call ncfopn('difest.nc','w')
cdiag call ncdims('x',itdm)
cdiag call ncdims('y',jtdm)
cdiag call ncdims('z',kdm)
cdiag call ncwrtr('coriop','x y',coriop,ip,1)
cdiag call ncwrtr('betafp','x y',betafp,ip,1)
cdiag call ncwrtr('kfpla','x y',tmp2d_1,ip,1)
cdiag call ncwrtr('bcrrd','x y',tmp2d_2,ip,1)
cdiag call ncwrtr('afeql','x y',tmp2d_3,ip,1)
cdiag call ncwrtr('dp','x y z',dp(1-nbdy,1-nbdy,k1n),ip,1)
cdiag call ncwrtr('temp','x y z',temp(1-nbdy,1-nbdy,k1n),ip,1)
cdiag call ncwrtr('saln','x y z',saln(1-nbdy,1-nbdy,k1n),ip,1)
cdiag call ncwrtr('u','x y z',u(1-nbdy,1-nbdy,k1n),iu,1)
cdiag call ncwrtr('v','x y z',v(1-nbdy,1-nbdy,k1n),iv,1)
cdiag call ncwrtr('bvfsq','x y z',tmp3d_1,ip,1)
cdiag call ncwrtr('rig','x y z',tmp3d_2,ip,1)
cdiag call ncwrtr('egr','x y z',tmp3d_3,ip,1)
cdiag call ncwrtr('rhisc','x y z',tmp3d_4,ip,1)
cdiag call ncwrtr('els','x y z',tmp3d_5,ip,1)
cdiag call ncwrtr('difint','x y z',difint,ip,1)
cdiag call ncwrtr('difdia','x y z',difdia,ip,1)
cdiag call ncfcls
cdiag call xcstop('(difest)')
cdiag        stop '(difest)'
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest:'
        endif
        call chksummsk(difint,ip,kk,'difint')
        call chksummsk(difiso,ip,kk,'difiso')
        call chksummsk(difdia,ip,kk,'difdia')
      endif
c
      return
      end
