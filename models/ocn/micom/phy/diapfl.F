      subroutine diapfl(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- diapycnal mixing
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
c --- parameters:
c ---   teps   - small temperature difference used in removing
c ---            instabilitites [K]
c ---   small  - small relative density increment []
c ---   restmx - if true, restore mixed layer base to its initial depth
c ---   dprmin - smallest layer thickness allowed in evaluating
c ---            local gradient richardson number [g/cm/s**2]
c ---   thrgbt - thin region near bottom that will be mixed prior to
c ---            more physical based mixing [g/cm/s**2]
c ---   dpdmin - layer thickness limit for what is considered a thin
c ---            interior layer [g/cm/s**2]
c ---   ri0    - critical gradient richardson number for shear
c ---            instability and gravity current mixing []
c ---   nus0   - maximum shear instability diffusivity [cm^2/s]
c ---   nug0   - maximum gravity current diffusivity [cm^2/s]
c ---   nuls0  - maximum diffusivity to avoid local weak stability
c ---            [cm^2/s]
c ---   dpgc   - thickness of region near the bottom where the maximum
c ---            diffusivity is increased due to gravity current mixing
c ---            processes [g/cm/s**2]
c ---   dsglc0 - critical local interface density difference [g/cm**3]
c ---   maxflx - maximum interface flux [g/cm/s**2]
c ---   dsig0  - apply a strong density restoration if the local density
c ---            error is greater than dsig0 [g/cm**3]
      real teps,small
      logical restmx
      real dprmin,thrgbt,dpdmin,ri0,nus0,nug0,nuls0,dpgc,dsglc0,maxflx,
     .     dsig0
      parameter (teps=1.e-9,small=1.e-22,restmx=.false.,
     .           dprmin=980.6,thrgbt=2.*98.06,dpdmin=98060.,ri0=1.,
     .           nus0=50.,nug0=500.,nuls0=500.,dpgc=300.*98060.,
     .           dsglc0=1.e-5,maxflx=5.*98060.,dsig0=1.e-6)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: rig
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  du2u,du2l,dv2u,dv2l
      real, dimension(0:kdm) :: flxl,flxnzl
      real, dimension(kdm+1) :: flxu,flxnzu,pn,po
      real, dimension(kdm) ::
     .  ttem,ssal,delp,dens,sigr,ri,dsglc,ttem_t,ssal_t,delp_t,dens_t,
     .  sigr_t,ri_t,dsglc_t,pres,dsgiu,dsgil,delpn,aa,cc,dd,gam,uo,un
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  mskuu,mskul,mskvu,mskvl
      integer, dimension(kdm) :: ktil,ktilu,ktill
      logical, dimension(kdm) :: istil
      real esci,dpsum,tdp,sdp,dsig,alfa,beta,dentrg,tnew,snew,z,q,
     .     a,b,c,d,e,f,delpml,ennsq,nus,flngth,flxmxu,flxmxl,
     .     flxdru,flxdrl,bei,udpn,rmax_niter
      integer i,j,k,l,kn,msum,kmin,kmax,ku,kl,kt,ntil,
     .        ntl,ko,niter
      logical changed,done
c
      real dsigdt,dsigds,sig,sofsig,tsroot,rho
      external dsigdt,dsigds,sig,sofsig,tsroot,rho
c
      real tsum0,tsum1,ssum0,ssum1
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc,ttrc_t
      real, dimension(ntr) :: trdp
      integer nt
#endif
c
      rmax_niter=0.
c
c --- inverse of e-folding time scale (in time steps) for restoring
c --- target density
      esci=delt1/86400.
c
c --- ------------------------------------------------------------------
c --- compute local gradient richardson number
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do j=1,jj+1
        do i=1,ii+1
          rig(i,j,1)=-99.
          du2u(i,j)=0.
          du2l(i,j)=0.
          dv2u(i,j)=0.
          dv2l(i,j)=0.
          mskuu(i,j)=0
          mskul(i,j)=0
          mskvu(i,j)=0
          mskvl(i,j)=0
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=klist(i,j)
        enddo
        enddo
      enddo
      call xctilr(util1, 1,1, 1,1, halo_ps)
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          klist(i,j)=nint(util1(i,j))
        enddo
        enddo
      enddo
      call xctilr(u(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
      call xctilr(v(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)
c
      do k=2,kk
c
c --- - find vertical velocity gradients across interfaces, only taking
c --- - into account layers thicker than a treshold thickness
        kn=k+nn
c
c$OMP PARALLEL DO PRIVATE(kmin,kt,q)
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            kmin=min(klist(i-1,j),klist(i,j))
            if (k.lt.kmin.or.dpu(i,j,kn).lt.dprmin) then
              du2u(i,j)=0.
              du2l(i,j)=0.
              mskuu(i,j)=0
              mskul(i,j)=0
            else
              kt=k-1
              do while (dpu(i,j,kt+nn).lt.dprmin.and.kt.ge.kmin)
                kt=kt-1
              enddo
              if (kt.lt.kmin) then
                du2u(i,j)=0.
                mskuu(i,j)=0
              else
                q=u(i,j,kn)-u(i,j,kt+nn)
                du2u(i,j)=q*q
                mskuu(i,j)=1
              endif
              kt=k+1
              do while (dpu(i,j,min(kk,kt)+nn).lt.dprmin.and.kt.le.kk)
                kt=kt+1
              enddo
              if (kt.gt.kk) then
                du2l(i,j)=0.
                mskul(i,j)=0
              else
                q=u(i,j,kn)-u(i,j,kt+nn)
                du2l(i,j)=q*q
                mskul(i,j)=1
              endif
            endif
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kmin,kt,q)
        do j=1,jj+1
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            kmin=min(klist(i,j-1),klist(i,j))
            if (k.lt.kmin.or.dpv(i,j,kn).lt.dprmin) then
              dv2u(i,j)=0.
              dv2l(i,j)=0.
              mskvu(i,j)=0
              mskvl(i,j)=0
            else
              kt=k-1
              do while (dpv(i,j,kt+nn).lt.dprmin.and.kt.ge.kmin)
                kt=kt-1
              enddo
              if (kt.lt.kmin) then
                dv2u(i,j)=0.
                mskvu(i,j)=0
              else
                q=v(i,j,kn)-v(i,j,kt+nn)
                dv2u(i,j)=q*q
                mskvu(i,j)=1
              endif
              kt=k+1
              do while (dpv(i,j,min(kk,kt)+nn).lt.dprmin.and.kt.le.kk)
                kt=kt+1
              enddo
              if (kt.gt.kk) then
                dv2l(i,j)=0.
                mskvl(i,j)=0
              else
                q=v(i,j,kn)-v(i,j,kt+nn)
                dv2l(i,j)=q*q
                mskvl(i,j)=1
              endif
            endif
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - compute local gradient richardson number
c$OMP PARALLEL DO PRIVATE(msum,dsig)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            rig(i,j,k)=-99.
            if (k.ge.klist(i,j)) then
              msum=mskuu(i  ,j)+mskul(i  ,j)
     .            +mskuu(i+1,j)+mskul(i+1,j)
     .            +mskvu(i,j  )+mskvl(i,j  )
     .            +mskvu(i,j+1)+mskvl(i,j+1)
              if (msum.gt.0) then
                if     (k.eq.klist(i,j)) then
                  if (k.eq.kk) then
                    dsig=sigma(i,j,kn)-sigma(i,j,k1n)
                  else
                    dsig=.5*(sigma(i,j,kn+1)-sigma(i,j,k1n))
                  endif
                elseif (k.eq.kk) then
                  dsig=sigma(i,j,kn)-sigma(i,j,kn-1)
                else
                  dsig=.5*(sigma(i,j,kn+1)-sigma(i,j,kn-1))
                endif
                rig(i,j,k)=alpha0*alpha0*dsig*dp(i,j,kn)*msum
     .                     /max(1.e-9,du2u(i  ,j)+du2l(i  ,j)
     .                               +du2u(i+1,j)+du2l(i+1,j)
     .                               +dv2u(i,j  )+dv2l(i,j  )
     .                               +dv2u(i,j+1)+dv2l(i,j+1))
              endif
            endif
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c --- extrapolate the richardson number to points with no value
c$OMP PARALLEL DO PRIVATE(kmin,kmax,ku,kt)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          kmin=klist(i,j)
          do while (rig(i,j,min(kk,kmin)).eq.-99..and.kmin.le.kk)
            kmin=kmin+1
          enddo
          if (kmin.le.kk) then
            do k=klist(i,j)+1,kmin-1
              rig(i,j,k)=rig(i,j,kmin)
            enddo
            kmax=kk
            do while (rig(i,j,kmax).eq.-99.)
              kmax=kmax-1
            enddo
            do k=kmax+1,kk
              rig(i,j,k)=rig(i,j,kmax)
            enddo
            ku=kmin
            do k=kmin+1,kmax-1
              if (rig(i,j,k).eq.-99.) then
                kt=k+1
                do while (rig(i,j,kt).eq.-99.)
                  kt=kt+1
                enddo
                rig(i,j,k)=min(rig(i,j,ku),rig(i,j,kt))
              else
                ku=k
              endif
            enddo
          endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ttem,ssal,delp,dens,sigr,ri,dsglc,ttem_t,
c$OMP+ ssal_t,delp_t,dens_t,sigr_t,ri_t,dsglc_t,pres,dsgiu,dsgil,flxu,
c$OMP+ flxl,flxnzu,flxnzl,delpn,aa,cc,dd,gam,dpsum,tdp,q,sdp,dsig,alfa,
c$OMP+ beta,dentrg,tnew,snew,z,a,b,c,d,e,f,delpml,ennsq,nus,flngth,
c$OMP+ flxmxu,flxmxl,flxdru,flxdrl,bei,ktil,ktilu,ktill,kn,kmin,kmax,
c$OMP+ ku,kl,kt,ntil,ntl,niter,istil,changed,done
#ifdef TRC
c$OMP+,ttrc,ttrc_t,trdp
#endif
c$OMP+ )
c$OMP+ REDUCTION(max:rmax_niter)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- copy variables into 1d arrays
          do k=1,kk
            kn=k+nn
            ttem(k)=temp(i,j,kn)
            ssal(k)=saln(i,j,kn)
            delp(k)=dp(i,j,kn)
            dens(k)=sigma(i,j,kn)
            sigr(k)=sigmar(k)
            ri(k)=rig(i,j,k)
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,k)=trc(i,j,k,nt)
            enddo
#endif
          enddo
c
c --- --- to avoid layers with near zero thickness near the bottom, mix
c --- --- layers within a thin region near the bottom into first
c --- --- physical layer
          dpsum=0.
          tdp=0.
          sdp=0.
#ifdef TRC
          do nt=1,ntrc  
            trdp(nt)=0.
          enddo
#endif
          k=kk+1
          do while (dpsum.lt.thrgbt)
            k=k-1
            dpsum=dpsum+delp(k)
            tdp=tdp+ttem(k)*delp(k)
            sdp=sdp+ssal(k)*delp(k)
#ifdef TRC
            do nt=1,ntrc  
              trdp(nt)=trdp(nt)+ttrc(nt,k)*delp(k)
            enddo
#endif
            delp(k)=0.
          enddo
          delp(k)=dpsum
          q=1./dpsum
          ttem(k)=tdp*q
          ssal(k)=sdp*q
#ifdef TRC
          do nt=1,ntrc  
            ttrc(nt,k)=trdp(nt)*q
          enddo
#endif
c
c --- --- locate range of physical layers
          kmin=klist(i,j)-1
          kmax=1
          do k=2,kk
            if (delp(k).gt.epsil) kmax=k
          enddo
c
          ttem(kmin)=ttem(1)
          ssal(kmin)=ssal(1)
          delp(kmin)=delp(1)
          dens(kmin)=dens(1)
#ifdef TRC
          do nt=1,ntrc  
            ttrc(nt,kmin)=ttrc(nt,1)
          enddo
#endif
c
c --- ------------------------------------------------------------------
c --- --- remove instability due to nonlinear equation of state
c --- ------------------------------------------------------------------
c
          changed=.true.
          niter=0
c
          do while (changed.and.niter.lt.10)
c
            changed=.false.
            niter=niter+1
c
            k=kmin
c
            do while (k.lt.kmax)
c
              k=k+1
c
              if (dens(k).lt.dens(k-1)) then
                write (lp,*) '(dens(k).lt.dens(k-1))',
     .                       i0+i,j0+j,k,kmin,kmax
              endif
              alfa=-alpha0*dsigdt(ttem(k),ssal(k))
              beta= alpha0*dsigds(ttem(k),ssal(k))
c
              if (beta*(ssal(k)-ssal(k-1))
     .           -alfa*(ttem(k)-ttem(k-1)).lt.0.) then
c
                changed=.true.
c
                done=.false.
                ku=k
c
                do while (.not.done)
                  done=.true.
                  ku=ku-1
c
                  dentrg=min(dens(k),max(sigr(k),dens(ku)+epsil))
                  if (abs(ttem(k)-ttem(ku)).lt.teps) then
                    tnew=ttem(k)
                    snew=sofsig(dentrg,tnew)
                    z=max(0.,(snew-ssal(k))*delp(k)/(ssal(ku)-snew))
                  else
                    a=ttem(ku)-ttem(k)
                    b=0.
                    c=0.
                    d=a
                    e=ssal(ku)-ssal(k)
                    f=ssal(k)*ttem(ku)-ssal(ku)*ttem(k)
                    tnew=tsroot(a,b,c,d,e,f,dentrg,ttem(ku),-999.)
                    if (abs(ttem(k)-tnew).ge.abs(ttem(k)-ttem(ku))) then
                      write (lp,*) 'strange temp',i0+i,j0+j,k,ku,
     .                             ttem(k),ttem(ku),tnew
                      if (delp(k).lt.epsil.or.
     .                    delp(k)+delp(ku).lt.epsil) then
                        z=0.
                        tnew=ttem(ku)
                        snew=sofsig(dentrg,tnew)
                      elseif (delp(ku).lt.epsil) then
                        z=epsil
                      else
                        z=.1*delp(ku)
                        tnew=(ttem(k)*delp(k)+ttem(ku)*z)/(delp(k)+z)
                        snew=(ssal(k)*delp(k)+ssal(ku)*z)/(delp(k)+z)
                      endif
                    else
                      snew=(ssal(k)*(ttem(ku)-tnew)
     .                     +ssal(ku)*(tnew-ttem(k)))
     .                     /(ttem(ku)-ttem(k))
                      z=max(0.,(tnew-ttem(k))*delp(k)/(ttem(ku)-tnew))
                    endif
                  endif
c
                  if (ku.eq.kmin.and.z.gt.delp(ku)-thkmin*onem) then
c
                    q=1./(delp(k)+delp(ku))
                    ttem(ku)=(ttem(k)*delp(k)+ttem(ku)*delp(ku))*q
                    ssal(ku)=(ssal(k)*delp(k)+ssal(ku)*delp(ku))*q
                    dens(ku)=sig(ttem(ku),ssal(ku))
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,ku)=(ttrc(nt,k )*delp(k )
     .                            +ttrc(nt,ku)*delp(ku))*q
                    enddo
#endif
                    delp(ku)=delp(k)+delp(ku)
                    delp(k)=0.
                    do while (dens(ku).gt.dens(min(kmax,kmin+1)).and.
     .                        kmin.lt.kmax)
                      kmin=kmin+1
                      if (kmin.gt.k) then
                        q=1./(delp(kmin)+delp(ku))
                        ttem(ku)=(ttem(kmin)*delp(kmin)
     .                           +ttem(ku)*delp(ku))*q
                        ssal(ku)=(ssal(kmin)*delp(kmin)
     .                           +ssal(ku)*delp(ku))*q
                        dens(ku)=sig(ttem(ku),ssal(ku))
#ifdef TRC
                        do nt=1,ntrc
                          ttrc(nt,ku)=(ttrc(nt,kmin)*delp(kmin)
     .                                +ttrc(nt,ku)*delp(ku))*q
                        enddo
#endif
                        delp(ku)=delp(kmin)+delp(ku)
                        delp(kmin)=0.
                        k=kmin
                      endif
                    enddo
                    if (kmin.gt.ku) then
                      ttem(kmin)=ttem(ku)
                      ssal(kmin)=ssal(ku)
                      dens(kmin)=dens(ku)
#ifdef TRC
                      do nt=1,ntrc
                        ttrc(nt,kmin)=ttrc(nt,ku)
                      enddo
#endif
                      delp(kmin)=delp(ku)
                    endif
                    ku=kmin
                    ttem(k)=ttem(ku)
                    ssal(k)=sofsig(dens(k),ttem(k))
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=ttrc(nt,ku)
                    enddo
#endif
c
                  elseif (z.gt.delp(ku)) then
c
                    q=1./(delp(k)+delp(ku))
                    ttem(k)=(ttem(k)*delp(k)+ttem(ku)*delp(ku))*q
                    ssal(k)=(ssal(k)*delp(k)+ssal(ku)*delp(ku))*q
#ifdef TRC
                    do nt=1,ntrc
                      ttrc(nt,k)=(ttrc(nt,k )*delp(k )
     .                           +ttrc(nt,ku)*delp(ku))*q
                    enddo
#endif
                    delp(k)=delp(k)+delp(ku)
                    delp(ku)=0.
                    done=.false.
c
                  else
c
                    ttem(k)=tnew
                    ssal(k)=snew
                    dens(k)=sig(ttem(k),ssal(k))
                    if (delp(k)+z.gt.epsil) then
                      q=1./(delp(k)+z)
#ifdef TRC
                      do nt=1,ntrc  
                        ttrc(nt,k)=(ttrc(nt,k)*delp(k)+ttrc(nt,ku)*z)*q
                      enddo
#endif
                    endif
                    delp(k)=delp(k)+z
                    delp(ku)=delp(ku)-z
c
                  endif
c
                enddo
c
                do kt=ku+1,k-1
                  q=real(kt-k)/real(ku-k)
                  ttem(kt)=ttem(k)+q*(ttem(ku)-ttem(k))
                  ssal(kt)=ssal(k)+q*(ssal(ku)-ssal(k))
                  dens(kt)=sig(ttem(kt),ssal(kt))
#ifdef TRC
                  do nt=1,ntrc  
                    ttrc(nt,kt)=ttrc(nt,k)+q*(ttrc(nt,ku)-ttrc(nt,k))
                  enddo
#endif
                enddo
c
              endif
            enddo
          enddo
c
          rmax_niter=max(rmax_niter,real(niter))
c
c --- ------------------------------------------------------------------
c --- --- diapycnal mixing
c --- ------------------------------------------------------------------
c
          if (kmin.lt.kmax) then
c
            delpml=delp(kmin)
c
c --- ----- save lower interface pressure
            pres(kmin)=delp(kmin)
            do k=kmin+1,kmax
              pres(k)=pres(k-1)+delp(k)
            enddo
c
c --- ----- compute for each layer the minimum local density difference
c --- ----- across the upper and lower layer interfaces
            do k=kmin+1,kmax
              dsglc(k)=rho(pres(k),ttem(k  ),ssal(k  ))
     .                -rho(pres(k),ttem(k-1),ssal(k-1))
            enddo
            do k=kmin+1,kmax-1
              dsglc(k)=min(dsglc(k),dsglc(k+1))
            enddo
c
c --- ----- exclude thin interior layers in the computation of regular
c --- ----- diapycnal fluxes
c           tsum0=0.
c           ssum0=0.
c           do k=kmin,kmax
c             tsum0=tsum0+ttem(k)*delp(k)
c             ssum0=ssum0+ssal(k)*delp(k)
c           enddo
c           if (ntrc.gt.0) then
c             tsum0=0.
c             do k=kmin,kmax
c               tsum0=tsum0+ttrc(1,k)*delp(k)
c             enddo
c           endif
c
c           if (i.eq.itest.and.j.eq.jtest) then
c             write (*,*) 'kmin,kmax',kmin,kmax
c             do k=kmin,kmax
c               write (*,'(i2,3f12.6,3e12.6)')
c    .            k,delp(k)/onem,ttem(k),ssal(k),ri(k),dens(k),
c    .            dens(k)-sigr(k)
c             enddo
c             write (*,*) 'tsum0,ssum0',
c    .                    tsum0/p(i,j,kk+1),ssum0/p(i,j,kk+1)
c             write (*,*)
c           endif
c
            kl=kmax
            do while (delp(kl).le.dpdmin)
              istil(kl)=.false.
              kl=kl-1
            enddo
            istil(kl)=.false.
            ntil=0
            do k=kl-1,kmin,-1
              if (delp(k).le.dpdmin) then
                ntil=ntil+1
                istil(k)=.true.
                ktil(ntil)=k
                ttem_t(ntil)=ttem(k)
                ssal_t(ntil)=ssal(k)
#ifdef TRC
                do nt=1,ntrc  
                  ttrc_t(nt,ntil)=ttrc(nt,k)
                enddo
#endif
                delp_t(ntil)=delp(k)
                dens_t(ntil)=dens(k)
                sigr_t(ntil)=sigr(k)
                ri_t(ntil)=ri(k)
                dsglc_t(ntil)=dsglc(k)
              else
                istil(k)=.false.
              endif
            enddo
            if (ntil.gt.0) then
              ku=kl
              do k=kl-1,kmin,-1
                if (delp(k).gt.dpdmin) then
                  ku=ku-1
                  ttem(ku)=ttem(k)
                  ssal(ku)=ssal(k)
#ifdef TRC
                  do nt=1,ntrc  
                    ttrc(nt,ku)=ttrc(nt,k)
                  enddo
#endif
                  delp(ku)=delp(k)
                  dens(ku)=dens(k)
                  sigr(ku)=sigr(k)
                  pres(ku)=pres(k)
                  ri(ku)=ri(k)
                  dsglc(ku)=dsglc(k)
                endif
              enddo
              kmin=ku
            endif
c
c           if (i.eq.itest.and.j.eq.jtest) then
c             write (*,*) 'kmin,kmax,ntil',kmin,kmax,ntil
c             do k=kmin,kmax
c               write (*,'(i2,3f12.6,3e12.6)')
c    .            k,delp(k)/onem,ttem(k),ssal(k),ri(k),dens(k),
c    .            dens(k)-sigr(k)
c             enddo
c             write (*,*)
c           endif
c
c --- ----- compute interface mass fluxes due to vertical turbulent
c --- ----- diffusion
c
            flxl(kmin-1)=0.
            flxu(kmin)=0.
            flxl(kmin)=0.
            flxu(kmax)=0.
            flxl(kmax)=0.
            flxu(kmax+1)=0.
            flxnzu(kmin)=0.
            flxnzl(kmin)=0.
            flxnzu(kmax)=0.
            flxnzl(kmax)=0.
c
            k=kmin+1
            do while (k.lt.kmax)
c
c --- ------- ennsq = buoy.freq.^2
              ennsq=.5*g*g*max(0.,dens(k+1)-dens(k-1))
     .              /max(delp(k),onecm)
c
              if (ri(k).lt.ri0) then
c
c --- --------- maximum diffusivity is increased near the bottom to
c --- --------- provide additional mixing of gravity currents
                q=(pres(kmax)-pres(k)+.5*delp(k))
     .            /min(dpgc,.5*pres(kmax))
                q=max(0.,1.-q*q)
                q=q*q*q
                nus=q*nug0+(1.-q)*nus0
c
c --- --------- parameterization of diffusivity as a function of local
c --- --------- gradient richardson number
                q=ri(k)/ri0
                q=max(0.,1.-q*q)
                nus=nus*q*q*q
              else
                nus=0.
              endif
c
c --- ------- enhance diffusivity when local stability is weak
              if (dsglc(k).lt.dsglc0) then
                q=max(0.,dsglc(k))/dsglc0
                q=max(0.,1.-q*q)
                nus=nus+nuls0*q*q*q
              endif
c
c --- ------- store (exch.coeff x buoy.freq.^2 / g x time step) in
c --- ------- -flngth- (dimensions of flngth: length in pressure units)
c --- 
c --- ------- use the following if exch.coeff. = diapyc / buoyancy
c --- ------- frequency
              flngth=(diapyc*sqrt(ennsq)+nus*ennsq)*delt1/alpha0
c --- 
c --- ------- use the following if exch.coeff. = diapyc
ccc           flngth=(diapyc+nus)*ennsq*delt1/alpha0
c
              alfa=-alpha0*dsigdt(ttem(k),ssal(k))
              beta= alpha0*dsigds(ttem(k),ssal(k))
c
              dsgiu(k)=1./max(beta*(ssal(k)-ssal(k-1))
     .                       -alfa*(ttem(k)-ttem(k-1)),small)
              dsgil(k)=1./max(beta*(ssal(k+1)-ssal(k))
     .                       -alfa*(ttem(k+1)-ttem(k)),small)
c
              flxmxu=flngth*dsgiu(k)
              flxmxl=flngth*dsgil(k)
              q=min(1.,maxflx/max(flxmxu,flxmxl),
     .              .25*delp(k-1)/flxmxu,.25*delp(k+1)/flxmxl)
              flxu(k)=q*flxmxu
              flxl(k)=q*flxmxl
              if (max(flxu(k),flxl(k)).lt.epsil) then
                q=1./max(flxmxu,flxmxl)
                flxnzu(k)=(flxmxu*q)*epsil
                flxnzl(k)=(flxmxl*q)*epsil
              else
                flxnzu(k)=flxu(k)
                flxnzl(k)=flxl(k)
              endif
              k=k+1
c
            enddo
c
c --- ----- modify mass fluxes to restore layer density
c
            do k=kmin,kmax
              delpn(k)=max(0.,delp(k)
     .                       +flxu(k)+flxl(k)-flxu(k+1)-flxl(k-1))
            enddo
c
            do k=kmin+1,kmax-1
c
c --- ------- store mixing length required for density restoration in
c --- ------- -flngth- (dimensions of flngth: length in pressure units)
              dsig=dens(k)-sigr(k)
              q=dsig/dsig0
              q=max(0.,1.-q*q)
              q=q*q*q
              flngth=(1.+(esci-1.)*q)*dsig*delp(k)
c
              if     (flngth.gt.0.) then
                flxdru=min(.1*delpn(k),flngth*dsgiu(k))
                q=min(1.,.25*delpn(k-1)/flxdru)
                if (q.lt.1.) then
                  flxu(k)=flxu(k)+q*flxdru
                  flxnzu(k)=flxnzu(k)+q*flxdru
                  flxdrl=max(-.25*delpn(k),(q-1.)*flngth*dsgil(k))
                  flxl(k)=flxl(k)+flxdrl
                  flxnzl(k)=flxnzl(k)+flxdrl
                else
                  flxu(k)=flxu(k)+flxdru
                  flxnzu(k)=flxnzu(k)+flxdru
                endif
              elseif (flngth.lt.0.) then
                flxdrl=min(.1*delpn(k),-flngth*dsgil(k))
                q=min(1.,.25*delpn(k+1)/flxdrl)
                if (q.lt.1.) then
                  flxl(k)=flxl(k)+q*flxdrl
                  flxnzl(k)=flxnzl(k)+q*flxdrl
                  flxdru=max(-.25*delpn(k),(1.-q)*flngth*dsgiu(k))
                  flxu(k)=flxu(k)+flxdru
                  flxnzu(k)=flxnzu(k)+flxdru
                else
                  flxl(k)=flxl(k)+flxdrl
                  flxnzl(k)=flxnzl(k)+flxdrl
                endif
              endif
c
            enddo
c
            k=kmax
            dsig=dens(k)-sigr(k)
            q=dsig/dsig0
            q=max(0.,1.-q*q)
            q=q*q*q
            flngth=(1.+(esci-1.)*q)*dsig*delp(k)
            alfa=-alpha0*dsigdt(ttem(k),ssal(k))
            beta= alpha0*dsigds(ttem(k),ssal(k))
            dsgiu(k)=1./max(beta*(ssal(k)-ssal(k-1))
     .                     -alfa*(ttem(k)-ttem(k-1)),small)
            flxu(k)=min(.1*delpn(k),.25*delpn(k-1),max(-.1*delpn(k),
     .                  flngth*dsgiu(k)))
            flxnzu(k)=flxu(k)
c
c --- ----- obtain coefficients for tridiagonal system for salinity and
c --- ----- temperature
            do k=kmin,kmax
              q=1./max(epsil,delp(k)+flxnzu(k)+flxnzl(k))
              aa(k)=-flxnzu(k)*q
              cc(k)=-flxnzl(k)*q
              dd(k)=delp(k)*q
            enddo
c
c --- ----- solve the tridiagonal system for temperature and salinity
            bei=1.
            do k=kmin+1,kmax
              gam(k)=cc(k-1)*bei
              bei=1./(1.-aa(k)*gam(k))
              ssal(k)=(dd(k)*ssal(k)-aa(k)*ssal(k-1))*bei
              ttem(k)=(dd(k)*ttem(k)-aa(k)*ttem(k-1))*bei
            enddo
            do k=kmax-1,kmin+1,-1
              ssal(k)=ssal(k)-gam(k+1)*ssal(k+1)
              ttem(k)=ttem(k)-gam(k+1)*ttem(k+1)
            enddo
c
#ifdef TRC
c --- ----- to ensure monotonicity of tracers, enforce positive
c --- ----- interface fluxes
            do k=kmin+1,kmax
              q=min(0.,flxnzl(k-1))
              flxnzl(k-1)=flxnzl(k-1)-q
              flxnzu(k  )=flxnzu(k  )-q
              q=min(0.,flxnzu(k  ))
              flxnzl(k-1)=flxnzl(k-1)-q
              flxnzu(k  )=flxnzu(k  )-q
            enddo
            do k=kmin,kmax
              q=1./max(epsil,delp(k)+flxnzu(k)+flxnzl(k))
              aa(k)=-flxnzu(k)*q
              cc(k)=-flxnzl(k)*q
              dd(k)=delp(k)*q
            enddo
c
c --- ----- solve the tridiagonal system for tracers
            k=kmin
            do nt=1,ntrc
              ttrc(nt,k)=dd(k)*ttrc(nt,k)
            enddo
            bei=1.
            do k=kmin+1,kmax
              gam(k)=cc(k-1)*bei
              bei=1./(1.-aa(k)*gam(k))
              do nt=1,ntrc
                ttrc(nt,k)=(dd(k)*ttrc(nt,k)-aa(k)*ttrc(nt,k-1))*bei
              enddo
            enddo
            do k=kmax-1,kmin,-1
              do nt=1,ntrc  
                ttrc(nt,k)=ttrc(nt,k)-gam(k+1)*ttrc(nt,k+1)
              enddo
            enddo
#endif
c
c --- ----- update layer thicknesses and density
            do k=kmin,kmax
              delp(k)=max(0.,delp(k)
     .                      +flxu(k)+flxl(k)-flxu(k+1)-flxl(k-1))
              dens(k)=sig(ttem(k),ssal(k))
            enddo
c
            if (ntil.gt.0) then
c
c --- ------- put thin interior layers back into the arrays and inflate
c --- ------- the thin layers by entraining water from the nearest thick
c --- ------- layers
c
              kmin=kmin-ntil
              ku=kmin
              do ntl=ntil,1,-1
                do k=ku,ktil(ntl)-1
                  kl=k+ntl
                  ttem(k)=ttem(kl)
                  ssal(k)=ssal(kl)
#ifdef TRC
                  do nt=1,ntrc  
                    ttrc(nt,k)=ttrc(nt,kl)
                  enddo
#endif
                  delp(k)=delp(kl)
                  dens(k)=dens(kl)
                  sigr(k)=sigr(kl)
                  ri(k)=ri(kl)
                  dsglc(k)=dsglc(kl)
                enddo
                k=ktil(ntl)
                ttem(k)=ttem_t(ntl)
                ssal(k)=ssal_t(ntl)
#ifdef TRC
                do nt=1,ntrc  
                  ttrc(nt,k)=ttrc_t(nt,ntl)
                enddo
#endif
                delp(k)=delp_t(ntl)
                dens(k)=dens_t(ntl)
                sigr(k)=sigr_t(ntl)
                ri(k)=ri_t(ntl)
                dsglc(k)=dsglc_t(ntl)
                ku=k+1
              enddo
c
c --- ------- save lower interface pressure
              pres(kmin)=delp(kmin)
              do k=kmin+1,kmax
                pres(k)=pres(k-1)+delp(k)
              enddo
c
              do ntl=1,ntil
c
                k=ktil(ntl)
                ku=k-1
                do while ((istil(ku).or.delp(ku).le.dpdmin)
     .                    .and.ku.gt.kmin)
                  ku=ku-1
                enddo
                ktilu(ntl)=ku
                kl=k+1
                do while ((istil(kl).or.delp(kl).le.dpdmin)
     .                    .and.kl.lt.kmax)
                  kl=kl+1
                enddo
                ktill(ntl)=kl
c
c --- --------- ennsq = buoy.freq.^2
                ennsq=.5*g*g*max(0.,dens(kl)-dens(ku))
     .                /max(delp(k),onecm)
c
                if (ri(k).lt.ri0) then
c
c --- ----------- maximum diffusivity is increased near the bottom to
c --- ----------- provide additional mixing of gravity currents
                  q=(pres(kmax)-pres(k)+.5*delp(k))
     .              /min(dpgc,.5*pres(kmax))
                  q=max(0.,1.-q*q)
                  q=q*q*q
                  nus=q*nug0+(1.-q)*nus0
c
c --- ----------- parameterization of diffusivity as a function of local
c --- ----------- gradient richardson number
                  q=ri(k)/ri0
                  q=max(0.,1.-q*q)
                  nus=nus*q*q*q
                else
                  nus=0.
                endif
c
c --- --------- enhance diffusivity when local stability is weak
                if (dsglc(k).lt.dsglc0) then
                  q=max(0.,dsglc(k))/dsglc0
                  q=max(0.,1.-q*q)
                  nus=nus+nuls0*q*q*q
                endif
c
c --- --------- store (exch.coeff x buoy.freq.^2 / g x time step) in
c --- --------- -flngth- (dimensions of flngth: length in pressure units)
c --- 
c --- --------- use the following if exch.coeff. = diapyc / buoyancy
c --- --------- frequency
                flngth=(diapyc*sqrt(ennsq)+nus*ennsq)*delt1/alpha0
c ---
c --- --------- use the following if exch.coeff. = diapyc
ccc             flngth=(diapyc+nus)*ennsq*delt1/alpha0
c
                alfa=-alpha0*dsigdt(ttem(k),ssal(k))
                beta= alpha0*dsigds(ttem(k),ssal(k))
c
                dsgiu(k)=1./max(beta*(ssal(k)-ssal(ku))
     .                         -alfa*(ttem(k)-ttem(ku)),small)
                dsgil(k)=1./max(beta*(ssal(kl)-ssal(k))
     .                         -alfa*(ttem(kl)-ttem(k)),small)
c
                flxmxu=flngth*dsgiu(k)
                flxmxl=flngth*dsgil(k)
                q=min(1.,maxflx/max(flxmxu,flxmxl),
     .                .25*delp(ku)/flxmxu,.25*delp(kl)/flxmxl)
                flxmxu=q*flxmxu
                flxmxl=q*flxmxl
                delp(ku)=delp(ku)-flxmxu
                delp(kl)=delp(kl)-flxmxl
                f=delp(k)
                delp(k)=delp(k)+flxmxu+flxmxl
                q=1./delp(k)
                ssal(k)=(f*ssal(k)+flxmxu*ssal(ku)+flxmxl*ssal(kl))*q
                ttem(k)=(f*ttem(k)+flxmxu*ttem(ku)+flxmxl*ttem(kl))*q
                dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                do nt=1,ntrc  
                  ttrc(nt,k)=(f*ttrc(nt,k)+flxmxu*ttrc(nt,ku)
     .                                    +flxmxl*ttrc(nt,kl))*q
                enddo
#endif
c
              enddo
c
              do ntl=1,ntil
c
                k=ktil(ntl)
                ku=ktilu(ntl)
                kl=ktill(ntl)
c
c --- --------- store mixing length required for density restoration in
c --- --------- -flngth- (dimensions of flngth: length in pressure units)
                flngth=(dens(k)-sigr(k))*delp(k)
c
                alfa=-alpha0*dsigdt(ttem(k),ssal(k))
                beta= alpha0*dsigds(ttem(k),ssal(k))
c
                dsgiu(k)=1./max(beta*(ssal(k)-ssal(ku))
     .                         -alfa*(ttem(k)-ttem(ku)),small)
                dsgil(k)=1./max(beta*(ssal(kl)-ssal(k))
     .                         -alfa*(ttem(kl)-ttem(k)),small)
c
                flxdru=0.
                flxdrl=0.
                if     (flngth.gt.0.) then
                  flxdru=flngth*dsgiu(k)
                  q=min(1.,.25*delp(ku)/flxdru)
                  if (q.lt.1.) then
                    flxdru=q*flxdru
                    flxdrl=max(-.25*delp(k),(q-1.)*flngth*dsgil(k))
                  endif
                elseif (flngth.lt.0.) then
                  flxdrl=-flngth*dsgil(k)
                  q=min(1.,.25*delp(kl)/flxdrl)
                  if (q.lt.1.) then
                    flxdrl=q*flxdrl
                    flxdru=max(-.25*delp(k),(1.-q)*flngth*dsgiu(k))
                  endif
                endif
c
                a=delp(ku)
                b=delp(k )
                c=delp(kl)
                delp(ku)=delp(ku)-flxdru
                delp(k )=delp(k )+flxdru+flxdrl
                delp(kl)=delp(kl)-flxdrl
                q=1./delp(k)
                ssal(k)=(b*ssal(k)+flxdru*ssal(ku)+flxdrl*ssal(kl))*q
                ttem(k)=(b*ttem(k)+flxdru*ttem(ku)+flxdrl*ttem(kl))*q
                dens(k)=sig(ttem(k),ssal(k))
#ifdef TRC
                if (flxdru.lt.0.) then
                  q=1./delp(ku)
                  do nt=1,ntrc  
                    ttrc(nt,ku)=(a*ttrc(nt,ku)-flxdru*ttrc(nt,k))*q
                  enddo
                endif
                if (flxdrl.lt.0.) then
                  q=1./delp(kl)
                  do nt=1,ntrc  
                    ttrc(nt,kl)=(c*ttrc(nt,kl)-flxdrl*ttrc(nt,k))*q
                  enddo
                endif
                q=1./delp(k)
                do nt=1,ntrc  
                  ttrc(nt,k)=
     .              ((b+min(0.,flxdru)+min(0.,flxdrl))*ttrc(nt,k)
     .              +max(0.,flxdru)*ttrc(nt,ku)
     .              +max(0.,flxdrl)*ttrc(nt,kl))*q
                enddo
#endif
c
              enddo
c
            endif
c
c           tsum1=0.
c           ssum1=0.
c           do k=kmin,kmax
c             tsum1=tsum1+ttem(k)*delp(k)
c             ssum1=ssum1+ssal(k)*delp(k)
c           enddo
c           if (abs(tsum1-tsum0)/p(i,j,kk+1).gt.1.e-12.or.
c    .          abs(ssum1-ssum0)/p(i,j,kk+1).gt.1.e-12) then
c             write (*,*) i,j,kmin,kmax
c             write (*,*) 'tsum',tsum0/p(i,j,kk+1),tsum1/p(i,j,kk+1)
c             write (*,*) 'ssum',ssum0/p(i,j,kk+1),ssum1/p(i,j,kk+1)
c           endif
c           if (ntrc.gt.0) then
c             tsum1=0.
c             do k=kmin,kmax
c               tsum1=tsum1+ttrc(1,k)*delp(k)
c             enddo
c             if (abs(tsum1-tsum0)/p(i,j,kk+1).gt.1.e-12) then
c               write (*,*) i,j,kmin,kmax
c               write (*,*) 'tsum',tsum0/p(i,j,kk+1),tsum1/p(i,j,kk+1)
c               stop
c             endif
c           endif
c
c           if (i.eq.itest.and.j.eq.jtest) then
c             write (*,*) 'kmin,kmax',kmin,kmax
c             do k=kmin,kmax
c               write (*,'(i2,3f12.6,3e12.6)')
c    .            k,delp(k)/onem,ttem(k),ssal(k),ri(k),dens(k),
c    .            dens(k)-sigr(k)
c             enddo
c             write (*,*) 'tsum1,ssum1',
c    .                    tsum1/p(i,j,kk+1),ssum1/p(i,j,kk+1)
c             write (*,*)
c             do nt=1,ntil
c               write (*,*) nt,ktil(nt),ktilu(nt),ktill(nt),
c    .                      ktil(max(1,nt-1))-ktil(nt)
c             enddo
c           endif
c
            if (restmx) then
c
c --- ------- restore mixed layer thickness
              if (delp(kmin).lt.delpml) then
                pres(kmin)=0.
                do k=kmin,kmax
                  pres(k+1)=pres(k)+delp(k)
                enddo
                tdp=ttem(kmin)*delp(kmin)
                sdp=ssal(kmin)*delp(kmin)
#ifdef TRC
                do nt=1,ntrc  
                  trdp(nt)=ttrc(nt,kmin)*delp(kmin)
                enddo
#endif
                k=kmin+1
                do while (pres(k).lt.delpml)
                  q=min(pres(k+1),delpml)-pres(k)
                  tdp=tdp+ttem(k)*q
                  sdp=sdp+ssal(k)*q
#ifdef TRC
                  do nt=1,ntrc  
                    trdp(nt)=trdp(nt)+ttrc(nt,k)*q
                  enddo
#endif
                  delp(k)=max(pres(k+1),delpml)-max(pres(k),delpml)
                  k=k+1
                enddo
                q=1./delpml
                delp(kmin)=delpml
                ttem(kmin)=tdp*q
                ssal(kmin)=sdp*q
#ifdef TRC
                do nt=1,ntrc  
                  ttrc(nt,kmin)=trdp(nt)*q
                enddo
#endif
              endif
c
            endif
c
          endif
c
c --- --- copy variables from 1d arrays back to main arrays and set t/s
c --- --- in remaining massless layers
          ttem(1)=ttem(kmin)
          ssal(1)=ssal(kmin)
          delp(1)=delp(kmin)
          dens(1)=dens(kmin)
#ifdef TRC
          do nt=1,ntrc
            ttrc(nt,1)=ttrc(nt,kmin)
          enddo
#endif
          if (kmin.ge.kmax) then
            do k=2,kk
              ttem(k)=max(ttem(1),temmin(k))
              dens(k)=sigr(k)
              ssal(k)=sofsig(dens(k),ttem(k))
              delp(k)=0.
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,1)
              enddo
#endif
            enddo
          else
            do k=2,kmin
              ttem(k)=ttem(kmin+1)
              dens(k)=sigr(k)
              ssal(k)=sofsig(dens(k),ttem(k))
              delp(k)=0.
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,kmin+1)
              enddo
#endif
            enddo
            do k=kmax+1,kk
              ttem(k)=ttem(kmax)
              dens(k)=sigr(k)
              ssal(k)=sofsig(dens(k),ttem(k))
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,kmax)
              enddo
#endif
            enddo
          endif
c
          do k=1,kk
            kn=k+nn
            temp(i,j,kn)=ttem(k)
            saln(i,j,kn)=ssal(k)
            diaflx(i,j,k)=diaflx(i,j,k)+delp(k)-dp(i,j,kn)
            dp(i,j,kn)=delp(k)
            sigma(i,j,kn)=dens(k)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
#ifdef TRC
            do nt=1,ntrc
              trc(i,j,k,nt)=ttrc(nt,k)
            enddo
#endif
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xcmaxr(rmax_niter)
      if (nint(rmax_niter).eq.10) then
        call xcstop('(diapfl: no convergence!)')
               stop '(diapfl: no convergence!)'
      endif
c
c --- ------------------------------------------------------------------
c --- Redistribute momentum
c --- ------------------------------------------------------------------
c
      call xctilr(p, 1,kk+1, 1,1, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(kn,uo,po,pn,ko,un,udpn)
      do j=1,jj
c
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k=1,kk
            kn=k+nn
            uo(k)=u(i,j,kn)
          enddo
          po(1)=0.
          pn(1)=0.
          do k=2,kk+1
            po(k)=pu(i,j,k)
            pn(k)=min(pbu(i,j,n),.5*(p(i  ,j,k)+p(i-1,j,k)),
     .                pu(i,j,kk+1))
          enddo
c
          ko=1
          do kn=1,kk
            if (pn(kn+1)-pn(kn).eq.0.) then
              un(kn)=0.
            else
              udpn=0.
              do while (pn(kn+1).gt.po(ko+1))
                udpn=udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko=ko+1
              enddo
              un(kn)=(udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))
     .               /(pn(kn+1)-pn(kn))
            endif
          enddo
          do k=1,kk
            kn=k+nn
            u(i,j,kn)=un(k)
          enddo
c
        enddo
        enddo
c
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          do k=1,kk
            kn=k+nn
            uo(k)=v(i,j,kn)
          enddo
          po(1)=0.
          pn(1)=0.
          do k=2,kk+1
            po(k)=pv(i,j,k)
            pn(k)=min(pbv(i,j,n),.5*(p(i,j  ,k)+p(i,j-1,k)),
     .                pv(i,j,kk+1))
          enddo
c
          ko=1
          do kn=1,kk
            if (pn(kn+1)-pn(kn).eq.0.) then
              un(kn)=0.
            else
              udpn=0.
              do while (pn(kn+1).gt.po(ko+1))
                udpn=udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko=ko+1
              enddo
              un(kn)=(udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))
     .               /(pn(kn+1)-pn(kn))
            endif
          enddo
          do k=1,kk
            kn=k+nn
            v(i,j,kn)=un(k)
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'diapfl:'
        endif
        call chksummsk(rig,ip,kk,'rig')
        call chksummsk(p,ip,kk+1,'p')
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(diaflx,ip,kk,'diaflx')
        call chksummsk(u,iu,2*kk,'u')
        call chksummsk(v,iv,2*kk,'v')
#ifdef TRC
        do nt=1,ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
