      subroutine diapfl(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Diapycnal mixing
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
c --- Parameters:
c ---   dsgmnr - minimum ratio of linearized density jump to target
c ---            density jump across a layer interface []
c ---   fcmxr  - maximum ratio of density restoration flux to the
c ---            estimated mean diapycnal flux []
c ---   dsgcr0 - gradually reduce the limiting of the density
c ---            restoration flux if the ratio of reference density
c ---            deviation to local vertical reference density
c ---            difference is greater than dsgcr0 []
c ---   dfeps  - small number used in the estimation of flux change
c ---            limit []
c ---   gbbl   - efficientcy factor for bottom boundary layer mixing []
c ---   kappa  - von Karman constant []
c ---   ustmin - minimum value of ustar used in computing the length
c ---            scale bottom boundary layer mixing [cm/s]
      real dsgmnr,fcmxr,dsgcr0,dfeps,gbbl,kappa,ustmin
      parameter (dsgmnr=.1,fcmxr=.25,dsgcr0=.25,dfeps=1.e-12,gbbl=.2,
     .           kappa=.4,ustmin=.01)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: fpug,fplg
      real, dimension(kdm) ::
     .  ttem,ssal,delp,dens,sigr,nu,fpu,fpl,fcu,fcl,dsgu,dsgl,dsghm,dsg,
     .  dsgui,dsgli,fmax,f,f0,fold,h,gtd,uc
      real, dimension(kdm+1) :: pres
      real c,delpu,delpl,nubbl,dsgdt,dsgds,fcmx,dsgc,q,dflim,ctd,bitd,r,
     .     s,t,dfdg,atd,maxdf,dtd
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: kming
      integer i,j,k,l,kn,kfpl,kmin,kmax,kl,kfmaxu
      logical, dimension(kdm) :: rstdns
      logical done,dwnwrd,remfmx
c
      real dsigdt,dsigds,sig,sofsig
      external dsigdt,dsigds,sig,sofsig
c
      real, dimension(kdm) :: ttem0,ssal0,delp0,dens0,sigr0,nu0
      integer niter
c
#if defined(TKE) && defined(TKEDDF)
#  include "param_tke.h"
#  include "common_tke.h"
      real, dimension(kdm) :: ttke
#  ifdef GLS
      real, dimension(kdm) :: tgls
#  endif
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc
      integer nt
#endif
c
c --- Constant in the diffusion equation
      c=g*g*delt1/(alpha0*alpha0)
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ ttem,ssal,delp,dens,sigr,nu,fpu,fpl,fcu,fcl,dsgu,dsgl,dsghm,dsg,
c$OMP+ dsgui,dsgli,fmax,f,f0,fold,h,gtd,pres,delpu,delpl,nubbl,dsgdt,
c$OMP+ dsgds,fcmx,dsgc,q,dflim,ctd,bitd,r,s,t,dfdg,atd,maxdf,dtd,k,kn,
c$OMP+ kfpl,kmin,kmax,kl,kfmaxu,rstdns,done,dwnwrd,remfmx,ttem0,ssal0,
c$OMP+ delp0,dens0,sigr0,nu0,niter
#if defined(TKE) && defined(TKEDDF)
c$OMP+,ttke
#  ifdef GLS
c$OMP+,tgls
#  endif
#endif
#ifdef TRC
c$OMP+,ttrc
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Copy variables into 1d arrays.
          do k=1,kk
            kn=k+nn
            ttem(k)=temp(i,j,kn)
            ssal(k)=saln(i,j,kn)
            delp(k)=dp(i,j,kn)
            dens(k)=sigma(i,j,kn)
            sigr(k)=sigmar(i,j,k)
            nu(k)=difdia(i,j,k)
            rstdns(k)=.true.
cdiag
            ttem0(k)=ttem(k)
            ssal0(k)=ssal(k)
            delp0(k)=delp(k)
            dens0(k)=dens(k)
            sigr0(k)=sigr(k)
            nu0(k)=nu(k)
cdiag
#if defined(TKE) && defined(TKEDDF)
            ttke(k)=tke(i,j,kn)
#  ifdef GLS
            tgls(k)=gls_psi(i,j,kn)
#  endif
#endif
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,k)=trc(i,j,k,nt)
            enddo
#endif
          enddo
c
c --- --- Locate range of physical layers.
          kfpl=kfpla(i,j,n)
          kmin=kfpl-2
          kmax=1
          do k=2,kk
            if (delp(k).gt.epsil) kmax=k
          enddo
c
          if (kmin.lt.kmax) then
c
c --- ----- Locate first layer where potential density should be kept
c --- ----- close to the reference potential density.
            rstdns(kfpl)=.false.
            if (kfpl.ne.kmax) then
              if (dens(kfpl).gt..5*(sigr(kfpl)+sigr(kfpl+1))) then
                rstdns(kfpl+1)=.false.
              endif
            endif
c
c --- ----- Copy mixed layer variables to the layers with indexes kmin
c --- ----- and kmin+1
            delp(kmin+1)=delp(2)
            delp(kmin  )=delp(1)
            ttem(kmin+1)=ttem(2)
            ttem(kmin  )=ttem(1)
            ssal(kmin+1)=ssal(2)
            ssal(kmin  )=ssal(1)
            nu(kmin+1)=nu(2)
            nu(kmin  )=nu(1)
#if defined(TKE) && defined(TKEDDF)
            ttke(kmin+1)=ttke(2)
            ttke(kmin  )=ttke(1)
#  ifdef GLS
            tgls(kmin+1)=tgls(2)
            tgls(kmin  )=tgls(1)
#  endif
#endif
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,kmin+1)=ttrc(nt,2)
              ttrc(nt,kmin  )=ttrc(nt,1)
            enddo
#endif
c
c --- ----- Find interface pressure
            pres(kmin)=0.
            do k=kmin,kmax
              pres(k+1)=pres(k)+delp(k)
            enddo
c
c --- ----- Compute mass flux between the layers in the mixed layers and
c --- ----- at the mixed layer base.  
            k=kmin
            fpu(k)=0.
            fpl(k)=min(pres(k+1),pres(kmax+1)-pres(k+1),
     .                 c*nu(k)*(delp(k)+delp(k+1))
     .                 /(2.*delp(k)*delp(k+1)))
            k=kmin+1
            fpu(k)=fpl(k-1)
            delpu=max(onem,delp(k))
            delpl=max(onem,delp(k+1))
            fpl(k)=min(pres(k+1),pres(kmax+1)-pres(k+1),
     .                 c*nu(k)*(delpu+delpl)/(2.*delpu*delpl))
            fpl(kmax)=0.
c
            if (kfpl.le.kmax) then
c
              if (kfpl.lt.kmax) then
c
c --- --------- Bottom boundary mixing is parameterized by assuming that
c --- --------- a part of the energy extracted from the mean flow by the
c --- --------- bottom drag drives diapycnal mixing.
                k=kmax-1
                nubbl=gbbl*ustarb(i,j)**3
     .                *exp(-(delp(k+1)+.5*delp(k))*abs(coriop(i,j))
     .                      *alpha0/(kappa*max(ustmin,ustarb(i,j))*g))
     .                /(alpha0*g*(sigr(k+1)-sigr(k)))
                nu(k)=max(nu(k),nubbl)
                difdia(i,j,k)=nu(k)
              endif
c
c --- ------- Compute linearized density jumps across upper and lower
c --- ------- interfaces of layers, average density jumps for layers,
c --- ------- and buoyancy flux corrections to nudge the layer density
c --- ------- towards the reference density. The flux corrections are
c --- ------- limited so their absolute value is guaranteed to be less
c --- ------- than the buoyancy flux due to diffusion.
              k=kfpl-1
              dsgli(k)=1.
              fcl(k)=-fpl(k)
              do k=kfpl,kmax-1
                if (rstdns(k)) then
                  dsgdt=dsigdt(ttem(k),ssal(k))
                  dsgds=dsigds(ttem(k),ssal(k))
                  dsgu(k)=max(dsgmnr*(sigr(k)-sigr(k-1)),
     .                        dsgdt*(ttem(k)-ttem(k-1))
     .                       +dsgds*(ssal(k)-ssal(k-1)))
                  dsgl(k)=max(dsgmnr*(sigr(k+1)-sigr(k)),
     .                        dsgdt*(ttem(k+1)-ttem(k))
     .                       +dsgds*(ssal(k+1)-ssal(k)))
                  dsghm(k)=2.*dsgu(k)*dsgl(k)/(dsgu(k)+dsgl(k))
                  dsg(k)=.5*(dsgu(k)+dsgl(k))
c                 dsg(k)=dsghm(k)
                  dsgui(k)=1./dsgu(k)
                  dsgli(k)=1./dsgl(k)
                  fcmx=.25*(sqrt(delp(k)*delp(k)
     .                          +4.*c*nu(k)*dsg(k)*(dsgui(k)+dsgli(k)))
     .                     -delp(k))*dsghm(k)*fcmxr
                  dsgc=dens(k)-sigr(k)
                  if (dsgc.gt.0.) then
                    fcl(k)=0.
                    if (dens(k-1).lt.sigr(k)) then
                      q=max(0.,(dens(k)-sigr(k+1))
     .                         /((sigr(k)-sigr(k+1))*(1.-dsgcr0)))
                      q=max(0.,1.-q*q)
                      q=q*q*q
                      fcu(k)=dsgc*delp(k)
                      fcu(k)=min(q*fcu(k)+(1.-q)*fcmx,fcu(k))
                    else
                      fcu(k)=0.
                    endif
                  else
                    fcu(k)=0.
                    if (dens(k+1).gt.sigr(k)) then
                      q=max(0.,(dens(k)-sigr(k-1))
     .                         /((sigr(k)-sigr(k-1))*(1.-dsgcr0)))
                      q=max(0.,1.-q*q)
                      q=q*q*q
                      fcl(k)=dsgc*delp(k)
                      fcl(k)=max(q*fcl(k)-(1.-q)*fcmx,fcl(k))
                    else
                      fcl(k)=0.
                    endif
                  endif
                else
                  dsgu(k)=1.
                  dsgl(k)=1.
                  dsghm(k)=1.
                  dsg(k)=1.
                  dsgui(k)=1.
                  dsgli(k)=1.
                  fcl(k)=0.
                  fcu(k)=0.
                endif
              enddo
              k=kmax
              dsgdt=dsigdt(ttem(k),ssal(k))
              dsgds=dsigds(ttem(k),ssal(k))
              dsgu(k)=max(dsgmnr*(sigr(k)-sigr(k-1)),
     .                    dsgdt*(ttem(k)-ttem(k-1))
     .                   +dsgds*(ssal(k)-ssal(k-1)))
              dsgui(k)=1./dsgu(k)
              if (dens(k).gt.sigr(k).and.dens(k-1).lt.sigr(k)) then
                fpu(k)=min(delp(k-1),
     .                     (dens(k)-sigr(k))*delp(k)*dsgui(k))
              else
                fpu(k)=0.
              endif
              fcu(k)=fpu(k)*dsgu(k)
c
c --- ------- Find maximum allowable buoyancy fluxes and further limit
c --- ------- the flux corrections to keep layer interfaces within the
c --- ------- fluid domain.
              fmax(kfpl-1)=0.
              fmax(kmax)=0.
              done=.false.
              niter=0
              do while (.not.done)
                done=.true.
                do k=kmax-1,kfpl,-1
                  q=((fmax(k+1)+fcu(k+1))*dsgui(k+1)
     .              +pres(kmax+1)-pres(k+1))*dsgl(k)
                  fcl(k)=max(-q,fcl(k))
                  fmax(k)=q+fcl(k)
                enddo
                kfmaxu=0
                do k=kfpl,kmax-1
                  q=((fmax(k-1)-fcl(k-1))*dsgli(k-1)
     .              +pres(k)-pres(kfpl))*dsgu(k)
                  if (fcu(k).gt.q) then
                    fcu(k)=q
                    done=.false.
                  endif
                  if (fmax(k).gt.q-fcu(k)) then
                    fmax(k)=q-fcu(k)
                    kfmaxu=k
                  endif
                enddo
                if (niter.eq.100) then
                  write (lp,*)
     .              'micom: diapfl: no convergence in flux limit!',
     .              i+i0,j+j0
                  open (10,file='diapfl.uf',form='unformatted')
                  write (10) kk,kfpl
                  write (10) g,alpha0,epsil,onem,delt1,dsgmnr,q,q
                  write (10) ttem0,ssal0,delp0,dens0,sigr0,nu0
                  close (10)
                  exit
                endif
              enddo
c
c --- ------- Make a first guess for buoyancy fluxes and set some
c --- ------- boundary conditions
              k=kfpl-1
              f0(k)=0.
              f(k)=0.
              gtd(k)=0.
              dflim=0.
              do k=kfpl,kmax-1
                f(k)=min(fmax(k),
     .                   .5*sqrt(c*nu(k)*dsg(k)
     .                           *(dsgui(k)+dsgli(k)))*dsghm(k),
     .                   c*nu(k)*dsg(k)/max(epsil,delp(k)))
                fold(k)=f(k)
                h(k)=fcu(k  )*dsgui(k  )-fcl(k  )*dsgli(k  )
     .              +fcl(k-1)*dsgli(k-1)-fcu(k+1)*dsgui(k+1)
                dflim=max(dflim,fmax(k))
              enddo
              k=kmax
              f0(k)=0.
              f(k)=0.
              gtd(k)=0.
              dflim=dflim*dfeps
c
c --- ------- Solve the diffusion equation for layer thickness using
c --- ------- backward time integration by an iterative algorithm.
              niter=0.
              dwnwrd=.false.
              do
c
c --- --------- Solve the equation by alternate downward and upward
c --- --------- propagation trough trough the layers
                dwnwrd=.not.dwnwrd
                if (dwnwrd) then
c
c --- ----------- Do a downward first pass.
                  ctd=0.
                  bitd=1.
                  remfmx=.false.
                  do k=kfpl,kmax-1
c
                    if (remfmx) then
                      gtd(k)=0.
                      f0(k)=fmax(k)
                      f(k)=fmax(k)
                    else
c
c --- --------------- Find the backward solution of the layer buoyancy
c --- --------------- flux, assuming no dependency on adjacent layer
c --- --------------- fluxes, and the sensitivity of the layer flux to
c --- --------------- changes in the fluxes of adjacent layers.
                      q=f0(k-1)*dsgli(k-1)+f(k+1)*dsgui(k+1)
     .                 -delp(k)-h(k)
                      r=4.*c*nu(k)*dsg(k)*(dsgui(k)+dsgli(k))
                      t=.25*dsghm(k)
                      if (q.lt.0.) then
                        s=r/(q*q)
                        if (s.lt.1.e-3) then
c
c --- ------------------- For certain parameters, use a taylor expansion
c --- ------------------- of the flux and flux sensitivity expressions
c --- ------------------- to avoid roundoff errors.
                          r=.00390625*s
                          q=-q*r*(128.-s*(32.-s*(16.-s*(10.-s*7. ))))
                          f0(k)=q*t
                          q=r*(128.-s*(96.-s*(80.-s*(70.-s*63.))))
                          dfdg=q*t
                        else
                          s=sqrt(q*q+r)
                          f0(k)=(q+s)*t
                          dfdg=(1.+q/s)*t
                        endif
                      else
                        s=sqrt(q*q+r)
                        f0(k)=(q+s)*t
                        dfdg=(1.+q/s)*t
                      endif
c
                      if (f0(k).ge.fmax(k)) then
c
c --- ----------------- If the maximum flux associated with the lower
c --- ----------------- fluid boundary has been reached, all subsequent
c --- ----------------- fluxes will be set to the maximum allowable
c --- ----------------- flux.
                        f0(k)=fmax(k)
                        dfdg=0.
                        if (k.gt.kfmaxu) remfmx=.true.
                      endif
c
c --- --------------- Modify the buoyancy fluxes by taking into account
c --- --------------- linearized contributions of adjacent layer fluxes.
c --- --------------- This linearization forms a tridiagonal set of
c --- --------------- equations.
                      gtd(k)=ctd*bitd
                      atd=-dfdg*dsgli(k-1)
                      ctd=-dfdg*dsgui(k+1)
                      bitd=1./(1.-atd*gtd(k))
                      f(k)=(f0(k)-atd*(f(k-1)-f0(k-1))+ctd*f(k+1))*bitd
                    endif
                  enddo
c
c --- ----------- Complete the solving of the tridiagonal set of
c --- ----------- equations and find the maximum flux change.
                  maxdf=0.
                  do k=kmax-1,kfpl,-1
                    f(k)=min(fmax(k),f(k)-gtd(k+1)*f(k+1))
                    maxdf=max(maxdf,abs(f(k)-fold(k)))
                    fold(k)=f(k)
                  enddo
                else
c
c --- ----------- Do a upward first pass.
                  atd=0.
                  bitd=1.
                  remfmx=.false.
                  do k=kmax-1,kfpl,-1
c
                    if (remfmx) then
                      gtd(k)=0.
                      f0(k)=fmax(k)
                      f(k)=fmax(k)
                    else
c
c --- --------------- Find the backward solution of the layer buoyancy
c --- --------------- flux, assuming no dependency on adjacent layer
c --- --------------- fluxes, and the sensitivity of the layer flux to
c --- --------------- changes in the fluxes of adjacent layers.
                      q=f(k-1)*dsgli(k-1)+f0(k+1)*dsgui(k+1)
     .                 -delp(k)-h(k)
                      r=4.*c*nu(k)*dsg(k)*(dsgui(k)+dsgli(k))
                      t=.25*dsghm(k)
                      if (q.lt.0.) then
                        s=r/(q*q)
                        if (s.lt.1.e-3) then
c
c --- ------------------- For certain parameters, use a taylor expansion
c --- ------------------- of the flux and flux sensitivity expressions
c --- ------------------- to avoid roundoff errors.
                          r=.00390625*s
                          q=-q*r*(128.-s*(32.-s*(16.-s*(10.-s*7. ))))
                          f0(k)=q*t
                          q=r*(128.-s*(96.-s*(80.-s*(70.-s*63.))))
                          dfdg=q*t
                        else
                          s=sqrt(q*q+r)
                          f0(k)=(q+s)*t
                          dfdg=(1.+q/s)*t
                        endif
                      else
                        s=sqrt(q*q+r)
                        f0(k)=(q+s)*t
                        dfdg=(1.+q/s)*t
                      endif
c
                      if (f0(k).ge.fmax(k)) then
c
c --- ----------------- If the maximum flux associated with the upper
c --- ----------------- fluid boundary has been reached, all subsequent
c --- ----------------- fluxes will be set to the maximum allowable
c --- ----------------- flux.
                        f0(k)=fmax(k)
                        dfdg=0.
                        if (k.le.kfmaxu) remfmx=.true.
                      endif
c
c --- --------------- Modify the buoyancy fluxes by taking into account
c --- --------------- linearized contributions of adjacent layer fluxes.
c --- --------------- This linearization forms a tridiagonal set of
c --- --------------- equations.
                      gtd(k)=atd*bitd
                      atd=-dfdg*dsgli(k-1)
                      ctd=-dfdg*dsgui(k+1)
                      bitd=1./(1.-ctd*gtd(k))
                      f(k)=(f0(k)+atd*f(k-1)-ctd*(f(k+1)-f0(k+1)))*bitd
                    endif
                  enddo
c
c --- ----------- Complete the solving of the tridiagonal set of
c --- ----------- equations and find the maximum flux change.
                  maxdf=0.
                  do k=kfpl,kmax-1
                    f(k)=min(fmax(k),f(k)-gtd(k-1)*f(k-1))
                    maxdf=max(maxdf,abs(f(k)-fold(k)))
                    fold(k)=f(k)
                  enddo
                endif
c
c --- --------- If the maximum flux change is below a treshold, stop the
c --- --------- iteration.
                niter=niter+1
                if (maxdf.le.dflim) exit
                if (niter.eq.100) then
                  write (lp,*)
     .          'micom: diapfl: no convergence in implicit diffusion!',
     .          i+i0,j+j0,maxdf,dflim
                  open (10,file='diapfl.uf',form='unformatted')
                  write (10) kk,kfpl
                  write (10) g,alpha0,epsil,onem,delt1,dsgmnr,q,q
                  write (10) ttem0,ssal0,delp0,dens0,sigr0,nu0
                  close (10)
                  exit
                endif
              enddo
c
c --- ------- Compute the mass fluxes
              do k=kfpl,kmax-1
                fpu(k)=(f(k)+fcu(k))*dsgui(k)
                fpl(k)=(f(k)-fcl(k))*dsgli(k)
              enddo
              fpu(kfpl)=fpl(kmin+1)
c
            endif
c
c --- ----- Solve the diffusion equation for temperature and salinity by
c --- ----- backward integration forming a tridiagonal set of equations.
            ctd=0.
            bitd=1.
            do k=kmin,kmax
              gtd(k)=ctd*bitd
              q=1./(delp(k)+fpu(k)+fpl(k))
              atd=-fpu(k)*q
              ctd=-fpl(k)*q
              dtd=delp(k)*q
              bitd=1./(1.-atd*gtd(k))
              ssal(k)=(dtd*ssal(k)-atd*ssal(max(1,k-1)))*bitd
              ttem(k)=(dtd*ttem(k)-atd*ttem(max(1,k-1)))*bitd
#if defined(TKE) && defined(TKEDDF)
              ttke(k)=(dtd*ttke(k)-atd*ttke(max(1,k-1)))*bitd
#  ifdef GLS
              tgls(k)=(dtd*tgls(k)-atd*tgls(max(1,k-1)))*bitd
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=(dtd*ttrc(nt,k)-atd*ttrc(nt,max(1,k-1)))*bitd
              enddo
#endif
            enddo
            do k=kmax-1,kmin,-1
              ssal(k)=ssal(k)-gtd(k+1)*ssal(k+1)
              ttem(k)=ttem(k)-gtd(k+1)*ttem(k+1)
              dens(k)=sig(ttem(k),ssal(k))
#if defined(TKE) && defined(TKEDDF)
              ttke(k)=ttke(k)-gtd(k+1)*ttke(k+1)
#  ifdef GLS
              tgls(k)=tgls(k)-gtd(k+1)*tgls(k+1)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,k)-gtd(k+1)*ttrc(nt,k+1)
              enddo
#endif
            enddo
            do k=kfpl,kmax-1
              delp(k)=max(0.,delp(k)+fpu(k)+fpl(k)-fpl(k-1)-fpu(k+1))
            enddo
            delp(kmax)=max(0.,delp(kmax)+fpu(kmax)-fpl(kmax-1))
c
c --- ----- Copy variables back to the mixed layers from the layers with
c --- ----- index kmin and kmin+1
            ttem(1)=ttem(kmin  )
            ttem(2)=ttem(kmin+1)
            ssal(1)=ssal(kmin  )
            ssal(2)=ssal(kmin+1)
            dens(1)=dens(kmin  )
            dens(2)=dens(kmin+1)
            if (kmin.gt.1) then
              if (kmin.eq.2) then
                delp(2)=delp(kmin+1)
                delp(kmin+1)=0.
              else
                delp(kmin  )=0.
              endif
            endif
#if defined(TKE) && defined(TKEDDF)
            ttke(1)=ttke(kmin  )
            ttke(2)=ttke(kmin+1)
#  ifdef GLS
            tgls(1)=tgls(kmin  )
            tgls(2)=tgls(kmin+1)
#  endif
#endif
#ifdef TRC
            do nt=1,ntrc
              ttrc(nt,1)=ttrc(nt,kmin  )
              ttrc(nt,2)=ttrc(nt,kmin+1)
            enddo
#endif
c
          endif
c
c --- --- Fill massless layers with resonable values of temperature,
c --- --- salinity, and tracers
          if (kfpl.gt.kmax) then
            do k=3,kk
              ttem(k)=max(ttem(2),temmin(i,j,k))
              dens(k)=sigr(k)
              ssal(k)=sofsig(dens(k),ttem(k))
              delp(k)=0.
#if defined(TKE) && defined(TKEDDF)
              ttke(k)=max(ttke(2),tke_min)
#  ifdef GLS
              tgls(k)=max(tgls(2),gls_psi_min)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,2)
              enddo
#endif
            enddo
          else
            do k=3,kfpl-1
              ttem(k)=ttem(kfpl)
              dens(k)=sigr(k)
              ssal(k)=sofsig(dens(k),ttem(k))
              delp(k)=0.
#if defined(TKE) && defined(TKEDDF)
              ttke(k)=ttke(kfpl)
#  ifdef GLS
              tgls(k)=tgls(kfpl)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,kfpl)
              enddo
#endif
            enddo
            do k=kmax+1,kk
              ttem(k)=ttem(kmax)
              dens(k)=sigr(k)
              ssal(k)=sofsig(dens(k),ttem(k))
#if defined(TKE) && defined(TKEDDF)
              ttke(k)=ttke(kmax)
#  ifdef GLS
              tgls(k)=tgls(kmax)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                ttrc(nt,k)=ttrc(nt,kmax)
              enddo
#endif
            enddo
          endif
c
c --- --- Copy 1d arrays to 3d arrays
          do k=1,kk
            kn=k+nn
            temp(i,j,kn)=ttem(k)
            saln(i,j,kn)=ssal(k)
            diaflx(i,j,k)=diaflx(i,j,k)+delp(k)-dp(i,j,kn)
            dp(i,j,kn)=delp(k)
            sigma(i,j,kn)=dens(k)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
#if defined(TKE) && defined(TKEDDF)
            tke(i,j,kn)=max(ttke(k),tke_min)
#  ifdef GLS
            gls_psi(i,j,kn)=max(tgls(k),gls_psi_min)
#  endif
#endif
#ifdef TRC
            do nt=1,ntrc
              trc(i,j,k,nt)=ttrc(nt,k)
            enddo
#endif
          enddo
c
c --- --- Save variables used for momentum mixing
          kming(i,j)=kmin
          if (kmin.lt.kmax) then
            do k=1,kmin
              fpug(i,j,k)=fpl(kmin)
              fplg(i,j,k)=fpl(kmin)
            enddo
            do k=kmin+1,kmax
              fpug(i,j,k)=fpu(k)
              fplg(i,j,k)=fpl(k)
            enddo
            do k=kmax+1,kk
              fpug(i,j,k)=0.
              fplg(i,j,k)=0.
            enddo
          else
            do k=1,kk
              fpug(i,j,k)=0.
              fplg(i,j,k)=0.
            enddo
          endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Diapycnal mixing of momentum.
c --- ------------------------------------------------------------------
c
      call xctilr(p, 1,kk+1, 1,1, halo_ps)
      call xctilr(fpug, 1,kk, 1,1, halo_ps)
      call xctilr(fplg, 1,kk, 1,1, halo_ps)
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kming(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kming(i,j)=nint(util1(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ kmin,kmax,uc,delp,kn,k,fpu,fpl,q,ctd,bitd,gtd,atd,ctd,dtd)
      do j=1,jj
c
c --- - Mixing of u-component.
c
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c
c --- --- Find range of mass containing layers.
          kmin=min(kming(i-1,j),kming(i,j))
          kmax=1
          do k=2,kk
            if (dpu(i,j,k+nn).gt.0.) kmax=k
          enddo
c
c --- --- If number of mass containing layers is less than 2, do not do
c --- --- anything in this water column.
          if (kmin.lt.kmax) then
c
c --- ----- Copy vertical column variables into 1d arrays.
            uc(kmin+1)=u(i,j,k1n+1)
            uc(kmin  )=u(i,j,k1n  )
            delp(kmin+1)=dpu(i,j,k1n+1)
            delp(kmin  )=dpu(i,j,k1n  )
            do k=kmin+2,kmax
              kn=k+nn
              uc(k)=u(i,j,kn)
              delp(k)=dpu(i,j,kn)
            enddo
c
c --- ----- Interpolate infterface mass fluxes to u-points. Limit mass
c --- ----- fluxes so the interfaces do not descend below the u-point
c --- ----- depth.
            k=kmin
            fpu(k)=0.
            fpl(k)=.5*(fplg(i-1,j,k)+fplg(i,j,k))
            q=0.
            do k=kmin+1,kmax
              fpu(k)=.5*(fpug(i-1,j,k)+fpug(i,j,k))
              fpl(k)=.5*(fplg(i-1,j,k)+fplg(i,j,k))
              q=q+delp(k-1)
              fpl(k-1)=min(pu(i,j,kk+1)-q+fpu(k),fpl(k-1))
            enddo
            fpl(kmax)=0.
c
c --- ----- Solve the diffusion equation for velocity by backward
c --- ----- integration forming a tridiagonal set of equations.
            ctd=0.
            bitd=1.
            do k=kmin,kmax
              gtd(k)=ctd*bitd
              q=1./(delp(k)+fpu(k)+fpl(k))
              atd=-fpu(k)*q
              ctd=-fpl(k)*q
              dtd=delp(k)*q
              bitd=1./(1.-atd*gtd(k))
              uc(k)=(dtd*uc(k)-atd*uc(max(kmin,k-1)))*bitd
            enddo
            do k=kmax-1,kmin,-1
              uc(k)=uc(k)-gtd(k+1)*uc(k+1)
            enddo
c
c --- ----- Put velocity back in main array.
            u(i,j,k1n  )=uc(kmin  )
            u(i,j,k1n+1)=uc(kmin+1)
            do k=kmin+2,kmax
              u(i,j,k+nn)=uc(k)
            enddo
c
c --- ----- If interfaces are lifted above the bottom because of
c --- ----- diapycnal mixing, give the newly opened layers the velocity
c --- ----- of the lowest initial mass containing layer.
            do k=kmax+1,kk
              if (.5*(p(i,j,k)+p(i-1,j,k)).lt.pu(i,j,kk+1))
     .          u(i,j,k+nn)=uc(kmax)
            enddo
c
          endif
c
        enddo
        enddo
c
c --- - Mixing of v-component.
c
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c
c --- --- Find range of mass containing layers.
          kmin=min(kming(i,j-1),kming(i,j))
          kmax=1
          do k=2,kk
            if (dpv(i,j,k+nn).gt.0.) kmax=k
          enddo
c
c --- --- If number of mass containing layers is less than 2, do not do
c --- --- anything in this water column.
          if (kmin.lt.kmax) then
c
c --- ----- Copy vertical column variables into 1d arrays.
            uc(kmin+1)=v(i,j,k1n+1)
            uc(kmin  )=v(i,j,k1n  )
            delp(kmin+1)=dpv(i,j,k1n+1)
            delp(kmin  )=dpv(i,j,k1n  )
            do k=kmin+2,kmax
              kn=k+nn
              uc(k)=v(i,j,kn)
              delp(k)=dpv(i,j,kn)
            enddo
c
c --- ----- Interpolate infterface mass fluxes to v-points. Limit mass
c --- ----- fluxes so the interfaces do not descend below the v-point
c --- ----- depth.
            k=kmin
            fpu(k)=0.
            fpl(k)=.5*(fplg(i,j-1,k)+fplg(i,j,k))
            q=0.
            do k=kmin+1,kmax
              fpu(k)=.5*(fpug(i,j-1,k)+fpug(i,j,k))
              fpl(k)=.5*(fplg(i,j-1,k)+fplg(i,j,k))
              q=q+delp(k-1)
              fpl(k-1)=min(pv(i,j,kk+1)-q+fpu(k),fpl(k-1))
            enddo
            fpl(kmax)=0.
c
c --- ----- Solve the diffusion equation for velocity by backward
c --- ----- integration forming a tridiagonal set of equations.
            ctd=0.
            bitd=1.
            do k=kmin,kmax
              gtd(k)=ctd*bitd
              q=1./(delp(k)+fpu(k)+fpl(k))
              atd=-fpu(k)*q
              ctd=-fpl(k)*q
              dtd=delp(k)*q
              bitd=1./(1.-atd*gtd(k))
              uc(k)=(dtd*uc(k)-atd*uc(max(kmin,k-1)))*bitd
            enddo
            do k=kmax-1,kmin,-1
              uc(k)=uc(k)-gtd(k+1)*uc(k+1)
            enddo
c
c --- ----- Put velocity back in main array.
            v(i,j,k1n  )=uc(kmin  )
            v(i,j,k1n+1)=uc(kmin+1)
            do k=kmin+2,kmax
              v(i,j,k+nn)=uc(k)
            enddo
c
c --- ----- If interfaces are lifted above the bottom because of
c --- ----- diapycnal mixing, give the newly opened layers the velocity
c --- ----- of the lowest initial mass containing layer.
            do k=kmax+1,kk
              if (.5*(p(i,j,k)+p(i,j-1,k)).lt.pv(i,j,kk+1))
     .          v(i,j,k+nn)=uc(kmax)
            enddo
c
          endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kn,q)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i-1,j,k+1))) 
     .              -min(q,.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)=
     .        max(0.,min(q,.5*(p(i,j,k+1)+p(i,j-1,k+1)))
     .              -min(q,.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c     do j=1,jj
c       do l=1,isu(j)
c       do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c         q=0.
c         do k=1,kk
c           kn=k+nn
c           q=q+u(i,j,kn)*dpu(i,j,kn)
c         enddo
c         if (abs(q).gt.1.e-4) then
c           write (lp,*) 'diapfl: u imbalance:',q,i,j
c         endif
c       enddo
c       enddo
c       do l=1,isv(j)
c       do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c         q=0.
c         do k=1,kk
c           kn=k+nn
c           q=q+v(i,j,kn)*dpv(i,j,kn)
c         enddo
c         if (abs(q).gt.1.e-4) then
c           write (lp,*) 'diapfl: v imbalance:',q,i,j
c         endif
c       enddo
c       enddo
c     enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'diapfl:'
        endif
        call chksummsk(p,ip,kk+1,'p')
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(diaflx,ip,kk,'diaflx')
        call chksummsk(u,iu,2*kk,'u')
        call chksummsk(v,iv,2*kk,'v')
#if defined(TKE) && defined(TKEDDF)
        call chksummsk(tke,ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(gls_psi,ip,2*kk,'gls_psi')
#  endif
#endif
#ifdef TRC
        do nt=1,ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
