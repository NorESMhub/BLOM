! ------------------------------------------------------------------------------
! Copyright (C) 2009-2022 Mats Bentsen, Mehmet Ilicak

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_diapfl

  use dimensions,    only: idm, jdm, kdm
  use mod_constants, only: g, alpha0, spval, epsilp, onem, L_mks2cgs
  use mod_time,      only: delt1
  use mod_xc,        only: xchalt, xctilr, ii, jj, kk, isp, ifp, ilp, &
                           i0, j0, lp, isu, ifu, ilu, isv, ifv, ilv, &
                           iu, iv, nbdy, halo_ps, mnproc, ip
  use mod_vcoord,    only: sigmar
  use mod_grid,      only: coriop
  use mod_eos,       only: sig, dsigdt, dsigds, sofsig
  use mod_temmin,    only: temmin
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                           p, pu, pv, kfpla
  use mod_diffusion, only: difdia
  use mod_forcing,   only: ustarb
  use mod_utility,   only: util1
  use mod_checksum,  only: csdiag, chksummsk
  use mod_tracers,   only: ntr, itrtke, itrgls, trc
  use mod_tke,       only: tke_min, gls_psi_min
  use mod_ifdefs,    only: use_TRC, use_TKE, use_GLS

  implicit none
  private

  public :: diapfl

contains

  subroutine diapfl(n,nn,k1n)

    ! --- ------------------------------------------------------------------
    ! --- Diapycnal mixing
    ! --- ------------------------------------------------------------------


    ! Arguments
    integer, intent(in) :: n,nn,k1n

    ! Local variables
    real :: dsgmnr,fcmxr,dsgcr0,dfeps,gbbl,kappa,ustmin
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), save :: fpug=spval,fplg = spval
    real, dimension(kdm) :: ttem,ssal,delp,dens,sigr,nu,fpu,fpl,fcu,fcl,dsgu,dsgl,dsghm,dsg
    real, dimension(kdm) :: dsgui,dsgli,fmax,f,f0,fold,h,gtd,uc
    real, dimension(kdm+1) :: pres
    real :: c,delpu,delpl,nubbl,dsgdt,dsgds,fcmx,dsgc,q,dflim,ctd,bitd,r
    real :: s,t,dfdg,atd,maxdf,dtd,pold,pnew,fpum,fpup,fplm,fplp
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: kming
    integer :: i,j,k,l,kn,kfpl,kmin,kmax,kfmaxu
    logical, dimension(kdm) :: rstdns
    logical :: done,dwnwrd,remfmx
    real, dimension(kdm) :: ttem0,ssal0,delp0,dens0,sigr0,nu0
    integer :: niter
    real, dimension(ntr,kdm) :: ttrc
    integer :: nt

    ! --- Parameters:
    ! ---   dsgmnr - minimum ratio of linearized density jump to target
    ! ---            density jump across a layer interface []
    ! ---   fcmxr  - maximum ratio of density restoration flux to the
    ! ---            estimated mean diapycnal flux []
    ! ---   dsgcr0 - gradually reduce the limiting of the density
    ! ---            restoration flux if the ratio of reference density
    ! ---            deviation to local vertical reference density
    ! ---            difference is greater than dsgcr0 []
    ! ---   dfeps  - small number used in the estimation of flux change
    ! ---            limit []
    ! ---   gbbl   - efficientcy factor for bottom boundary layer mixing []
    ! ---   kappa  - von Karman constant []
    ! ---   ustmin - minimum value of ustar used in computing the length
    ! ---            scale bottom boundary layer mixing [cm/s]
    parameter(dsgmnr=.1,fcmxr=.25,dsgcr0=.25,dfeps=1.e-12,gbbl=.2,kappa=.4,ustmin=.0001*L_mks2cgs)

    ! --- Constant in the diffusion equation
    c = g*g*delt1/(alpha0*alpha0)

    !$omp parallel do private( &
    !$omp l,i,k,kn,ttem,ssal,delp,dens,sigr,nu,rstdns,ttem0,ssal0,delp0, &
    !$omp dens0,sigr0,nu0,kfpl,kmin,kmax,pres,fpu,fpl,delpu,delpl,nubbl, &
    !$omp dsgli,fcl,dsgdt,dsgds,dsgu,dsgl,dsghm,dsg,dsgui,fcmx,dsgc,q,fcu, &
    !$omp fmax,done,niter,kfmaxu,f0,f,gtd,dflim,fold,h,dwnwrd,ctd,bitd, &
    !$omp remfmx,r,t,s,dfdg,atd,maxdf,dtd, &
    !$omp nt,ttrc)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! --- --- Copy variables into 1d arrays.
          do k = 1,kk
            kn = k+nn
            ttem(k) = temp(i,j,kn)
            ssal(k) = saln(i,j,kn)
            delp(k) = dp(i,j,kn)
            dens(k) = sigma(i,j,kn)
            sigr(k) = sigmar(i,j,k)
            nu(k) = difdia(i,j,k)
            rstdns(k) = .true.
            !diag
            ttem0(k) = ttem(k)
            ssal0(k) = ssal(k)
            delp0(k) = delp(k)
            dens0(k) = dens(k)
            sigr0(k) = sigr(k)
            nu0(k) = nu(k)
            !diag
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,k) = trc(i,j,kn,nt)
              end do
            end if
          end do

          ! --- --- Locate range of physical layers.
          kfpl = kfpla(i,j,n)
          kmin = kfpl-2
          kmax = 1
          do k = 2,kk
            if (delp(k) > epsilp) kmax = k
          end do

          if (kmin < kmax) then

            ! --- ----- Locate first layer where potential density should be kept
            ! --- ----- close to the reference potential density.
            rstdns(kfpl) = .false.
            if (kfpl /= kmax) then
              if (dens(kfpl) > .5*(sigr(kfpl)+sigr(kfpl+1))) then
                rstdns(kfpl+1) = .false.
              end if
            end if

            ! --- ----- Copy mixed layer variables to the layers with indexes kmin
            ! --- ----- and kmin+1
            delp(kmin+1) = delp(2)
            delp(kmin  ) = delp(1)
            ttem(kmin+1) = ttem(2)
            ttem(kmin  ) = ttem(1)
            ssal(kmin+1) = ssal(2)
            ssal(kmin  ) = ssal(1)
            nu(kmin+1) = nu(2)
            nu(kmin  ) = nu(1)
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,kmin+1) = ttrc(nt,2)
                ttrc(nt,kmin  ) = ttrc(nt,1)
              end do
            end if

            ! --- ----- Find interface pressure
            pres(kmin) = 0.
            do k = kmin,kmax
              pres(k+1) = pres(k)+delp(k)
            end do

            ! --- ----- Compute mass flux between the layers in the mixed layers and
            ! --- ----- at the mixed layer base.
            k = kmin
            fpu(k) = 0.
            fpl(k) = min(pres(k+1),pres(kmax+1)-pres(k+1), &
                 c*nu(k)*(delp(k)+delp(k+1)) &
                 /(2.*delp(k)*delp(k+1)))
            k = kmin+1
            fpu(k) = fpl(k-1)
            delpu = max(onem,delp(k))
            delpl = max(onem,delp(k+1))
            fpl(k) = min(pres(k+1),pres(kmax+1)-pres(k+1), &
                 c*nu(k)*(delpu+delpl)/(2.*delpu*delpl))
            fpl(kmax) = 0.

            if (kfpl <= kmax) then

              if (kfpl < kmax) then

                ! --- --------- Bottom boundary mixing is parameterized by assuming that
                ! --- --------- a part of the energy extracted from the mean flow by the
                ! --- --------- bottom drag drives diapycnal mixing.
                k = kmax-1
                nubbl = gbbl*ustarb(i,j)**3 &
                     *exp(-(delp(k+1)+.5*delp(k))*abs(coriop(i,j)) &
                     *alpha0/(kappa*max(ustmin,ustarb(i,j))*g)) &
                     /(alpha0*g*(sigr(k+1)-sigr(k)))
                nu(k) = max(nu(k),nubbl)
                difdia(i,j,k) = nu(k)
              end if

              ! --- ------- Compute linearized density jumps across upper and lower
              ! --- ------- interfaces of layers, average density jumps for layers,
              ! --- ------- and buoyancy flux corrections to nudge the layer density
              ! --- ------- towards the reference density. The flux corrections are
              ! --- ------- limited so their absolute value is guaranteed to be less
              ! --- ------- than the buoyancy flux due to diffusion.
              k = kfpl-1
              dsgli(k) = 1.
              fcl(k) = -fpl(k)
              do k = kfpl,kmax-1
                if (rstdns(k)) then
                  dsgdt = dsigdt(ttem(k),ssal(k))
                  dsgds = dsigds(ttem(k),ssal(k))
                  dsgu(k) = max(dsgmnr*(sigr(k)-sigr(k-1)), &
                       dsgdt*(ttem(k)-ttem(k-1)) &
                       +dsgds*(ssal(k)-ssal(k-1)))
                  dsgl(k) = max(dsgmnr*(sigr(k+1)-sigr(k)), &
                       dsgdt*(ttem(k+1)-ttem(k)) &
                       +dsgds*(ssal(k+1)-ssal(k)))
                  dsghm(k) = 2.*dsgu(k)*dsgl(k)/(dsgu(k)+dsgl(k))
                  dsg(k) = .5*(dsgu(k)+dsgl(k))
                  !                 dsg(k)=dsghm(k)
                  dsgui(k) = 1./dsgu(k)
                  dsgli(k) = 1./dsgl(k)
                  fcmx = .25*(sqrt(delp(k)*delp(k) &
                       +4.*c*nu(k)*dsg(k)*(dsgui(k)+dsgli(k))) &
                       -delp(k))*dsghm(k)*fcmxr
                  dsgc = dens(k)-sigr(k)
                  if (dsgc > 0.) then
                    fcl(k) = 0.
                    if (dens(k-1) < sigr(k)) then
                      q = max(0.,(dens(k)-sigr(k+1)) &
                           /((sigr(k)-sigr(k+1))*(1.-dsgcr0)))
                      q = max(0.,1.-q*q)
                      q = q*q*q
                      fcu(k) = dsgc*delp(k)
                      fcu(k) = min(q*fcu(k)+(1.-q)*fcmx,fcu(k))
                    else
                      fcu(k) = 0.
                    end if
                  else
                    fcu(k) = 0.
                    if (dens(k+1) > sigr(k)) then
                      q = max(0.,(dens(k)-sigr(k-1)) &
                           /((sigr(k)-sigr(k-1))*(1.-dsgcr0)))
                      q = max(0.,1.-q*q)
                      q = q*q*q
                      fcl(k) = dsgc*delp(k)
                      fcl(k) = max(q*fcl(k)-(1.-q)*fcmx,fcl(k))
                    else
                      fcl(k) = 0.
                    end if
                  end if
                else
                  dsgu(k) = 1.
                  dsgl(k) = 1.
                  dsghm(k) = 1.
                  dsg(k) = 1.
                  dsgui(k) = 1.
                  dsgli(k) = 1.
                  fcl(k) = 0.
                  fcu(k) = 0.
                end if
              end do
              k = kmax
              dsgdt = dsigdt(ttem(k),ssal(k))
              dsgds = dsigds(ttem(k),ssal(k))
              dsgu(k) = max(dsgmnr*(sigr(k)-sigr(k-1)), &
                   dsgdt*(ttem(k)-ttem(k-1)) &
                   +dsgds*(ssal(k)-ssal(k-1)))
              dsgui(k) = 1./dsgu(k)
              if (dens(k) > sigr(k).and.dens(k-1) < sigr(k)) then
                fpu(k) = min(delp(k-1), &
                     (dens(k)-sigr(k))*delp(k)*dsgui(k))
              else
                fpu(k) = 0.
              end if
              fcu(k) = fpu(k)*dsgu(k)

              ! --- ------- Find maximum allowable buoyancy fluxes and further limit
              ! --- ------- the flux corrections to keep layer interfaces within the
              ! --- ------- fluid domain.
              fmax(kfpl-1) = 0.
              fmax(kmax) = 0.
              done = .false.
              niter = 0
              do while (.not.done)
                done = .true.
                do k = kmax-1,kfpl,-1
                  q = ((fmax(k+1)+fcu(k+1))*dsgui(k+1) &
                       +pres(kmax+1)-pres(k+1))*dsgl(k)
                  fcl(k) = max(-q,fcl(k))
                  fmax(k) = q+fcl(k)
                end do
                kfmaxu = 0
                do k = kfpl,kmax-1
                  q = ((fmax(k-1)-fcl(k-1))*dsgli(k-1) &
                       +pres(k)-pres(kfpl))*dsgu(k)
                  if (fcu(k) > q) then
                    fcu(k) = q
                    done = .false.
                  end if
                  if (fmax(k) > q-fcu(k)) then
                    fmax(k) = q-fcu(k)
                    kfmaxu = k
                  end if
                end do
                if (niter == 100) then
                  write (lp,*) &
                       'blom: diapfl: no convergence in flux limit!', &
                       i+i0,j+j0
                  open (10,file='diapfl.uf',form = 'unformatted')
                  write (10) kk,kfpl
                  write (10) g,alpha0,epsilp,onem,delt1,dsgmnr,q,q
                  write (10) ttem0,ssal0,delp0,dens0,sigr0,nu0
                  close (10)
                  call xchalt('(diapfl)')
                  stop '(diapfl)'
                  exit
                end if
              end do

              ! --- ------- Make a first guess for buoyancy fluxes and set some
              ! --- ------- boundary conditions
              k = kfpl-1
              f0(k) = 0.
              f(k) = 0.
              gtd(k) = 0.
              dflim = 0.
              do k = kfpl,kmax-1
                f(k) = min(fmax(k), &
                     .5*sqrt(c*nu(k)*dsg(k) &
                     *(dsgui(k)+dsgli(k)))*dsghm(k), &
                     c*nu(k)*dsg(k)/max(epsilp,delp(k)))
                fold(k) = f(k)
                h(k) = fcu(k  )*dsgui(k  )-fcl(k  )*dsgli(k  ) &
                     +fcl(k-1)*dsgli(k-1)-fcu(k+1)*dsgui(k+1)
                dflim = max(dflim,fmax(k))
              end do
              k = kmax
              f0(k) = 0.
              f(k) = 0.
              gtd(k) = 0.
              dflim = dflim*dfeps

              ! --- ------- Solve the diffusion equation for layer thickness using
              ! --- ------- backward time integration by an iterative algorithm.
              niter = 0.
              dwnwrd = .false.
              do

                ! --- --------- Solve the equation by alternate downward and upward
                ! --- --------- propagation trough trough the layers
                dwnwrd = .not.dwnwrd
                if (dwnwrd) then

                  ! --- ----------- Do a downward first pass.
                  ctd = 0.
                  bitd = 1.
                  remfmx = .false.
                  do k = kfpl,kmax-1

                    if (remfmx) then
                      gtd(k) = 0.
                      f0(k) = fmax(k)
                      f(k) = fmax(k)
                    else

                      ! --- --------------- Find the backward solution of the layer buoyancy
                      ! --- --------------- flux, assuming no dependency on adjacent layer
                      ! --- --------------- fluxes, and the sensitivity of the layer flux to
                      ! --- --------------- changes in the fluxes of adjacent layers.
                      q = f0(k-1)*dsgli(k-1)+f(k+1)*dsgui(k+1) &
                           -delp(k)-h(k)
                      r = 4.*c*nu(k)*dsg(k)*(dsgui(k)+dsgli(k))
                      t = .25*dsghm(k)
                      if (q < 0.) then
                        s = r/(q*q)
                        if (s < 1.e-3) then

                          ! --- ------------------- For certain parameters, use a taylor expansion
                          ! --- ------------------- of the flux and flux sensitivity expressions
                          ! --- ------------------- to avoid roundoff errors.
                          r = .00390625*s
                          q = -q*r*(128.-s*(32.-s*(16.-s*(10.-s*7. ))))
                          f0(k) = q*t
                          q = r*(128.-s*(96.-s*(80.-s*(70.-s*63.))))
                          dfdg = q*t
                        else
                          s = sqrt(q*q+r)
                          f0(k) = (q+s)*t
                          dfdg = (1.+q/s)*t
                        end if
                      else
                        s = sqrt(q*q+r)
                        f0(k) = (q+s)*t
                        dfdg = (1.+q/s)*t
                      end if

                      if (f0(k) >= fmax(k)) then

                        ! --- ----------------- If the maximum flux associated with the lower
                        ! --- ----------------- fluid boundary has been reached, all subsequent
                        ! --- ----------------- fluxes will be set to the maximum allowable
                        ! --- ----------------- flux.
                        f0(k) = fmax(k)
                        dfdg = 0.
                        if (k > kfmaxu) remfmx = .true.
                      end if

                      ! --- --------------- Modify the buoyancy fluxes by taking into account
                      ! --- --------------- linearized contributions of adjacent layer fluxes.
                      ! --- --------------- This linearization forms a tridiagonal set of
                      ! --- --------------- equations.
                      gtd(k) = ctd*bitd
                      atd = -dfdg*dsgli(k-1)
                      ctd = -dfdg*dsgui(k+1)
                      bitd = 1./(1.-atd*gtd(k))
                      f(k) = (f0(k)-atd*(f(k-1)-f0(k-1))+ctd*f(k+1))*bitd
                    end if
                  end do

                  ! --- ----------- Complete the solving of the tridiagonal set of
                  ! --- ----------- equations and find the maximum flux change.
                  maxdf = 0.
                  do k = kmax-1,kfpl,-1
                    f(k) = min(fmax(k),f(k)-gtd(k+1)*f(k+1))
                    maxdf = max(maxdf,abs(f(k)-fold(k)))
                    fold(k) = f(k)
                  end do
                else

                  ! --- ----------- Do a upward first pass.
                  atd = 0.
                  bitd = 1.
                  remfmx = .false.
                  do k = kmax-1,kfpl,-1

                    if (remfmx) then
                      gtd(k) = 0.
                      f0(k) = fmax(k)
                      f(k) = fmax(k)
                    else

                      ! --- --------------- Find the backward solution of the layer buoyancy
                      ! --- --------------- flux, assuming no dependency on adjacent layer
                      ! --- --------------- fluxes, and the sensitivity of the layer flux to
                      ! --- --------------- changes in the fluxes of adjacent layers.
                      q = f(k-1)*dsgli(k-1)+f0(k+1)*dsgui(k+1) &
                           -delp(k)-h(k)
                      r = 4.*c*nu(k)*dsg(k)*(dsgui(k)+dsgli(k))
                      t = .25*dsghm(k)
                      if (q < 0.) then
                        s = r/(q*q)
                        if (s < 1.e-3) then

                          ! --- ------------------- For certain parameters, use a taylor expansion
                          ! --- ------------------- of the flux and flux sensitivity expressions
                          ! --- ------------------- to avoid roundoff errors.
                          r = .00390625*s
                          q = -q*r*(128.-s*(32.-s*(16.-s*(10.-s*7. ))))
                          f0(k) = q*t
                          q = r*(128.-s*(96.-s*(80.-s*(70.-s*63.))))
                          dfdg = q*t
                        else
                          s = sqrt(q*q+r)
                          f0(k) = (q+s)*t
                          dfdg = (1.+q/s)*t
                        end if
                      else
                        s = sqrt(q*q+r)
                        f0(k) = (q+s)*t
                        dfdg = (1.+q/s)*t
                      end if

                      if (f0(k) >= fmax(k)) then

                        ! --- ----------------- If the maximum flux associated with the upper
                        ! --- ----------------- fluid boundary has been reached, all subsequent
                        ! --- ----------------- fluxes will be set to the maximum allowable
                        ! --- ----------------- flux.
                        f0(k) = fmax(k)
                        dfdg = 0.
                        if (k <= kfmaxu) remfmx = .true.
                      end if

                      ! --- --------------- Modify the buoyancy fluxes by taking into account
                      ! --- --------------- linearized contributions of adjacent layer fluxes.
                      ! --- --------------- This linearization forms a tridiagonal set of
                      ! --- --------------- equations.
                      gtd(k) = atd*bitd
                      atd = -dfdg*dsgli(k-1)
                      ctd = -dfdg*dsgui(k+1)
                      bitd = 1./(1.-ctd*gtd(k))
                      f(k) = (f0(k)+atd*f(k-1)-ctd*(f(k+1)-f0(k+1)))*bitd
                    end if
                  end do

                  ! --- ----------- Complete the solving of the tridiagonal set of
                  ! --- ----------- equations and find the maximum flux change.
                  maxdf = 0.
                  do k = kfpl,kmax-1
                    f(k) = min(fmax(k),f(k)-gtd(k-1)*f(k-1))
                    maxdf = max(maxdf,abs(f(k)-fold(k)))
                    fold(k) = f(k)
                  end do
                end if

                ! --- --------- If the maximum flux change is below a treshold, stop the
                ! --- --------- iteration.
                niter = niter+1
                if (maxdf <= dflim) exit
                if (niter == 100) then
                  write (lp,*) &
                       'blom: diapfl: no convergence in implicit diffusion!', &
                       i+i0,j+j0,maxdf,dflim
                  open (10,file='diapfl.uf',form = 'unformatted')
                  write (10) kk,kfpl
                  write (10) g,alpha0,epsilp,onem,delt1,dsgmnr,q,q
                  write (10) ttem0,ssal0,delp0,dens0,sigr0,nu0
                  close (10)
                  call xchalt('(diapfl)')
                  stop '(diapfl)'
                  exit
                end if
              end do

              ! --- ------- Compute the mass fluxes
              do k = kfpl,kmax-1
                fpu(k) = (f(k)+fcu(k))*dsgui(k)
                fpl(k) = (f(k)-fcl(k))*dsgli(k)
              end do
              fpu(kfpl) = fpl(kmin+1)

            end if

            ! --- ----- Solve the diffusion equation for temperature and salinity by
            ! --- ----- backward integration forming a tridiagonal set of equations.
            ctd = 0.
            bitd = 1.
            do k = kmin,kmax
              gtd(k) = ctd*bitd
              q = 1./(delp(k)+fpu(k)+fpl(k))
              atd = -fpu(k)*q
              ctd = -fpl(k)*q
              dtd = delp(k)*q
              bitd = 1./(1.-atd*gtd(k))
              ssal(k) = (dtd*ssal(k)-atd*ssal(max(1,k-1)))*bitd
              ttem(k) = (dtd*ttem(k)-atd*ttem(max(1,k-1)))*bitd
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,k) = (dtd*ttrc(nt,k)-atd*ttrc(nt,max(1,k-1)))*bitd
                end do
              end if
            end do
            do k = kmax-1,kmin,-1
              ssal(k) = ssal(k)-gtd(k+1)*ssal(k+1)
              ttem(k) = ttem(k)-gtd(k+1)*ttem(k+1)
              dens(k) = sig(ttem(k),ssal(k))
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,k) = ttrc(nt,k)-gtd(k+1)*ttrc(nt,k+1)
                end do
              end if
            end do
            do k = kfpl,kmax-1
              delp(k) = max(0.,delp(k)+fpu(k)+fpl(k)-fpl(k-1)-fpu(k+1))
            end do
            delp(kmax) = max(0.,delp(kmax)+fpu(kmax)-fpl(kmax-1))

            ! --- ----- Copy variables back to the mixed layers from the layers with
            ! --- ----- index kmin and kmin+1
            ttem(1) = ttem(kmin  )
            ttem(2) = ttem(kmin+1)
            ssal(1) = ssal(kmin  )
            ssal(2) = ssal(kmin+1)
            dens(1) = dens(kmin  )
            dens(2) = dens(kmin+1)
            if (kmin > 1) then
              if (kmin == 2) then
                delp(2) = delp(kmin+1)
                delp(kmin+1) = 0.
              else
                delp(kmin  ) = 0.
              end if
            end if
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,1) = ttrc(nt,kmin  )
                ttrc(nt,2) = ttrc(nt,kmin+1)
              end do
            end if

          end if

          ! --- --- Fill massless layers with resonable values of temperature,
          ! --- --- salinity, and tracers
          if (kfpl > kmax) then
            do k = 3,kk
              ttem(k) = max(ttem(2),temmin(i,j,k))
              dens(k) = sigr(k)
              ssal(k) = sofsig(dens(k),ttem(k))
              delp(k) = 0.
              if (use_TRC) then
                do nt = 1,ntr
                  if (use_TKE) then
                    if (nt == itrtke) then
                      ttrc(nt,k) = max(ttrc(nt,2),tke_min)
                      cycle
                    end if
                    if (use_GLS) then
                      if (nt == itrgls) then
                        ttrc(nt,k) = max(ttrc(nt,2),gls_psi_min)
                        cycle
                      end if
                    end if
                  end if
                  ttrc(nt,k) = ttrc(nt,2)
                end do
              end if
            end do
          else
            do k = 3,kfpl-1
              ttem(k) = ttem(kfpl)
              dens(k) = sigr(k)
              ssal(k) = sofsig(dens(k),ttem(k))
              delp(k) = 0.
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,k) = ttrc(nt,kfpl)
                end do
              end if
            end do
            do k = kmax+1,kk
              ttem(k) = ttem(kmax)
              dens(k) = sigr(k)
              ssal(k) = sofsig(dens(k),ttem(k))
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,k) = ttrc(nt,kmax)
                end do
              end if
            end do
          end if

          ! --- --- Copy 1d arrays to 3d arrays
          do k = 1,kk
            kn = k+nn
            temp(i,j,kn) = ttem(k)
            saln(i,j,kn) = ssal(k)
            dp(i,j,kn) = delp(k)
            sigma(i,j,kn) = dens(k)
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
            if (use_TRC) then
              do nt = 1,ntr
                if (use_TKE) then
                  if (nt == itrtke) then
                    trc(i,j,kn,nt) = max(ttrc(nt,k),tke_min)
                    cycle
                  end if
                  if (use_GLS) then
                    if (nt == itrgls) then
                      trc(i,j,kn,nt) = max(ttrc(nt,k),gls_psi_min)
                      cycle
                    end if
                  end if
                end if
                trc(i,j,kn,nt) = ttrc(nt,k)
              end do
            end if
          end do

          ! --- --- Save variables used for momentum mixing
          kming(i,j) = kmin
          if (kmin < kmax) then
            do k = 1,kmin
              fpug(i,j,k) = fpl(kmin)
              fplg(i,j,k) = fpl(kmin)
            end do
            do k = kmin+1,kmax
              fpug(i,j,k) = fpu(k)
              fplg(i,j,k) = fpl(k)
            end do
            do k = kmax+1,kk
              fpug(i,j,k) = 0.
              fplg(i,j,k) = 0.
            end do
          else
            do k = 1,kk
              fpug(i,j,k) = 0.
              fplg(i,j,k) = 0.
            end do
          end if

        end do
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- Diapycnal mixing of momentum.
    ! --- ------------------------------------------------------------------

    call xctilr(p, 1,kk+1, 1,1, halo_ps)
    call xctilr(fpug, 1,kk, 1,1, halo_ps)
    call xctilr(fplg, 1,kk, 1,1, halo_ps)
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j) = kming(i,j)
        end do
      end do
    end do
    !$omp end parallel do
    call xctilr(util1, 1,1, 1,1, halo_ps)

    !$omp parallel do private(l,i)
    do j = 0,jj+1
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kming(i,j) = nint(util1(i,j))
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private( &
    !$omp l,i,kmin,kmax,k,uc,delp,kn,fpu,pold,pnew,fpum,fplm,fpup,fplp,fpl, &
    !$omp ctd,bitd,gtd,q,atd,dtd)
    do j = 1,jj

      ! --- - Mixing of u-component.

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))

          ! --- --- Find range of mass containing layers.
          kmin = min(kming(i-1,j),kming(i,j))
          kmax = 1
          do k = 2,kk
            if (dpu(i,j,k+nn) > 0.) kmax = k
          end do

          ! --- --- If number of mass containing layers is less than 2, do not do
          ! --- --- anything in this water column.
          if (kmin < kmax) then

            ! --- ----- Copy vertical column variables into 1d arrays.
            uc(kmin+1) = u(i,j,k1n+1)
            uc(kmin  ) = u(i,j,k1n  )
            delp(kmin+1) = dpu(i,j,k1n+1)
            delp(kmin  ) = dpu(i,j,k1n  )
            do k = kmin+2,kmax
              kn = k+nn
              uc(k) = u(i,j,kn)
              delp(k) = dpu(i,j,kn)
            end do

            ! --- ----- Interpolate interface mass fluxes to u-points. Limit mass
            ! --- ----- fluxes for correct reconstruction of layer thickness near
            ! --- ----- bathymetry.
            k = kmin
            fpu(k) = 0.
            do k = kmin+1,kmax
              pold = p(i-1,j,k)-fplg(i-1,j,k-1)+fpug(i-1,j,k)
              pnew = p(i-1,j,k)
              if (pold <= pu(i,j,kk+1)) then
                if (pnew <= pu(i,j,kk+1)) then
                  fpum = fpug(i-1,j,k  )
                  fplm = fplg(i-1,j,k-1)
                else
                  fpum = fpug(i-1,j,k  )
                  fplm = fplg(i-1,j,k-1)-pnew+pu(i,j,kk+1)
                end if
              else
                if (pnew <= pu(i,j,kk+1)) then
                  fpum = fpug(i-1,j,k  )-pold+pu(i,j,kk+1)
                  fplm = fplg(i-1,j,k-1)
                else
                  fpum = .5*(fpug(i-1,j,k  )+fplg(i-1,j,k-1))
                  fplm = fpum
                end if
              end if
              pold = p(i  ,j,k)-fplg(i  ,j,k-1)+fpug(i  ,j,k)
              pnew = p(i  ,j,k)
              if (pold <= pu(i,j,kk+1)) then
                if (pnew <= pu(i,j,kk+1)) then
                  fpup = fpug(i  ,j,k  )
                  fplp = fplg(i  ,j,k-1)
                else
                  fpup = fpug(i  ,j,k  )
                  fplp = fplg(i  ,j,k-1)-pnew+pu(i,j,kk+1)
                end if
              else
                if (pnew <= pu(i,j,kk+1)) then
                  fpup = fpug(i  ,j,k  )-pold+pu(i,j,kk+1)
                  fplp = fplg(i  ,j,k-1)
                else
                  fpup = .5*(fpug(i  ,j,k  )+fplg(i  ,j,k-1))
                  fplp = fpup
                end if
              end if
              fpu(k  ) = .5*(fpum+fpup)
              fpl(k-1) = .5*(fplm+fplp)
            end do
            fpl(kmax) = 0.

            ! --- ----- Solve the diffusion equation for velocity by backward
            ! --- ----- integration forming a tridiagonal set of equations.
            ctd = 0.
            bitd = 1.
            do k = kmin,kmax
              gtd(k) = ctd*bitd
              q = 1./(delp(k)+fpu(k)+fpl(k))
              atd = -fpu(k)*q
              ctd = -fpl(k)*q
              dtd = delp(k)*q
              bitd = 1./(1.-atd*gtd(k))
              uc(k) = (dtd*uc(k)-atd*uc(max(kmin,k-1)))*bitd
            end do
            do k = kmax-1,kmin,-1
              uc(k) = uc(k)-gtd(k+1)*uc(k+1)
            end do

            ! --- ----- Put velocity back in main array.
            u(i,j,k1n  ) = uc(kmin  )
            u(i,j,k1n+1) = uc(kmin+1)
            do k = kmin+2,kmax
              u(i,j,k+nn) = uc(k)
            end do

            ! --- ----- If interfaces are lifted above the bottom because of
            ! --- ----- diapycnal mixing, give the newly opened layers the velocity
            ! --- ----- of the lowest initial mass containing layer.
            do k = kmax+1,kk
              if (min(p(i-1,j,k),p(i,j,k)) < pu(i,j,kk+1)) &
                   u(i,j,k+nn) = uc(kmax)
            end do

          end if

        end do
      end do

      ! --- - Mixing of v-component.

      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))

          ! --- --- Find range of mass containing layers.
          kmin = min(kming(i,j-1),kming(i,j))
          kmax = 1
          do k = 2,kk
            if (dpv(i,j,k+nn) > 0.) kmax = k
          end do

          ! --- --- If number of mass containing layers is less than 2, do not do
          ! --- --- anything in this water column.
          if (kmin < kmax) then

            ! --- ----- Copy vertical column variables into 1d arrays.
            uc(kmin+1) = v(i,j,k1n+1)
            uc(kmin  ) = v(i,j,k1n  )
            delp(kmin+1) = dpv(i,j,k1n+1)
            delp(kmin  ) = dpv(i,j,k1n  )
            do k = kmin+2,kmax
              kn = k+nn
              uc(k) = v(i,j,kn)
              delp(k) = dpv(i,j,kn)
            end do

            ! --- ----- Interpolate interface mass fluxes to v-points. Limit mass
            ! --- ----- fluxes for correct reconstruction of layer thickness near
            ! --- ----- bathymetry.
            k = kmin
            fpu(k) = 0.
            do k = kmin+1,kmax
              pold = p(i,j-1,k)-fplg(i,j-1,k-1)+fpug(i,j-1,k)
              pnew = p(i,j-1,k)
              if (pold <= pv(i,j,kk+1)) then
                if (pnew <= pv(i,j,kk+1)) then
                  fpum = fpug(i,j-1,k  )
                  fplm = fplg(i,j-1,k-1)
                else
                  fpum = fpug(i,j-1,k  )
                  fplm = fplg(i,j-1,k-1)-pnew+pv(i,j,kk+1)
                end if
              else
                if (pnew <= pv(i,j,kk+1)) then
                  fpum = fpug(i,j-1,k  )-pold+pv(i,j,kk+1)
                  fplm = fplg(i,j-1,k-1)
                else
                  fpum = .5*(fpug(i,j-1,k  )+fplg(i,j-1,k-1))
                  fplm = fpum
                end if
              end if
              pold = p(i,j  ,k)-fplg(i,j  ,k-1)+fpug(i,j  ,k)
              pnew = p(i,j  ,k)
              if (pold <= pv(i,j,kk+1)) then
                if (pnew <= pv(i,j,kk+1)) then
                  fpup = fpug(i,j  ,k  )
                  fplp = fplg(i,j  ,k-1)
                else
                  fpup = fpug(i,j  ,k  )
                  fplp = fplg(i,j  ,k-1)-pnew+pv(i,j,kk+1)
                end if
              else
                if (pnew <= pv(i,j,kk+1)) then
                  fpup = fpug(i,j  ,k  )-pold+pv(i,j,kk+1)
                  fplp = fplg(i,j  ,k-1)
                else
                  fpup = .5*(fpug(i,j  ,k  )+fplg(i,j  ,k-1))
                  fplp = fpup
                end if
              end if
              fpu(k  ) = .5*(fpum+fpup)
              fpl(k-1) = .5*(fplm+fplp)
            end do
            fpl(kmax) = 0.

            ! --- ----- Solve the diffusion equation for velocity by backward
            ! --- ----- integration forming a tridiagonal set of equations.
            ctd = 0.
            bitd = 1.
            do k = kmin,kmax
              gtd(k) = ctd*bitd
              q = 1./(delp(k)+fpu(k)+fpl(k))
              atd = -fpu(k)*q
              ctd = -fpl(k)*q
              dtd = delp(k)*q
              bitd = 1./(1.-atd*gtd(k))
              uc(k) = (dtd*uc(k)-atd*uc(max(kmin,k-1)))*bitd
            end do
            do k = kmax-1,kmin,-1
              uc(k) = uc(k)-gtd(k+1)*uc(k+1)
            end do

            ! --- ----- Put velocity back in main array.
            v(i,j,k1n  ) = uc(kmin  )
            v(i,j,k1n+1) = uc(kmin+1)
            do k = kmin+2,kmax
              v(i,j,k+nn) = uc(k)
            end do

            ! --- ----- If interfaces are lifted above the bottom because of
            ! --- ----- diapycnal mixing, give the newly opened layers the velocity
            ! --- ----- of the lowest initial mass containing layer.
            do k = kmax+1,kk
              if (min(p(i,j-1,k),p(i,j,k)) < pv(i,j,kk+1)) &
                   v(i,j,k+nn) = uc(kmax)
            end do

          end if

        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,kn,l,i,q)
    do j = 1,jj
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            q = min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)= &
                 .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                 +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,kn,l,i,q)
    do j = 1,jj+1
      do k = 1,kk
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            q = min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)= &
                 .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                 +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !     do j=1,jj
    !       do l=1,isu(j)
    !       do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
    !         q=0.
    !         do k=1,kk
    !           kn=k+nn
    !           q=q+u(i,j,kn)*dpu(i,j,kn)
    !         enddo
    !         if (abs(q).gt.1.e-4) then
    !           write (lp,*) 'diapfl: u imbalance:',q,i,j
    !         endif
    !       enddo
    !       enddo
    !       do l=1,isv(j)
    !       do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
    !         q=0.
    !         do k=1,kk
    !           kn=k+nn
    !           q=q+v(i,j,kn)*dpv(i,j,kn)
    !         enddo
    !         if (abs(q).gt.1.e-4) then
    !           write (lp,*) 'diapfl: v imbalance:',q,i,j
    !         endif
    !       enddo
    !       enddo
    !     enddo

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'diapfl:'
      end if
      call chksummsk(p,ip,kk+1,'p')
      call chksummsk(dp,ip,2*kk,'dp')
      call chksummsk(temp,ip,2*kk,'temp')
      call chksummsk(saln,ip,2*kk,'saln')
      call chksummsk(sigma,ip,2*kk,'sigma')
      call chksummsk(u,iu,2*kk,'u')
      call chksummsk(v,iv,2*kk,'v')
      if (use_TRC) then
        do nt = 1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        end do
      end if
    end if

  end subroutine diapfl

end module mod_diapfl
