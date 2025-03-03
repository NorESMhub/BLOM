! ------------------------------------------------------------------------------
! Copyright (C) 2009-2025 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_mxlayr

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to the bulk
  ! mixed layer model.
  ! ------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm
  use mod_types,     only: r8
  use mod_constants, only: grav, spcifh, alpha0, epsilp, spval, onem, &
                           tencm, onecm, onemm, onemu
  use mod_time,      only: delt1
  use mod_xc,        only: xcstop, xctilr, isp, ifp, ilp, isu, &
                           ifu, ilu, isv, ifv, ip, ilv, nbdy, &
                           i0, j0, iu, iv, halo_ps, halo_uv, halo_vv
  use mod_vcoord,    only: sigmar
  use mod_grid,      only: scp2, scuxi, scvyi, coriop
  use mod_eos,       only: rho, sig, sig0, dsigdt, dsigdt0, &
                           dsigds, dsigds0, p_alpha,p_p_alpha
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                           p, pu, pv, kfpla
  use mod_swabs,     only: swbgal, swbgfc, swamxd
  use mod_forcing,   only: surflx, surrlx, sswflx, &
                           salflx, brnflx, salrlx, salt_corr, trc_corr, &
                           ustar, ustar3, buoyfl
  use mod_eddtra,    only: ce, lfmin, tau_mlr
  use mod_niw,       only: niwgf, niwbf, idkedt
  use mod_utility,   only: util1, util2, util3
  use mod_checksum,  only: csdiag, chksummsk
  use mod_tracers,   only: ntr, itrtke, itrgls, trc, trflx
  use mod_tke,       only: tke_min, gls_psi_min
  use mod_ifdefs,    only: use_TRC, use_TKE, use_GLS
  use mod_nctools

  implicit none
  private

  ! Variables to be set in namelist:
  real(r8) :: &
       rm0    ! Efficiency factor of wind TKE generation in
              ! the Oberhuber (1993) TKE closure [].
  real(r8) :: &
       rm5    ! Efficiency factor of TKE generation by
              ! momentum entrainment in the Oberhuber (1993)
              ! TKE closure [].
  character(len=80) :: &
       mlrttp ! Type of mixed layer restratification time
              ! scale. Valid types: 'variable', 'constant',
              ! 'limited'.

  ! Constants used in the bulk mixed layer model:
  real(r8) :: &
       mltmin = 5._r8  ! Minimum mixed-layer thickness [m].
  real(r8) :: &
       thktop = 10._r8 ! Thickness of top layer [m].

  ! Diagnostic variables:
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       mtkeus  ! Mixed layer TKE tendency related to friction
               ! velocity [cm3 s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       mtkeni  ! Mixed layer TKE tendency related to near
               ! inertial motions [cm3 s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       mtkebf  ! Mixed layer TKE tendency related to buoyancy
               ! forcing [cm3 s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       mtkers  ! Mixed layer TKE tendency related to eddy
               ! restratification [cm3 s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       mtkepe  ! Mixed layer TKE tendency related to pot.
               ! energy change [cm3 s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       mtkeke  ! Mixed layer TKE tendency related to kin.
               ! energy change [cm3 s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), protected :: &
       pbrnda  ! Brine plume pressure depth [g cm-1 s-2].

  ! Public module variables
  public :: rm0,rm5,mlrttp,mltmin
  public :: mtkeus,mtkeni,mtkebf,mtkers,mtkepe,mtkeke,pbrnda

  ! Public routines
  public :: inivar_mxlayr
  public :: mxlayr

contains

  subroutine inivar_mxlayr

    ! --- ------------------------------------------------------------------
    ! --- Initialize arrays.
    ! --- ------------------------------------------------------------------

    ! Local variables
    integer :: i,j

    !$omp parallel do private(i)
    do j = 1-nbdy,jj+nbdy
      do i = 1-nbdy,ii+nbdy
        mtkeus(i,j) = spval
        mtkeni(i,j) = spval
        mtkebf(i,j) = spval
        mtkers(i,j) = spval
        mtkepe(i,j) = spval
        mtkeke(i,j) = spval
        pbrnda(i,j) = spval
      end do
    end do
    !$omp end parallel do

  end subroutine inivar_mxlayr

  ! --- ----------------------------------------------------------------

  subroutine mxlayr(m,n,mm,nn,k1m,k1n)

    ! --- -----------------------------------------------------------------
    ! --- Modify mixed layer depth by enforcing a turbulent kinetic energy
    ! --- balance and apply surface forcing.
    ! --- -----------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(kdm+1) :: pres,po,pn
    real, dimension(kdm) :: ttem,ssal,delp,dens,densr,bc,uo,un
    real :: q,tmxl,smxl,rtau,rlf,alfa,beta,bfltot,bflpsw,pradd,pmxl, &
            lui,lbi,lei,cus,cni,cbftot,cbfpsw,crs,rm1,rm2,rm3,rm4,tkew, &
            dpmxl,tkeo,dtke,tfsl,sfsl,dpfsl,dptopl,dpt,pup,plo, &
            drhup,drhlo,pbrnd,bcwsum,bdpsum,tup,sup,dup,dsgdt,dsgds, &
            bpc,bpmldp,pswbas,pswup,pswlo,ttmp,stmp,sigtmp,sigfsl, &
            tmxl0,smxl0,dpe0,tdps,sdps,dpe,dps,udpn,um,vm,dke,dke0, &
            tkeu,tkel,uk,vk
    integer :: rtsflg,i,j,k,l,kn,kfpl,nitr,kmax,kfmax,ko
    logical :: chngd
    real, dimension(ntr,kdm) :: ttrc
    real, dimension(ntr) :: trfsl,trdps
    integer :: nt

    !  Parameters for Oberhuber (1993) TKE closure:
    !    mu     - parameter for the decay of TKE generated by surface
    !           - buoyancy flux when the buoyancy flux is destabilizing
    !           - [].
    !    kappa  - von Karman constant [].
    !    ustmin - minimum value of ustar used in computing the length
    !             scales for wind and buoyancy induced mixing [cm/s].
    !    mldjmp - minimum density jump at the mixed layer base used in
    !             the computation of potential energy change due to
    !             entrainment [g/cm**3].
    !    maxitr - maximum number of iterations allowed in the computation
    !             of TKE balance [].
    real :: kappa,mu,ustmin,mldjmp
    integer :: maxitr
    parameter (kappa=.4,mu=2.,ustmin = .001, &
               mldjmp=1.e-3,maxitr = 20)

    !  Parameters for the parameterization of restratification by mixed
    !  layer eddies by Fox-Kemper et al. (2008):
    !    cori20 - coriolis parameter at 20N [1/s].
    !    ci     - constant that appears when integrating the shape
    !             function over the mixed layer depth [].
    real :: cori20,ci,slbg0
    parameter (cori20 = 4.9745e-5,ci = 44./63.,slbg0 = 0.)

    !  Parameters for brine plume parameterization:
    !    bpdrho - density contrast between surface and brine plume depth
    !             [g/cm**3].
    !    bpmndp - minimum distribution thickness of salt from sea-ice
    !             freezing [g/cm/s**2].
    !    bpmxdp - maximum distribution depth below the mixed layer base
    !             of salt from sea-ice freezing [g/cm/s**2].
    !    bpdpmn - minimum layer thickness salt from sea-ice freezing
    !             is distributed over [g/cm/s**2].
    !    dsgmnr - minimum ratio of linearized density jump to target
    !             density jump across a layer interface [].
    real :: bpdrho,bpmndp,bpmxdp,bpdpmn,dsgmnr
    parameter (bpdrho=.4,bpmndp = 10.*onem, &
               bpmxdp=500.*onem,bpdpmn=1.*onem,dsgmnr = .1)

    ! ------------------------------------------------------------------
    ! Resolve type of mixed layer restratification time scale.
    ! ------------------------------------------------------------------

    if     (mlrttp == 'variable') then
      rtsflg = 1
    else if (mlrttp == 'constant') then
      rtsflg = 2
    else if (mlrttp == 'limited') then
      rtsflg = 3
    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' mlrttp = ',trim(mlrttp), &
             ' is unsupported!'
      end if
      call xcstop('(mxlayr)')
      stop '(mxlayr)'
    end if

    ! ------------------------------------------------------------------
    ! Compute squared lateral buoyancy gradient in the mixed layer and
    ! store it in -util1-.
    ! ------------------------------------------------------------------

    !$omp parallel do private(l,i,q,tmxl,smxl)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          q = 1./(dp(i,j,1+nn)+dp(i,j,2+nn))
          tmxl = (temp(i,j,1+nn)*dp(i,j,1+nn) &
                 +temp(i,j,2+nn)*dp(i,j,2+nn))*q
          smxl = (saln(i,j,1+nn)*dp(i,j,1+nn) &
                 +saln(i,j,2+nn)*dp(i,j,2+nn))*q
          util1(i,j) = grav*alpha0*sig0(tmxl,smxl)
        end do
      end do
    end do
    !$omp end parallel do
    call xctilr(util1, 1,1, 1,1, halo_ps)
    !$omp parallel do private(l,i,q)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          q = (util1(i,j)-util1(i-1,j))*scuxi(i,j)
          util2(i,j) = q*q
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i,q)
    do j = 1,jj+1
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          q = (util1(i,j)-util1(i,j-1))*scvyi(i,j)
          util3(i,j) = q*q
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          if     (ip(i-1,j)+ip(i+1,j) == 2) then
            util1(i,j) = .5*(util2(i,j)+util2(i+1,j))
          else if (ip(i-1,j) == 1) then
            util1(i,j) = util2(i  ,j)
          else if (ip(i+1,j) == 1) then
            util1(i,j) = util2(i+1,j)
          else
            util1(i,j) = 0.
          end if
          if     (ip(i,j-1)+ip(i,j+1) == 2) then
            util1(i,j) = util1(i,j)+.5*(util3(i,j)+util3(i,j+1))
          else if (ip(i,j-1) == 1) then
            util1(i,j) = util1(i,j)+util3(i,j)
          else if (ip(i,j+1) == 1) then
            util1(i,j) = util1(i,j)+util3(i,j+1)
          end if
          util1(i,j) = util1(i,j)+slbg0
        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------------
    ! Some mixed layer restratification parameterization constants.
    ! ------------------------------------------------------------------

    rtau = 1./tau_mlr
    rlf = 1./lfmin

    if (rm5 > 0.) then
      call xctilr(u(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
      call xctilr(v(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)
    end if

    !$omp parallel do private( &
    !$omp l,i,pres,k,kn,ttem,ssal,delp,dens,densr,q,tmxl,smxl,alfa,beta, &
    !$omp bfltot,bflpsw,pmxl,lui,lbi,lei,cus,cni,cbftot,cbfpsw,crs,rm1,rm2, &
    !$omp rm3,rm4,tkew,pradd,kfpl,dpmxl,tkeo,nitr,dtke,dpfsl,dptopl,tfsl, &
    !$omp sfsl,dpt,kmax,kfmax,bpmldp,pup,drhup,plo,drhlo,pbrnd,bcwsum, &
    !$omp bdpsum,tup,sup,dup,dsgdt,dsgds,bc,bpc,pswbas,pswup,pswlo,ttmp, &
    !$omp stmp,sigtmp,sigfsl,dps,tdps,sdps,tmxl0,smxl0,um,vm,dpe0,dke0, &
    !$omp tkeu,uk,vk,dpe,dke,tkel,chngd, &
    !$omp nt,ttrc,trfsl,trdps)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! Extract single column from 3-d fields.
          pres(1) = p(i,j,1)
          do k = 1,kk
            kn = k+nn
            ttem(k) = temp(i,j,kn)
            ssal(k) = saln(i,j,kn)
            delp(k) = dp(i,j,kn)
            pres(k+1) = pres(k)+delp(k)
            dens(k) = sigma(i,j,kn)
            densr(k) = sigmar(i,j,k)
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,k) = trc(i,j,kn,nt)
              end do
            end if
          end do

          ! ------------------------------------------------------------------
          ! Compute the turbulent kinetic energy (TKE) balance of the
          ! mixed layer.
          ! ------------------------------------------------------------------

          ! bfltot = total buoyancy flux [cm**2/sec**3]
          ! bflpsw = buoyancy flux due to penetrating short-wave radiation
          ! [cm**2/sec**3]
          ! note: surface density increases (column is destabilized) if
          ! bfltot > 0
          q = 1./(delp(1)+delp(2))
          tmxl = (ttem(1)*delp(1)+ttem(2)*delp(2))*q
          smxl = (ssal(1)*delp(1)+ssal(2)*delp(2))*q
          alfa = -alpha0*dsigdt0(tmxl,smxl)
          beta = alpha0*dsigds0(tmxl,smxl)
          bfltot = grav*alpha0*(alfa*surflx(i,j)/spcifh &
                            -beta*(salflx(i,j)-brnflx(i,j)))
          buoyfl(i,j,1) = bfltot
          bflpsw = grav*alpha0*alfa*swbgfc(i,j)*sswflx(i,j)/spcifh

          pmxl = pres(3)
          q = alpha0/grav
          lui = abs(coriop(i,j))*q/(kappa*max(ustmin,ustar(i,j)))
          lei = 1./(onem*swbgal(i,j))
          cus = rm0*ustar3(i,j)
          cni = niwgf*niwbf*idkedt(i,j)
          cbftot = .5*bfltot*q
          cbfpsw = .5*bflpsw*q
          if     (rtsflg == 1) then
            crs = ci*ce*util1(i,j)*q**3 &
                 *sqrt(scp2(i,j)/(coriop(i,j)*coriop(i,j)+rtau*rtau))*rlf
          else if (rtsflg == 2) then
            crs = ci*ce*util1(i,j)*q**3*sqrt(scp2(i,j))*rlf/cori20
          else if (rtsflg == 3) then
            crs = ci*ce*util1(i,j)*rlf*q**3 &
                 *sqrt(scp2(i,j)/min(cori20*cori20, &
                                     coriop(i,j)*coriop(i,j)+rtau*rtau))
          end if

          rm1 = exp(-lui*pmxl)
          q = lei*pmxl
          rm3 = exp(-q)
          rm4 = 2./q
          q = (cbftot-cbfpsw*(rm4*(1.-rm3)-rm3))
          if (q < 0.) then
            lbi = lui
            rm2 = rm1
          else
            lbi = lui*kappa/mu
            rm2 = exp(-lbi*pmxl)
          end if
          mtkeus(i,j) = cus*rm1
          mtkeni(i,j) = cni*rm1
          mtkebf(i,j) = q*rm2*pmxl
          mtkers(i,j) = -crs*pmxl*pmxl*pmxl
          mtkepe(i,j) = 0.
          mtkeke(i,j) = 0.
          tkew = mtkeus(i,j)+mtkeni(i,j)+mtkebf(i,j)+mtkers(i,j)

          pradd = swamxd*onem
          kfpl = kfpla(i,j,n)

          if (tkew < 0.and.pmxl > mltmin*onem) then

            ! -----------------------------------------------------------------
            ! If there is a TKE deficit in the mixed layer, reduce the
            ! mixed layer depth until the TKE balance is restored.
            ! -----------------------------------------------------------------

            if (pres(3)*lbi > 1.) then
              pmxl = 1./lbi
              dpmxl = min(pmxl-pres(1),pres(3)-pmxl,tencm)
              pmxl = pmxl-.5*dpmxl
            else
              dpmxl = -tencm
              pmxl = pres(3)+dpmxl
            end if
            tkeo = tkew
            nitr = 0
            do
              nitr = nitr+1
              rm1 = exp(-lui*pmxl)
              q = lei*max(tencm,pmxl)
              rm3 = exp(-q)
              rm4 = 2./q
              q = (cbftot-cbfpsw*(rm4*(1.-rm3)-rm3))
              if (q < 0.) then
                lbi = lui
                rm2 = rm1
              else
                lbi = lui*kappa/mu
                rm2 = exp(-lbi*pmxl)
              end if
              mtkeus(i,j) = cus*rm1
              mtkeni(i,j) = cni*rm1
              mtkebf(i,j) = q*rm2*pmxl
              mtkers(i,j) = -crs*pmxl*pmxl*pmxl
              mtkepe(i,j) = 0.
              mtkeke(i,j) = 0.
              tkew = mtkeus(i,j)+mtkeni(i,j)+mtkebf(i,j)+mtkers(i,j)
              if (.not.(nitr == 1.and.pres(3)*lbi > 1.)) then
                dtke = (tkew-tkeo)/dpmxl
                if (abs(dtke)<(abs(tkew)+1.e-22)/(pres(3)-pres(1))) then
                  if (tkew < 0.) then
                    dpmxl = .5*(pres(1)-pmxl)
                  else
                    dpmxl = .5*(pres(3)-pmxl)
                  end if
                else
                  dpmxl = max(pres(1)-pmxl,min(pres(3)-pmxl,-tkew/dtke))
                end if
              end if
              pmxl = pmxl+dpmxl
              tkeo = tkew
              if (abs(dpmxl) < onemm.or.nitr == maxitr) exit
            end do
            if (nitr == maxitr) then
              write (lp,*) 'reached maxitr when detraining',i+i0,j+j0
              write (lp,*) 'dpth = ',pres(3)/onem,';'
              write (lp,*) 'pmxl = ',pmxl/onem,';'
              write (lp,*) 'corio = ',coriop(i,j),';'
              write (lp,*) 'ustar = ',ustar(i,j),';'
              write (lp,*) 'bfltot = ',bfltot,';'
              write (lp,*) 'bflpsw = ',bflpsw,';'
              write (lp,*) 'bg2 = ',util1(i,j),';'
              write (lp,*) 'ce = ',ce*sqrt(scp2(i,j))*rlf,';'
              write (lp,*)
              ! call xchalt('(mxlayr)')
              ! stop '(mxlayr)'
            end if

            pmxl = max(mltmin*onem,pmxl)
            dpfsl = pres(3)-pmxl
            dptopl = min(thktop*onem,.5*(pmxl-pres(1)))

            if (pmxl < pres(2)) then

              q = 1./dpfsl
              tfsl = (ttem(2)*delp(2)+ttem(1)*(pres(2)-pmxl))*q
              sfsl = (ssal(2)*delp(2)+ssal(1)*(pres(2)-pmxl))*q
              ttem(2) = ttem(1)
              ssal(2) = ssal(1)
              if (use_TRC) then
                do nt = 1,ntr
                  trfsl(nt) = (ttrc(nt,2)*delp(2) &
                              +ttrc(nt,1)*(pres(2)-pmxl))*q
                  ttrc(nt,2) = ttrc(nt,1)
                end do
              end if
              delp(2) = pmxl-pres(1)-dptopl

            else

              tfsl = ttem(2)
              sfsl = ssal(2)
              if (use_TRC) then
                do nt = 1,ntr
                  trfsl(nt) = ttrc(nt,2)
                end do
              end if
              delp(2) = pmxl-pres(2)

              if (delp(1) > dptopl) then
                dpt = delp(1)-dptopl
                q = 1./(delp(2)+dpt)
                ttem(2) = (ttem(2)*delp(2)+ttem(1)*dpt)*q
                ssal(2) = (ssal(2)*delp(2)+ssal(1)*dpt)*q
                if (use_TRC) then
                  do nt = 1,ntr
                    ttrc(nt,2) = (ttrc(nt,2)*delp(2)+ttrc(nt,1)*dpt)*q
                  end do
                end if
                delp(2) = delp(2)+dpt
              else
                dpt = dptopl-delp(1)
                q = 1./(delp(1)+dpt)
                ttem(1) = (ttem(1)*delp(1)+ttem(2)*dpt)*q
                ssal(1) = (ssal(1)*delp(1)+ssal(2)*dpt)*q
                if (use_TRC) then
                  do nt = 1,ntr
                    ttrc(nt,1) = (ttrc(nt,1)*delp(1)+ttrc(nt,2)*dpt)*q
                  end do
                end if
                delp(2) = delp(2)-dpt
              end if

            end if

            delp(1) = dptopl

            ! ------------------------------------------------------------------
            ! Apply forcing to the water column
            ! ------------------------------------------------------------------

            kmax = 1
            do k = 2,kk
              if (delp(k) > epsilp) kmax = k
            end do
            kfmax = 0

            ! Apply brine forcing to the water column below surface layer
            pbrnda(i,j) = 0.
            if (brnflx(i,j) < 0.) then
              if (kfpl > kmax) then
                if (dpfsl > onemu) then
                  bpmldp = min(bpmndp,dpfsl+delp(2))
                  q = brnflx(i,j)*delt1*grav/bpmldp
                  ssal(2) = ssal(2)-q*max(0.,bpmldp-dpfsl)/delp(2)
                  sfsl = sfsl-q*min(dpfsl,bpmldp)/dpfsl
                else
                  ssal(2) = ssal(2)-brnflx(i,j)*delt1*grav/delp(2)
                end if
              else
                pup = pres(3)
                drhup = 0.
                k = kfpl
                do while (k <= kmax)
                  if (delp(k) > onemu) then
                    plo = pres(k)+.5*delp(k)
                    drhlo = rho(plo,ttem(k),ssal(k)) &
                         -rho(plo,ttem(1),ssal(1))
                    if (drhlo > bpdrho) exit
                    pup = plo
                    drhup = drhlo
                  end if
                  k = k+1
                end do
                if (k > kmax) then
                  pbrnd = pres(kmax+1)
                else
                  pbrnd = ((drhlo-bpdrho)*pup+(bpdrho-drhup)*plo) &
                       /(drhlo-drhup)
                end if
                pbrnd = min(pbrnd,pres(3)+bpmxdp)
                pbrnda(i,j) = pbrnd
                k = kfpl
                bcwsum = 0.
                bdpsum = 0.
                tup = tfsl
                sup = sfsl
                dup = sig(tfsl,sfsl)
                do while (k < kmax.and.pres(k+1) < pbrnd)
                  if (k == kfpl.or.dup < densr(k)) then
                    dsgdt = dsigdt(ttem(k),ssal(k))
                    dsgds = dsigds(ttem(k),ssal(k))
                    bc(k) = max(dsgmnr*(densr(k)-densr(k-1)), &
                                dsgdt*(ttem(k)-tup) &
                               +dsgds*(ssal(k)-sup)) &
                               /(dsgds*max(bpdpmn,delp(k)))
                    bcwsum = bcwsum+bc(k)*delp(k)
                    bdpsum = bdpsum+delp(k)
                  else
                    bc(k) = 0.
                  end if
                  tup = ttem(k)
                  sup = ssal(k)
                  dup = dens(k)
                  k = k+1
                end do
                if (k == kfpl.or.dup < densr(k)) then
                  dsgdt = dsigdt(ttem(k),ssal(k))
                  dsgds = dsigds(ttem(k),ssal(k))
                  bc(k) = max(dsgmnr*(densr(k)-densr(k-1)), &
                              dsgdt*(ttem(k)-tup)+dsgds*(ssal(k)-sup)) &
                         *max(bpdpmn,&
                              pbrnd-pres(k))/(dsgds*max(bpdpmn,delp(k))**2)
                  bcwsum = bcwsum+bc(k)*delp(k)
                  bdpsum = bdpsum+delp(k)
                else
                  bc(k) = 0.
                end if
                kfmax = k
                if (bdpsum <= epsilp) then
                  if (dpfsl > onemu) then
                    bpmldp = min(bpmndp,dpfsl+delp(2))
                    q = brnflx(i,j)*delt1*grav/bpmldp
                    ssal(2) = ssal(2)-q*max(0.,bpmldp-dpfsl)/delp(2)
                    sfsl = sfsl-q*min(dpfsl,bpmldp)/dpfsl
                  else
                    ssal(2) = ssal(2)-brnflx(i,j)*delt1*grav/delp(2)
                  end if
                else
                  if (bdpsum < bpmndp) then
                    bpmldp = min(bpmndp,bdpsum+dpfsl+delp(2))
                    q = brnflx(i,j)*delt1*grav/bpmldp
                    ssal(2) = ssal(2) &
                         -q*max(0.,bpmldp-bdpsum-dpfsl)/delp(2)
                    if (dpfsl > onemu) then
                      sfsl = sfsl-q*min(dpfsl,bpmldp-bdpsum)/dpfsl
                      bpc = q*bdpsum/bcwsum
                    else
                      bpc = q*(bdpsum+dpfsl)/bcwsum
                    end if
                  else
                    bpc = brnflx(i,j)*delt1*grav/bcwsum
                  end if
                  do k = kfpl,kfmax
                    ssal(k) = ssal(k)-bpc*bc(k)
                  end do
                end if
              end if
            end if

            ! Apply heat forcing to the water column below surface layer
            pswbas = swbgfc(i,j)*exp(-lei*delp(1))
            pswup = pswbas
            pswlo = swbgfc(i,j)*exp(-lei*min(pradd,pmxl))
            q = delt1*grav/delp(2)
            ttem(2) = ttem(2)-(pswup-pswlo)*sswflx(i,j)*q/spcifh
            pswup = pswlo
            pswlo = swbgfc(i,j)*exp(-lei*min(pradd,pres(3)))
            if (dpfsl > onemu) then
              tfsl = tfsl-(pswup-pswlo)*sswflx(i,j)*delt1*grav/(spcifh*dpfsl)
              pswup = pswlo
            end if
            k = kfpl
            do while (k < kmax)
              if (delp(k) > onemu) then
                pswlo = swbgfc(i,j)*exp(-lei*min(pradd,pres(k+1)))
                ttem(k) = ttem(k)-(pswup-pswlo)*sswflx(i,j)*delt1*grav/(spcifh*delp(k))
                pswup = pswlo
                kfmax = max(kfmax,k)
              end if
              k = k+1
              if (pres(k) > pradd) exit
            end do

            ! Apply heat and salt forcing to top layer
            q = delt1*grav/delp(1)
            ttem(1) = ttem(1) &
                 -(surflx(i,j)-(pswbas-pswup)*sswflx(i,j) &
                 +surrlx(i,j))*q/spcifh
            ssal(1) = ssal(1) &
                 -(salflx(i,j)-brnflx(i,j) &
                 +salrlx(i,j))*q
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,1) = ttrc(nt,1)-trflx(nt,i,j)*q
              end do
            end if

            ! Update density for layers where forcing has been applied
            dens(1) = sig(ttem(1),ssal(1))
            dens(2) = sig(ttem(2),ssal(2))
            do k = kfpl,kfmax
              dens(k) = sig(ttem(k),ssal(k))
            end do

            if (dpfsl <= onemu) then
              q = 1./(dpfsl+delp(2))
              ttem(2) = (tfsl*dpfsl+ttem(2)*delp(2))*q
              ssal(2) = (sfsl*dpfsl+ssal(2)*delp(2))*q
              dens(2) = sig(ttem(2),ssal(2))
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,2) = (trfsl(nt)*dpfsl+ttrc(nt,2)*delp(2))*q
                end do
              end if
              delp(2) = dpfsl+delp(2)
            else

              ! Place the content of the fossil mixed layer in isopycnic
              ! layers.
              k = min(kk,kfpl)

              if (k == 3) then
                q = 1./(dpfsl+delp(k))
                ttem(k) = (tfsl*dpfsl+ttem(k)*delp(k))*q
                ssal(k) = (sfsl*dpfsl+ssal(k)*delp(k))*q
                dens(k) = sig(ttem(k),ssal(k))
                if (use_TRC) then
                  do nt = 1,ntr
                    ttrc(nt,k) = (trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                  end do
                end if
                delp(k) = dpfsl+delp(k)
              else
                q = 1./(dpfsl+delp(k))
                ttmp = (tfsl*dpfsl+ttem(k)*delp(k))*q
                stmp = (sfsl*dpfsl+ssal(k)*delp(k))*q
                sigtmp = sig(ttmp,stmp)
                sigfsl = sig(tfsl,sfsl)
                if (sigtmp >= densr(k)) then
                  if (sigfsl > dens(k).and. &
                       dens(k) <= dens(min(kk,k+1)).and. &
                       rho(pmxl,tfsl   ,sfsl   ) < &
                       rho(pmxl,ttem(k),ssal(k))) then
                    k = k-1
                    q = 1./(dpfsl+delp(k))
                    ttem(k) = (tfsl*dpfsl+ttem(k)*delp(k))*q
                    ssal(k) = (sfsl*dpfsl+ssal(k)*delp(k))*q
                    dens(k) = sig(ttem(k),ssal(k))
                    if (use_TRC) then
                      do nt = 1,ntr
                        ttrc(nt,k) = (trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                      end do
                    end if
                    delp(k) = dpfsl+delp(k)
                  else
                    ttem(k) = ttmp
                    ssal(k) = stmp
                    dens(k) = sigtmp
                    if (use_TRC) then
                      do nt = 1,ntr
                        ttrc(nt,k) = (trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                      end do
                    end if
                    delp(k) = dpfsl+delp(k)
                  end if
                else
                  if (delp(k) > onemu.and.dens(k) > densr(k).and. &
                      sigfsl < densr(k)-(1.e-6)) then
                    dps = min(dpfsl,&
                              delp(k)*(dens(k)-densr(k))/(densr(k)-sigfsl))
                    q = 1./(dps+delp(k))
                    ttem(k) = (tfsl*dps+ttem(k)*delp(k))*q
                    ssal(k) = (sfsl*dps+ssal(k)*delp(k))*q
                    dens(k) = sig(ttem(k),ssal(k))
                    if (use_TRC) then
                      do nt = 1,ntr
                        ttrc(nt,k) = (trfsl(nt)*dps+ttrc(nt,k)*delp(k))*q
                      end do
                    end if
                    delp(k) = dps+delp(k)
                    dpfsl = dpfsl-dps
                    if (dpfsl <= onemu) then
                      q = 1./(dpfsl+delp(2))
                      ttem(2) = (tfsl*dpfsl+ttem(2)*delp(2))*q
                      ssal(2) = (sfsl*dpfsl+ssal(2)*delp(2))*q
                      dens(2) = sig(ttem(2),ssal(2))
                      if (use_TRC) then
                        do nt = 1,ntr
                          ttrc(nt,2) = (trfsl(nt)*dpfsl + ttrc(nt,2)*delp(2))*q
                        end do
                      end if
                      delp(2) = dpfsl+delp(2)
                    else
                      k = k-1
                      do while (sigfsl < densr(k))
                        if (k == 3) exit
                        k = k-1
                      end do
                      q = 1./(dpfsl+delp(k))
                      ttem(k) = (tfsl*dpfsl+ttem(k)*delp(k))*q
                      ssal(k) = (sfsl*dpfsl+ssal(k)*delp(k))*q
                      dens(k) = sig(ttem(k),ssal(k))
                      if (use_TRC) then
                        do nt = 1,ntr
                          ttrc(nt,k) = (trfsl(nt)*dpfsl + ttrc(nt,k)*delp(k))*q
                        end do
                      end if
                      delp(k) = dpfsl+delp(k)
                    end if
                  else
                    k = k-1
                    do while (sigfsl < densr(k))
                      if (k == 3) exit
                      k = k-1
                    end do
                    q = 1./(dpfsl+delp(k))
                    ttem(k) = (tfsl*dpfsl+ttem(k)*delp(k))*q
                    ssal(k) = (sfsl*dpfsl+ssal(k)*delp(k))*q
                    dens(k) = sig(ttem(k),ssal(k))
                    if (use_TRC) then
                      do nt = 1,ntr
                        ttrc(nt,k) = (trfsl(nt)*dpfsl+ttrc(nt,k)*delp(k))*q
                      end do
                    end if
                    delp(k) = dpfsl+delp(k)
                  end if
                end if
              end if

            end if

          else

            if (tkew < 0.) then
              pmxl = mltmin*onem
              tdps = ttem(2)*delp(2)
              sdps = ssal(2)*delp(2)
              if (use_TRC) then
                do nt = 1,ntr
                  trdps(nt) = ttrc(nt,2)*delp(2)
                end do
              end if
              k = kfpl
              do while (k <= kk)
                q = min(pmxl,pres(k+1))-pres(k)
                tdps = tdps+ttem(k)*q
                sdps = sdps+ssal(k)*q
                if (use_TRC) then
                  do nt = 1,ntr
                    trdps(nt) = trdps(nt)+ttrc(nt,k)*q
                  end do
                end if
                delp(k) = pres(k+1)-min(pmxl,pres(k+1))
                if (pres(k+1) > pmxl) exit
                k = k+1
              end do
            else

              ! -----------------------------------------------------------------
              ! If there is a TKE surplus in the mixed layer, increase
              ! the mixed layer depth until the TKE balance is restored.
              ! -----------------------------------------------------------------

              q = 1./(delp(1)+delp(2))
              tmxl0 = (ttem(1)*delp(1)+ttem(2)*delp(2))*q
              smxl0 = (ssal(1)*delp(1)+ssal(2)*delp(2))*q
              um = (u(i  ,j,1+nn)*dpu(i  ,j,1+nn) &
                   +u(i+1,j,1+nn)*dpu(i+1,j,1+nn) &
                   +u(i  ,j,2+nn)*dpu(i  ,j,2+nn) &
                   +u(i+1,j,2+nn)*dpu(i+1,j,2+nn)) &
                   /max(onecm,&
                        dpu(i  ,j,1+nn)+dpu(i+1,j,1+nn) &
                       +dpu(i  ,j,2+nn)+dpu(i+1,j,2+nn))
              vm = (v(i,j  ,1+nn)*dpv(i,j  ,1+nn) &
                   +v(i,j+1,1+nn)*dpv(i,j+1,1+nn) &
                   +v(i,j  ,2+nn)*dpv(i,j  ,2+nn) &
                   +v(i,j+1,2+nn)*dpv(i,j+1,2+nn)) &
                   /max(onecm,&
                        dpv(i,j,1+nn)+dpv(i,j+1,1+nn) &
                       +dpv(i,j,2+nn)+dpv(i,j+1,2+nn))
              dpe0 = 0.
              dke0 = 0.
              tkeu = tkew
              k = kfpl
              tdps = ttem(2)*delp(2)
              sdps = ssal(2)*delp(2)
              if (use_TRC) then
                do nt = 1,ntr
                  trdps(nt) = ttrc(nt,2)*delp(2)
                end do
              end if
              do
                if (k > kk) then
                  exit
                else if (delp(k) < epsilp) then
                  k = k+1
                else
                  pmxl = pres(k+1)
                  uk = (u(i  ,j,k+nn)*dpu(i  ,j,k+nn) &
                       +u(i+1,j,k+nn)*dpu(i+1,j,k+nn)) &
                       /max(onecm,dpu(i  ,j,k+nn)+dpu(i+1,j,k+nn))
                  vk = (v(i,j  ,k+nn)*dpv(i,j  ,k+nn) &
                       +v(i,j+1,k+nn)*dpv(i,j+1,k+nn)) &
                       /max(onecm,dpv(i,j  ,k+nn)+dpv(i,j+1,k+nn))
                  nitr = 0
                  do
                    nitr = nitr+1
                    tmxl = (tmxl0*(pres(k)-pres(1)) &
                         +ttem(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    smxl = (smxl0*(pres(k)-pres(1)) &
                         +ssal(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    dpe = dpe0 &
                         +max(.5*alpha0*alpha0*mldjmp &
                        *(pres(k)-pres(1))*(pmxl-pres(k)), &
                          p_p_alpha(pmxl,pres(1),tmxl,smxl) &
                         -p_p_alpha(pmxl,pres(k),ttem(k),ssal(k)) &
                         -p_p_alpha(pres(k),pres(1),tmxl0,smxl0) &
                         -(pres(1)-pres(k)) &
                         *p_alpha(pmxl,pres(k),ttem(k),ssal(k))) &
                         *alpha0/(delt1*grav)
                    dke = dke0 &
                         +.5*rm5*(pres(k)-pres(1))*(pmxl-pres(k)) &
                         *((uk-um)**2+(vk-vm)**2)*alpha0 &
                         /((pmxl-pres(1))*delt1*grav)
                    rm1 = exp(-lui*pmxl)
                    q = lei*pmxl
                    rm3 = exp(-q)
                    rm4 = 2./q
                    q = (cbftot-cbfpsw*(rm4*(1.-rm3)-rm3))
                    if (q < 0.) then
                      lbi = lui
                      rm2 = rm1
                    else
                      lbi = lui*kappa/mu
                      rm2 = exp(-lbi*pmxl)
                    end if
                    mtkeus(i,j) = cus*rm1
                    mtkeni(i,j) = cni*rm1
                    mtkebf(i,j) = q*rm2*pmxl
                    mtkers(i,j) = -crs*pmxl*pmxl*pmxl
                    mtkepe(i,j) = -dpe
                    mtkeke(i,j) = dke
                    tkew = mtkeus(i,j)+mtkeni(i,j)+mtkebf(i,j)+mtkers(i,j) &
                          +mtkepe(i,j)+mtkeke(i,j)
                    if (nitr == 1) then
                      if (tkew > 0.) then
                        exit
                      else
                        pmxl = pres(k)
                        dpmxl = min(tencm,.5*delp(k))
                        tkel = tkew
                        tkew = tkeu
                      end if
                    else
                      dtke = (tkew-tkeo)/dpmxl
                      chngd = .false.
                      if (nitr == 2) then
                        if (dtke > -tkew/(pres(k+1)-pmxl)) then
                          pmxl = pres(k+1)
                          dpmxl = -min(tencm,.5*delp(k))
                          tkew = tkel
                          chngd = .true.
                        end if
                      end if
                      if (.not.chngd) then
                        if (abs(dtke) < &
                             (abs(tkew)+1.e-22)/delp(k)) then
                          if (tkew < 0.) then
                            dpmxl = .5*(pres(k)-pmxl)
                          else
                            dpmxl = pres(k+1)-pmxl
                          end if
                        else
                          dpmxl = max(pres(k)-pmxl, &
                                  min(pres(k+1)-pmxl,-tkew/dtke))
                        end if
                        dpmxl = max(max(mltmin*onem,pres(k))-pmxl,dpmxl)
                      end if
                    end if
                    pmxl = pmxl+dpmxl
                    tkeo = tkew
                    if (abs(dpmxl) < onemm.or.nitr == maxitr) exit
                  end do
                  if (nitr == maxitr) then
                    write (lp,*) 'reached maxitr when entraining',i+i0,j+j0
                    write (lp,*) 'dpth = ',pres(3)/onem,';'
                    write (lp,*) 'pmxl = ',pmxl/onem,';'
                    write (lp,*) 'corio = ',coriop(i,j),';'
                    write (lp,*) 'ustar = ',ustar(i,j),';'
                    write (lp,*) 'bfltot = ',bfltot,';'
                    write (lp,*) 'bflpsw = ',bflpsw,';'
                    write (lp,*) 'bg2 = ',util1(i,j),';'
                    write (lp,*) 'ce = ',ce*sqrt(scp2(i,j))*rlf,';'
                    write (lp,*) 'pres(3) = ',pres(3)/onem,';'
                    do kn = kfpla(i,j,n),k
                      write (lp,*) 'pres(',kn+1,') = ',pres(kn+1)/onem,';'
                    end do
                    write (lp,*) 'ttem(1) = ',ttem(1),';'
                    do kn = kfpla(i,j,n),k
                      write (lp,*) 'ttem(',kn,') = ',ttem(kn),';'
                    end do
                    write (lp,*) 'ssal(1) = ',ssal(1),';'
                    do kn = kfpla(i,j,n),k
                      write (lp,*) 'ssal(',kn,') = ',ssal(kn),';'
                    end do
                    ! call xchalt('(mxlayr)')
                    ! stop '(mxlayr)'
                  end if
                  if (pmxl < pres(k+1)-epsilp.and.nitr < maxitr) then
                    tdps = tdps+ttem(k)*(pmxl-pres(k))
                    sdps = sdps+ssal(k)*(pmxl-pres(k))
                    if (use_TRC) then
                      do nt = 1,ntr
                        trdps(nt) = trdps(nt)+ttrc(nt,k)*(pmxl-pres(k))
                      end do
                    end if
                    delp(k) = pres(k+1)-pmxl
                    exit
                  else
                    tdps = tdps+ttem(k)*delp(k)
                    sdps = sdps+ssal(k)*delp(k)
                    if (use_TRC) then
                      do nt = 1,ntr
                        trdps(nt) = trdps(nt)+ttrc(nt,k)*delp(k)
                      end do
                    end if
                    pmxl = pres(k+1)
                    tmxl = (tmxl0*(pres(k)-pres(1)) &
                         +ttem(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    smxl = (smxl0*(pres(k)-pres(1)) &
                         +ssal(k)*(pmxl-pres(k)))/(pmxl-pres(1))
                    dpe = dpe0 &
                         +max(.5*alpha0*alpha0*mldjmp &
                         *(pres(k)-pres(1))*(pmxl-pres(k)), &
                          p_p_alpha(pmxl,pres(1),tmxl,smxl) &
                         -p_p_alpha(pmxl,pres(k),ttem(k),ssal(k)) &
                         -p_p_alpha(pres(k),pres(1),tmxl0,smxl0) &
                         -(pres(1)-pres(k)) &
                         *p_alpha(pmxl,pres(k),ttem(k),ssal(k)))*alpha0/(delt1*grav)
                    dpe0 = dpe
                    dke = dke0 &
                         +.5*rm5*(pres(k)-pres(1))*(pmxl-pres(k)) &
                         *((uk-um)**2+(vk-vm)**2)*alpha0/((pmxl-pres(1))*delt1*grav)
                    dke0 = dke
                    tmxl0 = tmxl
                    smxl0 = smxl
                    um = (um*(pres(k)-pres(1))+uk*(pmxl-pres(k)))/(pmxl-pres(1))
                    vm = (vm*(pres(k)-pres(1))+vk*(pmxl-pres(k)))/(pmxl-pres(1))
                    delp(k) = 0.
                    k = k+1
                  end if
                end if
              end do
            end if

            pres(3) = min(pres(kk+1),pmxl)
            delp(2) = pres(3)-pres(2)
            q = 1./delp(2)
            ttem(2) = tdps*q
            ssal(2) = sdps*q
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,2) = trdps(nt)*q
              end do
            end if
            kfpl = k
            do k = 4,kfpl
              pres(k) = pres(3)
            end do

            ! Restore top layer to its reference pressure thickness
            dptopl = min(thktop*onem,.5*(pres(3)-pres(1)))
            if (delp(1) > dptopl) then
              dpt = delp(1)-dptopl
              q = 1./(delp(2)+dpt)
              ttem(2) = (ttem(2)*delp(2)+ttem(1)*dpt)*q
              ssal(2) = (ssal(2)*delp(2)+ssal(1)*dpt)*q
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,2) = (ttrc(nt,2)*delp(2)+ttrc(nt,1)*dpt)*q
                end do
              end if
              delp(2) = delp(2)+dpt
            else
              dpt = dptopl-delp(1)
              q = 1./(delp(1)+dpt)
              ttem(1) = (ttem(1)*delp(1)+ttem(2)*dpt)*q
              ssal(1) = (ssal(1)*delp(1)+ssal(2)*dpt)*q
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,1) = (ttrc(nt,1)*delp(1)+ttrc(nt,2)*dpt)*q
                end do
              end if
              delp(2) = delp(2)-dpt
            end if
            delp(1) = dptopl
            pres(2) = pres(1)+delp(1)

            ! ------------------------------------------------------------------
            ! Apply forcing to the water column
            ! ------------------------------------------------------------------

            kmax = 1
            do k = 2,kk
              if (delp(k) > epsilp) kmax = k
            end do
            kfmax = 0

            ! Apply brine forcing to the water column below surface layer
            pbrnda(i,j) = 0.
            if (brnflx(i,j) < 0.) then
              if (kfpl > kmax) then
                ssal(2) = ssal(2)-brnflx(i,j)*delt1*grav/delp(2)
              else
                pup = pres(3)
                drhup = 0.
                k = kfpl
                do while (k <= kmax)
                  if (delp(k) > onemu) then
                    plo = pres(k)+.5*delp(k)
                    drhlo = rho(plo,ttem(k),ssal(k)) &
                           -rho(plo,ttem(1),ssal(1))
                    if (drhlo > bpdrho) exit
                    pup = plo
                    drhup = drhlo
                  end if
                  k = k+1
                end do
                if (k > kmax) then
                  pbrnd = pres(kmax+1)
                else
                  pbrnd = ((drhlo-bpdrho)*pup+(bpdrho-drhup)*plo) &
                          /(drhlo-drhup)
                end if
                pbrnd = min(pbrnd,pres(3)+bpmxdp)
                pbrnda(i,j) = pbrnd
                k = kfpl
                bcwsum = 0.
                bdpsum = 0.
                tup = ttem(2)
                sup = ssal(2)
                dup = sig(ttem(2),ssal(2))
                do while (k < kmax.and.pres(k+1) < pbrnd)
                  if (k == kfpl.or.dup < densr(k)) then
                    dsgdt = dsigdt(ttem(k),ssal(k))
                    dsgds = dsigds(ttem(k),ssal(k))
                    bc(k) = max(dsgmnr*(densr(k)-densr(k-1)), &
                                dsgdt*(ttem(k)-tup)+dsgds*(ssal(k)-sup)) &
                               /(dsgds*max(bpdpmn,delp(k)))
                    bcwsum = bcwsum+bc(k)*delp(k)
                    bdpsum = bdpsum+delp(k)
                  else
                    bc(k) = 0.
                  end if
                  tup = ttem(k)
                  sup = ssal(k)
                  dup = dens(k)
                  k = k+1
                end do
                if (k == kfpl.or.dup < densr(k)) then
                  dsgdt = dsigdt(ttem(k),ssal(k))
                  dsgds = dsigds(ttem(k),ssal(k))
                  bc(k) = max(dsgmnr*(densr(k)-densr(k-1)), &
                              dsgdt*(ttem(k)-tup) &
                              +dsgds*(ssal(k)-sup)) &
                              *max(bpdpmn,pbrnd-pres(k)) &
                                   /(dsgds*max(bpdpmn,delp(k))**2)
                  bcwsum = bcwsum+bc(k)*delp(k)
                  bdpsum = bdpsum+delp(k)
                else
                  bc(k) = 0.
                end if
                kfmax = k
                if (bdpsum <= epsilp) then
                  ssal(2) = ssal(2)-brnflx(i,j)*delt1*grav/delp(2)
                else
                  if (bdpsum < bpmndp) then
                    bpmldp = min(bpmndp,bdpsum+delp(2))
                    q = brnflx(i,j)*delt1*grav/bpmldp
                    ssal(2) = ssal(2)-q*(bpmldp-bdpsum)/delp(2)
                    bpc = q*bdpsum/bcwsum
                  else
                    bpc = brnflx(i,j)*delt1*grav/bcwsum
                  end if
                  do k = kfpl,kfmax
                    ssal(k) = ssal(k)-bpc*bc(k)
                  end do
                end if
              end if
            end if

            ! Apply heat forcing to the water column below surface layer
            pswbas = swbgfc(i,j)*exp(-lei*delp(1))
            pswup = pswbas
            pswlo = swbgfc(i,j)*exp(-lei*min(pradd,pres(3)))
            q = delt1*grav/delp(2)
            ttem(2) = ttem(2)-(pswup-pswlo)*sswflx(i,j)*q/spcifh
            pswup = pswlo
            k = kfpl
            do while (k < kmax)
              if (delp(k) > onemu) then
                pswlo = swbgfc(i,j)*exp(-lei*min(pradd,pres(k+1)))
                ttem(k) = ttem(k)-(pswup-pswlo)*sswflx(i,j)*delt1*grav/(spcifh*delp(k))
                pswup = pswlo
                kfmax = max(kfmax,k)
              end if
              k = k+1
              if (pres(k) > pradd) exit
            end do

            ! Apply heat and salt forcing to top layer
            q = delt1*grav/delp(1)
            ttem(1) = ttem(1) &
                 -(surflx(i,j)-(pswbas-pswup)*sswflx(i,j)+surrlx(i,j))*q/spcifh
            ssal(1) = ssal(1) - (salflx(i,j)-brnflx(i,j) + salrlx(i,j))*q
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,1) = ttrc(nt,1)-trflx(nt,i,j)*q
              end do
            end if

            ! Update density for layers where forcing has been applied
            dens(1) = sig(ttem(1),ssal(1))
            dens(2) = sig(ttem(2),ssal(2))
            do k = kfpl,kfmax
              dens(k) = sig(ttem(k),ssal(k))
            end do

          end if

          ! Define first physical layer.
          k = 3
          dps = 0.
          do while (delp(k) < epsilp)
            dps = dps+delp(k)
            delp(k) = 0.
            k = k+1
            if (k > kk) exit
          end do
          if (k > kk) then
            delp(2) = delp(2)+dps
          else
            delp(k) = delp(k)+dps
          end if
          kfpla(i,j,n) = k

          ! Put single column back into 3-d fields.
          do k = 1,kk
            kn = k+nn
            temp(i,j,kn) = ttem(k)
            salt_corr(i,j) = salt_corr(i,j) - min(0._r8, ssal(k))*delp(k)/grav
            saln(i,j,kn) = max(0._r8, ssal(k))
            sigma(i,j,kn) = dens(k)
            dp(i,j,kn) = delp(k)
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
                trc_corr(i,j,nt) = trc_corr(i,j,nt) &
                                 - min(0._r8, ttrc(nt,k))*delp(k)/grav
                trc(i,j,kn,nt) = max(0._r8, ttrc(nt,k))
              end do
            end if
          end do

        end do
      end do
    end do
    !$omp end parallel do

    ! store 'old' interface pressures in -pu,pv-

    !$omp parallel do private(k,kn,l,i)
    do j = 1,jj
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            pu(i,j,k+1) = pu(i,j,k)+dpu(i,j,kn)
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            pv(i,j,k+1) = pv(i,j,k)+dpv(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! store 'new' layer thicknesses in -dpu,dpv-

    call xctilr(dp(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)

    !$omp parallel do private(k,kn,l,i)
    do j = -2,jj+2
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,kn,l,i,q)
    do j = -1,jj+2
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q = min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)= &
                 .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                 +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          end do
        end do
        do l = 1,isv(j)
          do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q = min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)= &
                 .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                    +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,k,kn,uo,po,pn,ko,un,udpn)
    do j = 1,jj

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k = 1,kk
            kn = k+nn
            uo(k) = u(i,j,kn)
          end do
          po(1) = 0.
          pn(1) = 0.
          do k = 2,kk+1
            po(k) = pu(i,j,k)
            pn(k) = .5*(min(pu(i,j,kk+1),p(i  ,j,k)) &
                      +min(pu(i,j,kk+1),p(i-1,j,k)))
          end do

          ko = 1
          do kn = 1,kk
            if (pn(kn+1)-pn(kn) == 0.) then
              un(kn) = 0.
            else
              udpn = 0.
              do while (pn(kn+1) > po(ko+1))
                udpn = udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko = ko+1
              end do
              un(kn) = (udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))&
                      /(pn(kn+1)-pn(kn))
            end if
          end do
          do k = 1,kk
            kn = k+nn
            u(i,j,kn) = un(k)
          end do

        end do
      end do

      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          do k = 1,kk
            kn = k+nn
            uo(k) = v(i,j,kn)
          end do
          po(1) = 0.
          pn(1) = 0.
          do k = 2,kk+1
            po(k) = pv(i,j,k)
            pn(k) = .5*(min(pv(i,j,kk+1),p(i,j  ,k)) &
                       +min(pv(i,j,kk+1),p(i,j-1,k)))
          end do

          ko = 1
          do kn = 1,kk
            if (pn(kn+1)-pn(kn) == 0.) then
              un(kn) = 0.
            else
              udpn = 0.
              do while (pn(kn+1) > po(ko+1))
                udpn = udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko = ko+1
              end do
              un(kn) = (udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn)))) &
                      /(pn(kn+1)-pn(kn))
            end if
          end do
          do k = 1,kk
            kn = k+nn
            v(i,j,kn) = un(k)
          end do

        end do
      end do
    end do
    !$omp end parallel do

    ! do j=1,jj
    !   do l=1,isu(j)
    !     do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
    !       q=0.
    !       do k=1,kk
    !         kn=k+nn
    !         q=q+u(i,j,kn)*dpu(i,j,kn)
    !       enddo
    !       if (abs(q).gt.1.e-4) then
    !         write (lp,*) 'mxlayr: u imbalance:',q,i,j
    !       endif
    !     enddo
    !   enddo
    !   do l=1,isv(j)
    !     do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
    !       q=0.
    !       do k=1,kk
    !         kn=k+nn
    !         q=q+v(i,j,kn)*dpv(i,j,kn)
    !       enddo
    !       if (abs(q).gt.1.e-4) then
    !         write (lp,*) 'mxlayr: v imbalance:',q,i,j
    !       endif
    !     enddo
    !   enddo
    ! enddo

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'mxlayr:'
      end if
      call chksummsk(dp,ip,2*kk,'dp')
      call chksummsk(temp,ip,2*kk,'temp')
      call chksummsk(saln,ip,2*kk,'saln')
      call chksummsk(sigma,ip,2*kk,'sigma')
      call chksummsk(u,iu,2*kk,'u')
      call chksummsk(v,iv,2*kk,'v')
      call chksummsk(dpu,iu,2*kk,'dpu')
      call chksummsk(dpv,iv,2*kk,'dpv')
      if (use_TRC) then
        do nt = 1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        end do
      end if
    end if

  end subroutine mxlayr

end module mod_mxlayr
