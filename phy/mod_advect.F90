! ------------------------------------------------------------------------------
! Copyright (C) 2007-2020 Mats Bentsen, Mehmet Ilicak

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

module mod_advect

  ! --- ------------------------------------------------------------------
  ! --- This module contains variables and procedures related to advection
  ! --- of layer pressure thickness and tracers by calling incremental
  ! --- remapping routines.
  ! --- ------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm
  use mod_types,     only: r8
  use mod_constants, only: onemm
  use mod_time,      only: delt1, dlt
  use mod_xc,        only: xctilr, xcstop, ii, jj, kk, isp, ifp, ilp, &
                           iu, iv, ip, isu, ifu, ilu, isv, ifv, ilv, &
                           ifp, lp, halo_ps, mnproc, nbdy
  use mod_grid,      only: scuy, scvx, scp2i, scp2
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                           uflx, vflx, utflx, vtflx, usflx, vsflx, &
                           p, pbu, pbv, ubflxs_p, vbflxs_p
  use mod_diffusion, only: umfltd, vmfltd, umflsm, vmflsm
  use mod_remap,     only: remap_eitvel, remap_eitflx
  use mod_utility,   only: utotm, vtotm, umax, vmax
  use mod_checksum,  only: csdiag, chksummsk
  use mod_tracers,   only: ntr, itrtke, itrgls, trc, uflxtr, vflxtr
  use mod_ifdefs,    only: use_TRC, use_TKE, use_TKEADV

  implicit none
  private

  ! Variables to be set in namelist:
  ! Method of applying eddy-induced transport in the remap transport algorithm.
  ! Valid methods: 'eitvel', 'eitflx'.
  character(len = 80), public ::  rmpmth

  ! Public routines
  public :: advect

contains

  subroutine advect(m,n,mm,nn,k1m,k1n)

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    integer :: i,j,k,l,km,kn,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin,umflei,vmflei
    integer :: nt

    !$omp parallel do private( &
    !$omp l,i,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne)
    do j = -1,jj+2
      do l = 1,isp(j)
        do i = max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
          iw = i-iu(i  ,j)
          ie = i+iu(i+1,j)
          js = j-iv(i,j  )
          jn = j+iv(i,j+1)
          isw = i*(1-ip(iw,js))+iw*ip(iw,js)
          jsw = j*(1-ip(iw,js))+js*ip(iw,js)
          ise = i*(1-ip(ie,js))+ie*ip(ie,js)
          jse = j*(1-ip(ie,js))+js*ip(ie,js)
          inw = i*(1-ip(iw,jn))+iw*ip(iw,jn)
          jnw = j*(1-ip(iw,jn))+jn*ip(iw,jn)
          ine = i*(1-ip(ie,jn))+ie*ip(ie,jn)
          jne = j*(1-ip(ie,jn))+jn*ip(ie,jn)
          pbmin(i,j)= &
               min(p(isw,jsw,kk+1),p(i  ,js ,kk+1),p(ise,jse,kk+1), &
               p(iw ,j  ,kk+1),p(i  ,j  ,kk+1),p(ie ,j  ,kk+1), &
               p(inw,jnw,kk+1),p(i  ,jn ,kk+1),p(ine,jne,kk+1))
        end do
      end do
    end do
    !$omp end parallel do

    if (use_TRC) then
      do nt = 1,ntr
        if (use_TKE .and. .not. use_TKEADV) then
          if (nt == itrtke.or.nt == itrgls) cycle
        end if
        call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1,kk, 3,3, halo_ps)
      end do
    end if

    if (rmpmth == 'eitvel') then
      !$omp parallel do private(km,kn,j,l,i) &
      !$omp firstprivate(utotm,vtotm)
      do k = 1,kk
        km = k+mm
        kn = k+nn

        ! --- --- advective and diffusive velocity at mid time level

        do j = -1,jj+2
          do l = 1,isu(j)
            do i = max(0,ifu(j,l)),min(ii+2,ilu(j,l))
              utotm(i,j) = u(i,j,km) &
                   +(ubflxs_p(i,j,m)*dlt/pbu(i,j,m) &
                   +(umfltd(i,j,km)+umflsm(i,j,km)) &
                   /max(onemm,dpu(i,j,kn))) &
                   /(delt1*scuy(i,j))
              utotm(i,j) = max(-umax(i,j),min(umax(i,j),utotm(i,j)))
            end do
          end do
        end do
        do j = 0,jj+2
          do l = 1,isv(j)
            do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              vtotm(i,j) = v(i,j,km) &
                   +(vbflxs_p(i,j,m)*dlt/pbv(i,j,m) &
                   +(vmfltd(i,j,km)+vmflsm(i,j,km)) &
                   /max(onemm,dpv(i,j,kn))) &
                   /(delt1*scvx(i,j))
              vtotm(i,j) = max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
            end do
          end do
        end do

        call remap_eitvel(scuy,scvx,scp2i,scp2,pbmin, &
             pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n), &
             p(1-nbdy,1-nbdy,k+1),utotm,vtotm,delt1,1, &
             dp(1-nbdy,1-nbdy,kn), &
             temp(1-nbdy,1-nbdy,kn), &
             saln(1-nbdy,1-nbdy,kn), &
             uflx(1-nbdy,1-nbdy,km), &
             vflx(1-nbdy,1-nbdy,km), &
             utflx(1-nbdy,1-nbdy,km), &
             vtflx(1-nbdy,1-nbdy,km), &
             usflx(1-nbdy,1-nbdy,km), &
             vsflx(1-nbdy,1-nbdy,km), &
             kn,trc)
      end do
      !$omp end parallel do

    else if (rmpmth == 'eitflx') then

      !$omp parallel do private(km,kn,j,l,i) &
      !$omp firstprivate(utotm,vtotm)
      do k = 1,kk
        km = k+mm
        kn = k+nn

        ! --- --- advective velocity and total eddy-induced mass flux at mid
        ! --- --- time level
        do j = -1,jj+2
          do l = 1,isu(j)
            do i = max(0,ifu(j,l)),min(ii+2,ilu(j,l))
              utotm(i,j) = u(i,j,km) &
                   +dlt*ubflxs_p(i,j,m) &
                   /(delt1*pbu(i,j,m)*scuy(i,j))
              utotm(i,j) = max(-umax(i,j),min(umax(i,j),utotm(i,j)))
              umflei(i,j) = umfltd(i,j,km)+umflsm(i,j,km)
            end do
          end do
        end do
        do j = 0,jj+2
          do l = 1,isv(j)
            do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              vtotm(i,j) = v(i,j,km) &
                   +dlt*vbflxs_p(i,j,m) &
                   /(delt1*pbv(i,j,m)*scvx(i,j))
              vtotm(i,j) = max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
              vmflei(i,j) = vmfltd(i,j,km)+vmflsm(i,j,km)
            end do
          end do
        end do

        call remap_eitflx(scuy,scvx,scp2i,scp2,pbmin, &
             pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n), &
             p(1-nbdy,1-nbdy,k+1),utotm,vtotm, &
             umflei,vmflei, &
             delt1,1, &
             dp(1-nbdy,1-nbdy,kn), &
             temp(1-nbdy,1-nbdy,kn), &
             saln(1-nbdy,1-nbdy,kn), &
             uflx(1-nbdy,1-nbdy,km), &
             vflx(1-nbdy,1-nbdy,km), &
             utflx(1-nbdy,1-nbdy,km), &
             vtflx(1-nbdy,1-nbdy,km), &
             usflx(1-nbdy,1-nbdy,km), &
             vsflx(1-nbdy,1-nbdy,km), &
             kn,trc)
      end do
      !$omp end parallel do
    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' rmpmth = ',trim(rmpmth),' is unsupported!'
      end if
      call xcstop('(advect)')
      stop '(advect)'
    end if

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'advect:'
      end if
      call chksummsk(dp,ip,2*kk,'dp')
      call chksummsk(temp,ip,2*kk,'temp')
      call chksummsk(saln,ip,2*kk,'saln')
      call chksummsk(uflx,iu,2*kk,'uflx')
      call chksummsk(vflx,iv,2*kk,'vflx')
      if (use_TRC) then
        do nt = 1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        end do
      end if
    end if

  end subroutine advect

end module mod_advect
