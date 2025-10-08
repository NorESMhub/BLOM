! ------------------------------------------------------------------------------
! Copyright (C) 2005-2025 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_pbcor

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to the
  ! application of corrected layer-wise mass fluxes so that the
  ! vertical sum of layer thicknesses match the bottom pressure found
  ! by the integration of the barotropic equations.
  ! ------------------------------------------------------------------

  use dimensions,    only: idm, jdm
  use mod_types,     only: r8
  use mod_constants, only: epsilp
  use mod_time,      only: dlt
  use mod_xc,        only: xctilr, xcstop, lp, ip, iu, iv, ii, jj, kk, &
                           isp, ifp, lp, isu, ifu, ilu, isv, ifv, ilv, &
                           ilp, nbdy, mnproc, halo_ps, halo_uv, halo_vv
  use mod_grid,      only: scp2i
  use mod_eos,       only: sig
  use mod_state,     only: dp, temp, saln, sigma, uflx, vflx, &
                           utflx, vtflx, usflx, vsflx, p, &
                           ubflxs, vbflxs, pb, ubflxs_p, vbflxs_p, pb_p
  use mod_utility,   only: utotm, vtotm, utotn, vtotn, uflux, vflux, &
                           uflux2, vflux2, uflux3, vflux3
  use mod_checksum,  only: csdiag, chksum
  use mod_tracers,   only: ntr, trc, uflxtr, vflxtr, itrtke, itrgls
  use mod_ifdefs,    only: use_TRC, use_TKE, use_TKEADV

  implicit none
  private

  ! Variables to be set in namelist:
  character(len = 80), public :: &
       bmcmth ! Baroclinic mass flux correction method. Valid methods:
              ! 'uc' (upstream column), 'dluc' (depth limited upstream
              ! column).

  ! Parameters:
  real(r8), parameter :: &
       dpeps1 = 1.e-5_r8, & ! Small layer pressure thickness [kg m-1 s-2].
       dpeps2 = 1.e-7_r8    ! Small layer pressure thickness [kg m-1 s-2].

  ! Public routines
  public :: pbcor1, pbcor2

contains

  subroutine pbcor1(m,n,mm,nn,k1m,k1n)

    ! ------------------------------------------------------------------
    ! Correct the layer thicknesses to match the predictive bottom
    ! pressure from the barotropic solution.
    ! ------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu_t,pbv_t
    real, dimension(1-nbdy:idm+nbdy) :: pbfac
    real :: dpo,dpni
    integer :: i,j,k,l,kn,km
    integer :: nt
    character(len = 2) cnt

    !$omp parallel do private(k,kn,l,i)
    do j = 0,jj+1
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,k,km)
    do j = 1,jj
      if     (bmcmth == 'uc') then
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j) = dlt*ubflxs_p(i,j,m)
          end do
        end do
      else if (bmcmth == 'dluc') then
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j) = dlt*ubflxs_p(i,j,m)
            pbu_t(i,j) = min(p(i,j,kk+1),p(i-1,j,kk+1))
          end do
        end do
      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' bmcmth = ',trim(bmcmth),' is unsupported!'
        end if
        call xcstop('(pbcor1)')
        stop '(pbcor1)'
      end if
      do k = 1,kk
        km = k+mm
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j) = utotm(i,j)-uflx(i,j,km)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,k,km)
    do j = 1,jj+1
      if     (bmcmth == 'uc') then
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j) = dlt*vbflxs_p(i,j,m)
          end do
        end do
      else if (bmcmth == 'dluc') then
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j) = dlt*vbflxs_p(i,j,m)
            pbv_t(i,j) = min(p(i,j,kk+1),p(i,j-1,kk+1))
          end do
        end do
      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' bmcmth = ',trim(bmcmth),' is unsupported!'
        end if
        call xcstop('(pbcor1)')
        stop '(pbcor1)'
      end if
      do k = 1,kk
        km = k+mm
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j) = vtotm(i,j)-vflx(i,j,km)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    do k = 1,kk
      km = k+mm
      kn = k+nn

      if     (bmcmth == 'uc') then
        !$omp parallel do private(l,i,nt)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotm(i,j) > 0.) then
                  uflux(i,j) = utotm(i,j)*dp(i-1,j,kn)/p(i-1,j,kk+1)
                 uflux2(i,j) = uflux(i,j)*saln(i-1,j,kn)
                 uflux3(i,j) = uflux(i,j)*temp(i-1,j,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i-1,j,kn,nt)
                  end do
                end if
              else
                 uflux(i,j) = utotm(i,j)*dp(i  ,j,kn)/p(i  ,j,kk+1)
                uflux2(i,j) = uflux(i,j)*saln(i  ,j,kn)
                uflux3(i,j) = uflux(i,j)*temp(i  ,j,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i  ,j,kn,nt)
                  end do
                end if
              end if
               uflx(i,j,km) = uflx(i,j,km)+uflux(i,j)
              usflx(i,j,km) = usflx(i,j,km)+uflux2(i,j)
              utflx(i,j,km) = utflx(i,j,km)+uflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do

        !$omp parallel do private(l,i,nt)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotm(i,j) > 0.) then
                 vflux(i,j) = vtotm(i,j)*dp(i,j-1,kn)/p(i,j-1,kk+1)
                vflux2(i,j) = vflux(i,j)*saln(i,j-1,kn)
                vflux3(i,j) = vflux(i,j)*temp(i,j-1,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j-1,kn,nt)
                  end do
                end if
              else
                 vflux(i,j) = vtotm(i,j)*dp(i,j  ,kn)/p(i,j  ,kk+1)
                vflux2(i,j) = vflux(i,j)*saln(i,j  ,kn)
                vflux3(i,j) = vflux(i,j)*temp(i,j  ,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j  ,kn,nt)
                  end do
                end if
              end if
               vflx(i,j,km) = vflx(i,j,km)+vflux(i,j)
              vsflx(i,j,km) = vsflx(i,j,km)+vflux2(i,j)
              vtflx(i,j,km) = vtflx(i,j,km)+vflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (bmcmth == 'dluc') then

        !$omp parallel do private(l,i,nt)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotm(i,j) > 0.) then
                uflux(i,j)= &
                     utotm(i,j) &
                     *max(0.,min(pbu_t(i,j),&
                                 p(i-1,j,k+1))-p(i-1,j,k))/pbu_t(i,j)
                uflux2(i,j) = uflux(i,j)*saln(i-1,j,kn)
                uflux3(i,j) = uflux(i,j)*temp(i-1,j,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i-1,j,kn,nt)
                  end do
                end if
              else
                uflux(i,j)= &
                     utotm(i,j) &
                     *max(0.,min(pbu_t(i,j),&
                                 p(i,j,k+1))-p(i,j,k))/pbu_t(i,j)
                uflux2(i,j) = uflux(i,j)*saln(i  ,j,kn)
                uflux3(i,j) = uflux(i,j)*temp(i  ,j,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i  ,j,kn,nt)
                  end do
                end if
              end if
              uflx(i,j,km) = uflx(i,j,km)+uflux(i,j)
              usflx(i,j,km) = usflx(i,j,km)+uflux2(i,j)
              utflx(i,j,km) = utflx(i,j,km)+uflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do

        !$omp parallel do private(l,i,nt)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotm(i,j) > 0.) then
                vflux(i,j)= &
                     vtotm(i,j) &
                     *max(0.,min(pbv_t(i,j),&
                                 p(i,j-1,k+1))-p(i,j-1,k))/pbv_t(i,j)
                vflux2(i,j) = vflux(i,j)*saln(i,j-1,kn)
                vflux3(i,j) = vflux(i,j)*temp(i,j-1,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j-1,kn,nt)
                  end do
                end if
              else
                vflux(i,j)= &
                     vtotm(i,j) &
                     *max(0.,min(pbv_t(i,j),&
                                 p(i,j,k+1))-p(i,j,k))/pbv_t(i,j)
                vflux2(i,j) = vflux(i,j)*saln(i,j  ,kn)
                vflux3(i,j) = vflux(i,j)*temp(i,j  ,kn)
                if (use_TRC) then
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j  ,kn,nt)
                  end do
                end if
              end if
               vflx(i,j,km) = vflx(i,j,km)+vflux(i,j)
              vsflx(i,j,km) = vsflx(i,j,km)+vflux2(i,j)
              vtflx(i,j,km) = vtflx(i,j,km)+vflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' bmcmth = ',trim(bmcmth),' is unsupported!'
        end if
        call xcstop('(pbcor1)')
        stop '(pbcor1)'
      end if

      !$omp parallel do private(l,i,dpo,dpni,nt)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpo = dp(i,j,kn)
            dp(i,j,kn) = max(0.,dpo-(uflux(i+1,j)-uflux(i,j) &
                                    +vflux(i,j+1)-vflux(i,j))*scp2i(i,j))
            dpo = dpo+dpeps1
            dpni = 1./(dp(i,j,kn)+dpeps1)
            saln(i,j,kn) = (dpo*saln(i,j,kn) &
                           -(uflux2(i+1,j)-uflux2(i,j) &
                            +vflux2(i,j+1)-vflux2(i,j))*scp2i(i,j))*dpni
            temp(i,j,kn) = (dpo*temp(i,j,kn) &
                           -(uflux3(i+1,j)-uflux3(i,j) &
                            +vflux3(i,j+1)-vflux3(i,j))*scp2i(i,j))*dpni
            if (use_TRC) then
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trc(i,j,kn,nt) = (dpo*trc(i,j,kn,nt) &
                                 -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j) &
                                  +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
              end do
            end if
            if (dp(i,j,kn) < dpeps2) dp(i,j,kn) = 0.
          end do
        end do
      end do
      !$omp end parallel do

    end do

    !$omp parallel do private(k,kn,l,i,pbfac)
    do j = 1,jj
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfac(i) = pb_p(i,j)/p(i,j,kk+1)
        end do
      end do
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,kn) = dp(i,j,kn)*pbfac(i)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'pbcor1:'
      end if
      call chksum(dp  , 2*kk, halo_ps, 'dp'  )
      call chksum(temp, 2*kk, halo_ps, 'temp')
      call chksum(saln, 2*kk, halo_ps, 'saln')
      call chksum(uflx, 2*kk, halo_uv, 'uflx')
      call chksum(vflx, 2*kk, halo_vv, 'vflx')
      if (use_TRC) then
        do nt = 1,ntr
          write(cnt, '(i2.2)') nt
          call chksum(trc(1-nbdy,1-nbdy,1,nt), 2*kk, halo_ps, 'trc'//cnt)
        end do
      end if
    end if

  end subroutine pbcor1

  ! ------------------------------------------------------------------

  subroutine pbcor2(m,n,mm,nn,k1m,k1n)

    ! ------------------------------------------------------------------
    ! Correct the layer thicknesses to better match the corrected bottom
    ! pressure from the barotropic solution.
    ! ------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu_t,pbv_t
    real, dimension(1-nbdy:idm+nbdy) :: pbfac
    real :: dpo,dpni
    integer :: i,j,k,l,kn,km
    integer :: nt
    character(len = 2) cnt

    call xctilr(ubflxs(1-nbdy,1-nbdy,n), 1,1, 1,1, halo_uv)
    call xctilr(vbflxs(1-nbdy,1-nbdy,n), 1,1, 1,1, halo_vv)
    if (use_TRC) then
      do nt = 1,ntr
        call xctilr(trc(1-nbdy,1-nbdy,k1m,nt), 1,kk, 1,1, halo_ps)
      end do
    end if

    !$omp parallel do private(k,km,l,i)
    do j = 0,jj+1
      do k = 1,kk
        km = k+mm
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,km) = max(0.,dp(i,j,km))+epsilp
            p(i,j,k+1) = p(i,j,k)+dp(i,j,km)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,k,kn)
    do j = 1,jj
      if     (bmcmth == 'uc') then
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j) = dlt*ubflxs(i,j,n)
          end do
        end do
      else if (bmcmth == 'dluc') then
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j) = dlt*ubflxs(i,j,n)
            pbu_t(i,j) = min(p(i,j,kk+1),p(i-1,j,kk+1))
          end do
        end do
      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' bmcmth = ',trim(bmcmth),' is unsupported!'
        end if
        call xcstop('(pbcor2)')
        stop '(pbcor2)'
      end if
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j) = utotn(i,j)-uflx(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i,k,kn)
    do j = 1,jj+1
      if     (bmcmth == 'uc') then
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j) = dlt*vbflxs(i,j,n)
          end do
        end do
      else if (bmcmth == 'dluc') then
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j) = dlt*vbflxs(i,j,n)
            pbv_t(i,j) = min(p(i,j,kk+1),p(i,j-1,kk+1))
          end do
        end do
      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' bmcmth = ',trim(bmcmth),' is unsupported!'
        end if
        call xcstop('(pbcor2)')
        stop '(pbcor2)'
      end if
      do k = 1,kk
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j) = vtotn(i,j)-vflx(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    do k = 1,kk
      kn = k+nn
      km = k+mm

      if     (bmcmth == 'uc') then
        !$omp parallel do private(l,i,nt)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotn(i,j) > 0.) then
                 uflux(i,j) = utotn(i,j)*dp(i-1,j,km)/p(i-1,j,kk+1)
                uflux2(i,j) = uflux(i,j)*saln(i-1,j,km)
                uflux3(i,j) = uflux(i,j)*temp(i-1,j,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i-1,j,km,nt)
                  end do
                end if
              else
                 uflux(i,j) = utotn(i,j)*dp(i  ,j,km)/p(i  ,j,kk+1)
                uflux2(i,j) = uflux(i,j)*saln(i  ,j,km)
                uflux3(i,j) = uflux(i,j)*temp(i  ,j,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i  ,j,km,nt)
                  end do
                end if
              end if
               uflx(i,j,kn) = uflx(i,j,kn)+uflux(i,j)
              usflx(i,j,kn) = usflx(i,j,kn)+uflux2(i,j)
              utflx(i,j,kn) = utflx(i,j,kn)+uflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do

        !$omp parallel do private(l,i,nt)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotn(i,j) > 0.) then
                 vflux(i,j) = vtotn(i,j)*dp(i,j-1,km)/p(i,j-1,kk+1)
                vflux2(i,j) = vflux(i,j)*saln(i,j-1,km)
                vflux3(i,j) = vflux(i,j)*temp(i,j-1,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j-1,km,nt)
                  end do
                end if
              else
                 vflux(i,j) = vtotn(i,j)*dp(i,j  ,km)/p(i,j  ,kk+1)
                vflux2(i,j) = vflux(i,j)*saln(i,j  ,km)
                vflux3(i,j) = vflux(i,j)*temp(i,j  ,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j  ,km,nt)
                  end do
                end if
              end if
               vflx(i,j,kn) = vflx(i,j,kn)+vflux(i,j)
              vsflx(i,j,kn) = vsflx(i,j,kn)+vflux2(i,j)
              vtflx(i,j,kn) = vtflx(i,j,kn)+vflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (bmcmth == 'dluc') then
        !$omp parallel do private(l,i,nt)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotn(i,j) > 0.) then
                uflux(i,j)= &
                     utotn(i,j) &
                     *max(0.,min(pbu_t(i,j),&
                                 p(i-1,j,k+1))-p(i-1,j,k))/pbu_t(i,j)
                uflux2(i,j) = uflux(i,j)*saln(i-1,j,km)
                uflux3(i,j) = uflux(i,j)*temp(i-1,j,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i-1,j,km,nt)
                  end do
                end if
              else
                uflux(i,j)= &
                     utotn(i,j) &
                     *max(0.,min(pbu_t(i,j),&
                                 p(i,j,k+1))-p(i,j,k))/pbu_t(i,j)
                uflux2(i,j) = uflux(i,j)*saln(i  ,j,km)
                uflux3(i,j) = uflux(i,j)*temp(i  ,j,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    uflxtr(nt,i,j) = uflux(i,j)*trc(i  ,j,km,nt)
                  end do
                end if
              end if
               uflx(i,j,kn) = uflx(i,j,kn)+uflux(i,j)
              usflx(i,j,kn) = usflx(i,j,kn)+uflux2(i,j)
              utflx(i,j,kn) = utflx(i,j,kn)+uflux3(i,j)
            end do
          end do
        end do
        !$omp end parallel do

        !$omp parallel do private(l,i,nt)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotn(i,j) > 0.) then
                vflux(i,j)= &
                     vtotn(i,j) &
                     *max(0.,min(pbv_t(i,j),&
                                 p(i,j-1,k+1))-p(i,j-1,k))/pbv_t(i,j)
                vflux2(i,j) = vflux(i,j)*saln(i,j-1,km)
                vflux3(i,j) = vflux(i,j)*temp(i,j-1,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j-1,km,nt)
                  end do
                end if
              else
                vflux(i,j)= &
                     vtotn(i,j) &
                     *max(0.,min(pbv_t(i,j),&
                                 p(i,j,k+1))-p(i,j,k))/pbv_t(i,j)
                vflux2(i,j) = vflux(i,j)*saln(i,j  ,km)
                vflux3(i,j) = vflux(i,j)*temp(i,j  ,km)
                if (use_TRC) then
                  do nt = 1,ntr
                    vflxtr(nt,i,j) = vflux(i,j)*trc(i,j  ,km,nt)
                  end do
                end if
              end if
               vflx(i,j,kn) = vflx(i,j,kn)+vflux(i,j)
              vsflx(i,j,kn) = vsflx(i,j,kn)+vflux2(i,j)
              vtflx(i,j,kn) = vtflx(i,j,kn)+vflux3(i,j)
            end do
          end do
        end do
      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' bmcmth = ',trim(bmcmth),' is unsupported!'
        end if
        call xcstop('(pbcor2)')
        stop '(pbcor2)'
      end if

      !$omp parallel do private(l,i,dpo,dpni,nt)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpo = dp(i,j,km)
            dp(i,j,km) = dpo-scp2i(i,j)*(uflux(i+1,j)-uflux(i,j) &
                                        +vflux(i,j+1)-vflux(i,j))
            dpni = 1./dp(i,j,km)
            saln(i,j,km) = (dpo*saln(i,j,km) &
                 -scp2i(i,j)*(uflux2(i+1,j)-uflux2(i,j) &
                             +vflux2(i,j+1)-vflux2(i,j)))*dpni
            temp(i,j,km) = (dpo*temp(i,j,km) &
                 -scp2i(i,j)*(uflux3(i+1,j)-uflux3(i,j) &
                             +vflux3(i,j+1)-vflux3(i,j)))*dpni
            if (use_TRC) then
              do nt = 1,ntr
                trc(i,j,km,nt) = (dpo*trc(i,j,km,nt) &
                     -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j) &
                      +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
              end do
            end if
            sigma(i,j,km) = sig(temp(i,j,km),saln(i,j,km))
            dp(i,j,km) = dp(i,j,km)-epsilp
            if (dp(i,j,km) < dpeps2) dp(i,j,km) = 0.
          end do
        end do
      end do
      !$omp end parallel do

    end do

    !$omp parallel do private(k,km,l,i,pbfac)
    do j = 1,jj
      do k = 1,kk
        km = k+mm
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,km)
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfac(i) = pb(i,j,m)/p(i,j,kk+1)
        end do
      end do
      do k = 1,kk
        km = k+mm
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,km) = dp(i,j,km)*pbfac(i)
            p(i,j,k+1) = p(i,j,k)+dp(i,j,km)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'pbcor2:'
      end if
      call chksum(dp   , 2*kk, halo_ps, 'dp'   )
      call chksum(temp , 2*kk, halo_ps, 'temp' )
      call chksum(saln , 2*kk, halo_ps, 'saln' )
      call chksum(p    , kk+1, halo_ps, 'p'    )
      call chksum(sigma, 2*kk, halo_ps, 'sigma')
      call chksum(uflx , 2*kk, halo_uv, 'uflx' )
      call chksum(vflx , 2*kk, halo_vv, 'vflx' )
      do nt = 1,ntr
        write(cnt, '(i2.2)') nt
        call chksum(trc(1-nbdy,1-nbdy,1,nt), 2*kk, halo_ps, 'trc'//cnt)
      end do
    end if

  end subroutine pbcor2

end module mod_pbcor
