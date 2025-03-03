! ------------------------------------------------------------------------------
! Copyright (C) 2006-2024 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_diffus

  use mod_time,      only: delt1
  use mod_xc
  use mod_grid,      only: scuy, scvx, scp2, scuxi, scvyi
  use mod_eos,       only: sig
  use mod_state,     only: dp, temp, saln, sigma, &
                           utflx, vtflx, usflx, vsflx
  use mod_diffusion, only: ltedtp_opt, ltedtp_neutral, difiso, &
                           utflld, vtflld, usflld, vsflld
  use mod_checksum,  only: csdiag, chksummsk
  use mod_tracers,   only: ntr, itrtke, itrgls, trc, uflxtr, vflxtr
  use mod_ifdefs,    only: use_TRC, use_TKE, use_TKEIDF

  implicit none
  private

  public :: diffus

contains

  subroutine diffus(m,n,mm,nn,k1m,k1n)

    ! ------------------------------------------------------------------
    ! diffusion of tracers
    ! ------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    integer :: i,j,k,l,kn,km
    real :: q
    integer :: nt
    real :: dpeps
    parameter (dpeps = 1.e-5)

    call xctilr(dp(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
    if (ltedtp_opt == ltedtp_neutral) then
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_ps)
      if (use_TRC) then
        do nt = 1,ntr
          if (use_TKE .and. .not. use_TKEIDF) then
            if (nt == itrtke.or.nt == itrgls) cycle
          end if
          call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1,kk, 1,1, halo_ps)
        end do
      end if
      return
    else
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      if (use_TRC) then
        do nt = 1,ntr
          if (use_TKE .and. .not. use_TKEIDF) then
            if (nt == itrtke.or.nt == itrgls) cycle
          end if
          call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1,kk, 2,2, halo_ps)
        end do
      end if
      call xctilr(difiso, 1,kk, 2,2, halo_ps)
    end if

    do k = 1,kk
      kn = k+nn
      km = k+mm

      !$omp parallel do private(l,i,q,nt)
      do j = 0,jj+1
        do l = 1,isu(j)
          do i = max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            q = delt1*.5*(difiso(i-1,j,k)+difiso(i,j,k)) &
                *scuy(i,j)*scuxi(i,j) &
                *max(min(dp(i-1,j,kn),dp(i,j,kn)),dpeps)
            usflld(i,j,km) = q*(saln(i-1,j,kn)-saln(i,j,kn))
            utflld(i,j,km) = q*(temp(i-1,j,kn)-temp(i,j,kn))
            if (use_TRC) then
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEIDF) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                uflxtr(nt,i,j) = q*(trc(i-1,j,kn,nt)-trc(i,j,kn,nt))
              end do
            end if
            usflx(i,j,km) = usflx(i,j,km)+usflld(i,j,km)
            utflx(i,j,km) = utflx(i,j,km)+utflld(i,j,km)
          end do
        end do
      end do
      !$omp end parallel do

      !$omp parallel do private(l,i,q,nt)
      do j = 0,jj+2
        do l = 1,isv(j)
          do i = max(0,ifv(j,l)),min(ii+1,ilv(j,l))
            q = delt1*.5*(difiso(i,j-1,k)+difiso(i,j,k)) &
                *scvx(i,j)*scvyi(i,j) &
                *max(min(dp(i,j-1,kn),dp(i,j,kn)),dpeps)
            vsflld(i,j,km) = q*(saln(i,j-1,kn)-saln(i,j,kn))
            vtflld(i,j,km) = q*(temp(i,j-1,kn)-temp(i,j,kn))
            if (use_TRC) then
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEIDF) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                vflxtr(nt,i,j) = q*(trc(i,j-1,kn,nt)-trc(i,j,kn,nt))
              end do
            end if
            vsflx(i,j,km) = vsflx(i,j,km)+vsflld(i,j,km)
            vtflx(i,j,km) = vtflx(i,j,km)+vtflld(i,j,km)
          end do
        end do
      end do
      !$omp end parallel do

      !$omp parallel do private(l,i,q, nt)
      do j = 0,jj+1
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            q = 1./(scp2(i,j)*max(dp(i,j,kn),dpeps))
            saln(i,j,kn) = saln(i,j,kn) &
                           -q*(usflld(i+1,j,km)-usflld(i,j,km) &
                              +vsflld(i,j+1,km)-vsflld(i,j,km))
            temp(i,j,kn) = temp(i,j,kn) &
                          -q*(utflld(i+1,j,km)-utflld(i,j,km) &
                             +vtflld(i,j+1,km)-vtflld(i,j,km))
            if (use_TRC) then
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEIDF) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trc(i,j,kn,nt) = trc(i,j,kn,nt) &
                     -q*(uflxtr(nt,i+1,j)-uflxtr(nt,i,j) &
                        +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))
              end do
            end if
            sigma(i,j,kn) = sig(temp(i,j,kn),saln(i,j,kn))
          end do
        end do
      end do
      !$omp end parallel do

    end do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'diffus:'
      end if
      call chksummsk(temp,ip,2*kk,'temp')
      call chksummsk(saln,ip,2*kk,'saln')
      if (use_TRC) then
        do nt = 1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        end do
      end if
    end if

  end subroutine diffus

end module mod_diffus
