! ------------------------------------------------------------------------------
! Copyright (C) 2006-2024 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin

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

module mod_temmin

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to defining
  ! the lower bound of physical temperature values in isopycnic
  ! layers.
  ! ------------------------------------------------------------------

  use dimensions,    only: idm,jdm, kdm
  use mod_types,     only: r8
  use mod_config,    only: expcnf
  use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml, sigmar
  use mod_xc,        only: xcstop, ii, jj, kk, isp, ifp, ilp, lp, &
                           mnproc, nbdy
  use mod_eos,       only: ap11, ap12, ap13, ap14, ap15, ap16, &
                           ap21, ap22, ap23, ap24, ap25, ap26, atf
  use mod_pointtest, only: itest, jtest, ptest

  implicit none
  private

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), public :: &
       temmin ! Minimum temperature allowed in an isopycnic layer
              ! [deg C].

  public :: settemmin

contains

  ! --- ------------------------------------------------------------------

  subroutine settemmin()

    ! ------------------------------------------------------------------
    ! Set minimum physical temperature values in isopycnic layers
    ! ------------------------------------------------------------------

    integer :: i,j,k,l
    real :: salfrz,a,b,c

    if     (vcoord_tag /= vcoord_isopyc_bulkml .or. &
         expcnf == 'cesm' .or. expcnf == 'single_column') then

      ! Set temmin to a constant freezing temperature for all layers
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              temmin(i,j,k) = -3.
            end do
          end do
        end do
      end do
      !$omp end parallel do

    else if (expcnf == 'ben02clim'.or.expcnf == 'ben02syn'.or. &
             expcnf == 'fuk95'.or.expcnf == 'channel') then

      ! Let temmin be the freezing temperature of a given potential
      ! density. This can be achieved by using potential density given
      ! in the function sig and the salinity dependent freezing
      ! temperature given in the function swtfrz.

      !$omp parallel do private(k,l,i,a,b,c,salfrz)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              a = ((ap14-ap24*sigmar(i,j,k))*atf &
                   + ap15-ap25*sigmar(i,j,k) )*atf &
                   +ap16-ap26*sigmar(i,j,k)
              b = (ap12-ap22*sigmar(i,j,k))*atf+ap13-ap23*sigmar(i,j,k)
              c = ap11-ap21*sigmar(i,j,k)
              salfrz = (-b+sqrt(b*b-4.*a*c))/(2.*a)
              temmin(i,j,k) = atf*salfrz
            end do
          end do
        end do
      end do
      !$omp end parallel do

    else if (expcnf == 'isomip1'.or.expcnf == 'isomip2') then

      ! Set temmin to a low value.
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              temmin(i,j,k) = -10.
            end do
          end do
        end do
      end do
      !$omp end parallel do

    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' expcnf = ',trim(expcnf),' is unsupported!'
      end if
      call xcstop('(settemmin)')
      stop '(settemmin)'
    end if

    if (mnproc == ptest) then
      write (lp,'(a/(6(i5,f8.3)))') 'minimum temperature values:', &
           (k,temmin(itest,jtest,k),k = 2,kk)
    end if

  end subroutine settemmin

end module mod_temmin
