! ------------------------------------------------------------------------------
! Copyright (C) 2015-2024 Mats Bentsen, Mariana Vertenstein
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

module mod_idarlx

  use dimensions,  only: itdm, jtdm
  use mod_xc,      only: xchalt, xcaput, xctilr, mnproc, lp, nbdy, &
                         ii, jj, halo_ps
  use mod_forcing, only: tflxap, sflxap, tflxdi, sflxdi, nflxdi, &
                         aptflx, apsflx, ditflx, disflx

  implicit none
  private

  public :: idarlx

contains

  subroutine idarlx()

    ! ------------------------------------------------------------------
    ! Initialize diagnosing/application of relaxation fluxes
    ! ------------------------------------------------------------------

    ! Local varaibles
    real, dimension(itdm,jtdm) :: tmp2d
    integer :: i,j,k
    integer :: nfu

    if (aptflx) then
      if (mnproc == 1) then
        open (newunit=nfu,file='tflxdi.uf',form = 'unformatted')
        read (nfu) i,j
        if (i /= itdm.or.j /= jtdm) then
          write (lp,*) 'wrong dimensions in tflxdi.uf'
          call xchalt('(idarlx)')
          stop '(idarlx)'
        end if
      end if
      do k = 1,48
        if (mnproc == 1) then
          read (nfu) tmp2d
        end if
        call xcaput(tmp2d,tflxap(1-nbdy,1-nbdy,k),1)
      end do
      if (mnproc == 1) then
        close (unit = nfu)
      end if
      call xctilr(tflxap, 1,48, nbdy,nbdy, halo_ps)
    end if
    if (apsflx) then
      if (mnproc == 1) then
        open (newunit=nfu,file='sflxdi.uf',form = 'unformatted')
        read (nfu) i,j
        if (i /= itdm.or.j /= jtdm) then
          write (lp,*) 'wrong dimensions in sflxdi.uf'
          call xchalt('(idarlx)')
          stop '(idarlx)'
        end if
      end if
      do k = 1,48
        if (mnproc == 1) then
          read (nfu) tmp2d
        end if
        call xcaput(tmp2d,sflxap(1-nbdy,1-nbdy,k),1)
      end do
      if (mnproc == 1) then
        close (unit = nfu)
      end if
      call xctilr(sflxap, 1,48, nbdy,nbdy, halo_ps)
    end if

    if (ditflx.or.disflx) then
      do k = 1,48
        nflxdi(k) = 0
      end do
      if (ditflx) then
        !$omp parallel do private(k,i)
        do j = 1-nbdy,jj+nbdy
          do k = 1,48
            do i = 1-nbdy,ii+nbdy
              tflxdi(i,j,k) = 0.
            end do
          end do
        end do
        !$omp end parallel do
      end if
      if (disflx) then
        !$omp parallel do private(k,i)
        do j = 1-nbdy,jj+nbdy
          do k = 1,48
            do i = 1-nbdy,ii+nbdy
              sflxdi(i,j,k) = 0.
            end do
          end do
        end do
        !$omp end parallel do
      end if
    end if

  end subroutine idarlx

end module mod_idarlx
