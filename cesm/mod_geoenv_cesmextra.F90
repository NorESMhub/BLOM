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

module mod_geoenv_cesmextra

  use dimensions,  only: itdm, jtdm
  use mod_xc,      only: cplmsk, mnproc, lp, ii, jj, xchalt, xcbcst, xcaput
  use mod_grid,    only: grfile
  use mod_utility, only: util1
  use netcdf

  implicit none
  private

  public :: geoenv_cesmextra

contains

  subroutine geoenv_cesmextra()

    ! ------------------------------------------------------------------
    ! Read additional grid parameters when configured for coupling to
    ! CESM
    ! ------------------------------------------------------------------

    ! Local variables
    real, dimension(itdm,jtdm) :: tmpg
    integer :: i,j,status,ncid,dimid,varid

    ! ------------------------------------------------------------------
    ! read grid information from grid file
    ! ------------------------------------------------------------------

    if (mnproc == 1) then
      write (lp,'(2a)') ' reading additional grid information from ',trim(grfile)

      ! - open netcdf file
      status = nf90_open(grfile,nf90_nowrite,ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_open: ',trim(grfile),': ', nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if

      ! check dimensions
      status = nf90_inq_dimid(ncid,'x',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = i)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: x: ', nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
      status = nf90_inq_dimid(ncid,'y',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = j)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: y: ', nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
      if (i /= itdm.or.j /= jtdm) then
        write (lp,'(2a)') ' wrong dimensions in ',trim(grfile)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
    end if

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'cplmask',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: cplmask: ', nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
    end if
    if (mnproc == 1) then
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: cplmask: ', nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
    end if
    call xcaput(tmpg,util1,1)
    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        cplmsk(i,j) = nint(util1(i,j))
      end do
    end do
    !$omp end parallel do

    ! close grid information file

    if (mnproc == 1) then
      status = nf90_close(ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_close: ',trim(grfile),': ', nf90_strerror(status)
        call xchalt('(geoenv_cesmextra)')
        stop '(geoenv_cesmextra)'
      end if
    end if

  end subroutine geoenv_cesmextra

end module mod_geoenv_cesmextra
