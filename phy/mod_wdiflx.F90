! ------------------------------------------------------------------------------
! Copyright (C) 2005-2024 Mats Bentsen, Alok Kumar Gupta, Mariana Vertenstein
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

module mod_wdiflx

  use mod_config,  only: runid
  use mod_time,    only: date,time
  use mod_xc
  use mod_forcing, only: tflxdi, sflxdi, nflxdi, ditflx, disflx
  use mod_nctools
  use mod_dia,     only : iotype

  implicit none
  private

  public :: wdiflx

contains

  subroutine wdiflx()

    ! --- Write accumulated diagnosed heat and salt fluxes

    ! Local variables
    character(len=256) :: fname
    integer :: i,j,k

    if (ditflx) then

      write (fname,'(2a,i4.4,a)') &
           trim(runid),'_tflxdi_',date%year-1,'.nc'
      if (mnproc == 1) then
        write (lp,'(2a)') &
             'Writing diagnostic heat flux to ',trim(fname)
      end if
      call ncfopn(fname,'w','c',1,iotype)

#ifdef COMPDIA
      call ncdimc('pcomp',ip,0)
#else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#endif
      call ncdims('week',48)

      call ncputr('time',time)
      call ncputi('nflxdi',nflxdi)

#ifdef COMPDIA
      call ncdefvar('tflxdi','pcomp week',ndouble,8)
      call ncedef
      call nccomp('tflxdi','pcomp week',tflxdi,ip,1.,0.,8)
#else
      call ncdefvar('tflxdi','x y week',ndouble,8)
      call ncedef
      call ncwrtr('tflxdi','x y week',tflxdi,ip,1,1.,0.,8)
#endif

      call ncfcls

      !$omp parallel do private(j,i)
      do k = 1,48
        nflxdi(k) = 0
        do j = 1,jj
          do i = 1,ii
            tflxdi(i,j,k) = 0.
          end do
        end do
      end do
      !$omp end parallel do

    end if

    if (disflx) then

      write (fname,'(2a,i4.4,a)') &
           trim(runid),'_sflxdi_',date%year-1,'.nc'
      if (mnproc == 1) then
        write (lp,'(2a)') &
             'Writing diagnostic salt flux to ',trim(fname)
      end if
      call ncfopn(fname,'w','c',1,iotype)

#ifdef COMPDIA
      call ncdimc('pcomp',ip,0)
#else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#endif
      call ncdims('week',48)

      call ncputr('time',time)
      call ncputi('nflxdi',nflxdi)

#ifdef COMPDIA
      call nccomp('sflxdi','pcomp week',sflxdi,ip,1.,0.,8)
#else
      call ncwrtr('sflxdi','x y week',sflxdi,ip,1,1.,0.,8)
#endif

      call ncfcls

      !$omp parallel do private(j,i)
      do k = 1,48
        nflxdi(k) = 0
        do j = 1,jj
          do i = 1,ii
            sflxdi(i,j,k) = 0.
          end do
        end do
      end do
      !$omp end parallel do

    end if

  end subroutine wdiflx

end module mod_wdiflx
