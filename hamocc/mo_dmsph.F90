! Copyright (C) 2020  J. Tjiputra
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.

module mo_dmsph

  implicit none
  private
  public :: get_pi_ph,get_dmsph,pi_ph_path,with_dmsph

  ! Activate/deactivate calculation of DMS as a function of pH.
  logical :: with_dmsph = .false.

  ! Path to input data, set through namelist in hamocc_init.F
  character(len=256),save    :: pi_ph_path = ''

  ! Length of surface PI pH record from file
  ! Assume monthly record by default (otherwise change get_dmsph)
  integer, parameter :: pi_ph_record = 12

  ! dms_ph scaling factor
  real, parameter :: dms_gamma = 0.87

  ! surface PI pH field
  real, dimension(:,:,:), allocatable :: pi_ph

CONTAINS


  SUBROUTINE GET_PI_PH(kpie,kpje,omask)
  !**********************************************************************

    use mo_control_bgc, only: io_stdo_bgc
    use netcdf,         only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mod_xc,         only: mnproc,xchalt

    implicit none
    INTEGER, INTENT(in) :: kpie,kpje
    INTEGER ::i,j,l

    REAL,intent(in) ::omask(kpie,kpje)

    ! define the fields
    REAL :: pi_ph_in(kpie,kpje,pi_ph_record)
    INTEGER ncid,ncstat

    ! allocate pi_ph field
    if(.not. allocated(pi_ph)) call alloc_pi_ph(kpie,kpje)

    !
    ! Open netCDF data file
    !
    IF(mnproc==1) THEN
       ncstat = NF90_OPEN(trim(pi_ph_path)//'MONTHLY_PI_PH.nc',        &
     &                   NF90_NOWRITE, ncid)
       write(io_stdo_bgc,*) 'HAMOCC: opening MONTHLY_PI_PH file'
       IF (ncstat.NE.NF90_NOERR ) THEN
          CALL xchalt('(get_pi_ph: Problem with netCDF1)')
                 stop '(get_pi_ph: Problem with netCDF1)'
       END IF
    END IF
    !
    ! Read  data
    call read_netcdf_var(ncid,'pH',pi_ph_in(1,1,1),pi_ph_record,0,0)
    !
    ! Close file
    IF(mnproc==1) THEN
       ncstat = NF90_CLOSE(ncid)
       IF ( ncstat .NE. NF90_NOERR ) THEN
          CALL xchalt('(get_pi_ph: Problem with netCDF200)')
                 stop '(get_pi_ph: Problem with netCDF200)'
       END IF
    END IF

    ! set missings over land
    do l=1,pi_ph_record
       do j=1,kpje
          do i=1,kpie
             if(omask(i,j).gt.0.5) then
                pi_ph(i,j,l) = pi_ph_in(i,j,l)
             else
                pi_ph(i,j,l) = 0.
             endif
          enddo
       enddo
    enddo

  END SUBROUTINE GET_PI_PH


  !**********************************************************************
  ! Allocate the PI_PH field.
  !**********************************************************************
  subroutine alloc_pi_ph(kpie,kpje)
    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc

    implicit none
    integer, intent(in) :: kpie,kpje
    integer             :: errstat

    IF (mnproc.eq.1) THEN
       WRITE(io_stdo_bgc,*)'Memory allocation for variable pi_ph ...'
       WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
       WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
       WRITE(io_stdo_bgc,*)'Third dimension    : ',pi_ph_record
    ENDIF

    ALLOCATE (pi_ph(kpie,kpje,pi_ph_record),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pi_ph'
    pi_ph(:,:,:) = 0.0

  end subroutine alloc_pi_ph


  !**********************************************************************
  ! Get DMS_PH as a function of pH
  !**********************************************************************
  function get_dmsph(i,j) result(dms_ph)
    use mo_carbch, only: hi
    use mod_time,  only: date

    implicit none

    integer, intent(in) :: i,j
    real :: dms_ph

    dms_ph  = 1. + (-log10(hi(i,j,1)) - pi_ph(i,j,date%month))*dms_gamma

  end function get_dmsph


end module mo_dmsph
