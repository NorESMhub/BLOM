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

module mo_read_pi_ph

  implicit none
  private
  public :: ini_pi_ph,get_pi_ph,pi_ph_file,pi_ph

  ! Path to input data, set through namelist in hamocc_init.F
  character(len=256),save    :: pi_ph_file = ''

  ! Length of surface PI pH record from file
  !  - Current implementation only support monthly records.
  integer, parameter :: pi_ph_record = 12

  ! surface PI pH climatology
  real, dimension(:,:,:), allocatable :: pi_ph_clim

  ! surface PI pH monthly data
  real, dimension(:,:), allocatable :: pi_ph

CONTAINS


  !**********************************************************************
  ! PUBLIC SUBROUTINE : INI_PI_PH
  !
  ! Initialise the PI_PH field from climatology.
  !**********************************************************************
  subroutine ini_pi_ph(kpie,kpje,omask)
    use mo_control_bgc, only: io_stdo_bgc,with_dmsph
    use netcdf,         only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mod_xc,         only: mnproc,xchalt

    implicit none
    INTEGER, INTENT(in) :: kpie,kpje
    INTEGER ::i,j,l

    REAL,intent(in) ::omask(kpie,kpje)

    ! define the fields
    REAL :: pi_ph_in(kpie,kpje,pi_ph_record)
    INTEGER ncid,ncstat

    ! Allocate pi_ph field (required argument for hmaocc4bcm)
    if(.not. allocated(pi_ph)) call alloc_pi_ph(kpie,kpje)

    ! Only read PI pH climatology if needed for DMS
    if(with_dmsph) then

      ! allocate pi_ph_clim field
      if(.not. allocated(pi_ph_clim)) call alloc_pi_ph_clim(kpie,kpje)

      ! Open netCDF data file
      IF(mnproc==1) THEN
        ncstat = NF90_OPEN(trim(pi_ph_file), NF90_NOWRITE, ncid)
        write(io_stdo_bgc,*) 'HAMOCC: opening PI_PH climatology file'
        IF (ncstat.NE.NF90_NOERR ) THEN
          CALL xchalt('(ini_pi_ph: Problem with netCDF1)')
          stop '(ini_pi_ph: Problem with netCDF1)'
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
          CALL xchalt('(ini_pi_ph: Problem with netCDF200)')
          stop '(ini_pi_ph: Problem with netCDF200)'
        END IF
      END IF

      ! set missings over land
      do l=1,pi_ph_record
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j).gt.0.5) then
              pi_ph_clim(i,j,l) = pi_ph_in(i,j,l)
            else
              pi_ph_clim(i,j,l) = 0.
            endif
          enddo
        enddo
      enddo
    endif

  end subroutine ini_pi_ph


  !**********************************************************************
  ! PUBLIC SUBROUTINE : GET_PI_PH
  !
  ! Return PI_PH field for a given month.
  !**********************************************************************
  subroutine get_pi_ph(kpie,kpje,kplmon)
    use mo_control_bgc, only: with_dmsph

    implicit none

    integer, intent(in) :: kpie,kpje,kplmon
    integer, save :: oldmonth=0

    ! Ensure that pi_ph is allocated
    if(.not. allocated(pi_ph)) call alloc_pi_ph(kpie,kpje)

    ! Update only if PI pH climatology is used for DMS
    if(with_dmsph) then
      if(kplmon /= oldmonth) then
        pi_ph = reshape(pi_ph_clim(:,:,kplmon), [kpie,kpje])
        oldmonth = kplmon
      endif
    endif

  end subroutine get_pi_ph


  !**********************************************************************
  ! PRIVATE SUBROUTINE : ALLOC_PI_PH
  !
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
    ENDIF

    ALLOCATE (pi_ph(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pi_ph'
    pi_ph(:,:) = 0.0

  end subroutine alloc_pi_ph


  !**********************************************************************
  ! PRIVATE SUBROUTINE : ALLOC_PI_PH_CLIM
  !
  ! Allocate the PI_PH_CLIM field.
  !**********************************************************************
  subroutine alloc_pi_ph_clim(kpie,kpje)
    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc

    implicit none
    integer, intent(in) :: kpie,kpje
    integer             :: errstat

    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable pi_ph_clim ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',pi_ph_record
    ENDIF

    ALLOCATE (pi_ph_clim(kpie,kpje,pi_ph_record),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pi_ph_clim'
    pi_ph_clim(:,:,:) = 0.0

  end subroutine alloc_pi_ph_clim


end module mo_read_pi_ph
