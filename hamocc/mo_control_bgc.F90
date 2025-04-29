! Copyright (C) 2002  S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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

module mo_control_bgc

  !*************************************************************************************************
  ! Control variables for iHAMOCC.
  !
  ! S.Legutke,        *MPI-MaD, HH*    28.02.02
  !
  ! Modified
  ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  ! - removed unused variables
  !*************************************************************************************************

  implicit none
  public

  ! Routines
  public :: get_bgc_namelist

  ! Logical unit number for I/O.
  integer :: io_stdo_bgc              !  standard out.

  ! File containing namelists
  character(len=:), allocatable, protected :: bgc_namelist

  ! Control variables
  real    :: dtbgc                    !  time step length [sec].
  real    :: dtb                      !  time step length [days].
  integer :: ndtdaybgc                !  time steps per day.

  integer :: ldtbgc                   !  time step number from bgc restart file
  integer :: ldtrunbgc                !  actual time steps of run.

  real    :: rmasks = 0.0             !  value at wet cells in sediment.
  real    :: rmasko = 99999.00        !  value at wet cells in ocean.

  ! Variables set via namelist bgcnml
  logical           :: l_3Dvarsedpor          = .false. ! apply spatially variable sediment porosity
  logical           :: do_ndep                = .true.  ! apply n-deposition
  logical           :: do_n2onh3_coupled      = .false. ! for coupled simulations, use field provided by atmosphere
  logical           :: do_rivinpt             = .true.  ! apply riverine input
  logical           :: do_sedspinup           = .false. ! apply sediment spin-up
  logical           :: do_oalk                = .false. ! apply ocean alkalinization
  logical           :: with_dmsph             = .false. ! apply DMS with pH dependence
  logical           :: use_M4AGO              = .false. ! run with M4AGO settling scheme
  logical           :: lkwrbioz_off           = .false. ! if true, allow remin and primary prod throughout full water column
  logical           :: lTO2depremin           = .false. ! Temperature- and O2-dependent remineralization of POM
  logical           :: ldyn_sed_age           = .false. ! switch for dynamic sediment age in combination with use_sediment_quality
  integer           :: sedspin_yr_s           = -1      ! start year for sediment spin-up
  integer           :: sedspin_yr_e           = -1      ! end   year for sediment spin-up
  integer           :: sedspin_ncyc           = -1      ! sediment spin-up sub-cycles
  character(len=64) :: ocn_co2_type                     ! indicates co2 coupling to an active atm
                                                        ! model if set to 'diagnostic'
                                                        ! or 'prognostic'

  ! Logical switches set via namelist config_bgc
  logical           :: use_BROMO              = .false.
  logical           :: use_AGG                = .false.
  logical           :: use_WLIN               = .true.
  logical           :: use_natDIC             = .false.
  logical           :: use_CFC                = .false.
  logical           :: use_cisonew            = .false.
  logical           :: use_PBGC_OCNP_TIMESTEP = .false.
  logical           :: use_PBGC_CK_TIMESTEP   = .false.
  logical           :: use_FB_BGC_OCE         = .false.
  logical           :: use_BOXATM             = .false.
  logical           :: use_sedbypass          = .false.
  logical           :: use_extNcycle          = .false.
  logical           :: use_coupler_ndep       = .false.
  logical           :: use_pref_tracers       = .true.
  logical           :: use_shelfsea_res_time  = .false.
  logical           :: use_sediment_quality   = .false.

contains

  subroutine get_bgc_namelist
    !-------------------------------------------------------------------------
    ! Get filename for namelist file
    !-------------------------------------------------------------------------
    use mod_config, only: inst_suffix
    use mod_xc,     only: xchalt

    logical :: exists

    if (.not. allocated(bgc_namelist)) then
      inquire (file='ocn_in'//trim(inst_suffix), exist=exists)
      if (exists) then
        allocate(character(len=len('ocn_in'//trim(inst_suffix))) :: bgc_namelist)
        bgc_namelist = 'ocn_in'//trim(inst_suffix)
      else
        inquire (file='limits', exist=exists)
        if (exists) then
          allocate(character(len=len('limits')) :: bgc_namelist)
          bgc_namelist = 'limits'
        else
          call xchalt('cannot find limits file')
          stop        'cannot find limits file'
        endif
      endif
    endif
  end subroutine get_bgc_namelist

end module mo_control_bgc
