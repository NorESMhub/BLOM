! Copyright (C) 2020  J. Schwinger
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


module mo_read_fedep

  !*************************************************************************************************
  ! Declaration, memory allocation, and routines related to reading iron deposition input data
  !
  !  J.Schwinger,      *NORCE Climate, Bergen*   2020-05-27
  !
  !  Modified
  !  J. Schwinger,     *NORCE climate, Bergen*   2022-06-02
  !   -revise structure of this module, split into a module for reading the data (mo_read_fedep)
  !    and a module that applies the fluxes in core hamocc (mo_apply_fedep)
  !*************************************************************************************************

  use mo_kind, only: bgc_fnmlen

  implicit none
  private

  public :: ini_read_fedep ! Initialise the module for reading iron deposition data
  public :: get_fedep      ! Get the iron (dust) deposition for a given month

  ! File name (incl. full path) for input data, set through namelist in hamocc_init
  character(len=bgc_fnmlen), public :: fedepfile=''
  character(len=64),         public :: fedep_source=''

  ! Array to store iron deposition fluxes after reading from file
  real, allocatable,  private :: dustflx_tot(:,:,:)
  real, allocatable,  private :: dustflx_sfe(:,:,:)

contains

  subroutine ini_read_fedep(kpie,kpje,omask)

    !***********************************************************************************************
    ! Initialise the iron deposition module, read in the iron (dust) data set.
    !
    ! J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
    !***********************************************************************************************

    use netcdf,             only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mod_xc,             only: mnproc,xchalt
    use mo_control_bgc,     only: io_stdo_bgc
    use mo_param_bgc,       only: sec_per_day,frac_ironindust,frac_soliron,fetune
    use mo_chemcon,         only: mw_fe
    use mo_netcdf_bgcrw,    only: read_netcdf_var

    ! Arguments
    integer,          intent(in) :: kpie              ! 1st dimension of model grid.
    integer,          intent(in) :: kpje              ! 2nd dimension of model grid.
    real,             intent(in) :: omask(kpie,kpje)  ! land/ocean mask (1=ocean)

    ! Local variables
    integer             :: i,j,l
    integer             :: ncid,ncstat,ncvarid,errstat

    ! allocate field to hold iron deposition fluxes
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_fedep:'
      write(io_stdo_bgc,*)' '
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable dustflx_tot ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    :  12'
    endif
    allocate (dustflx_tot(kpie,kpje,12),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory dustflx_tot'
    dustflx_tot(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable dustflx_sfe ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    :  12'
    endif
    allocate (dustflx_sfe(kpie,kpje,12),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory dustflx_sfe'
    dustflx_sfe(:,:,:) = 0.0

    ! Open netCDF data file
    if (mnproc==1) then
      ncstat = NF90_OPEN(trim(fedepfile),NF90_NOWRITE, ncid)
      if (ncstat /= NF90_NOERR ) then
        call xchalt('(ini_read_fedep: Problem with netCDF1)')
        stop        '(ini_read_fedep: Problem with netCDF1)'
      end if
    end if

    ! Read  data
    select case(trim(fedep_source))

      case('mahw2006')
        ! Dust/iron deposition using the Mahowald et al. 2006 model data.
        ! The variable 'DUST' contains the total dust-flux in kg/m2/month; this
        ! is converted to kmol fe/m2/s of bioavailable (soluble) iron using the molar
        ! weight of iron (mw_fe=55.85) and assuming that dust contains a fraction
        ! 'frac_ironindust' of iron, and that soluble iron is a fraction 'frac_soliron' of
        ! total iron. A tuning factor 'fetune' is applied to the soluble fraction of iron.
        ! Note the conversion kg/m2/month -> kg/m2/s is assuming 30 days per month.
        call read_netcdf_var(ncid,'DUST',dustflx_tot(1,1,1),12,0,0)
        dustflx_tot(:,:,:) = dustflx_tot(:,:,:)/30.0/sec_per_day
        dustflx_sfe(:,:,:) = dustflx_tot(:,:,:)*frac_ironindust*frac_soliron/mw_fe*fetune

      case('GESAMP2018')
        ! Iron deposition using data from the GESAMP atmospheric iron deposition
        ! model intercomparison study (Myriokefalitakis et al. 2018,
        ! doi:10.5194/bg-15-6659-2018).
        ! The variables 'TFe' and 'LFe' contain the total and soluble (labile) iron
        ! fluxes in kg/m2/s; soluble iron is converted to kmol Fe/m2/s using the molar weight of
        ! iron (mw_fe=55.85). A tuning factor 'fetune' is applied to the soluble fraction of iron.
        ! The total dust-flux is 'back-calculated' from total iron assuming a fixed fraction
        ! of iron in dust.
        call read_netcdf_var(ncid,'TFe',dustflx_tot(1,1,1),12,0,0)
        call read_netcdf_var(ncid,'LFe',dustflx_sfe(1,1,1),12,0,0)
        dustflx_tot(:,:,:) = dustflx_tot(:,:,:)/frac_ironindust
        dustflx_sfe(:,:,:) = dustflx_sfe(:,:,:)/mw_fe*fetune

      case default
        if (mnproc==1) then
          call xchalt('(ini_read_fedep: invalid fedep_source)')
          stop        '(ini_read_fedep: invalid fedep_source)'
        end if

    end select

    ! Close file
    if (mnproc==1) then
      ncstat = NF90_CLOSE(ncid)
      if ( ncstat /=  NF90_NOERR ) then
        call xchalt('(ini_read_fedep: Problem with netCDF2)')
        stop        '(ini_read_fedep: Problem with netCDF2)'
      end if
    end if

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_fedep: Using dust deposition file '//trim(fedepfile)
    endif

    ! set flux to zero over land
    do l=1,12
      do j=1,kpje
        do i=1,kpie

          if(omask(i,j).lt.0.5) then
            dustflx_tot(i,j,l) = 0.0
            dustflx_sfe(i,j,l) = 0.0
          endif

        enddo
      enddo
    enddo

  end subroutine ini_read_fedep


  subroutine get_fedep(kpie,kpje,kbnd,kplmon,dust,dust_stream)

    !***********************************************************************************************
    ! Get iron (dust) deposition for current month
    !
    !  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
    !***********************************************************************************************

    use mod_xc             , only: mnproc,xchalt,nbdy
    use mo_output_forcing  , only: output_forcing
    use mo_param1_bgc      , only: itdust,isfe,ndust
    use mo_param_bgc       , only: sec_per_day,frac_ironindust,frac_soliron,fetune
    use mo_chemcon         , only: mw_fe

    integer,        intent(in)  :: kpie                   ! 1st dimension of model grid
    integer,        intent(in)  :: kpje                   ! 2nd dimension of model grid
    integer,        intent(in)  :: kbnd                   ! nb of halo-rows
    integer,        intent(in)  :: kplmon                 ! current month.
    real,           intent(out) :: dust(kpie,kpje,ndust)  ! dust/sFe flux for current month
    real, optional, intent(in)  :: dust_stream(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,ndust)

    integer :: i,j,n
    logical :: debug = .true.
    logical :: first_time = .true.

    if (present(dust_stream)) then

      select case(trim(fedep_source))

        case('mahw2006')
          do j=1,kpje
            do i=1,kpie
              ! note that the stream is read from a file where units already are in kg/m2/s
              dust(i,j,itdust) = dust_stream(i,j,1)
              dust(i,j,isfe)   = dust_stream(i,j,1)*frac_ironindust*frac_soliron/mw_fe*fetune
            end do
          end do

        case('GESAMP2018')
          do j=1,kpje
            do i=1,kpie
              ! note that the stream is read from a file where units already are in kg/m2/s
              dust(i,j,itdust) = dust_stream(i,j,1)/frac_ironindust
              dust(i,j,isfe)   = dust_stream(i,j,2)/mw_fe*fetune
            end do
          end do

        case default
          if (mnproc==1) then
            call xchalt('(get_fedep: invalid fedep_source for dust stream)')
            stop        '(get_fedep: invalid fedep_source for dust stream)'
          end if

      end select

    else
      dust(:,:,itdust) = dustflx_tot(:,:,kplmon)
      dust(:,:,isfe)   = dustflx_sfe(:,:,kplmon)
    end if

    if (debug) then
       if (first_time) then
          call output_forcing('dustflx_tot.nc', 'dustflx_tot', kpie, kpje, dust(:,:,itdust))
          call output_forcing('dustflx_sfe.nc', 'dustflx_sfe', kpie, kpje, dust(:,:,isfe))
          first_time = .false.
       end if
    end if

  end subroutine get_fedep

end module mo_read_fedep
