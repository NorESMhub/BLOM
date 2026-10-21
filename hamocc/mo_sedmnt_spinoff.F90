! Copyright (C) 2021-2024  j. maerz
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


module mo_sedmnt_spinoff

  !*************************************************************************************************
  ! Routine to perform an offline sediment spinup
  !
  ! j.maerz *UiB, Bergen* 2025-04-11
  !
  !*************************************************************************************************


  implicit none
  private

  public :: sed_offline_spinup

contains

  subroutine sed_offline_spinup(kpie,kpje,kpke,kbnd,omask,kplyear,kldtday,pddpo,  prho,psao,ptho)
    use mod_xc,             only: mnproc
    use mo_control_bgc,     only: io_stdo_bgc,ldtbgc,dtbgc
    use mo_powach,          only: powach
    use mo_sedshi,          only: sedshi
    use mo_vgrid,           only: dp_min,kbo,ptiestu
    use mo_sedmnt,          only: prorca,prcaca,produs,silpro
    use mo_read_sedspinoff, only: clim_prorca,clim_prcaca,clim_produs,clim_silpro
    use mo_carchm,          only: carchm_kequi
    use mo_carbch,          only: keqb
    use mo_param_bgc,       only: rcar

    ! Arguments
    integer, intent(in)  :: kpie                                            ! 1st dimension of model grid.
    integer, intent(in)  :: kpje                                            ! 2nd dimension of model grid.
    integer, intent(in)  :: kpke                                            ! 3rd (vertical) dimension of model grid
    integer, intent(in)  :: kbnd                                            ! number of halo grid points.
    integer, intent(in)  :: kplyear                                         ! current year.
    integer, intent(in)  :: kldtday                                         ! number of time step in current day.
    real,    intent(in)  :: omask(kpie,kpje)                                ! land/ocean mask.
    real,    intent(in)  :: pddpo(kpie,kpje,kpke)                           ! size of grid cell (3rd dimension) [m]
    ! The following SH/COULD BE CHANGED TO MODEL INPUT VALUES (yr-mean climatology)!?
    ! Can be provided through the BLOM restart file...
    real,    intent(in)  :: prho   (kpie,kpje,kpke)                         ! density [g/cm^3].
    real,    intent(in)  :: psao   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) ! salinity [psu.].
    real,    intent(in)  :: ptho   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) ! potential temperature [deg C].

    ! Minimum and maximum SST/SSS set for carbon chemistry and gas exchange calculations
    ! parameters as in mo_carchm
    real, parameter :: temp_min = -1.0
    real, parameter :: temp_max = 40.0
    real, parameter :: saln_min =  5.0
    real, parameter :: saln_max = 40.0

    ! Local variables
    integer       :: i,j
    real          :: Kh0,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa,prb,t,s
    logical       :: lspin = .true.
    logical, save :: init_offl_spinup = .true.

    !--------------------------------------------------------------------
    ! Increment bgc time step counter of experiment (it seems required for restart writing)
    ldtbgc = ldtbgc + 1

    if (init_offl_spinup) then
      ! At first call, initialize the forcing fields and provide bottom water chemistry constants
      prorca = 0.
      prcaca = 0.
      produs = 0.
      silpro = 0.
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'Initializing bottom equilibrium constants ...'
      endif
      !$OMP PARALLEL DO PRIVATE(i,j,Kh0,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa,prb,t,s)
      do i=1,kpie
        do j=1,kpje
          if (omask(i,j) > 0.5 .and. pddpo(i,j,kbo(i,j)) > dp_min) then
            t    = min(temp_max,max(temp_min,ptho(i,j,kbo(i,j))))
            s    = min(saln_max,max(saln_min,psao(i,j,kbo(i,j))))
            prb  = ptiestu(i,j,kbo(i,j))*98060.0*1.027e-6 ! pressure in unit bars, 98060 = onem
            call carchm_kequi(t,s,prb,Kh0,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa)

            keqb( 1,i,j)  = K1
            keqb( 2,i,j)  = K2
            keqb( 3,i,j)  = Kb
            keqb( 4,i,j)  = Kw
            keqb( 5,i,j)  = Ks1
            keqb( 6,i,j)  = Kf
            keqb( 7,i,j)  = Ksi
            keqb( 8,i,j)  = K1p
            keqb( 9,i,j)  = K2p
            keqb(10,i,j)  = K3p
            keqb(11,i,j)  = Kspc

            ! Initializing the climatological deposition fluxes
            prorca(i,j) = clim_prorca(i,j)/rcar*1e-3*dtbgc ! mol C m-2 s-1  -> kmol P m-2/time step
            prcaca(i,j) = clim_prcaca(i,j)     *1e-3*dtbgc ! mol Ca m-2 s-1 -> kmol Ca m-2/time step
            produs(i,j) = clim_produs(i,j)     *1e-3*dtbgc ! g m-2 s-1      -> kg m-2/time step
            silpro(i,j) = clim_silpro(i,j)     *1e-3*dtbgc ! mol Si m-2 s-1 -> kmol Si m-2/time step
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO

      ! set logical to false after initialization
      init_offl_spinup = .false.
    endif

    call powach(kpie,kpje,kpke,kbnd,prho,omask,psao,ptho,lspin)

    if (kldtday == 1 .or. kldtday == 2) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*) 'Sediment shifting ...'
      endif
      call sedshi(kpie,kpje,omask,kplyear)
    endif

  end subroutine sed_offline_spinup

end module mo_sedmnt_spinoff
