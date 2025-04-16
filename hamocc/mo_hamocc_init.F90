! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen,
!                     P.-G. Chiu
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

module mo_hamocc_init

  implicit none
  private

  public :: hamocc_init

contains

  subroutine hamocc_init(read_rest,rstfnm_hamocc)

    !***********************************************************************************************
    ! Initialize HAMOCC and its interface to BLOM.
    ! Interface to ocean model (parameter list):
    ! - HAMOCC intialization when coupled to BLOM.
    !
    !  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-25
    !***********************************************************************************************

    use mod_time,       only: date,baclin
    use mod_xc,         only: ii,jj,kk,idm,jdm,kdm,nbdy,isp,ifp,ilp,mnproc,lp,xchalt
    use mod_grid,       only: plon,plat,depths
    use mod_tracers,    only: ntrbgc,ntr,itrbgc,trc
    use mo_control_bgc, only: bgc_namelist,get_bgc_namelist,do_ndep,do_rivinpt,do_oalk,            &
                              do_sedspinup,sedspin_yr_s,sedspin_yr_e,sedspin_ncyc,                 &
                              dtb,dtbgc,io_stdo_bgc,ldtbgc,                                        &
                              ldtrunbgc,ndtdaybgc,with_dmsph,l_3Dvarsedpor,use_M4AGO,              &
                              lkwrbioz_off,do_n2o_coupled,do_nh3_coupled,                          &
                              ocn_co2_type, use_sedbypass, use_BOXATM, use_BROMO,use_extNcycle,    &
                              use_coupler_ndep,lTO2depremin,use_sediment_quality,ldyn_sed_age
    use mo_param1_bgc,  only: ks,init_por2octra_mapping
    use mo_param_bgc,   only: ini_parambgc,claydens,calcdens,calcwei,opaldens,opalwei,ropal
    use mo_carbch,      only: alloc_mem_carbch,ocetra,atm,atm_co2
    use mo_biomod,      only: alloc_mem_biomod
    use mo_sedmnt,      only: alloc_mem_sedmnt,sedlay,powtra,burial,ini_sedmnt,prorca_mavg
    use mo_vgrid,       only: alloc_mem_vgrid,set_vgrid
    use mo_bgcmean,     only: alloc_mem_bgcmean
    use mo_read_rivin,  only: ini_read_rivin,rivinfile
    use mo_read_fedep,  only: ini_read_fedep,fedepfile
    use mo_read_ndep,   only: ini_read_ndep,ndepfile
    use mo_read_oafx,   only: ini_read_oafx
    use mo_read_pi_ph,  only: ini_pi_ph,pi_ph_file
    use mo_read_sedpor, only: read_sedpor,sedporfile
    use mo_read_sedqual,only: read_sedqual,sedqualfile
    use mo_clim_swa,    only: ini_swa_clim,swaclimfile
    use mo_Gdata_read,  only: inidic,inialk,inipo4,inioxy,inino3,inisil,inid13c,inid14c
    use mo_intfcblom,   only: alloc_mem_intfcblom,nphys,bgc_dx,bgc_dy,bgc_dp,bgc_rho,omask,        &
                              sedlay2,powtra2,burial2,blom2hamocc,atm2,prorca_mavg2
    use mo_ini_fields,  only: ini_fields_ocean,ini_fields_atm
    use mo_aufr_bgc,    only: aufr_bgc
    use mo_extNsediment,only: alloc_mem_extNsediment_diag
    use mo_ihamocc4m4ago, only: alloc_mem_m4ago
    use mo_m4ago_HAMOCCinit,only: init_m4ago_nml_params, init_m4ago_derived_params
    use mo_read_shelfmask,only: ini_read_shelfmask,shelfsea_maskfile


    ! Arguments
    integer,          intent(in) :: read_rest     ! flag indicating whether to read restart files.
    character(len=*), intent(in) :: rstfnm_hamocc ! restart filename.

    ! Local variables
    integer :: i,j,k,l,nt
    integer :: iounit
    real    :: sed_por(idm,jdm,ks)         = 0.
    real    :: sed_POCage_init(idm,jdm,ks) = 0.
    real    :: prorca_mavg_init(idm,jdm)   = 0.

    namelist /bgcnml/ atm_co2,fedepfile,do_rivinpt,rivinfile,do_ndep,ndepfile,do_oalk,             &
         &            do_sedspinup,sedspin_yr_s,sedspin_yr_e,sedspin_ncyc,                         &
         &            inidic,inialk,inipo4,inioxy,inino3,inisil,inid13c,inid14c,swaclimfile,       &
         &            with_dmsph,pi_ph_file,l_3Dvarsedpor,sedporfile,ocn_co2_type,use_M4AGO,       &
         &            do_n2o_coupled,lkwrbioz_off,lTO2depremin,shelfsea_maskfile,sedqualfile,      &
         &            ldyn_sed_age,do_nh3_coupled
    !
    ! --- Set io units and some control parameters
    !
    io_stdo_bgc = lp              !  standard out.
    dtbgc = nphys*baclin          !  time step length [sec].
    ndtdaybgc=NINT(86400./dtbgc)  !  time steps per day [No].
    dtb=1./ndtdaybgc              !  time step length [days].
    ldtbgc = 0
    ldtrunbgc = 0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*) 'iHAMOCC: initialisation'
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) 'restart',read_rest
      write(io_stdo_bgc,*) 'dims',idm,jdm,kdm
      write(io_stdo_bgc,*) 'date',date
      write(io_stdo_bgc,*) 'time step',dtbgc
    endif
    !
    ! --- Read the HAMOCC BGCNML namelist and check the value of some variables.
    !
    if(.not. allocated(bgc_namelist)) call get_bgc_namelist
    open (newunit=iounit, file=bgc_namelist, status='old', action='read')
    read (unit=iounit, nml=BGCNML)
    close (unit=iounit)

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*) 'iHAMOCC: reading namelist BGCNML'
      write(io_stdo_bgc,nml=BGCNML)

      if(do_sedspinup) then
        if(sedspin_yr_s<0 .or. sedspin_yr_e<0 .or. sedspin_yr_s>sedspin_yr_e) then
          call xchalt('(hamocc_init: invalid sediment spinup start/end year)')
          stop        '(hamocc_init: invalid sediment spinup start/end year)'
        endif
        if(sedspin_ncyc < 2) then
          call xchalt('(hamocc_init: invalid nb. of sediment spinup subcycles)')
          stop        '(hamocc_init: invalid nb. of sediment spinup subcycles)'
        endif
      endif
    endif

    ! init the index-mapping between pore water and ocean tracers
    call init_por2octra_mapping()
    !
    ! --- Memory allocation
    !
    call alloc_mem_intfcblom(idm,jdm,kdm)
    call alloc_mem_bgcmean(idm,jdm,kdm)
    call alloc_mem_vgrid(idm,jdm,kdm)
    call alloc_mem_biomod(idm,jdm,kdm)
    call alloc_mem_sedmnt(idm,jdm)
    call alloc_mem_carbch(idm,jdm,kdm)
    if (use_M4AGO) then
      call alloc_mem_M4AGO(idm,jdm,kdm)
    endif
    if (use_extNcycle .and. .not. use_sedbypass) then
      call alloc_mem_extNsediment_diag(idm,jdm,ks)
    endif
    !
    ! --- initialise trc array (two time levels)
    !
    do nt=itrbgc,itrbgc+ntrbgc-1
      do k=1,2*kk
        do j=1,jj
          do i=1,ii
            trc(i,j,k,nt)=0.0
          enddo
        enddo
      enddo
    enddo
    !
    ! --- initialise HAMOCC land/ocean mask
    !
    do j=1,jj
      do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          omask(i,j)=1.
        enddo
      enddo
    enddo
    !
    ! --- BLOM to HAMOCC interface
    !
    call blom2hamocc(2,1,kk,0)
    !
    ! --- Calculate variables related to the vertical grid
    !
    call set_vgrid(idm,jdm,kdm,bgc_dp)
    !
    ! --- Initialize parameters
    !
    call ini_parambgc()
    if (use_M4AGO) then
      call init_m4ago_nml_params(claydens,calcdens,calcwei,opaldens,opalwei)
      call init_m4ago_derived_params(ropal)
    endif

    ! --- Initialize atmospheric fields with (updated) parameter values
    call ini_fields_atm(idm,jdm)

    ! --- Initialize sediment and ocean tracers
    call ini_fields_ocean(read_rest,idm,jdm,kdm,nbdy,bgc_dp,bgc_rho,omask,plon,plat)

    ! --- Initialize sediment layering
    !     First, read the porosity and potentially apply it in ini_sedmnt
    call read_sedpor(idm,jdm,ks,omask,sed_por)
    !     Second, read the sediment POC age and climatological prorca and pot. apply it in ini_sedmnt
    call read_sedqual(idm,jdm,ks,omask,sed_POCage_init,prorca_mavg_init)
    call ini_sedmnt(idm,jdm,omask,sed_por,sed_POCage_init,prorca_mavg_init)

    !
    ! --- Initialise reading of input data (dust, n-deposition, river, etc.)
    !
    call ini_read_fedep(idm,jdm,omask)
    if (.not. use_coupler_ndep) then
       call ini_read_ndep(idm,jdm)
    end if
    call ini_read_rivin(idm,jdm,omask)
    call ini_read_shelfmask(idm,jdm,nbdy,depths,omask)
    call ini_read_oafx(idm,jdm,bgc_dx,bgc_dy,plat,omask)
    if (use_BROMO) then
      call ini_swa_clim(idm,jdm,omask)
    endif
    call ini_pi_ph(idm,jdm,omask)
    !
    ! --- Read restart fields from restart file if requested, otherwise
    !     (at first start-up) copy ocetra and sediment arrays (which are
    !     initialised in BELEG_VARS) to both timelevels of their respective
    !     two-time-level counterpart
    !
    if (read_rest.eq.1) then
      call AUFR_BGC(idm,jdm,kdm,ntr,ntrbgc,itrbgc,trc,                           &
           &   date%year,date%month,date%day,omask,rstfnm_hamocc)
    else
      trc(1:idm,1:jdm,1    :kdm,  itrbgc:itrbgc+ntrbgc-1) = ocetra(:,:,:,:)
      trc(1:idm,1:jdm,kdm+1:2*kdm,itrbgc:itrbgc+ntrbgc-1) = ocetra(:,:,:,:)
      if (.not. use_sedbypass) then
        sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
        sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:)
        powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
        powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:)
        burial2(:,:,1,:)         = burial(:,:,:)
        burial2(:,:,2,:)         = burial(:,:,:)
        if (use_sediment_quality) then
          prorca_mavg2(:,:,1)       = prorca_mavg(:,:)
          prorca_mavg2(:,:,2)       = prorca_mavg(:,:)
        endif
      endif
      if (use_BOXATM) then
        atm2(:,:,1,:)            = atm(:,:,:)
        atm2(:,:,2,:)            = atm(:,:,:)
      endif
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) 'iHAMOCC: finished initialisation'
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*)
    endif

  end subroutine hamocc_init

end module mo_hamocc_init
