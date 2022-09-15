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


subroutine hamocc_init(read_rest,rstfnm_hamocc)
!******************************************************************************
!
!  HAMOCC_INIT - initialize HAMOCC and its interface to BLOM.
!
!
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-25
!
!
!  Purpose
!  -------
!  - HAMOCC intialization when coupled to BLOM.
!
!
!  Interface to ocean model (parameter list):
!  -----------------------------------------
!  *INTEGER*   *read_rest*     - flag indicating whether to read restart files.
!  *INTEGER*   *rstfnm_hamocc* - restart filename.
!
!******************************************************************************
  use mod_time,       only: date,baclin
  use mod_xc,         only: ii,jj,kk,idm,jdm,kdm,nbdy,isp,ifp,ilp,              &
       &                    mnproc,lp,nfu,xchalt
  use mod_grid,       only: plon,plat
  use mod_tracers,    only: ntrbgc,ntr,itrbgc,trc
  use mo_control_bgc, only: bgc_namelist,get_bgc_namelist,                      &
       &                    do_ndep,do_rivinpt,do_oalk,do_sedspinup,            &
       &                    sedspin_yr_s,sedspin_yr_e,sedspin_ncyc,             &
       &                    dtb,dtbgc,io_stdo_bgc,ldtbgc,                       &
       &                    ldtrunbgc,ndtdaybgc,with_dmsph,l_3Dvarsedpor
  use mo_param1_bgc,  only: ks
  use mo_carbch,      only: alloc_mem_carbch,ocetra,atm,atm_co2
  use mo_biomod,      only: alloc_mem_biomod
  use mo_sedmnt,      only: alloc_mem_sedmnt,sedlay,powtra,burial
  use mo_vgrid,       only: alloc_mem_vgrid,set_vgrid
  use mo_bgcmean,     only: alloc_mem_bgcmean
  use mo_read_rivin,  only: ini_read_rivin,rivinfile
  use mo_read_fedep,  only: ini_read_fedep,fedepfile
  use mo_read_ndep,   only: ini_read_ndep,ndepfile
  use mo_read_oafx,   only: ini_read_oafx,oalkfile,oalkscen
  use mo_read_pi_ph,  only: ini_pi_ph,pi_ph_file
  use mo_read_sedpor, only: read_sedpor,sedporfile
  use mo_clim_swa,    only: ini_swa_clim,swaclimfile
  use mo_Gdata_read,  only: inidic,inialk,inipo4,inioxy,inino3,                 &
       &                    inisil,inid13c,inid14c
  use mo_intfcblom,   only: alloc_mem_intfcblom,nphys,                          &
       &                    bgc_dx,bgc_dy,bgc_dp,bgc_rho,                       &
       &                    omask,sedlay2,powtra2,burial2,                      &
       &                    blom2hamocc
#ifdef BOXATM
  use mo_intfcblom,   only: atm2
#endif

  implicit none

  integer,          intent(in) :: read_rest
  character(len=*), intent(in) :: rstfnm_hamocc

  integer :: i,j,k,l,nt
  integer :: iounit
  real    :: sed_por(idm,jdm,ks) = 0.

  namelist /bgcnml/ atm_co2,fedepfile,do_rivinpt,rivinfile,do_ndep,ndepfile,    &
       &   do_oalk,oalkscen,oalkfile,do_sedspinup,sedspin_yr_s,                 &
       &   sedspin_yr_e,sedspin_ncyc,                                           &
       &   inidic,inialk,inipo4,inioxy,inino3,inisil,                           &
       &   inid13c,inid14c,swaclimfile,                                         &
       &   with_dmsph,pi_ph_file,l_3Dvarsedpor,sedporfile
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
     WRITE(io_stdo_bgc,*)'********************************************'
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
  open (newunit=iounit, file=bgc_namelist, status='old'                         &
       &   ,action='read')
  read (unit=iounit, nml=BGCNML)
  close (unit=iounit)
  IF (mnproc.eq.1) THEN

     write(io_stdo_bgc,*)
     write(io_stdo_bgc,*) 'iHAMOCC: reading namelist BGCNML'
     write(io_stdo_bgc,nml=BGCNML)

     if(do_sedspinup) then
        if(sedspin_yr_s<0 .or. sedspin_yr_e<0 .or.                              &
             &   sedspin_yr_s>sedspin_yr_e) then
           call xchalt('(invalid sediment spinup start/end year)')
           stop        '(invalid sediment spinup start/end year)'
        endif
        if(sedspin_ncyc < 2) then
           call xchalt('(invalid nb. of sediment spinup subcycles)')
           stop        '(invalid nb. of sediment spinup subcycles)'
        endif
     endif

  ENDIF
  !
  ! --- Memory allocation
  !
  CALL ALLOC_MEM_INTFCBLOM(idm,jdm,kdm)
  CALL ALLOC_MEM_BGCMEAN(idm,jdm,kdm)
  CALL ALLOC_MEM_VGRID(idm,jdm,kdm)
  CALL ALLOC_MEM_BIOMOD(idm,jdm,kdm)
  CALL ALLOC_MEM_SEDMNT(idm,jdm)
  CALL ALLOC_MEM_CARBCH(idm,jdm,kdm)
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
  ! --- Initialize sediment layering
  ! First raed the porosity, then apply it in bodensed
  CALL read_sedpor(idm,jdm,ks,omask,sed_por)
  CALL BODENSED(idm,jdm,kdm,bgc_dp,omask,sed_por)
  !
  ! --- Initialize parameters, sediment and ocean tracer.
  !
  CALL BELEG_PARM(idm,jdm)
  CALL BELEG_VARS(read_rest,idm,jdm,kdm,nbdy,bgc_dp,bgc_rho,omask,              &
       &   plon,plat)
  !
  ! --- Initialise reading of input data (dust, n-deposition, river, etc.)
  !
  CALL ini_read_fedep(idm,jdm,omask)

  CALL ini_read_ndep(idm,jdm)

  CALL ini_read_rivin(idm,jdm,omask)

  CALL ini_read_oafx(idm,jdm,bgc_dx,bgc_dy,plat,omask)

#ifdef BROMO
  CALL ini_swa_clim(idm,jdm,omask)
#endif

  call ini_pi_ph(idm,jdm,omask)
  !
  ! --- Read restart fields from restart file if requested, otherwise
  !     (at first start-up) copy ocetra and sediment arrays (which are
  !     initialised in BELEG_VARS) to both timelevels of their respective
  !     two-time-level counterpart
  !
  IF(read_rest.eq.1) THEN
     CALL AUFR_BGC(idm,jdm,kdm,ntr,ntrbgc,itrbgc,trc,                           &
          &   date%year,date%month,date%day,omask,rstfnm_hamocc)
  ELSE
     trc(1:idm,1:jdm,1:kdm,      itrbgc:itrbgc+ntrbgc-1) =                      &
          &   ocetra(:,:,:,:)
     trc(1:idm,1:jdm,kdm+1:2*kdm,itrbgc:itrbgc+ntrbgc-1) =                      &
          &   ocetra(:,:,:,:)
#ifndef sedbypass
     sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
     sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:)
     powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
     powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:)
     burial2(:,:,1,:)         = burial(:,:,:)
     burial2(:,:,2,:)         = burial(:,:,:)
#endif
#if defined(BOXATM)
     atm2(:,:,1,:)            = atm(:,:,:)
     atm2(:,:,2,:)            = atm(:,:,:)
#endif
  ENDIF

  if (mnproc.eq.1) then
     write(io_stdo_bgc,*)
     WRITE(io_stdo_bgc,*)'********************************************'
     write(io_stdo_bgc,*) 'iHAMOCC: finished initialisation'
     write(io_stdo_bgc,*)
  endif

!******************************************************************************
end subroutine hamocc_init
