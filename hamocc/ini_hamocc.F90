! Copyright (C) 2001  S. Legutke
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


      SUBROUTINE INI_HAMOCC(kpaufr,kpie,kpje,kpke,kbnd,                       &
                            kplyear,kplmonth,kplday,                          &
                            pdt,pdlxp,pdlyp,pddpo,prho,omask,pglon,pglat,     &
                            ntr,ntrbgc,itrbgc,trc,                            &
#ifndef sedbypass
                            sedlay2,powtra2,burial2,                          &    
#endif
                            rstfnm_ocn,path)

!****************************************************************
!
!**** *INI_BGC* - initialize marine bio-geo-chemistry module.
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!     J.Schwinger       *GFI, Bergen*    2013-10-21
!     - added GNEWS2 option for riverine input of carbon and nutrients
!     - code cleanup
!
!     J.Schwinger,      *GFI, Bergen*    2014-05-21
!     - adapted code for use with two time level tracer field in MICOM
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved reading of namelist and initialisation of dust to 
!       ini_hamocc.F90
!     - added sediment bypass preprocessor option
!     
!     Purpose
!     -------
!     - allocate and initialize bgc variables
!     - allocate and initialize sediment layering
!     - initialize bgc constants (beleg.F90)
!     - read restart fields if kpaufr=1
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpaufr*     - 1/0 for read / do not read restart file
!     *INTEGER* *kpie*       - zonal dimension of model grid.
!     *INTEGER* *kpje*       - meridional dimension of model grid.
!     *INTEGER* *kpke*       - vertical dimension of model grid.
!     *INTEGER* *kbnd*       - nb of halo grid points
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmonth*   - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *REAL*    *pdt*        - ocean model time step [sec].
!     *REAL*    *pdlxp*      - size of scalar grid cell (zonal) [m].
!     *REAL*    *pdlyp*      - size of scalar grid cell (meridional) [m].
!     *REAL*    *pddpo*      - size of scalar grid cell (vertical dimension) [m].
!     *REAL*    *prho*       - density [g/cm^3].
!     *REAL*    *omask*      - land/ocean mask
!     *REAL*    *pglon*      - geographical longitude of grid points [degree E].
!     *REAL*    *pglat*      - geographical latitude  of grid points [degree N].
!     *INTEGER* *ntr*        - number of tracers in tracer field
!     *INTEGER* *ntrbgc*     - number of biogechemical tracers in tracer field
!     *INTEGER* *itrbgc*     - start index for biogeochemical tracers in tracer field
!     *REAL*    *trc*        - initial/restart tracer field to be passed to the 
!                              ocean model [mol/kg]
!     *REAL*    *sedlay2*    - initial/restart sediment (two time levels) field
!     *REAL*    *powtra2*    - initial/restart pore water tracer (two time levels) field
!     *REAL*    *burial2*    - initial/restart sediment burial (two time levels) field
!     *CHAR*    *rstfnm_ocn* - restart file name-informations
!     *CHAR*    *path*       - path to data files
!
!**********************************************************************
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc
      use mo_vgrid, only: alloc_mem_vgrid,set_vgrid
      use mod_xc, only: mnproc,lp,nfu
      use mo_bgcmean
      use mo_ndep, only: ini_ndep,ndepfname
      use mo_riverinpt, only: ini_riverinpt

 
      implicit none

      INTEGER, intent(in)  :: kpaufr
      INTEGER, intent(in)  :: kpie,kpje,kpke,kbnd,ntr,ntrbgc,itrbgc
      INTEGER, intent(in)  :: kplyear,kplmonth,kplday      
      REAL,    intent(in)  :: pdt
      REAL,    intent(in)  :: pdlxp(kpie,kpje)
      REAL,    intent(in)  :: pdlyp(kpie,kpje)
      REAL,    intent(in)  :: pddpo(kpie,kpje,kpke)
      REAL,    intent(in)  :: prho (kpie,kpje,kpke)
      REAL,    intent(in)  :: omask(kpie,kpje)
      REAL,    intent(in)  :: pglon(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: pglat(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: trc(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,2*kpke,ntr)
#ifndef sedbypass
      REAL,    intent(out) :: sedlay2(kpie,kpje,2*ks,nsedtra)
      REAL,    intent(out) :: powtra2(kpie,kpje,2*ks,npowtra)
      REAL,    intent(out) :: burial2(kpie,kpje,2,   nsedtra)
#endif
      character(len=*), intent(in) :: rstfnm_ocn,path

      ! Local variables
      INTEGER :: i,j,k,l

      namelist /bgcnml/ atm_co2,do_rivinpt,do_ndep,ndepfname


!
! Define io units
!
      io_stdo_bgc = lp      !  standard out.
      io_nml = nfu          !  namelist

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) 'HAMOCC initialisation'
      write(io_stdo_bgc,*) 'restart',kpaufr
      write(io_stdo_bgc,*) 'dims',kpie,kpje,kpke
      write(io_stdo_bgc,*) 'time',kplyear,kplmonth,kplday
      write(io_stdo_bgc,*) 'time step',pdt
      endif

!
! Read the HAMOCC BGCNML namelist.
!
      open (unit=io_nml,file='ocn_in',status='old',action='read')
      read (unit=io_nml,nml=BGCNML)
      close (unit=io_nml)
      IF (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)
        write(io_stdo_bgc,*) 'HAMOCC: reading namelist BGCNML'
        write(io_stdo_bgc,nml=BGCNML)
      ENDIF

!                    
! Set control constants ( mo_control_bgc )
!
      dtbgc = pdt                   !  time step length [sec].
      ndtdaybgc=NINT(86400./dtbgc)  !  time steps per day [No].
      dtb=1./ndtdaybgc              !  time step length [days].

!
! Initialize some namelist parameters
!
      isac = 1

!
! Initialize time step counter of run.
!
      ldtrunbgc = 0

!                        
! Allocate memory : model output
!
      CALL ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

!                        
! Allocate memory : vgrid
!
      CALL ALLOC_MEM_VGRID(kpie,kpje,kpke)

!                        
! Allocate memory : biology
!
      CALL ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

!                        
! Allocate memory : sediment
!
      CALL ALLOC_MEM_SEDMNT(kpie,kpje)

!                        
! Allocate memory : inorganic carbon cycle
!
      CALL ALLOC_MEM_CARBCH(kpie,kpje,kpke)

!
! Calculate variables related to the vertical grid
!
      call set_vgrid(kpie,kpje,kpke,pddpo)

!                        
! Initialize sediment layering
!
      CALL BODENSED(kpie,kpje,kpke,pddpo)

!                        
! Initialize sediment and ocean tracer.
! 
      CALL BELEG_BGC(kpaufr,kpie,kpje,kpke,kbnd,pddpo,prho,omask,        &
                     pglon,pglat,path)

!
! Initialise dust input, n-deposition and river input
!
      CALL GET_DUST(kpie,kpje,kpke,omask,path)

      CALL ini_ndep(kpie,kpje,path)

      CALL INI_RIVERINPT(path)
#ifdef DMSPH
      CALL GET_PI_PH(kpie,kpje,kpke,omask,path)
#endif
!                        
! Read restart fields from restart file if requested, otherwise 
! (at first start-up) copy ocetra and sediment arrays (which are
! initialised in BELEG) to both timelevels of their respective
! two-time-level counterpart
!
      IF(kpaufr.eq.1) THEN
         CALL AUFR_BGC(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,trc,             &
#ifndef sedbypass
     &                 sedlay2,powtra2,burial2,                          &
#endif
     &                 kplyear,kplmonth,kplday,omask,rstfnm_ocn)
      ELSE
         trc(1:kpie,1:kpje,1:kpke,       itrbgc:itrbgc+ntrbgc-1) =       &
     &     ocetra(:,:,:,:)
         trc(1:kpie,1:kpje,kpke+1:2*kpke,itrbgc:itrbgc+ntrbgc-1) =       &
     &     ocetra(:,:,:,:)
#ifndef sedbypass
         sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
         sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:) 
         powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
         powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:) 
         burial2(:,:,1,:)         = burial(:,:,:)
         burial2(:,:,2,:)         = burial(:,:,:) 
#endif
#if defined(BOXATM) || defined(DIFFAT)
         atm2(:,:,1,:)            = atm(:,:,:)
         atm2(:,:,2,:)            = atm(:,:,:)
#endif
      ENDIF


!
! Global inventory of all tracers
!      
!      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,1)


      RETURN
      END
