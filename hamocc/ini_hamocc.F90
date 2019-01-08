      SUBROUTINE INI_HAMOCC(kpaufr,pdt,kpndtrun,kpie,kpje,kpke,kpbe      &
     &            ,pddpo,ptho,psao,prho,pdlxp,pdlyp,ptiestu,ptiestw      &
     &            ,kplyear,kplmonth,kplday,kpldtoce                      &
     &            ,pglon,pglat,omask,ntr,ntrbgc,itrbgc,trc               &
#ifndef sedbypass
     &            ,sedlay2,powtra2,burial2                               &    
#endif
     &            ,rstfnm_ocn,path)

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
!     *REAL*    *pdt*        - ocean model time step [sec].
!     *INTEGER* *kpndtrun*   - total no. of time steps of run.
!     *INTEGER* *kpie*       - zonal dimension of model grid.
!     *INTEGER* *kpje*       - meridional dimension of model grid.
!     *INTEGER* *kpke*       - vertical dimension of model grid.
!     *REAL*    *pddpo*      - size of scalar grid cell (vertical dimension) [m].
!     *REAL*    *ptho*       - potential temperature [deg C].
!     *REAL*    *psao*       - salinity [psu].
!     *REAL*    *prho*       - density [g/cm^3].
!     *REAL*    *pdlxp*      - size of scalar grid cell (zonal) [m].
!     *REAL*    *pdlyp*      - size of scalar grid cell (meridional) [m].
!     *REAL*    *ptiestu*    - depth of level [m].
!     *REAL*    *ptiestw*    - depth of level interface [m].
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmonth*   - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *INTEGER* *kpldtoce*   - step  in ocean restart date
!     *REAL*    *pglon*      - geographical longitude of grid points [degree E].
!     *REAL*    *pglat*      - geographical latitude  of grid points [degree N].
!     *REAL*    *omask*      - land/ocean mask
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
      use mod_xc, only: mnproc,lp,nfu
      use mo_bgcmean
      use mo_ndep, only: ini_ndep,ndepfname
      use mo_riverinpt, only: ini_riverinpt

 
      implicit none
      INTEGER :: kpie,kpje,kpke,kpbe,ntr,ntrbgc,itrbgc
      INTEGER :: kplyear,kplmonth,kplday,kpldtoce
      INTEGER :: kpaufr,kpndtrun,k,l
      INTEGER :: i,j
      
      REAL    :: pddpo(kpie,kpje,kpke)
      REAL    :: ptho (kpie,kpje,kpke)
      REAL    :: psao (kpie,kpje,kpke)
      REAL    :: prho (kpie,kpje,kpke)
      REAL    :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL    :: pglon(kpie,kpje),pglat(kpie,kpje)
      REAL    :: ptiestu(kpie,kpje,kpke+1),ptiestw(kpie,kpje,kpke+1)
      REAL    :: omask(kpie,kpje)
      REAL    :: trc(1-kpbe:kpie+kpbe,1-kpbe:kpje+kpbe,2*kpke,ntr)
#ifndef sedbypass
      REAL    :: sedlay2(kpie,kpje,2*ks,nsedtra)
      REAL    :: powtra2(kpie,kpje,2*ks,npowtra)
      REAL    :: burial2(kpie,kpje,2,   nsedtra)
#endif
      REAL    :: pdt
      character(len=*) :: rstfnm_ocn,path

      namelist /bgcnml/ atm_co2,do_rivinpt,do_ndep,ndepfname


!
! Define io units
!
      io_stdo_bgc = lp      !  standard out.
      io_nml = nfu          !  namelist

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) 'HAMOCC initialisation'
      write(io_stdo_bgc,*) 'restart',kpaufr,kpndtrun
      write(io_stdo_bgc,*) 'dims',kpie,kpje,kpke
      write(io_stdo_bgc,*) 'time',kplyear,kplmonth,kplday,kpldtoce
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
      
      ndtrunbgc = kpndtrun

!
! Initialize some namelist parameters
!
      isac = 1

!
! Initialize time step counter of run.
!
      ldtrunbgc = 0

      CALL ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

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
! Initialize sediment layering
!
      CALL BODENSED(kpie,kpje,kpke,pddpo)

!                        
! Initialize sediment and ocean tracer.
! 
      CALL BELEG_BGC(kpaufr,kpie,kpje,kpke,pddpo,ptiestw,prho,           &
     &               omask,pglon,pglat,path)
     
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
     &                 kplyear,kplmonth,kplday,kpldtoce,omask,           &
     &                 rstfnm_ocn)
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
      ENDIF


!
! Global inventory of all tracers
!      
!      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,1)


      RETURN
      END
