      MODULE mo_carbch


!***********************************************************************
!
!**** *MODULE mo_carbch* - Variables for inorganic carbon cycle.
!
!     S.Legutke,        *MPI-MaD, HH*     31.10.01
!
!     Modified
!     --------
!  
!     Patrick Wetzel    *MPI-Met, HH*     16.04.02
!     - new: atm, atdifv, suppco2
!     - changed: chemc(:,:,:) to chemcm(:,:,:,:)
!     - new: bgcmean(:,:,:,:)
!     J. Schwinger      *UiB-GfI, Bergen* 04.05.12
!     - added initialisation of all vars after allocation
!
!     Purpose
!     -------
!     - declaration and memory allocation
!
!     *ocetra*       *REAL*  - .
!     *hi*           *REAL*  - .
!     *co3*          *REAL*  - .
!     *chemcm*       *REAL*  - .
!     *co2flu*       *REAL*  - .
!     *co2fluacc     *REAL*  - .
!     *akw3*         *REAL*  - .
!     *akb3*         *REAL*  - .
!     *ak13*         *REAL*  - .
!     *ak23*         *REAL*  - .
!     *aksp*         *REAL*  - .
!
!**********************************************************************
      
      implicit none
      
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ocetra
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atm      
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atmflx
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: co3
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: hi
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: OmegaC 
      REAL, DIMENSION (:,:),     ALLOCATABLE :: kwb
      REAL, DIMENSION (:,:),     ALLOCATABLE :: kbb
      REAL, DIMENSION (:,:),     ALLOCATABLE :: k1b
      REAL, DIMENSION (:,:),     ALLOCATABLE :: k2b
      REAL, DIMENSION (:,:),     ALLOCATABLE :: kspb
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: satoxy
      REAL, DIMENSION (:,:),     ALLOCATABLE :: satn2o
      REAL, DIMENSION (:,:),     ALLOCATABLE :: atdifv
      REAL, DIMENSION (:,:),     ALLOCATABLE :: suppco2
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: sedfluxo
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: dusty
      
      REAL :: dmspar(6)
            
! decay coefficient for sco214 plus more for 13C/14C
      REAL :: c14dec, D14Catm, Ratm, D14Cocn, Roc14
      REAL :: roc13, atcoa
      REAL :: ozkoa, c14prod, c14inva, eins, c14ret, c14prosta

#ifndef DIFFAT            
      REAL :: atm_co2, atm_o2, atm_n2 
      REAL :: atm_c13, atm_c14
#endif  
#ifdef CFC
      REAL :: atm_cfc11,atm_cfc12,atm_sf6
#endif

#ifdef EMS_CO2
      REAL :: emission, ems_per_step
#endif
#ifdef ANTC14
      REAL :: D14C_north, D14C_south, D14C_equator
      REAL, DIMENSION (:,:), ALLOCATABLE :: Rbomb
#endif


      CONTAINS

      SUBROUTINE ALLOC_MEM_CARBCH(kpie,kpje,kpke)

      use mod_xc
      use mo_control_bgc
      use mo_param1_bgc 

      INTEGER :: kpie,kpje,kpke
      INTEGER :: errstat
      
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable ocetra ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nocetra
        ENDIF

        ALLOCATE (ocetra(kpie,kpje,kpke,nocetra),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory ocetra'
        ocetra(:,:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable hi ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (hi(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory hi'
        hi(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (co3(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory co3'
        co3(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable OmegaC ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (OmegaC(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory OmegaC'
        OmegaC(:,:,:) = 0.0


        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedfluxo ..'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',npowtra
        ENDIF

        ALLOCATE (sedfluxo(kpie,kpje,npowtra),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sedfluxo'
        sedfluxo(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (satn2o(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory satn2o'
        satn2o(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable kspb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (kspb(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory kspb'
        kspb(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (satoxy(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory satoxy'
        satoxy(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable k1b ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (k1b(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory k1b'
        k1b(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable k2b ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (k2b(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory k2b'
        k2b(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable kbb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (kbb(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory kbb'
        kbb(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable kwb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (kwb(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory kwb'
        kwb(:,:) = 0.0

!#if defined(DIFFAT) || defined(CCSMCOUPLED)

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atm ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',natm
        ENDIF

        ALLOCATE (atm(kpie,kpje,natm),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory atm'
        atm(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atmflx ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',natm
        ENDIF

        ALLOCATE (atmflx(kpie,kpje,natm),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory atmflx'
        atmflx(:,:,:) = 0.0

!#endif	

#if defined(DIFFAT)
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atdifv ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (atdifv(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory atdifv'
        atdifv(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable suppco2 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (suppco2(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory suppco2'
        suppco2(:,:) = 0.0
#endif	

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable dusty ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',12
        ENDIF

        ALLOCATE (dusty(kpie,kpje,12),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory dusty'
        dusty(:,:,:) = 0.0

#ifdef ANTC14
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable Rbomb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF
        
        ALLOCATE (Rbomb(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory Rbomb'
        Rbomb(:,:) = 0.0
#endif	

      END SUBROUTINE ALLOC_MEM_CARBCH

      END MODULE mo_carbch
