      MODULE mo_carbch


!***********************************************************************
!
!**** *MODULE mo_carbch* - Variables for inorganic carbon cycle.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!  
!     Patrick Wetzel    *MPI-Met, HH*    16.04.02
!     - new: atm, atdifv, suppco2
!     - changed: chemc(:,:,:) to chemcm(:,:,:,:)
!     - new: bgcmean(:,:,:,:)
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
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: chemcm
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atm      
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atmflx
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: co3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: hi
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: OmegaC 
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: akw3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: akb3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ak13
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ak23
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: aksp
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: satoxy
      REAL, DIMENSION (:,:),   ALLOCATABLE :: satn2o
      REAL, DIMENSION (:,:),   ALLOCATABLE :: atdifv
      REAL, DIMENSION (:,:),   ALLOCATABLE :: suppco2
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: sedfluxi
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: sedfluxo
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: dusty
      
      REAL :: dmspar(6)
            
! decay coefficient for sco214 plus more for 13C/14C
      REAL :: c14dec, D14Catm, Ratm, D14Cocn, Roc14
      REAL :: roc13, atcoa
      REAL :: ozkoa, c14prod, c14inva, eins, c14ret, c14prosta

#ifndef DIFFAT            
      REAL :: atm_co2, atm_o2, atm_n2 
      REAL :: atm_c13, atm_c14, atmacon,atmacmol
#endif  

#ifdef EMS_CO2
!#ifdef DIFFAT
      REAL :: emission, ems_per_step
!#else
      REAL :: co2_atm_1,co2_atm_2,co2_atm_3
!#endif      
#endif
#ifdef PCFC
      REAL ::           cfc11_atm_1s,cfc11_atm_2s,cfc11_atm_3s        &
     &                 ,cfc11_atm_1n,cfc11_atm_2n,cfc11_atm_3n        &
     &                 ,cfc12_atm_1s,cfc12_atm_2s,cfc12_atm_3s        &
     &                 ,cfc12_atm_1n,cfc12_atm_2n,cfc12_atm_3n        

      REAL, DIMENSION (:,:), ALLOCATABLE :: cfc_int
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

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable hi ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (hi(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory hi'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (co3(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory co3'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable OmegaC ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (OmegaC(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory OmegaC'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable chemcm ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',8
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',12
        ENDIF

        ALLOCATE (chemcm(kpie,kpje,8,12),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory chemcm'

!        IF (mnproc.eq.1) THEN
!        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedfluxi ...'
!        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
!        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
!	 WRITE(io_stdo_bgc,*)'Third dimension    : ',3
!        ENDIF
!
!        ALLOCATE (sedfluxi(kpie,kpje,3))

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedfluxo ..'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',npowtra
        ENDIF

        ALLOCATE (sedfluxo(kpie,kpje,npowtra),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sedfluxo'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (satn2o(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory satn2o'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable aksp ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (aksp(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory aksp'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (satoxy(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory satoxy'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable ak23 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (ak23(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory ak23'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable ak13 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (ak13(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory ak13'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable akb3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (akb3(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory akb3'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable akw3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (akw3(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory akw3'
!#if defined(DIFFAT) || defined(CCSMCOUPLED)

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atm ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',natm
        ENDIF

        ALLOCATE (atm(kpie,kpje,natm),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory atm'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atmflx ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',natm
        ENDIF

        ALLOCATE (atmflx(kpie,kpje,natm),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory atmflx'
!#endif	
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atdifv ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (atdifv(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory atdifv'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable suppco2 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (suppco2(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory suppco2'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable dusty ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',12
        ENDIF

        ALLOCATE (dusty(kpie,kpje,12),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory dusty'

#ifdef PCFC
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable cfc_int ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF
        
        ALLOCATE (cfc_int(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory cfc_int'
#endif
#ifdef ANTC14
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable Rbomb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF
        
        ALLOCATE (Rbomb(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory Rbomb'
#endif	

      END SUBROUTINE ALLOC_MEM_CARBCH

      END MODULE mo_carbch
