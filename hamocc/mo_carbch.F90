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
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: keqb

      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: satoxy
      REAL, DIMENSION (:,:),     ALLOCATABLE :: satn2o
      REAL, DIMENSION (:,:),     ALLOCATABLE :: atdifv
      REAL, DIMENSION (:,:),     ALLOCATABLE :: suppco2
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: sedfluxo
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: dusty
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: phyto_growth
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: pi_ph
#ifdef natDIC
      REAL                                   :: atm_co2_nat
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: nathi
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: natco3
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: natOmegaC
#endif
      
      REAL :: dmspar(6)

#ifdef __c_isotopes
      REAL :: c14_t_half, c14dec
      REAL :: Roc14, Roc13
      REAL :: frac_airsea13_prom, frac_airsea13, frac_photo13_prom, frac_photo13, frac_caco313
      REAL :: frac_photo14, frac_caco314, frac_airsea14
#endif

#ifndef DIFFAT            
      REAL :: atm_co2, atm_o2, atm_n2 
      REAL :: atm_c13, atm_c14
#endif  
#ifdef CFC
      REAL :: atm_cfc11_nh,atm_cfc11_sh
      REAL :: atm_cfc12_nh,atm_cfc12_sh
      REAL :: atm_sf6_nh,atm_sf6_sh
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
        WRITE(io_stdo_bgc,*)'Memory allocation for variable pi_ph ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (pi_ph(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory pi_ph'
        pi_ph(:,:,:) = 0.0

#ifdef natDIC
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable nathi ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (nathi(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory nathi'
        nathi(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable natco3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (natco3(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory natco3'
        natco3(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable natOmegaC ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (natOmegaC(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory natOmegaC'
        natOmegaC(:,:,:) = 0.0
#endif

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
        WRITE(io_stdo_bgc,*)'Memory allocation for variable keqb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',11
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpie
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpje
        ENDIF

        ALLOCATE (keqb(11,kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory keqb'
        keqb(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (satoxy(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory satoxy'
        satoxy(:,:,:) = 0.0

#ifdef __c_isotopes
        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable phyto_growth ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        ENDIF

        ALLOCATE (phyto_growth(kpie,kpje,kpke),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory phyto_growth'
        phyto_growth(:,:,:) = 0.0
#endif

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
