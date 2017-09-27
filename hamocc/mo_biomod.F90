      MODULE mo_biomod

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/mo_biomod.f90,v $\\
!$Revision: 1.2 $\\
!$Date: 2004/11/12 15:37:21 $\\

!***********************************************************************
!
!**** *MODULE mo_biomod* - Variables for marine biology.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!     note: kbo,bolay shall be moved to mo_control_bgc
!     
!     I. Kriest, GEOMAR, 11.08.2016
!     - included T-dependence of cyanobacteria growth
!     - modified stoichiometry for denitrification
! 
!     Purpose
!     -------
!     - declaration and memory allocation.
!
!     *kbo*         *INTEGER*  - number of wet cells in column.
!     *bolay*          *REAL*  - height of bottom cell.
!
!**********************************************************************
      implicit none

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: kbo
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: kwrbioz
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: k0100,k0500,k1000,k2000,k4000
      REAL, DIMENSION (:,:), ALLOCATABLE :: bolay

      REAL, DIMENSION (:,:), ALLOCATABLE :: expoor
      REAL, DIMENSION (:,:), ALLOCATABLE :: expoca
      REAL, DIMENSION (:,:), ALLOCATABLE :: exposi

      REAL, DIMENSION (:,:), ALLOCATABLE :: strahl

#ifdef FB_BGC_OCE     
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: abs_oce
#endif 

      REAL :: phytomi,grami,grazra,pi_alpha
      REAL :: remido,dyphy,zinges,epsher,spemor,gammap,gammaz,ecan
      REAL :: ro2ut,rcar,rnit,rnoi,rdnit0,rdnit1,rdnit2,rdn2o1,rdn2o2,rcalc,ropal
      REAL :: bluefix,tf2,tf1,tf0,tff  
      REAL :: bkphy,bkzoo,bkopal,bifr13,bifr14
      REAL :: wpoc,wcal,wopal,drempoc,dremdoc,dremn2o
      REAL :: dphymor,dzoomor,dremopal
      REAL :: dremsul
      REAL :: psedi,csedi,ssedi
      REAL :: perc_diron, riron, fesoly, relaxfe, fetune, wdust,bolaymin 
      REAL :: perc_disil
      REAL :: ctochl, atten_w, atten_c, atten_f
      REAL :: vol0
#ifdef AGG
      REAL :: SinkExp, FractDim, Stick, cellmass, cellsink, fsh, fse
      REAL :: alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
      REAL :: vsmall,safe,pupper,plower,zdis,nmldmin
      REAL :: dustd1,dustd2,dustd3,dustsink,calmax
#elif defined(WLIN)
      REAL :: wmin,wmax,wlin
#endif
#ifdef __c_isotopes
      REAL :: factor_13c, factor_14c, atm_c14_cal, atm_dc14_cal
      REAL :: PDB, ref14c, prei_d13C_atm, prei_dd14C_atm, atm_c13_cal
#endif

      CONTAINS

      SUBROUTINE ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

      use mod_xc
      use mo_control_bgc
      use mo_param1_bgc 

      INTEGER :: kpie,kpje,kpke
      INTEGER :: errstat
      
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable expoor ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (expoor(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory expoor'
         expoor(:,:) = 0.0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable expoca ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (expoca(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory expoca'
         expoca(:,:) = 0.0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable exposi ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (exposi(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory exposi'
         exposi(:,:) = 0.0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable kbo ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (kbo(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory kbo'
         kbo(:,:) = 0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable kwrbioz...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (kwrbioz(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory kwrbioz'
         kwrbioz(:,:) = 0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variables k0100, k0500, k1000, k2000 ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (k0100(kpie,kpje),stat=errstat)
         ALLOCATE (k0500(kpie,kpje),stat=errstat)
         ALLOCATE (k1000(kpie,kpje),stat=errstat)
         ALLOCATE (k2000(kpie,kpje),stat=errstat)
         ALLOCATE (k4000(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory k0100, k0500, k1000, k2000'
         k0100(:,:) = 0
         k0500(:,:) = 0
         k1000(:,:) = 0
         k2000(:,:) = 0
         k4000(:,:) = 0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable strahl ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (strahl(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory starhl'
         strahl(:,:) = 0.0

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable bolay ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (bolay(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory bolay'
         bolay(:,:) = 0.0


#ifdef FB_BGC_OCE 
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable abs_oce'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (abs_oce(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory abs_oce'
         abs_oce(:,:,:) = 0.0
#endif 

      END SUBROUTINE ALLOC_MEM_BIOMOD

      END MODULE mo_biomod
