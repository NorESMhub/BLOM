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
      REAL, DIMENSION (:,:), ALLOCATABLE :: bolay

      REAL, DIMENSION (:,:), ALLOCATABLE :: expoor
      REAL, DIMENSION (:,:), ALLOCATABLE :: expoca
      REAL, DIMENSION (:,:), ALLOCATABLE :: exposi

      REAL, DIMENSION (:,:), ALLOCATABLE :: strahl

      REAL, DIMENSION (:,:), ALLOCATABLE :: alar1max
      REAL, DIMENSION (:,:), ALLOCATABLE :: TSFmax
      REAL, DIMENSION (:,:), ALLOCATABLE :: TMFmax

#ifdef FB_BGC_OCE     
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: abs_oce
#endif 

      REAL :: phytomi,grami,grazra,rrrcl,pi_alpha
      REAL :: remido,dyphy,zinges,epsher,spemor,gammap,gammaz,ecan
      REAL :: ro2ut,rcar,rnit,rnoi,rnit23,rnit13,rcalc,ropal,bluefix
      REAL :: bkphy,bkzoo,bkopal,bifr13,bifr14,plafr13,plafr14
      REAL :: wpoc,wcal,wopal,drempoc,dremdoc,dremcalc,dremn2o
      REAL :: dphymor,dzoomor,dremopal,calmax, gutc
      REAL :: dremsul
      REAL :: psedi,csedi,ssedi
      REAL :: perc_diron, riron, fesoly, relaxfe, wdust,bolaymin  
      REAL :: ctochl, atten_w, atten_c, atten_f
      REAL :: sco212_sfc,alkali_sfc,phosph_sfc,ano3_sfc,silica_sfc
      REAL :: vol0

#ifdef AGG      
      REAL :: SinkExp, FractDim, Stick, cellmass, cellsink, fsh, fse
      REAL :: alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
      REAL :: vsmall,safe,pupper,plower,zdis
      REAL :: dustd1,dustd2,dustd3,dustsink
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

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable expoca ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (expoca(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory expoca'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable exposi ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (exposi(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory exposi'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable kbo ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (kbo(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory kbo'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable kwrbioz...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (kwrbioz(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory kwrbioz'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable strahl ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (strahl(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory starhl'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable bolay ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (bolay(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory bolay'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable alar1max'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (alar1max(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory alar1max'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable TSFmax'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (TSFmax(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory TSFmax'

         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable TMFmax'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (TMFmax(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory TMFmax'

#ifdef FB_BGC_OCE 
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable abs_oce'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (abs_oce(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory abs_oce'
#endif 

      END SUBROUTINE ALLOC_MEM_BIOMOD

      END MODULE mo_biomod
