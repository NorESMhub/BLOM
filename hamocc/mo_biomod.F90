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
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved accumulation of all output fields to seperate subroutine,
!       new global fields for output defined here
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

      REAL, DIMENSION (:,:),   ALLOCATABLE :: bolay
      REAL, DIMENSION (:,:),   ALLOCATABLE :: strahl
#ifdef FB_BGC_OCE     
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: abs_oce
#endif 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: expoor
      REAL, DIMENSION (:,:),   ALLOCATABLE :: expoca
      REAL, DIMENSION (:,:),   ALLOCATABLE :: exposi
      REAL, DIMENSION (:,:),   ALLOCATABLE :: intphosy 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: intdnit 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: intnfix 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: intdmsprod 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: intdms_bac
      REAL, DIMENSION (:,:),   ALLOCATABLE :: intdms_uv 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: carflx0100 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: carflx0500 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: carflx1000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: carflx2000 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: carflx4000 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: carflx_bot 
      REAL, DIMENSION (:,:),   ALLOCATABLE :: bsiflx0100
      REAL, DIMENSION (:,:),   ALLOCATABLE :: bsiflx0500
      REAL, DIMENSION (:,:),   ALLOCATABLE :: bsiflx1000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: bsiflx2000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: bsiflx4000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: bsiflx_bot
      REAL, DIMENSION (:,:),   ALLOCATABLE :: calflx0100
      REAL, DIMENSION (:,:),   ALLOCATABLE :: calflx0500
      REAL, DIMENSION (:,:),   ALLOCATABLE :: calflx1000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: calflx2000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: calflx4000
      REAL, DIMENSION (:,:),   ALLOCATABLE :: calflx_bot
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: phosy3d
#ifdef AGG
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: wmass
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: wnumb
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: eps3d
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: asize3d
#endif

      REAL :: phytomi,grami,grazra,pi_alpha
      REAL :: remido,dyphy,zinges,epsher,spemor,gammap,gammaz,ecan
      REAL :: ro2ut,rcar,rnit,rnoi,rdnit0,rdnit1,rdnit2,rdn2o1,rdn2o2,rcalc,ropal
      REAL :: bluefix,tf2,tf1,tf0,tff  
      REAL :: bkphy,bkzoo,bkopal
      REAL :: wpoc,wcal,wopal
      REAL :: drempoc,dremopal,dremn2o,dremsul
      REAL :: perc_diron, riron, fesoly, relaxfe, fetune, wdust
      REAL :: ctochl, atten_w, atten_c, atten_f
#ifdef cisonew
      REAL :: c13fac,c14fac
      REAL :: re1312,re14to,prei13,prei14
      REAL :: bifr13,bifr14,growth_co2,bifr13_perm
#endif
#ifdef AGG
      REAL :: SinkExp, FractDim, Stick, cellmass, cellsink, fsh, fse
      REAL :: alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
      REAL :: vsmall,safe,pupper,plower,zdis,nmldmin
      REAL :: dustd1,dustd2,dustd3,dustsink,calmax
#elif defined(WLIN)
      REAL :: wmin,wmax,wlin
#endif
      CONTAINS

      SUBROUTINE ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

      use mod_xc
      use mo_control_bgc
      use mo_param1_bgc 

      INTEGER :: kpie,kpje,kpke
      INTEGER :: errstat
      
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
         WRITE(io_stdo_bgc,*)'Memory allocation for variable bolay ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (bolay(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory bolay'
         bolay(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable strahl ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (strahl(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory strahl'
         strahl(:,:) = 0.0


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
         WRITE(io_stdo_bgc,*)'Memory allocation for variable intphosy ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (intphosy(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory intphosy'
         intphosy(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable intdnit ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (intdnit(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory intdnit'
         intdnit(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable intnfix ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (intnfix(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory intnfix'
         intnfix(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable intdmsprod, intdms_bac, intdms_uv ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (intdmsprod(kpie,kpje),stat=errstat)
         ALLOCATE (intdms_bac(kpie,kpje),stat=errstat)
         ALLOCATE (intdms_uv(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory intdmsprod, intdms_bac, intdms_uv'
         intdmsprod(:,:) = 0.0
         intdms_bac(:,:) = 0.0
         intdms_uv(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable carflx* ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (carflx0100(kpie,kpje),stat=errstat)
         ALLOCATE (carflx0500(kpie,kpje),stat=errstat)
         ALLOCATE (carflx1000(kpie,kpje),stat=errstat)
         ALLOCATE (carflx2000(kpie,kpje),stat=errstat)
         ALLOCATE (carflx4000(kpie,kpje),stat=errstat)
         ALLOCATE (carflx_bot(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory carflx*'
         carflx0100(:,:) = 0.0
         carflx0500(:,:) = 0.0
         carflx1000(:,:) = 0.0
         carflx2000(:,:) = 0.0
         carflx4000(:,:) = 0.0
         carflx_bot(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable bsiflx* ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (bsiflx0100(kpie,kpje),stat=errstat)
         ALLOCATE (bsiflx0500(kpie,kpje),stat=errstat)
         ALLOCATE (bsiflx1000(kpie,kpje),stat=errstat)
         ALLOCATE (bsiflx2000(kpie,kpje),stat=errstat)
         ALLOCATE (bsiflx4000(kpie,kpje),stat=errstat)
         ALLOCATE (bsiflx_bot(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory bsiflx*'
         bsiflx0100(:,:) = 0.0
         bsiflx0500(:,:) = 0.0
         bsiflx1000(:,:) = 0.0
         bsiflx2000(:,:) = 0.0
         bsiflx4000(:,:) = 0.0
         bsiflx_bot(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable calflx* ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         ENDIF

         ALLOCATE (calflx0100(kpie,kpje),stat=errstat)
         ALLOCATE (calflx0500(kpie,kpje),stat=errstat)
         ALLOCATE (calflx1000(kpie,kpje),stat=errstat)
         ALLOCATE (calflx2000(kpie,kpje),stat=errstat)
         ALLOCATE (calflx4000(kpie,kpje),stat=errstat)
         ALLOCATE (calflx_bot(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory bsiflx*'
         calflx0100(:,:) = 0.0
         calflx0500(:,:) = 0.0
         calflx1000(:,:) = 0.0
         calflx2000(:,:) = 0.0
         calflx4000(:,:) = 0.0
         calflx_bot(:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable phosy3d ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (phosy3d(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory phosy3d'
         phosy3d(:,:,:) = 0.0
         

#ifdef AGG
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable wmass ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (wmass(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory eps3d'
         wmass(:,:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable wnumb ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (wnumb(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory eps3d'
         wnumb(:,:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable eps3d ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (eps3d(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory eps3d'
         eps3d(:,:,:) = 0.0


         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)'Memory allocation for variable asize3d ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
         ENDIF

         ALLOCATE (asize3d(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory asize3d'
         asize3d(:,:,:) = 0.0
#endif



      END SUBROUTINE ALLOC_MEM_BIOMOD

      END MODULE mo_biomod
