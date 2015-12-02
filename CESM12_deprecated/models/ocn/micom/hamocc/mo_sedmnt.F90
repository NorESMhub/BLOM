       MODULE mo_sedmnt

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/mo_sedmnt.f90,v $\\
!$Revision: 1.2 $\\
!$Date: 2004/11/12 15:37:21 $\\

!***********************************************************************
!
!**** *MODULE mo_sedmnt* - Variables for sediment modules.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - declaration and memory allocation
!
!     *sedlay*         *REAL*  - .
!     *sedla1*         *REAL*  - .
!     *sedtot*         *REAL*  - .
!     *sedtoa*         *REAL*  - .
!     *seffel*         *REAL*  - .
!     *sedhpl*         *REAL*  - .
!     *powtra*         *REAL*  - .
!     *prorca*         *REAL*  - .
!     *prcaca*         *REAL*  - .
!     *silpro*         *REAL*  - .
!     *porwat*         *REAL*  - .
!     *porsol*         *REAL*  - .
!     *seddzi*         *REAL*  - .
!     *dzs*            *REAL*  - .
!     *porwah*         *REAL*  - .
!     *seddw*          *REAL*  - .
!     *sedict*         *REAL*  - .
!     *rno3*           *REAL*  - .
!     *calcon*         *REAL*  - .
!     *ansed*          *REAL*  - .
!     *o2ut*           *REAL*  - .
!
!**********************************************************************
      implicit none

      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: sedlay
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: powtra

      REAL, DIMENSION (:,:,:), ALLOCATABLE :: sedhpl

      REAL, DIMENSION (:), ALLOCATABLE :: seddw
      REAL, DIMENSION (:), ALLOCATABLE :: porsol
      REAL, DIMENSION (:), ALLOCATABLE :: porwah
      REAL, DIMENSION (:), ALLOCATABLE :: porwat

      REAL, DIMENSION (:), ALLOCATABLE :: dzs
      REAL, DIMENSION (:), ALLOCATABLE :: seddzi

      REAL, DIMENSION (:,:), ALLOCATABLE :: silpro
      REAL, DIMENSION (:,:), ALLOCATABLE :: prorca
      REAL, DIMENSION (:,:), ALLOCATABLE :: pror13
      REAL, DIMENSION (:,:), ALLOCATABLE :: prca13
      REAL, DIMENSION (:,:), ALLOCATABLE :: pror14
      REAL, DIMENSION (:,:), ALLOCATABLE :: prca14
      REAL, DIMENSION (:,:), ALLOCATABLE :: prcaca
      REAL, DIMENSION (:,:), ALLOCATABLE :: produs
      
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: burial


      REAL :: sedict,rno3,o2ut,ansed,sedac,sedifti
      REAL :: calcwei, opalwei, orgwei
      REAL :: calcdens, opaldens, orgdens, claydens
      REAL :: calfa, oplfa, orgfa, clafa, solfu

      CONTAINS

      SUBROUTINE ALLOC_MEM_SEDMNT(kpie,kpje)

      use mod_xc
      use mo_control_bgc
      use mo_param1_bgc 
      
      INTEGER :: kpie,kpje
      INTEGER :: errstat

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)' '
        WRITE(io_stdo_bgc,*)'******************************************'
        WRITE(io_stdo_bgc,*)' '
        WRITE(io_stdo_bgc,*)'Memory allocation for sediment modules :'
        ENDIF

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedlay ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nsedtra
        ENDIF

        ALLOCATE (sedlay(kpie,kpje,ks,nsedtra),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sedlay'
        sedlay(:,:,:,:) = 0.0


        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedhpl ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        ENDIF

        ALLOCATE (sedhpl(kpie,kpje,ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sedhpl'
        sedhpl(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable burial ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nsedtra
        ENDIF

        ALLOCATE (burial(kpie,kpje,nsedtra),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory burial'
        burial(:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable powtra ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',npowtra
        ENDIF

        ALLOCATE (powtra(kpie,kpje,ks,npowtra),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory powtra'
        powtra(:,:,:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable silpro ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (silpro(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory silpro'
        silpro(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable prorca ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (prorca(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory prorca'
        prorca(:,:) = 0.0
        ALLOCATE (pror13(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory pror13'
        pror13(:,:) = 0.0
        ALLOCATE (pror14(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory pror14'
        pror14(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable prcaca ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (prcaca(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory prcaca'
        prcaca(:,:) = 0.0
        ALLOCATE (prca13(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory prca13'
        prca13(:,:) = 0.0
        ALLOCATE (prca14(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory prca14'
        prca14(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable produs ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        ENDIF

        ALLOCATE (produs(kpie,kpje),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory produs'
        produs(:,:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable dzs ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ksp
        ENDIF

        ALLOCATE (dzs(ksp),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory dzs'
        dzs(:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable seddzi ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ksp
        ENDIF

        ALLOCATE (seddzi(ksp),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory seddzi'
        seddzi(:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable seddw ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks
        ENDIF

        ALLOCATE (seddw(ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory seddw'
        seddw(:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable porsol ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks
        ENDIF

        ALLOCATE (porsol(ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory porsol'
        porsol(:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable porwah ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks
        ENDIF

        ALLOCATE (porwah(ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory porwah'
        porwah(:) = 0.0

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable porwat ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks
        ENDIF

        ALLOCATE (porwat(ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory porwat'
        porwat(:) = 0.0



      END SUBROUTINE ALLOC_MEM_SEDMNT

      END MODULE mo_sedmnt
