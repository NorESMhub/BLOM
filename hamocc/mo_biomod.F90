! Copyright (C) 2001  S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, I. Kriest
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


      MODULE mo_biomod
!******************************************************************************
!
! MODULE mo_biomod - Variables for marine biology.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
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
!******************************************************************************
      implicit none

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
      REAL :: c14fac
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
!******************************************************************************
! ALLOC_MEM_BIOMOD - Allocate variables in this module
!******************************************************************************
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc

      INTEGER, intent(in) :: kpie,kpje,kpke
      INTEGER             :: errstat
      

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'***************************************************'
      WRITE(io_stdo_bgc,*)'Memory allocation for marine biology module :'
      WRITE(io_stdo_bgc,*)' '
      ENDIF


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


!******************************************************************************
     END SUBROUTINE ALLOC_MEM_BIOMOD

      END MODULE mo_biomod
