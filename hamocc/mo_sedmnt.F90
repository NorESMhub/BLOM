! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
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


       MODULE mo_sedmnt
!******************************************************************************
!
! MODULE mo_sedmnt - Variables for sediment modules.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - added sediment bypass preprocessor option
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
!******************************************************************************
      use mo_param1_bgc, only: ks,ksp,nsedtra,npowtra

      implicit none

      REAL, save :: dzs(ksp)    = 0.0
      REAL, save :: seddzi(ksp) = 0.0
      REAL, save :: seddw(ks)   = 0.0
      REAL, save :: porsol(ks)  = 0.0
      REAL, save :: porwah(ks)  = 0.0
      REAL, save :: porwat(ks)  = 0.0

      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: sedlay
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: powtra
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: sedhpl

      REAL, DIMENSION (:,:),     ALLOCATABLE :: silpro
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prorca
      REAL, DIMENSION (:,:),     ALLOCATABLE :: pror13
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prca13
      REAL, DIMENSION (:,:),     ALLOCATABLE :: pror14
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prca14
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prcaca
      REAL, DIMENSION (:,:),     ALLOCATABLE :: produs
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: burial

      REAL :: sedict,rno3,o2ut,ansed,sedac,sedifti
      REAL :: calcwei, opalwei, orgwei
      REAL :: calcdens, opaldens, orgdens, claydens
      REAL :: calfa, oplfa, orgfa, clafa, solfu

      CONTAINS


      SUBROUTINE ALLOC_MEM_SEDMNT(kpie,kpje)
!******************************************************************************
! ALLOC_MEM_SEDMNT - Allocate variables in this module
!******************************************************************************
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc
      
      INTEGER, intent(in) :: kpie,kpje
      INTEGER             :: errstat

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'***************************************************'
      WRITE(io_stdo_bgc,*)'Memory allocation for sediment module :'
      WRITE(io_stdo_bgc,*)' '
      ENDIF

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
#ifdef cisonew
      ALLOCATE (pror13(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pror13'
      pror13(:,:) = 0.0
      ALLOCATE (pror14(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pror14'
      pror14(:,:) = 0.0
#endif

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable prcaca ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (prcaca(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prcaca'
      prcaca(:,:) = 0.0
#ifdef cisonew
      ALLOCATE (prca13(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prca13'
      prca13(:,:) = 0.0
      ALLOCATE (prca14(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prca14'
      prca14(:,:) = 0.0
#endif

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable produs ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (produs(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory produs'
      produs(:,:) = 0.0


#ifndef sedbypass
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
#endif


!******************************************************************************
      END SUBROUTINE ALLOC_MEM_SEDMNT

      END MODULE mo_sedmnt
