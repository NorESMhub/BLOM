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
!     - initialization of sediment
!
!     Description:
!     ------------
!     Public routines and variable of this module:
!
!     -subroutine alloc_mem_sedmnt
!        Allocate memory for sediment variables
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
!     -subroutine ini_sedmnt
!         Initialize sediment parameters (some are also used in water column)
!     -subroutine ini_sedmnt_fields
!         Initialize 2D and 3D sediment fields
!
!******************************************************************************
 use mo_param1_bgc, only: ks,ksp,nsedtra,npowtra
 use mo_control_bgc, only: io_stdo_bgc
 use mod_xc,         only: mnproc

      implicit none

      REAL, save :: dzs(ksp)    = 0.0
      REAL, save :: seddzi(ksp) = 0.0
      REAL, save :: seddw(ks)   = 0.0

      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: sedlay
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: powtra
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: sedhpl
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: porsol
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: porwah
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: porwat
      REAL, DIMENSION (:,:),     ALLOCATABLE :: solfu
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: zcoefsu
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: zcoeflo

      REAL, DIMENSION (:,:),     ALLOCATABLE :: silpro
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prorca
      REAL, DIMENSION (:,:),     ALLOCATABLE :: pror13
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prca13
      REAL, DIMENSION (:,:),     ALLOCATABLE :: pror14
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prca14
      REAL, DIMENSION (:,:),     ALLOCATABLE :: prcaca
      REAL, DIMENSION (:,:),     ALLOCATABLE :: produs
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: burial

      REAL :: sedict,rno3,o2ut,ansed
      REAL :: calcwei, opalwei, orgwei
      REAL :: calcdens, opaldens, orgdens, claydens
      REAL :: calfa, oplfa, orgfa, clafa
      REAL :: disso_sil,silsat,disso_poc,sed_denit,disso_caco3

 CONTAINS
     
      !========================================================================
 SUBROUTINE ini_sedmnt(kpie,kpje,kpke,omask,sed_por)
      use mo_control_bgc, only: dtbgc

      implicit none

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)
      real,    intent(in) :: sed_por(kpie,kpje,ks)

      integer :: k

      sedict = 1.e-9 * dtbgc ! Molecular diffusion coefficient
      ! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]
      ! THIS NEEDS TO BE CHANGED TO disso=3.e-8! THIS IS ONLY KEPT FOR THE MOMENT
      ! FOR BACKWARDS COMPATIBILITY
      !disso_sil = 3.e-8*dtbgc  ! (2011-01-04) EMR
      !disso_sil = 1.e-6*dtbgc  ! test vom 03.03.04 half live sil ca. 20.000 yr 
      disso_sil = 1.e-6*dtbgc
      ! Silicate saturation concentration is 1 mol/m3
      silsat    = 0.001

      ! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]
      disso_poc = 0.01 / 86400. * dtbgc  !  disso=3.e-5 was quite high

      ! Denitrification rate constant of POP (disso) [1/sec]
      sed_denit =  0.01/86400. * dtbgc 

      ! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
      disso_caco3 = 1.e-7 * dtbgc

      ! ******************************************************************
      ! densities etc. for SEDIMENT SHIFTING

      ! define weight of calcium carbonate, opal, and poc [kg/kmol]
      calcwei = 100.           ! 40+12+3*16 kg/kmol C
      opalwei = 60.            ! 28 + 2*16  kg/kmol Si
      orgwei  = 30.            ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                           ! after Alldredge, 1998:
                           ! POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3

      ! define densities of opal, caco3, poc [kg/m3]
      calcdens = 2600.
      opaldens = 2200.
      orgdens  = 1000.
      claydens = 2600.         !quartz

      ! define volumes occupied by solid constituents [m3/kmol]
      calfa = calcwei / calcdens
      oplfa = opalwei / opaldens
      orgfa = orgwei / orgdens
      clafa = 1. / claydens    !clay is calculated in kg/m3

      ! sediment layer thickness
      dzs(1) = 0.001
      dzs(2) = 0.003
      dzs(3) = 0.005
      dzs(4) = 0.007
      dzs(5) = 0.009
      dzs(6) = 0.011
      dzs(7) = 0.013
      dzs(8) = 0.015
      dzs(9) = 0.017
      dzs(10) = 0.019
      dzs(11) = 0.021
      dzs(12) = 0.023
      dzs(13) = 0.025

      if (mnproc == 1) then
        write(io_stdo_bgc,*)  ' '
        write(io_stdo_bgc,*)  'Sediment layer thickness [m] : '
        write(io_stdo_bgc,'(5F9.3)') dzs
        write(io_stdo_bgc,*)  ' '
      endif

      seddzi(1) = 500.
      do k = 1, ks
        seddzi(k+1) = 1. / dzs(k+1)          ! inverse of grid cell size
        seddw(k) = 0.5 * (dzs(k) + dzs(k+1)) ! distance between grid cell centers (pressure points)
      enddo

#ifndef sedbypass
      ! 2d and 3d fields are not allocated in case of sedbypass
      ! so only initialize them if we are using the sediment
      CALL ini_sedmnt_por(kpie,kpje,kpke,omask,sed_por)
#endif     
 END SUBROUTINE ini_sedmnt

      !========================================================================
 SUBROUTINE ini_sedmnt_por(kpie,kpje,kpke,omask,sed_por)
      !
      ! Initialization of:
      !   - 3D porosity field (cell center and cell boundaries)
      !   - solid volume fraction at cell center
      !   - vertical molecular diffusion coefficients scaled with porosity
      !
      use mo_control_bgc, only: l_3Dvarsedpor

      implicit none

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)
      real,    intent(in) :: sed_por(kpie,kpje,ks)

      ! local
      integer :: i,j,k

      ! this initialization can be done via reading a porosity map
      ! porwat is the poroisty at the (pressure point) center of the grid cell
      if (l_3Dvarsedpor)then
       ! lon-lat variable sediment porosity from input file
       do k=1,ks
       do j=1,kpje
       do i=1,kpie
        if(omask(i,j).gt. 0.5)then
          porwat(i,j,k) = sed_por(i,j,k)
        endif
       enddo
       enddo
       enddo
      else
        porwat(:,:,1) = 0.85
        porwat(:,:,2) = 0.83
        porwat(:,:,3) = 0.8
        porwat(:,:,4) = 0.79
        porwat(:,:,5) = 0.77
        porwat(:,:,6) = 0.75
        porwat(:,:,7) = 0.73
        porwat(:,:,8) = 0.7
        porwat(:,:,9) = 0.68
        porwat(:,:,10) = 0.66
        porwat(:,:,11) = 0.64
        porwat(:,:,12) = 0.62
      endif

      if (mnproc == 1) then
        write(io_stdo_bgc,*)  'Pore water in sediment initialized'
      endif

      do k = 1, ks
      do j = 1, kpje
      do i = 1, kpie
        porsol(i,j,k) = 1. - porwat(i,j,k)                                  ! solid volume fraction at grid center
        if(k >= 2) porwah(i,j,k) = 0.5 * (porwat(i,j,k) + porwat(i,j,k-1))  ! porosity at cell interfaces
        if(k == 1) porwah(i,j,k) = 0.5 * (1. + porwat(i,j,1))
      enddo
      enddo
      enddo

      ! determine total solid sediment volume
      solfu = 0.
      do i = 1, kpie
      do j = 1, kpje
      do k = 1, ks
        solfu(i,j) = solfu(i,j) + seddw(k) * porsol(i,j,k)
      enddo
      enddo
      enddo

      ! Initialize porosity-dependent diffusion coefficients of sediment
      zcoefsu(:,:,0) = 0.0
      do k = 1,ks
      do j = 1, kpje
      do i = 1, kpie
        ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
        zcoefsu(i,j,k  ) = -sedict * seddzi(k) * porwah(i,j,k)
        zcoeflo(i,j,k-1) = -sedict * seddzi(k) * porwah(i,j,k)    ! why the same ?
      enddo
      enddo
      enddo
      zcoeflo(:,:,ks) = 0.0                    ! diffusion coefficient for bottom sediment layer
      
      if (mnproc == 1) then
        write(io_stdo_bgc,*)  'Pore water diffusion coefficients in sediment initialized'
      endif

 END SUBROUTINE ini_sedmnt_por


      !========================================================================
 SUBROUTINE ALLOC_MEM_SEDMNT(kpie,kpje)
 !******************************************************************************
 ! ALLOC_MEM_SEDMNT - Allocate variables in this module
 !******************************************************************************
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
      WRITE(io_stdo_bgc,*)'Memory allocation for variable porsol ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
      ENDIF

      ALLOCATE (porsol(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porsol'
      porsol(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable porwah ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
      ENDIF

      ALLOCATE (porwah(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porwah'
      porwah(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable porwat ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
      ENDIF

      ALLOCATE (porwat(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porwat'
      porwat(:,:,:) = 0.0
      
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable solfu ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (solfu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory solfu'
      solfu(:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable zcoefsu ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
      ENDIF

      ALLOCATE (zcoefsu(kpie,kpje,0:ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory zcoefsu'
      zcoefsu(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable zcoeflo ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
      ENDIF

      ALLOCATE (zcoeflo(kpie,kpje,0:ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory zcoeflo'
      zcoeflo(:,:,:) = 0.0

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
