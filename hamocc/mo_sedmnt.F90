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


module mo_sedmnt

  !******************************************************************************
  ! Variables for sediment modules.
  ! - declaration and memory allocation
  ! - initialization of sediment
  !
  ! S.Legutke,        *MPI-MaD, HH*    31.10.01
  ! Modified:
  ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  ! - added sediment bypass preprocessor option
  !******************************************************************************

  use mod_xc,         only: mnproc
  use mo_param1_bgc,  only: ks,ksp,nsedtra,npowtra
  use mo_control_bgc, only: io_stdo_bgc,use_sedbypass,use_cisonew

  implicit none
  private

  ! Routines
  public  :: ini_sedmnt        ! Initialize sediment parameters and sediment vertical grid
  public  :: alloc_mem_sedmnt  ! Allocate memory for sediment variables
  private :: ini_sedmnt_por    ! Initialize 2D and 3D sediment fields

  ! Module variables
  real, protected, public :: dzs(ksp)    = 0.0
  real, protected, public :: seddzi(ksp) = 0.0
  real, protected, public :: seddw(ks)   = 0.0

  real, dimension (:,:,:,:), allocatable, public :: sedlay
  real, dimension (:,:,:,:), allocatable, public :: powtra
  real, dimension (:,:,:),   allocatable, public :: sedhpl
  real, dimension (:,:,:),   allocatable, public :: porsol
  real, dimension (:,:,:),   allocatable, public :: porwah
  real, dimension (:,:,:),   allocatable, public :: porwat
  real, dimension (:,:),     allocatable, public :: solfu
  real, dimension (:,:,:),   allocatable, public :: zcoefsu
  real, dimension (:,:,:),   allocatable, public :: zcoeflo

  real, dimension (:,:),     allocatable, public :: silpro
  real, dimension (:,:),     allocatable, public :: prorca
  real, dimension (:,:),     allocatable, public :: pror13
  real, dimension (:,:),     allocatable, public :: prca13
  real, dimension (:,:),     allocatable, public :: pror14
  real, dimension (:,:),     allocatable, public :: prca14
  real, dimension (:,:),     allocatable, public :: prcaca
  real, dimension (:,:),     allocatable, public :: produs
  real, dimension (:,:,:),   allocatable, public :: burial

  real, protected, public :: calfa, oplfa, orgfa, clafa

CONTAINS

  !******************************************************************************
  subroutine ini_sedmnt(kpie,kpje,kpke,omask,sed_por)

    use mo_param_bgc, only: claydens,calcwei,calcdens,opalwei,opaldens,orgwei,orgdens,sedict

    ! Arguments
    integer, intent(in) :: kpie,kpje,kpke
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: sed_por(kpie,kpje,ks)

    ! Local variables
    integer :: k

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

    if (.not. use_sedbypass) then
      ! 2d and 3d fields are not allocated in case of sedbypass
      ! so only initialize them if we are using the sediment
      call ini_sedmnt_por(kpie,kpje,kpke,omask,sed_por)
    endif

  end subroutine ini_sedmnt

  !******************************************************************************
  subroutine ini_sedmnt_por(kpie,kpje,kpke,omask,sed_por)
    !
    ! Initialization of:
    ! - 3D porosity field (cell center and cell boundaries)
    ! - solid volume fraction at cell center
    ! - vertical molecular diffusion coefficients scaled with porosity
    !
    use mo_control_bgc, only: l_3Dvarsedpor
    use mo_param_bgc,   only: sedict

    ! Arguments
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
            if(omask(i,j).gt. 0.5) then
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

  end subroutine ini_sedmnt_por

  !******************************************************************************
  subroutine alloc_mem_sedmnt(kpie,kpje)

    ! ------------------------------------------------------
    !  Allocate variables in this module
    ! ------------------------------------------------------

    ! Arguments
    integer, intent(in) :: kpie,kpje

    ! Local variables
    integer :: errstat

    if (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'Memory allocation for sediment module :'
      write(io_stdo_bgc,*)' '
    endif

    if (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)'Memory allocation for variable silpro ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (silpro(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory silpro'
    silpro(:,:) = 0.0

    if (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)'Memory allocation for variable prorca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (prorca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory prorca'
    prorca(:,:) = 0.0
    if (use_cisonew) then
      allocate (pror13(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pror13'
      pror13(:,:) = 0.0
      allocate (pror14(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pror14'
      pror14(:,:) = 0.0
    endif

    if (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)'Memory allocation for variable prcaca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (prcaca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory prcaca'
    prcaca(:,:) = 0.0
    if (use_cisonew) then
      allocate (prca13(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prca13'
      prca13(:,:) = 0.0
      allocate (prca14(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prca14'
      prca14(:,:) = 0.0
    endif

    if (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)'Memory allocation for variable produs ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (produs(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory produs'
    produs(:,:) = 0.0

    if (.not. use_sedbypass) then
      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable sedlay ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
        write(io_stdo_bgc,*)'Forth dimension    : ',nsedtra
      endif
      allocate (sedlay(kpie,kpje,ks,nsedtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedlay'
      sedlay(:,:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable sedhpl ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (sedhpl(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedhpl'
      sedhpl(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable porsol ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (porsol(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porsol'
      porsol(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable porwah ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (porwah(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porwah'
      porwah(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable porwat ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (porwat(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porwat'
      porwat(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable solfu ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (solfu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory solfu'
      solfu(:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable zcoefsu ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (zcoefsu(kpie,kpje,0:ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory zcoefsu'
      zcoefsu(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable zcoeflo ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (zcoeflo(kpie,kpje,0:ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory zcoeflo'
      zcoeflo(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable burial ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',nsedtra
      endif
      allocate (burial(kpie,kpje,nsedtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory burial'
      burial(:,:,:) = 0.0

      if (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable powtra ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
        write(io_stdo_bgc,*)'Forth dimension    : ',npowtra
      endif
      allocate (powtra(kpie,kpje,ks,npowtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory powtra'
      powtra(:,:,:,:) = 0.0
    endif

  end subroutine alloc_mem_sedmnt

end module mo_sedmnt
