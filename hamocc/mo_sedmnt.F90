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

  !*************************************************************************************************
  ! Sediment variables, declaration and memory allocation, and initialization
  !
  ! S.Legutke,        *MPI-MaD, HH*    31.10.01
  !
  ! Modified:
  ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  !  - added sediment bypass preprocessor option
  !*************************************************************************************************

  use mod_xc,         only: mnproc
  use mo_kind,        only: rp
  use mo_param1_bgc,  only: ks,ksp,nsedtra,npowtra
  use mo_control_bgc, only: io_stdo_bgc,use_sedbypass,use_cisonew,use_sediment_quality

  implicit none
  private

  ! Routines
  public  :: ini_sedmnt        ! Initialize sediment parameters and sediment vertical grid
  public  :: alloc_mem_sedmnt  ! Allocate memory for sediment variables
  private :: ini_sedmnt_por    ! Initialize 2D and 3D sediment fields

  ! Module variables
  real(rp), protected, public :: dzs(ksp)    = 0.0_rp
  real(rp), protected, public :: seddzi(ksp) = 0.0_rp
  real(rp), protected, public :: seddw(ks)   = 0.0_rp

  real(rp), dimension (:,:,:,:), allocatable, public :: sedlay
  real(rp), dimension (:,:,:,:), allocatable, public :: powtra
  real(rp), dimension (:,:,:),   allocatable, public :: sedhpl
  real(rp), dimension (:,:,:),   allocatable, public :: porsol
  real(rp), dimension (:,:,:),   allocatable, public :: porwah
  real(rp), dimension (:,:,:),   allocatable, public :: porwat
  real(rp), dimension (:,:),     allocatable, public :: solfu
  real(rp), dimension (:,:,:),   allocatable, public :: zcoefsu
  real(rp), dimension (:,:,:),   allocatable, public :: zcoeflo

  real(rp), dimension (:,:),     allocatable, public :: silpro
  real(rp), dimension (:,:),     allocatable, public :: prorca
  real(rp), dimension (:,:),     allocatable, public :: prorca_mavg
  real(rp), dimension (:,:),     allocatable, public :: pror13
  real(rp), dimension (:,:),     allocatable, public :: prca13
  real(rp), dimension (:,:),     allocatable, public :: pror14
  real(rp), dimension (:,:),     allocatable, public :: prca14
  real(rp), dimension (:,:),     allocatable, public :: prcaca
  real(rp), dimension (:,:),     allocatable, public :: produs
  real(rp), dimension (:,:,:),   allocatable, public :: burial

  ! Output diagnostics
  real(rp), dimension (:,:,:),   allocatable, public :: sed_rem_aerob
  real(rp), dimension (:,:,:),   allocatable, public :: sed_rem_denit
  real(rp), dimension (:,:,:),   allocatable, public :: sed_rem_sulf

  ! values for sediment quality-driven remineralization
  real(rp), dimension(:,:,:),    allocatable, public :: sed_reactivity_a
  real(rp), dimension(:,:,:),    allocatable, public :: sed_reactivity_k
  real(rp), dimension(:,:,:),    allocatable, public :: sed_applied_reminrate

  real(rp), protected, public :: calfa, oplfa, orgfa, clafa

CONTAINS

  subroutine ini_sedmnt(kpie,kpje,omask,sed_por,sed_POCage_init,prorca_mavg_init)
    !***********************************************************************************************

    use mo_param_bgc, only: claydens,calcwei,calcdens,opalwei,opaldens,orgwei,orgdens,sedict

    ! Arguments
    integer, intent(in) :: kpie,kpje
    real(rp),intent(in) :: omask(kpie,kpje)
    real(rp),intent(in) :: sed_por(kpie,kpje,ks)
    real(rp),intent(in) :: sed_POCage_init(kpie,kpje,ks)
    real(rp),intent(in) :: prorca_mavg_init(kpie,kpje)

    ! Local variables
    integer :: k

    ! define volumes occupied by solid constituents [m3/kmol]
    calfa = calcwei / calcdens
    oplfa = opalwei / opaldens
    orgfa = orgwei / orgdens
    clafa = 1._rp / claydens    !clay is calculated in kg/m3

    ! sediment layer thickness
    dzs(1) = 0.001_rp
    dzs(2) = 0.003_rp
    dzs(3) = 0.005_rp
    dzs(4) = 0.007_rp
    dzs(5) = 0.009_rp
    dzs(6) = 0.011_rp
    dzs(7) = 0.013_rp
    dzs(8) = 0.015_rp
    dzs(9) = 0.017_rp
    dzs(10) = 0.019_rp
    dzs(11) = 0.021_rp
    dzs(12) = 0.023_rp
    dzs(13) = 0.025_rp

    if (mnproc == 1) then
      write(io_stdo_bgc,*)  ' '
      write(io_stdo_bgc,*)  'Sediment layer thickness [m] : '
      write(io_stdo_bgc,'(5F9.3)') dzs
      write(io_stdo_bgc,*)  ' '
    endif

    seddzi(1) = 500._rp
    do k = 1, ks
      seddzi(k+1) = 1._rp / dzs(k+1)          ! inverse of grid cell size
      seddw(k) = 0.5_rp * (dzs(k) + dzs(k+1)) ! distance between grid cell centers (pressure points)
    enddo

    if (.not. use_sedbypass) then
      ! 2d and 3d fields are not allocated in case of sedbypass
      ! so only initialize them if we are using the sediment
      call ini_sedmnt_por(kpie,kpje,omask,sed_por)
      if (use_sediment_quality) then
        call ini_sed_qual(kpie,kpje,omask,sed_POCage_init,prorca_mavg_init)
      endif
    endif


  end subroutine ini_sedmnt

  subroutine ini_sedmnt_por(kpie,kpje,omask,sed_por)
    !***********************************************************************************************
    ! Initialization of:
    ! - 3D porosity field (cell center and cell boundaries)
    ! - solid volume fraction at cell center
    ! - vertical molecular diffusion coefficients scaled with porosity
    !***********************************************************************************************
    use mo_control_bgc, only: l_3Dvarsedpor
    use mo_param_bgc,   only: sedict

    ! Arguments
    integer, intent(in) :: kpie,kpje
    real(rp),intent(in) :: omask(kpie,kpje)
    real(rp),intent(in) :: sed_por(kpie,kpje,ks)

    ! local
    integer :: i,j,k

    ! this initialization can be done via reading a porosity map
    ! porwat is the porosity at the (pressure point) center of the grid cell
    if (l_3Dvarsedpor)then
      ! lon-lat variable sediment porosity from input file
      do k=1,ks
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j).gt. 0.5_rp) then
              porwat(i,j,k) = sed_por(i,j,k)
            endif
          enddo
        enddo
      enddo
    else
      porwat(:,:,1) = 0.85_rp
      porwat(:,:,2) = 0.83_rp
      porwat(:,:,3) = 0.8_rp
      porwat(:,:,4) = 0.79_rp
      porwat(:,:,5) = 0.77_rp
      porwat(:,:,6) = 0.75_rp
      porwat(:,:,7) = 0.73_rp
      porwat(:,:,8) = 0.7_rp
      porwat(:,:,9) = 0.68_rp
      porwat(:,:,10) = 0.66_rp
      porwat(:,:,11) = 0.64_rp
      porwat(:,:,12) = 0.62_rp
    endif

    if (mnproc == 1) then
      write(io_stdo_bgc,*)  'Pore water in sediment initialized'
    endif

    do k = 1, ks
      do j = 1, kpje
        do i = 1, kpie
          porsol(i,j,k) = 1._rp - porwat(i,j,k)                                  ! solid volume fraction at grid center
          if(k >= 2) porwah(i,j,k) = 0.5_rp * (porwat(i,j,k) + porwat(i,j,k-1))  ! porosity at cell interfaces
          if(k == 1) porwah(i,j,k) = 0.5_rp * (1._rp + porwat(i,j,1))
        enddo
      enddo
    enddo

    ! determine total solid sediment volume
    solfu = 0._rp
    do i = 1, kpie
      do j = 1, kpje
        do k = 1, ks
          solfu(i,j) = solfu(i,j) + seddw(k) * porsol(i,j,k)
        enddo
      enddo
    enddo

    ! Initialize porosity-dependent diffusion coefficients of sediment
    zcoefsu(:,:,0) = 0.0_rp
    do k = 1,ks
      do j = 1, kpje
        do i = 1, kpie
          ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
          zcoefsu(i,j,k  ) = -sedict * seddzi(k) * porwah(i,j,k)
          zcoeflo(i,j,k-1) = -sedict * seddzi(k) * porwah(i,j,k)    ! why the same ?
        enddo
      enddo
    enddo
    zcoeflo(:,:,ks) = 0.0_rp                    ! diffusion coefficient for bottom sediment layer

    if (mnproc == 1) then
      write(io_stdo_bgc,*)  'Pore water diffusion coefficients in sediment initialized'
    endif

  end subroutine ini_sedmnt_por

  subroutine ini_sed_qual(kpie,kpje,omask,sed_POCage_init,prorca_mavg_init)
    !-----------------------------------------------------------------------------------------------
    ! Initialize moving average prorca and sediment POC age
    ! use burial age equiv. to oldest sed layer
    !-----------------------------------------------------------------------------------------------
    use mo_param1_bgc,   only: issso12_age
    use mo_param_bgc,    only: sec_per_day

    implicit none
    ! Arguments
    integer, intent(in) :: kpie,kpje
    real(rp),intent(in) :: omask(kpie,kpje)
    real(rp),intent(in) :: sed_POCage_init(kpie,kpje,ks)
    real(rp),intent(in) :: prorca_mavg_init(kpie,kpje)

    ! Local variables
    integer :: i,j,k

    if (mnproc == 1) then
      write(io_stdo_bgc,*)  ' '
      write(io_stdo_bgc,*)  'Initializing sediment quality: age and moving average prorca'
      write(io_stdo_bgc,*)  ' '
    endif

    do i = 1,kpie
      do j = 1,kpje
          ! Units: prorca_mavg_init expected to be in [kmol P m-2 s-1]
          !        - needs to be converted to [mmol P m-2 d-1]
          prorca_mavg(i,j)        = prorca_mavg_init(i,j)*1.0e6_rp/sec_per_day
          burial(i,j,issso12_age) = sed_POCage_init(i,j,ks)
        do k = 1,ks
          sedlay(i,j,k,issso12_age) = sed_POCage_init(i,j,k)
        enddo
      enddo
    enddo

  end subroutine ini_sed_qual

  subroutine alloc_mem_sedmnt(kpie,kpje)
    !***********************************************************************************************
    !  Allocate variables in this module
    !***********************************************************************************************
    use mo_control_bgc, only: use_extNcycle

    ! Arguments
    integer, intent(in) :: kpie,kpje

    ! Local variables
    integer :: errstat

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'Memory allocation for sediment module :'
      write(io_stdo_bgc,*)' '
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable silpro ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (silpro(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory silpro'
    silpro(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable prorca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (prorca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory prorca'
    prorca(:,:) = 0.0_rp
    if (use_cisonew) then
      allocate (pror13(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pror13'
      pror13(:,:) = 0.0_rp
      allocate (pror14(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pror14'
      pror14(:,:) = 0.0_rp
    endif

    if (use_sediment_quality) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable prorca_mavg ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (prorca_mavg(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prorca_mavg'
      prorca_mavg(:,:) = 0.0_rp


      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable sed_reactivity_a ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (sed_reactivity_a(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sed_reactivity_a'
      sed_reactivity_a(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable sed_reactivity_k ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (sed_reactivity_k(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sed_reactivity_k'
      sed_reactivity_k(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable sed_applied_reminrate ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (sed_applied_reminrate(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sed_applied_reminrate'
      sed_applied_reminrate(:,:,:) = 0.0_rp
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable prcaca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (prcaca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory prcaca'
    prcaca(:,:) = 0.0_rp
    if (use_cisonew) then
      allocate (prca13(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prca13'
      prca13(:,:) = 0.0_rp
      allocate (prca14(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory prca14'
      prca14(:,:) = 0.0_rp
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable produs ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (produs(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory produs'
    produs(:,:) = 0.0_rp

    if (.not. use_sedbypass) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable sedlay ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
        write(io_stdo_bgc,*)'Forth dimension    : ',nsedtra
      endif
      allocate (sedlay(kpie,kpje,ks,nsedtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedlay'
      sedlay(:,:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable sedhpl ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (sedhpl(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedhpl'
      sedhpl(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable porsol ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (porsol(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porsol'
      porsol(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable porwah ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (porwah(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porwah'
      porwah(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable porwat ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (porwat(kpie,kpje,ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory porwat'
      porwat(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable solfu ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (solfu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory solfu'
      solfu(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable zcoefsu ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (zcoefsu(kpie,kpje,0:ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory zcoefsu'
      zcoefsu(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable zcoeflo ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
      endif
      allocate (zcoeflo(kpie,kpje,0:ks),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory zcoeflo'
      zcoeflo(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable burial ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',nsedtra
      endif
      allocate (burial(kpie,kpje,nsedtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory burial'
      burial(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable powtra ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
        write(io_stdo_bgc,*)'Forth dimension    : ',npowtra
      endif
      allocate (powtra(kpie,kpje,ks,npowtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory powtra'
      powtra(:,:,:,:) = 0.0_rp

      if (.not. use_extNcycle) then
        if (mnproc.eq.1) then
          write(io_stdo_bgc,*)'Memory allocation for variable sed_rem_aerob ..'
          write(io_stdo_bgc,*)'First dimension    : ',kpie
          write(io_stdo_bgc,*)'Second dimension   : ',kpje
          write(io_stdo_bgc,*)'Third dimension    : ',ks
        endif
        allocate (sed_rem_aerob(kpie,kpje,ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sed_rem_aerob'
        sed_rem_aerob(:,:,:) = 0.0_rp

        if (mnproc.eq.1) then
          write(io_stdo_bgc,*)'Memory allocation for variable sed_rem_denit ..'
          write(io_stdo_bgc,*)'First dimension    : ',kpie
          write(io_stdo_bgc,*)'Second dimension   : ',kpje
          write(io_stdo_bgc,*)'Third dimension    : ',ks
        endif
        allocate (sed_rem_denit(kpie,kpje,ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sed_rem_denit'
        sed_rem_denit(:,:,:) = 0.0_rp

        if (mnproc.eq.1) then
          write(io_stdo_bgc,*)'Memory allocation for variable sed_rem_sulf ..'
          write(io_stdo_bgc,*)'First dimension    : ',kpie
          write(io_stdo_bgc,*)'Second dimension   : ',kpje
          write(io_stdo_bgc,*)'Third dimension    : ',ks
        endif
        allocate (sed_rem_sulf(kpie,kpje,ks),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory sed_rem_sulf'
        sed_rem_sulf(:,:,:) = 0.0_rp
      endif
    endif ! use_sedbypass

  end subroutine alloc_mem_sedmnt

end module mo_sedmnt
