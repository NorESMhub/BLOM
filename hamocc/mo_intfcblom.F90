! Copyright (C) 2020  J. Schwinger
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


module mo_intfcblom

  !******************************************************************************
  ! Variables for BLOM-iHAMOCC interface
  ! - Declaration and memory allocation related to the BLOM-iHAMOCC interface.
  ! - This includes 2-time-level copies of sediment and amospheric fields.
  !
  !  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
  !******************************************************************************

  use mo_control_bgc, only: use_sedbypass,use_BOXATM

  implicit none
  private

  ! Routines

  public :: alloc_mem_intfcblom ! Allocate memory for BLOM interface variables
  public :: blom2hamocc         ! Transfer fields from BLOM to HAMOCC
  public :: hamocc2blom         ! Transfer fields from HAMOCC to BLOM

  ! Module variables

  integer, parameter, public :: nphys=2       ! number of bgc timesteps per ocean timestep.

  real, allocatable, public  :: bgc_dx(:,:)    ! size of grid cell (longitudinal) [m].
  real, allocatable, public  :: bgc_dy(:,:)    ! size of grid cell (latitudinal) [m].
  real, allocatable, public  :: bgc_dp(:,:,:)  ! size of grid cell (depth) [m].
  real, allocatable, public  :: bgc_rho(:,:,:) ! sea water density [kg/m^3].

  real, allocatable, public :: omask(:,:)     ! land ocean mask.

  ! The following arrays are used to keep a two time-level copy of sediment
  ! and prognostic atmosphere fields. These arrays are copied back and forth
  ! in blom2hamocc.F and hamocc2blom.F in the same manner as the tracer field.
  ! Also, they written/read to and from restart files:

  ! Two time-level copy of sediment fields
  real, allocatable, public :: sedlay2(:,:,:,:)
  real, allocatable, public :: powtra2(:,:,:,:)
  real, allocatable, public :: burial2(:,:,:,:)

  ! Two time level copy of prognostic atmosphere field used if BOXATM is activated
  real, allocatable, public :: atm2(:,:,:,:)

contains

  subroutine alloc_mem_intfcblom(kpie,kpje,kpke)

    !******************************************************************************
    ! Allocate variables in this module
    !
    ! J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
    !******************************************************************************

    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc
    use mo_param1_bgc,  only: ks,nsedtra,npowtra,natm

    ! Arguments
    integer, intent(in) :: kpie,kpje,kpke

    ! Local variables
    integer :: errstat

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'Memory allocation for module mo_intfcblom :'
      write(io_stdo_bgc,*)' '
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable bgc_dx, bgc_dy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (bgc_dx(kpie,kpje),stat=errstat)
    allocate (bgc_dy(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory bgc_dx, bgc_dy'
    bgc_dx(:,:) = 0.0
    bgc_dy(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable bgc_dp ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (bgc_dp(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory bgc_dp'
    bgc_dp(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable bgc_rho ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (bgc_rho(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory bgc_dp'
    bgc_rho(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable omask ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate(omask(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory omask'
    omask(:,:) = 0.0

    if (.not. use_sedbypass) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable sedlay2 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',2*ks
        write(io_stdo_bgc,*)'Fourth dimension   : ',nsedtra
      endif
      allocate (sedlay2(kpie,kpje,2*ks,nsedtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedlay2'
      sedlay2(:,:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable powtra2 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',2*ks
        write(io_stdo_bgc,*)'Fourth dimension   : ',npowtra
      endif
      allocate (powtra2(kpie,kpje,2*ks,npowtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory powtra2'
      powtra2(:,:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable burial2 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',2
        write(io_stdo_bgc,*)'Fourth dimension   : ',nsedtra
      endif
      allocate (burial2(kpie,kpje,2,nsedtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory burial2'
      burial2(:,:,:,:) = 0.0
    endif

    if (use_BOXATM) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable atm2 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',2
        write(io_stdo_bgc,*)'Fourth dimension   : ',natm
      endif
      allocate (atm2(kpie,kpje,2,natm),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory atm2'
      atm2(:,:,:,:) = 0.0
    endif

  end subroutine alloc_mem_intfcblom


  subroutine blom2hamocc(m,n,mm,nn)

    !******************************************************************************
    ! Interface between BLOM and HAMOCC.
    !
    ! K. Assmann        *GFI, UiB        initial version
    ! J. Schwinger      *GFI, UiB        2013-04-22
    ! Modified:
    ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    ! - removed inverse of layer thickness
    ! - added sediment bypass preprocessor option
    ! M. Bentsen,       *NORCE, Bergen*          2020-05-03
    ! - changed ocean model from MICOM to BLOM
    ! T. Torsvik,       *University of Bergen*   2021-08-26
    ! - integrate subroutine into module mo_intfcblom
    !******************************************************************************

    use mod_constants, only: onem
    use mod_xc,        only: ii,jdm,jj,kdm,kk,ifp,isp,ilp,idm
    use mod_grid,      only: scpx,scpy
    use mod_state,     only: dp,temp,saln
    use mod_eos,       only: rho,p_alpha
    use mod_difest,    only: hOBL
    use mod_tracers,   only: ntrbgc,itrbgc,trc
    use mo_param1_bgc, only: ks,nsedtra,npowtra,natm
    use mo_carbch,     only: ocetra,atm
    use mo_sedmnt,     only: sedlay,powtra,sedhpl,burial
    use mo_vgrid,      only: kmle, kmle_static

    ! Arguments
    integer, intent(in) :: m,n,mm,nn

    ! Local variables
    integer :: i,j,k,l,nns,kn
    real    :: p1,p2,ldp,th,s,pa
    real    :: rp(idm,jdm,kdm+1)

    nns=(n-1)*ks

    rp(:,:,:) = 0.0

    ! --- calculate pressure at interfaces (necesarry since p has
    ! --- not been calculated at restart)

    !$OMP PARALLEL DO PRIVATE(k,kn,l,i)
    do k=1,kk
      kn=k+nn
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            rp(i,j,k+1) = rp(i,j,k) + dp(i,j,kn)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- 2D fields
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,jj
      do i=1,ii
        !
        ! --- - dimension of grid box in meters
        bgc_dx(i,j) = scpx(i,j)/1.e2
        bgc_dy(i,j) = scpy(i,j)/1.e2
        !
        ! --- - index of level above OBL depth
        ! ---   isopycninc coords: hOBL(i,j) = hOBL_static = 3.  =>  kmle(i,j) = 2
        ! ---   hybrid coords: hOBL defined according to cvmix_kpp_compute_kOBL_depth
        kmle(i,j) = nint(hOBL(i,j))-1
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- 3D fields
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(k,kn,l,i,th,s,p1,p2,ldp,pa)
    do k=1,kk
      kn=k+nn
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            !
            ! --- - integrated specific volume
            th = temp(i,j,kn)
            s  = saln(i,j,kn)
            p1 = rp(i,j,k)
            if(dp(i,j,kn) == 0.0) then
              ldp = 1.0
              pa  = ldp/rho(p1,th,s)
            else if(dp(i,j,kn) <  1.0e-2) then
              ldp = dp(i,j,kn)
              pa  = ldp/rho(p1,th,s)
            else
              ldp = dp(i,j,kn)
              p2  = p1+ldp
              pa  = p_alpha(p1,p2,th,s)
            endif
            !
            ! --- - density in g/cm^3
            bgc_rho(i,j,k)=ldp/pa
            !
            ! --- - layer thickness in meters
            bgc_dp(i,j,k) = 0.0
            if(dp(i,j,kn).ne.0.0) bgc_dp(i,j,k) = pa / onem
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- - return if restart (HAMOCC fields are not allocated yet)
    ! --- ------------------------------------------------------------------
    if( .not. allocated(ocetra) ) return

    ! --- ------------------------------------------------------------------
    ! --- pass tracer fields from ocean model; convert mol/kg -> kmol/m^3
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(k,kn,l,i)
    do k=1,kk
      kn=k+nn
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            ocetra(i,j,k,:) = trc(i,j,kn,itrbgc:itrbgc+ntrbgc-1) * bgc_rho(i,j,k)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- pass sediments fields (a two time-level copy of sediment fields
    ! --- is kept outside HAMOCC)
    ! --- ------------------------------------------------------------------

    if (.not. use_sedbypass) then
      nns=(n-1)*ks

      !$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do k=1,ks
        kn=k+nns
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              sedlay(i,j,k,:) = sedlay2(i,j,kn,:)
              powtra(i,j,k,:) = powtra2(i,j,kn,:)
              burial(i,j,:)   = burial2(i,j,n,:)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    endif

    ! --- ------------------------------------------------------------------
    ! --- pass atmosphere fields if required (a two time-level copy of
    ! --- atmosphere fields is kept outside HAMOCC)
    ! --- ------------------------------------------------------------------

    if (use_BOXATM) then
      !$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          atm(i,j,:) = atm2(i,j,n,:)
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif

  end subroutine blom2hamocc


  subroutine hamocc2blom(m,n,mm,nn)

    !******************************************************************************
    ! Interface between BLOM and HAMOCC.
    ! Pass flux and tracer fields back from HAMOCC to BLOM.
    ! The local HAMOCC arrays are copied back in the appropriate
    ! time-level of the tracer field. Note that also sediment fields
    ! are copied back, since a two time-level copy of sediment fields
    ! is kept outside HAMOCC. For the sediment fields the same time-
    ! smothing as for the tracer field (i.e. analog to tmsmt2.F) is
    ! performed to avoid a seperation of the two time levels.
    !
    ! J. Schwinger      *GFI, UiB        2014-05-21 initial version
    ! Modified:
    ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    ! - added sediment bypass preprocessor option
    ! M. Bentsen,       *NORCE, Bergen*          2020-05-03
    ! - changed ocean model from MICOM to BLOM
    ! T. Torsvik,       *University of Bergen*   2021-08-26
    ! - integrate subroutine into module mo_intfcblom
    !******************************************************************************

    use mod_xc,        only: ii,jj,kk,ifp,ilp,isp
    use mod_tracers,   only: ntrbgc,itrbgc,trc
    use mod_tmsmt,     only: wts1, wts2
    use mo_carbch,     only: ocetra,atm
    use mo_param1_bgc, only: ks,nsedtra,npowtra,natm
    use mo_sedmnt,     only: sedlay,powtra,sedhpl,burial

    ! Arguments
    integer, intent(in) :: m,n,mm,nn

    ! Local variables
    integer :: i,j,k,l,nns,mms,kn,km

    ! --- ------------------------------------------------------------------
    ! --- pass tracer fields to ocean model; convert kmol/m^3 -> mol/kg
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(k,kn,l,i)
    do k=1,kk
      kn=k+nn
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            trc(i,j,kn,itrbgc:itrbgc+ntrbgc-1) = ocetra(i,j,k,:)/bgc_rho(i,j,k)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- apply time smoothing for sediment fields and pass them back
    ! --- ------------------------------------------------------------------

    if (.not. use_sedbypass) then
      nns=(n-1)*ks
      mms=(m-1)*ks

      !$OMP PARALLEL DO PRIVATE(k,km,kn,l,i)
      do k=1,ks
        km=k+mms
        kn=k+nns
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))              ! time smoothing (analog to tmsmt2.F)
              sedlay2(i,j,km,:) = wts1*sedlay2(i,j,km,:)   &  ! mid timelevel
                   + wts2*sedlay2(i,j,kn,:)   &  ! old timelevel
                   + wts2*sedlay(i,j,k,:)        ! new timelevel
              powtra2(i,j,km,:) = wts1*powtra2(i,j,km,:)   &
                   + wts2*powtra2(i,j,kn,:)   &
                   + wts2*powtra(i,j,k,:)
              burial2(i,j,m,:)  = wts1*burial2(i,j,m,:)    &
                   + wts2*burial2(i,j,n,:)    &
                   + wts2*burial(i,j,:)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do k=1,ks
        kn=k+nns
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              sedlay2(i,j,kn,:) = sedlay(i,j,k,:)  ! new time level replaces old time level here
              powtra2(i,j,kn,:) = powtra(i,j,k,:)
              burial2(i,j,n,:)  = burial(i,j,:)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    endif ! .not. use_sedbypass

    ! --- ------------------------------------------------------------------
    ! --- apply time smoothing for atmosphere fields if required
    ! --- ------------------------------------------------------------------

    if (use_BOXATM) then
      !$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii                                  ! time smoothing (analog to tmsmt2.F)
          atm2(i,j,m,:) = wts1*atm2(i,j,m,:)   &   ! mid timelevel
                        + wts2*atm2(i,j,n,:)   &   ! old timelevel
                        + wts2*atm(i,j,:)          ! new timelevel
        enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          atm2(i,j,n,:) = atm(i,j,:)  ! new time level replaces old time level here
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif

  end subroutine hamocc2blom

end module mo_intfcblom
