! ------------------------------------------------------------------------------
! Copyright (C) 2022 Mats Bentsen
!
! This file is part of BLOM.
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
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_nuopc_methods
! ------------------------------------------------------------------------------
! This module contains routines operating on BLOM data structures needed by the
! NUOPC cap.
! ------------------------------------------------------------------------------

   use mod_types,      only: r8
   use mod_constants,  only: rearth, onem, L_mks2cgs
   use mod_time,       only: nstep, baclin, delt1, dlt
   use mod_xc
   use mod_grid,       only: scuy, scvx, scp2, scuxi, scvyi, plon, plat, cosang, sinang
   use mod_state,      only: u, v, dp, temp, saln, pbu, pbv, ubflxs, vbflxs, sealv
   use mod_forcing,    only: wavsrc_opt, wavsrc_extern, sprfac, prfac, flxco2,flxdms, flxbrf
   use mod_difest,     only: obldepth
   use mod_vcoord,     only: vcoord_type_tag, isopyc_bulkml, cntiso_hybrid
   use mod_cesm,       only: frzpot, mltpot, &
                             swa_da, nsf_da, hmlt_da, lip_da, sop_da, eva_da, &
                             rnf_da, rfi_da, fmltfz_da, sfl_da, ztx_da, mty_da, &
                             ustarw_da, slp_da, abswnd_da, ficem_da, lamult_da, &
                             lasl_da, ustokes_da, vstokes_da, atmco2_da, atmbrf_da, &
                             flxdms_da, flxbrf_da, l1ci, l2ci
   use mod_utility,    only: util1, util2
   use mod_checksum,   only: csdiag, chksummsk
   use mod_fill_global,only: fill_global
   use shr_const_mod,  only: SHR_CONST_RHOSW, SHR_CONST_LATICE, SHR_CONST_TKFRZ
#ifdef HAMOCC
   use mo_carbch,      only: ocetra
   use mo_control_bgc, only: use_BROMO, ocn_co2_type
#endif

   implicit none
   private

   ! Parameters.
   character(len=*), parameter :: modname = '(mod_nuopc_methods)'
   real(r8), parameter :: iL_mks2cgs = 1./L_mks2cgs

   type :: fldlist_type
      character(len=128) :: stdname
      integer :: ungridded_lbound = 0
      integer :: ungridded_ubound = 0
      real(r8), dimension(:), pointer :: dataptr
   end type fldlist_type
   integer, parameter :: fldsMax = 100

   real(r8), dimension(:), allocatable :: mod2med_areacor, med2mod_areacor
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      acc_u, acc_v, acc_dhdx, acc_dhdy, acc_t, acc_s, acc_frzpot, acc_bld, &
      acc_fco2, acc_fdms, acc_fbrf
   real(r8) :: tlast_coupled
   integer :: jjcpl

   public :: fldlist_type, fldsmax, tlast_coupled, &
             blom_logwrite, blom_getgindex, blom_checkmesh, blom_setareacor, &
             blom_getglobdim, blom_getprecipfact, blom_accflds, &
             blom_advertise_imports, blom_advertise_exports, &
             blom_importflds, blom_exportflds

   ! Indices for import fields
   integer :: &
        index_Si_ifrac    = -1, &
        index_So_duu10n   = -1, &
        index_Fioi_melth  = -1, &
        index_Fioi_meltw  = -1, &
        index_Fioi_salt   = -1, &
        index_Fioi_bcpho  = -1, &
        index_Fioi_bcphi  = -1, &
        index_Fioi_flxdst = -1, &
        index_Foxx_rofl   = -1, &
        index_Foxx_rofi   = -1, &
        index_Foxx_tauy   = -1, &
        index_Foxx_taux   = -1, &
        index_Foxx_lat    = -1, &
        index_Foxx_sen    = -1, &
        index_Foxx_lwup   = -1, &
        index_Foxx_evap   = -1, &
        index_Foxx_swnet  = -1, &
        index_Sw_lamult   = -1, &
        index_Sw_ustokes  = -1, &
        index_Sw_vstokes  = -1, &
        index_Sw_hstokes  = -1, &
        index_Faxa_lwdn   = -1, &
        index_Faxa_snow   = -1, &
        index_Faxa_rain   = -1, &
        index_Sa_pslv     = -1, &
        index_Sa_co2diag  = -1, &
        index_Sa_co2prog  = -1

   ! Indices for export fields
   integer  :: &
        index_So_omask   = -1, &
        index_So_u       = -1, &
        index_So_v       = -1, &
        index_So_dhdx    = -1, &
        index_So_dhdy    = -1, &
        index_So_t       = -1, &
        index_So_s       = -1, &
        index_So_bldepth = -1, &
        index_Fioo_q     = -1, &
        index_Faoo_dms   = -1, &
        index_Faoo_brf   = -1, &
        index_Faoo_fco2_ocn = -1

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

  subroutine fldlist_add(num, fldlist, stdname, index, ungridded_lbound, ungridded_ubound)

    ! ---------------------------------------------------------------------------
    ! Add to list of field information.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    integer           , intent(inout) :: num
    type(fldlist_type), intent(inout) :: fldlist(:)
    character(len=*)  , intent(in)    :: stdname
    integer           , intent(out)   :: index
    integer, optional , intent(in)    :: ungridded_lbound, ungridded_ubound

    ! Local parameters.
    character(len=*), parameter :: &
         subname = modname//':(fldlist_add)'

    ! Local variables.
    integer :: rc

    num = num + 1
    if (num > fldsMax) then
      write(lp,'(a,3i6,2(f21.13,3x),d21.5)') subname// &
           ': BLOM ERROR: number of fields exceeds fldsMax for '//trim(stdname)
      call xchalt(subname)
      stop subname
    endif
    fldlist(num)%stdname = trim(stdname)

    index = num

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
      fldlist(num)%ungridded_lbound = ungridded_lbound
      fldlist(num)%ungridded_ubound = ungridded_ubound
    endif

  end subroutine fldlist_add

  ! ---------------------------------------------------------------------------
  ! Public procedures.
  ! ---------------------------------------------------------------------------

  subroutine blom_advertise_imports(flds_scalar_name, fldsToOcn_num, fldsToOcn, &
       flds_co2a, flds_co2c)

    ! -------------------------------------------------------------------
    ! Determine fldsToOcn for import fields
    ! -------------------------------------------------------------------

    character(len=*)   , intent(in)                  :: flds_scalar_name
    integer            , intent(inout)               :: fldsToOcn_num
    type(fldlist_type) , intent(inout), dimension(:) :: fldsToOcn
    logical            , intent(in)                  :: flds_co2a
    logical            , intent(in)                  :: flds_co2c

    integer :: index_scalar

    call fldlist_add(fldsToOcn_num, fldsToOcn, trim(flds_scalar_name), index_scalar)

    ! From ice:
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_ifrac'   , index_Si_ifrac )
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_melth' , index_Fioi_melth)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_meltw' , index_Fioi_meltw)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_salt'  , index_Fioi_salt)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_bcpho' , index_Fioi_bcpho)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_bcphi' , index_Fioi_bcphi)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_flxdst', index_Fioi_flxdst)

    ! From river:
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofl', index_Foxx_rofl)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofi', index_Foxx_rofi)

    ! From fields computed mediator:
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'So_duu10n'  , index_So_duu10n)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_tauy'  , index_Foxx_tauy)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_taux'  , index_Foxx_taux)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lat'   , index_Foxx_lat)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_sen'   , index_Foxx_sen)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lwup'  , index_Foxx_lwup)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_evap'  , index_Foxx_evap)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_swnet' , index_Foxx_swnet)

    ! From wave:
    if (wavsrc_opt == wavsrc_extern) then
      call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_lamult'  , index_Sw_lamult)
      call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_ustokes' , index_Sw_ustokes)
      call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_vstokes' , index_Sw_vstokes)
      call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_hstokes' , index_Sw_hstokes)
    end if

    ! From atmosphere:
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_pslv'   , index_Sa_pslv  )
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_lwdn' , index_Faxa_lwdn)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_snow' , index_Faxa_snow)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_rain' , index_Faxa_rain)

    ! From atm co2 fields:
    if (flds_co2a .or. flds_co2c) then
      call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2diag' ,index_Sa_co2diag)
      call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2prog', index_Sa_co2prog)
    endif

  end subroutine blom_advertise_imports

  subroutine blom_advertise_exports(flds_scalar_name, fldsFrOcn_num, fldsFrOcn)
    ! -------------------------------------------------------------------
    ! Determine fldsToOcn for export fields
    ! -------------------------------------------------------------------

    character(len=*)                 , intent(in)    :: flds_scalar_name
    integer                          , intent(inout) :: fldsFrOcn_num
    type(fldlist_type), dimension(:) , intent(inout) :: fldsFrOcn

    integer :: index_scalar

    call fldlist_add(fldsFrOcn_num, fldsFrOcn, trim(flds_scalar_name), index_scalar)

    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_omask'      , index_So_omask)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_t'          , index_So_t)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_u'          , index_So_u)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_v'          , index_So_v)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_s'          , index_So_s)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdx'       , index_So_dhdx)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdy'       , index_So_dhdy)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_bldepth'    , index_So_bldepth)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Fioo_q'        , index_Fioo_q)
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Faoo_fco2_ocn' , index_Faoo_fco2_ocn)
#ifdef HAMOCC
    ! call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faoo_dms', index_Faoo_dms)
    ! if (use_BROMO) then
    !    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faoo_brf', index_Faoo_brf)
    ! end if
#endif
  end subroutine blom_advertise_exports

  subroutine blom_logwrite(msg)
    ! ---------------------------------------------------------------------------
    ! Write message string to standard out from master PE.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    character(len=*), intent(in) :: msg

    if (mnproc == 1) write(lp,'(a)') trim(msg)

  end subroutine blom_logwrite

  subroutine blom_getgindex(gindex)
    ! ---------------------------------------------------------------------------
    ! Get global index space for the computational domain.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    integer, allocatable, dimension(:), intent(out) :: gindex

    ! Local variables.
    integer :: mproc_next, i, j, n

    ! Set the j-extent of the local ocean domain to be exchanged. Needed
    ! because of duplication of the last global domain row when using a
    ! tripolar grid.
    if (nreg == 2 .and. nproc == jpr) then
      jjcpl = jj - 1
    else
      jjcpl = jj
    endif

    ! Create the global index space for the computational domain. Also append
    ! indices of eliminated grid cells adjacent to the domain and with larger
    ! global i-index.
    mproc_next = mod(mproc, ipr) + 1
    do while (ii_pe(mproc_next,nproc) == 0)
      mproc_next = mod(mproc_next, ipr) + 1
    enddo
    allocate(gindex(mod(i0_pe(mproc_next,nproc) - i0 + itdm, itdm)*jjcpl))
    n = 0
    do j = 1, jjcpl
      do i = 1, ii
        n = n + 1
        gindex(n) = (j0 + j - 1)*itdm + i0 + i
      enddo
    enddo
    do j = 1, jjcpl
      do i = ii + 1, mod(i0_pe(mproc_next,nproc) - i0 + itdm, itdm)
        n = n + 1
        gindex(n) = (j0 + j - 1)*itdm + mod(i0 + i - 1, itdm) + 1
      enddo
    enddo

  end subroutine blom_getgindex

  subroutine blom_checkmesh(lonmesh, latmesh, maskmesh)
    ! ---------------------------------------------------------------------------
    ! Check for consistency of lat, lon and mask between mediator mesh and model
    ! grid.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    real(r8), dimension(:), pointer, intent(in) :: lonmesh, latmesh
    integer, dimension(:), pointer, intent(in) :: maskmesh

    ! Local parameters.
    character(len=*), parameter :: &
         subname = modname//':(blom_checkmesh)'

    ! Local variables.
    real(r8) :: diff_lon, diff_lat
    integer :: mproc_next, i, j, n

    do j = 1, jjcpl
      do i = 1, ii
        n = (j - 1)*ii + i
        diff_lon = abs(mod(lonmesh(n) - plon(i,j),360._r8))
        if (diff_lon > 1.e-3_r8) then
          write(lp,'(a,3i6,2(f21.13,3x),d21.5)') subname// &
               ': BLOM ERROR: n, i, j, lonmesh(n), plon(i,j), diff_lon = ', &
               n, i, j, lonmesh(n), plon(i,j), diff_lon
          call xchalt(subname)
          stop subname
        endif
        diff_lat = abs(latmesh(n) - plat(i,j))
        if (diff_lat > 1.e-3_r8) then
          write(lp,'(a,3i6,2(f21.13,3x),d21.5)') subname// &
               ': BLOM ERROR: n, i, j, latmesh(n), plat(i,j), diff_lat = ', &
               n, i, j, latmesh(n), plat(i,j), diff_lat
          call xchalt(subname)
          stop subname
        endif
        if (maskmesh(n) /= ip(i,j)) then
          write(lp,'(a,3i6,2(f21.13,3x),d21.5)') subname// &
               ': BLOM ERROR: n, i, j, maskmesh(n), ip(i,j) = ', &
               n, i, j, maskmesh(n), ip(i,j)
          call xchalt(subname)
          stop subname
        endif
      enddo
    enddo

  end subroutine blom_checkmesh

  subroutine blom_getprecipfact(precip_fact_provided, precip_fact)
    ! ---------------------------------------------------------------------------
    ! Get precipitation factor.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    logical, intent(out)  :: precip_fact_provided
    real(r8), intent(out) :: precip_fact

    precip_fact_provided = sprfac
    precip_fact = prfac

  end subroutine blom_getprecipfact

  subroutine blom_getglobdim(nx_global, ny_global)
    ! ---------------------------------------------------------------------------
    ! Get global dimensions of export/import domain.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    integer, intent(out) :: nx_global, ny_global

    nx_global = itdm
    if (nreg == 2) then
      ny_global = jtdm - 1
    else
      ny_global = jtdm
    endif

  end subroutine blom_getglobdim

  subroutine blom_setareacor(areamesh, maskmesh)
    ! ---------------------------------------------------------------------------
    ! Set flux area correction factors.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    real(r8), dimension(:), pointer, intent(in) :: areamesh
    integer, dimension(:), pointer, intent(in) :: maskmesh

    ! Local parameters.
    character(len=*), parameter :: &
         subname = modname//':(blom_setareacor)'

    ! Local variables.
    real(r8) :: areamodel, &
         max_mod2med_areacor, max_med2mod_areacor, &
         min_mod2med_areacor, min_med2mod_areacor
    integer :: i, j, n

    allocate(mod2med_areacor(size(areamesh)), &
         med2mod_areacor(size(areamesh)))
    mod2med_areacor(:) = 1._r8
    med2mod_areacor(:) = 1._r8

    !$omp parallel do private(i, n)
    do j = 1, jjcpl
      do i = 1, ii
        n = (j - 1)*ii + i
        if (maskmesh(n) /= 0) then
          areamodel = scp2(i,j)/(rearth*rearth)
          mod2med_areacor(n) = areamodel/areamesh(n)
          med2mod_areacor(n) = areamesh(n)/areamodel
        endif
      enddo
    enddo
    !$omp end parallel do

    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call xcmax(max_mod2med_areacor)
    call xcmin(min_mod2med_areacor)
    call xcmax(max_med2mod_areacor)
    call xcmin(min_med2mod_areacor)
    if (mnproc == 1) then
      write(lp,'(a,2g23.15)') &
           subname//': min_mod2med_areacor, max_mod2med_areacor ', &
           min_mod2med_areacor, max_mod2med_areacor
      write(lp,'(a,2g23.15)') &
           subname//': min_med2mod_areacor, max_med2mod_areacor ', &
           min_med2mod_areacor, max_med2mod_areacor
    endif

  end subroutine blom_setareacor

  subroutine blom_accflds
    ! ---------------------------------------------------------------------------
    ! Accumulate export fields to be averaged before sent to the mediator.
    ! ---------------------------------------------------------------------------

    ! Local parameters.
    character(len=*), parameter :: &
         subname = modname//':(blom_accflds)'

    ! Local variables.
    real(r8) :: q
    integer m, n, mm, nn, k1m, k1n, i, j, l
    logical :: first_call = .true.

    ! ------------------------------------------------------------------------
    ! Set accumulation arrays to zero if this is the first call after a
    ! coupling interval.
    ! ------------------------------------------------------------------------

    if (tlast_coupled == 0._r8) then
      acc_u     (:,:) = 0._r8
      acc_v     (:,:) = 0._r8
      acc_dhdx  (:,:) = 0._r8
      acc_dhdy  (:,:) = 0._r8
      acc_t     (:,:) = 0._r8
      acc_s     (:,:) = 0._r8
      acc_frzpot(:,:) = 0._r8
      acc_bld   (:,:) = 0._r8
      acc_fco2  (:,:) = 0._r8
      acc_fdms  (:,:) = 0._r8
      acc_fbrf  (:,:) = 0._r8
    endif

    ! ------------------------------------------------------------------------
    ! Accumulate fields in send buffer
    ! ------------------------------------------------------------------------

    m = mod(nstep + 1, 2) + 1
    n = mod(nstep    , 2) + 1
    mm = (m - 1)*kk
    nn = (n - 1)*kk
    k1m = 1 + mm
    k1n = 1 + nn

    call xctilr(sealv, 1,1, 1,1, halo_ps)

    !$omp parallel do private(l, i)
    do j = 1, jj
      do l = 1, isu(j)
        do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
          acc_u(i,j) = acc_u(i,j) &
               + ( u(i,j,k1n) &
               + (ubflxs(i,j,m) + ubflxs(i,j,n))*dlt &
               /(pbu(i,j,n)*scuy(i,j)*delt1))*baclin
          acc_dhdx(i,j) = acc_dhdx(i,j) &
               + (sealv(i,j) - sealv(i-1,j))*scuxi(i,j)*baclin
        enddo
      enddo
      do l = 1, isv(j)
        do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
          acc_v(i,j) = acc_v(i,j) &
               + ( v(i,j,k1n)  &
               + (vbflxs(i,j,m) + vbflxs(i,j,n))*dlt &
               /(pbv(i,j,n)*scvx(i,j)*delt1))*baclin
          acc_dhdy(i,j) = acc_dhdy(i,j) &
               + (sealv(i,j) - sealv(i,j-1))*scvyi(i,j)*baclin
        enddo
      enddo
      do l = 1, isp(j)
        do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
          acc_t(i,j) = acc_t(i,j) + temp(i,j,k1n)*baclin
          acc_s(i,j) = acc_s(i,j) + saln(i,j,k1n)*baclin
          acc_frzpot(i,j) = acc_frzpot(i,j) + frzpot(i,j)
        enddo
      enddo
    enddo
    !$omp end parallel do

    select case (vcoord_type_tag)
    case (isopyc_bulkml)
      q = baclin/onem
      !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            acc_bld(i,j) = (dp(i,j,1+nn) + dp(i,j,2+nn))*q
          enddo
        enddo
      enddo
      !$omp end parallel do
    case (cntiso_hybrid)
      !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            acc_bld(i,j) = OBLdepth(i,j)*baclin
          enddo
        enddo
      enddo
      !$omp end parallel do
    case default
      if (mnproc == 1.and. first_call) then
        write(lp,*) subname//': unsupported vertical coordinate!'
      end if
      call xcstop(subname)
      stop subname
    end select

    if (index_Faoo_fco2_ocn > 0) then
      !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            acc_fco2(i,j) = acc_fco2(i,j) + flxco2(i,j)*baclin
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif

    if (index_Faoo_dms > 0) then
      !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            acc_fdms(i,j) = acc_fdms(i,j) + flxdms(i,j)*baclin
          enddo
        enddo
      enddo
      !$omp end parallel do
    end if

    if (index_Faoo_brf > 0) then
      !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            acc_fbrf(i,j) = acc_fbrf(i,j) + flxbrf(i,j)*baclin
          enddo
        enddo
      enddo
      !$omp end parallel do
    end if

    ! ------------------------------------------------------------------------
    ! Increment time since last coupling.
    ! ------------------------------------------------------------------------

    tlast_coupled = tlast_coupled + baclin

    if (first_call) then
      first_call = .false.
    end if

  end subroutine blom_accflds

  subroutine blom_importflds(fldlist_num, fldlist)
    ! ---------------------------------------------------------------------------
    ! Import fields from mediator to BLOM arrays.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    integer, intent(in) :: fldlist_num
    type(fldlist_type), dimension(:), intent(in) :: fldlist

    ! Local parameters.
    character(len=*), parameter :: &
         subname = modname//':(blom_importflds)'
    real(r8), parameter :: &
         mval = - 1.e12_r8, &
         fval = - 1.e13_r8
    logical :: first_call = .true.

    ! Local variables.
    real(r8) :: afac, utmp, vtmp
    integer :: n, i, j, l
    integer :: index_co2

    ! Update time level indices.
    if (l1ci == 1 .and. l2ci == 1) then
      l1ci = 2
      l2ci = 2
    else
      l1ci = l2ci
      l2ci = 3 - l2ci
    endif

    !$omp parallel do private(i, n, afac, utmp, vtmp)
    do j = 1, jjcpl
      do i = 1, ii
        if     (ip(i,j) == 0) then
          util1(i,j) = mval
          util2(i,j) = mval
          ustarw_da(i,j,l2ci) = mval
        elseif (cplmsk(i,j) == 0) then
          util1(i,j) = fval
          util2(i,j) = fval
          ustarw_da(i,j,l2ci) = fval
        else
          n = (j - 1)*ii + i
          afac = med2mod_areacor(n)

          utmp = fldlist(index_Foxx_taux)%dataptr(n)*afac
          vtmp = fldlist(index_Foxx_tauy)%dataptr(n)*afac
          util1(i,j) =   utmp*cosang(i,j) + vtmp*sinang(i,j)
          util2(i,j) = - utmp*sinang(i,j) + vtmp*cosang(i,j)

          ! Friction velocity [m s-1].
          ustarw_da(i,j,l2ci) = sqrt(sqrt(utmp*utmp + vtmp*vtmp) &
               /SHR_CONST_RHOSW)
        endif
      enddo
    enddo
    !$omp end parallel do

    call fill_global(mval, fval, halo_pv, util1)
    call fill_global(mval, fval, halo_pv, util2)
    call fill_global(mval, fval, halo_ps, ustarw_da(1-nbdy,1-nbdy,l2ci))

    call xctilr(util1, 1,1, 1,1, halo_pv)
    call xctilr(util2, 1,1, 1,1, halo_pv)

    !$omp parallel do private(l, i)
    do j = 1, jj
      do l = 1, isu(j)
        do i = max(1,ifu(j,l)), min(ii,ilu(j,l))
          ! x-component of momentum flux [kg m-1 s-2].
          ztx_da(i,j,l2ci) = .5_r8*(util1(i-1,j) + util1(i,j))
        enddo
      enddo
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)), min(ii,ilv(j,l))
          ! y-component of momentum flux [kg m-1 s-2].
          mty_da(i,j,l2ci) = .5_r8*(util2(i,j-1) + util2(i,j))
        enddo
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(i, n, afac)
    do j = 1, jjcpl
      do i = 1, ii

        if (ip(i,j) == 0) then
          lip_da(i,j,l2ci) = mval
          sop_da(i,j,l2ci) = mval
          eva_da(i,j,l2ci) = mval
          rnf_da(i,j,l2ci) = mval
          rfi_da(i,j,l2ci) = mval
          fmltfz_da(i,j,l2ci) = mval
          sfl_da(i,j,l2ci) = mval
          swa_da(i,j,l2ci) = mval
          nsf_da(i,j,l2ci) = mval
          hmlt_da(i,j,l2ci) = mval
          slp_da(i,j,l2ci) = mval
          abswnd_da(i,j,l2ci) = mval
          ficem_da(i,j,l2ci) = mval
        elseif (cplmsk(i,j) == 0) then
          lip_da(i,j,l2ci) = 0._r8
          sop_da(i,j,l2ci) = 0._r8
          eva_da(i,j,l2ci) = 0._r8
          rnf_da(i,j,l2ci) = 0._r8
          rfi_da(i,j,l2ci) = 0._r8
          fmltfz_da(i,j,l2ci) = 0._r8
          sfl_da(i,j,l2ci) = 0._r8
          swa_da(i,j,l2ci) = 0._r8
          nsf_da(i,j,l2ci) = 0._r8
          hmlt_da(i,j,l2ci) = 0._r8
          slp_da(i,j,l2ci) = fval
          abswnd_da(i,j,l2ci) = fval
          ficem_da(i,j,l2ci) = fval
        else
          n = (j - 1)*ii + i
          afac = med2mod_areacor(n)

          ! Liquid water flux, positive downwards [kg m-2 s-1].
          lip_da(i,j,l2ci) = fldlist(index_Faxa_rain)%dataptr(n)*afac

          ! Solid precipitation, positive downwards [kg m-2 s-1].
          sop_da(i,j,l2ci) = fldlist(index_Faxa_snow)%dataptr(n)*afac

          ! Evaporation, positive downwards [kg m-2 s-1].
          eva_da(i,j,l2ci) = fldlist(index_Foxx_evap)%dataptr(n)*afac

          ! Liquid runoff, positive downwards [kg m-2 s-1].
          rnf_da(i,j,l2ci) = fldlist(index_Foxx_rofl)%dataptr(n)*afac

          ! Frozen runoff, positive downwards [kg m-2 s-1].
          rfi_da(i,j,l2ci) = fldlist(index_Foxx_rofi)%dataptr(n)*afac

          ! Fresh water due to melting/freezing, positive downwards
          ! [kg m-2 s-1].
          fmltfz_da(i,j,l2ci) = fldlist(index_Fioi_meltw)%dataptr(n)*afac

          ! Salt flux, positive downwards [kg m-2 s-1].
          sfl_da(i,j,l2ci) = fldlist(index_Fioi_salt)%dataptr(n)*afac

          ! Shortwave heat flux, positive downwards [W m-2].
          swa_da(i,j,l2ci) = fldlist(index_Foxx_swnet)%dataptr(n)*afac

          ! Non-solar heat flux, positive downwards [W m-2].
          nsf_da(i,j,l2ci) = ( fldlist(index_Foxx_lat)%dataptr(n) &
               + fldlist(index_Foxx_sen)%dataptr(n) &
               + fldlist(index_Foxx_lwup)%dataptr(n) &
               + fldlist(index_Faxa_lwdn)%dataptr(n) &
               - ( fldlist(index_Faxa_snow)%dataptr(n) &
               + fldlist(index_Foxx_rofi)%dataptr(n)) &
               *SHR_CONST_LATICE)*afac

          ! Heat flux due to melting, positive downwards [W m-2].
          hmlt_da(i,j,l2ci) = fldlist(index_Fioi_melth)%dataptr(n)*afac

          ! Sea level pressure [kg m-1 s-2].
          slp_da(i,j,l2ci) = fldlist(index_Sa_pslv)%dataptr(n)

          ! 10m wind speed [m s-1].
          abswnd_da(i,j,l2ci) = sqrt(fldlist(index_So_duu10n)%dataptr(n))

          ! Ice fraction [].
          ficem_da(i,j,l2ci) = fldlist(index_Si_ifrac)%dataptr(n)

        endif

      enddo
    enddo
    !$omp end parallel do

    if (nreg == 2) then
      call xctilr(lip_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(sop_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(eva_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(rnf_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(rfi_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(fmltfz_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(sfl_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(swa_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(nsf_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
      call xctilr(hmlt_da(1-nbdy,1-nbdy,l2ci), 1,1, 0,0, halo_ps)
    endif

    call fill_global(mval, fval, halo_ps, slp_da(1-nbdy,1-nbdy,l2ci))
    call fill_global(mval, fval, halo_ps, abswnd_da(1-nbdy,1-nbdy,l2ci))
    call fill_global(mval, fval, halo_ps, ficem_da(1-nbdy,1-nbdy,l2ci))

    if (wavsrc_opt == wavsrc_extern) then
      !$omp parallel do private(i, n, utmp, vtmp)
      do j = 1, jjcpl
        do i = 1, ii
          if     (ip(i,j) == 0) then
            util1(i,j) = mval
            util2(i,j) = mval
            lamult_da(i,j,l2ci) = mval
            lasl_da(i,j,l2ci) = mval
          elseif (cplmsk(i,j) == 0) then
            util1(i,j) = fval
            util2(i,j) = fval
            lamult_da(i,j,l2ci) = fval
            lasl_da(i,j,l2ci) = fval
          else
            n = (j - 1)*ii + i

            utmp = fldlist(index_Sw_ustokes)%dataptr(n)
            vtmp = fldlist(index_Sw_vstokes)%dataptr(n)
            util1(i,j) =   utmp*cosang(i,j) + vtmp*sinang(i,j)
            util2(i,j) = - utmp*sinang(i,j) + vtmp*cosang(i,j)

            ! Langmuir enhancement factor [].
            lamult_da(i,j,l2ci) = fldlist(index_Sw_lamult)%dataptr(n)

            ! Surface layer averaged Langmuir number [].
            lasl_da(i,j,l2ci) = fldlist(index_Sw_hstokes)%dataptr(n)

          endif
        enddo
      enddo
      !$omp end parallel do

      call fill_global(mval, fval, halo_pv, util1)
      call fill_global(mval, fval, halo_pv, util2)
      call fill_global(mval, fval, halo_ps, lamult_da(1-nbdy,1-nbdy,l2ci))
      call fill_global(mval, fval, halo_ps, lasl_da(1-nbdy,1-nbdy,l2ci))

      call xctilr(util1, 1,1, 1,1, halo_pv)
      call xctilr(util2, 1,1, 1,1, halo_pv)

      !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isu(j)
          do i = max(1,ifu(j,l)), min(ii,ilu(j,l))
            ! x-component of surface Stokes drift [m s-1].
            ustokes_da(i,j,l2ci) = .5_r8*(util1(i-1,j) + util1(i,j))
          enddo
        enddo
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)), min(ii,ilv(j,l))
            ! y-component of surface Stokes drift [m s-1].
            vstokes_da(i,j,l2ci) = .5_r8*(util2(i,j-1) + util2(i,j))
          enddo
        enddo
      enddo
      !$omp end parallel do

    end if

    ! CO2 flux

    index_co2 = -1
#ifdef HAMOCC
    if (ocn_co2_type == 'diagnostic' .and. index_Sa_co2diag > 0) then
      index_co2 = index_Sa_co2diag
    else if (ocn_co2_type == 'prognostic' .and. index_Sa_co2prog > 0) then
      index_co2 = index_Sa_co2prog
    end if
#endif
    if (index_co2 > 0) then
      !$omp parallel do private(i, n)
      do j = 1, jjcpl
        do i = 1, ii
          if     (ip(i,j) == 0) then
            atmco2_da(i,j,l2ci) = mval
          elseif (cplmsk(i,j) == 0) then
            atmco2_da(i,j,l2ci) = fval
          else
            n = (j - 1)*ii + i
            ! Atmospheric co2 concentration [ppmv?]
            atmco2_da(i,j,l2ci) = fldlist(index_co2)%dataptr(n)
          endif
        enddo
      enddo
      !$omp end parallel do
      call fill_global(mval, fval, halo_ps, atmco2_da(1-nbdy,1-nbdy,l2ci))
      if (mnproc == 1 .and. first_call) then
        write(lp,*) subname//': atmospheric co2 obtained from mediator'
      end if
    else
      !$omp parallel do private(i)
      do j = 1, jj
        do i = 1, ii
          if (ip(i,j) == 0) then
            atmco2_da(i,j,l2ci) = mval
          else
            atmco2_da(i,j,l2ci) = -1
          endif
        enddo
      enddo
      !$omp end parallel do
      if (mnproc == 1 .and. first_call)  then
        write(lp,*) subname//': atmospheric co2 not obtained from mediator'
      endif
    end if

    if (csdiag) then
      if (mnproc == 1 .and. first_call) then
        write(lp,*) subname//':'
      endif
      call chksummsk(ustarw_da(1-nbdy,1-nbdy,l2ci),ip,1,'ustarw')
      call chksummsk(ztx_da(1-nbdy,1-nbdy,l2ci),iu,1,'ztx')
      call chksummsk(mty_da(1-nbdy,1-nbdy,l2ci),iv,1,'mty')
      call chksummsk(lip_da(1-nbdy,1-nbdy,l2ci),ip,1,'lip')
      call chksummsk(sop_da(1-nbdy,1-nbdy,l2ci),ip,1,'sop')
      call chksummsk(eva_da(1-nbdy,1-nbdy,l2ci),ip,1,'eva')
      call chksummsk(rnf_da(1-nbdy,1-nbdy,l2ci),ip,1,'rnf')
      call chksummsk(rfi_da(1-nbdy,1-nbdy,l2ci),ip,1,'rfi')
      call chksummsk(fmltfz_da(1-nbdy,1-nbdy,l2ci),ip,1,'fmltfz')
      call chksummsk(sfl_da(1-nbdy,1-nbdy,l2ci),ip,1,'sfl')
      call chksummsk(swa_da(1-nbdy,1-nbdy,l2ci),ip,1,'swa')
      call chksummsk(nsf_da(1-nbdy,1-nbdy,l2ci),ip,1,'nsf')
      call chksummsk(hmlt_da(1-nbdy,1-nbdy,l2ci),ip,1,'hmlt')
      call chksummsk(slp_da(1-nbdy,1-nbdy,l2ci),ip,1,'slp')
      call chksummsk(abswnd_da(1-nbdy,1-nbdy,l2ci),ip,1,'abswnd')
      call chksummsk(ficem_da(1-nbdy,1-nbdy,l2ci),ip,1,'ficem')
      call chksummsk(atmco2_da(1-nbdy,1-nbdy,l2ci),ip,1,'atmco2')
    endif

    if (first_call) then
      first_call = .false.
    end if

  end subroutine blom_importflds

  subroutine blom_exportflds(fldlist_num, fldlist)
    ! ---------------------------------------------------------------------------
    ! Export from BLOM arrays to mediator fields.
    ! ---------------------------------------------------------------------------

    ! Input/output arguments.
    integer, intent(in) :: fldlist_num
    type(fldlist_type), dimension(:), intent(in) :: fldlist

    ! Local parameters.
    character(len=*), parameter :: &
         subname = modname//':(blom_exportflds)'

    ! Local variables.
    real(r8) :: tfac, utmp, vtmp
    integer  :: n, l, i, j
    logical, save :: first_call = .true.

    tfac = 1._r8/tlast_coupled

    ! ------------------------------------------------------------------------
    ! Provide standard export fields.
    ! ------------------------------------------------------------------------

    call xctilr(acc_u,    1,1, 1,1, halo_uv)
    call xctilr(acc_v,    1,1, 1,1, halo_vv)
    call xctilr(acc_dhdx, 1,1, 1,1, halo_uv)
    call xctilr(acc_dhdy, 1,1, 1,1, halo_vv)

    fldlist(index_So_omask)%dataptr(:) = 0._r8
    fldlist(index_So_u)%dataptr(:) = 0._r8
    fldlist(index_So_v)%dataptr(:) = 0._r8
    fldlist(index_So_dhdx)%dataptr(:) = 0._r8
    fldlist(index_So_dhdy)%dataptr(:) = 0._r8
    fldlist(index_So_t)%dataptr(:) = 0._r8
    fldlist(index_So_s)%dataptr(:) = 0._r8
    fldlist(index_So_bldepth)%dataptr(:) = 0._r8
    fldlist(index_Fioo_q)%dataptr(:) = 0._r8

    !$omp parallel do private(l, i, n, utmp, vtmp)
    do j = 1, jjcpl
      do l = 1, isp(j)
        do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
          n = (j - 1)*ii + i

          ! Ocean mask [].
          fldlist(index_So_omask)%dataptr(n) = 1._r8

          ! Surface velocity, interpolated onto scalar points and rotated [m s-1].
          utmp = .5_r8*(acc_u(i,j) + acc_u(i+1,j))*tfac*iL_mks2cgs
          vtmp = .5_r8*(acc_v(i,j) + acc_v(i,j+1))*tfac*iL_mks2cgs
          fldlist(index_So_u)%dataptr(n) = utmp*cosang(i,j) &
               - vtmp*sinang(i,j)
          fldlist(index_So_v)%dataptr(n) = utmp*sinang(i,j) &
               + vtmp*cosang(i,j)

          ! Surface gradient, interpolated onto scalar points and rotated [].
          utmp = (acc_dhdx(i,j)*iu(i,j) + acc_dhdx(i+1,j)*iu(i+1,j))*tfac &
               /max(1, iu(i,j) + iu(i+1,j))
          vtmp = (acc_dhdy(i,j)*iv(i,j) + acc_dhdy(i,j+1)*iv(i,j+1))*tfac &
               /max(1, iv(i,j) + iv(i,j+1))
          fldlist(index_So_dhdx)%dataptr(n) = utmp*cosang(i,j) &
               - vtmp*sinang(i,j)
          fldlist(index_So_dhdy)%dataptr(n) = utmp*sinang(i,j) &
               + vtmp*cosang(i,j)

          ! Surface temperature [K].
          fldlist(index_So_t)%dataptr(n) = acc_t(i,j)*tfac + SHR_CONST_TKFRZ

          ! Surface salinity [g kg-1].
          fldlist(index_So_s)%dataptr(n) = acc_s(i,j)*tfac

          ! Boundary layer depth [m].
          fldlist(index_So_bldepth)%dataptr(n) = acc_bld(i,j)*tfac

          ! Freezing/melting potential [W m-2].
          if (acc_frzpot(i,j) > 0._r8) then
            fldlist(index_Fioo_q)%dataptr(n) = acc_frzpot(i,j)*tfac*mod2med_areacor(n)
          else
            fldlist(index_Fioo_q)%dataptr(n) = mltpot(i,j)*tfac*mod2med_areacor(n)
          endif

        enddo
      enddo
    enddo
    !$omp end parallel do

    if (index_Faoo_fco2_ocn > 0) then
      ! CO2 flux [kg CO2 m-2 s-1]
      if (associated(fldlist(index_Faoo_fco2_ocn)%dataptr)) then
        fldlist(index_Faoo_fco2_ocn)%dataptr(:) = 0._r8
        !$omp parallel do private(l, i, n)
        do j = 1, jjcpl
          do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
              n = (j - 1)*ii + i
              fldlist(index_Faoo_fco2_ocn)%dataptr(n) = acc_fco2(i,j)*tfac*mod2med_areacor(n)
            enddo
          enddo
        enddo
        !$omp end parallel do
      end if
    else
      if (first_call) then
        if (mnproc == 1 .and. first_call) then
          write(lp,*) subname//': co2 flux not sent to coupler'
        end if
      end if
    end if

    if (index_Faoo_dms > 0) then
      ! dms flux (kmol DMS/m^2/s)
      fldlist(index_Faoo_dms)%dataptr(:) = 0._r8
      !$omp parallel do private(l, i, n)
      do j = 1, jjcpl
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            n = (j - 1)*ii + i
            fldlist(index_Faoo_dms)%dataptr(n) = acc_fdms(i,j)*tfac*mod2med_areacor(n)
          enddo
        enddo
      enddo
      !$omp end parallel do
    else
      if (first_call) then
        if (mnproc == 1 .and. first_call) then
          write(lp,*) subname//': dms flux not sent to coupler'
        end if
      end if
    end if

    if (index_Faoo_brf > 0) then
      ! brf flux (kmol BRF/m^2/s)
      fldlist(index_Faoo_brf)%dataptr(:) = 0._r8
      !$omp parallel do private(l, i, n)
      do j = 1, jjcpl
        do l = 1, isp(j)
          do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            n = (j - 1)*ii + i
            fldlist(index_Faoo_brf)%dataptr(n) = acc_fbrf(i,j)*tfac*mod2med_areacor(n)
          enddo
        enddo
      enddo
      !$omp end parallel do
    else
      if (first_call) then
        if (mnproc == 1 .and. first_call) then
          write(lp,*) subname//': brf flux not sent to coupler'
        end if
      end if
    end if

    if (first_call) then
      first_call = .false.
    end if

    tlast_coupled = 0._r8

  end subroutine blom_exportflds

end module mod_nuopc_methods
