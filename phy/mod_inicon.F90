! ------------------------------------------------------------------------------
! Copyright (C) 2008-2025 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin,
!                         Mariana Vertenstein
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

module mod_inicon
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to advection of layer
! pressure thickness and tracers by calling incremental remapping routines.
! ------------------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm, itdm, jtdm
  use mod_types,     only: r8
  use mod_config,    only: expcnf
  use mod_constants, only: grav, epsilp, epsilz, onem
  use mod_time,      only: nstep, delt1, dlt
  use mod_xc,        only: xchalt, xcbcst, xcaput, xcstop, xctilr, &
                           mnproc, lp, ii, jj, kk, isp, ifp, ilp, &
                           isu, ifu, ilu, isv, ifv, ilv, isq, ifq, ilq, &
                           i0, j0, ip, iu, iv, iq, halo_ps, nbdy, nreg
  use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml, &
                           vcoord_cntiso_hybrid, sigref_spec, sigmar
  use mod_ale_regrid_remap, only: regrid_method_tag, regrid_method_direct, &
                           ale_regrid_remap
  use mod_grid,      only: scuy, scvx, scuyi, scvxi, plon, plat, depths, &
                           corioq
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, p, pu, &
                           pv, phi, ubflxs, vbflxs, ub, vb, pb, pbu, &
                           pbv, ubflxs_p, vbflxs_p, pb_p, pbu_p, pbv_p, &
                           ubcors_p, vbcors_p, kfpla
  use mod_pgforc,    only: pgfx, pgfy, pgfxm, pgfym, &
                           xixp, xixm, xiyp, xiym, pgforc
  use mod_barotp,    only: ubflx, vbflx, pb_mn, ubflx_mn, vbflx_mn, &
                           pvtrop
  use mod_tmsmt,     only: dpold
  use mod_temmin,    only: settemmin
  use mod_cesm,      only: inicon_cesm
  use mod_fuk95,     only: inicon_fuk95
  use mod_channel,   only: inicon_channel
  use mod_eos,       only: rho, sig, sofsig, delphi
  use mod_swtfrz,    only: swtfrz
  use mod_pointtest, only: itest, jtest, ptest
  use mod_checksum,  only: csdiag, chksummsk
  use mod_inicon_ben02, only: inicon_ben02
  use mod_utility,   only: fnmlen
  use mod_fill_global, only: fill_global
  use mod_hor3map,   only: recon_grd_struct, recon_src_struct, remap_struct, &
                           hor3map_ppm, hor3map_monotonic, &
                           hor3map_non_oscillatory, &
                           hor3map_non_oscillatory_posdef, &
                           initialize_rcgs, initialize_rcss, initialize_rms, &
                           prepare_reconstruction, reconstruct, &
                           prepare_remapping, remap, &
                           hor3map_noerr, hor3map_errstr
  use gsw_mod_toolbox, only: gsw_p_from_z, gsw_sa_from_sp, gsw_pt0_from_t
  use netcdf

  implicit none

  private

  ! Arrays and variables that are either allocated and specified in the NUOPC
  ! interface or internally in this module.
  real(r8), allocatable, dimension(:,:,:) :: t_woa, s_woa
  real(r8), allocatable, dimension(:,:) :: depth_bnds_woa
  real(r8), allocatable, dimension(:) :: depth_woa
  real(r8) :: t_woa_fval, s_woa_fval
  integer :: kdm_woa

  ! Variables to be set in namelist:
  character(len = fnmlen) :: &
    icfile                       ! Name of file containing initial conditions,
                                 ! that is either a valid restart file or file
                                 ! with climatological based initial conditions.
  logical :: &
    woa_nuopc_provided = .false. ! If true, WOA climatology has been provided
                                 ! via NUOPC.

  public :: icfile, woa_nuopc_provided, t_woa, s_woa, &
            depth_bnds_woa, depth_woa, t_woa_fval, s_woa_fval, kdm_woa, &
            inicon

contains

  ! ----------------------------------------------------------------------------
  ! Private procedures.
  ! ----------------------------------------------------------------------------

  function getpl(th,s,phiu,phil,pup) result(plo)
  ! ----------------------------------------------------------------------------
  ! Get lower pressure interface of a layer knowing the temperature, salinity of
  ! the layer and the geopotential at upper and lower interface.
  ! ----------------------------------------------------------------------------

    ! Arguments
    real(r8), intent(in) :: &
         th, &   ! Layer potential temperature [deg C].
         s, &    ! Layer salinity [g kg-1].
         phiu, & ! Geopotential at upper interface [m2 s-2].
         phil, & ! Geopotential at lower interface [m2 s-2].
         pup     ! Pressure at upper interface [kg m-1 s-2].

    ! Function output
    real(r8) :: plo ! Pressure at lower interface [kg m-1 s-2].

    ! Local variables
    real(r8) :: q,dphi,alpu,alpl

    ! first guess on pressure interface
    plo = pup-rho(pup,th,s)*(phil-phiu)

    ! improve the accuracy of the pressure interface by an
    ! iterative procedure
    q = 1._r8
    do while (abs(q) > 1.e-5_r8)
      call delphi(pup,plo,th,s,dphi,alpu,alpl)
      q = (phil-phiu-dphi)/alpl
      plo = plo-q
    end do

  end function getpl

  subroutine read_regridded_woa
  ! ----------------------------------------------------------------------------
  ! Read WOA climatology horizontally regridded onto the model grid.
  ! ----------------------------------------------------------------------------

    ! Local variables.
    real(r8), dimension(itdm,jtdm) :: tmp2d
    integer, dimension(3) :: istart, icount
    integer :: errstat, ncid, dimid, varid, i, j, k

    if (mnproc == 1) then

      write(lp,*) 'reading WOA climatology from '//trim(icfile)

      ! Open netCDF file.
      errstat = nf90_open(icfile, nf90_nowrite, ncid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_open: '//trim(icfile)//': '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif

      ! Check dimensions.
      errstat = nf90_inq_dimid(ncid, 'x', dimid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_dimid: x: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_inquire_dimension(ncid, dimid, len = i)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inquire_dimension: x: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_inq_dimid(ncid, 'y', dimid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_dimid: y: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_inquire_dimension(ncid, dimid, len = j)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inquire_dimension: y: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      if ((nreg == 2 .and. (i /= itdm .or. j /= jtdm-1)) .or. &
          (nreg /= 2 .and. (i /= itdm .or. j /= jtdm  ))) then
        write (lp,*) 'wrong dimensions in '//trim(icfile)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif

      ! Get depth dimension.
      errstat = nf90_inq_dimid(ncid, 'depth', dimid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_dimid: depth: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_inquire_dimension(ncid, dimid, len = kdm_woa)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inquire_dimension: depth: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif

      ! Set start and count vectors for reading horizontal slices of data.
      istart(1) = 1
      istart(2) = 1
      icount(1) = itdm
      if (nreg == 2) then
        icount(2) = jtdm - 1
      else
        icount(2) = jtdm
      endif
      icount(3) = 1

    endif

    ! Allocate arrays to receive data from file.
    call xcbcst(kdm_woa)
    allocate(t_woa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm_woa), &
             s_woa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm_woa), &
             depth_woa(kdm_woa), &
             depth_bnds_woa(2,kdm_woa), &
             stat = errstat)
    if (errstat /= 0) then
      write(lp,*) 'Failed to allocate WOA arrays!'
      call xchalt('(read_regridded_woa)')
             stop '(read_regridded_woa)'
    endif

    ! Read in situ temperature climatology.

    if (mnproc == 1) then
      errstat = nf90_inq_varid(ncid, 't_an', varid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_varid: t_an: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_get_att(ncid, varid, '_FillValue', t_woa_fval)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_get_att: t_an: _FillValue: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
    endif
    call xcbcst(t_woa_fval)

    do k = 1, kdm_woa
      if (mnproc == 1) then
        istart(3) = k
        errstat = nf90_get_var(ncid, varid, tmp2d, istart, icount)
        if (errstat /= nf90_noerr) then
           write(lp,*) 'nf90_get_var: t_an: '//nf90_strerror(errstat)
           call xchalt('(read_regridded_woa)')
                  stop '(read_regridded_woa)'
        endif
      endif
      call xcaput(tmp2d, t_woa(1-nbdy,1-nbdy,k), 1)
    enddo

    ! Read salinity climatology.

    if (mnproc == 1) then
      errstat = nf90_inq_varid(ncid, 's_an', varid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_varid: s_an: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_get_att(ncid, varid, '_FillValue', s_woa_fval)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_get_att: s_an: _FillValue: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
    endif
    call xcbcst(s_woa_fval)

    do k = 1, kdm_woa
      if (mnproc == 1) then
        istart(3) = k
        errstat = nf90_get_var(ncid, varid, tmp2d, istart, icount)
        if (errstat /= nf90_noerr) then
          write(lp,*) 'nf90_get_var: s_an: '//nf90_strerror(errstat)
          call xchalt('(read_regridded_woa)')
                 stop '(read_regridded_woa)'
        endif
      endif
      call xcaput(tmp2d, s_woa(1-nbdy,1-nbdy,k), 1)
    enddo

    ! Read depth and depth_bnds arrays.
    if (mnproc == 1) then
      errstat = nf90_inq_varid(ncid, 'depth', varid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_varid: depth: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_get_var(ncid, varid, depth_woa)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_get_var: depth: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_inq_varid(ncid, 'depth_bnds', varid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_inq_varid: depth_bnds: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
      errstat = nf90_get_var(ncid, varid, depth_bnds_woa)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_get_var: depth_bnds: '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
    endif
    call xcbcst(depth_woa)
    call xcbcst(depth_bnds_woa(1,:))
    call xcbcst(depth_bnds_woa(2,:))

    ! Close file.
    if (mnproc == 1) then
      errstat = nf90_close(ncid)
      if (errstat /= nf90_noerr) then
        write(lp,*) 'nf90_close: '//trim(icfile)//': '//nf90_strerror(errstat)
        call xchalt('(read_regridded_woa)')
               stop '(read_regridded_woa)'
      endif
    endif

    if (nreg == 2) then
      call xctilr(t_woa, 1, kdm_woa, 0,0, halo_ps)
      call xctilr(s_woa, 1, kdm_woa, 0,0, halo_ps)
    endif

  end subroutine read_regridded_woa

  subroutine inicon_woa_file
  ! ----------------------------------------------------------------------------
  ! Create initial interface geopotential and layer potential temperature and
  ! salinity directly from World Ocean Atlas (WOA) climatology at standard depth
  ! levels.
  ! ----------------------------------------------------------------------------

    ! Local variables.
    real(r8), allocatable, dimension(:) :: z_src_ref, z_src, sp_src, pt_src
    real(r8), dimension(kdm+1) :: z_dst_ref, z_dst
    real(r8), dimension(kdm) :: pt_dst, sp_dst
    real(r8) :: t_woa_mval, s_woa_mval, rk_woa, dk_woa, p_midlev, sa
    type(recon_grd_struct) :: rcgs
    type(recon_src_struct) :: pt_rcss, sp_rcss
    type(remap_struct) :: rms
    integer :: errstat, k_woa, i, j, k, l
    logical :: filling_failed

    if (trim(sigref_spec) /= 'namelist') then
      if (mnproc == 1) &
        write(lp,*) 'Initial conditions from WOA require sigref_spec == ''namelist''!'
      call xcstop('(inicon_woa_file)')
             stop '(inicon_woa_file)'
    endif

    if (vcoord_tag == vcoord_isopyc_bulkml) then
      if (mnproc == 1) &
        write(lp,*) 'Initial conditions from WOA is not supported for vcoord_type = ''isopyc_bulkml''!'
      call xcstop('(inicon_woa_file)')
             stop '(inicon_woa_file)'
    endif

    ! If WOA climatology is not provided via NUOPC, read climatology that has
    ! been horizontally regridded onto the model grid.
    if (.not. woa_nuopc_provided) call read_regridded_woa

    ! Allocate 1D WOA arrays.
    allocate(z_src_ref(kdm_woa+1), &
             z_src(kdm_woa+1), &
             sp_src(kdm_woa), &
             pt_src(kdm_woa), &
             stat = errstat)
    if (errstat /= 0) then
      write(lp,*) 'Failed to allocate WOA 1D arrays!'
      call xchalt('(inicon_woa_file)')
             stop '(inicon_woa_file)'
    endif

    ! Fill missing values within the model grid depth range by laterally
    ! extrapolating values from neighbouring points.
    t_woa_mval = 2._r8*t_woa_fval
    s_woa_mval = 2._r8*s_woa_fval
    do k = 1, kdm_woa
      do j = 1, jj
        do i = 1, ii
          if (ip(i,j) == 0 .or. depths(i,j) < depth_bnds_woa(1,k)) then
            t_woa(i,j,k) = t_woa_mval
            s_woa(i,j,k) = s_woa_mval
          endif
        enddo
      enddo
      if (mnproc == 1) then
         write(lp,*)'calling fill_global at k = ',k
      end if
      call fill_global(t_woa_mval, t_woa_fval, halo_ps, &
                       t_woa(1-nbdy,1-nbdy,k))
      call fill_global(s_woa_mval, s_woa_fval, halo_ps, &
                       s_woa(1-nbdy,1-nbdy,k))
    enddo

    ! Check that there are no missing values in the shallowest depth level.
    filling_failed = .false.
    do j = 1, jj
      do l = 1, isp(j)
      do i = ifp(j,l), ilp(j,l)
        if (t_woa(i,j,1) == t_woa_fval .or. &
            s_woa(i,j,1) == s_woa_fval) filling_failed = .true.
      enddo
      enddo
    enddo
    call xcbcst(filling_failed)
    if (filling_failed) then
      if (mnproc == 1) &
        write(lp,*) 'Failed to fill missing values in shallowest WOA depth level!'
      call xcstop('(inicon_woa_file)')
             stop '(inicon_woa_file)'
    endif

    ! Create source and destination depth interface arrays, where the latter has
    ! length matching the number of model interfaces. The source depths are
    ! mapped to destination depths in index space.
    z_src_ref(1) = - depth_bnds_woa(1,1)
    do k = 1, kdm_woa
      z_src_ref(k+1) = - depth_bnds_woa(2,k)
    enddo
    z_dst_ref(1) = z_src_ref(1)
    do k = 2, kk
      rk_woa = real(kdm_woa*(k - 1), r8)/real(kk, r8) + 1._r8
      k_woa = int(rk_woa)
      dk_woa = rk_woa - k_woa
      z_dst_ref(k) = z_src_ref(k_woa  )*(1._r8 - dk_woa) &
                   + z_src_ref(k_woa+1)*dk_woa
    enddo
    z_dst_ref(kk+1) = z_src_ref(kdm_woa+1)

    ! Prepare vertical remapping.

    rcgs%n_src = kdm_woa
    rcgs%method = hor3map_ppm

    pt_rcss%limiting = hor3map_non_oscillatory
    pt_rcss%pc_left_bndr = .true.
    pt_rcss%pc_right_bndr = .true.

    sp_rcss%limiting = hor3map_non_oscillatory_posdef
    sp_rcss%pc_left_bndr = .true.
    sp_rcss%pc_right_bndr = .true.

    ! Process all columns to convert from in situ to potential temperature and
    ! vertically remap to array with model dimensions.

    do j = 1, jj
      do l = 1, isp(j)
      do i = ifp(j,l), ilp(j,l)

        ! Create source arrays of practical salinity (sp_src), potential
        ! temperature (pt_src) and source interface depths bounded by model
        ! depth (z_src). Fill remaining missing values from above.
        do k = 1, kdm_woa
          if (t_woa(i,j,k) /= t_woa_fval .and. &
              t_woa(i,j,k) /= t_woa_mval .and. &
              s_woa(i,j,k) /= s_woa_fval .and. &
              s_woa(i,j,k) /= s_woa_mval) then
            sp_src(k) = s_woa(i,j,k)
            p_midlev = gsw_p_from_z(- depth_woa(k), plat(i,j))
            sa = gsw_sa_from_sp(sp_src(k), p_midlev, plon(i,j), plat(i,j))
            pt_src(k) = gsw_pt0_from_t(sa, t_woa(i,j,k), p_midlev)
          else
            sp_src(k) = sp_src(k-1)
            pt_src(k) = pt_src(k-1)
          endif
          z_src(k) = max(z_src_ref(k), - depths(i,j))
        enddo
        z_src(kdm_woa+1) = - depths(i,j)

        ! Bound destination interface depths by model depth.
        do k = 1, kk
          z_dst(k) = max(z_dst_ref(k), - depths(i,j))
        enddo
         z_dst(kk+1) = - depths(i,j)


        ! Prepare vertical reconstruction and remapping.

        errstat = prepare_reconstruction(rcgs, z_src)
        if (errstat /= hor3map_noerr) then
          write(lp,*) trim(hor3map_errstr(errstat))
          call xchalt('(inicon_woa_file)')
                 stop '(inicon_woa_file)'
        endif

        errstat = prepare_remapping(rcgs, rms, z_dst)
        if (errstat /= hor3map_noerr) then
          write(lp,*) trim(hor3map_errstr(errstat))
          call xchalt('(inicon_woa_file)')
                 stop '(inicon_woa_file)'
        endif

        ! Reconstruct and remap potential temperature.

        errstat = reconstruct(rcgs, pt_rcss, pt_src)
        if (errstat /= hor3map_noerr) then
          write(*,*) trim(hor3map_errstr(errstat))
          call xchalt('(inicon_woa_file)')
                 stop '(inicon_woa_file)'
        endif

        errstat = remap(pt_rcss, rms, pt_dst)
        if (errstat /= hor3map_noerr) then
          write(*,*) trim(hor3map_errstr(errstat))
          call xchalt('(inicon_woa_file)')
                 stop '(inicon_woa_file)'
        endif

        ! Reconstruct and remap practical salinity.

        errstat = reconstruct(rcgs, sp_rcss, sp_src)
        if (errstat /= hor3map_noerr) then
          write(*,*) trim(hor3map_errstr(errstat))
          call xchalt('(inicon_woa_file)')
                 stop '(inicon_woa_file)'
        endif

        errstat = remap(sp_rcss, rms, sp_dst)
        if (errstat /= hor3map_noerr) then
          write(*,*) trim(hor3map_errstr(errstat))
          call xchalt('(inicon_woa_file)')
                 stop '(inicon_woa_file)'
        endif

        ! Copy potential temperature and practical salinity to model arrays and
        ! compute geopotential at layer interfaces.
        do k = 1, kk
          temp(i,j,k) = pt_dst(k)
          saln(i,j,k) = sp_dst(k)
          phi(i,j,k) = grav*z_dst(k)
        enddo
        phi(i,j,kk+1) = grav*z_dst(kk+1)

      enddo
      enddo
    enddo

    ! Deallocate arrays.
    deallocate(t_woa, s_woa, depth_woa, depth_bnds_woa, &
               z_src_ref, z_src, sp_src, pt_src)

  end subroutine inicon_woa_file

  subroutine inicon_layer_file
  ! ----------------------------------------------------------------------------
  ! Read initial conditions on model horizontal grid from file containing layer
  ! variables of thickness, potential temperature, salinity and reference
  ! potential density.
  ! ----------------------------------------------------------------------------

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dz
    real, dimension(itdm,jtdm) :: tmp2d
    real :: dsig,a0,a1,a2
    integer, dimension(3) :: start,count
    integer :: i,j,kdmic,k,l,status,ncid,dimid,varid,kb

    if (mnproc == 1) then

      write(lp,*) 'reading layer initial condition from '//trim(icfile)

      ! Open netcdf file.
      status = nf90_open(icfile,nf90_nowrite,ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_open: ',trim(icfile),': ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if

      ! Check dimensions.
      status = nf90_inq_dimid(ncid,'x',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = i)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: x: ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
      status = nf90_inq_dimid(ncid,'y',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = j)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: y: ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
      status = nf90_inq_dimid(ncid,'z',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: z: ',nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = kdmic)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: z: ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
      if (i /= itdm .or. j /= jtdm .or. &
          (kdmic /= kdm .and. vcoord_tag == vcoord_isopyc_bulkml) .or. &
          (kdmic >  kdm .and. vcoord_tag == vcoord_cntiso_hybrid .and. &
           trim(sigref_spec) == 'inicon')) then
        write (lp,*) 'wrong dimensions in '//trim(icfile)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if

    end if

    call xcbcst(kdmic)

    start(1) = 1
    start(2) = 1
    count(1) = itdm
    count(2) = jtdm
    count(3) = 1

    ! Read reference potential density.
    if ( vcoord_tag == vcoord_isopyc_bulkml .or. &
        (vcoord_tag == vcoord_cntiso_hybrid .and. &
         trim(sigref_spec) == 'inicon')) then
      if (mnproc == 1) then
        status = nf90_inq_varid(ncid,'sigma',varid)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: sigma: ', &
               nf90_strerror(status)
          call xchalt('(inicon_layer_file)')
          stop '(inicon_layer_file)'
        end if
      end if
      do k = 1,kdmic
        if (mnproc == 1) then
          start(3) = k
          status = nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: sigma: ', &
                 nf90_strerror(status)
            call xchalt('(inicon_layer_file)')
            stop '(inicon_layer_file)'
          end if
        end if
        call xcaput(tmp2d,sigmar(1-nbdy,1-nbdy,k),1)
      end do
    end if

    ! Read potential temperature.
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'temp',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: temp: ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: temp: ', &
               nf90_strerror(status)
          call xchalt('(inicon_layer_file)')
          stop '(inicon_layer_file)'
        end if
      end if
      call xcaput(tmp2d,temp(1-nbdy,1-nbdy,k),1)
    end do

    ! Read salinity.
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'saln',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: saln: ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: saln: ', &
               nf90_strerror(status)
          call xchalt('(inicon_layer_file)')
          stop '(inicon_layer_file)'
        end if
      end if
      call xcaput(tmp2d,saln(1-nbdy,1-nbdy,k),1)
    end do

    ! Read layer thickness.
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'dz',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: dz: ',nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: dz: ',nf90_strerror(status)
          call xchalt('(inicon_layer_file)')
          stop '(inicon_layer_file)'
        end if
      end if
      call xcaput(tmp2d,dz(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_close(ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_close: ',trim(icfile),': ', &
             nf90_strerror(status)
        call xchalt('(inicon_layer_file)')
        stop '(inicon_layer_file)'
      end if
    end if

    if (vcoord_tag == vcoord_cntiso_hybrid .and. &
        trim(sigref_spec) == 'inicon') then
      !$omp parallel do private(l,i,k,kb,dsig,a0,a1,a2)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (mnproc == ptest.and.i == itest.and.j == jtest) then
              write(lp,*) &
                   'Layer reference potential densities from initial condition:'
              do k = 1,kdmic
                write(lp,*) k,sigmar(i,j,k)
              end do
            end if
            do k = kdmic,2,-1
              sigmar(i,j,k+kk-kdmic) = .5_r8*(sigmar(i,j,k-1) &
                                             +sigmar(i,j,k  ))
            end do
            kb = kk-kdmic+2
            dsig = sigmar(i,j,kb+1)-sigmar(i,j,kb)
            a0 = 1./(kb-1)**2
            a1 = (dsig+(2.*sigmar(i,j,kb)-dsig*kb)*kb)*a0
            a2 = (dsig*(kb-1)-sigmar(i,j,kb))*a0
            a0 = -a1-a2
            do k = 1,kb-1
              sigmar(i,j,k) = a0+(a1+a2*k)*k
            end do
            if (mnproc == ptest.and.i == itest.and.j == jtest) then
              write(lp,*) &
                   'Generated interface reference potential densities:'
              do k = 1,kk
                write(lp,*) k,sigmar(i,j,k)
              end do
            end if
          end do
        end do
      end do
      !$omp end parallel do

    end if

    ! Construct interface depths [m] from layer thicknesses [m].
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          ! z(i,j,1)=z(i,j,1)
          z(i,j,1) = 0.
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,kdmic
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            z(i,j,k+1) = min(depths(i,j), &
                             z(i,j,k)+dz(i,j,k))
          end do
        end do
      end do
      do k = kdmic+1,kk
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            z(i,j,k+1) = z(i,j,kdmic+1)
          end do
        end do
      end do
      do k = 2,kk
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (z(i,j,kk+1)-z(i,j,k) < 1.e-6) then
              z(i,j,k) = depths(i,j)
            end if
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          z(i,j,kk+1) = depths(i,j)
        end do
      end do
    end do
    !$omp end parallel do

    ! Compute layer interface geopotential.
    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,kk+1
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            phi(i,j,k) = -grav*z(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine inicon_layer_file

  subroutine inicon_file
  ! ----------------------------------------------------------------------------
  ! Obtain initial conditions from file.
  ! ----------------------------------------------------------------------------

    ! Local variables.
    integer :: errstat, ncid, varid
    logical :: woa_icfile, layer_icfile

    if (woa_nuopc_provided) then
       woa_icfile = .true.
    end if

    if (.not. woa_nuopc_provided) then
      ! Check type of initial condition file.
      woa_icfile = .false.
      layer_icfile = .false.
      if (mnproc == 1) then
        errstat = nf90_open(icfile, nf90_nowrite, ncid)
        if (errstat == nf90_noerr) then
          woa_icfile = &
            (nf90_inq_varid(ncid, 'depth'     , varid) == nf90_noerr) .and. &
            (nf90_inq_varid(ncid, 'depth_bnds', varid) == nf90_noerr) .and. &
            (nf90_inq_varid(ncid, 't_an'      , varid) == nf90_noerr) .and. &
            (nf90_inq_varid(ncid, 's_an'      , varid) == nf90_noerr)
          layer_icfile = &
            (nf90_inq_varid(ncid, 'sigma', varid) == nf90_noerr) .and. &
            (nf90_inq_varid(ncid, 'temp' , varid) == nf90_noerr) .and. &
            (nf90_inq_varid(ncid, 'saln' , varid) == nf90_noerr) .and. &
            (nf90_inq_varid(ncid, 'dz'   , varid) == nf90_noerr)
          errstat = nf90_close(ncid)
          if (errstat /= nf90_noerr) then
            write(lp,*) 'nf90_close: '//trim(icfile)//': '// &
                        nf90_strerror(errstat)
            call xchalt('(inicon_file)')
                   stop '(inicon_file)'
          endif
        else
          write(lp,*) 'nf90_open: '//trim(icfile)//': '//nf90_strerror(errstat)
          call xchalt('(inicon_file)')
                 stop '(inicon_file)'
        endif
      endif
      call xcbcst(woa_icfile)
      call xcbcst(layer_icfile)
      if (.not. (woa_icfile .or. layer_icfile)) then
        if (mnproc == 1) then
          write(lp,*) 'inicon_file: cannot recognize type of initial '// &
                      'condition file '//trim(icfile)
        endif
        call xcstop('(inicon_file)')
               stop '(inicon_file)'
      endif
    endif

    ! Obtain initial conditions based on the source type.
    if     (woa_nuopc_provided .or. woa_icfile) then
      call inicon_woa_file
    elseif (layer_icfile) then
      call inicon_layer_file
    else
      if (mnproc == 1) then
        write(lp,*) 'inicon_file: no valid source of initial conditions'
      endif
      call xcstop('(inicon_file)')
             stop '(inicon_file)'
    endif

  end subroutine inicon_file

  ! ----------------------------------------------------------------------------
  ! Public procedures.
  ! ----------------------------------------------------------------------------

  subroutine inicon
  ! ----------------------------------------------------------------------------
  ! Define initial conditions.
  ! ----------------------------------------------------------------------------

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tfrz
    integer :: i,j,k,l,regrid_method_tag_orig
    real :: q,tsfac,dps

    ! --------------------------------------------------------------------------
    ! Define layer interface heights and layer temperature and salinity.
    ! --------------------------------------------------------------------------

    select case (trim(expcnf))
      case ('cesm', 'ben02clim', 'ben02syn', 'noforcing', 'single_column')
        call inicon_file
      case ('fuk95')
        call inicon_fuk95
      case ('channel')
        call inicon_channel
      case ('isomip1')
        ! call inicon_isomip1
      case ('isomip2')
        ! call inicon_isomip2
      case default
        if (mnproc == 1) then
          write (lp,'(3a)') ' inicon: expcnf = ', trim(expcnf), &
               ' is unsupported!'
        end if
        call xcstop('(inicon)')
        stop '(inicon)'
    end select

    ! --------------------------------------------------------------------------
    ! Set minimum physical temperature for each isopycnic layer.
    ! --------------------------------------------------------------------------

    call settemmin

    ! --------------------------------------------------------------------------
    ! Initialize configuration specific variables.
    ! --------------------------------------------------------------------------

    select case (trim(expcnf))
    case ('cesm')
      call inicon_cesm
    case ('ben02clim', 'ben02syn', 'single_column')
      call inicon_ben02
    end select

    ! --------------------------------------------------------------------------
    ! Make sure layer temperature is greater than the lower physical bound and
    ! make temperature, salinity, and potential density variables consistent.
    ! --------------------------------------------------------------------------

    select case (vcoord_tag)

      case (vcoord_isopyc_bulkml)

        do k = 1,2
          tfrz(1:ii,1:jj) = swtfrz(p(1:ii,1:jj,1),saln(1:ii,1:jj,k))
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                temp(i,j,k) = max(tfrz(i,j),temp(i,j,k))
                sigma(i,j,k) = sig(temp(i,j,k),saln(i,j,k))
              end do
            end do
          end do
          !$omp end parallel do
        end do
        do k = 3,kk
          tfrz(1:ii,1:jj) = swtfrz(p(1:ii,1:jj,1),saln(1:ii,1:jj,k))
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                temp(i,j,k) = max(tfrz(i,j),temp(i,j,k))
                saln(i,j,k) = sofsig(sigmar(i,j,k),temp(i,j,k))
                sigma(i,j,k) = sig(temp(i,j,k),saln(i,j,k))
              end do
            end do
          end do
          !$omp end parallel do
        end do

      case default

        do k = 1,kk
          tfrz(1:ii,1:jj) = swtfrz(p(1:ii,1:jj,1),saln(1:ii,1:jj,k))
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                temp(i,j,k) = max(tfrz(i,j),temp(i,j,k))
                sigma(i,j,k) = sig(temp(i,j,k),saln(i,j,k))
              end do
            end do
          end do
          !$omp end parallel do
        end do

    end select

    if (mnproc == ptest) then
      write (lp,'('' sigmar(k)    :'',7f9.5/(15x,7f9.5))') &
           (sigmar(itest,jtest,k),k = 1,kk)
    end if

    ! --------------------------------------------------------------------------
    ! Find layer interface pressure.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(l,i,k)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          p(i,j,1) = getpl(temp(i,j,1),saln(i,j,1),0.,phi(i,j,1),0.)
        end do
      end do
      do k = 1,kk
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1) = getpl(temp(i,j,k), saln(i,j,k), &
                               phi(i,j,k), phi(i,j,k+1), p(i,j,k))
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(p, 1,kk+1, 2,2, halo_ps)
    call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 1,1, halo_ps)

    ! --------------------------------------------------------------------------
    ! Set layer thickness and bottom pressure.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(k,l,i)
    do j = 0,jj+1
      do k = 1,kk
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,k) = p(i,j,k+1)-p(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,l,i)
    do j = 0,jj+1
      do k = 1,kk
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i)
    do j = 0,jj+1
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          pb(i,j,1) = p(i,j,kk+1)
          pb(i,j,2) = pb(i,j,1)
          pb_mn(i,j,1) = pb(i,j,1)
          pb_mn(i,j,2) = pb(i,j,1)
          pb_p(i,j) = pb(i,j,1)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          pbu(i,j,1) = min(pb(i,j,1),pb(i-1,j,1))
          pbu(i,j,2) = pbu(i,j,1)
          pbu_p(i,j) = pbu(i,j,1)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          pbv(i,j,1) = min(pb(i,j,1),pb(i,j-1,1))
          pbv(i,j,2) = pbv(i,j,1)
          pbv_p(i,j) = pbv(i,j,1)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,l,i,q)
    do j = -1,jj+2
      do k = 1,kk
        do l = 1,isu(j)
          do i = max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q = min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,k)= &
                 .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                    +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
            pu(i,j,k+1) = pu(i,j,k)+dpu(i,j,k)
          end do
        end do
        do l = 1,isv(j)
          do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q = min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,k)= &
                 .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                    +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
            pv(i,j,k+1) = pv(i,j,k)+dpv(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      regrid_method_tag_orig = regrid_method_tag
      regrid_method_tag = regrid_method_direct
      call ale_regrid_remap(2,1,kk,0,kk+1,1)
      regrid_method_tag = regrid_method_tag_orig
    end if
    call xctilr(temp, 1,kk, 1,1, halo_ps)
    call xctilr(saln, 1,kk, 1,1, halo_ps)

    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,kk
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            temp(i,j,k+kk) = temp(i,j,k)
            saln(i,j,k+kk) = saln(i,j,k)
            sigma(i,j,k+kk) = sigma(i,j,k)
            dp(i,j,k+kk) = dp(i,j,k)
          end do
        end do
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            dpu(i,j,k+kk) = dpu(i,j,k)
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            dpv(i,j,k+kk) = dpv(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Initialize potential vorticity of barotropic flow.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(l,i,q)
    do j = 0,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          q = 2./(pb_p(i,j)+pb_p(i-1,j))
          pvtrop(i,j  ,1) = corioq(i,j  )*q
          pvtrop(i,j+1,1) = corioq(i,j+1)*q
          pvtrop(i,j  ,2) = pvtrop(i,j  ,1)
          pvtrop(i,j+1,2) = pvtrop(i,j+1,1)
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i,q)
    do j = 1,jj
      do l = 1,isv(j)
        do i = max(0,ifv(j,l)),min(ii,ilv(j,l))
          q = 2./(pb_p(i,j)+pb_p(i,j-1))
          pvtrop(i  ,j,1) = corioq(i  ,j)*q
          pvtrop(i+1,j,1) = corioq(i+1,j)*q
          pvtrop(i  ,j,2) = pvtrop(i  ,j,1)
          pvtrop(i+1,j,2) = pvtrop(i+1,j,1)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isq(j)
        do i = max(1,ifq(j,l)),min(ii,ilq(j,l))
          pvtrop(i,j,1) = corioq(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  ) &
                                         +pb_p(i,j-1)+pb_p(i-1,j-1))
          pvtrop(i,j,2) = pvtrop(i,j,1)
        end do
      end do
    end do
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Separate baroclinic and barotropic velocity components.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(l,i,k)
    do j = 1,jj

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          ub(i,j,1) = 0.
        end do
        do k = 1,kk
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            ub(i,j,1) = ub(i,j,1)+u(i,j,k)*dpu(i,j,k)
          end do
        end do
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          ub(i,j,1) = ub(i,j,1)/pbu_p(i,j)
          ub(i,j,2) = ub(i,j,1)
        end do
      end do

      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          vb(i,j,1) = 0.
        end do
        do k = 1,kk
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vb(i,j,1) = vb(i,j,1)+v(i,j,k)*dpv(i,j,k)
          end do
        end do
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          vb(i,j,1) = vb(i,j,1)/pbv_p(i,j)
          vb(i,j,2) = vb(i,j,1)
        end do
      end do
    end do
    !$omp end parallel do

    do k = 1,kk
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            u(i,j,k) = u(i,j,k)-ub(i,j,1)
            u(i,j,k+kk) = u(i,j,k)
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            v(i,j,k) = v(i,j,k)-vb(i,j,1)
            v(i,j,k+kk) = v(i,j,k)
          end do
        end do
      end do
      !$omp end parallel do
    end do

    tsfac = delt1/dlt

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          ubflx_mn(i,j,1) = ub(i,j,1)*pbu_p(i,j)*scuy(i,j)
          ubflx_mn(i,j,2) = ubflx_mn(i,j,1)
          ubflx(i,j,1) = ubflx_mn(i,j,1)
          ubflx(i,j,2) = ubflx_mn(i,j,1)
          ubflxs(i,j,1) = ubflx_mn(i,j,1)*tsfac
          ubflxs(i,j,2) = ubflxs(i,j,1)
          ubflxs(i,j,3) = ubflxs(i,j,1)
          ubflxs_p(i,j,1) = ubflxs(i,j,1)
          ubflxs_p(i,j,2) = ubflxs(i,j,1)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          vbflx_mn(i,j,1) = vb(i,j,1)*pbv_p(i,j)*scvx(i,j)
          vbflx_mn(i,j,2) = vbflx_mn(i,j,1)
          vbflx_mn(i,j,2) = vbflx_mn(i,j,1)
          vbflx(i,j,1) = vbflx_mn(i,j,1)
          vbflx(i,j,2) = vbflx_mn(i,j,1)
          vbflxs(i,j,1) = vbflx_mn(i,j,1)*tsfac
          vbflxs(i,j,2) = vbflxs(i,j,1)
          vbflxs(i,j,3) = vbflxs(i,j,1)
          vbflxs_p(i,j,1) = vbflxs(i,j,1)
          vbflxs_p(i,j,2) = vbflxs(i,j,1)
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          ubcors_p(i,j)= (vbflx_mn(i  ,j  ,1)*scvxi(i  ,j  ) &
                         +vbflx_mn(i  ,j+1,1)*scvxi(i  ,j+1) &
                         +vbflx_mn(i-1,j  ,1)*scvxi(i-1,j  ) &
                         +vbflx_mn(i-1,j+1,1)*scvxi(i-1,j+1)) &
                         *(pvtrop(i,j,1)+pvtrop(i,j+1,1))*.125*tsfac
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          ubcors_p(i,j) = -(ubflx_mn(i  ,j  ,1)*scuyi(i  ,j  ) &
                           +ubflx_mn(i+1,j  ,1)*scuyi(i+1,j  ) &
                           +ubflx_mn(i  ,j-1,1)*scuyi(i  ,j-1) &
                           +ubflx_mn(i+1,j-1,1)*scuyi(i+1,j-1)) &
                           *(pvtrop(i,j,1)+pvtrop(i+1,j,1))*.125*tsfac
        end do
      end do
    end do
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Initialize fields related to the pressure gradient force.
    ! --------------------------------------------------------------------------

    call pgforc(2,1,kk,0,kk+1,1)

    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,kk
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            pgfx(i,j,k+kk) = pgfx(i,j,k)
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            pgfy(i,j,k+kk) = pgfy(i,j,k)
          end do
        end do
      end do
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          pgfxm(i,j,2) = pgfxm(i,j,1)
          xixp(i,j,2) = xixp(i,j,1)
          xixm(i,j,2) = xixm(i,j,1)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          pgfym(i,j,2) = pgfym(i,j,1)
          xiyp(i,j,2) = xiyp(i,j,1)
          xiym(i,j,2) = xiym(i,j,1)
        end do
      end do
    end do
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Define first physical interior layer.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(l,i,k,dps)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          k = 3
          dps = 0.
          do while (dp(i,j,k) < epsilp)
            dps = dps+dp(i,j,k)
            dp(i,j,k) = 0.
            dp(i,j,k+kk) = 0.
            k = k+1
            if (k > kk) exit
          end do
          if (k > kk) then
            dp(i,j,2) = dp(i,j,2)+dps
            dp(i,j,2+kk) = dp(i,j,2)
          else
            dp(i,j,k) = dp(i,j,k)+dps
            dp(i,j,k+kk) = dp(i,j,k)
          end if
          kfpla(i,j,1) = k
          kfpla(i,j,2) = k
        end do
      end do
    end do
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Set other time level layer thicknesses.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,kk
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,k+kk) = dp(i,j,k)
            dpold(i,j,k) = dp(i,j,k)
            dpold(i,j,k+kk) = dp(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (mnproc == ptest) then
      i = itest
      j = jtest
      write (lp,103) nstep,i0+i,j0+j, &
           '  init.profile  temp    saln    dens   thkns    dpth', &
           (k,temp(i,j,k),saln(i,j,k), &
           sig(temp(i,j,k),saln(i,j,k)), &
           dp(i,j,k)/onem,p(i,j,k+1)/onem,k = 1,kk)
103   format (i9,2i5,a/(28x,i3,3f8.2,2f8.1))
    end if

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'inicon:'
      end if
      call chksummsk(p,ip,kk+1,'p')
      call chksummsk(dp,ip,2*kk,'dp')
      call chksummsk(temp,ip,2*kk,'temp')
      call chksummsk(saln,ip,2*kk,'saln')
      call chksummsk(sigma,ip,2*kk,'sigma')
      call chksummsk(pb,ip,3,'pb')
      call chksummsk(pbu,iu,2,'pbu')
      call chksummsk(pbv,iv,2,'pbv')
      call chksummsk(pvtrop,iq,2,'pvtrop')
      call chksummsk(pu,iu,kk+1,'pu')
      call chksummsk(pv,iv,kk+1,'pv')
    end if

  end subroutine inicon

end module mod_inicon
