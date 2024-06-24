! ------------------------------------------------------------------------------
! Copyright (C) 2008-2024 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin,
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

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to advection
  ! of layer pressure thickness and tracers by calling incremental
  ! remapping routines.
  ! ------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm, itdm, jtdm
  use mod_types,     only: r8
  use mod_config,    only: expcnf
  use mod_constants, only: g, epsilp, onem, &
                           L_mks2cgs, M_mks2cgs, P_mks2cgs
  use mod_time,      only: nstep, delt1, dlt
  use mod_xc,        only: xchalt, xcbcst, xcaput, xcstop, xctilr, &
                           mnproc, lp, ii, jj, kk, isp, ifp, ilp, &
                           isu, ifu, ilu, isv, ifv, ilv, isq, ifq, ilq, &
                           i0, j0, ip, iu, iv, iq, halo_ps, nbdy
  use mod_vcoord,    only: vcoord_type_tag, isopyc_bulkml, &
                           cntiso_hybrid, sigmar, &
                           cntiso_hybrid_regrid_direct_remap, &
                           remap_velocity
  use mod_grid,      only: scuy, scvx, scuyi, scvxi, depths, &
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
  use mod_fuk95,     only: ictsz_fuk95
  use mod_channel,   only: ictsz_channel
  use mod_eos,       only: rho, sig, sofsig, delphi
  use mod_swtfrz,    only: swtfrz
  use mod_pointtest, only: itest, jtest, ptest
  use mod_checksum,  only: csdiag, chksummsk
  use mod_inicon_ben02, only: inicon_ben02
  use netcdf

  implicit none
  private

  ! Variables to be set in namelist:
  character(len = 256), public :: &
       icfile ! Name of file containing initial conditions, that is
  ! either a valid restart file or 'inicon.nc' if
  ! climatological based initial conditions are desired.

  ! Public routines
  public :: inicon

  ! Private routines
  private :: getpl, ictsz_file

contains

  ! ------------------------------------------------------------------
  ! Private procedures.
  ! ------------------------------------------------------------------

  function getpl(th,s,phiu,phil,pup) result(plo)

    ! ------------------------------------------------------------------
    ! Get lower pressure interface of a layer knowing the temperature,
    ! salinity of the layer and the geopotential at upper and lower
    ! interface.
    ! ------------------------------------------------------------------

    ! Arguments
    real(r8), intent(in) :: &
         th    ! Layer potential temperature [deg C]. &
    real(r8), intent(in) :: &
         s     ! Layer salinity [g kg-1]. &
    real(r8), intent(in) :: &
         phiu  ! Geopotential at upper interface [cm2 s-2]. &
    real(r8), intent(in) :: &
         phil  ! Geopotential at lower interface [cm2 s-2]. &
    real(r8), intent(in) :: &
         pup   ! Pressure at upper interface [g cm-1 s-2].

    ! Function output
    real(r8) :: plo ! Pressure at lower interface [g cm-1 s-2].

    ! Local variables
    real(r8) :: q,dphi,alpu,alpl

    ! first guess on pressure interface
    plo = pup-rho(pup,th,s)*(phil-phiu)

    ! improve the accuracy of the pressure interface by an
    ! iterative procedure
    q = 1._r8
    do while (abs(q) > 1.e-5_r8*p_mks2cgs)
      call delphi(pup,plo,th,s,dphi,alpu,alpl)
      q = (phil-phiu-dphi)/alpl
      plo = plo-q
    end do

  end function getpl

  ! ------------------------------------------------------------------

  subroutine ictsz_file

    ! ------------------------------------------------------------------
    ! Read initial conditions from file to define layer temperature and
    ! salinity and geopotential at layer interfaces.
    ! ------------------------------------------------------------------

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dz
    real, dimension(itdm,jtdm) :: tmp2d
    real :: dsig,a0,a1,a2
    integer, dimension(3) :: start,count
    integer :: i,j,kdmic,k,l,status,ncid,dimid,varid,kb
    real :: im_mks2cgs

    iM_mks2cgs = 1.0 / M_mks2cgs

    if (mnproc == 1) then
      write (lp,'(2a)') ' reading initial condition from ', &
           trim(icfile)
      call flush(lp)

      ! Open netcdf file
      status = nf90_open(icfile,nf90_nowrite,ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_open: ',trim(icfile),': ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if

      ! Check dimensions
      status = nf90_inq_dimid(ncid,'x',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = i)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: x: ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
      status = nf90_inq_dimid(ncid,'y',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = j)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: y: ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
      status = nf90_inq_dimid(ncid,'z',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: z: ',nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = kdmic)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: z: ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
      if (i /= itdm.or.j /= jtdm.or. &
           (kdmic /= kdm.and.vcoord_type_tag /= cntiso_hybrid).or. &
           (kdmic > kdm.and.vcoord_type_tag == cntiso_hybrid)) then
        write (lp,*) 'wrong dimensions in '//trim(icfile)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if

    end if

    call xcbcst(kdmic)

    start(1) = 1
    start(2) = 1
    count(1) = itdm
    count(2) = jtdm
    count(3) = 1

    ! Read reference potential density
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'sigma',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: sigma: ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: sigma: ', &
               nf90_strerror(status)
          call xchalt('(ictsz_file)')
          stop '(ictsz_file)'
        end if
      end if
      call xcaput(tmp2d,sigmar(1-nbdy,1-nbdy,k),1)
    end do

    ! Read potential temperature
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'temp',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: temp: ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: temp: ', &
               nf90_strerror(status)
          call xchalt('(ictsz_file)')
          stop '(ictsz_file)'
        end if
      end if
      call xcaput(tmp2d,temp(1-nbdy,1-nbdy,k),1)
    end do

    ! Read salinity
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'saln',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: saln: ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: saln: ', &
               nf90_strerror(status)
          call xchalt('(ictsz_file)')
          stop '(ictsz_file)'
        end if
      end if
      call xcaput(tmp2d,saln(1-nbdy,1-nbdy,k),1)
    end do

    ! Read layer thickness
    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'dz',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: dz: ',nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
    end if
    do k = 1,kdmic
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmp2d,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: dz: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
          stop '(ictsz_file)'
        end if
      end if
      call xcaput(tmp2d,dz(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_close(ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_close: ',trim(icfile),': ', &
             nf90_strerror(status)
        call xchalt('(ictsz_file)')
        stop '(ictsz_file)'
      end if
    end if

    if (vcoord_type_tag == cntiso_hybrid) then
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

    ! Construct interface depths [cm] from layer thicknesses [m] and
    ! convert unit of reference potential density from [kg/m^3] to
    ! [g/cm^3]
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          ! z(i,j,1)=z(i,j,1)*L_mks2cgs
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
            z(i,j,k+1) = min(depths(i,j)*L_mks2cgs, &
                             z(i,j,k)+dz(i,j,k)*L_mks2cgs)
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
            if (z(i,j,kk+1)-z(i,j,k) < 1.e-6*l_mks2cgs) then
              z(i,j,k) = depths(i,j)*L_mks2cgs
            end if
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          z(i,j,kk+1) = depths(i,j)*L_mks2cgs
        end do
      end do
      do k = 1,kk
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            sigmar(i,j,k) = sigmar(i,j,k)*iM_mks2cgs
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! compute layer interface geopotential
    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,kk+1
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            phi(i,j,k) = -g*z(i,j,k)
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine ictsz_file

  ! ------------------------------------------------------------------
  ! Public procedures.
  ! ------------------------------------------------------------------

  subroutine inicon()
  ! ------------------------------------------------------------------
  ! Define initial conditions
  ! ------------------------------------------------------------------

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tfrz
    integer :: i,j,k,l
    real :: q,tsfac,dps

    ! ------------------------------------------------------------------
    ! Define layer interface heights and layer temperature and salinity
    ! ------------------------------------------------------------------

    select case (trim(expcnf))
      case ('cesm', 'ben02clim', 'ben02syn', 'single_column')
        call ictsz_file
      case ('fuk95')
        call ictsz_fuk95
      case ('channel')
        call ictsz_channel
      case ('isomip1')
        ! call ictsz_isomip1
      case ('isomip2')
        ! call ictsz_isomip2
      case default
        if (mnproc == 1) then
          write (lp,'(3a)') ' inicon: expcnf = ', trim(expcnf), &
               ' is unsupported!'
        end if
        call xcstop('(inicon)')
        stop '(inicon)'
    end select

    ! ------------------------------------------------------------------
    ! Set minimum physical temperature for each isopycnic layer
    ! ------------------------------------------------------------------

    call settemmin

    ! ------------------------------------------------------------------
    ! Initialize configuration specific variables
    ! ------------------------------------------------------------------

    select case (trim(expcnf))
    case ('cesm')
      call inicon_cesm
    case ('ben02clim', 'ben02syn', 'single_column')
      call inicon_ben02
    end select

    ! ------------------------------------------------------------------
    ! Make sure layer temperature is greater than the lower physical
    ! bound and make temperature, salinity, and potential density
    ! variables consistent.
    ! ------------------------------------------------------------------

    select case (vcoord_type_tag)

    case (isopyc_bulkml)

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

    case (cntiso_hybrid)

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

    case default

      if (mnproc == 1) then
        write (lp,*) 'inicon: unsupported vertical coordinate!'
      end if
      call xcstop('(inicon)')
      stop '(inicon)'

    end select

    if (mnproc == ptest) then
      write (lp,'('' sigmar(k)    :'',7f9.5/(15x,7f9.5))') &
           (sigmar(itest,jtest,k),k = 1,kk)
    end if

    ! ------------------------------------------------------------------
    ! Find layer interface pressure
    ! ------------------------------------------------------------------

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

    ! ------------------------------------------------------------------
    ! Set layer thickness and bottom pressure
    ! ------------------------------------------------------------------

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

    if (vcoord_type_tag == cntiso_hybrid) then
      call cntiso_hybrid_regrid_direct_remap(2,1,kk,0,kk+1,1)
      call remap_velocity(2,1,kk,0,kk+1,1)
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

    ! ------------------------------------------------------------------
    ! Initialize potential vorticity of barotropic flow
    ! ------------------------------------------------------------------

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

    ! ------------------------------------------------------------------
    ! separate baroclinic and barotropic velocity components
    ! ------------------------------------------------------------------

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

    ! ------------------------------------------------------------------
    ! Initialize fields related to the pressure gradient force
    ! ------------------------------------------------------------------

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

    ! ------------------------------------------------------------------
    ! Define first physical interior layer
    ! ------------------------------------------------------------------

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

    ! ------------------------------------------------------------------
    ! Set other time level layer thicknesses
    ! ------------------------------------------------------------------

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
           M_mks2cgs*sig(temp(i,j,k),saln(i,j,k)), &
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
