! ------------------------------------------------------------------------------
! Copyright (C) 2015-2025 Mats Bentsen, Ping-Gin Chiu, Mehmet Ilicak,
!                         Aleksi Nummelin, Mariana Vertenstein
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

module mod_geoenv

  use mod_config,    only: inst_suffix
  use mod_constants, only: rearth, pi, radian
  use mod_xc,        only: xchalt, xcaput, xcbcst, itdm, jtdm, &
                           lp, nbdy, i0, ii, j0, jj,  mnproc
  use mod_diffusion, only: rhsctp, tbfile
  use mod_grid,      only: grfile, qclon, qclat, pclon, pclat, uclon, &
                           uclat, vclon, vclat, scqx, scqy, scpx, scpy, &
                           scux, scuy, scvx, scvy, scq2, scp2, scu2, &
                           scv2, qlon, qlat, plon, plat, ulon, ulat, &
                           vlon, vlat, depths, corioq, coriop, betafp, &
                           betatp, angle, cosang, sinang, hangle, nwp
  use mod_utility,   only: fnmlen
  use netcdf

  implicit none
  private

  public :: geoenv_file
  public :: geoenv_test

contains

  subroutine geoenv_file

    ! ------------------------------------------------------------------
    ! Get bathymetry and grid specification from file and compute
    ! Coriolis parameter
    ! ------------------------------------------------------------------

    ! Local variables
    integer, parameter :: cwmlen = 100
    character (len = fnmlen) :: nlfnm
    character (len = 80), dimension(cwmlen) :: cwmtag
    character (len = 1), dimension(cwmlen) :: cwmedg
    real, dimension(itdm,jtdm) :: tmpg
    real, dimension(cwmlen) :: cwmwth
    integer, dimension(cwmlen) :: cwmi,cwmj
    integer, dimension(3) :: start,count
    integer :: i,j,k,status,ncid,dimid,varid,nfu,ios,ncwm,l
    logical :: fexist

    namelist /cwmod/ cwmtag,cwmedg,cwmi,cwmj,cwmwth

    ! ------------------------------------------------------------------
    ! read grid information from grid file
    ! ------------------------------------------------------------------

    if (mnproc == 1) then
      write (lp,'(2a)') ' reading grid information from ',trim(grfile)
      call flush(lp)

      ! open netcdf file
      status = nf90_open(grfile,nf90_nowrite,ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_open: ',trim(grfile),': ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if

      ! check dimensions
      status = nf90_inq_dimid(ncid,'x',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = i)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: x: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_inq_dimid(ncid,'y',dimid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status=nf90_inquire_dimension(ncid,dimid,len = j)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inquire_dimension: y: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      if (i /= itdm.or.j /= jtdm) then
        write (lp,'(2a)') ' wrong dimensions in ',trim(grfile)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if

      ! read bathymetry
      status = nf90_inq_varid(ncid,'pdepth',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: pdepth: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: pdepth: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if

      ! count number of wet points for subsequent xcsum testing
      nwp = 0
      do j = 1,jtdm
        do i = 1,itdm
          if (tmpg(i,j) > 0.) nwp = nwp+1
        end do
      end do
    end if

    call xcaput(tmpg,depths,1)

    ! read grid coordinates

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qlon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qlon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: qlon: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,qlon,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qlat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qlat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: qlat: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,qlat,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'plon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: plon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: plon: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,plon,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'plat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: plat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: plat: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,plat,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'ulon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: ulon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: ulon: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,ulon,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'ulat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: ulat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: ulat: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,ulat,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'vlon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: vlon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: vlon: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,vlon,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'vlat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: vlat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: vlat: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,vlat,1)

    start(1) = 1
    start(2) = 1
    count(1) = itdm
    count(2) = jtdm
    count(3) = 1

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qclon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qclon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: qclon: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,qclon(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qclat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qclat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: qclat: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,qclat(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'pclon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: pclon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: pclon: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,pclon(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'pclat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: pclat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: pclat: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,pclat(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'uclon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: uclon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: uclon: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,uclon(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'uclat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: uclat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: uclat: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,uclat(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'vclon',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: vclon: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: vclon: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,vclon(1-nbdy,1-nbdy,k),1)
    end do

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'vclat',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: vclat: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    do k = 1,4
      if (mnproc == 1) then
        start(3) = k
        status = nf90_get_var(ncid,varid,tmpg,start,count)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: vclat: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,vclat(1-nbdy,1-nbdy,k),1)
    end do

    ! read scale factors

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qdx',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qdx: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: qdx: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scqx,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qdy',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qdy: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: qdy: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scqy,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'pdx',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: pdx: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: pdx: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scpx,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'pdy',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: pdy: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: pdy: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scpy,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'udx',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: udx: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: udx: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scux,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'udy',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: udy: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: udy: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scuy,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'vdx',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: vdx: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: vdx: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scvx,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'vdy',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: vdy: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: vdy: ',nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scvy,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'qarea',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: qarea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: qarea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scq2,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'parea',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: parea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: parea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scp2,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'uarea',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: uarea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: uarea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scu2,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'varea',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: varea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: varea: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,scv2,1)

    if (mnproc == 1) then
      status = nf90_inq_varid(ncid,'angle',varid)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_inq_varid: angle: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    if (mnproc == 1) then
      status = nf90_get_var(ncid,varid,tmpg)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') ' nf90_get_var: angle: ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if
    call xcaput(tmpg,angle,1)

    ! close grid information file

    if (mnproc == 1) then
      status = nf90_close(ncid)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') ' nf90_close: ',trim(grfile),': ', &
             nf90_strerror(status)
        call xchalt('(geoenv_file)')
        stop '(geoenv_file)'
      end if
    end if

    ! ------------------------------------------------------------------
    ! read topographic beta if needed
    ! ------------------------------------------------------------------

    if (rhsctp) then

      if (mnproc == 1) then
        write (lp,'(2a)') ' reading topographic beta from ', &
             trim(tbfile)
        call flush(lp)
        status = nf90_open(tbfile,nf90_nowrite,ncid)
        if (status /= nf90_noerr) then
          write(lp,'(4a)') ' nf90_open: ',trim(tbfile),': ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
        status = nf90_inq_varid(ncid,'topo_beta',varid)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: topo_beta: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
        status = nf90_get_var(ncid,varid,tmpg)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: topo_beta: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,betatp,1)

      if (mnproc == 1) then
        status = nf90_inq_varid(ncid,'hangle',varid)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: hangle: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
        status = nf90_get_var(ncid,varid,tmpg)
        if (status /= nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: hangle: ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      call xcaput(tmpg,hangle,1)

      if (mnproc == 1) then
        status = nf90_close(ncid)
        if (status /= nf90_noerr) then
          write(lp,'(4a)') ' nf90_close: ',trim(tbfile),': ', &
               nf90_strerror(status)
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if

    end if
    ! ------------------------------------------------------------------
    ! Apply channel width modifications if specified in namelist
    ! ------------------------------------------------------------------

    if (mnproc == 1) then
      nlfnm = 'ocn_in'//trim(inst_suffix)
      inquire(file=nlfnm,exist = fexist)
      if (fexist) then
        open (newunit=nfu,file=nlfnm,status='old',action = 'read')
      else
        nlfnm = 'limits'//trim(inst_suffix)
        inquire(file=nlfnm,exist = fexist)
        if (fexist) then
          open (newunit=nfu,file=nlfnm,status='old',action = 'read')
        else
          write (lp,*) 'geoenv_file: could not find namelist file!'
          call xchalt('(geoenv_file)')
          stop '(geoenv_file)'
        end if
      end if
      cwmtag = ''
      read (unit=nfu,nml=cwmod,iostat = ios)
      close (unit = nfu)
    end if
    call xcbcst(ios)
    if (ios /= 0) then
      if (mnproc == 1) then
        write (lp,*) 'geoenv_file: no valid channel width '// &
             'modifications found in namelist.'
      end if
    else
      if (mnproc == 1) then
        ncwm = 1
        do while (cwmtag(ncwm) /= '')
          ncwm = ncwm+1
          if (ncwm > cwmlen) exit
        end do
        ncwm = ncwm-1
      end if
      call xcbcst(ncwm)
      do l = 1,ncwm
        if (mnproc == 1) then
          write (lp,'(a,i3,a)') ' channel width modification number ', &
               l,':'
          write (lp,'(3a)') '   cwmtag = ''',trim(cwmtag(l)),''''
          write (lp,'(3a)') '   cwmedg = ''',cwmedg(l),''''
          write (lp,'(a,i5)') '   cwmi   =',cwmi(l)
          write (lp,'(a,i5)') '   cwmj   =',cwmj(l)
          write (lp,'(a,e10.3)') '   cwmwth = ',cwmwth(l)
          if (cwmedg(l) /= 'u'.and.cwmedg(l) /= 'v') then
            write (lp,*) &
                 'geoenv_file: channel width modification grid cell edge '// &
                 'specification must be either ''u'' or ''v''!'
            call xchalt('(geoenv_file)')
            stop '(geoenv_file)'
          end if
          if (cwmi(l) < 1.or.cwmi(l) > itdm.or. &
               cwmj(l) < 1.or.cwmj(l) > jtdm) then
            write (lp,*) &
                 'geoenv_file: channel width modification indices out of bounds!'
            call xchalt('(geoenv_file)')
            stop '(geoenv_file)'
          end if
        end if
        call xcbcst(cwmtag(l))
        call xcbcst(cwmedg(l))
        call xcbcst(cwmi(l))
        call xcbcst(cwmj(l))
        call xcbcst(cwmwth(l))
        if (cwmi(l) > i0.and.cwmi(l) <= i0+ii.and. &
             cwmj(l) > j0.and.cwmj(l) <= j0+jj) then
          i = cwmi(l)-i0
          j = cwmj(l)-j0
          if (cwmedg(l) == 'u') then
            write (lp,*) 'blom: geoenv_file: ',trim(cwmtag(l)), &
                 i+i0,j+j0,'scuy:',scuy(i,j),'->',cwmwth(l)
            scu2(i,j) = scu2(i,j)*cwmwth(l)/scuy(i,j)
            scuy(i,j) = cwmwth(l)
          else
            write (lp,*) 'blom: geoenv_file: ',trim(cwmtag(l)), &
                 i+i0,j+j0,'scvx:',scvx(i,j),'->',cwmwth(l)
            scv2(i,j) = scv2(i,j)*cwmwth(l)/scvx(i,j)
            scvx(i,j) = cwmwth(l)
          end if
        end if
      end do
    end if

    ! ------------------------------------------------------------------
    ! Precompute cosine and sine of local angle of i-direction with
    ! eastward direction, and compute Coriolis and beta plane parameter
    ! ------------------------------------------------------------------

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii

        cosang(i,j) = cos(angle(i,j))
        sinang(i,j) = sin(angle(i,j))

        corioq(i,j) = sin(qlat(i,j)/radian)*4.*pi/86164.
        coriop(i,j) = sin(plat(i,j)/radian)*4.*pi/86164.
        betafp(i,j) = cos(plat(i,j)/radian)*4.*pi/(86164.*rearth)

      end do
    end do
    !$omp end parallel do

  end subroutine geoenv_file

  subroutine geoenv_test

    ! ------------------------------------------------------------------
    ! Define bathymetry, grid specification and Coriolis parameter for
    ! test case
    ! ------------------------------------------------------------------

    depths = 0.
    nwp = 0
    qlon = 0.
    qlat = 0.
    plon = 0.
    plat = 0.
    ulon = 0.
    ulat = 0.
    vlon = 0.
    vlat = 0.
    qclon = 0.
    qclat = 0.
    pclon = 0.
    pclat = 0.
    uclon = 0.
    uclat = 0.
    vclon = 0.
    vclat = 0.
    scqx = 0.
    scqy = 0.
    scpx = 0.
    scpy = 0.
    scux = 0.
    scuy = 0.
    scvx = 0.
    scvy = 0.
    scq2 = 0.
    scp2 = 0.
    scu2 = 0.
    scv2 = 0.
    angle = 0.
    corioq = 0.
    coriop = 0.
    betafp = 0.
    betatp = 0.
    cosang = 0.
    sinang = 0.
    hangle = 0.

  end subroutine geoenv_test

end module mod_geoenv
