! ------------------------------------------------------------------------------
! Copyright (C) 2004-2024 Ingo Bethke, Mats Bentsen, Alok Kumar Gupta

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_nctools
  ! ----------------------------------------------------------------------
  ! --- NetCDF tools -----------------------------------------------------
  ! ----------------------------------------------------------------------

  ! Description:
  !    nctools.f is a boundle of fortran subroutines designed to write
  !    output from the mpi version of BLOM to NetCDF data-format.

  !    Prerequisites are the linkage of the f90 netcdf library and the
  !    access to the module netcdf.mod.

  ! Comments:
  !    The order of actions is relevant: i.e. dimensions should be
  !    declared iniflag, then fields not containing the unlimited (time)
  !    dimension, and finally the fields containing the time dimension.
  !    Attributes may be defined any time.

  ! Contents:

  !    ncfopn        - creates a file or opens existing file for reading
  !    ncfcls        - closes an opened nc-file

  !    ncdims        - defines a simple axis-dimension
  !    ncdimc        - defines a dimension for compressed data storage
  !                    and optionally stores the index information that
  !                    is necessary to de-compress data
  !    nctime        - defines a time dimension and stores a time value

  !    ncattr        - adds a text attribute

  !    ncputr        - writes a skalar or vector in real format
  !    ncputi        - writes a skalar or vector in integer format
  !    ncgetr        - reads a skalar or vector in real format
  !    ncgeti        - reads a skalar or vector in integer format

  !    ncread        - reads a 2d or 3d field

  !    ncwrtc        - writes a string array
  !    ncwrti        - writes a 2d or 3d field as int4
  !    ncwrtr        - writes a 2d or 3d field as real8
  !    ncpack        - writes a 2d or 3d field in packed format as int2
  !                    with scale factor and offset
  !    nccomp        - writes a 2d or 3d field in compressed (skipping
  !                    land points) in real8 format
  !    nccopa        - writes a 2d or 3d field in compressed/packed format
  !                    as int2 with scale factor and offset

  !    ncerro        - displays error massage
  !    ncsevl        - evaluates strings with deliminators ' ',':' and '-'
  !    ncinqv        - inquires if variable exits

  ! Revision history:
  !    may2008       - switched from F77-API to F90-API
  !    may2008       - removed precision argument from ncwrti and ncwrtr
  !    apr2008       - added header padding
  !    apr2008       - added ncputr, ncputi, ncgetr, ncgeti
  !    mar2008       - added nccomp, nccopa, ncread, ncdimc
  !    mar2008       - added MPI support

  ! Contact:
  !    Ingo Bethke (ingo.bethke@nersc.no)

  ! ----------------------------------------------------------------------

  use dimensions,   only: itdm, jtdm, idm, jdm, iqr, jqr, ijqr
  use mod_xc,       only: xcstop, xchalt, xcaget, xcaput, xcmin, xcmax, &
                          xcbcst, ii, jj, kk, nbdy, lp, i0, j0, &
                          mpe_1, mnproc, nproc, mproc
#ifdef PNETCDF
  use mod_xc,       only: xcgetrow, xcgetrowint2, xcgetrow4
  use mod_xc,       only: mpicomm,mpierr,mpireq,mpistat
#endif
  use mod_calendar, only: date_type, daynum_diff, calendar_noerr, &
                          calendar_errstr
  use netcdf

  implicit none
  public

#ifdef PNETCDF
#  include <pnetcdf.inc>
#  include <mpif.h>
#endif

  interface ncputr
    module procedure ncputrs,ncputrv
  end interface ncputr

  interface ncputi
    module procedure ncputis,ncputiv
  end interface ncputi

  interface ncgetr
    module procedure ncgetrs,ncgetrv
  end interface ncgetr

  interface ncgeti
    module procedure ncgetis,ncgetiv
  end interface ncgeti

  ! --- Declare global variables
  logical, private :: flgpad
  integer, private :: ncid,rhid,status,rec,io_type
#ifdef PNETCDF
  integer(kind = mpi_offset_kind) :: clen,istart(5),icount(5),tkd
#endif
  integer*2, parameter :: i2fill=-32768,i2max = 32767
  real, parameter :: fillr8 = 9.9692099683868690e+36
  real(kind=4), parameter :: fillr4 = 9.9692099683868690e+36
  integer :: ndouble,nchar,nfint

contains


  subroutine ncfopn(fnm,faccess,frmt,irec,iotype)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Opens a NetCDF file for reading or writing.
    ! --- Arguments:
    !       char(*) fnm      (in) -  file name
    !       char(*) faccess  (in) -  'w' for write and 'r' for read
    !       char(*) frmt (in)     -  'CLASSIC','64BIT','HDF5'
    !       integer irec          -  time step to be written
    ! ----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: fnm
    character(len=*), intent(in) :: faccess
    character(len=*), intent(in) :: frmt
    integer,          intent(in) :: irec
    integer,          intent(in) :: iotype

    ! Local variables
    integer :: oldmode
    integer, parameter :: nf90__64bit_offset = 512
    integer, parameter :: nf90__hdf5 = 4096
#ifdef PNETCDF
    integer*4, save :: info = MPI_INFO_NULL
    character(len = 3) :: stripestr
    character(len = 9) :: stripestr2
    integer :: ierr
#endif
    integer :: testio

    ! --- Open file
    IO_TYPE = iotype
    testio = 0
    if(io_type  ==  1) then
#ifdef PNETCDF
      testio = 1
      ndouble = nf_double
      nchar = nf_char
      nfint = nf_int
      rec = 1
      write(stripestr,('(i3)')) 8
      write(stripestr2,('(i9)')) 1024*1024
      if (faccess(1:1) == 'r') then
        call ncerro(nfmpi_open(mpicomm,fnm,nf_nowrite,MPI_INFO_NULL,ncid))
        rec = irec
      else if (irec /= 1) then
        status = nfmpi_open(mpicomm,fnm,nf_write,INFO,ncid)
        if (status /= nf_noerr) then
          call mpi_info_create(info,ierr)
          call mpi_info_set(info,'romio_ds_read','disable',ierr)
          call mpi_info_set(info,'romio_ds_write','disable',ierr)
          call mpi_info_set(info,"striping_factor",stripestr,ierr)
          call mpi_info_set(info,"striping_unit",stripestr2,ierr)
          if (frmt(1:1) == '6') then
            call ncerro(nfmpi_create(mpicomm,fnm,IOR(nf_clobber,nf_64bit_offset),INFO,ncid))
          else if (frmt(1:1) == 'h'.or.frmt(1:1) == 'h') then
            call ncerro(nfmpi_create(mpicomm,fnm,IOR(nf_clobber,nf_64bit_data),INFO,ncid))
          else
            call ncerro(nfmpi_create(mpicomm,fnm,IOR(NF_CLOBBER,nf_format_classic),INFO,ncid))
          end if
        else
          call ncerro(nfmpi_redef(ncid))
          rec = irec
        end if
      else
        call mpi_info_create(info,ierr)
        call mpi_info_set(info,'romio_ds_read','disable',ierr)
        call mpi_info_set(info,'romio_ds_write','disable',ierr)
        call mpi_info_set(info,'striping_factor',stripestr,ierr)
        call mpi_info_set(info,'striping_unit',stripestr2,ierr)
        call mpi_info_set(info,'striping_unit',stripestr2,ierr)
        call mpi_info_set(info,'nc_header_align_size',stripestr2,ierr)
        if (frmt(1:1) == '6') then
          call ncerro(nfmpi_create(mpicomm,fnm,IOR(nf_clobber,nf_64bit_offset),INFO,ncid))
        else if (frmt(1:1) == 'h'.or.frmt(1:1) == 'h') then
          call ncerro(nfmpi_create(mpicomm,fnm,IOR(nf_clobber,nf_64bit_data),INFO,ncid))
        else
          call ncerro(nfmpi_create(mpicomm,fnm,IOR(nf_clobber,nf_format_classic),INFO,ncid))
        end if
      end if
      ! --- Initialise header padding
      if (rec == 1) then
        flgpad = .false.
      else
        flgpad = .true.
      end if
#endif
      if(testio  ==  0) then
        write(lp,*) 'check iotype in namelist'
        call xchalt('(ncerro)')
        stop '(ncerro)'
      end if
    else if (io_type  ==  0) then
      ndouble = nf90_double
      nchar = nf90_char
      nfint = nf90_int
      if (mnproc == 1) then

        ! --- - Open file
        rec = 1
        if (faccess(1:1) == 'r') then
          call ncerro(nf90_open(fnm,nf90_nowrite,ncid))
          rec = irec
        else if (irec /= 1) then
          status = nf90_open(fnm,nf90_write,ncid)
          if (status /= nf90_noerr) then
            if (frmt(1:1) == '6') then
              call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__64bit_offset),ncid))
            else if (frmt(1:1) == 'h'.or.frmt(1:1) == 'h') then
              call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__hdf5),ncid))
            else
              call ncerro(nf90_create(fnm,nf90_clobber,ncid))
            end if
          else
            call ncerro(nf90_redef(ncid))
            rec = irec
          end if
          call ncerro(nf90_set_fill(ncid,nf90_nofill,oldmode))
        else
          if (frmt(1:1) == '6') then
            call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__64bit_offset),ncid))
          else if (frmt(1:1) == 'h'.or.frmt(1:1) == 'h') then
            call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__hdf5),ncid))
          else
            call ncerro(nf90_create(fnm,nf90_clobber,ncid))
          end if
          call ncerro(nf90_set_fill(ncid,nf90_nofill,oldmode))
        end if

        ! --- - Set variable id to global (needed for global attributes)
        rhid = nf90_global

        ! --- - Initialise header padding
        if (rec == 1) then
          flgpad = .false.
        else
          flgpad = .true.
        end if

      end if
    end if

  end subroutine ncfopn



  subroutine ncfcls

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Closes NetCDF file
    ! ----------------------------------------------------------------------

    if (io_type  ==  1) then
#ifdef PNETCDF
      call ncerro(nfmpi_close(ncid))
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1) then
        call ncerro(nf90_close(ncid))
      end if
    end if

  end subroutine ncfcls


  subroutine ncdims(dnm,dim)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Creates a simple dimension
    ! --- Arguments:
    !       char(*) dnm  (in) -  name of the dimension
    !       integer dim  (in) -  number of values on axis
    ! ----------------------------------------------------------------------

    character(len=*), intent(in) :: dnm
    integer, intent(in) :: dim

    integer :: dimid

    ! --- define dimension
    if (io_type  ==  1) then
#ifdef PNETCDF
      clen = dim
      if (rec == 1) then
         call ncerro(nfmpi_def_dim(ncid,trim(dnm),clen,dimid))
      end if
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1.and.rec == 1) then
         call ncerro(nf90_def_dim(ncid,dnm,dim,dimid))
      end if
    end if

  end subroutine ncdims



  subroutine ncdimc(dnm,msk,flg)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Creates a dimension for compressed data storage and optionally
    !       stores the index information that is neccessary to de-compress
    !       data.
    ! --- Arguments:
    !       char(*)      dnm (in) -  name of dimension, e.g. 'pcomp'
    !       integer(...) msk (in) -  2d mask with dimensions
    !                                (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    !       integer      flg (in) -  1 if index variable should be written
    ! ----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: dnm
    integer,          intent(in), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: msk
    integer,          intent(in) :: flg

    ! Local variables
    integer :: n,i,j,dimid
    integer, dimension(itdm*jtdm) :: ivec
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rmsk
    real, dimension(itdm,jtdm) :: rmskt

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        rmsk(i,j) = msk(i,j)
      end do
    end do
    !$omp end parallel do
    call xcaget(rmskt,rmsk,1)
    if (mnproc == 1.and.rec == 1) then
      n = 0
      if (flg == 1) then
        do j = 1,jtdm
          do i = 1,itdm
            if (rmskt(i,j) > 0.5) then
              n = n+1
              ivec(n) = i+(j-1)*itdm
            end if
          end do
        end do
      else
        do j = 1,jtdm
          do i = 1,itdm
            if (rmskt(i,j) > 0.5) n = n+1
          end do
        end do
      end if
    end if
    if (io_type  ==  0) then
      if (mnproc == 1.and.rec == 1) then
        call ncerro(nf90_def_dim(ncid,dnm,n,dimid))
        if (flg == 1) then
          call ncerro(nf90_def_var(ncid,dnm,nf90_int,dimid,rhid))
          if (flgpad) then
            call ncerro(nf90_enddef(ncid))
          else
            call ncerro(nf90_enddef(ncid,81920,4,40960,4))
            flgpad = .true.
          end if
          call ncerro(nf90_put_var(ncid,rhid,ivec(1:n)))
          call ncerro(nf90_redef(ncid))
        end if
      end if
    else if (io_type  ==  1) then
#ifdef PNETCDF
      call xcbcst(n)
      clen = n
      if (rec == 1) then
        call ncerro(nfmpi_def_dim(ncid,dnm,clen,dimid))
        if (flg == 1) then
          call ncerro(nfmpi_def_var(ncid,dnm,nf_int,1,dimid,rhid))
          call ncerro(nfmpi_enddef(ncid))
          istart(1) = 1
          icount(1) = 0
          if (mnproc == 1) icount(1) = n
          call ncerro(nfmpi_put_vara_int_all(ncid,rhid,istart,icount,ivec(1:n)))
          call ncerro(nfmpi_redef(ncid))
        end if
      end if
#endif
    end if

  end subroutine ncdimc



  subroutine ncattr(attname,attrib)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Adds a single text attribute to the variable that has been
    !       created last. If no variable has been created then a global
    !       attribute is added.
    ! --- Arguments:
    !       char(*) attname (in) -  attribute name
    !       char(*) attrib  (in) -  attribute text
    ! ----------------------------------------------------------------------

    character(len=*), intent(in) ::  attname,attrib

    if (io_type  ==  1) then
#ifdef PNETCDF
      clen = len(trim(attrib))
      if (rec == 1) then
         call ncerro(nfmpi_put_att_text(ncid,rhid,attname,clen,attrib))
      end if
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1.and.rec == 1) then
         call ncerro(nf90_put_att(ncid,rhid,attname,attrib))
      end if
    end if

  end subroutine ncattr



  subroutine ncputrs(vnm,rval)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Write a skalar as global attribute in real8 format.
    ! --- Arguments:
    !       char(*) vnm    (in) -  variable name
    !       real rval      (in) -  real skalar
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    real, intent(in) :: rval

    if (io_type  ==  1) then
#ifdef PNETCDF
      clen = 1
      if (rec == 1) then
         call ncerro(nfmpi_put_att_double(ncid,nf_global,vnm,nf_double,clen,rval))
      end if
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1.and.rec == 1) then
         call ncerro(nf90_put_att(ncid,nf90_global,vnm,rval))
      end if
    end if

  end subroutine ncputrs



  subroutine ncputrv(vnm,rval)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Write a vector as global attribute in real8 format.
    ! --- Arguments:
    !       char(*) vnm    (in) -  variable name
    !       real rval(:)   (in) -  real vector
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    real, intent(in) :: rval(:)

    if (io_type  ==  1) then
#ifdef PNETCDF
      clen = size(rval)
      if (rec == 1) then
         call ncerro(nfmpi_put_att_double(ncid,nf_global,vnm,nf_double,clen,rval))
      end if
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1.and.rec == 1) then
         call ncerro(nf90_put_att(ncid,nf90_global,vnm,rval))
      end if
    end if

  end subroutine ncputrv



  subroutine ncputis(vnm,ival)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes a scalar as global attribute in int4 format.
    ! --- Arguments:
    !       char(*) vnm       (in) -  variable name
    !       integer ival      (in) -  integer skalar
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    integer, intent(in) :: ival

    if (io_type  ==  1) then
#ifdef PNETCDF
      clen = 1
      if (rec == 1) then
         call ncerro(nfmpi_put_att_int(ncid,nf_global,vnm,nf_int,clen,ival))
      end if
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1.and.rec == 1) then
         call ncerro(nf90_put_att(ncid,nf90_global,vnm,ival))
      end if
    end if

  end subroutine ncputis



  subroutine ncputiv(vnm,ival)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes a vector as global attribute in int4 format.
    ! --- Arguments:
    !       char(*) vnm       (in) -  variable name
    !       integer ival(:)   (in) -  integer vector
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    integer, intent(in) :: ival(:)

    if (io_type  ==  1) then
#ifdef PNETCDF
      clen = size(ival)
      if (rec == 1) then
         call ncerro(nfmpi_put_att_int(ncid,nf_global,vnm,nf_int,clen,ival))
      end if
#endif
    else if (io_type  ==  0) then

      if (mnproc == 1.and.rec == 1) then
         call ncerro(nf90_put_att(ncid,nf90_global,vnm,ival))
      end if
    end if

  end subroutine ncputiv



  subroutine ncgetrs(vnm,rval)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !      Reads a skalar in real8 format.
    ! --- Arguments:
    !       char(*) vnm  (in)  -  variable name
    !       real    rval (out) -  real skalar
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    real, intent(out) :: rval

    if (io_type  ==  1) then
#ifdef PNETCDF
      call ncerro(nfmpi_get_att_double(ncid,nf_global,vnm,rval))
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1) then
         call ncerro(nf90_get_att(ncid,nf90_global,vnm,rval))
      end if
      call xcbcst(rval)
    end if

  end subroutine ncgetrs



  subroutine ncgetrv(vnm,rval)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !      Reads a vector in real8 format.
    ! --- Arguments:
    !       char(*) vnm  (in)  -  variable name
    !       real(:) rval (out) -  real vector
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    real, intent(out) :: rval(:)

    if (io_type  ==  1) then
#ifdef PNETCDF
      call ncerro(nfmpi_get_att_double(ncid,nf_global,vnm,rval))
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1) then
         call ncerro(nf90_get_att(ncid,nf90_global,vnm,rval))
      end if
      call xcbcst(rval)
    end if

  end subroutine ncgetrv



  subroutine ncgetis(vnm,ival)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !      Reads a skalar in integer format.
    ! --- Arguments:
    !       char(*)    vnm  (in)  -  variable name
    !       integer    ival (out) -  integer skalar
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    integer, intent(out) :: ival

    if (io_type  ==  1) then
#ifdef PNETCDF
      call ncerro(nfmpi_get_att_int(ncid,nf_global,vnm,ival))
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1) then
         call ncerro(nf90_get_att(ncid,nf90_global,vnm,ival))
      end if
      call xcbcst(ival)
    end if

  end subroutine ncgetis



  subroutine ncgetiv(vnm,ival)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !      Read a vector in integer format.
    ! --- Arguments:
    !       char(*)    vnm  (in)  -  variable name
    !       integer(:) ival (out) -  integer vector
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    integer, intent(out) :: ival(:)

    if (io_type  ==  1) then
#ifdef PNETCDF
      call ncerro(nfmpi_get_att_int(ncid,nf_global,vnm,ival))
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1) then
         call ncerro(nf90_get_att(ncid,nf90_global,vnm,ival))
      end if
      call xcbcst(ival)
    end if

  end subroutine ncgetiv



  subroutine nctime(datenum,calendar,units,startdate)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Create time dimension and store time value.
    !       Valid calendars are '360_day', 'noleap' = '365_day',
    !       'all_leap' = '366_day', 'julian', 'proleptic_gregorian' and
    !       'standard' = 'mixed' = 'gregorian'.
    ! --- Arguments:
    !       real    datenum   (in) -  time value
    !       char(*) calendar  (in) -  calendar choice
    !       char(*) units     (in) -  time units,
    !                                 e.g. 'hours since 0001-01-01'
    !       char(*) startdate (in) -  start date, e.g. '1997-01-01'
    ! --- Comment:
    !       The script assumes that the precision of startdate is not
    !	    higher than the unit specified in the units string. E.g. if
    !       units='years since 0001-01-01' and startdate='1997-06-01' then
    !       the resulting offset will be rounded to 1996 yrs.
    ! ----------------------------------------------------------------------

    real :: datenum
    character(len=*), intent(in) :: calendar,units,startdate

    integer :: bstrn,sstrn,strind(2,10),by,bm,bd,bh,sy,sm,sd,sh,dndiff
    real :: caloffset
    data by/1/,bm/1/,bd/1/,bh/0/,sy/1/,sm/1/,sd/1/,sh/0/

    integer :: errstat


    if (io_type  ==  1) then
#ifdef PNETCDF
      ! --- - define time dimension
      istart(1) = 1
      icount(1) = 0

      ! --- - define time coordinate variable
      call ncerro(nfmpi_inq_varid(ncid,'time',rhid))

      if (mnproc == 1) then
        ! --- --- analyse units string
        call ncsevl(units,bstrn,strind)
        read(units(strind(1,3):strind(2,3)),*) by
        if (bstrn >= 4) read(units(strind(1,4):strind(2,4)),*) bm
        if (bstrn >= 5) read(units(strind(1,5):strind(2,5)),*) bd
        if (bstrn >= 6) read(units(strind(1,6):strind(2,6)),*) bh

        ! --- --- analyse startdate string
        call ncsevl(startdate,sstrn,strind)
        read(startdate(strind(1,1):strind(2,1)),*) sy
        if (sstrn >= 2) read(startdate(strind(1,2):strind(2,2)),*) sm
        if (sstrn >= 3) read(startdate(strind(1,3):strind(2,3)),*) sd
        if (sstrn >= 4) read(startdate(strind(1,4):strind(2,4)),*) sh

        ! --- --- calculate calendar offset
        if     (units(1:1) == 'y') then
          caloffset = real(sy-by)
        else if (units(1:1) == 'm') then
          caloffset = real(sy-by)*12.+real(sm-bm)
        else
          errstat = daynum_diff(calendar,date_type(by,bm,bd), &
               date_type(sy,sm,sd),dndiff)
          if (errstat /= calendar_noerr) then
            write (lp, '(2a)') ' nctime: daynum_diff error: ',trim(calendar_errstr(errstat))
            call xchalt('(nctime)')
            stop '(nctime)'
          end if
          if     (units(1:1) == 'd') then
            caloffset = real(dndiff)
          else if (units(1:1) == 'h') then
            caloffset = real(dndiff)*24.+real(sh-bh)
          end if
        end if
        datenum = datenum+caloffset
        istart(1) = rec
        icount(1) = 1
      end if

      ! --- - Change to data mode and write time
      call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,datenum))
#endif
    else if (io_type  ==  0) then
      if (mnproc == 1) then
        call ncerro(nf90_inq_varid(ncid,'time',rhid))

        ! --- --- analyse units string
        call ncsevl(units,bstrn,strind)
        read(units(strind(1,3):strind(2,3)),*) by
        if (bstrn >= 4) read(units(strind(1,4):strind(2,4)),*) bm
        if (bstrn >= 5) read(units(strind(1,5):strind(2,5)),*) bd
        if (bstrn >= 6) read(units(strind(1,6):strind(2,6)),*) bh

        ! --- --- analyse startdate string
        call ncsevl(startdate,sstrn,strind)
        read(startdate(strind(1,1):strind(2,1)),*) sy
        if (sstrn >= 2) read(startdate(strind(1,2):strind(2,2)),*) sm
        if (sstrn >= 3) read(startdate(strind(1,3):strind(2,3)),*) sd
        if (sstrn >= 4) read(startdate(strind(1,4):strind(2,4)),*) sh

        ! --- --- calculate calendar offset
        if     (units(1:1) == 'y') then
          caloffset = real(sy-by)
        else if (units(1:1) == 'm') then
          caloffset = real(sy-by)*12.+real(sm-bm)
        else
          errstat = daynum_diff(calendar,date_type(by,bm,bd),date_type(sy,sm,sd),dndiff)
          if (errstat /= calendar_noerr) then
            write (lp, '(2a)') ' nctime: daynum_diff error: ',trim(calendar_errstr(errstat))
            call xchalt('(nctime)')
            stop '(nctime)'
          end if
          if     (units(1:1) == 'd') then
            caloffset = real(dndiff)
          else if (units(1:1) == 'h') then
            caloffset = real(dndiff)*24.+real(sh-bh)
          end if
        end if
        datenum = datenum+caloffset
        call ncerro(nf90_put_var(ncid,rhid,datenum,(/rec/)))
      end if
    end if

  end subroutine nctime



  subroutine ncread(vnm, fld, msk, mskflg, fill, scf_arg)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Reads 2d or 3d field.
    ! --- Arguments:
    !       char(len=*)  vnm    (in)  -  variable name
    !       real(...)    fld    (out) -  field with dimension
    !                                    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,:)
    !       integer(...) msk    (in)  -  field mask, if mskflg=1 then the
    !                                    dimensions assumed to be
    !                                    (idm+2*nbdy)*(jdm+2*nbdy)
    !       integer      mskflg (in)  -  1 if mask exists
    !       real         fill   (in)  -  fill value to be set on land
    ! ----------------------------------------------------------------------

    character(len=*), intent(in)           :: vnm
    real,             intent(out)          :: fld(*)
    integer,          intent(in)           :: msk(*)
    integer,          intent(in)           :: mskflg
    real,             intent(in)           :: fill
    real,             intent(in), optional :: scf_arg

    character(len=100) :: dimname
    integer :: i,j,ij,k,kd,n,ndm
    integer, parameter :: maxdm=5
    integer, parameter :: ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
    integer :: ndims,dimids(maxdm),dimlen
#ifdef PNETCDF
    integer(kind = mpi_offset_kind) :: tdimlen
#endif
    integer, dimension(maxdm) :: start,count
    real, allocatable, dimension(:,:,:) :: rfldt
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
    real, dimension(itdm,jtdm) :: rmskt
    real, dimension(itdm*jtdm) :: rfldtcmp
    real :: ofs,scf
    logical :: cmpflg

    ! --- Initialise fields
    cmpflg = .false.
    kd = 1
    if(io_type  ==  1) then
#ifdef PNETCDF
      do n = 1,maxdm
        istart(n) = 1
        icount(n) = 1
      end do

      status = nfmpi_inq_varid(ncid,vnm,rhid)
      if (status /= nf_noerr) then
        if (mnproc == 1) then
          write(lp,*) 'WARNING: Problems reading variable ',trim(vnm)
          call flush(lp)
        end if
        call ncerro(status)
      end if
      call ncerro(nfmpi_inq_varndims(ncid,rhid,ndims))
      call ncerro(nfmpi_inq_vardimid(ncid,rhid,dimids))
      ndm = 1
      do n = 1,ndims
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(n),tdimlen))

        call ncerro(nfmpi_inq_dimname(ncid,dimids(n),dimname))
        if (dimname(2:5) == 'comp') cmpflg = .true.
        if ((n == 2.and.cmpflg).or.(n == 3.and..not.cmpflg)) kd = tdimlen
        ndm = ndm*tdimlen
      end do

      ! --- Get attributes
      status = nfmpi_get_att_double(ncid,rhid,'scale_factor',scf)
      if (status /= nf_noerr) scf = 1.
      if (present(scf_arg)) scf = scf*scf_arg
      status = nfmpi_get_att_double(ncid,rhid,'add_offset',ofs)
      if (status /= nf_noerr) ofs = 0.

      ! --- get mask for the full domain
      if (mskflg == 1) then
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          end do
        end do
        !$omp end parallel do
        call xcaget(rmskt,rmsk,1)
      end if

      ! --- Get data
      if (cmpflg) then
        allocate(rfldt(itdm,jtdm,1))
        rfldt = 0.0
        do k = 1,kd
          if (mnproc == 1) then
            istart(2) = k
            icount(1) = nint(dble(ndm)/dble(kd))
          else
            icount(1) = 0
            icount(2) = 0
          end if
          call ncerro(nfmpi_get_vara_double_all(ncid,rhid,istart,icount,rfldtcmp))
          if (mnproc == 1) then
            n = 0
            do j = 1,jtdm
              do i = 1,itdm
                if (rmskt(i,j) > 0.5) then
                  n = n+1
                  rfldt(i,j,1) = rfldtcmp(n)
                end if
              end do
            end do
          end if
          call xcaput(rfldt,rfld,1)
          if (mskflg == 1) then
            if (scf == 1..and.ofs == 0) then
              !$omp parallel do private(i,ij)
              do j = 1,jj
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij) == 1) then
                    fld(ij+(k-1)*ijdm) = rfld(i,j)
                  else
                    fld(ij+(k-1)*ijdm) = fill
                  end if
                end do
              end do
              !$omp end parallel do
            else
              !$omp parallel do private(i,ij)
              do j = 1,jj
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij) == 1) then
                    fld(ij+(k-1)*ijdm) = rfld(i,j)*scf+ofs
                  else
                    fld(ij+(k-1)*ijdm) = fill
                  end if
                end do
              end do
              !$omp end parallel do
            end if
          else
            if (scf == 1..and.ofs == 0) then
              !$omp parallel do private(i,ij)
              do j = 1,jj
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm) = rfld(i,j)
                end do
              end do
              !$omp end parallel do
            else
              !$omp parallel do private(i,ij)
              do j = 1,jj
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm) = rfld(i,j)*scf+ofs
                end do
              end do
              !$omp end parallel do
            end if
          end if
        end do
      else
        istart(1) = i0+1
        istart(2) = j0+1
        istart(3) = 1
        icount(1) = ii
        icount(2) = jj
        icount(3) = kd
        allocate(rfldt(ii,jj,kd))
        call ncerro(nfmpi_get_vara_double_all(ncid,rhid,istart,icount,rfldt))

        if (mskflg == 1) then
          if (scf == 1..and.ofs == 0) then
            !$omp parallel do private(k,i,ij)
            do j = 1,jj
              do k = 1,kd
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij) == 1) then
                    fld(ij+(k-1)*ijdm) = rfldt(i,j,k)
                  else
                    fld(ij+(k-1)*ijdm) = fill
                  end if
                end do
              end do
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(k,i,ij)
            do j = 1,jj
              do k = 1,kd
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij) == 1) then
                    fld(ij+(k-1)*ijdm) = rfldt(i,j,k)*scf+ofs
                  else
                    fld(ij+(k-1)*ijdm) = fill
                  end if
                end do
              end do
            end do
            !$omp end parallel do
          end if
        else
          if (scf == 1..and.ofs == 0) then
            !$omp parallel do private(k,i,ij)
            do j = 1,jj
              do k = 1,kd
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm) = rfldt(i,j,k)
                end do
              end do
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(k,i,ij)
            do j = 1,jj
              do k = 1,kd
                do i = 1,ii
                  ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm) = rfldt(i,j,k)*scf+ofs
                end do
              end do
            end do
            !$omp end parallel do
          end if
        end if
      end if
      deallocate(rfldt)
#endif
    else if(io_type  ==  0) then
      do n = 1,maxdm
        start(n) = 1
        count(n) = 1
      end do

      ! --- Inquire dimensions, ect.
      if (mnproc == 1) then
        status = nf90_inq_varid(ncid,vnm,rhid)
        if (status /= nf90_noerr) then
          write(lp,*) 'WARNING: Problems reading variable ',trim(vnm)
          call flush(lp)
          call ncerro(status)
        end if
        call ncerro(nf90_inquire_variable(ncid,rhid,ndims = ndims))
        call ncerro(nf90_inquire_variable(ncid,rhid,dimids = dimids))
        ndm = 1
        do n = 1,ndims
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),len = dimlen))
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),name = dimname))
          if (dimname(2:5) == 'comp') cmpflg = .true.
          if ((n == 2.and.cmpflg).or.(n == 3.and..not.cmpflg)) kd = dimlen
          ndm = ndm*dimlen
        end do
      end if
      call xcbcst(cmpflg)
      call xcbcst(kd)

      ! --- Get attributes
      if (mnproc == 1) then
        status = nf90_get_att(ncid,rhid,'scale_factor',scf)
        if (status /= nf90_noerr) scf = 1.
        status = nf90_get_att(ncid,rhid,'add_offset',ofs)
        if (status /= nf90_noerr) ofs = 0.
      end if
      call xcbcst(scf)
      call xcbcst(ofs)

      ! --- get mask for the full domain
      if (mskflg == 1) then
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          end do
        end do
        !$omp end parallel do
        call xcaget(rmskt,rmsk,1)
      end if
      allocate(rfldt(itdm,jtdm,1))

      ! --- Get data
      do k = 1,kd
        if (mnproc == 1) then
          if (cmpflg) then
            start(2) = k
            count(1) = nint(dble(ndm)/dble(kd))
          else
            start(3) = k
            count(1) = itdm
            count(2) = jtdm
          end if
          if (cmpflg) then
            call ncerro(nf90_get_var(ncid,rhid,rfldtcmp,start,count))
            n = 0
            do j = 1,jtdm
              do i = 1,itdm
                if (rmskt(i,j) > 0.5) then
                  n = n+1
                  rfldt(i,j,1) = rfldtcmp(n)
                end if
              end do
            end do
          else
            call ncerro(nf90_get_var(ncid,rhid,rfldt,start,count))
          end if
        end if
        call xcaput(rfldt,rfld,1)
        if (mskflg == 1) then
          if (scf == 1..and.ofs == 0) then
            !$omp parallel do private(i,ij)
            do j = 1,jj
              do i = 1,ii
                ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                if (msk(ij) == 1) then
                  fld(ij+(k-1)*ijdm) = rfld(i,j)
                else
                  fld(ij+(k-1)*ijdm) = fill
                end if
              end do
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(i,ij)
            do j = 1,jj
              do i = 1,ii
                ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                if (msk(ij) == 1) then
                  fld(ij+(k-1)*ijdm) = rfld(i,j)*scf+ofs
                else
                  fld(ij+(k-1)*ijdm) = fill
                end if
              end do
            end do
            !$omp end parallel do
          end if
        else
          if (scf == 1..and.ofs == 0) then
            !$omp parallel do private(i,ij)
            do j = 1,jj
              do i = 1,ii
                ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                fld(ij+(k-1)*ijdm) = rfld(i,j)
              end do
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(i,ij)
            do j = 1,jj
              do i = 1,ii
                ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                fld(ij+(k-1)*ijdm) = rfld(i,j)*scf+ofs
              end do
            end do
            !$omp end parallel do
          end if
        end if
      end do
      deallocate(rfldt)
    end if

  end subroutine ncread



  subroutine ncpack(vnm,dims,fld,msk,mskflg,sfac,offs)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes real field to nc-file in packed format using scale factor
    !       and offset.
    ! --- Arguments:
    !       char(*)  vnm    (in) -  variable name
    !       char(*)  dims   (in) -  axes string, e.g. 'x y z time'
    !       real(*)  fld    (in) -  field with dimension
    !       int(*)   msk    (in) -  field mask, if mskflg=1 the dimension
    !                               is assumed (idm+2*nbdy)*(jdm+2*nbdy)
    !       integer  mskflg (in) -  flag indicating the presence of a mask
    !       real     sfac   (in) -  additional scale factor
    !       real     offs   (in) -  additional offset
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    character(len=*) , intent(in) :: dims
    real             , intent(in) :: fld(*)
    integer          , intent(in) :: msk(*)
    integer          , intent(in) :: mskflg
    real             , intent(in) :: sfac
    real             , intent(in) :: offs

    character(len=4) :: c4
    integer :: i,j,ij,ijk,k,n,kd
    integer, parameter :: maxdm=5, ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
    real :: scf,ofs,arng(2),fldmin,fldmax
    logical :: uvflg
    integer, dimension(maxdm) :: start,count
    integer*2, allocatable, dimension(:,:,:) :: fldout,fld_out
    real, allocatable, dimension(:,:,:) ::rfld
    integer :: dimid,dimids(maxdm),strn,strind(2,maxdm)
    real, dimension(itdm,jtdm) :: rfldt

    ! --- Initialise fields
    uvflg = .false.
    kd = 1

    if(io_type  ==  1) then
#ifdef PNETCDF
      do i = 1,5
        istart(i) = 1
        icount(i) = 1
      end do


      ! --- define variable
      call ncsevl(dims,strn,strind)
      c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))

      if (c4 == 'time') then
        istart(strn) = rec
        icount(strn) = 1
      end if
      tkd = 1

      do n = 1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
        dimids(n) = dimid
        if ((n == strn.and.c4 /= 'time').or. &
            (n > 2.and.n == strn-1.and.c4 == 'time')) then
           call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        end if
        if (dims(1:1) == 'u'.or.dims(1:1) == 'v') uvflg = .true.
      end do

      kd = tkd
      istart(1) = 1
      istart(2) = j0+1
      if(mproc  ==  mpe_1(nproc) ) then
        icount(1) = itdm
        icount(2) = jj
        if(kd  >  1) icount(3) = kd
      else
        do n = 1,strn
          icount(n) = 0
        end do
      end if

      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      ! --- compute scale factor and offset
      fldmin = abs(fillr8)
      fldmax = -abs(fillr8)

      if (mskflg == 1) then
        !$omp parallel do reduction(min:fldmin) reduction(max:fldmax)
        !$omp private(j,i,ij,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij) == 1) then
                ijk = ij+(k-1)*ijdm
                fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
              end if
            end do
          end do
        end do
        !$omp end parallel do
      else if (mskflg == 2) then
        !$omp parallel do reduction(min:fldmin) reduction(max:fldmax)
        !$omp private(j,i,ij,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk = ij+(k-1)*ijdm
              if (msk(ij) == 1.and.fld(ijk) /= fillr8) then
                fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
              end if
            end do
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do reduction(min:fldmin) reduction(max:fldmax)
        !$omp private(j,i,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ijk = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
            end do
          end do
        end do
        !$omp end parallel do
      end if
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = max(abs(fldmax),abs(fldmin))/dble(i2max)
        end if
        ofs = 0.d0
      else
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = (fldmax-fldmin)/dble(2*i2max)
        end if
        ofs = 0.5*(fldmin+fldmax)
      end if
      arng(1) = fldmin
      arng(2) = fldmax

      ! --- Define attributes
      if (rec == 1) then
        call ncerro(nfmpi_redef(ncid))
        clen = 2
        call ncerro(nfmpi_put_att_double(ncid,rhid,'actual_range',nf_double,clen,arng))
        clen = 1
        call ncerro(nfmpi_put_att_double(ncid,rhid,'scale_factor',nf_double,clen,scf))
        call ncerro(nfmpi_put_att_double(ncid,rhid,'add_offset',nf_double,clen,ofs))
        call ncerro(nfmpi_enddef(ncid))
      end if
      allocate(rfld(ii,jj,kd))
      allocate(fldout(ii,jj,kd))
      allocate(fld_out(itdm,jj,kd))
      fld_out = i2fill

      ! --- Prepare and write output field
      if (mskflg == 1) then
        !$omp parallel do
        !$omp private(j,i,ij)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij) == 1) then
                rfld(i,j,k) = nint(((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf) - i2fill
              else
                rfld(i,j,k) = 0
              end if
            end do
          end do
        end do
        !$omp end parallel do
      else if (mskflg == 2) then
        !$omp parallel do
        !$omp private(j,i,ij,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk = ij+(k-1)*ijdm
              if (msk(ij) == 1.and.fld(ijk) /= fillr8) then
                rfld(i,j,k) = nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
              else
                rfld(i,j,k) = 0
              end if
            end do
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do
        !$omp private(j,i,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ijk = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              rfld(i,j,k) = nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            end do
          end do
        end do
        !$omp end parallel do
      end if

      !$omp parallel do private(j,i)
      do k = 1,kd
        do j = 1,jj
          do i = 1,ii
            fldout(i,j,k) = rfld(i,j,k)+i2fill
          end do
        end do
      end do
      !$omp end parallel do
      call xcgetrowint2(fld_out,fldout,kd)
      call ncerro(nfmpi_put_vara_int2_all(ncid,rhid,istart,icount,fld_out))

      if (mnproc == 1) write(lp,*) trim(vnm),arng
      deallocate(rfld,fldout,fld_out)
#endif
    else if(io_type  ==  0) then
      do n = 1,maxdm
        start(n) = 1
        count(n) = 1
      end do
      count(1) = itdm
      count(2) = jtdm

      ! --- define variable
      if (mnproc == 1) then
        call ncsevl(dims,strn,strind)
        c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4 == 'time') then
          start(strn) = rec
        else
          start(strn) = 1
        end if
        do n = 1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
          dimids(n) = dimid
          if ((n == strn.and.c4 /= 'time').or. &
              (n > 2.and.n == strn-1.and.c4 == 'time')) then
             call ncerro(nf90_inquire_dimension(ncid,dimid,len = kd))
          end if
          if (dims(1:1) == 'u'.or.dims(1:1) == 'v') uvflg = .true.
        end do
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      end if
      call xcbcst(uvflg)
      call xcbcst(kd)

      ! --- compute scale factor and offset
      fldmin = abs(fillr8)
      fldmax = -abs(fillr8)
      if (mskflg == 1) then
        !$omp parallel do reduction(min:fldmin) reduction(max:fldmax) &
        !$omp private(j,i,ij,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij) == 1) then
                ijk = ij+(k-1)*ijdm
                fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
              end if
            end do
          end do
        end do
        !$omp end parallel do
      else if (mskflg == 2) then
        !$omp parallel do reduction(min:fldmin) reduction(max:fldmax) &
        !$omp private(j,i,ij,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk = ij+(k-1)*ijdm
              if (msk(ij) == 1.and.fld(ijk) /= fillr8) then
                fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
              end if
            end do
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do reduction(min:fldmin) reduction(max:fldmax) &
        !$omp private(j,i,ijk)
        do k = 1,kd
          do j = 1,jj
            do i = 1,ii
              ijk = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
            end do
          end do
        end do
        !$omp end parallel do
      end if
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = max(abs(fldmax),abs(fldmin))/dble(i2max)
        end if
        ofs = 0.d0
      else
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = (fldmax-fldmin)/dble(2*i2max)
        end if
        ofs = 0.5*(fldmin+fldmax)
      end if
      arng(1) = fldmin
      arng(2) = fldmax

      ! --- Define attributes

      if (mnproc == 1 .and. rec == 1) then
        call ncerro(nf90_redef(ncid))
        call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
        call ncerro(nf90_put_att(ncid,rhid,'scale_factor',scf))
        call ncerro(nf90_put_att(ncid,rhid,'add_offset',ofs))
        call ncerro(nf90_enddef(ncid))
      end if
      allocate(rfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      allocate(fldout(itdm,jtdm,1))

      ! --- Prepare and write output field
      do k = 1,kd
        if (k > 1) start(3) = k
        if (mskflg == 1) then
          !$omp parallel do private(i,ij)
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij) == 1) then
                rfld(i,j,1) = nint(((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf) - i2fill
              else
                rfld(i,j,1) = 0
              end if
            end do
          end do
          !$omp end parallel do
        else if (mskflg == 2) then
          !$omp parallel do private(i,ij,ijk)
          do j = 1,jj
            do i = 1,ii
              ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk = ij+(k-1)*ijdm
              if (msk(ij) == 1.and.fld(ijk) /= fillr8) then
                rfld(i,j,1) = nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
              else
                rfld(i,j,1) = 0
              end if
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(i,ijk)
          do j = 1,jj
            do i = 1,ii
              ijk = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              rfld(i,j,1) = nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            end do
          end do
          !$omp end parallel do
        end if
        call xcaget(rfldt,rfld,1)
        if (mnproc == 1) then
          !$omp parallel do
          do j = 1,jtdm
            do i = 1,itdm
              fldout(i,j,1) = rfldt(i,j)+i2fill
            end do
          end do
          !$omp end parallel do
          call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
        end if
      end do

      ! --- Put file back to define mode
      if (mnproc == 1) then
        write(lp,*) trim(vnm),arng
      end if
      deallocate(rfld,fldout)

    end if
  end subroutine ncpack



  subroutine nccomp(vnm,dims,fld,msk,sfac,offs,prec)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes real field to nc-file in compressed format.
    !       The horizontal dimensions are replaced with a single, compressed
    !       dimension and only ocean points are written.
    ! --- Arguments:
    !       char(*)  vnm  (in) -  variable name
    !       char(*)  dims (in) -  axes string, e.g. 'pcomp z time'
    !       real(*)  fld  (in) -  field with dimension
    !                             (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
    !       int(*)   msk  (in) -  field mask with dimensions
    !                             (idm+2*nbdy)*(jdm+2*nbdy)
    !       real     sfac (in) -  additional scale factor
    !       real     offs (in) -  additional offset
    !       integer  prec (in) -  precision: 4=real4, 8=real8
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    character(len=*) , intent(in) :: dims
    real             , intent(in) :: fld(*)
    integer          , intent(in) :: msk(*)
    real             , intent(in) :: sfac
    real             , intent(in) :: offs
    integer          , intent(in) :: prec

    integer, parameter :: maxdm=5, ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
    integer :: dimid,dimids(maxdm),strn,strind(2,maxdm)
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
    real, dimension(itdm,jtdm) :: rfldt,rmskt
    real, dimension(itdm*jtdm) :: fldout
    real(kind = 4), dimension(itdm*jtdm) :: fldoutr4
    character(len=4) :: c4
    integer :: i,j,ij,ijk,k,n,kd
    integer, dimension(maxdm) :: start,count

    ! --- Initialise fields
    kd = 1
    if(io_type  ==  1) then
#ifdef PNETCDF
      do n = 1,maxdm
        istart(n) = 1
        icount(n) = 1
      end do

      ! --- define variable
      call ncsevl(dims,strn,strind)
      c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
      if (c4 == 'time') then
        istart(strn) = rec
        icount(strn) = 1
      end if
      tkd = 1
      do n = 1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
        dimids(n) = dimid
        if (n == 2.and.dims(strind(1,n):strind(2,n)) /= 'time') then
           call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        end if
      end do
      kd = tkd
      clen = 1
      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      if (mnproc /= 1) then
        do n = 1,maxdm
          icount(n) = 0
        end do
      end if
#endif
    else if(io_type  ==  0) then
      do n = 1,maxdm
        start(n) = 1
        count(n) = 1
      end do

      ! --- define variable
      if (mnproc == 1) then
        call ncsevl(dims,strn,strind)
        c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4 == 'time') then
          start(strn) = rec
        else
          start(strn) = 1
        end if
        do n = 1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
          dimids(n) = dimid
          if (n == 2.and.dims(strind(1,n):strind(2,n)) /= 'time') then
             call ncerro(nf90_inquire_dimension(ncid,dimid,len = kd))
          end if
        end do
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      end if
      call xcbcst(kd)
    end if

    ! --- get mask for the full domain
    !$omp parallel do private(i,ij)
    do j = 1,jj
      do i = 1,ii
        ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
        rmsk(i,j) = msk(ij)
      end do
    end do
    !$omp end parallel do
    call xcaget(rmskt,rmsk,1)

    ! --- Prepare output field
    do k = 1,kd
      !$omp parallel do private(i,ij,ijk)
      do j = 1,jj
        do i = 1,ii
          ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          ijk = ij+(k-1)*ijdm
          if (msk(ij) /= 0.and.fld(ijk) /= fillr8) then
            rfld(i,j) = fld(ijk)*sfac+offs
          else
            rfld(i,j) = fillr8
          end if
        end do
      end do
      !$omp end parallel do
      call xcaget(rfldt,rfld,1)
      if(io_type  ==  1) then
#ifdef PNETCDF
        n = 0
        if (prec /= 4) then
          if (mnproc == 1) then
            do j = 1,jtdm
              do i = 1,itdm
                if (rmskt(i,j) > .5) then
                  n = n+1
                  fldout(n) = rfldt(i,j)
                end if
              end do
            end do
            istart(2) = max(start(2),k)
            icount(1) = n
          end if
          call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,fldout))
        else
          if (mnproc == 1) then
            do j = 1,jtdm
              do i = 1,itdm
                if (rmskt(i,j) > .5) then
                  n = n+1
                  fldoutr4(n) = rfldt(i,j)
                end if
              end do
            end do
            istart(2) = max(start(2),k)
            icount(1) = n
          end if
          call ncerro(nfmpi_put_vara_real_all(ncid,rhid,istart,icount,fldoutr4))
        end if
#endif
      else if(io_type  ==  0) then
        if (mnproc == 1) then
          n = 0
          if (prec /= 4) then
            do j = 1,jtdm
              do i = 1,itdm
                if (rmskt(i,j) > .5) then
                  n = n+1
                  fldout(n) = rfldt(i,j)
                end if
              end do
            end do
            start(2) = max(start(2),k)
            count(1) = n
            call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
          else
            do j = 1,jtdm
              do i = 1,itdm
                if (rmskt(i,j) > .5) then
                  n = n+1
                  fldoutr4(n) = rfldt(i,j)
                end if
              end do
            end do
            start(2) = max(start(2),k)
            count(1) = n
            call ncerro(nf90_put_var(ncid,rhid,fldoutr4,start,count))
          end if
        end if
      end if
    end do
  end subroutine nccomp



  subroutine nccopa(vnm,dims,fld,msk,sfac,offs)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes real field to nc-file in compressed+packed format using
    !       scale factor and offset.
    !       The horizontal dimensions are replaced with a single, compressed
    !       dimension and only ocean points are written.
    ! --- Arguments:
    !       char(*)  vnm  (in) -  variable name
    !       char(*)  dims (in) -  axes string, e.g. 'pcomp z time'
    !       real(*)  fld  (in) -  field with dimension
    !                             (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
    !       int(*)   msk  (in) -  field mask with dimensions
    !                             (idm+2*nbdy)*(jdm+2*nbdy)
    !       real     sfac (in) -  additional scale factor
    !       real     offs (in) -  additional offset
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    character(len=*) , intent(in) :: dims
    real             , intent(in) :: fld(*)
    integer          , intent(in) :: msk(*)
    real             , intent(in) :: sfac
    real             , intent(in) :: offs

    integer :: i,j,ij,ijk,ijdm,k,kd,n
    integer, parameter :: maxdm = 5
    integer :: dimid,dimids(maxdm),strn,strind(2,maxdm)
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
    real, dimension(itdm,jtdm) :: rfldt,rmskt
    real :: scf,ofs,arng(2),fldmin,fldmax
    integer*2, dimension(itdm*jtdm) :: fldout
    logical :: uvflg
    integer, dimension(maxdm) :: start,count

    ! --- Initialise fields
    uvflg = .false.
    kd = 1

    if (io_type  ==  1) then
#ifdef PNETCDF
      do n = 1,maxdm
        istart(n) = 1
        icount(n) = 1
      end do
      call ncsevl(dims,strn,strind)
      if (dims(strind(1,strn):strind(2,strn)) == 'time') then
        istart(strn) = rec
        icount(strn) = 1
      end if
      tkd = 1
      do n = 1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
        dimids(n) = dimid
        if (n == 2.and.dims(strind(1,n):strind(2,n)) /= 'time') &
             call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        if (dims(1:1) == 'u'.or.dims(1:1) == 'v') uvflg = .true.
      end do
      kd = tkd

      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      ! --- compute scale factor and offset
      fldmin = abs(fillr8)
      fldmax = -abs(fillr8)

      ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
      !$omp parallel do private(j,i,ij,ijk) &
      !$omp reduction(min:fldmin) reduction(max:fldmax)
      do k = 1,kd
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
            if (msk(ij) == 1.and.fld(ijk) < fillr8) then
              fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
            end if
          end do
        end do
      end do
      !$omp end parallel do
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = max(abs(fldmax),abs(fldmin))/dble(i2max)
        end if
        ofs = 0.d0
      else
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = (fldmax-fldmin)/dble(2*i2max)
        end if
        ofs = 0.5*(fldmin+fldmax)
      end if
      arng(1) = fldmin
      arng(2) = fldmax

      ! --- Define attributes
      call ncerro(nfmpi_redef(ncid))
      clen = 2
      call ncerro(nfmpi_put_att_double(ncid,rhid,'actual_range', &
           nf_double,clen,arng))
      clen = 1
      call ncerro(nfmpi_put_att_double(ncid,rhid,'scale_factor', &
           nf_double,clen,scf))
      call ncerro(nfmpi_put_att_double(ncid,rhid,'add_offset',nf_double, &
           clen,ofs))
      call ncerro(nfmpi_enddef(ncid))

      if (mnproc /= 1) then
        do n = 1,maxdm
          icount(n) = 0
        end do
      end if

      ! --- get mask for the full domain
      !$omp parallel do private(i,ij)
      do j = 1,jj
        do i = 1,ii
          ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          rmsk(i,j) = msk(ij)
        end do
      end do
      !$omp end parallel do
      call xcaget(rmskt,rmsk,1)

      ! --- Prepare and write output field
      do k = 1,kd
        !$omp parallel do
        !$omp private(i,ij,ijk)
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk = ij+(k-1)*ijdm
            if (msk(ij) == 1.and.fld(ijk) < fillr8) then
              rfld(i,j) = nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            else
              rfld(i,j) = 0
            end if
          end do
        end do
        !$omp end parallel do
        call xcaget(rfldt,rfld,1)
        if (mnproc == 1) then
          n = 0
          do j = 1,jtdm
            do i = 1,itdm
              if (rmskt(i,j) > 0.5) then
                n = n+1
                fldout(n) = rfldt(i,j)+i2fill
              end if
            end do
          end do
          istart(2) = max(istart(2),k)
          icount(1) = n
        end if
        call ncerro(nfmpi_put_vara_int2_all(ncid,rhid,istart,icount,fldout))
      end do

      ! --- Put file back to define mode
      write(lp,*) trim(vnm),arng
#endif
    else if(io_type  ==  0) then
      do n = 1,maxdm
        start(n) = 1
        count(n) = 1
      end do
      ! --- define variable
      if (mnproc == 1) then
        call ncsevl(dims,strn,strind)
        if (dims(strind(1,strn):strind(2,strn)) == 'time') then
          start(strn) = rec
        else
          start(strn) = 1
        end if
        do n = 1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
          dimids(n) = dimid
          if (n == 2.and.dims(strind(1,n):strind(2,n)) /= 'time') then
             call ncerro(nf90_inquire_dimension(ncid,dimid,len = kd))
          end if
          if (dims(1:1) == 'u'.or.dims(1:1) == 'v') uvflg = .true.
        end do
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      end if
      call xcbcst(uvflg)
      call xcbcst(kd)

      ! --- compute scale factor and offset
      fldmin = abs(fillr8)
      fldmax = -abs(fillr8)

      ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
      !$omp parallel do private(j,i,ij,ijk) &
      !$omp reduction(min:fldmin) reduction(max:fldmax)
      do k = 1,kd
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
            if (msk(ij) == 1.and.fld(ijk) < fillr8) then
              fldmin = min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax = max(fldmax,(fld(ijk)*sfac)+offs)
            end if
          end do
        end do
      end do
      !$omp end parallel do
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = max(abs(fldmax),abs(fldmin))/dble(i2max)
        end if
        ofs = 0.d0
      else
        if (fldmin >= fldmax) then
          scf = 1.d0
        else
          scf = (fldmax-fldmin)/dble(2*i2max)
        end if
        ofs = 0.5*(fldmin+fldmax)
      end if
      arng(1) = fldmin
      arng(2) = fldmax

      ! --- Define attributes
      if (mnproc == 1) then
        call ncerro(nf90_redef(ncid))
        call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
        call ncerro(nf90_put_att(ncid,rhid,'scale_factor',scf))
        call ncerro(nf90_put_att(ncid,rhid,'add_offset',ofs))
        call ncerro(nf90_enddef(ncid))
      end if

      ! --- get mask for the full domain
      !$omp parallel do private(i,ij)
      do j = 1,jj
        do i = 1,ii
          ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          rmsk(i,j) = msk(ij)
        end do
      end do
      !$omp end parallel do
      call xcaget(rmskt,rmsk,1)

      ! --- Prepare and write output field
      do k = 1,kd
        !$omp parallel do private(i,ij,ijk)
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk = ij+(k-1)*ijdm
            if (msk(ij) == 1.and.fld(ijk) < fillr8) then
              rfld(i,j) = nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            else
              rfld(i,j) = 0
            end if
          end do
        end do
        !$omp end parallel do
        call xcaget(rfldt,rfld,1)
        if (mnproc == 1) then
          n = 0
          do j = 1,jtdm
            do i = 1,itdm
              if (rmskt(i,j) > 0.5) then
                n = n+1
                fldout(n) = rfldt(i,j)+i2fill
              end if
            end do
          end do
          start(2) = max(start(2),k)
          count(1) = n
          call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
        end if
      end do
    end if

  end subroutine nccopa



  subroutine ncwrtc(vnm,dims,fld)
    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes a string array to nc-file.
    ! --- Arguments:
    !       char(*)  vnm       (in) -  variable name
    !       char(*)  dims      (in) -  string with dimension names
    !       char(*)  fld(*)    (in) -  input field
    ! ----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: vnm
    character(len=*), intent(in) :: dims
    character(len=*), intent(in) :: fld(*)

    ! Local variables
    integer :: n
    integer, parameter :: maxdm = 5
    integer :: dimids(maxdm),strn,strind(2,maxdm)
#ifdef PNETCDF
    integer(kind = MPI_OFFSET_KIND), dimension(maxdm) :: start,count
    integer(kind = mpi_offset_kind) :: nsp,slenmaxp
#endif
    integer :: ns,slenmax

    if(io_type  ==  1) then
#ifdef PNETCDF
      call ncsevl(dims,strn,strind)
      if (strn > 2) then
        if (mnproc == 1) then
          write(lp,*) 'ncwrtc: number of dimensions has to be <=2'
          call flush(lp)
        end if
        call xcstop('(ncerro)')
        stop '(ncerro)'
      end if
      do n = 1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimids(n)))
      end do
      if (strn == 1) then
        nsp = 1
      else
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(1),slenmaxp))
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(2),nsp))
      end if

      ! --- write data
      if (mnproc == 1) then
        count(2) = 1
        count(1) = slenmaxp
        start(1) = 1
        start(2) = 1
      else
        count(2) = 0
        count(1) = 0
        start(1) = 1
        start(2) = 1
      end if
      ! --- - define variable
      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      ! --- - leave define mode

      do n = 1,nsp
        if (mnproc == 1) start(2) = n
        call ncerro(nfmpi_put_vara_text_all(ncid,rhid,start,count,fld(n)//'X'))
      end do
#endif
    else if(io_type  ==  0) then
      if (mnproc == 1) then
        ! --- - inquire dimensions
        call ncsevl(dims,strn,strind)
        if (strn > 2) then
          write(lp,*) 'ncwrtc: number of dimensions has to be <=2'
          call flush(lp)
          call xchalt('(ncerro)')
          stop '(ncerro)'
        end if
        do n = 1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimids(n)))
        end do
        if (strn == 1) then
          ns = 1
        else
          call ncerro(nf90_inquire_dimension(ncid,dimids(1),len = slenmax))
          call ncerro(nf90_inquire_dimension(ncid,dimids(2),len = ns))
        end if
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))
        ! --- - write data
        do n = 1,ns
          call ncerro(nf90_put_var(ncid,rhid,fld(n)(1:slenmax)//'X',(/1,n/),(/slenmax,1/)))
        end do
      end if
    end if

  end subroutine ncwrtc



  subroutine ncwrtr(vnm,dims,fld,msk,mskflg,sfac,offs,prec)
    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes field to nc-file as real4 or real8.
    ! --- Arguments:
    !       char(*)  vnm       (in) -  variable name
    !       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
    !       real(*)  fld       (in) -  field with dimension
    !                                  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
    !       int(*)   msk       (in) -  field mask with dimensions
    !                                  (idm+2*nbdy)*(jdm+2*nbdy)
    !       integer  mskflg    (in) -  set to 1 if mask is used, 0 else
    !       real     sfac      (in) -  additional scale factor
    !       real     offs      (in) -  additional offset
    !       integer  prec      (in) -  precision: 4=real4, 8=real8
    ! ----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: vnm
    character(len=*), intent(in) :: dims
    real,             intent(in) :: fld(*)
    integer,          intent(in) :: msk(*)
    integer,          intent(in) :: mskflg
    real,             intent(in) :: sfac
    real,             intent(in) :: offs
    integer,          intent(in) :: prec

    ! Local variables
    integer :: i,j,ij,ijk,k,kd,n
    integer, parameter :: maxdm=5
    integer, parameter :: ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
    integer :: dimid,dimids(maxdm),strn,strind(2,maxdm)
    real(kind = 4), allocatable, dimension(:,:,:) :: r4fldt
    real(kind = 4), allocatable, dimension(:,:,:) :: wr4fldt
    real, allocatable, dimension(:,:,:) :: rfld,wrfld
    real, allocatable, dimension(:,:) :: rmsk,rfldt,rmskt
    integer, dimension(maxdm) :: start,count
    real :: fldmin,fldmax
    character :: c4*4

    ! --- Initialise fields
    kd = 1
    if(io_type  ==  1) then
#ifdef PNETCDF
      do i = 1,5
        istart(i) = 1
        icount(i) = 1
      end do
      tkd = 1
      allocate(rmsk(ii,jj))
      clen = 1
      call ncsevl(dims,strn,strind)
      c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))

      if (c4  == 'time') then
        istart(strn) = rec
        icount(strn) = 1
      end if

      tkd = 1
      do n = 1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
        dimids(n) = dimid
        if (n > 2.and. &
          ((n == strn.and.c4 /= 'time') .or. (n == strn-1.and.c4 == 'time'))) then
           call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        end if
      end do

      kd = tkd

      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      istart(1) = 1
      istart(2) = j0+1

      if(mproc  ==  mpe_1(nproc) ) then
        icount(1) = itdm
        icount(2) = jj
        if (kd > 1) icount(3) = kd
      else
        icount(1) = 0
        icount(2) = 0
      end if


      ! --- get mask for the full domain
      if (mskflg == 1) then
        !$omp parallel do private(i,ij)
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            rmsk(i,j) = msk(ij)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = 1.
          end do
        end do
        !$omp end parallel do
      end if

      ! --- Prepare output field
      fldmin = abs(fillr8)
      fldmax = -abs(fillr8)
      allocate(rfld(ii,jj,kd))
      rfld = fillr8
      do k = 1,kd
        !$omp parallel do private(i,ij,ijk) &
        !$omp reduction(min:fldmin) reduction(max:fldmax)
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk = ij+(k-1)*ijdm
            rfld(i,j,k) = fld(ijk)
            if (msk(ij) /= 0.and.fld(ijk) /= fillr8) then
              fldmin = min(rfld(i,j,k),fldmin)
              fldmax = max(rfld(i,j,k),fldmax)
            end if
          end do
        end do
        !$omp end parallel do
      end do
      if (prec /= 4) then
        allocate(wrfld(itdm,jj,kd))
        wrfld = fillr8
        if (sfac /= 1..or.offs /= 0.) then
          !$omp parallel do private(j,i)
          do k = 1,kd
            do j = 1,jj
              do i = 1,ii
                if (rmsk(i,j) < 0.5.or.rfld(i,j,k) == fillr8) then
                  rfld(i,j,k) = fillr8
                else
                  rfld(i,j,k) = rfld(i,j,k)*sfac+offs
                end if
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(j,i)
          do k = 1,kd
            do j = 1,jj
              do i = 1,ii
                if (rmsk(i,j) < 0.5.or.rfld(i,j,k) == fillr8) &
                     rfld(i,j,k) = fillr8
              end do
            end do
          end do
          !$omp end parallel do
        end if
        call xcgetrow(wrfld, rfld, kd)
        call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,wrfld))
        deallocate(wrfld,rfld,rmsk)
      else
        allocate(wr4fldt(itdm,jj,kd))
        allocate(r4fldt(ii,jj,kd))
        r4fldt = fillr8
        wr4fldt = fillr8
        if (sfac /= 1..or.offs /= 0.) then
          !$omp parallel do private(j,i)
          do k = 1,kd
            do j = 1,jj
              do i = 1,ii
                if (rmsk(i,j) < 0.5.or.rfld(i,j,k) == fillr8) then
                  r4fldt(i,j,k) = fillr8
                else
                  r4fldt(i,j,k) = rfld(i,j,k)*sfac+offs
                end if
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(j,i)
          do k = 1,kd
            do j = 1,jj
              do i = 1,ii
                if (rmsk(i,j) < 0.5.or.rfld(i,j,k) == fillr8) then
                  r4fldt(i,j,k) = fillr8
                else
                  r4fldt(i,j,k) = rfld(i,j,k)
                end if
              end do
            end do
          end do
          !$omp end parallel do
        end if
        call xcgetrow4(wr4fldt, r4fldt, kd)
        call ncerro(nfmpi_put_vara_real_all(ncid,rhid,istart,icount,wr4fldt))
        deallocate(wr4fldt,r4fldt,rfld,rmsk)
      end if
#endif
    else if(io_type  ==  0) then
      do n = 1,maxdm
        start(n) = 1
        count(n) = 1
      end do
      count(1) = itdm
      count(2) = jtdm

      allocate(rfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      allocate(rmsk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(rmskt(itdm,jtdm))
      allocate(rfldt(itdm,jtdm))
      if(prec  ==  4)allocate(r4fldt(itdm,jtdm,1))
      if (mnproc == 1) then
        call ncsevl(dims,strn,strind)
        c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4 == 'time') then
          start(strn) = rec
        else
          start(strn) = 1
        end if
        do n = 1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
          dimids(n) = dimid
          if (n > 2.and. &
               ((n == strn.and.c4 /= 'time').or. &
               (n == strn-1.and.c4 == 'time'))) &
               call ncerro(nf90_inquire_dimension(ncid,dimid,len = kd))
        end do
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      end if
      call xcbcst(kd)

      ! --- get mask for the full domain
      if (mskflg == 1) then
        !$omp parallel do private(i,ij)
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            rmsk(i,j) = msk(ij)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = 1.
          end do
        end do
        !$omp end parallel do
      end if
      call xcaget(rmskt,rmsk,1)

      ! --- Prepare output field
      fldmin = abs(fillr8)
      fldmax = -abs(fillr8)

      do k = 1,kd
        if (k > 1) start(3) = k
        !$omp parallel do private(i,ij,ijk) &
        !$omp reduction(min:fldmin) reduction(max:fldmax)
        do j = 1,jj
          do i = 1,ii
            ij = i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk = ij+(k-1)*ijdm
            rfld(i,j,1) = fld(ijk)
            if (msk(ij) /= 0 .and. fld(ijk) /= fillr8) then
              fldmin = min(rfld(i,j,1),fldmin)
              fldmax = max(rfld(i,j,1),fldmax)
            end if
          end do
        end do
        !$omp end parallel do
        call xcaget(rfldt,rfld,1)
        if (mnproc == 1) then
          if (prec /= 4) then
            if (sfac /= 1..or.offs /= 0.) then
              !$omp parallel do private(i)
              do j = 1,jtdm
                do i = 1,itdm
                  if (rmskt(i,j) < 0.5.or.rfldt(i,j) == fillr8) then
                    rfldt(i,j) = fillr8
                  else
                    rfldt(i,j) = rfldt(i,j)*sfac+offs
                  end if
                end do
              end do
              !$omp end parallel do
            else
              !$omp parallel do private(i)
              do j = 1,jtdm
                do i = 1,itdm
                  if (rmskt(i,j) < 0.5.or.rfldt(i,j) == fillr8) &
                       rfldt(i,j) = fillr8
                end do
              end do
              !$omp end parallel do
            end if
            call ncerro(nf90_put_var(ncid,rhid,rfldt,start,count))
          else
            if (sfac /= 1..or.offs /= 0.) then
              !$omp parallel do private(i)
              do j = 1,jtdm
                do i = 1,itdm
                  if (rmskt(i,j) < 0.5.or.rfldt(i,j) == fillr8) then
                    r4fldt(i,j,1) = fillr8
                  else
                    r4fldt(i,j,1) = rfldt(i,j)*sfac+offs
                  end if
                end do
              end do
              !$omp end parallel do
            else
              !$omp parallel do private(i)
              do j = 1,jtdm
                do i = 1,itdm
                  if (rmskt(i,j) < 0.5.or.rfldt(i,j) == fillr8) then
                    r4fldt(i,j,1) = fillr8
                  else
                    r4fldt(i,j,1) = rfldt(i,j)
                  end if
                end do
              end do
              !$omp end parallel do
            end if
            call ncerro(nf90_put_var(ncid,rhid,r4fldt,start,count))
          end if
        end if
      end do
      if (prec /= 4) then
        deallocate(rmsk,rfld,rmskt,rfldt)
      else
        deallocate(r4fldt,rmsk,rfld,rmskt,rfldt)
      end if
    end if
    call xcmin(fldmin)
    call xcmax(fldmax)

    ! --- Put file back to define mode

    if (mnproc == 1 ) then
      if (sfac > 0.) then
        write(lp,*) trim(vnm),(fldmin*sfac+offs),(fldmax*sfac+offs)
      else
        write(lp,*) trim(vnm),(fldmax*sfac+offs),(fldmin*sfac+offs)
      end if
    end if
  end subroutine ncwrtr



  subroutine ncwrti(vnm,dims,fld,msk,mskflg)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes field to nc-file as int4.

    ! --- Arguments:
    !       char(*)    vnm       (in) -  variable name
    !       char(*)    dims      (in) -  axes string, e.g. 'pcomp z time'
    !       integer(*) fld       (in) -  field with dimension
    !                                   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
    !       integer(*) msk       (in) -  field mask with dimensions
    !                                   (idm+2*nbdy)*(jdm+2*nbdy)
    !       integer    mskflg    (in) -  set to 1 if mask is used, 0 else
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    character(len=*) , intent(in) :: dims
    integer          , intent(in) :: fld(*)
    integer          , intent(in) :: msk(*)
    integer          , intent(in) :: mskflg

    integer :: i,j,k,n,kd
    integer, parameter :: maxdm=5,ijdm = (idm+2*nbdy)*(jdm+2*nbdy)
    integer :: dimid,dimids(maxdm),strn,strind(2,maxdm)
    real, allocatable, dimension(:,:) :: rmsk,rfldt,rmskt
    real, allocatable, dimension(:,:,:) :: rfld,wrfld
    integer, allocatable, dimension(:,:,:) :: irfld
    integer, dimension(maxdm) :: start,count
    integer  :: fill
    character(len=4) :: c4

    ! --- Initialise fields
    kd = 1

    ! --- define variable
    if(io_type  ==  1) then
#ifdef PNETCDF
      do i = 1,5
        istart(i) = 1
        icount(i) = 1
      end do

      tkd = 1
      fill = nf_fill_int
      call ncsevl(dims,strn,strind)
      c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
      if (c4 == 'time') then
        istart(strn) = rec
        icount(strn) = 1
      end if
      tkd = 1
      do n = 1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
        dimids(n) = dimid
        if (n > 2.and. &
          ((n == strn.and.c4 /= 'time') .or. (n == strn-1.and.c4 == 'time'))) then
           call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        end if
      end do
      kd = tkd

      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      istart(1) = 1
      istart(2) = j0+1
      if(mproc  ==  mpe_1(nproc) ) then
        icount(1) = itdm
        icount(2) = jj
        if (kd > 1) icount(3) = kd
      else
        icount(1) = 0
        icount(2) = 0
      end if

      allocate(rmsk(ii,jj))
      allocate(rfld(ii,jj,kd))

      ! --- get mask for the full domain
      if (mskflg == 1) then
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = 1.
          end do
        end do
        !$omp end parallel do
      end if
      !$omp parallel do private(k,i)
      do j = 1,jj
        do k = 1,kd
          do i = 1,ii
            rfld(i,j,k) = fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(k,i)
      do j = 1,jj
        do k = 1,kd
          do i = 1,ii
            if (rmsk(i,j) < 0.5) rfld(i,j,k) = fill
          end do
        end do
      end do
      !$omp end parallel do

      allocate(wrfld(itdm,jj,kd))
      allocate(irfld(itdm,jj,kd))
      wrfld = fill
      call xcgetrow(wrfld, rfld, kd)
      !$omp parallel do private(k,i)
      do j = 1,jj
        do k = 1,kd
          do i = 1,itdm
            irfld(i,j,k) = wrfld(i,j,k)
          end do
        end do
      end do
      !$omp end parallel do

      call ncerro(nfmpi_put_vara_int_all(ncid,rhid,istart,icount,irfld))

      deallocate(wrfld,rfld,rmsk,irfld)
#endif

    else if(io_type  ==  0) then
      do n = 1,maxdm
        start(n) = 1
        count(n) = 1
      end do
      count(1) = itdm
      count(2) = jtdm
      allocate(rfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      allocate(rmsk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(rfldt(itdm,jtdm))
      allocate(rmskt(itdm,jtdm))
      fill = nf90_fill_int
      if (mnproc == 1) then
        call ncsevl(dims,strn,strind)
        c4 = dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4 == 'time') then
          start(strn) = rec
        else
          start(strn) = 1
        end if
        do n = 1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimid))
          dimids(n) = dimid
          if (n > 2.and. &
               ((n == strn.and.c4 /= 'time').or. &
               (n == strn-1.and.c4 == 'time'))) &
               call ncerro(nf90_inquire_dimension(ncid,dimid,len = kd))
        end do
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      end if
      call xcbcst(kd)

      ! --- get mask for the full domain
      if (mskflg == 1) then
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rmsk(i,j) = 1.
          end do
        end do
        !$omp end parallel do
      end if
      call xcaget(rmskt,rmsk,1)

      ! --- Prepare output field
      do k = 1,kd
        if (k > 1) start(3) = k
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            rfld(i,j,1) = fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
          end do
        end do
        !$omp end parallel do
        call xcaget(rfldt,rfld,1)
        if (mnproc == 1) then
          !$omp parallel do private(i)
          do j = 1,jtdm
            do i = 1,itdm
              if (rmskt(i,j) < 0.5) rfldt(i,j) = fill
            end do
          end do
          !$omp end parallel do
          call ncerro(nf90_put_var(ncid,rhid,rfldt,start,count))
        end if
      end do
    end if

  end subroutine ncwrti



  subroutine ncwrt1(vnm,dims,fld)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes real8 field from master process to nc-file.

    ! --- Arguments:
    !       char(*)  vnm       (in) -  variable name
    !       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
    !       real(*)  fld       (in) -  input field
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    character(len=*) , intent(in) :: dims
    real             , intent(in) :: fld(*)

    integer, parameter :: maxdm = 5
    integer :: n,strind(2,maxdm),ndims
#ifdef PNETCDF
    integer (kind = MPI_OFFSET_KIND), dimension(maxdm) :: dimlenp
#endif
    integer, dimension(maxdm) :: dimids
    integer, dimension(maxdm) :: start,count,dimlen


    if(io_type  ==  1) then
#ifdef PNETCDF
      call ncsevl(dims,ndims,strind)
      do n = 1,ndims
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimids(n)))
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(n),dimlenp(n)))
      end do
      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))
      do n = 1,maxdm
        istart(n) = 1
        icount(n) = 0
      end do
      if(mnproc  ==  1)then
        do n = 1,ndims-1
          icount(n) = dimlenp(n)
        end do

        if (dims(strind(1,ndims):strind(2,ndims)) == 'time') then
          istart(ndims) = rec
          icount(ndims) = 1
        else
          istart(ndims) = 1
          icount(ndims) = dimlenp(ndims)
        end if
      end if

      ! --- Write data to file
      call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,fld))

      ! --- Put file back to define mode
#endif
    else if(io_type  ==  0) then
      if (mnproc == 1) then

        ! --- - Define variable
        call ncsevl(dims,ndims,strind)
        do n = 1,ndims
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimids(n)))
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),len = dimlen(n)))
        end do
        call ncerro(nf90_inq_varid(ncid,vnm,rhid))

        ! --- - Prepare start and count
        do n = 1,maxdm
          start(n) = 1
          count(n) = 1
        end do
        do n = 1,ndims-1
          count(n) = dimlen(n)
        end do
        if (dims(strind(1,ndims):strind(2,ndims)) == 'time') then
          start(ndims) = rec
          count(ndims) = 1
        else
          start(ndims) = 1
          count(ndims) = dimlen(ndims)
        end if

        ! --- - Write data to file
        call ncerro(nf90_put_var(ncid,rhid,fld(1:product(count)),start,count))
      end if
    end if

  end subroutine ncwrt1



  ! ----------------------------------------------------------------------
  ! --- auxilary routines ------------------------------------------------
  ! ----------------------------------------------------------------------



  subroutine ncerro(ncstatus)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Display netcdf error massages

    ! --- Arguments:
    !       int ncstatus (in) -  netcdf status
    ! ----------------------------------------------------------------------

    integer :: ncstatus

    if(io_type  ==  1) then
#ifdef PNETCDF
      if (ncstatus /= nf_noerr) then
        write(lp,*) 'NetCDF error:',nfmpi_strerror(ncstatus)
        call flush(lp)
        call xchalt('(ncerro)')
        stop '(ncerro)'
      end if
#endif
    else if(io_type  ==  0) then
      if (ncstatus /= nf90_noerr) then
        write(lp,*) 'NetCDF error: ',nf90_strerror(ncstatus)
        call flush(lp)
        call xchalt('(ncerro)')
        stop '(ncerro)'
      end if
    end if

  end subroutine ncerro



  subroutine ncsevl(strg,strgn,strgind)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Finds the number and the locations of sub-strings.
    !       Valid deliminators are ' ', '-' and ':'.

    ! --- Arguments:
    !       char(*) strg   (in)  -  input string
    !       int strgn      (out) -  number of sub-strings
    !       int strgind(*) (out) -  start/end locations of sub-strings (the
    !                              dimension must at least equal to strgn*2)
    ! ----------------------------------------------------------------------

    character(len=*), intent(in)  :: strg
    integer,          intent(out) :: strgn
    integer,          intent(out) :: strgind(*)

    character :: charold,charnew
    integer :: i

    charold = ' '
    strgn = 0
    do i = 1,len(strg)
      charnew = strg(i:i)
      if ((charold == ' '.or.charold == '-'.or.charold == ':').and. &
          (charnew /= ' '.and.charnew /= '-'.and.charnew /= ':')) then
        strgn = strgn+1
        strgind(strgn) = i
      else if ((charnew == ' '.or.charnew == '-'.or.charnew == ':') .and. &
               (charold /= ' '.and.charold /= '-'.and.charold /= ':')) then
        strgn = strgn+1
        strgind(strgn) = i-1
      end if
      charold = charnew
    end do
    if (mod(strgn,2) == 1) then
      strgn = strgn+1
      strgind(strgn) = len(strg)
    end if
    strgn = int(strgn/2.)

  end subroutine ncsevl

  logical function ncinqa(nm)
    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Inquires if attribute exists.
    ! --- Arguments:
    !       char(*)       nm    (in)  -  attribute name
    !       logical   ncinqa   (out)  -  return value
    ! ----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) ::  nm
#ifdef PNETCDF
    integer(kind = mpi_offset_kind) :: attlen
#endif
    integer :: attlen1

    if(io_type  ==  1) then
#ifdef PNETCDF
      status = nfmpi_inq_attlen(ncid,nf_global,trim(nm),attlen)
      if (status /= nf_noerr) then
        ncinqa = .false.
      else
        if (attlen /= 0) then
          ncinqa = .true.
        else
          ncinqa = .false.
        end if
      end if
#endif
    else if(io_type  ==  0) then
      if(mnproc  ==  1) then
        status = nf90_inquire_attribute(ncid,nf90_global,trim(nm), &
             len = attlen1)
        if (status /= nf90_noerr) then

          ncinqa = .false.
        else
          if (attlen1 /= 0) then
            ncinqa = .true.
          else
            ncinqa = .false.
          end if
        end if
      end if
      call xcbcst(ncinqa)
    end if

  end function ncinqa


  logical function ncinqv(vnm)
    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Inquires if variable exists.
    ! --- Arguments:
    !       char(*)      vnm    (in)  -  variable name
    !       logical   ncinqv   (out)  -  return value
    ! ----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in) :: vnm
    if(io_type  ==  1) then
#ifdef PNETCDF
      status = nfmpi_inq_varid(ncid,trim(vnm),rhid)
      if (status /= nf_noerr) then
        ncinqv = .false.
      else
        ncinqv = .true.
      end if
#endif
    else if(io_type  ==  0) then
      if(mnproc  ==  1) then
        status = nf90_inq_varid(ncid,trim(vnm),rhid)
        if (status /= nf90_noerr) then
          ncinqv = .false.
        else
          ncinqv = .true.
        end if
      end if
      call xcbcst(ncinqv)
    end if

  end function ncinqv

  subroutine ncdefvar(vnm,dims,itype,iatt)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes real8 field from master process to nc-file.

    ! --- Arguments:
    !       char(*)  vnm       (in) -  variable name
    !       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
    !       real(*)  fld       (in) -  input field
    ! ----------------------------------------------------------------------

    character(len=*) , intent(in) :: vnm
    character(len=*) , intent(in) :: dims
    integer          , intent(in) :: itype
    integer          , intent(in) :: iatt

    integer, parameter :: maxdm = 5
    integer :: n,strind(2,maxdm),ndims
    integer, dimension(maxdm) :: dimids

    if (rec == 1) then
      call ncsevl(dims,ndims,strind)
      if( io_type  ==  1) then
#ifdef PNETCDF
        do n = 1,ndims
          call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimids(n)))
        end do

        call ncerro(nfmpi_def_var(ncid,vnm,itype,ndims,dimids,rhid))
        clen = 1
        if(iatt  ==  8) then
          call ncerro(nfmpi_put_att_double(ncid,rhid,'_FillValue', &
               nf_double,clen,fillr8))
        else if(iatt  ==  4) then
          call ncerro(nfmpi_put_att_real(ncid,rhid,'_FillValue', &
               nf_real,clen,fillr4))
        else if(iatt  ==  2) then
          call ncerro(nfmpi_put_att_int(ncid,rhid,'_FillValue', &
               nf_int,clen,nf_fill_int))
        end if
#endif
      else if( io_type  ==  0) then
        do n = 1,ndims
          call ncerro(nf90_inq_dimid(ncid, dims(strind(1,n):strind(2,n)),dimids(n)))
        end do
        call ncerro(nf90_def_var(ncid,vnm,itype,dimids(1:ndims),rhid))
        if(iatt  ==  8) then
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr8))
        else if(iatt  ==  4) then
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr4))
        else if(iatt  ==  2) then
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',nf90_fill_int))
        end if
      end if
    end if

  end subroutine ncdefvar

  subroutine ncdefvar3d(frmt,cmpflg,gridid, &
       vnm,vlngnm,vstdnm,vunits,isize)

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Writes real8 field from master process to nc-file.

    ! --- Arguments:
    !       char(*)  vnm       (in) -  variable name
    !       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
    !       real(*)  fld       (in) -  input field
    ! ----------------------------------------------------------------------

    integer, parameter :: maxdm = 5
    integer :: frmt,cmpflg,n,isize
    character(len = *) :: gridid,vnm,vlngnm,vstdnm,vunits

    character(len = 100) :: dims
    integer :: strind(2,maxdm),ndims
    integer, dimension(maxdm) :: dimids

    real :: arng(2)
    integer*2 i2min,vrng(2)
    if (rec == 1) then
      i2min = -i2max
      vrng(1) = i2min
      vrng(2) = i2max
      arng(1) = fillr8
      arng(2) = fillr8


      ! --- Check whether field should be written
      if (frmt == 0) return

      ! --- Create dimension string
      if(isize  ==  0) then
        if (cmpflg == 1) then
          dims = gridid(1:1)//'comp time'
        else
          dims = 'x y time'
        end if
      else if (isize  ==  1) then
        if (cmpflg == 1) then
          dims = gridid(1:1)//'comp sigma time'
        else
          dims = 'x y sigma time'
        end if
      else if (isize  ==  2) then
        if (cmpflg == 1) then
          dims = gridid(1:1)//'comp depth time'
        else
          dims = 'x y depth time'
        end if
      else if (isize  ==  3) then
        if (cmpflg == 1) then
          dims = 'pcomp ks time'
        else
          dims = 'x y ks time'
        end if
      else if (isize  ==  4) then
        if (cmpflg == 1) then
          dims = 'pcomp time'
        else
          dims = 'x y time'
        end if
      else
        write (lp,*) 'unknown data size!'
        call xchalt('(ncdefvar3d)')
        stop '(ncdefvar3d)'
      end if

      call ncsevl(dims,ndims,strind)
      if( io_type  ==  1) then
#ifdef PNETCDF
        do n = 1,ndims
          call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),dimids(n)))
        end do


        ! --- Check output format
        if (frmt == 2) then
          call ncerro(nfmpi_def_var(ncid,vnm,nf_int2,ndims,dimids,rhid))
          clen = 1
          call ncerro(nfmpi_put_att_int2(ncid,rhid,'_FillValue',nf_int2,clen,i2fill))
          clen = 2
          call ncerro(nfmpi_put_att_int2(ncid,rhid,'valid_range',nf_int2,clen,vrng))

          call ncerro(nfmpi_put_att_double(ncid,rhid,'actual_range',nf_double,clen,arng))
          clen = 1
          call ncerro(nfmpi_put_att_double(ncid,rhid,'scale_factor',nf_double,clen,fillr8))
          call ncerro(nfmpi_put_att_double(ncid,rhid,'add_offset',nf_double,clen,fillr8))

          if (cmpflg == 1) then
            clen = 3
            call ncerro(nfmpi_put_att_text(ncid,rhid,'compress',clen,'x y'))
          end if
        else if (frmt == 4) then
          if (cmpflg == 1) then
            call ncerro(nfmpi_def_var(ncid,vnm,nf_real,ndims,dimids,rhid))
            clen = 3
            call ncerro(nfmpi_put_att_text(ncid,rhid,'compress',clen,'x y'))
          else
            call ncerro(nfmpi_def_var(ncid,vnm,nf_real,ndims,dimids,rhid))
          end if
          clen = 1
          call ncerro(nfmpi_put_att_real(ncid,rhid,'_FillValue', &
               nf_real,clen,fillr4))
        else if (frmt == 8) then
          if (cmpflg == 1) then
            call ncerro(nfmpi_def_var(ncid,vnm,nf_double,ndims,dimids,rhid))
            clen = 3
            call ncerro(nfmpi_put_att_text(ncid,rhid,'compress',clen,'x y'))
          else
            call ncerro(nfmpi_def_var(ncid,vnm,nf_double,ndims,dimids,rhid))
          end if
          clen = 1
          call ncerro(nfmpi_put_att_double(ncid,rhid,'_FillValue',nf_double,clen,fillr8))
        else
          write (lp,*) 'unknown output format!'
          call xchalt('(ncdefvar3d)')
          stop '(ncdefvar3d)'
        end if
#endif
      else if( io_type  ==  0) then
        do n = 1,ndims
          call ncerro(nf90_inq_dimid(ncid, dims(strind(1,n):strind(2,n)),dimids(n)))
        end do

        if (frmt == 2) then
          call ncerro(nf90_def_var(ncid,vnm,nf90_int2,dimids(1:ndims),rhid))
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',i2fill))
          call ncerro(nf90_put_att(ncid,rhid,'valid_range',vrng))
          call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
          call ncerro(nf90_put_att(ncid,rhid,'scale_factor',fillr8))
          call ncerro(nf90_put_att(ncid,rhid,'add_offset',fillr8))
          if (cmpflg == 1) then
             call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))
          end if
        else if (frmt == 4) then
          call ncerro(nf90_def_var(ncid,vnm,nf90_real,dimids(1:ndims),rhid))
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr4))
          if (cmpflg == 1) then
             call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))
          end if
        else if (frmt == 8) then
          call ncerro(nf90_def_var(ncid,vnm,nf90_double, &
               dimids(1:ndims),rhid))
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr8))
          if (cmpflg == 1) then
             call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))
          end if
        else
          write (lp,*) 'unknown output format!'
          call xchalt('(ncdefvar3d)')
          stop '(ncdefvar3d)'
        end if

      end if

      ! --- Define attributes
      if (len(trim(vunits)) /= 0) call ncattr('units',vunits)
      if (len(trim(vlngnm)) /= 0) call ncattr('long_name',vlngnm)
      if (len(trim(vstdnm)) /= 0) call ncattr('standard_name',vstdnm)
      call ncattr('coordinates',gridid(1:1)//'lon '//gridid(1:1)//'lat')
      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
    end if
  end subroutine ncdefvar3d

  subroutine ncedef

    ! ----------------------------------------------------------------------
    ! --- Description:
    !       Closes NetCDF file
    ! ----------------------------------------------------------------------

    if(io_type  ==  1) then
#ifdef PNETCDF
      call ncerro(nfmpi_enddef(ncid))
#endif
    else if(io_type  ==  0) then
      if (mnproc == 1) then
        if (flgpad) then
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,81920,4,40960,4))
          flgpad = .true.
        end if

      end if
    end if
  end subroutine ncedef

end module mod_nctools
