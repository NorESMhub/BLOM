! ------------------------------------------------------------------------------
! Copyright (C) 2004-2020 Ingo Bethke, Mats Bentsen, Alok Kumar Gupta
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

      module mod_nctools
c ----------------------------------------------------------------------
c --- NetCDF tools -----------------------------------------------------
c ----------------------------------------------------------------------
c
c Description:
c    nctools.f is a boundle of fortran subroutines designed to write
c    output from the mpi version of BLOM to NetCDF data-format. 
c
c    Prerequisites are the linkage of the f90 netcdf library and the
c    access to the module netcdf.mod.
c
c Comments:
c    The order of actions is relevant: i.e. dimensions should be
c    declared iniflag, then fields not containing the unlimited (time)
c    dimension, and finally the fields containing the time dimension.
c    Attributes may be defined any time.
c
c Contents:
c
c    ncfopn        - creates a file or opens existing file for reading
c    ncfcls        - closes an opened nc-file
c
c    ncdims        - defines a simple axis-dimension
c    ncdimc        - defines a dimension for compressed data storage 
c                    and optionally stores the index information that 
c                    is necessary to de-compress data
c    nctime        - defines a time dimension and stores a time value
c
c    ncattr        - adds a text attribute
c
c    ncputr        - writes a skalar or vector in real format 
c    ncputi        - writes a skalar or vector in integer format  
c    ncgetr        - reads a skalar or vector in real format
c    ncgeti        - reads a skalar or vector in integer format 
c
c    ncread        - reads a 2d or 3d field 
c
c    ncwrtc        - writes a string array 
c    ncwrti        - writes a 2d or 3d field as int4 
c    ncwrtr        - writes a 2d or 3d field as real8 
c    ncpack        - writes a 2d or 3d field in packed format as int2 
c                    with scale factor and offset
c    nccomp        - writes a 2d or 3d field in compressed (skipping  
c                    land points) in real8 format 
c    nccopa        - writes a 2d or 3d field in compressed/packed format   
c                    as int2 with scale factor and offset  
c
c    ncerro        - displays error massage
c    ncsevl        - evaluates strings with deliminators ' ',':' and '-'
c    ncdnum        - computes a reference number for a given date
c    ncinqv        - inquires if variable exits 
c
c Revision history:
c    may2008       - switched from F77-API to F90-API
c    may2008       - removed precision argument from ncwrti and ncwrtr 
c    apr2008       - added header padding
c    apr2008       - added ncputr, ncputi, ncgetr, ncgeti 
c    mar2008       - added nccomp, nccopa, ncread, ncdimc  
c    mar2008       - added MPI support 
c
c Contact:
c    Ingo Bethke (ingo.bethke@nersc.no)
c
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
      implicit none
#ifdef PNETCDF
#  include <pnetcdf.inc>
#  include <mpif.h>
#endif
c
c 
      private :: ii,jj,kk,idm,jdm,kdm,nbdy
c
      interface ncputr 
      module procedure ncputrs,ncputrv 
      end interface ncputr
c
      interface ncputi 
      module procedure ncputis,ncputiv 
      end interface ncputi
c
      interface ncgetr 
      module procedure ncgetrs,ncgetrv 
      end interface ncgetr
c
      interface ncgeti 
      module procedure ncgetis,ncgetiv 
      end interface ncgeti
c
c --- Declare global variables  
      logical, private, save :: flgpad 
      integer, private, save :: ncid,rhid,status,rec,io_type
#ifdef PNETCDF
      integer(kind=MPI_OFFSET_KIND) clen,istart(5),icount(5),tkd
#endif
      integer*2, parameter :: i2fill=-32768,i2max=32767
      real, parameter :: fillr8=9.9692099683868690e+36
      real(kind=4), parameter :: fillr4=9.9692099683868690e+36
      integer ndouble,nchar,nfint
c
      contains             


      subroutine ncfopn(fnm,faccess,frmt,irec,iotype)
c
c ----------------------------------------------------------------------
c --- Description:
c       Opens a NetCDF file for reading or writing. 
c
c --- Arguments:
c       char(*) fnm      (in)      -  file name
c       char(*) faccess  (in)      -  'w' for write and 'r' for read 
c       char(*) frmt (in,optional) -  'CLASSIC','64BIT','HDF5'
c       integer irec               -  time step to be written 
c ----------------------------------------------------------------------
c
      implicit none
c
      integer :: irec
      character*(*), intent(in) :: fnm,faccess,frmt
c
      integer :: oldmode
      integer, parameter :: nf90__64bit_offset=512
      integer, parameter :: nf90__hdf5=4096
c
#ifdef PNETCDF
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
      integer*4 ,save :: info=MPI_INFO_NULL
      character(len=3) :: stripestr
      character(len=9) :: stripestr2
#endif
      integer ierr,iotype,testio

c
c --- Open file        
      IO_TYPE=iotype
      testio=0
      IF(IO_TYPE .eq. 1) then 
#ifdef PNETCDF
      testio=1 
      ndouble=nf_double
      nchar=nf_char
      nfint=nf_int
      rec=1
      write(stripestr,('(i3)')) 8 
      write(stripestr2,('(i9)')) 1024*1024
      if (faccess(1:1).eq.'r') then
        call ncerro(nfmpi_open(mpicomm,fnm,nf_nowrite,MPI_INFO_NULL,
     .              ncid))
        rec=irec
      elseif (irec.ne.1) then
        status=nfmpi_open(mpicomm,fnm,nf_write,INFO,ncid)
        if (status.ne.nf_noerr) then
         call mpi_info_create(info,ierr)
         call mpi_info_set(info,'romio_ds_read','disable',ierr)
         call mpi_info_set(info,'romio_ds_write','disable',ierr)
         call mpi_info_set(info,"striping_factor",stripestr,ierr)
         call mpi_info_set(info,"striping_unit",stripestr2,ierr)
          if (frmt(1:1).eq.'6') then
            call ncerro(nfmpi_create(mpicomm,fnm,
     .               IOR(nf_clobber,nf_64bit_offset),INFO,ncid))
          elseif (frmt(1:1).eq.'h'.or.frmt(1:1).eq.'H') then
            call ncerro(nfmpi_create(mpicomm,fnm,
     .               IOR(nf_clobber,nf_64bit_data),INFO,ncid))
          else
            call ncerro(nfmpi_create(mpicomm,fnm,
     +               IOR(NF_CLOBBER,nf_format_classic),INFO,ncid))
          endif
        else
          call ncerro(nfmpi_redef(ncid))
          rec=irec
        endif
      else
        call mpi_info_create(info,ierr)
        call mpi_info_set(info,'romio_ds_read','disable',ierr)
        call mpi_info_set(info,'romio_ds_write','disable',ierr)
        call mpi_info_set(info,'striping_factor',stripestr,ierr)
        call mpi_info_set(info,'striping_unit',stripestr2,ierr)
        call mpi_info_set(info,'striping_unit',stripestr2,ierr)
        call mpi_info_set(info,'nc_header_align_size',stripestr2,ierr)
        if (frmt(1:1).eq.'6') then
          call ncerro(nfmpi_create(mpicomm,fnm,
     .             IOR(nf_clobber,nf_64bit_offset),INFO,ncid))
        elseif (frmt(1:1).eq.'h'.or.frmt(1:1).eq.'H') then
          call ncerro(nfmpi_create(mpicomm,fnm,
     .             IOR(nf_clobber,nf_64bit_data),INFO,ncid))
        else
          call ncerro(nfmpi_create(mpicomm,fnm,
     .             IOR(nf_clobber,nf_format_classic),INFO,ncid))
        endif
      endif
c --- Initialise header padding
      if (rec.eq.1) then
        flgpad=.false.
      else
        flgpad=.true.
      endif
#endif
      if(testio .eq. 0) then
      write(lp,*) 'check iotype in namelist'
      call xchalt('(ncerro)')
      stop '(ncerro)'
      endif 
      elseif (io_type .eq. 0) then
      ndouble=nf90_double
      nchar=nf90_char
      nfint=nf90_int
      if (mnproc.eq.1) then 
c 
c --- - Open file        
        rec=1
        if (faccess(1:1).eq.'r') then 
          call ncerro(nf90_open(fnm,nf90_nowrite,ncid))
          rec=irec
        elseif (irec.ne.1) then
          status=nf90_open(fnm,nf90_write,ncid)
          if (status.ne.nf90_noerr) then
            if (frmt(1:1).eq.'6') then 
              call ncerro(nf90_create(fnm,OR(nf90_clobber,
     .                                       nf90__64bit_offset),ncid))
            elseif (frmt(1:1).eq.'h'.or.frmt(1:1).eq.'H') then
              call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__hdf5),
     .                                ncid))
            else 
              call ncerro(nf90_create(fnm,nf90_clobber,ncid))
            endif
          else
            call ncerro(nf90_redef(ncid))
            rec=irec
          endif              
          call ncerro(nf90_set_fill(ncid,nf90_nofill,oldmode))
        else 
          if (frmt(1:1).eq.'6') then 
            call ncerro(nf90_create(fnm,OR(nf90_clobber,
     .                                     nf90__64bit_offset),ncid))
          elseif (frmt(1:1).eq.'h'.or.frmt(1:1).eq.'H') then
            call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__hdf5),
     .                              ncid))
          else 
            call ncerro(nf90_create(fnm,nf90_clobber,ncid))
          endif
          call ncerro(nf90_set_fill(ncid,nf90_nofill,oldmode))
        endif
c
c --- - Set variable id to global (needed for global attributes)
        rhid=nf90_global
c
c --- - Initialise header padding
        if (rec.eq.1) then 
          flgpad=.false. 
        else
          flgpad=.true.
        endif
c
      endif 
      endif
c
      end subroutine ncfopn



      subroutine ncfcls
c
c ----------------------------------------------------------------------
c --- Description:
c       Closes NetCDF file
c ----------------------------------------------------------------------
c
      implicit none
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
      call ncerro(nfmpi_close(ncid))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1) then 
        call ncerro(nf90_close(ncid))
      endif 
      endif
c
      end subroutine ncfcls


      subroutine ncdims(dnm,dim)
c
c ----------------------------------------------------------------------
c --- Description:
c       Creates a simple dimension
c
c --- Arguments:
c       char(*) dnm  (in) -  name of the dimension
c       integer dim  (in) -  number of values on axis
c ----------------------------------------------------------------------
c
      implicit none
c
      integer :: dimid,dim
      character*(*) :: dnm
c
c --- define dimension
      if (io_type .eq. 1) then
#ifdef PNETCDF
      clen = dim
      if (rec.eq.1)
     .  call ncerro(nfmpi_def_dim(ncid,trim(dnm),clen,dimid))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1.and.rec.eq.1)  
     .  call ncerro(nf90_def_dim(ncid,dnm,dim,dimid))
      endif
c
      end subroutine ncdims



      subroutine ncdimc(dnm,msk,flg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Creates a dimension for compressed data storage and optionally
c       stores the index information that is neccessary to de-compress
c       data.
c
c --- Arguments:
c       char(*)      dnm (in) -  name of dimension, e.g. 'pcomp'
c       integer(...) msk (in) -  2d mask with dimensions
c                                (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) 
c       integer      flg (in) -  1 if index variable should be written 
c ----------------------------------------------------------------------
c
      implicit none
c
      integer :: n,i,j,dimid,flg 
      integer, dimension(itdm*jtdm) :: ivec    
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: msk   
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rmsk   
      real, dimension(itdm,jtdm) :: rmskt
      character*(*) dnm 
c 
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          rmsk(i,j)=msk(i,j)
        enddo
      enddo 
c$OMP END PARALLEL DO
      call xcaget(rmskt,rmsk,1)
      if (mnproc.eq.1.and.rec.eq.1) then 
        n=0
        if (flg.eq.1) then 
          do j=1,jtdm
            do i=1,itdm
              if (rmskt(i,j).gt.0.5) then
                n=n+1
                ivec(n)=i+(j-1)*itdm
              endif
            enddo 
          enddo  
        else 
          do j=1,jtdm
            do i=1,itdm
              if (rmskt(i,j).gt.0.5) n=n+1
            enddo 
          enddo  
        endif
        endif
        if (io_type .eq. 0) then
        if (mnproc.eq.1.and.rec.eq.1) then
        call ncerro(nf90_def_dim(ncid,dnm,n,dimid))
        if (flg.eq.1) then
          call ncerro(nf90_def_var(ncid,dnm,nf90_int,dimid,rhid))
          if (flgpad) then 
            call ncerro(nf90_enddef(ncid))
          else
            call ncerro(nf90_enddef(ncid,81920,4,40960,4))
            flgpad=.true.
          endif 
          call ncerro(nf90_put_var(ncid,rhid,ivec(1:n)))
          call ncerro(nf90_redef(ncid))
        endif 
        endif
      elseif (io_type .eq. 1) then
#ifdef PNETCDF
      call xcbcst(n)
      clen=n
      if (rec.eq.1) then
      call ncerro(nfmpi_def_dim(ncid,dnm,clen,dimid))
      if (flg.eq.1) then
        call ncerro(nfmpi_def_var(ncid,dnm,nf_int,1,dimid,rhid))
          call ncerro(nfmpi_enddef(ncid))
        istart(1)=1
        icount(1)=0
        if (mnproc.eq.1) icount(1)=n
        call ncerro(nfmpi_put_vara_int_all(ncid,rhid,istart,icount,
     .                                     ivec(1:n)))
        call ncerro(nfmpi_redef(ncid))
      endif
      endif
#endif
      endif
c
      end subroutine ncdimc



      subroutine ncattr(attname,attrib)
c
c ----------------------------------------------------------------------
c --- Description:
c       Adds a single text attribute to the variable that has been 
c       created last. If no variable has been created then a global 
c       attribute is added.
c
c --- Arguments:
c       char(*) attname (in) -  attribute name
c       char(*) attrib  (in) -  attribute text
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) attname,attrib
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
      clen=len(trim(attrib))
      if (rec.eq.1)
     .  call ncerro(nfmpi_put_att_text(ncid,rhid,attname,clen,attrib))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1.and.rec.eq.1)  
     .  call ncerro(nf90_put_att(ncid,rhid,attname,attrib))
      endif
c
      end subroutine ncattr



      subroutine ncputrs(vnm,rval)
c
c ----------------------------------------------------------------------
c --- Description:
c       Write a skalar as global attribute in real8 format. 
c
c --- Arguments:
c       char(*) vnm    (in) -  variable name
c       real rval      (in) -  real skalar
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) vnm
      real rval 
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
      clen=1
      if (rec.eq.1)
     .  call ncerro(nfmpi_put_att_double(ncid,nf_global,vnm,nf_double,
     .                                   clen,rval))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,rval))
      endif
c
      end subroutine ncputrs



      subroutine ncputrv(vnm,rval)
c
c ----------------------------------------------------------------------
c --- Description:
c       Write a vector as global attribute in real8 format. 
c
c --- Arguments:
c       char(*) vnm    (in) -  variable name
c       real rval(:)   (in) -  real vector
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) vnm
      real rval(:)
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
      clen = size(rval)
      if (rec.eq.1)
     .  call ncerro(nfmpi_put_att_double(ncid,nf_global,vnm,nf_double,
     .                                   clen,rval))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,rval))
      endif
c
      end subroutine ncputrv



      subroutine ncputis(vnm,ival)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes a skalar as global attribute in int4 format.
c
c --- Arguments:
c       char(*) vnm       (in) -  variable name
c       integer ival      (in) -  integer skalar 
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) vnm
      integer ival
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
      clen=1
      if (rec.eq.1)
     .  call ncerro(nfmpi_put_att_int(ncid,nf_global,vnm,nf_int,
     .                                clen,ival))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,ival))
      endif
c
      end subroutine ncputis



      subroutine ncputiv(vnm,ival)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes a vector as global attribute in int4 format.
c
c --- Arguments:
c       char(*) vnm       (in) -  variable name
c       integer ival(:)   (in) -  integer vector
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) vnm
      integer ival(:)
      if (io_type .eq. 1) then
#ifdef PNETCDF
      clen=size(ival)
      if (rec.eq.1)
     .  call ncerro(nfmpi_put_att_int(ncid,nf_global,vnm,nf_int,
     .                                clen,ival))
#endif
      elseif (io_type .eq. 0) then
c
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,ival))
      endif
c
      end subroutine ncputiv



      subroutine ncgetrs(vnm,rval)
c
c ----------------------------------------------------------------------
c --- Description:
c      Reads a skalar in real8 format.
c
c --- Arguments:
c       char(*) vnm  (in)  -  variable name
c       real    rval (out) -  real skalar
c ----------------------------------------------------------------------
c 
      implicit none
c
      character*(*) :: vnm
      real :: rval
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
        call ncerro(nfmpi_get_att_double(ncid,nf_global,vnm,rval))
#endif
      elseif (io_type .eq. 0) then
        if (mnproc.eq.1) 
     .    call ncerro(nf90_get_att(ncid,nf90_global,vnm,rval))
        call xcbcst(rval)
      endif
c
      end subroutine ncgetrs



      subroutine ncgetrv(vnm,rval)
c
c ----------------------------------------------------------------------
c --- Description:
c      Reads a vector in real8 format.
c
c --- Arguments:
c       char(*) vnm  (in)  -  variable name
c       real(:) rval (out) -  real vector
c ----------------------------------------------------------------------
c 
      implicit none
c
      character*(*) :: vnm
      real :: rval(:)
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
        call ncerro(nfmpi_get_att_double(ncid,nf_global,vnm,rval))
#endif
      elseif (io_type .eq. 0) then
        if (mnproc.eq.1) 
     .    call ncerro(nf90_get_att(ncid,nf90_global,vnm,rval))
        call xcbcst(rval)
      endif
c
      end subroutine ncgetrv



      subroutine ncgetis(vnm,ival)
c
c ----------------------------------------------------------------------
c --- Description:
c      Reads a skalar in integer format.
c 
c --- Arguments:
c       char(*)    vnm  (in)  -  variable name
c       integer    ival (out) -  integer skalar
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm
      integer :: ival 
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
        call ncerro(nfmpi_get_att_int(ncid,nf_global,vnm,ival))
#endif
      elseif (io_type .eq. 0) then
        if (mnproc.eq.1)
     .    call ncerro(nf90_get_att(ncid,nf90_global,vnm,ival))
        call xcbcst(ival)
      endif
c
      end subroutine ncgetis



      subroutine ncgetiv(vnm,ival)
c
c ----------------------------------------------------------------------
c --- Description:
c      Read a vector in integer format.
c 
c --- Arguments:
c       char(*)    vnm  (in)  -  variable name
c       integer(:) ival (out) -  integer vector
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm
      integer :: ival(:)
c
      if (io_type .eq. 1) then
#ifdef PNETCDF
        call ncerro(nfmpi_get_att_int(ncid,nf_global,vnm,ival))
#endif
      elseif (io_type .eq. 0) then
        if (mnproc.eq.1)
     .    call ncerro(nf90_get_att(ncid,nf90_global,vnm,ival))
        call xcbcst(ival)
      endif
c
      end subroutine ncgetiv



      subroutine nctime(datenum,calendar,units,startdate)
c
c ----------------------------------------------------------------------
c --- Description:
c       Create time dimension and store time value.
c       Valid calendars are '360_day', 'noleap' = '365_day',
c       'all_leap' = '366_day', 'julian', 'proleptic_gregorian' and
c       'standard' = 'mixed' = 'gregorian'.
c
c --- Arguments:
c       real    datenum    (in) -  time value
c       char(*) calendar   (in) -  calendar choice
c       char(*) units      (in) -  time units, 
c                                  e.g. 'hours since 0001-01-01'
c       char(*) startdatde (in) -  start date, e.g. '1997-01-01'
c
c --- Comment:
c       The script assumes that the precision of startdate is not
c	higher than the unit specified in the units string. E.g. if
c       units='years since 0001-01-01' and startdate='1997-06-01' then
c       the resulting offset will be rounded to 1996 yrs.
c ----------------------------------------------------------------------
c
      implicit none
c
      integer :: bstrn,sstrn,strind(2,10),days,by,bm,bd,bh,
     .   sy,sm,sd,sh
      character*(*) :: calendar,units,startdate
      real :: datenum,caloffset
      data by/1/,bm/1/,bd/1/,bh/0/,sy/1/,sm/1/,sd/1/,sh/0/
c

      if (io_type .eq. 1) then
#ifdef PNETCDF
c --- define time dimension
      istart(1)=1
      icount(1)=0
c
c --- define time coordinate variable
         call ncerro(nfmpi_inq_varid(ncid,'time',rhid))
c
      if (mnproc.eq.1) then
c --- - analyse units string
        call ncsevl(units,bstrn,strind)
        read(units(strind(1,3):strind(2,3)),*) by
        if (bstrn.ge.4) read(units(strind(1,4):strind(2,4)),*) bm
        if (bstrn.ge.5) read(units(strind(1,5):strind(2,5)),*) bd
        if (bstrn.ge.6) read(units(strind(1,6):strind(2,6)),*) bh
c
c --- - analyse startdate string
        call ncsevl(startdate,sstrn,strind)
        read(startdate(strind(1,1):strind(2,1)),*) sy
        if (sstrn.ge.2) read(startdate(strind(1,2):strind(2,2)),*) sm
        if (sstrn.ge.3) read(startdate(strind(1,3):strind(2,3)),*) sd
        if (sstrn.ge.4) read(startdate(strind(1,4):strind(2,4)),*) sh
c
c --- - calculate calendar offset
        if     (units(1:1).eq.'y') then
          caloffset=real(sy-by)
        elseif (units(1:1).eq.'m') then
          caloffset=real(sy-by)*12.+real(sm-bm)
        elseif (units(1:1).eq.'d') then
          call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
          caloffset=real(days)
        elseif (units(1:1).eq.'h') then
          call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
          caloffset=real(days)*24.+real(sh-bh)
        endif
        datenum=datenum+caloffset
        istart(1)=rec
        icount(1)=1
      endif
c
c --- Change to data mode and write time
      call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,
     .            datenum))
#endif
      elseif (io_type .eq. 0) then
      if (mnproc.eq.1) then 
          call ncerro(nf90_inq_varid(ncid,'time',rhid))
c
c --- - analyse units string
        call ncsevl(units,bstrn,strind)
        read(units(strind(1,3):strind(2,3)),*) by
        if (bstrn.ge.4) read(units(strind(1,4):strind(2,4)),*) bm
        if (bstrn.ge.5) read(units(strind(1,5):strind(2,5)),*) bd
        if (bstrn.ge.6) read(units(strind(1,6):strind(2,6)),*) bh
c
c --- - analyse startdate string
        call ncsevl(startdate,sstrn,strind)
        read(startdate(strind(1,1):strind(2,1)),*) sy
        if (sstrn.ge.2) read(startdate(strind(1,2):strind(2,2)),*) sm
        if (sstrn.ge.3) read(startdate(strind(1,3):strind(2,3)),*) sd
        if (sstrn.ge.4) read(startdate(strind(1,4):strind(2,4)),*) sh
c
c --- - calculate calendar offset
        if     (units(1:1).eq.'y') then
          caloffset=real(sy-by)
        elseif (units(1:1).eq.'m') then
          caloffset=real(sy-by)*12.+real(sm-bm)
        elseif (units(1:1).eq.'d') then
          call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
          caloffset=real(days)
        elseif (units(1:1).eq.'h') then
          call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
          caloffset=real(days)*24.+real(sh-bh)
        endif
        datenum=datenum+caloffset
        call ncerro(nf90_put_var(ncid,rhid,datenum,(/rec/)))
      endif
      endif 
c
      end subroutine nctime
 


      subroutine ncread(vnm,fld,msk,mskflg,fill) 
c
c ----------------------------------------------------------------------
c --- Description:
c       Reads 2d or 3d field.
c
c --- Arguments:
c       char(*)      vnm    (in)  -  variable name
c       real(...)    fld    (out) -  field with dimension
c                                    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,:)
c       integer(...) msk    (in)  -  field mask, if mskflg=1 then the 
c                                    dimensions assumed to be 
c                                    (idm+2*nbdy)*(jdm+2*nbdy)
c       integer      mskflg (in)  -  1 if mask exists  
c       real         fill   (in)  -  fill value to be set on land   
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm 
      character*100 :: dimname 
      integer :: i,j,ij,k,kd,n,ndm
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer ndims,dimids(maxdm),dimlen,strind(2,maxdm),
     .  msk(*),mskflg
#ifdef PNETCDF
      integer(kind=MPI_OFFSET_KIND) tdimlen 
#endif
      integer, dimension(maxdm) :: start,count
      real, allocatable, dimension(:,:,:) :: rfldt
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
      real, dimension(itdm,jtdm) :: rmskt
      real, dimension(itdm*jtdm) :: rfldtcmp
      real :: fill,fld(*),ofs,scf
      logical cmpflg
c
c --- Initialise fields
      cmpflg=.false.
      kd=1
      IF(IO_TYPE .eq. 1) then
#ifdef PNETCDF
      do n=1,maxdm
        istart(n)=1
        icount(n)=1
      enddo  

      status=nfmpi_inq_varid(ncid,vnm,rhid)
      if (status.ne.nf_noerr) then
        if (mnproc.eq.1) then
          write(lp,*) 'WARNING: Problems reading variable ',trim(vnm)
          call flush(lp)
        endif
        call ncerro(status)
      endif
      call ncerro(nfmpi_inq_varndims(ncid,rhid,ndims))
      call ncerro(nfmpi_inq_vardimid(ncid,rhid,dimids))
      ndm=1
      do n=1,ndims
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(n),tdimlen))
         
        call ncerro(nfmpi_inq_dimname(ncid,dimids(n),dimname))
        if (dimname(2:5).eq.'comp') cmpflg=.true.
        if ((n.eq.2.and.cmpflg).or.(n.eq.3.and..not.cmpflg)) kd=tdimlen
        ndm=ndm*tdimlen
      enddo
c
c --- Get attributes
      status=nfmpi_get_att_double(ncid,rhid,'scale_factor',scf)
      if (status.ne.nf_noerr) scf=1.
      status=nfmpi_get_att_double(ncid,rhid,'add_offset',ofs)
      if (status.ne.nf_noerr) ofs=0.
c
c --- get mask for the full domain  
      if (mskflg.eq.1) then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcaget(rmskt,rmsk,1)
      endif
c
c --- Get data
      if (cmpflg) then
        allocate(rfldt(itdm,jtdm,1))
        rfldt=0.0
        do k=1,kd
          if (mnproc.eq.1) then
            istart(2)=k
            icount(1)=nint(dble(ndm)/dble(kd))
          else
            icount(1)=0
            icount(2)=0
          endif
            call ncerro(nfmpi_get_vara_double_all(ncid,rhid,istart,
     .                                            icount,rfldtcmp))
            if (mnproc.eq.1) then
              n=0
              do j=1,jtdm
                do i=1,itdm
                  if (rmskt(i,j).gt.0.5) then
                    n=n+1
                    rfldt(i,j,1)=rfldtcmp(n)
                  endif
                enddo
              enddo
            endif
            call xcaput(rfldt,rfld,1)
            if (mskflg.eq.1) then
              if (scf.eq.1..and.ofs.eq.0) then
c$OMP PARALLEL DO PRIVATE(i,ij)
                do j=1,jj
                  do i=1,ii
                    ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                    if (msk(ij).eq.1) then
                      fld(ij+(k-1)*ijdm)=rfld(i,j)
                    else
                      fld(ij+(k-1)*ijdm)=fill
                    endif
                  enddo
                enddo
c$OMP END PARALLEL DO
              else
c$OMP PARALLEL DO PRIVATE(i,ij)
              do j=1,jj
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij).eq.1) then
                    fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs
                  else
                    fld(ij+(k-1)*ijdm)=fill
                  endif
                enddo
              enddo
c$OMP END PARALLEL DO
            endif
          else
            if (scf.eq.1..and.ofs.eq.0) then
c$OMP PARALLEL DO PRIVATE(i,ij)
              do j=1,jj
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs
                enddo
              enddo
c$OMP END PARALLEL DO
            else
c$OMP PARALLEL DO PRIVATE(i,ij)
              do j=1,jj
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm)=rfld(i,j)
                enddo
              enddo
c$OMP END PARALLEL DO
            endif
          endif
        enddo
      else
        istart(1)=i0+1
        istart(2)=j0+1
        istart(3)=1
        icount(1)=ii
        icount(2)=jj
        icount(3)=kd
        allocate(rfldt(ii,jj,kd))
        call ncerro(nfmpi_get_vara_double_all(ncid,rhid,istart,
     .                                        icount,rfldt))

        if (mskflg.eq.1) then
          if (scf.eq.1..and.ofs.eq.0) then
c$OMP PARALLEL DO PRIVATE(k,i,ij)
            do j=1,jj
              do k=1,kd
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij).eq.1) then
                    fld(ij+(k-1)*ijdm)=rfldt(i,j,k)
                  else
                    fld(ij+(k-1)*ijdm)=fill
                  endif
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(k,i,ij)
            do j=1,jj
              do k=1,kd
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if (msk(ij).eq.1) then
                    fld(ij+(k-1)*ijdm)=rfldt(i,j,k)*scf+ofs
                  else
                    fld(ij+(k-1)*ijdm)=fill
                  endif
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif
        else
          if (scf.eq.1..and.ofs.eq.0) then
c$OMP PARALLEL DO PRIVATE(k,i,ij)
            do j=1,jj
              do k=1,kd
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm)=rfldt(i,j,k)*scf+ofs
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(k,i,ij)
            do j=1,jj
              do k=1,kd
                do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fld(ij+(k-1)*ijdm)=rfldt(i,j,k)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif
        endif
      endif
      deallocate(rfldt)
#endif
      else IF(IO_TYPE .eq. 0) then
      do n=1,maxdm
        start(n)=1
        count(n)=1
      enddo
c       
c --- Inquire dimensions, ect.
      if (mnproc.eq.1) then      
        status=nf90_inq_varid(ncid,vnm,rhid)
        if (status.ne.nf90_noerr) then 
          write(lp,*) 'WARNING: Problems reading variable ',trim(vnm)
          call flush(lp)
          call ncerro(status)
        endif
        call ncerro(nf90_inquire_variable(ncid,rhid,ndims=ndims))
        call ncerro(nf90_inquire_variable(ncid,rhid,dimids=dimids))
        ndm=1
        do n=1,ndims
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),len=dimlen))
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),
     .      name=dimname))
          if (dimname(2:5).eq.'comp') cmpflg=.true.
          if ((n.eq.2.and.cmpflg).or.(n.eq.3.and..not.cmpflg)) kd=dimlen
          ndm=ndm*dimlen
        enddo
      endif 
      call xcbcst(cmpflg)
      call xcbcst(kd)
c
c --- Get attributes
      if (mnproc.eq.1) then      
        status=nf90_get_att(ncid,rhid,'scale_factor',scf)
        if (status.ne.nf90_noerr) scf=1.
        status=nf90_get_att(ncid,rhid,'add_offset',ofs)
        if (status.ne.nf90_noerr) ofs=0.
      endif 
      call xcbcst(scf)
      call xcbcst(ofs)
c
c --- get mask for the full domain  
      if (mskflg.eq.1) then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcaget(rmskt,rmsk,1)
      endif
      allocate(rfldt(itdm,jtdm,1))
c
c --- Get data
      do k=1,kd
        if (mnproc.eq.1) then 
          if (cmpflg) then 
            start(2)=k
            count(1)=nint(dble(ndm)/dble(kd)) 
          else 
            start(3)=k 
            count(1)=itdm
            count(2)=jtdm
          endif
          if (cmpflg) then 
            call ncerro(nf90_get_var(ncid,rhid,rfldtcmp,start,count))
            n=0 
            do j=1,jtdm
              do i=1,itdm
                if (rmskt(i,j).gt.0.5) then
                  n=n+1
                  rfldt(i,j,1)=rfldtcmp(n) 
                endif 
              enddo
            enddo  
          else
            call ncerro(nf90_get_var(ncid,rhid,rfldt,start,count))
          endif  
        endif
        call xcaput(rfldt,rfld,1)
        if (mskflg.eq.1) then 
          if (scf.eq.1..and.ofs.eq.0) then 
c$OMP PARALLEL DO PRIVATE(i,ij)
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                if (msk(ij).eq.1) then
                  fld(ij+(k-1)*ijdm)=rfld(i,j)
                else
                  fld(ij+(k-1)*ijdm)=fill
                endif
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(i,ij)
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                if (msk(ij).eq.1) then
                  fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs
                else
                  fld(ij+(k-1)*ijdm)=fill
                endif
              enddo
            enddo
c$OMP END PARALLEL DO
          endif  
        else
          if (scf.eq.1..and.ofs.eq.0) then
c$OMP PARALLEL DO PRIVATE(i,ij)
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs 
              enddo
            enddo     
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(i,ij)
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                fld(ij+(k-1)*ijdm)=rfld(i,j)
              enddo
            enddo
c$OMP END PARALLEL DO
          endif
        endif        
      enddo 
      deallocate(rfldt)
      endif
c
      end subroutine ncread


 
      subroutine ncpack(vnm,dims,fld,msk,mskflg,sfac,offs)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real field to nc-file in packed format using scale factor
c       and offset.
c
c --- Arguments:
c       char(*)  vnm    (in) -  variable name
c       char(*)  dims   (in) -  axes string, e.g. 'x y z time'
c       real(*)  fld    (in) -  field with dimension
c       int(*)   msk    (in) -  field mask, if mskflg=1 the dimension
c                               is assumed (idm+2*nbdy)*(jdm+2*nbdy)
c       integer  mskflg (in) -  flag indicating the presence of a mask
c       real     sfac   (in) -  additional scale factor
c       real     offs   (in) -  additional offset
c ----------------------------------------------------------------------
c
      implicit none
c
      character*4 :: c4
      character*(*) :: vnm,dims
      integer :: i,j,ij,ijk,k,n,kd
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      real scf,sfac,ofs,offs,arng(2),fldmin,fldmax,fld(*)
      logical uvflg
      integer, dimension(maxdm) :: start,count
      integer*2, allocatable, dimension(:,:,:) :: fldout,fld_out
      real, allocatable, dimension(:,:,:) ::rfld
      integer dimid,dimids(maxdm),strn,strind(2,maxdm),msk(*),mskflg    
      real, dimension(itdm,jtdm) :: rfldt
c
c --- Initialise fields
      uvflg=.false.
      kd=1

      IF(IO_TYPE .eq. 1) then
#ifdef PNETCDF
      do i=1,5
      istart(i)=1
      icount(i)=1
      enddo

c
c --- define variable
      call ncsevl(dims,strn,strind)
      c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))

      if (c4.eq.'time') then
        istart(strn)=rec
        icount(strn)=1
      endif
      tkd=1 

      do n=1,strn
        call ncerro(nfmpi_inq_dimid(ncid,
     .    dims(strind(1,n):strind(2,n)),dimid))
        dimids(n)=dimid
        if ((n.eq.strn.and.c4.ne.'time').or.
     .      (n.gt.2.and.n.eq.strn-1.and.c4.eq.'time'))
     .    call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        if (dims(1:1).eq.'u'.or.dims(1:1).eq.'v') uvflg=.true.
      enddo

      kd=tkd
      istart(1) = 1
      istart(2) = j0+1
      if(mproc .eq. mpe_1(nproc) ) then
      icount(1) = itdm
      icount(2) = jj
      if(kd .gt. 1) icount(3)=kd
      else
      do n=1,strn
      icount(n)=0
      enddo
      endif

      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))
c
c --- compute scale factor and offset
      fldmin=abs(fillr8)
      fldmax=-abs(fillr8)
      
      if (mskflg.eq.1) then
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  PRIVATE(j,i,ij,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij).eq.1) then
                ijk=ij+(k-1)*ijdm
                fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (mskflg.eq.2) then
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  PRIVATE(j,i,ij,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk=ij+(k-1)*ijdm
              if (msk(ij).eq.1.and.fld(ijk).ne.fillr8) then
                fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+   PRIVATE(j,i,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
        endif
        ofs=0.d0
      else
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=(fldmax-fldmin)/dble(2*i2max)
        endif
        ofs=0.5*(fldmin+fldmax)
      endif
      arng(1)=fldmin
      arng(2)=fldmax
c
c --- Define attributes
      if (rec.eq.1) then
        call ncerro(nfmpi_redef(ncid))
        clen=2
        call ncerro(nfmpi_put_att_double(ncid,rhid,'actual_range',
     .                                   nf_double,clen,arng))
        clen=1
        call ncerro(nfmpi_put_att_double(ncid,rhid,'scale_factor',
     .                                   nf_double,clen,scf))
        call ncerro(nfmpi_put_att_double(ncid,rhid,'add_offset',
     .                                   nf_double,clen,ofs))
        call ncerro(nfmpi_enddef(ncid))
      endif
      allocate(rfld(ii,jj,kd))
      allocate(fldout(ii,jj,kd))
      allocate(fld_out(itdm,jj,kd))
      fld_out=i2fill
c
c --- Prepare and write output field
      if (mskflg.eq.1) then
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(j,i,ij) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij).eq.1) then
              rfld(i,j,k)=nint(((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf)
     .                   -i2fill
              else
                rfld(i,j,k)=0
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (mskflg.eq.2) then
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(j,i,ij,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk=ij+(k-1)*ijdm
              if (msk(ij).eq.1.and.fld(ijk).ne.fillr8) then
                rfld(i,j,k)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
              else
                rfld(i,j,k)=0
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(j,i,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              rfld(i,j,k)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c$OMP PARALLEL DO PRIVATE(j,i)
        do k=1,kd
          do j=1,jj
            do i=1,ii
              fldout(i,j,k)=rfld(i,j,k)+i2fill
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcgetrowint2(fld_out,fldout,kd)
          call ncerro(nfmpi_put_vara_int2_all(ncid,rhid,istart,icount,
     .                                        fld_out))
c
        if (mnproc.eq.1) write(lp,*) trim(vnm),arng
        deallocate(rfld,fldout,fld_out)
c
#endif
      ELSEIF(IO_TYPE .eq. 0) then
      do n=1,maxdm
        start(n)=1
        count(n)=1
      enddo
      count(1)=itdm
      count(2)=jtdm
c
c --- define variable 
      if (mnproc.eq.1) then
        call ncsevl(dims,strn,strind)
        c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4.eq.'time') then 
          start(strn)=rec 
        else 
          start(strn)=1
        endif 
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,
     .    dims(strind(1,n):strind(2,n)),dimid))
          dimids(n)=dimid
          if ((n.eq.strn.and.c4.ne.'time').or.
     .        (n.gt.2.and.n.eq.strn-1.and.c4.eq.'time')) 
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
          if (dims(1:1).eq.'u'.or.dims(1:1).eq.'v') uvflg=.true.
        enddo
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      endif
      call xcbcst(uvflg)
      call xcbcst(kd)
c
c --- compute scale factor and offset
      fldmin=abs(fillr8)
      fldmax=-abs(fillr8)
      if (mskflg.eq.1) then 
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  PRIVATE(j,i,ij,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij).eq.1) then
                ijk=ij+(k-1)*ijdm
                fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (mskflg.eq.2) then 
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  PRIVATE(j,i,ij,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk=ij+(k-1)*ijdm
              if (msk(ij).eq.1.and.fld(ijk).ne.fillr8) then
                fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+   PRIVATE(j,i,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then 
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
        endif
        ofs=0.d0
      else
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=(fldmax-fldmin)/dble(2*i2max)
        endif
        ofs=0.5*(fldmin+fldmax)
      endif 
      arng(1)=fldmin
      arng(2)=fldmax
c
c --- Define attributes

      if (mnproc.eq.1 .and. rec.eq.1) then 
          call ncerro(nf90_redef(ncid))
          call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
          call ncerro(nf90_put_att(ncid,rhid,'scale_factor',scf))
          call ncerro(nf90_put_att(ncid,rhid,'add_offset',ofs))
          call ncerro(nf90_enddef(ncid))
      endif
      allocate(rfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      allocate(fldout(itdm,jtdm,1))
c
c --- Prepare and write output field
      do k=1,kd
        if (k.gt.1) start(3)=k
        if (mskflg.eq.1) then
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(i,ij) 
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij).eq.1) then
              rfld(i,j,1)=nint(((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf)
     .                   -i2fill
              else
                rfld(i,j,1)=0 
              endif
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (mskflg.eq.2) then
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(i,ij,ijk) 
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk=ij+(k-1)*ijdm
              if (msk(ij).eq.1.and.fld(ijk).ne.fillr8) then
                rfld(i,j,1)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
              else
                rfld(i,j,1)=0
              endif
            enddo
          enddo
c$OMP END PARALLEL DO
        else
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(i,ijk) 
          do j=1,jj
            do i=1,ii 
              ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              rfld(i,j,1)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            enddo 
          enddo   
c$OMP END PARALLEL DO
        endif  
        call xcaget(rfldt,rfld,1)
        if (mnproc.eq.1) then
c$OMP PARALLEL DO 
          do j=1,jtdm
            do i=1,itdm
              fldout(i,j,1)=rfldt(i,j)+i2fill
            enddo 
          enddo 
c$OMP END PARALLEL DO
          call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
        endif 
      enddo 
c
c --- Put file back to define mode
      if (mnproc.eq.1) then
        write(lp,*) trim(vnm),arng
      endif
      deallocate(rfld,fldout)
c
      endif
      end subroutine ncpack



      subroutine nccomp(vnm,dims,fld,msk,sfac,offs,prec)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real field to nc-file in compressed format.
c       The horizontal dimensions are replaced with a single, compressed
c       dimension and only ocean points are written.  
c
c --- Arguments:
c       char(*)  vnm  (in) -  variable name
c       char(*)  dims (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld  (in) -  field with dimension
c                             (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       int(*)   msk  (in) -  field mask with dimensions 
c                             (idm+2*nbdy)*(jdm+2*nbdy)  
c       real     sfac (in) -  additional scale factor
c       real     offs (in) -  additional offset
c       integer  prec (in) -  precision: 4=real4, 8=real8  
c ----------------------------------------------------------------------
c
      implicit none
c
      character*4 :: c4
      character*(*) :: vnm,dims
      integer :: i,j,ij,ijk,k,n,prec,kd
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),strn,strind(2,maxdm),
     .  msk(*)     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real, dimension(itdm*jtdm) :: fldout
      real(kind=4), dimension(itdm*jtdm) :: fldoutr4
      real :: fld(*),sfac,offs
c 
c --- Initialise fields
      kd=1
      IF(IO_TYPE .eq. 1) then
#ifdef PNETCDF
      do n=1,maxdm
        istart(n)=1
        icount(n)=1
      enddo
c  
c --- define variable
      call ncsevl(dims,strn,strind)
      c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
      if (c4.eq.'time') then
        istart(strn)=rec
        icount(strn)=1
      endif
      tkd=1 
      do n=1,strn
       call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                             dimid))
        dimids(n)=dimid
        if (n.eq.2.and.dims(strind(1,n):strind(2,n)).ne.'time')
     .    call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
      enddo
      kd=tkd
      clen=1
      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      if (mnproc.ne.1) then
        do n=1,maxdm
          icount(n)=0
        enddo
      endif
#endif
      ELSEIF(IO_TYPE .eq. 0) then
      do n=1,maxdm
        start(n)=1
        count(n)=1
      enddo
c  
c --- define variable 
      if (mnproc.eq.1) then
        call ncsevl(dims,strn,strind)
        c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4.eq.'time') then
          start(strn)=rec
        else
          start(strn)=1
        endif
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                               dimid))
          dimids(n)=dimid
          if (n.eq.2.and.dims(strind(1,n):strind(2,n)).ne.'time') 
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
        enddo
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      endif
      call xcbcst(kd)
      endif
c
c --- get mask for the full domain
c$OMP PARALLEL DO PRIVATE(i,ij)
      do j=1,jj
        do i=1,ii
          ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          rmsk(i,j)=msk(ij)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
c$OMP PARALLEL DO PRIVATE(i,ij,ijk)
        do j=1,jj
          do i=1,ii 
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            if (msk(ij).ne.0.and.fld(ijk).ne.fillr8) then 
              rfld(i,j)=fld(ijk)*sfac+offs
            else
              rfld(i,j)=fillr8
            endif
          enddo 
        enddo       
c$OMP END PARALLEL DO
        call xcaget(rfldt,rfld,1) 
       IF(IO_TYPE .eq. 1) then
#ifdef PNETCDF
        n=0
        if (prec.ne.4) then
          if (mnproc.eq.1) then
            do j=1,jtdm
              do i=1,itdm
                if (rmskt(i,j).gt..5) then
                  n=n+1
                  fldout(n)=rfldt(i,j)
                endif
              enddo
            enddo
            istart(2)=max(start(2),k)
            icount(1)=n
          endif
          call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,
     .                                          fldout))
        else
          if (mnproc.eq.1) then
            do j=1,jtdm
              do i=1,itdm
                if (rmskt(i,j).gt..5) then
                  n=n+1
                  fldoutr4(n)=rfldt(i,j)
                endif
              enddo
            enddo
            istart(2)=max(start(2),k)
            icount(1)=n
          endif
          call ncerro(nfmpi_put_vara_real_all(ncid,rhid,istart,icount,
     .                                        fldoutr4))
        endif
#endif
        ELSEIF(IO_TYPE .eq. 0) then
        if (mnproc.eq.1) then
          n=0
          if (prec.ne.4) then 
            do j=1,jtdm
              do i=1,itdm
                if (rmskt(i,j).gt..5) then
                  n=n+1
                  fldout(n)=rfldt(i,j)
                endif 
              enddo 
            enddo 
            start(2)=max(start(2),k)
            count(1)=n
            call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
          else
            do j=1,jtdm
              do i=1,itdm
                if (rmskt(i,j).gt..5) then
                  n=n+1
                  fldoutr4(n)=rfldt(i,j)
                endif 
              enddo 
            enddo 
            start(2)=max(start(2),k)
            count(1)=n
            call ncerro(nf90_put_var(ncid,rhid,fldoutr4,start,count))
          endif  
        endif 
      endif
      enddo
      end subroutine nccomp



      subroutine nccopa(vnm,dims,fld,msk,sfac,offs)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real field to nc-file in compressed+packed format using 
c       scale factor and offset.
c       The horizontal dimensions are replaced with a single, compressed
c       dimension and only ocean points are written.  
c
c --- Arguments:
c       char(*)  vnm  (in) -  variable name
c       char(*)  dims (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld  (in) -  field with dimension
c                             (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       int(*)   msk  (in) -  field mask with dimensions 
c                             (idm+2*nbdy)*(jdm+2*nbdy)  
c       real     sfac (in) -  additional scale factor
c       real     offs (in) -  additional offset
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm,dims
      integer :: i,j,ij,ijk,ijdm,k,kd,n
      integer, parameter :: maxdm=5
      integer dimid,dimids(maxdm),strn,strind(2,maxdm),msk(*)     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real scf,sfac,ofs,offs,arng(2),fldmin,fldmax,fld(*)
      integer*2, dimension(itdm*jtdm) :: fldout
      logical uvflg
      integer, dimension(maxdm) :: start,count
c
c --- Initialise fields
      uvflg=.false.
      kd=1
c 
      if (io_type .eq. 1) then
#ifdef PNETCDF
      do n=1,maxdm
        istart(n)=1
        icount(n)=1
      enddo
      call ncsevl(dims,strn,strind)
      if (dims(strind(1,strn):strind(2,strn)).eq.'time') then
        istart(strn)=rec
        icount(strn)=1
      endif
      tkd=1
      do n=1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                              dimid))
        dimids(n)=dimid
        if (n.eq.2.and.dims(strind(1,n):strind(2,n)).ne.'time')
     .    call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
        if (dims(1:1).eq.'u'.or.dims(1:1).eq.'v') uvflg=.true.
      enddo
      kd=tkd

      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))
c
c --- compute scale factor and offset
      fldmin=abs(fillr8)
      fldmax=-abs(fillr8)
c
      ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
c$OMP PARALLEL DO PRIVATE(j,i,ij,ijk) 
c$OMP+  REDUCTION(min:fldmin) REDUCTION(max:fldmax)
      do k=1,kd
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
            if (msk(ij).eq.1.and.fld(ijk).lt.fillr8) then
              fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
            endif
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
        endif
        ofs=0.d0
      else
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=(fldmax-fldmin)/dble(2*i2max)
        endif
        ofs=0.5*(fldmin+fldmax)
      endif
      arng(1)=fldmin
      arng(2)=fldmax
c
c --- Define attributes
      call ncerro(nfmpi_redef(ncid))
      clen=2
      call ncerro(nfmpi_put_att_double(ncid,rhid,'actual_range',
     .                  nf_double,clen,arng))
      clen=1
      call ncerro(nfmpi_put_att_double(ncid,rhid,'scale_factor',
     . nf_double,clen,scf))
      call ncerro(nfmpi_put_att_double(ncid,rhid,'add_offset',nf_double,
     .                               clen,ofs))
      call ncerro(nfmpi_enddef(ncid))
c
      if (mnproc.ne.1) then
        do n=1,maxdm
          icount(n)=0
        enddo
      endif
c
c --- get mask for the full domain
c$OMP PARALLEL DO PRIVATE(i,ij)
      do j=1,jj
        do i=1,ii
          ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          rmsk(i,j)=msk(ij)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare and write output field    
      do k=1,kd
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(i,ij,ijk) 
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            if (msk(ij).eq.1.and.fld(ijk).lt.fillr8) then
              rfld(i,j)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            else
              rfld(i,j)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcaget(rfldt,rfld,1)
        if (mnproc.eq.1) then
          n=0
          do j=1,jtdm
            do i=1,itdm
              if (rmskt(i,j).gt.0.5) then
                n=n+1
                fldout(n)=rfldt(i,j)+i2fill
              endif
            enddo
          enddo
          istart(2)=max(istart(2),k)
          icount(1)=n
        endif
        call ncerro(nfmpi_put_vara_int2_all(ncid,rhid,istart,
     .                                      icount,fldout))
      enddo
c
c --- Put file back to define mode
      write(lp,*) trim(vnm),arng
#endif
      ELSEIF(IO_TYPE .eq. 0) then
      do n=1,maxdm
        start(n)=1
        count(n)=1
      enddo
c --- define variable 
      if (mnproc.eq.1) then
        call ncsevl(dims,strn,strind)
        if (dims(strind(1,strn):strind(2,strn)).eq.'time') then
          start(strn)=rec
        else
          start(strn)=1
        endif
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .      dimid))
          dimids(n)=dimid
          if (n.eq.2.and.dims(strind(1,n):strind(2,n)).ne.'time') 
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
          if (dims(1:1).eq.'u'.or.dims(1:1).eq.'v') uvflg=.true.
        enddo
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      endif
      call xcbcst(uvflg)
      call xcbcst(kd)
c
c --- compute scale factor and offset
      fldmin=abs(fillr8)
      fldmax=-abs(fillr8)
c
      ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
c$OMP PARALLEL DO PRIVATE(j,i,ij,ijk) 
c$OMP+  REDUCTION(min:fldmin) REDUCTION(max:fldmax)
      do k=1,kd
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1) 
            ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
            if (msk(ij).eq.1.and.fld(ijk).lt.fillr8) then
              fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
            endif
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcmin(fldmin)
      call xcmax(fldmax)
      if (uvflg) then 
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
        endif
        ofs=0.d0
      else
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=(fldmax-fldmin)/dble(2*i2max)
        endif
        ofs=0.5*(fldmin+fldmax)
      endif 
      arng(1)=fldmin
      arng(2)=fldmax
c
c --- Define attributes
      if (mnproc.eq.1) then
        call ncerro(nf90_redef(ncid))
        call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
        call ncerro(nf90_put_att(ncid,rhid,'scale_factor',scf))
        call ncerro(nf90_put_att(ncid,rhid,'add_offset',ofs))
        call ncerro(nf90_enddef(ncid))
      endif
c
c --- get mask for the full domain
c$OMP PARALLEL DO PRIVATE(i,ij)
      do j=1,jj
        do i=1,ii
          ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          rmsk(i,j)=msk(ij)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare and write output field    
      do k=1,kd
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(i,ij,ijk) 
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            if (msk(ij).eq.1.and.fld(ijk).lt.fillr8) then
              rfld(i,j)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-i2fill
            else
              rfld(i,j)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcaget(rfldt,rfld,1)
        if (mnproc.eq.1) then
          n=0
          do j=1,jtdm
            do i=1,itdm
              if (rmskt(i,j).gt.0.5) then
                n=n+1
                fldout(n)=rfldt(i,j)+i2fill
              endif 
            enddo 
          enddo 
          start(2)=max(start(2),k) 
          count(1)=n 
          call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
        endif 
      enddo 
      endif
c
      end subroutine nccopa



      subroutine ncwrtc(vnm,dims,fld) 
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes a string array to nc-file. 
c
c --- Arguments:
c       char(*)  vnm       (in) -  variable name
c       char(*)  dims      (in) -  string with dimension names
c       char(*)  fld(*)    (in) -  input field 
c
      implicit none 
c
      character(len=*) :: vnm,dims,fld(*) 
c
      integer :: n 
      integer, parameter :: maxdm=5 
      integer :: dimids(maxdm),strn,strind(2,maxdm)
#ifdef PNETCDF
      integer(kind=MPI_OFFSET_KIND), dimension(maxdm) :: start,count
      integer(kind=MPI_OFFSET_KIND) nsp,slenmaxp
#endif
      integer ns,slenmax
c
      IF(IO_TYPE .EQ. 1) THEN
#ifdef PNETCDF
      call ncsevl(dims,strn,strind)
      if (strn.gt.2) then
        if (mnproc.eq.1) then
          write(lp,*) 'ncwrtc: number of dimensions has to be <=2'
          call flush(lp)
        endif
        call xcstop('(ncerro)')
               stop '(ncerro)'
      endif
      do n=1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                              dimids(n)))
      enddo
      if (strn.eq.1) then
        nsp=1
      else
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(1),slenmaxp))
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(2),nsp))
      endif
c
c --- write data 
      if (mnproc.eq.1) then
        count(2)=1
        count(1)=slenmaxp
        start(1)=1
        start(2)=1
      else
        count(2)=0
        count(1)=0
        start(1)=1
        start(2)=1
      endif
c --- - define variable
        call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))
c
c --- - leave define mode 

      do n=1,nsp
         if (mnproc.eq.1) start(2)=n
         call ncerro(nfmpi_put_vara_text_all(ncid,rhid,
     .                                       start,count,fld(n)//'X'))
      enddo
#endif
      ELSEIF(IO_TYPE .EQ. 0) THEN
      if (mnproc.eq.1) then
c --- - inquire dimensions
        call ncsevl(dims,strn,strind) 
        if (strn.gt.2) then  
          write(lp,*) 'ncwrtc: number of dimensions has to be <=2'  
          call flush(lp)
          call xchalt('(ncerro)')
          stop '(ncerro)'
        endif
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                               dimids(n)))
        enddo 
        if (strn.eq.1) then 
          ns=1
        else 
          call ncerro(nf90_inquire_dimension(ncid,dimids(1),
     .                                       len=slenmax))
          call ncerro(nf90_inquire_dimension(ncid,dimids(2),len=ns))
        endif
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
c --- - write data 
        do n=1,ns
          call ncerro(nf90_put_var(ncid,rhid,fld(n)(1:slenmax)//'X',
     .                             (/1,n/),(/slenmax,1/))) 
        enddo 
      endif
      ENDIF
c
      end subroutine ncwrtc



      subroutine ncwrtr(vnm,dims,fld,msk,mskflg,sfac,offs,prec)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes field to nc-file as real4 or real8. 
c
c --- Arguments:
c       char(*)  vnm       (in) -  variable name
c       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld       (in) -  field with dimension
c                                  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       int(*)   msk       (in) -  field mask with dimensions 
c                                  (idm+2*nbdy)*(jdm+2*nbdy) 
c       integer  mskflg    (in) -  set to 1 if mask is used, 0 else   
c       real     sfac      (in) -  additional scale factor
c       real     offs      (in) -  additional offset
c       integer  prec      (in) -  precision: 4=real4, 8=real8  
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm,dims 
      integer :: i,j,ij,ijk,k,kd,n
      integer, parameter :: maxdm=5,ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer dimid,dimids(maxdm),strn,prec,strind(2,maxdm),msk(*),
     .  mskflg     
      real(kind=4), allocatable, dimension(:,:,:) :: r4fldt,wr4fldt
      real, allocatable, dimension(:,:,:) :: rfld,wrfld
      real, allocatable, dimension(:,:) :: rmsk,rfldt,rmskt
      integer, dimension(maxdm) :: start,count
      real :: fld(*),sfac,offs,fldmin,fldmax
      character c4*4 
c
c --- Initialise fields
      kd=1
      IF(IO_TYPE .eq. 1) then
#ifdef PNETCDF
      do i=1,5
      istart(i)=1
      icount(i)=1
      enddo
      tkd=1
      allocate(rmsk(ii,jj))
      clen=1
      call ncsevl(dims,strn,strind)
      c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))

      if (c4 .eq.'time') then
        istart(strn)=rec
        icount(strn)=1
      endif

      tkd=1
      do n=1,strn
       call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                             dimid))
        dimids(n)=dimid
        if (n.gt.2.and.
     .      ((n.eq.strn.and.c4.ne.'time').or.
     .       (n.eq.strn-1.and.c4.eq.'time')))
     .    call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
      enddo

      kd=tkd

             call ncerro(nfmpi_inq_varid(ncid,vnm,rhid)) 

      istart(1)=1
      istart(2)=j0+1

      if(mproc .eq. mpe_1(nproc) ) then
      icount(1)=itdm
      icount(2)=jj
      if (kd.gt.1) icount(3)=kd
      else
      icount(1)=0
      icount(2)=0
      endif

c
c --- get mask for the full domain
      if (mskflg.eq.1) then
c$OMP PARALLEL DO PRIVATE(i,ij)
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            rmsk(i,j)=msk(ij)
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=1.
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- Prepare output field   
      fldmin=abs(fillr8)
      fldmax=-abs(fillr8)
      allocate(rfld(ii,jj,kd))
      rfld=fillr8
      do k=1,kd
c$OMP PARALLEL DO PRIVATE(i,ij,ijk)
c$OMP+  REDUCTION(min:fldmin) REDUCTION(max:fldmax)
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            rfld(i,j,k)=fld(ijk)
            if (msk(ij).ne.0.and.fld(ijk).ne.fillr8) then
              fldmin=min(rfld(i,j,k),fldmin)
              fldmax=max(rfld(i,j,k),fldmax)
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
      if (prec.ne.4) then
       allocate(wrfld(itdm,jj,kd))
       wrfld=fillr8
        if (sfac.ne.1..or.offs.ne.0.) then
c$OMP PARALLEL DO PRIVATE(j,i)
          do k=1,kd
            do j=1,jj
              do i=1,ii
                if (rmsk(i,j).lt.0.5.or.rfld(i,j,k).eq.fillr8) then
                  rfld(i,j,k)=fillr8
                else
                  rfld(i,j,k)=rfld(i,j,k)*sfac+offs
                endif
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
c$OMP PARALLEL DO PRIVATE(j,i)
          do k=1,kd
            do j=1,jj
              do i=1,ii
                if (rmsk(i,j).lt.0.5.or.rfld(i,j,k).eq.fillr8)
     .            rfld(i,j,k)=fillr8
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        call xcgetrow(wrfld, rfld, kd)
        call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,icount,
     .                                        wrfld))
        deallocate(wrfld,rfld,rmsk)
      else
      allocate(wr4fldt(itdm,jj,kd))
      allocate(r4fldt(ii,jj,kd))
      r4fldt=fillr8
      wr4fldt=fillr8
        if (sfac.ne.1..or.offs.ne.0.) then
c$OMP PARALLEL DO PRIVATE(j,i)
          do k=1,kd
            do j=1,jj
              do i=1,ii
                if (rmsk(i,j).lt.0.5.or.rfld(i,j,k).eq.fillr8) then
                  r4fldt(i,j,k)=fillr8
                else
                  r4fldt(i,j,k)=rfld(i,j,k)*sfac+offs
                endif
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
c$OMP PARALLEL DO PRIVATE(j,i)
          do k=1,kd
            do j=1,jj
              do i=1,ii
                if (rmsk(i,j).lt.0.5.or.rfld(i,j,k).eq.fillr8) then
                  r4fldt(i,j,k)=fillr8
                else
                  r4fldt(i,j,k)=rfld(i,j,k)
                endif
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        call xcgetrow4(wr4fldt, r4fldt, kd)
        call ncerro(nfmpi_put_vara_real_all(ncid,rhid,istart,icount,
     .                                      wr4fldt))
        deallocate(wr4fldt,r4fldt,rfld,rmsk)
      endif
#endif
      ELSEIF(IO_TYPE .eq. 0) then 
      do n=1,maxdm
        start(n)=1
        count(n)=1
      enddo
      count(1)=itdm
      count(2)=jtdm

      allocate(rfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      allocate(rmsk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(rmskt(itdm,jtdm))
      allocate(rfldt(itdm,jtdm))
      if(prec .eq. 4)allocate(r4fldt(itdm,jtdm,1))
      if (mnproc.eq.1) then
        call ncsevl(dims,strn,strind)
        c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4.eq.'time') then
          start(strn)=rec
        else
          start(strn)=1
        endif
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                               dimid))
          dimids(n)=dimid
          if (n.gt.2.and.
     .        ((n.eq.strn.and.c4.ne.'time').or.
     .         (n.eq.strn-1.and.c4.eq.'time')))
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
        enddo
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      endif
      call xcbcst(kd)
c
c --- get mask for the full domain
      if (mskflg.eq.1) then 
c$OMP PARALLEL DO PRIVATE(i,ij)
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            rmsk(i,j)=msk(ij)
          enddo
        enddo
c$OMP END PARALLEL DO
      else 
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=1.
          enddo
        enddo
c$OMP END PARALLEL DO
      endif 
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field   
      fldmin=abs(fillr8)
      fldmax=-abs(fillr8)

      do k=1,kd
        if (k.gt.1) start(3)=k 
c$OMP PARALLEL DO PRIVATE(i,ij,ijk)
c$OMP+  REDUCTION(min:fldmin) REDUCTION(max:fldmax)
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            rfld(i,j,1)=fld(ijk)
            if (msk(ij).ne.0.and.fld(ijk).ne.fillr8) then
              fldmin=min(rfld(i,j,1),fldmin)
              fldmax=max(rfld(i,j,1),fldmax)
            endif
          enddo
        enddo  
c$OMP END PARALLEL DO
        call xcaget(rfldt,rfld,1) 
        if (mnproc.eq.1) then
          if (prec.ne.4) then 
            if (sfac.ne.1..or.offs.ne.0.) then   
c$OMP PARALLEL DO PRIVATE(i)
              do j=1,jtdm
                do i=1,itdm
                  if (rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fillr8) then
                    rfldt(i,j)=fillr8
                  else
                    rfldt(i,j)=rfldt(i,j)*sfac+offs
                  endif
                enddo
              enddo
c$OMP END PARALLEL DO
            else
c$OMP PARALLEL DO PRIVATE(i)
              do j=1,jtdm
                do i=1,itdm
                  if (rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fillr8)
     .              rfldt(i,j)=fillr8
                enddo
              enddo
c$OMP END PARALLEL DO
            endif
            call ncerro(nf90_put_var(ncid,rhid,rfldt,start,count))
          else
            if (sfac.ne.1..or.offs.ne.0.) then   
c$OMP PARALLEL DO PRIVATE(i)
              do j=1,jtdm
                do i=1,itdm                
                  if (rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fillr8) then 
                    r4fldt(i,j,1)=fillr8
                  else 
                    r4fldt(i,j,1)=rfldt(i,j)*sfac+offs
                  endif
                enddo 
              enddo 
c$OMP END PARALLEL DO
            else 
c$OMP PARALLEL DO PRIVATE(i)
              do j=1,jtdm
                do i=1,itdm
                  if (rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fillr8) then
                    r4fldt(i,j,1)=fillr8
                  else
                    r4fldt(i,j,1)=rfldt(i,j)
                  endif 
                enddo
              enddo
c$OMP END PARALLEL DO
            endif
            call ncerro(nf90_put_var(ncid,rhid,r4fldt,start,count))
          endif
        endif 
      enddo
      if (prec.ne.4) then
      deallocate(rmsk,rfld,rmskt,rfldt)
      else
      deallocate(r4fldt,rmsk,rfld,rmskt,rfldt)
      endif
      endif 
      call xcmin(fldmin)
      call xcmax(fldmax)
c
c --- Put file back to define mode

      if (mnproc.eq.1 ) then
        if (sfac.gt.0.) then
          write(lp,*) trim(vnm),(fldmin*sfac+offs),(fldmax*sfac+offs)
        else
          write(lp,*) trim(vnm),(fldmax*sfac+offs),(fldmin*sfac+offs)
        endif
      endif
      end subroutine ncwrtr



      subroutine ncwrti(vnm,dims,fld,msk,mskflg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes field to nc-file as int4. 
c
c --- Arguments:
c       char(*)    vnm       (in) -  variable name
c       char(*)    dims      (in) -  axes string, e.g. 'pcomp z time'
c       integer(*) fld       (in) -  field with dimension
c                                   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       integer(*) msk       (in) -  field mask with dimensions 
c                                   (idm+2*nbdy)*(jdm+2*nbdy) 
c       integer    mskflg    (in) -  set to 1 if mask is used, 0 else   
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm,dims 
      integer :: i,j,k,n,kd
      integer, parameter :: maxdm=5,ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer :: dimid,dimids(maxdm),strn,strind(2,maxdm),msk(*),mskflg,
     .  fld(*)     
      real, allocatable, dimension(:,:) :: rmsk,rfldt,rmskt
      real, allocatable, dimension(:,:,:) :: rfld,wrfld
      integer, allocatable, dimension(:,:,:) :: irfld
      integer, dimension(maxdm) :: start,count
      integer  :: fill
      character c4*4
c
c --- Initialise fields
      kd=1
c  
c --- define variable 
      IF(IO_TYPE .EQ. 1) THEN
#ifdef PNETCDF
      do i=1,5
      istart(i)=1
      icount(i)=1
      enddo

      tkd=1      
      fill=nf_fill_int
      call ncsevl(dims,strn,strind)
      c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
      if (c4.eq.'time') then
        istart(strn)=rec
        icount(strn)=1
      endif
      tkd=1
      do n=1,strn
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                              dimid))
        dimids(n)=dimid
        if (n.gt.2.and.
     .      ((n.eq.strn.and.c4.ne.'time').or.
     .       (n.eq.strn-1.and.c4.eq.'time')))
     .    call ncerro(nfmpi_inq_dimlen(ncid,dimid,tkd))
      enddo
      kd=tkd
           
      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))

      istart(1)=1
      istart(2)=j0+1
      if(mproc .eq. mpe_1(nproc) ) then
      icount(1)=itdm
      icount(2)=jj
      if (kd.gt.1) icount(3)=kd
      else
      icount(1)=0
      icount(2)=0
      endif

      allocate(rmsk(ii,jj))
      allocate(rfld(ii,jj,kd))
c
c --- get mask for the full domain
      if (mskflg.eq.1) then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=1.
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c$OMP PARALLEL DO PRIVATE(k,i)
      do j=1,jj
        do k=1,kd
          do i=1,ii
            rfld(i,j,k)=fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,i)
      do j=1,jj
        do k=1,kd
          do i=1,ii
            if (rmsk(i,j).lt.0.5) rfld(i,j,k)=fill
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c 
      allocate(wrfld(itdm,jj,kd))
      allocate(irfld(itdm,jj,kd))
      wrfld=fill
      call xcgetrow(wrfld, rfld, kd)
c$OMP PARALLEL DO PRIVATE(k,i)
      do j=1,jj
        do k=1,kd
          do i=1,itdm
          irfld(i,j,k)=wrfld(i,j,k)
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO

      call ncerro(nfmpi_put_vara_int_all(ncid,rhid,istart,icount,irfld))
c
      deallocate(wrfld,rfld,rmsk,irfld)
#endif
c
      ELSEIF(IO_TYPE .EQ. 0) THEN
      do n=1,maxdm
      start(n)=1
      count(n)=1
      enddo
      count(1)=itdm
      count(2)=jtdm
      allocate(rfld(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      allocate(rmsk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(rfldt(itdm,jtdm))
      allocate(rmskt(itdm,jtdm))
      fill=nf90_fill_int
      if (mnproc.eq.1) then
        call ncsevl(dims,strn,strind)
        c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4.eq.'time') then
          start(strn)=rec
        else
          start(strn)=1
        endif
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                               dimid))
          dimids(n)=dimid
          if (n.gt.2.and.
     .        ((n.eq.strn.and.c4.ne.'time').or.
     .         (n.eq.strn-1.and.c4.eq.'time')))
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
        enddo
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
      endif
      call xcbcst(kd)
c
c --- get mask for the full domain
      if (mskflg.eq.1) then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=1.
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
        if (k.gt.1) start(3)=k
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            rfld(i,j,1)=fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcaget(rfldt,rfld,1)
        if (mnproc.eq.1) then
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,jtdm
            do i=1,itdm
              if (rmskt(i,j).lt.0.5) rfldt(i,j)=fill
            enddo
          enddo
c$OMP END PARALLEL DO
          call ncerro(nf90_put_var(ncid,rhid,rfldt,start,count))
        endif
      enddo
      ENDIF
c
      end subroutine ncwrti


      
      subroutine ncwrt1(vnm,dims,fld)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real8 field from master process to nc-file. 
c
c --- Arguments:
c       char(*)  vnm       (in) -  variable name
c       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld       (in) -  input field  
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm,dims
c
      integer, parameter :: maxdm=5
      integer :: n,strind(2,maxdm),ndims
#ifdef PNETCDF
      integer (kind=MPI_OFFSET_KIND), dimension(maxdm) :: 
     . dimlenp
#endif
      integer , dimension(maxdm) :: dimids
      integer, dimension(maxdm) :: start,count,dimlen
c
      real :: fld(*)
c
      IF(IO_TYPE .EQ. 1) THEN
#ifdef PNETCDF
      call ncsevl(dims,ndims,strind)
      do n=1,ndims
        call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                              dimids(n)))
        call ncerro(nfmpi_inq_dimlen(ncid,dimids(n),dimlenp(n)))
      enddo
      call ncerro(nfmpi_inq_varid(ncid,vnm,rhid))
      do n=1,maxdm
        istart(n)=1
        icount(n)=0
      enddo
      if(mnproc .eq. 1)then
      do n=1,ndims-1
        icount(n)=dimlenp(n)
      enddo
      
      if (dims(strind(1,ndims):strind(2,ndims)).eq.'time') then
        istart(ndims)=rec
        icount(ndims)=1
      else
        istart(ndims)=1
        icount(ndims)=dimlenp(ndims)
      endif
      endif
c
c --- Write data to file 
      call ncerro(nfmpi_put_vara_double_all(ncid,rhid,istart,
     .                                      icount,fld))
c
c --- Put file back to define mode
#endif
      ELSEIF(IO_TYPE .EQ. 0) THEN
      if (mnproc.eq.1) then
c
c --- - Define variable 
        call ncsevl(dims,ndims,strind)
        do n=1,ndims
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .                               dimids(n)))
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),
     .                                       len=dimlen(n)))
        enddo
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
c
c --- - Prepare start and count
        do n=1,maxdm
          start(n)=1
          count(n)=1
        enddo
        do n=1,ndims-1 
          count(n)=dimlen(n) 
        enddo
        if (dims(strind(1,ndims):strind(2,ndims)).eq.'time') then
          start(ndims)=rec
          count(ndims)=1 
        else
          start(ndims)=1
          count(ndims)=dimlen(ndims) 
        endif
c
c --- - Write data to file 
        call ncerro(nf90_put_var(ncid,rhid,fld(1:product(count)),start,
     .                           count))
      endif
      ENDIF
c
      end subroutine ncwrt1



c ----------------------------------------------------------------------
c --- auxilary routines ------------------------------------------------
c ----------------------------------------------------------------------

 

      subroutine ncerro(ncstatus)
c
c ----------------------------------------------------------------------
c --- Description:
c       Display netcdf error massages  
c
c --- Arguments:
c       int ncstatus (in) -  netcdf status 
c ----------------------------------------------------------------------
c
      implicit none
c
      integer ncstatus
c
      IF(IO_TYPE .EQ. 1) THEN
#ifdef PNETCDF
      if (ncstatus.ne.nf_noerr) then
        write(lp,*) 'NetCDF error:',nfmpi_strerror(ncstatus)
        call flush(lp)
        call xchalt('(ncerro)')
               stop '(ncerro)'
      endif
#endif
      ELSEIF(IO_TYPE .EQ. 0) THEN
      if (ncstatus.ne.nf90_noerr) then 
        write(lp,*) 'NetCDF error: ',nf90_strerror(ncstatus)
        call flush(lp)
        call xchalt('(ncerro)')
               stop '(ncerro)'
      endif 
      ENDIF
c
      end subroutine ncerro



      subroutine ncsevl(strg,strgn,strgind)
c
c ----------------------------------------------------------------------
c --- Description:
c       Finds the number and the locations of sub-strings. 
c       Valid deliminators are ' ', '-' and ':'.  
c
c --- Arguments:
c       char(*) strg   (in)  -  input string
c       int strgn      (out) -  number of sub-strings
c       int strgind(*) (out) -  start/end locations of sub-strings (the
c                              dimension must at least equal to strgn*2)
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) strg
      character charold,charnew
      integer strgn,strgind(*),i
c
      charold=' '
      strgn=0
      do i=1,len(strg)
         charnew=strg(i:i)
         if ((charold.eq.' '.or.charold.eq.'-'.or.charold.eq.':').and.
     .       (charnew.ne.' '.and.charnew.ne.'-'.and.charnew.ne.':'))
     .   then
            strgn=strgn+1
            strgind(strgn)=i
         elseif ((charnew.eq.' '.or.charnew.eq.'-'.or.charnew.eq.':')
     .     .and.(charold.ne.' '.and.charold.ne.'-'.and.charold.ne.':'))
     .   then
            strgn=strgn+1
            strgind(strgn)=i-1
         endif
         charold=charnew
      enddo
      if (mod(strgn,2).eq.1) then
         strgn=strgn+1
         strgind(strgn)=len(strg)
      endif
      strgn=int(strgn/2.)
c
      end subroutine ncsevl


      subroutine nccaln(days,calendar,y1,m1,d1,y2,m2,d2)
c
c ----------------------------------------------------------------------
c --- Description:
c       nccaln calculates the number of days between two
c       dates. Valid calendars are '360_day', 'noleap' = '365_day',
c       'all_leap' = '366_day', 'julian', 'proleptic_gregorian' and
c       'standard' = 'mixed' = 'gregorian'.
c
c --- Arguments:
c       int days         (out) -  number of days
c       char(*) calendar (in)  -  choice of calendar
c       int y1           (in)  -  start year
c       int m1           (in)  -  start month
c       int d1           (in)  -  start day
c       int y2           (in)  -  last year
c       int m2           (in)  -  last month
c       int d2           (in)  -  last day
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) calendar
      integer n1,y1,m1,d1,n2,y2,m2,d2,days,nd1,nd2,nd3
c
      if (calendar(1:3).eq.'sta'.or.calendar(1:3).eq.'mix'.or.
     .    calendar(1:3).eq.'gre') then
        if (y1.lt.1582.or.(y1.eq.1582.and.m1.lt.10).or.
     .      (y1.eq.1582.and.m1.eq.10.and.d1.le.15)) then
          call ncdnum('julian',y1,m1,d1,n1)
        else
          call ncdnum('julian',1582,10,15,nd1)
          call ncdnum('proleptic',1582,10,15,nd2)
          call ncdnum('proleptic',y1,m1,d1,nd3)
          n1=nd1-nd2+nd3-10
        endif
        if (y2.lt.1582.or.(y2.eq.1582.and.m2.lt.10).or.
     .      (y2.eq.1582.and.m2.eq.10.and.d2.le.15)) then
          call ncdnum('julian',y2,m2,d2,n2)
        else
          call ncdnum('julian',1582,10,15,nd1)
          call ncdnum('proleptic',1582,10,15,nd2)
          call ncdnum('proleptic',y2,m2,d2,nd3)
          n2=nd1-nd2+nd3-10
        endif
      else
        call ncdnum(calendar,y1,m1,d1,n1)
        call ncdnum(calendar,y2,m2,d2,n2)
      endif
      days=n2-n1
c
      end subroutine nccaln



      subroutine ncdnum(calendar,year,month,day,ndays)
c
c ----------------------------------------------------------------------
c --- Description:
c       Auxilary subroutine for nccaln which calculates the number
c       of days between a given date and a reference date. Thereby the
c       reference date varies from calendar to calendar and is not
c       standerized. The purpose is more of differential nature rather
c       than to create an absolute time axis.
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) calendar
      integer y,year,month,day,ndays,yr1,yr4,yr100,yr400
      integer acc_leap(12),acc_noleap(12),yoffset
      parameter(yoffset=1000000,yr1=365,yr4=1461,yr100=36524,
     .   yr400=146097)
      data acc_leap   /0,31,60,91,121,152,182,213,244,274,305,335/
      data acc_noleap /0,31,59,90,120,151,181,212,243,273,304,334/
c
      y=year+yoffset
      if (calendar(1:3).eq.'pro') then
        ndays=int((y-1)/400.)*yr400
     .       +int(mod((y-1),400)/100.)*yr100
     .       +int(mod(mod((y-1),400),100)/4.)*yr4
     .       +mod(mod(mod((y-1),400),100),4)*yr1
        if (month.gt.2.and.((mod(y,4).eq.0.and.mod(y,100).ne.0).or.
     .     mod(y,400).eq.0)) then
          ndays=ndays+acc_leap(month)+day
        else
          ndays=ndays+acc_noleap(month)+day
        endif
      elseif (calendar(1:3).eq.'jul') then
        ndays=int((y-1)/4.)*yr4+mod(y,4)*yr1
        if (month.gt.2.and.mod(y,4).eq.0) then
          ndays=ndays+acc_leap(month)+day
        else
          ndays=ndays+acc_noleap(month)+day
        endif
      elseif (calendar(1:3).eq.'360') then
        ndays=(y-1)*360+30*(month-1)+day
      elseif (calendar(1:3).eq.'nol'.or.calendar(1:3).eq.'365') then
        ndays=(y-1)*365+acc_noleap(month)+day
      elseif (calendar(1:3).eq.'all'.or.calendar(1:3).eq.'366') then
        ndays=(y-1)*366+acc_leap(month)+day
      endif
c
      end subroutine ncdnum
   

      logical function ncinqa(nm)
c
c ----------------------------------------------------------------------
c --- Description:
c       Inquires if attribute exists.  
c
c --- Arguments:
c       char(*)       nm    (in)  -  attribute name
c       logical   ncinqa   (out)  -  return value
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) nm
#ifdef PNETCDF
      integer(kind=MPI_OFFSET_KIND) attlen
#endif
      integer attlen1   
c     
      IF(IO_TYPE .EQ. 1) THEN
#ifdef PNETCDF
      status=nfmpi_inq_attlen(ncid,nf_global,trim(nm),attlen)
      if (status.ne.nf_noerr) then
      ncinqa=.false.
      else
         if (attlen.ne.0) then
           ncinqa=.true.
         else
           ncinqa=.false.
         endif
      endif
#endif
      ELSEIF(IO_TYPE .EQ. 0) THEN
      IF(mnproc .eq. 1) THEN
      status=nf90_inquire_attribute(ncid,nf90_global,trim(nm),
     .  len=attlen1)
      if (status.ne.nf90_noerr) then

         ncinqa=.false.
      else
         if (attlen1.ne.0) then
           ncinqa=.true.
         else 
           ncinqa=.false.
         endif 
      endif
      endif
      call xcbcst(ncinqa)
      ENDIF
c
      end function ncinqa



      logical function ncinqv(vnm)
c
c ----------------------------------------------------------------------
c --- Description:
c       Inquires if variable exists.  
c
c --- Arguments:
c       char(*)      vnm    (in)  -  variable name
c       logical   ncinqv   (out)  -  return value
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) vnm
c     
      IF(IO_TYPE .EQ. 1) THEN
#ifdef PNETCDF
      status=nfmpi_inq_varid(ncid,trim(vnm),rhid)
      if (status.ne.nf_noerr) then
      ncinqv=.false.
      else
        ncinqv=.true.
      endif
#endif
      ELSEIF(IO_TYPE .EQ. 0) THEN
      IF(mnproc .eq. 1) THEN
      status=nf90_inq_varid(ncid,trim(vnm),rhid)
      if (status.ne.nf90_noerr) then 
        ncinqv=.false.
      else
        ncinqv=.true. 
      endif
      endif
      call xcbcst(ncinqv)
      ENDIF 
c
      end function ncinqv
      
      subroutine ncdefvar(vnm,dims,itype,iatt)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real8 field from master process to nc-file. 
c
c --- Arguments:
c       char(*)  vnm       (in) -  variable name
c       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld       (in) -  input field  
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) :: vnm,dims
c
      integer, parameter :: maxdm=5
      integer :: n,strind(2,maxdm),ndims,itype,iatt
      integer , dimension(maxdm) :: dimids
c
      if (rec.eq.1) then
      call ncsevl(dims,ndims,strind)
      if( IO_TYPE .EQ. 1) then
#ifdef PNETCDF
      do n=1,ndims
         call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .      dimids(n)))
      enddo

        call ncerro(nfmpi_def_var(ncid,vnm,itype,ndims,
     .        dimids,rhid))
        clen=1
        if(iatt .eq. 8) then
        call ncerro(nfmpi_put_att_double(ncid,rhid,'_FillValue',
     .                                  nf_double,clen,fillr8))
        elseif(iatt .eq. 4) then
        call ncerro(nfmpi_put_att_real(ncid,rhid,'_FillValue',
     .                                  nf_real,clen,fillr4))
        elseif(iatt .eq. 2) then
        call ncerro(nfmpi_put_att_int(ncid,rhid,'_FillValue',
     .                                  nf_int,clen,nf_fill_int))
        endif
#endif
        elseif( IO_TYPE .EQ. 0) then
        do n=1,ndims
          call ncerro(nf90_inq_dimid(ncid, dims(strind(1,n):strind(2,n))
     .  ,dimids(n)))
        enddo
        call ncerro(nf90_def_var(ncid,vnm,itype,dimids(1:ndims),rhid))
        if(iatt .eq. 8) then
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr8))
        elseif(iatt .eq. 4) then
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr4))
        elseif(iatt .eq. 2) then
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',nf90_fill_int))
        endif
        endif
        endif

      end subroutine ncdefvar

      subroutine ncdefvar3d(frmt,cmpflg,gridid,
     .  vnm,vlngnm,vstdnm,vunits,isize)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real8 field from master process to nc-file. 
c
c --- Arguments:
c       char(*)  vnm       (in) -  variable name
c       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld       (in) -  input field  
c ----------------------------------------------------------------------
c
      implicit none
c 
      integer, parameter :: maxdm=5
      integer :: frmt,cmpflg,n,isize
      character(len=*) :: gridid,vnm,vlngnm,vstdnm,vunits
c
      character(len=100) :: dims
      integer :: strind(2,maxdm),ndims
      integer , dimension(maxdm) :: dimids
c
      real :: arng(2)
      integer*2 i2min,vrng(2)
      if (rec.eq.1) then
      i2min=-i2max
      vrng(1)=i2min
      vrng(2)=i2max
      arng(1)=fillr8
      arng(2)=fillr8
c
c
c --- Check whether field should be written
      if (frmt.eq.0) return
c
c --- Create dimension string 
      if(isize .eq. 0) then
      if (cmpflg.eq.1) then
        dims=gridid(1:1)//'comp time'
      else
        dims='x y time'
      endif
      elseif (isize .eq. 1) then
      if (cmpflg.eq.1) then
        dims=gridid(1:1)//'comp sigma time'
      else
        dims='x y sigma time'
      endif
      elseif (isize .eq. 2) then
      if (cmpflg.eq.1) then
        dims=gridid(1:1)//'comp depth time'
      else
        dims='x y depth time'
      endif
      elseif (isize .eq. 3) then
      if (cmpflg.eq.1) then
        dims='pcomp ks time'
      else
        dims='x y ks time'
      endif
      elseif (isize .eq. 4) then
      if (cmpflg.eq.1) then
        dims='pcomp time'
      else
        dims='x y time'
      endif
      else
      write (lp,*) 'unknown data size!'
        call xchalt('(ncdefvar3d)')
               stop '(ncdefvar3d)'
      endif

      call ncsevl(dims,ndims,strind)
      if( IO_TYPE .EQ. 1) then
#ifdef PNETCDF
      do n=1,ndims
         call ncerro(nfmpi_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .      dimids(n)))
      enddo

c
c --- Check output format
      if (frmt.eq.2) then
        call ncerro(nfmpi_def_var(ncid,vnm,nf_int2,ndims,
     .        dimids,rhid))
        clen=1
          call ncerro(nfmpi_put_att_int2(ncid,rhid,'_FillValue',nf_int2,
     .              clen,i2fill))
        clen=2
        call ncerro(nfmpi_put_att_int2(ncid,rhid,'valid_range',nf_int2,
     .              clen,vrng))

        call ncerro(nfmpi_put_att_double(ncid,rhid,'actual_range',
     .                nf_double,clen,arng))
        clen=1
        call ncerro(nfmpi_put_att_double(ncid,rhid,'scale_factor',
     .                                   nf_double,clen,fillr8))
        call ncerro(nfmpi_put_att_double(ncid,rhid,'add_offset',
     .                nf_double,clen,fillr8))

        if (cmpflg.eq.1) then
        clen=3
        call ncerro(nfmpi_put_att_text(ncid,rhid,'compress',clen,'x y'))
        endif
      elseif (frmt.eq.4) then
        if (cmpflg.eq.1) then
          call ncerro(nfmpi_def_var(ncid,vnm,nf_real,ndims,
     .        dimids,rhid))
          clen=3
          call ncerro(nfmpi_put_att_text(ncid,rhid,'compress',clen,
     .                                   'x y'))
        else
          call ncerro(nfmpi_def_var(ncid,vnm,nf_real,ndims,
     .        dimids,rhid))
        endif
         clen=1
          call ncerro(nfmpi_put_att_real(ncid,rhid,'_FillValue',
     .          nf_real,clen,fillr4))
      elseif (frmt.eq.8) then
        if (cmpflg.eq.1) then
          call ncerro(nfmpi_def_var(ncid,vnm,nf_double,ndims,
     .        dimids,rhid))
          clen = 3
          call ncerro(nfmpi_put_att_text(ncid,rhid,'compress',clen,
     .                                   'x y'))
        else
        call ncerro(nfmpi_def_var(ncid,vnm,nf_double,ndims,
     .        dimids,rhid))
        endif
        clen=1
        call ncerro(nfmpi_put_att_double(ncid,rhid,'_FillValue',
     .          nf_double,clen,fillr8))
      else
        write (lp,*) 'unknown output format!'
        call xchalt('(ncdefvar3d)')
               stop '(ncdefvar3d)'
      endif
#endif
      elseif( IO_TYPE .EQ. 0) then
      do n=1,ndims
          call ncerro(nf90_inq_dimid(ncid, dims(strind(1,n):strind(2,n))
     .  ,dimids(n)))
      enddo

        if (frmt.eq.2) then
        call ncerro(nf90_def_var(ncid,vnm,nf90_int2,dimids(1:ndims),
     .                           rhid))
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',i2fill))
        call ncerro(nf90_put_att(ncid,rhid,'valid_range',vrng))
        call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
        call ncerro(nf90_put_att(ncid,rhid,'scale_factor',fillr8))
        call ncerro(nf90_put_att(ncid,rhid,'add_offset',fillr8))
        if (cmpflg.eq.1) 
     .  call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))

        elseif (frmt.eq.4) then
        call ncerro(nf90_def_var(ncid,vnm,nf90_real,
     .        dimids(1:ndims),rhid))
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr4))
        if (cmpflg.eq.1)
     .  call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))

        elseif (frmt.eq.8) then
        call ncerro(nf90_def_var(ncid,vnm,nf90_double,
     .        dimids(1:ndims),rhid))
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr8))
        if (cmpflg.eq.1)
     .  call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))

      else
        write (lp,*) 'unknown output format!'
        call xchalt('(ncdefvar3d)')
               stop '(ncdefvar3d)'
      endif

      endif

c --- Define attributes
      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
      call ncattr('coordinates',
     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
      endif
      end subroutine ncdefvar3d

      subroutine ncedef
c
c ----------------------------------------------------------------------
c --- Description:
c       Closes NetCDF file
c ----------------------------------------------------------------------
c
      implicit none
c
      IF(IO_TYPE .eq. 1) then 
#ifdef PNETCDF
      call ncerro(nfmpi_enddef(ncid))
#endif
      ELSEIF(IO_TYPE .eq. 0) then
      if (mnproc.eq.1) then
      if (flgpad) then
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,81920,4,40960,4))
          flgpad=.true.
        endif

      endif
      ENDIF
      end subroutine ncedef 
      end module mod_nctools
