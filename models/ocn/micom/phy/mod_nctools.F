      module mod_nctools
c ----------------------------------------------------------------------
c --- NetCDF tools -----------------------------------------------------
c ----------------------------------------------------------------------
c
c Description:
c    nctools.f is a boundle of fortran subroutines designed to write
c    output from the mpi version of MICOM to NetCDF data-format. 
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
c
      implicit none 
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
      integer, private, save :: ncid,rhid,status,rec
c
      contains             


      subroutine ncfopn(fnm,faccess,frmt,irec)
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
     .         nf90__64bit_offset),ncid))
            elseif (frmt(1:1).eq.'h'.or.frmt(1:1).eq.'H') then
              call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__hdf5),
     .          ncid))
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
     .        nf90__64bit_offset),ncid))
          elseif (frmt(1:1).eq.'h'.or.frmt(1:1).eq.'H') then
            call ncerro(nf90_create(fnm,OR(nf90_clobber,nf90__hdf5),
     .        ncid))
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
      if (mnproc.eq.1) then 
        call ncerro(nf90_close(ncid))
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
      integer :: dim,dimid
      character*(*) :: dnm 
c
c --- define dimension
      if (mnproc.eq.1.and.rec.eq.1)  
     .  call ncerro(nf90_def_dim(ncid,dnm,dim,dimid))
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
      do j=1,jj
        do i=1,ii
          rmsk(i,j)=msk(i,j)
        enddo
      enddo 
      call xcaget(rmskt,rmsk,1)
      if (mnproc.eq.1.and.rec.eq.1) then 
        n=0
        if(flg.eq.1)then 
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
              if(rmskt(i,j).gt.0.5) n=n+1
            enddo 
          enddo  
        endif  
        call ncerro(nf90_def_dim(ncid,dnm,n,dimid))
        if(flg.eq.1)then
          call ncerro(nf90_def_var(ncid,dnm,nf90_int,dimid,rhid))
          if(flgpad)then 
            call ncerro(nf90_enddef(ncid))
          else
            call ncerro(nf90_enddef(ncid,40960,4,0,4))
            flgpad=.true.
          endif 
          call ncerro(nf90_put_var(ncid,rhid,ivec(1:n)))
          call ncerro(nf90_redef(ncid))
        endif 
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
      if (mnproc.eq.1.and.rec.eq.1)  
     .  call ncerro(nf90_put_att(ncid,rhid,attname,attrib))
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
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,rval))
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
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,rval))
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
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,ival))
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
c
      if (mnproc.eq.1.and.rec.eq.1) 
     .  call ncerro(nf90_put_att(ncid,nf90_global,vnm,ival))
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
      if (mnproc.eq.1) 
     .  call ncerro(nf90_get_att(ncid,nf90_global,vnm,rval))
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
      if (mnproc.eq.1) 
     .  call ncerro(nf90_get_att(ncid,nf90_global,vnm,rval))
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
      if(mnproc.eq.1)
     .  call ncerro(nf90_get_att(ncid,nf90_global,vnm,ival))
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
      if(mnproc.eq.1)
     .  call ncerro(nf90_get_att(ncid,nf90_global,vnm,ival))
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
      integer :: dimid,bstrn,sstrn,strind(2,10),days,by,bm,bd,bh,
     .   sy,sm,sd,sh
      character*(*) :: calendar,units,startdate
      real :: datenum,caloffset
      data by/1/,bm/1/,bd/1/,bh/0/,sy/1/,sm/1/,sd/1/,sh/0/
c
      if (mnproc.eq.1) then 
c --- define time dimension
      if (rec.eq.1) then 
        call ncerro(nf90_def_dim(ncid,'time',nf90_unlimited,dimid))
      else 
        call ncerro(nf90_inq_dimid(ncid,'time',dimid))
      endif 
c
c --- define time coordinate variable
      if (rec.eq.1) then 
        call ncerro(nf90_def_var(ncid,'time',nf90_double,dimid,rhid))
        call ncerro(nf90_put_att(ncid,rhid,'long_name','time'))
        call ncerro(nf90_put_att(ncid,rhid,'units',units))
        call ncerro(nf90_put_att(ncid,rhid,'calendar',calendar))
      else 
        call ncerro(nf90_inq_varid(ncid,'time',rhid))
      endif 
c
c --- analyse units string
      call ncsevl(units,bstrn,strind)
      read(units(strind(1,3):strind(2,3)),*) by
      if (bstrn.ge.4) read(units(strind(1,4):strind(2,4)),*) bm
      if (bstrn.ge.5) read(units(strind(1,5):strind(2,5)),*) bd
      if (bstrn.ge.6) read(units(strind(1,6):strind(2,6)),*) bh
c
c --- analyse startdate string
      call ncsevl(startdate,sstrn,strind)
      read(startdate(strind(1,1):strind(2,1)),*) sy
      if (sstrn.ge.2) read(startdate(strind(1,2):strind(2,2)),*) sm
      if (sstrn.ge.3) read(startdate(strind(1,3):strind(2,3)),*) sd
      if (sstrn.ge.4) read(startdate(strind(1,4):strind(2,4)),*) sh
c
c --- calculate calendar offset
      if (units(1:1).eq.'y') then
        caloffset=real(sy-by)
      elseif(units(1:1).eq.'m')then
        caloffset=real(sy-by)*12.+real(sm-bm)
      elseif(units(1:1).eq.'d')then
        call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
        caloffset=real(days)
      elseif(units(1:1).eq.'h')then
        call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
        caloffset=real(days)*24.+real(sh-bh)
      endif
      datenum=datenum+caloffset
c
c --- Change to data mode and write time
      if(flgpad)then 
        call ncerro(nf90_enddef(ncid))
      else
        call ncerro(nf90_enddef(ncid,40960,4,0,4))
        flgpad=.true.
      endif 
      call ncerro(nf90_put_var(ncid,rhid,datenum,(/rec/)))
c
c --- Switch back to define mode
      call ncerro(nf90_redef(ncid))
c
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
      integer ndims,dimid,dimids(maxdm),dimlen,strn,strind(2,maxdm),
     .  msk(*),mskflg
      integer, dimension(maxdm) :: start,count
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real, dimension(itdm*jtdm) :: rfldtcmp
      real :: rkd,fill,fld(*),ofs,scf,rcmpflg
c
c --- Initialise fields
      rcmpflg=0.
      rkd=1.
      do n=1,maxdm
        start(n)=1
        count(n)=1
      enddo  
c       
c --- Inquire dimensions, ect.
      if(mnproc.eq.1)then      
        status=nf90_inq_varid(ncid,vnm,rhid)
        if(status.ne.nf90_noerr) then 
          write(lp,*) 'WARNINIG: Problems reading variable ',trim(vnm)
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
          if(dimname(2:5).eq.'comp') rcmpflg=1. 
          if (nint(n+rcmpflg).eq.3) rkd=dimlen
          ndm=ndm*dimlen
        enddo
      endif 
      call xcmaxr(rcmpflg)
      call xcmaxr(rkd)
      kd=nint(rkd)
c
c --- Get attributes
      scf=-abs(nf90_fill_double)
      ofs=-abs(nf90_fill_double)
      if(mnproc.eq.1)then      
        status=nf90_get_att(ncid,rhid,'scale_factor',scf)
        if(status.ne.nf90_noerr) scf=1.
        status=nf90_get_att(ncid,rhid,'add_offset',ofs)
        if(status.ne.nf90_noerr) ofs=0.
      endif 
      call xcmaxr(scf)
      call xcmaxr(ofs)
c
c --- get mask for the full domain  
      if(mskflg.eq.1)then
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          enddo
        enddo
        call xcaget(rmskt,rmsk,1)
      endif
c
c --- Get data
      do k=1,kd
        if(mnproc.eq.1)then 
          if(rcmpflg.gt.0.5)then 
            start(2)=k
            count(1)=nint(dble(ndm)/dble(kd)) 
          else 
            start(3)=k 
            count(1)=itdm
            count(2)=jtdm
          endif
          if(rcmpflg.gt.0.5)then 
            call ncerro(nf90_get_var(ncid,rhid,rfldtcmp,start,count))
            n=0 
            do j=1,jtdm
              do i=1,itdm
                if(rmskt(i,j).gt.0.5) then
                  n=n+1
                  rfldt(i,j)=rfldtcmp(n) 
                endif 
              enddo
            enddo  
          else
            call ncerro(nf90_get_var(ncid,rhid,rfldt,start,count))
          endif  
        endif
        call xcaput(rfldt,rfld,1)
        if(mskflg.eq.1)then 
          if(scf.eq.1..and.ofs.eq.0)then 
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                if(msk(ij).eq.1)then
                  fld(ij+(k-1)*ijdm)=rfld(i,j)
                else
                  fld(ij+(k-1)*ijdm)=fill
                endif
              enddo
            enddo
          else
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                if(msk(ij).eq.1)then
                  fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs
                else
                  fld(ij+(k-1)*ijdm)=fill
                endif
              enddo
            enddo
          endif  
        else
          if(scf.eq.1..and.ofs.eq.0)then
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs 
              enddo
            enddo     
          else
            do j=1,jj
              do i=1,ii
                ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                fld(ij+(k-1)*ijdm)=rfld(i,j)
              enddo
            enddo
          endif
        endif        
      enddo 
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
      integer :: i,j,ij,ijk,k,kd,n
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),strn,strind(2,maxdm),msk(*),mskflg    
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld
      real, dimension(itdm,jtdm) :: rfldt
      real scf,sfac,ofs,offs,arng(2),fldmin,fldmax,fld(*),rkd,ruvflg
     .  ,rfill
      integer*2, dimension(itdm,jtdm) :: fldout
      integer*2 fill,i2max,i2min,vrng(2)
      parameter(rfill=nf90_fill_double,fill=-32768,i2max=32767, 
     .  i2min=-i2max)
c
c --- Initialise fields
      ruvflg=0.
      rkd=1.
      kd=1
      vrng(1)=i2min
      vrng(2)=i2max
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
     .        dims(strind(1,n):strind(2,n)),dimid))
            dimids(n)=dimid
            if ((n.eq.strn.and.c4.ne.'time').or.
     .        (n.gt.2.and.n.eq.strn-1.and.c4.eq.'time')) 
     .        call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
            if (dims(1:1).eq.'u'.or.dims(1:1).eq.'v') ruvflg=1.
         enddo
         rkd=kd
         if (rec.eq.1) then 
           call ncerro(nf90_def_var(ncid,vnm,nf90_int2,dimids(1:strn),
     .       rhid))
         else 
           call ncerro(nf90_inq_varid(ncid,vnm,rhid))
         endif
      endif
      call xcmaxr(ruvflg)
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- compute scale factor and offset
      fldmin=abs(nf90_fill_double)
      fldmax=-abs(nf90_fill_double)
      if (mskflg.eq.1) then 
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  PRIVATE(ij,ijk) 
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
c$OMP+  PRIVATE(ij,ijk) 
        do k=1,kd
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk=ij+(k-1)*ijdm
              if (msk(ij).eq.1.and.fld(ijk).ne.rfill) then
                fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
                fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+   PRIVATE(ijk) 
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
      call xcminr(fldmin)
      call xcmaxr(fldmax)
      if (ruvflg.gt.0.5) then 
        if (fldmin.ge.fldmax )then
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
        if (rec.eq.1) then 
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fill))
          call ncerro(nf90_put_att(ncid,rhid,'valid_range',vrng))
          call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
          call ncerro(nf90_put_att(ncid,rhid,'scale_factor',scf))
          call ncerro(nf90_put_att(ncid,rhid,'add_offset',ofs))
        endif
        if(flgpad)then 
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif 
      endif
c
c --- Prepare and write output field
      do k=1,kd
        if (k.gt.1) start(3)=k
        if (mskflg.eq.1) then
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(ij) 
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              if (msk(ij).eq.1) then
                rfld(i,j)=nint(((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf)
     .            -fill
              else
                rfld(i,j)=0 
              endif
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (mskflg.eq.2) then
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(ij,ijk) 
          do j=1,jj
            do i=1,ii
              ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
              ijk=ij+(k-1)*ijdm
              if (msk(ij).eq.1.and.fld(ijk).ne.rfill) then
                rfld(i,j)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-fill
              else
                rfld(i,j)=0
              endif
            enddo
          enddo
c$OMP END PARALLEL DO
        else
c$OMP PARALLEL DO 
c$OMP+  PRIVATE(ijk) 
          do j=1,jj
            do i=1,ii 
              ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
              rfld(i,j)=nint(((fld(ijk)*sfac)+offs-ofs)/scf)-fill
            enddo 
          enddo   
c$OMP END PARALLEL DO
        endif  
        call xcaget(rfldt,rfld,1)
        if(mnproc.eq.1)then
c$OMP PARALLEL DO 
          do j=1,jtdm
            do i=1,itdm
              fldout(i,j)=rfldt(i,j)+fill
            enddo 
          enddo 
c$OMP END PARALLEL DO
          call ncerro(nf90_put_var(ncid,rhid,fldout,start,count))
        endif 
      enddo 
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        call ncerro(nf90_redef(ncid))
        write(lp,*) trim(vnm),arng
      endif
c
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
      integer :: i,j,ij,ijk,k,kd,n,prec
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),dimlen(maxdm),strn,strind(2,maxdm),
     .  msk(*)     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real, dimension(itdm*jtdm) :: fldout
      real(kind=4), dimension(itdm*jtdm) :: fldoutr4
      real :: fld(*),rkd,sfac,offs
      real, parameter :: fill=nf90_fill_double 
      real(kind=4) :: fillr4
c 
c --- Initialise fields
      rkd=1.
      kd=1
      do n=1,maxdm
        start(n)=1
        count(n)=1 
      enddo
c  
c --- define variable 
      if(mnproc.eq.1)then
        call ncsevl(dims,strn,strind)
        c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4.eq.'time') then
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
        enddo
        rkd=kd
        if (rec.eq.1) then 
          if (prec.ne.4) then         
            call ncerro(nf90_def_var(ncid,vnm,nf90_double,
     .       dimids(1:strn),rhid))
            call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fill))
          else 
            call ncerro(nf90_def_var(ncid,vnm,nf90_real,dimids(1:strn),
     .        rhid))
            fillr4=fill
            call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fillr4))
          endif 
          call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))
        else
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
        endif
        if(flgpad)then 
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif 
      endif
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- get mask for the full domain
      do j=1,jj
        do i=1,ii
          ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
          rmsk(i,j)=msk(ij)
        enddo
      enddo
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
        do j=1,jj
          do i=1,ii 
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            if (msk(ij).ne.0.and.fld(ijk).ne.fill) then 
              rfld(i,j)=fld(ijk)*sfac+offs
            else
              rfld(i,j)=fill
            endif
          enddo 
        enddo       
        call xcaget(rfldt,rfld,1) 
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
      enddo 
c
c --- Put file back to define mode
      if(mnproc.eq.1) call ncerro(nf90_redef(ncid))
c
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
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),strn,strind(2,maxdm),msk(*)     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real scf,sfac,ofs,offs,arng(2),fldmin,fldmax,fld(*),rkd,ruvflg
     .  ,rfill
      integer*2, dimension(itdm*jtdm) :: fldout
      integer*2 i2fill,i2max,i2min,vrng(2)
      parameter(rfill=nf90_fill_double,i2fill=-32768,i2max=32767,
     .  i2min=-i2max)
c
c --- Initialise fields
      ruvflg=0.
      rkd=1.
      kd=1
      vrng(1)=i2min
      vrng(2)=i2max
      do n=1,maxdm
        start(n)=1  
        count(n)=1 
      enddo
c 
c --- define variable 
      if(mnproc.eq.1)then
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
          if (dims(1:1).eq.'u'.or.dims(1:1).eq.'v') ruvflg=1.
        enddo
        rkd=kd
        if (rec.eq.1) then 
          call ncerro(nf90_def_var(ncid,vnm,nf90_int2,dimids(1:strn),
     .      rhid))
        else
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
        endif 
      endif
      call xcmaxr(ruvflg)
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- compute scale factor and offset
      fldmin=abs(nf90_fill_double)
      fldmax=-abs(nf90_fill_double)
c
      ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
c$OMP PARALLEL DO PRIVATE(ij,ijk) 
c$OMP+  REDUCTION(min:fldmin) REDUCTION(max:fldmax)
      do k=1,kd
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1) 
            ijk=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm
            if (msk(ij).eq.1.and.fld(ijk).lt.rfill) then
              fldmin=min(fldmin,(fld(ijk)*sfac)+offs)
              fldmax=max(fldmax,(fld(ijk)*sfac)+offs)
            endif
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcminr(fldmin)
      call xcmaxr(fldmax)
      if (ruvflg.gt.0.5) then 
        if (fldmin.ge.fldmax) then
          scf=1.d0
        else
          scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
        endif
        ofs=0.d0
      else
        if(fldmin.ge.fldmax)then
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
        call ncerro(nf90_put_att(ncid,rhid,'_FillValue',i2fill))
        call ncerro(nf90_put_att(ncid,rhid,'valid_range',vrng))
        call ncerro(nf90_put_att(ncid,rhid,'actual_range',arng))
        call ncerro(nf90_put_att(ncid,rhid,'scale_factor',scf))
        call ncerro(nf90_put_att(ncid,rhid,'add_offset',ofs))
        call ncerro(nf90_put_att(ncid,rhid,'compress','x y'))
        if(flgpad)then 
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif 
      endif
c
c --- get mask for the full domain
c$OMP PARALLEL DO PRIVATE(ij)
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
c$OMP+  PRIVATE(ij,ijk) 
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            if (msk(ij).eq.1.and.fld(ijk).lt.rfill) then
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
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        call ncerro(nf90_redef(ncid))
        write(lp,*) trim(vnm),arng
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
      integer :: n,ns,slenmax  
      integer, parameter :: maxdm=5 
      integer :: dimids(maxdm),strn,strind(2,maxdm)
c
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
     .      dimids(n)))
        enddo 
        if (strn.eq.1) then 
          ns=1
        else 
          call ncerro(nf90_inquire_dimension(ncid,dimids(1),
     .      len=slenmax))
          call ncerro(nf90_inquire_dimension(ncid,dimids(2),len=ns))
        endif
c
c --- - define variable
        call ncerro(nf90_def_var(ncid,vnm,nf90_char,dimids(1:strn),
     .    rhid))
c
c --- - leave define mode 
        if(flgpad)then 
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif 
c
c --- - write data 
        do n=1,ns
           call ncerro(nf90_put_var(ncid,rhid,fld(n)(1:slenmax)//'X',
     .       (/1,n/),(/slenmax,1/))) 
        end do 
          
c --- - switch back to define mode
        call ncerro(nf90_redef(ncid))
c
      endif
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
      integer, dimension(maxdm) :: istart,icount
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real(kind=4), dimension(itdm,jtdm) :: r4fldt
      real :: rkd,fld(*),sfac,offs,fldmin,fldmax
      real, parameter :: fill=nf90_fill_double
      character c4*4 
c
c --- Initialise fields
      rkd=1.
      kd=1
      do n=1,maxdm
        istart(n)=1  
        icount(n)=1  
      enddo
      icount(1)=itdm
      icount(2)=jtdm
c  
c --- define variable 
      if (mnproc.eq.1) then
        call ncsevl(dims,strn,strind)
        c4=dims(strind(1,strn):min(strind(1,strn)+3,strind(2,strn)))
        if (c4.eq.'time') then
          istart(strn)=rec
        else
          istart(strn)=1
        endif
        do n=1,strn
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .      dimid))
          dimids(n)=dimid
          if (n.gt.2.and.
     .        ((n.eq.strn.and.c4.ne.'time').or.
     .         (n.eq.strn-1.and.c4.eq.'time')))
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
        enddo
        rkd=kd
        if (rec.eq.1) then 
          if (prec.ne.4) then 
            call ncerro(nf90_def_var(ncid,vnm,nf90_double,
     .        dimids(1:strn),rhid))
          else 
            call ncerro(nf90_def_var(ncid,vnm,nf90_real,dimids(1:strn),
     .        rhid))
          endif 
          if(mskflg.eq.1)then
            if (prec.ne.4) then
              call ncerro(nf90_put_att(ncid,rhid,'_FillValue',
     .          nf90_fill_double))
            else
              call ncerro(nf90_put_att(ncid,rhid,'_FillValue',
     .          nf90_fill_real))
            endif
          endif
        else 
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
        endif 
        if(flgpad)then 
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif 
      endif
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- get mask for the full domain
      if(mskflg.eq.1)then 
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            rmsk(i,j)=msk(ij)
          enddo
        enddo
      else 
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=1.
          enddo
        enddo
      endif 
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field   
      fldmin=abs(nf90_fill_double)
      fldmax=-abs(nf90_fill_double)
      do k=1,kd
        if (k.gt.1) istart(3)=k 
c$OMP PARALLEL DO PRIVATE(ij,ijk)
c$OMP+  REDUCTION(min:fldmin) REDUCTION(max:fldmax)
        do j=1,jj
          do i=1,ii
            ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
            ijk=ij+(k-1)*ijdm
            rfld(i,j)=fld(ijk)
            if (msk(ij).ne.0.and.fld(ijk).ne.fill) then
              fldmin=min(rfld(i,j),fldmin)
              fldmax=max(rfld(i,j),fldmax)
            endif
          enddo
        enddo  
c$OMP END PARALLEL DO
        call xcaget(rfldt,rfld,1)  
        if (mnproc.eq.1) then
          if (prec.ne.4) then 
            if (sfac.ne.1..or.offs.ne.0.) then   
c$OMP PARALLEL DO 
              do j=1,jtdm
                do i=1,itdm
                  if(rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fill) then
                    rfldt(i,j)=fill
                  else
                    rfldt(i,j)=rfldt(i,j)*sfac+offs
                  endif
                enddo
              enddo
c$OMP END PARALLEL DO
            else
c$OMP PARALLEL DO 
              do j=1,jtdm
                do i=1,itdm
                  if(rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fill)
     .              rfldt(i,j)=fill
                enddo
              enddo
c$OMP END PARALLEL DO
            endif
            call ncerro(nf90_put_var(ncid,rhid,rfldt,istart,icount))
          else
            if (sfac.ne.1..or.offs.ne.0.) then   
c$OMP PARALLEL DO 
              do j=1,jtdm
                do i=1,itdm                
                  if(rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fill) then 
                    r4fldt(i,j)=fill
                  else 
                    r4fldt(i,j)=rfldt(i,j)*sfac+offs
                  endif
                enddo 
              enddo 
c$OMP END PARALLEL DO
            else 
c$OMP PARALLEL DO 
              do j=1,jtdm
                do i=1,itdm
                  if(rmskt(i,j).lt.0.5.or.rfldt(i,j).eq.fill) then
                    r4fldt(i,j)=fill
                  else
                    r4fldt(i,j)=rfldt(i,j)
                  endif 
                enddo
              enddo
c$OMP END PARALLEL DO
            endif
            call ncerro(nf90_put_var(ncid,rhid,r4fldt,istart,icount))
          endif
        endif 
      enddo 
      call xcminr(fldmin)
      call xcmaxr(fldmax)
c
c --- Put file back to define mode
      if (mnproc.eq.1) then
        call ncerro(nf90_redef(ncid))
        write(lp,*) trim(vnm),(fldmin*sfac+offs),(fldmax*sfac+offs)
      endif
c
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
      integer :: i,j,ij,k,kd,n
      integer, parameter :: maxdm=5,ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer :: dimid,dimids(maxdm),strn,strind(2,maxdm),msk(*),mskflg,
     .  fld(*)     
      integer, dimension(maxdm) :: start,count
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real :: rkd,fill
      character c4*4
c
c --- Initialise fields
      rkd=1.
      kd=1
      do n=1,maxdm
        start(n)=1  
        count(n)=1  
      enddo
      count(1)=itdm
      count(2)=jtdm
c  
c --- define variable 
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
     .      dimid))
          dimids(n)=dimid
          if (n.gt.2.and.
     .        ((n.eq.strn.and.c4.ne.'time').or.
     .         (n.eq.strn-1.and.c4.eq.'time')))
     .      call ncerro(nf90_inquire_dimension(ncid,dimid,len=kd))
        enddo
        rkd=kd
        if (rec.eq.1) then 
          call ncerro(nf90_def_var(ncid,vnm,nf90_int,dimids(1:strn),
     .      rhid))
        else 
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
        endif 
        if (mskflg.eq.1) then
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fill))
        endif
        if (flgpad) then 
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif 
      endif
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- get mask for the full domain
      if (mskflg.eq.1)then
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
          enddo
        enddo
      else
        do j=1,jj
          do i=1,ii
            rmsk(i,j)=1.
          enddo
        enddo
      endif
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
        if (k.gt.1) start(3)=k
        do j=1,jj
          do i=1,ii
            rfld(i,j)=fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
          enddo
        enddo
        call xcaget(rfldt,rfld,1)
        if (mnproc.eq.1) then
          do j=1,jtdm
            do i=1,itdm
              if(rmskt(i,j).lt.0.5) rfldt(i,j)=fill
            enddo
          enddo
          call ncerro(nf90_put_var(ncid,rhid,rfldt,start,count))
        endif
      enddo
c
c --- Put file back to define mode
      if (mnproc.eq.1) call ncerro(nf90_redef(ncid))
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
      integer :: n,nmax,strind(2,maxdm),ndims
      integer, dimension(maxdm) :: start,count,dimids,dimlen
c
      real :: fld(*),fill
c
      if(mnproc.eq.1)then
c
c --- - Define variable 
        fill=nf90_fill_double 
        call ncsevl(dims,ndims,strind)
        do n=1,ndims
          call ncerro(nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .      dimids(n)))
          call ncerro(nf90_inquire_dimension(ncid,dimids(n),
     .      len=dimlen(n)))
        enddo
        if (rec.eq.1) then
          call ncerro(nf90_def_var(ncid,vnm,nf90_double,dimids(1:ndims),
     .      rhid))
          call ncerro(nf90_put_att(ncid,rhid,'_FillValue',fill))
        else 
          call ncerro(nf90_inq_varid(ncid,vnm,rhid))
        endif
        if(flgpad)then
          call ncerro(nf90_enddef(ncid))
        else
          call ncerro(nf90_enddef(ncid,40960,4,0,4))
          flgpad=.true.
        endif
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
     .    count))
c
c --- - Put file back to define mode
        call ncerro(nf90_redef(ncid))
c
      endif
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
      if (ncstatus.ne.nf90_noerr) then 
        write(lp,*) 'NetCDF error: ',nf90_strerror(ncstatus)
        call flush(lp)
        call xchalt('(ncerro)')
        stop '(ncerro)'
      endif 
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
         if((charold.eq.' '.or.charold.eq.'-'.or.charold.eq.':').and.
     .      (charnew.ne.' '.and.charnew.ne.'-'.and.charnew.ne.':'))then
            strgn=strgn+1
            strgind(strgn)=i
         elseif((charnew.eq.' '.or.charnew.eq.'-'.or.charnew.eq.':')
     .    .and.(charold.ne.' '.and.charold.ne.'-'.and.charold.ne.':'))
     .   then
            strgn=strgn+1
            strgind(strgn)=i-1
         endif
         charold=charnew
      enddo
      if(mod(strgn,2).eq.1)then
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
      if(calendar(1:3).eq.'sta'.or.calendar(1:3).eq.'mix'.or.
     .   calendar(1:3).eq.'gre')then
         if(y1.lt.1582.or.(y1.eq.1582.and.m1.lt.10).or.
     .      (y1.eq.1582.and.m1.eq.10.and.d1.le.15))then
            call ncdnum('julian',y1,m1,d1,n1)
         else
            call ncdnum('julian',1582,10,15,nd1)
            call ncdnum('proleptic',1582,10,15,nd2)
            call ncdnum('proleptic',y1,m1,d1,nd3)
            n1=nd1-nd2+nd3-10
         endif
         if(y2.lt.1582.or.(y2.eq.1582.and.m2.lt.10).or.
     .      (y2.eq.1582.and.m2.eq.10.and.d2.le.15))then
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
      integer y,m,year,month,day,ndays,yr1,yr4,yr100,yr400
      integer acc_leap(12),acc_noleap(12),yoffset
      parameter(yoffset=1000000,yr1=365,yr4=1461,yr100=36524,
     .   yr400=146097)
      data acc_leap   /0,31,60,91,121,152,182,213,244,274,305,335/
      data acc_noleap /0,31,59,90,120,151,181,212,243,273,304,334/
c
      y=year+yoffset
      if(calendar(1:3).eq.'pro')then
         ndays=int((y-1)/400.)*yr400
     .      + int(mod((y-1),400)/100.)*yr100
     .      + int(mod(mod((y-1),400),100)/4.)*yr4
     .      + mod(mod(mod((y-1),400),100),4)*yr1
         if(month.gt.2.and.((mod(y,4).eq.0.and.mod(y,100).ne.0).or.
     .      mod(y,400).eq.0))then
            ndays=ndays+acc_leap(month)+day
         else
            ndays=ndays+acc_noleap(month)+day
         endif
      elseif(calendar(1:3).eq.'jul')then
         ndays=int((y-1)/4.)*yr4+mod(y,4)*yr1
         if(month.gt.2.and.mod(y,4).eq.0)then
            ndays=ndays+acc_leap(month)+day
         else
            ndays=ndays+acc_noleap(month)+day
         endif
      elseif(calendar(1:3).eq.'360')then
         ndays=(y-1)*360+30*(month-1)+day
      elseif(calendar(1:3).eq.'nol'.or.calendar(1:3).eq.'365')then
         ndays=(y-1)*365+acc_noleap(month)+day
      elseif(calendar(1:3).eq.'all'.or.calendar(1:3).eq.'366')then
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
      integer attlen     
c     
      status=nf90_inquire_attribute(ncid,nf90_global,trim(nm),
     .  len=attlen)
      if (status.ne.nf90_noerr) then
         ncinqa=.false.
      else
         if (attlen.ne.0) then
           ncinqa=.true.
         else 
           ncinqa=.false.
         endif 
      endif
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
      status=nf90_inq_varid(ncid,trim(vnm),rhid)
      if (status.ne.nf90_noerr) then 
         ncinqv=.false.
      else
         ncinqv=.true. 
      endif 
c
      end function ncinqv

      end module mod_nctools
