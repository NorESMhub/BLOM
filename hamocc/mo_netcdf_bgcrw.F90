! Copyright (C) 2001  S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen
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

module mo_netcdf_bgcrw

  implicit none
  private

  public :: netcdf_def_vardb
  public :: write_netcdf_var
  public :: read_netcdf_var

contains

  subroutine netcdf_def_vardb (kcid,kshort,yshort,kdims,kcdims,kcvarid,                            &
                               kunitl,yunit,klong,ylong,pmissing,klabel,kunit)

    ! **********************************************************************************************
    !  Interface to NETCDF routines - define NetCDF variable.
    !
    !  S.Legutke,        *MPI-MaD, HH*    10.10.01
    ! **********************************************************************************************

    use netcdf,  only: nf90_double,nf90_noerr,nf90_put_att,nf90_def_var
    use mod_xc,  only: mnproc,xchalt
    use mod_dia, only: iotype
    use mo_kind, only: rp
#ifdef PNETCDF
#   include <pnetcdf.inc>
#   include <mpif.h>
#endif

    ! Arguments
    integer,          intent(in)  :: kcid          ! file ID.
    integer,          intent(in)  :: kshort        ! length of short name.
    integer,          intent(in)  :: kdims         ! number of dimensions.
    integer,          intent(in)  :: kcdims(kdims) ! dimensions.
    integer,          intent(out) :: kcvarid       ! variable ID.
    integer,          intent(in)  :: kunitl        ! length of unit string.
    integer,          intent(in)  :: klong         ! length of long name.
    integer,          intent(in)  :: klabel        ! label for abort identification.
    integer,          intent(in)  :: kunit         ! stdout unit.
    character(len=*), intent(in)  :: yshort        ! short name.
    character(len=*), intent(in)  :: yunit         ! unit string.
    character(len=*), intent(in)  :: ylong         ! long name.

    ! Local variables
    integer                       :: k
    real(rp)                      :: pmissing
    character(len=24)             :: ystring
    integer                       :: ncstat
#ifdef PNETCDF
    integer(kind=MPI_OFFSET_KIND) :: clen
#endif

    ystring(1:21)='NETCDF stop at label '

    !
    !  Define variable
    !
    if (mnproc==1 .and. IOTYPE==0) then
      ncstat = NF90_DEF_VAR(kcid,yshort(1:kshort),NF90_DOUBLE,kcdims,kcvarid)
      if ( ncstat /=  NF90_NOERR ) then
        write(kunit,*) 'Problems with definition of NetCDF variable:'
        write(kunit,*) 'kcid           : ',kcid
        write(kunit,*) 'kshort         : ',kshort
        write(kunit,*) 'yshort(kshort) : ',yshort(1:kshort),'---'
        write(kunit,*) 'kdims          : ',kdims
        write(kunit,*) 'kcdims         : ',(kcdims(k),k=1,kdims)
        write(kunit,*) 'kcvarid        : ',kcvarid
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(netcdf_def_vardb)')
        stop '(netcdf_def_vardb)'
      endif
      !
      !  Set unit
      !
      ncstat = NF90_PUT_ATT(kcid,kcvarid,'units',yunit(1:kunitl))
      if ( ncstat /=  NF90_NOERR ) then
        write(kunit,*) 'Problems with definition of unit:'
        write(kunit,*) 'kcid          : ',kcid
        write(kunit,*) 'kcvarid       : ',kcvarid
        write(kunit,*) 'kunitl        : ',kunitl
        write(kunit,*) 'yunit(kunitl) : ',yunit(1:kunitl),'---'
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(netcdf_def_vardb)')
        stop '(netcdf_def_vardb)'
      endif
      !
      !  Set long name
      !
      ncstat = NF90_PUT_ATT(kcid,kcvarid,'long_name',ylong(1:klong))
      if ( ncstat /=  NF90_NOERR ) then
        write(kunit,*) 'Problems with definition of long name:'
        write(kunit,*) 'kcid         : ',kcid
        write(kunit,*) 'kcvarid      : ',kcvarid
        write(kunit,*) 'klong        : ',klong
        write(kunit,*) 'ylong(klong) : ',ylong(1:klong),'---'
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(netcdf_def_vardb)')
        stop '(netcdf_def_vardb)'
      endif
      !
      !  Set missing value
      !
      ncstat = NF90_PUT_ATT(kcid,kcvarid,'missing_value',pmissing)
      if ( ncstat /=  NF90_NOERR ) then
        write(kunit,*) 'Problems with definition of missing value:'
        write(kunit,*) 'kcid     : ',kcid
        write(kunit,*) 'kcvarid  : ',kcvarid
        write(kunit,*) 'pmissing : ',pmissing
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(netcdf_def_vardb)')
        stop '(netcdf_def_vardb)'
      endif
    else if (IOTYPE==1) then
#ifdef PNETCDF
      !
      !  Define variable
      !
      ncstat = nfmpi_def_var(kcid,yshort(1:kshort),nf_double,kdims,kcdims,kcvarid)
      if ( ncstat /=  NF_NOERR ) then
        write(kunit,*) 'Problems with definition of NetCDF variable:'
        write(kunit,*) 'kcid           : ',kcid
        write(kunit,*) 'kshort         : ',kshort
        write(kunit,*) 'yshort(kshort) : ',yshort(1:kshort),'---'
        write(kunit,*) 'kdims          : ',kdims
        write(kunit,*) 'kcdims         : ',(kcdims(k),k=1,kdims)
        write(kunit,*) 'kcvarid        : ',kcvarid
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(pnetcdf_def_vardb)')
        stop '(pnetcdf_def_vardb)'
      endif
      !
      !  Set unit
      !
      clen=len(trim(yunit(1:kunitl)))
      ncstat = NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'units',clen,yunit(1:kunitl))
      if ( ncstat /=  NF_NOERR ) then
        write(kunit,*) 'Problems with definition of unit:'
        write(kunit,*) 'kcid          : ',kcid
        write(kunit,*) 'kcvarid       : ',kcvarid
        write(kunit,*) 'kunitl        : ',kunitl
        write(kunit,*) 'yunit(kunitl) : ',yunit(1:kunitl),'---'
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(pnetcdf_def_vardb)')
        stop '(pnetcdf_def_vardb)'
      endif
      !
      !  Set long name
      !
      clen=len(trim(ylong(1:klong)))
      ncstat = NFMPI_PUT_ATT_TEXT(kcid,kcvarid,'long_name',clen,ylong(1:klong))
      if ( ncstat /=  NF_NOERR ) then
        write(kunit,*) 'Problems with definition of long name:'
        write(kunit,*) 'kcid         : ',kcid
        write(kunit,*) 'kcvarid      : ',kcvarid
        write(kunit,*) 'klong        : ',klong
        write(kunit,*) 'ylong(klong) : ',ylong(1:klong),'---'
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(pnetcdf_def_vardb)')
        stop '(pnetcdf_def_vardb)'
      endif
      !
      !  Set missing value
      !
      clen=1
      ncstat = NFMPI_PUT_ATT_DOUBLE(kcid,kcvarid,'missing_value',NF_DOUBLE,clen,pmissing)
      if ( ncstat /=  NF_NOERR ) then
        write(kunit,*) 'Problems with definition of missing value:'
        write(kunit,*) 'kcid     : ',kcid
        write(kunit,*) 'kcvarid  : ',kcvarid
        write(kunit,*) 'pmissing : ',pmissing
        write(ystring(22:24),'(I3)') klabel
        write(kunit,*) ystring
        call xchalt('(pnetcdf_def_vardb)')
        stop '(pnetcdf_def_vardb)'
      endif
#endif
    endif

  end subroutine netcdf_def_vardb


  subroutine write_netcdf_var(ncid,desc,arr,klev,time)

    !***********************************************************************************************
    ! Gathers a global variable from all PEs and writes it to a NETCDF file
    ! The NETCDF File is only accessed by mnproc=1
    !***********************************************************************************************

    use netcdf,  only: nf90_noerr,nf90_inq_varid,nf90_strerror,nf90_put_var
    use mod_xc,  only: itdm,jtdm,jdm,lp,mnproc,nbdy,idm,xchalt,xcaget
    use mod_dia, only: iotype
    use mo_kind, only: rp
#ifdef PNETCDF
    use mod_xc,  only: i0,ii,jj,j0,mproc,mpe_1,nproc,xcgetrow
#   include <pnetcdf.inc>
#   include <mpif.h>
#endif

    ! Arguments
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: desc
    integer,          intent(in)  :: klev
    integer,          intent(in)  :: time
    real(rp),         intent(in)  :: arr(idm,jdm,klev)

    ! Local variables
    integer                       :: k,i,j
    integer                       :: ndims
    real(rp)                      :: arr_g(itdm,jtdm)
    integer                       :: ncstat
    integer                       :: ncvarid
    integer, allocatable          :: start(:),count(:)
    real(rp),    allocatable      :: arr_l(:,:,:)
#ifdef PNETCDF
    real(rp),    allocatable      :: arr_g1(:,:,:)
    integer (kind=MPI_OFFSET_KIND), allocatable :: istart(:),icount(:)
#endif

    ! Write NETCDF data

    if (klev.eq.1.and.time.eq.0) then
      ndims=2
    elseif (klev.eq.1.or.time.eq.0) then
      ndims=3
    else
      ndims=4
    endif
    if (IOTYPE==0) then

      allocate(start(ndims))
      allocate(count(ndims))
      allocate(arr_l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      arr_l=0.0_rp
      start(1)=1
      count(1)=itdm
      start(2)=1
      count(2)=jtdm
      if (klev.gt.1.or.time.gt.0) then
        if (klev.gt.1.and.time.gt.0) then
          start(3)=1
          count(3)=1
          start(4)=time
          count(4)=1
        else if (klev.gt.1.and.time.eq.0) then
          start(3)=1
          count(3)=1
        else
          start(3)=time
          count(3)=1
        endif
      endif

      do k=1,klev
        do j=1,jdm
          do i=1,idm
            arr_l(i,j,1)=arr(i,j,k)
          enddo
        enddo
        call xcaget(arr_g,arr_l,1)
        if (mnproc.eq.1) then
          if (k.gt.1) then
            start(3)=k
            count(3)=1
          endif

          ncstat=nf90_inq_varid(ncid,desc,ncvarid)
          if (ncstat.ne.nf90_noerr) then
            write(lp,'(4a)') 'nf90_inq_varid: ',trim(desc),': ',nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
            stop        '(write_netcdf_var)'
          endif

          ncstat=nf90_put_var(ncid,ncvarid,arr_g,start,count)
          if (ncstat.ne.nf90_noerr) then
            write(lp,'(4a)') 'nf90_put_var: ',trim(desc),': ',nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
            stop        '(write_netcdf_var)'
          endif

          ! ncstat=nf90_sync(ncid)
          ! if (ncstat.ne.nf90_noerr) then
          !   write(lp,'(4a)') 'nf90_sync: ',trim(desc),': ',nf90_strerror(ncstat)
          !   call xchalt('(write_netcdf_var)')
          !   stop '(write_netcdf_var)'
          ! endif
        endif
      enddo
      deallocate(start,count,arr_l)

    else if (IOTYPE==1) then

#ifdef PNETCDF
      allocate(istart(ndims))
      allocate(icount(ndims))
      allocate(arr_l(ii,jj,klev))
      allocate(arr_g1(itdm,jj,klev))

      arr_l=0.0_rp
      arr_g1=0.0_rp

      if (klev.gt.1.or.time.gt.0) then
        if (klev.gt.1.and.time.gt.0) then
          istart(3)=1
          icount(3)=klev
          istart(4)=time
          icount(4)=1
        else if (klev.gt.1.and.time.eq.0) then
          istart(3)=1
          icount(3)=klev
        else
          istart(3)=time
          icount(3)=1
        endif
      endif

      istart(1)=1
      istart(2)=j0+1

      if(mproc .eq. mpe_1(nproc) ) then
        icount(1)=itdm
        icount(2)=jj
      else
        do i=1,ndims
          icount(i)=0
        enddo
      endif

      do k=1,klev
        do j=1,jj
          do i=1,ii
            arr_l(i,j,k)=arr(i,j,k)
          enddo
        enddo
      enddo

      call xcgetrow(arr_g1, arr_l, klev)

      ncstat=nfmpi_inq_varid(ncid,desc,ncvarid)
      if (ncstat.ne.nf_noerr) then
        write(lp,'(4a)') 'nfmpi_inq_varid: ',trim(desc),': ',nfmpi_strerror(ncstat)
        call xchalt('(write_pnetcdf_var)')
        stop        '(write_pnetcdf_var)'
      endif

      ncstat=nfmpi_put_vara_double_all(ncid,ncvarid,istart,icount,arr_g1)
      if (ncstat.ne.nf_noerr) then
        write(lp,'(4a)') 'nfmpi_put_var: ',trim(desc),': ',nfmpi_strerror(ncstat)
        call xchalt('(write_pnetcdf_var)')
        stop        '(write_pnetcdf_var)'
      endif

      ! ncstat=nfmpi_sync(ncid)
      ! if (ncstat.ne.nf_noerr) then
      !   write(lp,'(4a)') 'nfmpi_sync: ',trim(desc),': ',nfmpi_strerror(ncstat)
      !   call xchalt('(write_pnetcdf_var)')
      !   stop '(write_pnetcdf_var)'
      ! endif

      deallocate(istart,icount,arr_l,arr_g1)
#endif

    end if

  end subroutine write_netcdf_var


  subroutine read_netcdf_var(ncid,desc,arr,klev,time,typeio)

    !***********************************************************************************************
    ! Reads a variable from a NETCDF file and distributes it to all PEs
    ! The NETCDF File is only accessed by mnproc=1
    !***********************************************************************************************

    use netcdf, only: nf90_noerr,nf90_inq_varid,nf90_strerror,nf90_get_var
    use mod_xc, only: idm,itdm,jtdm,jdm,lp,mnproc,nbdy,xchalt,xcaput
    use mo_kind,only: rp
#ifdef PNETCDF
    use mod_xc, only: i0,ii,jj,j0
#   include <pnetcdf.inc>
#   include <mpif.h>
#endif

    ! Arguments
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: desc
    integer,          intent(in)  :: klev
    integer,          intent(in)  :: time
    integer,          intent(in)  :: typeio
    real(rp),         intent(out) :: arr(idm,jdm,klev)

    ! Local variables
    integer                       :: i,j,k
    integer                       :: ncstat
    integer                       :: ncvarid
    integer                       :: start(4),count(4)
    real(rp)                      :: arr_g(itdm,jtdm)
    real(rp), allocatable         :: arr_l(:,:,:)
#ifdef PNETCDF
    integer(kind=MPI_OFFSET_KIND) :: istart(4),icount(4)
#endif

    ! Read NETCDF data

    if (TYPEIO==0) then

      start=1
      count=0
      start(1)=1
      count(1)=itdm
      start(2)=1
      count(2)=jtdm
      if (klev.gt.1.or.time.gt.0) then
        if (klev.gt.1.and.time.gt.0) then
          start(3)=1
          count(3)=1
          start(4)=time
          count(4)=1
        else if (klev.gt.1.and.time.eq.0) then
          start(3)=1
          count(3)=1
        else
          start(3)=time
          count(3)=1
        endif
      endif
      allocate(arr_l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))

      if (mnproc.eq.1) then
        ncstat=nf90_inq_varid(ncid,desc,ncvarid)
        if (ncstat.ne.nf90_noerr) then
          write(lp,'(4a)') 'nf90_inq_varid: ',trim(desc),': ',nf90_strerror(ncstat)
          call xchalt('(read_netcdf_var)')
          stop        '(read_netcdf_var)'
        endif
      endif
      do k=1,klev
        if (mnproc.eq.1) then
          if (k.gt.1) then
            start(3)=k
            count(3)=1
          endif
          ncstat=nf90_get_var(ncid,ncvarid,arr_g,start,count)
          if (ncstat.ne.nf90_noerr) then
            write(lp,'(4a)') 'nf90_get_vara_double: ',trim(desc),': ',nf90_strerror(ncstat)
            call xchalt('(read_netcdf_var)')
            stop        '(read_netcdf_var)'
          endif
        endif
        call xcaput(arr_g,arr_l,1)
        do j=1,jdm
          do i=1,idm
            arr(i,j,k)=arr_l(i,j,1)
          enddo
        enddo
      enddo

    else if (TYPEIO==1) then

#ifdef PNETCDF
      allocate(arr_l(ii,jj,klev))
      arr=0.0_rp
      istart=1
      icount=0
      istart(1)=i0+1
      icount(1)=ii
      istart(2)=j0+1
      icount(2)=jj
      if (klev.gt.1.or.time.gt.0) then
        if (klev.gt.1.and.time.gt.0) then
          istart(3)=1
          icount(3)=klev
          istart(4)=time
          icount(4)=1
        else if (klev.gt.1.and.time.eq.0) then
          istart(3)=1
          icount(3)=klev
        else
          istart(3)=time
          icount(3)=1
        endif
      endif

      ncstat=nfmpi_inq_varid(ncid,desc,ncvarid)
      if (ncstat.ne.nf_noerr) then
        write(lp,'(4a)') 'nfmpi_inq_varid: ',trim(desc),': ',nfmpi_strerror(ncstat)
        call xchalt('(read_pnetcdf_var)')
        stop        '(read_pnetcdf_var)'
      endif

      ncstat=nfmpi_get_vara_double_all(ncid,ncvarid,istart,icount,arr_l)
      if (ncstat.ne.nf_noerr) then
        write(lp,'(4a)') 'nfmpi_get_vara_double: ',trim(desc),': ',nfmpi_strerror(ncstat)
        call xchalt('(read_pnetcdf_var)')
        stop        '(read_pnetcdf_var)'
      endif
      do k=1,klev
        do j=1,jj
          do i=1,ii
            arr(i,j,k)=arr_l(i,j,k)
          enddo
        enddo
      enddo
#endif
    else
      call xchalt('(read_pnetcdf_var) WRONG IOTYPE')
    endif

    deallocate(arr_l)

  end subroutine read_netcdf_var

end module mo_netcdf_bgcrw
