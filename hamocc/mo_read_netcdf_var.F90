! Copyright (C) 2020  I. Bethke, M. Bentsen
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

module mo_read_netcdf_var

  implicit none
  private

  public :: read_netcdf_var

contains

  subroutine read_netcdf_var(ncid,desc,arr,klev,time,typeio)

    !**************************************************************************
    ! Reads a variable from a NETCDF file and distributes it to all PEs
    ! The NETCDF File is only accessed by mnproc=1
    !**************************************************************************

    use netcdf, only: nf90_noerr,nf90_inq_varid,nf90_strerror,nf90_get_var
    use mod_xc, only: idm,itdm,jtdm,jdm,lp,mnproc,nbdy,xchalt,xcaput
#ifdef PNETCDF
    use mod_xc, only: i0,ii,jj,j0
#endif
#ifdef PNETCDF
#include <pnetcdf.inc>
#include <mpif.h>
#endif

    ! Arguments
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: desc
    integer,          intent(in)  :: klev
    integer,          intent(in)  :: time
    integer,          intent(in)  :: typeio
    real,             intent(out) :: arr(idm,jdm,klev)

    ! Local variables
    integer :: i,j,k
    integer :: ncstat
    integer :: ncvarid
    integer :: start(4),count(4)
    real    :: arr_g(itdm,jtdm)
    real, allocatable :: arr_l(:,:,:)
#ifdef PNETCDF
    integer (kind=MPI_OFFSET_KIND) :: istart(4),icount(4)
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
          stop '(read_netcdf_var)'
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
            stop '(read_netcdf_var)'
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
      arr=0.0
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
        stop '(read_pnetcdf_var)'
      endif

      ncstat=nfmpi_get_vara_double_all(ncid,ncvarid,istart,icount,arr_l)
      if (ncstat.ne.nf_noerr) then
        write(lp,'(4a)') 'nfmpi_get_vara_double: ',trim(desc),': ',nfmpi_strerror(ncstat)
        call xchalt('(read_pnetcdf_var)')
        stop '(read_pnetcdf_var)'
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

  end subroutine read_netcdf_var

end module mo_read_netcdf_var
