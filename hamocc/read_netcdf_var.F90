! Copyright (C) 2020  I. Bethke
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


      SUBROUTINE READ_NETCDF_VAR(ncid,desc,arr,klev,time,typeio)
!**************************************************************************
!
! Reads a variable from a NETCDF file and distributes it to all PEs
!
! The NETCDF File is only accessed by mnproc=1
!
!**************************************************************************
      USE netcdf
      USE mod_xc
      implicit none
#ifdef PNETCDF
#include <pnetcdf.inc>
#endif
#include <mpif.h>
      integer ncid, klev, time, ndims
      character (len=*) desc
      real arr(idm,jdm,klev),arr_g(itdm,jtdm)

      real, allocatable :: arr_l(:,:,:) 

      integer ncstat,ncvarid,i,j,k,typeio
      integer :: start(4),count(4)
      integer (kind=MPI_OFFSET_KIND) :: istart(4),icount(4)

! Read NETCDF data

      IF(TYPEIO==0) THEN
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
          write(lp,'(4a)') 'nf90_inq_varid: ',trim(desc),': ', &
     &                     nf90_strerror(ncstat)
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
            write(lp,'(4a)') 'nf90_get_vara_double: ',trim(desc),': ', &
     &                       nf90_strerror(ncstat)
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
      ELSE IF(TYPEIO==1) THEN
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
          write(lp,'(4a)') 'nfmpi_inq_varid: ',trim(desc),': ', &
     &                     nfmpi_strerror(ncstat)
          call xchalt('(read_pnetcdf_var)')
                 stop '(read_pnetcdf_var)'
        endif

      ncstat=nfmpi_get_vara_double_all(ncid,ncvarid,istart,icount,arr_l)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nfmpi_get_vara_double: ',trim(desc),': ', &
     &                       nfmpi_strerror(ncstat)
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
      ELSE
      call xchalt('(read_pnetcdf_var) WRONG IOTYPE')
      ENDIF
      END
