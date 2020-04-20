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


      SUBROUTINE WRITE_NETCDF_VAR(ncid,desc,arr,klev,time)
!**************************************************************************
!
! Gathers a global variable from all PEs and writes it to a NETCDF file
!
! The NETCDF File is only accessed by mnproc=1
!
!**************************************************************************
      use netcdf
      use mod_xc
      use mod_dia, only : iotype      
      implicit none
#ifdef PNETCDF
#  include <pnetcdf.inc>
#endif
#  include <mpif.h>
      integer ncid, klev, time, ndims
      character (len=*) desc
      real arr(idm,jdm,klev)

      real arr_g(itdm,jtdm)
      real , allocatable :: arr_g1(:,:,:),arr_l(:,:,:)
      integer ncstat,ncvarid,k,i,j
      integer, allocatable :: start(:),count(:) 
      integer (kind=MPI_OFFSET_KIND), allocatable :: istart(:),icount(:)
! Write NETCDF data
       
      if (klev.eq.1.and.time.eq.0) then
        ndims=2
      elseif (klev.eq.1.or.time.eq.0) then
        ndims=3
      else
        ndims=4
      endif
      IF(IOTYPE==0) THEN
      allocate(start(ndims))
      allocate(count(ndims))
      allocate(arr_l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      arr_l=0.0
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
            write(lp,'(4a)') 'nf90_inq_varid: ',trim(desc),': ', &
     &                       nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
          endif
          ncstat=nf90_put_var(ncid,ncvarid,arr_g,start,count)
          if (ncstat.ne.nf90_noerr) then
            write(lp,'(4a)') 'nf90_put_var: ',trim(desc),': ', &
     &                       nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
          endif
!          ncstat=nf90_sync(ncid)
!          if (ncstat.ne.nf90_noerr) then
!            write(lp,'(4a)') 'nf90_sync: ',trim(desc),': ', &
!     &                       nf90_strerror(ncstat)
!            call xchalt('(write_netcdf_var)')
!                   stop '(write_netcdf_var)'
!          endif
        endif
      enddo
      deallocate(start,count)
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      allocate(istart(ndims))
      allocate(icount(ndims))
      allocate(arr_l(ii,jj,klev))
      arr_l=0.0
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
       allocate(arr_g1(itdm,jj,klev))
       arr_g1=0.0
       call xcgetrow(arr_g1, arr_l, klev)

          ncstat=nfmpi_inq_varid(ncid,desc,ncvarid)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nfmpi_inq_varid: ',trim(desc),': ', &
     &                       nfmpi_strerror(ncstat)
            call xchalt('(write_pnetcdf_var)')
                   stop '(write_pnetcdf_var)'
          endif
        ncstat=nfmpi_put_vara_double_all(ncid,ncvarid,istart,     &
     &           icount,arr_g1)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nfmpi_put_var: ',trim(desc),': ', &
     &                       nfmpi_strerror(ncstat)
            call xchalt('(write_pnetcdf_var)')
                   stop '(write_pnetcdf_var)'
          endif
!          ncstat=nfmpi_sync(ncid)
!          if (ncstat.ne.nf_noerr) then
!            write(lp,'(4a)') 'nfmpi_sync: ',trim(desc),': ', &
!     &                       nfmpi_strerror(ncstat)
!            call xchalt('(write_pnetcdf_var)')
!                   stop '(write_pnetcdf_var)'
!          endif

          deallocate(istart,icount,arr_g1)
#endif
      ENDIF

      END
