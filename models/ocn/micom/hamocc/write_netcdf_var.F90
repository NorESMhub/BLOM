      SUBROUTINE WRITE_NETCDF_VAR(ncid,desc,arr,klev,time)

!**************************************************************************
!
! Gathers a global variable from all PEs and writes it to a NETCDF file
!
! The NETCDF File is only accessed by mnproc=1
!
!**************************************************************************

      use mod_xc
      
      implicit none

      integer ncid, klev, time, ndims
      character (len=*) desc
      real arr(idm,jdm,klev)

      real arr_g(itdm,jtdm),arr_l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

      include 'netcdf.inc'
      integer ncstat,ncvarid,k,i,j
      integer, allocatable :: start(:),count(:) 

! Write NETCDF data

      if (klev.eq.1.and.time.eq.0) then
        ndims=2
      elseif (klev.eq.1.or.time.eq.0) then
        ndims=3
      else
        ndims=4
      endif
      allocate(start(ndims))
      allocate(count(ndims))

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
            arr_l(i,j)=arr(i,j,k)
          enddo
        enddo
        call xcaget(arr_g,arr_l,1)
        if (mnproc.eq.1) then
          if (k.gt.1) then
            start(3)=k
            count(3)=1
          endif
          ncstat=nf_inq_varid(ncid,desc,ncvarid)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nf_inq_varid: ',trim(desc),': ', &
     &                       nf_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
          endif
          ncstat=nf_put_vara_double(ncid,ncvarid,start,count,arr_g)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nf_put_vara_double: ',trim(desc),': ', &
     &                       nf_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
          endif
          ncstat=nf_sync(ncid)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nf_sync: ',trim(desc),': ', &
     &                       nf_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
          endif
        endif
      enddo

      deallocate(start,count)

      END
