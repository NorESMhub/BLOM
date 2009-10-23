      SUBROUTINE READ_NETCDF_VAR(ncid,desc,arr,klev,time)

!**************************************************************************
!
! Reads a variable from a NETCDF file and distributes it to all PEs
!
! The NETCDF File is only accessed by mnproc=1
!
!**************************************************************************

      USE mod_xc
      
      implicit none

      integer ncid, klev, time, ndims
      character (len=*) desc
      real arr(idm,jdm,klev)

      real arr_g(itdm,jtdm),arr_l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

      include 'netcdf.inc'
      integer ncstat,ncvarid,i,j,k
      integer, allocatable :: start(:),count(:)

! Read NETCDF data

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
 
      if (mnproc.eq.1) then
        ncstat=nf_inq_varid(ncid,desc,ncvarid)
        if (ncstat.ne.nf_noerr) then
          write(lp,'(4a)') 'nf_inq_varid: ',trim(desc),': ', &
     &                     nf_strerror(ncstat)
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
          ncstat=nf_get_vara_double(ncid,ncvarid,start,count,arr_g)
          if (ncstat.ne.nf_noerr) then
            write(lp,'(4a)') 'nf_get_vara_double: ',trim(desc),': ', &
     &                       nf_strerror(ncstat)
            call xchalt('(read_netcdf_var)')
                   stop '(read_netcdf_var)'
          endif
        endif
        call xcaput(arr_g,arr_l,1)
        do j=1,jdm
          do i=1,idm
            arr(i,j,k)=arr_l(i,j)
          enddo
        enddo
      enddo

      END
