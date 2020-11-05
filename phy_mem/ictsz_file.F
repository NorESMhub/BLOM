! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen
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

      subroutine ictsz_file
c
c --- ------------------------------------------------------------------
c --- Read initial conditions from file to define layer temperature and
c --- salinity and geopotential at layer interfaces.
c --- ------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dz
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(3) :: start,count
      integer i,j,k,l,status,ncid,dimid,varid
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading initial condition from ',
     .                    trim(icfile)
        call flush(lp)
c
c --- - Open netcdf file
        status=nf90_open(icfile,nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_open: ',trim(icfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
c
c --- - Check dimensions
        status=nf90_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=i)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: x: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        status=nf90_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=j)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: y: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in '//trim(icfile)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
c
      endif
c
      start(1)=1
      start(2)=1
      count(1)=itdm
      count(2)=jtdm
      count(3)=1
c
c --- Read reference potential density
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'sigma',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: sigma: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: sigma: ',
     .                       nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,sigmar(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- Read potential temperature
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'temp',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: temp: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: temp: ',
     .                       nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,temp(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- Read salinity
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'saln',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: saln: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: saln: ',
     .                       nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,saln(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- Read layer thickness
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'dz',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: dz: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: dz: ',nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,dz(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_close(ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_close: ',trim(icfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
c
c --- Construct interface depths [cm] from layer thicknesses [m] and
c --- convert unit of reference potential density from [kg/m^3] to
c --- [g/cm^3]
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c         z(i,j,1)=z(i,j,1)*1.e2
          z(i,j,1)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      do k=1,kk
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            z(i,j,k+1)=min(depths(i,j)*1.e2,z(i,j,k)+dz(i,j,k)*1.e2)
            sigmar(i,j,k)=sigmar(i,j,k)*1.e-3
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=2,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (z(i,j,kk+1)-z(i,j,k).lt.1.e-4)
     .        z(i,j,k)=depths(i,j)*1.e2
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          z(i,j,kk+1)=depths(i,j)*1.e2
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute layer interface geopotential
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk+1
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            phi(i,j,k)=-g*z(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
