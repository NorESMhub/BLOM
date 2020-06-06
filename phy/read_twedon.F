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

      subroutine read_twedon
c
c --- ------------------------------------------------------------------
c --- Read tidal wave energy dissipation divided by bottom buoyancy
c --- frequency
c --- ------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
#include "common_blocks.h"
c
      real, dimension(itdm,jtdm) :: tmpg
      integer i,j,l,status,ncid,dimid,varid
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading tidal dissipation data from ',
     .                    trim(tdfile)
        call flush(lp)
c
c --- - open netcdf file
        status=nf90_open(tdfile,nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_open: ',trim(tdfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
c
c --- - check dimensions
        status=nf90_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=i)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: x: ',
     .                     nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
        status=nf90_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=j)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: y: ',
     .                     nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,'(2a)') ' wrong dimensions in ',trim(tdfile)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
c
c --- - read tidal dissipation data
        status=nf90_inq_varid(ncid,'twedon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: twedon: ',
     .                     nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: twedon: ',
     .                     nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
c
        status=nf90_close(ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_close: ',trim(tdfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(read_twedon)')
                 stop '(read_twedon)'
        endif
      endif
c
      call xcaput(tmpg,twedon,1)
c
c --- change units from W*s/m**2=kg/s**2 to g/s**2
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          twedon(i,j)=twedon(i,j)*1.e3
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(twedon, 1,1, nbdy,nbdy, halo_ps)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'read_twedon:'
        endif
        call chksummsk(twedon,ip,1,'twedon')
      endif
c
      return
      end
