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

      subroutine rdcsss
c
c --- ------------------------------------------------------------------
c --- Read monthly climatological sea surface salinity
c --- ------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(3) :: start,count
      integer i,j,k,status,ncid,dimid,varid
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading monthly climatological SSS from ',
     .                    trim(scfile)
        call flush(lp)
c
c --- - open netcdf file
        status=nf90_open(trim(scfile),nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_open: ',trim(scfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
c
c --- - check dimensions
        status=nf90_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=i)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: x: ',
     .                     nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
        status=nf90_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=j)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimlen: y: ',nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in '//trim(scfile)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
c
        status=nf90_inq_varid(ncid,'sss',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: sss: ',
     .                     nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
c
        start(1)=1
        start(2)=1
        count(1)=itdm
        count(2)=jtdm
        count(3)=1
c
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: sss: ',
     .                       nf90_strerror(status)
            call xchalt('(rdcsss)')
                   stop '(rdcsss)'
          endif
        endif
        call xcaput(tmp2d,sssclm(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_close(ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_close: ',trim(scfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(rdcsss)')
                 stop '(rdcsss)'
        endif
      endif
c
      call xctilr(sssclm, 1,12, nbdy,nbdy, halo_ps)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'rdcsss:'
        endif
        call chksummsk(sssclm,ip,12,'sssclm')
      endif
c
      return
      end
