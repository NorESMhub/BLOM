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

      subroutine geoenv_cesmextra
c
c --- ------------------------------------------------------------------
c --- Read additional grid parameters when configured for coupling to
c --- CESM
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
      integer i,j,status,ncid,dimid,varid
c
c --- ------------------------------------------------------------------
c --- read grid information from grid file
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading additional grid information from ',
     .                    trim(grfile)
        call flush(lp)
c
c --- - open netcdf file
        status=nf90_open(grfile,nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_open: ',trim(grfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
c
c --- - check dimensions
        status=nf90_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=i)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: x: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
        status=nf90_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=j)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: y: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,'(2a)') ' wrong dimensions in ',trim(grfile)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
      endif
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'cplmask',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: cplmask: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_get_var: cplmask: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
      endif
      call xcaput(tmpg,util1,1)
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          cplmsk(i,j)=nint(util1(i,j))
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- close grid information file
c
      if (mnproc.eq.1) then
        status=nf90_close(ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_close: ',trim(grfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_cesmextra)')
                 stop '(geoenv_cesmextra)'
        endif
      endif
c
      return
      end
