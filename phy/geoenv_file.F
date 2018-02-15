      subroutine geoenv_file
c
c --- ------------------------------------------------------------------
c --- Get bathymetry and grid specification from file and compute
c --- Coriolis parameter
c --- ------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
c
      real, dimension(itdm,jtdm) :: tmpg
      integer, dimension(3) :: start,count
      integer i,j,k,status,ncid,dimid,varid
c
c --- ------------------------------------------------------------------
c --- read grid information from grid.nc
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading grid information from ',
     .                    path(1:path_len)//'grid.nc'
        call flush(lp)
c
c --- - open netcdf file
        status=nf90_open(path(1:path_len)//'grid.nc',nf90_nowrite,
     .                   ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_open: grid.nc: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
c
c --- - check dimensions
        status=nf90_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_dimid: x: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=i)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inquire_dimension: x: ',
     .                      nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_dimid: y: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=j)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inquire_dimension: y: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in grid.nc'
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
c
c --- - read bathymetry
        status=nf90_inq_varid(ncid,'pdepth',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: pdepth: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: pdepth: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
c
c --- - count number of wet points for subsequent xcsum testing
        nwp=0
        do j=1,jtdm
          do i=1,itdm
            if (tmpg(i,j).gt.0.) nwp=nwp+1
          enddo
        enddo
      endif
c
      call xcaput(tmpg,depths,1)
c
c --- read grid coordinates
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qlon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qlon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: qlon: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,qlon,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qlat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qlat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: qlat: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,qlat,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'plon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: plon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: plon: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,plon,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'plat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: plat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: plat: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,plat,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'ulon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: ulon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: ulon: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,ulon,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'ulat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: ulat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: ulat: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,ulat,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'vlon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: vlon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: vlon: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,vlon,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'vlat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: vlat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: vlat: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,vlat,1)
c
      start(1)=1
      start(2)=1
      count(1)=itdm
      count(2)=jtdm
      count(3)=1
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qclon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qclon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: qclon: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,qclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qclat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qclat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: qclat: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,qclat(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'pclon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: pclon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: pclon: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,pclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'pclat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: pclat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: pclat: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,pclat(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'uclon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: uclon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: uclon: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,uclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'uclat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: uclat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: uclat: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,uclat(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'vclon',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: vclon: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: vclon: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,vclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'vclat',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: vclat: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmpg,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') 'nf90_get_var: vclat: ',
     .                       nf90_strerror(status)
            call xchalt('(geoenv_file)')
                   stop '(geoenv_file)'
          endif
        endif
        call xcaput(tmpg,vclat(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read scale factors
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qdx',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qdx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: qdx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scqx,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qdy',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qdy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: qdy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scqy,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'pdx',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: pdx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: pdx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scpx,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'pdy',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: pdy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: pdy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scpy,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'udx',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: udx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: udx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scux,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'udy',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: udy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: udy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scuy,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'vdx',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: vdx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: vdx: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scvx,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'vdy',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: vdy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: vdy: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scvy,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'qarea',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: qarea: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: qarea: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scq2,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'parea',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: parea: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: parea: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scp2,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'uarea',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: uarea: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: uarea: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scu2,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'varea',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: varea: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: varea: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,scv2,1)
c
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'angle',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_inq_varid: angle: ',
     .                     nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf90_get_var(ncid,varid,tmpg)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_var: angle: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
      call xcaput(tmpg,angle,1)
c
c --- close grid information file
c
      if (mnproc.eq.1) then
        status=nf90_close(ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_close: grid.nc: ',nf90_strerror(status)
          call xchalt('(geoenv_file)')
                 stop '(geoenv_file)'
        endif
      endif
c
c --- ------------------------------------------------------------------
c --- Get correct units of scale factors, precompute cosine and sine of
c --- local angle of i-direction and with eastward direction, and
c --- compute Coriolis and beta plane parameter
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do i=1,ii
c
          scqx(i,j)=scqx(i,j)*1.e2
          scqy(i,j)=scqy(i,j)*1.e2
          scpx(i,j)=scpx(i,j)*1.e2
          scpy(i,j)=scpy(i,j)*1.e2
          scux(i,j)=scux(i,j)*1.e2
          scuy(i,j)=scuy(i,j)*1.e2
          scvx(i,j)=scvx(i,j)*1.e2
          scvy(i,j)=scvy(i,j)*1.e2
          scq2(i,j)=scq2(i,j)*1.e4
          scp2(i,j)=scp2(i,j)*1.e4
          scu2(i,j)=scu2(i,j)*1.e4
          scv2(i,j)=scv2(i,j)*1.e4
c
          cosang(i,j)=cos(angle(i,j))
          sinang(i,j)=sin(angle(i,j))
c
          corioq(i,j)=sin(qlat(i,j)/radian)*4.*pi/86164.
          coriop(i,j)=sin(plat(i,j)/radian)*4.*pi/86164.
          betafp(i,j)=cos(plat(i,j)/radian)*4.*pi/(86164.*rearth)
c
        enddo
      enddo
c
      return
      end