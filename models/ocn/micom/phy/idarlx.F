      subroutine idarlx
c
c --- ------------------------------------------------------------------
c --- Initialize diagnosing/application of relaxation fluxes
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c
      real, dimension(itdm,jtdm) :: tmp2d
      integer i,j,k
c
      if (aptflx) then
        if (mnproc.eq.1) then
          open (unit=nfu,file=path(1:path_len)//'tflxdi.uf',
     .          form='unformatted')
          read (nfu) i,j
          if (i.ne.itdm.or.j.ne.jtdm) then
            write (lp,*) 'wrong dimensions in tflxdi.uf'
            call xchalt('(idarlx)')
                   stop '(idarlx)'
          endif
        endif
        do k=1,48
          if (mnproc.eq.1) then
            read (nfu) tmp2d
          endif
          call xcaput(tmp2d,tflxap(1-nbdy,1-nbdy,k),1)
        enddo
        if (mnproc.eq.1) then
          close (unit=nfu)
        endif
        call xctilr(tflxap, 1,48, nbdy,nbdy, halo_ps)
      endif
      if (apsflx) then
        if (mnproc.eq.1) then
          open (unit=nfu,file=path(1:path_len)//'sflxdi.uf',
     .          form='unformatted')
          read (nfu) i,j
          if (i.ne.itdm.or.j.ne.jtdm) then
            write (lp,*) 'wrong dimensions in sflxdi.uf'
            call xchalt('(idarlx)')
                   stop '(idarlx)'
          endif
        endif
        do k=1,48
          if (mnproc.eq.1) then
            read (nfu) tmp2d
          endif
          call xcaput(tmp2d,sflxap(1-nbdy,1-nbdy,k),1)
        enddo
        if (mnproc.eq.1) then
          close (unit=nfu)
        endif
        call xctilr(sflxap, 1,48, nbdy,nbdy, halo_ps)
      endif
c
      if (ditflx.or.disflx) then
        do k=1,48
          nflxdi(k)=0
        enddo
        if (ditflx) then
c$OMP PARALLEL DO
          do k=1,48
            do j=1-nbdy,jj+nbdy
              do i=1-nbdy,ii+nbdy
                tflxdi(i,j,k)=0.
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        if (disflx) then
c$OMP PARALLEL DO
          do k=1,48
            do j=1-nbdy,jj+nbdy
              do i=1-nbdy,ii+nbdy
                sflxdi(i,j,k)=0.
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
      endif
c
      return
      end
