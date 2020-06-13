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
          open (unit=nfu,file='tflxdi.uf',form='unformatted')
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
          open (unit=nfu,file='sflxdi.uf',form='unformatted')
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
c$OMP PARALLEL DO PRIVATE(k,i)
          do j=1-nbdy,jj+nbdy
            do k=1,48
              do i=1-nbdy,ii+nbdy
                tflxdi(i,j,k)=0.
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        if (disflx) then
c$OMP PARALLEL DO PRIVATE(k,i)
          do j=1-nbdy,jj+nbdy
            do k=1,48
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
