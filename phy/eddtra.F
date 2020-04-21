! ------------------------------------------------------------------------------
! Copyright (C) 2015-2018 Mats Bentsen
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

      subroutine eddtra(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Compute eddy-induced transport.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      integer i,j,k,l,km
c
c --- ------------------------------------------------------------------
c --- Compute eddy-induced transport of mass.
c --- ------------------------------------------------------------------
c
      if     (eitmth.eq.'intdif') then
        call eddtra_intdif(m,n,mm,nn,k1m,k1n)
      elseif (eitmth.eq.'gm') then
        call eddtra_gm(m,n,mm,nn,k1m,k1n)
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' eitmth=',trim(eitmth),' is unsupported!'
        endif
        call xcstop('(eddtra)')
               stop '(eddtra)'
      endif
c
c --- ------------------------------------------------------------------
c --- Diagnose eddy-induced transport components of heat and salt.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
      do j=1,jj
        do k=1,kk
          km=k+mm
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            utfltd(i,j,km)=.5*umfltd(i,j,km)
     .                     *(temp(i-1,j,km)+temp(i,j,km))
            usfltd(i,j,km)=.5*umfltd(i,j,km)
     .                     *(saln(i-1,j,km)+saln(i,j,km))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtfltd(i,j,km)=.5*vmfltd(i,j,km)
     .                     *(temp(i,j-1,km)+temp(i,j,km))
            vsfltd(i,j,km)=.5*vmfltd(i,j,km)
     .                     *(saln(i,j-1,km)+saln(i,j,km))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'eddtra:'
        endif
        call chksummsk(umfltd(1-nbdy,1-nbdy,k1m),iu,kk,'umfltd')
        call chksummsk(vmfltd(1-nbdy,1-nbdy,k1m),iv,kk,'vmfltd')
        call chksummsk(utfltd(1-nbdy,1-nbdy,k1m),iu,kk,'utfltd')
        call chksummsk(vtfltd(1-nbdy,1-nbdy,k1m),iv,kk,'vtfltd')
        call chksummsk(usfltd(1-nbdy,1-nbdy,k1m),iu,kk,'usfltd')
        call chksummsk(vsfltd(1-nbdy,1-nbdy,k1m),iv,kk,'vsfltd')
      endif
c
      return
      end
