! ------------------------------------------------------------------------------
! Copyright (C) 2012-2018 Mats Bentsen

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

      subroutine fwbbal(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Balance the freshwater budget by computing a correcting factor
c --- to be applied to precipitation and runoff. The correction factor
c --- is based on the fresh water budget for the previous year.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
c
      real totei,totpr
      integer i,j,l
c
c --- Accumulate two groups of fresh water fluxes. One is evaporation
c --- and sea-ice melting/freezing and the other is precipitation and
c --- runoff. The fresh water fluxes are weighted with the time step in
c --- case it varies during the accumulation period.
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          eiacc(i,j)=eiacc(i,j)
     .              +(eva(i,j)+fmltfz(i,j))*baclin
          pracc(i,j)=pracc(i,j)
     .              +(lip(i,j)+sop(i,j)+rnf(i,j)+rfi(i,j))*baclin
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute new correction factor at the end of a year and reset
c --- accumulation arrays
      if (nday_of_year.eq.1.and.mod(nstep,nstep_in_day).eq.0) then
c
c --- - Weight the accumulated fluxes with grid cell area and do global
c --- - sums, but only including grid cells connected to the world ocean.
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            eiacc(i,j)=eiacc(i,j)*scp2(i,j)
            pracc(i,j)=pracc(i,j)*scp2(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(totei,eiacc,ipwocn)
        call xcsum(totpr,pracc,ipwocn)
c
c --- - Update correction factor
        prfac=-prfac*totei/totpr
        if (mnproc.eq.1) then
          write (lp,*)
     .      'new correction factor for precipitation/runoff:',prfac
          call flush(lp)
        endif
c
c --- - Reset accumulation arrays
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            eiacc(i,j)=0.
            pracc(i,j)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'fwbbal:'
        endif
        call chksummsk(eiacc,ip,1,'eiacc')
        call chksummsk(pracc,ip,1,'pracc')
      endif
c
      return
      end
