! ------------------------------------------------------------------------------
! Copyright (C) 2005-2018 Mats Bentsen, Alok Kumar Gupta
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

      subroutine wdiflx
c
c --- Write accumulated diagnosed heat and salt fluxes
c
      use mod_xc
      use mod_nctools
      use mod_dia, only : iotype
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
c
      character fname*256
      integer i,j,k
c
      if (ditflx) then
c
        write (fname,'(2a,i4.4,a)')
     .    runid(1:runid_len),'_tflxdi_',nyear-1,'.nc'
        if (mnproc.eq.1) then
          write (lp,'(2a)')
     .      'Writing diagnostic heat flux to ',trim(fname)
        endif
        call ncfopn(fname,'w','c',1,iotype)
c
#ifdef COMPDIA
        call ncdimc('pcomp',ip,0)
#else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
#endif 
        call ncdims('week',48)
c
        call ncputr('time',time)
        call ncputi('nflxdi',nflxdi)
c
#ifdef COMPDIA
        call ncdefvar('tflxdi','pcomp week',ndouble,8)
        call ncedef 
        call nccomp('tflxdi','pcomp week',tflxdi,ip,1.,0.,8)
#else
        call ncdefvar('tflxdi','x y week',ndouble,8)
        call ncedef
        call ncwrtr('tflxdi','x y week',tflxdi,ip,1,1.,0.,8)
#endif 
c       
        call ncfcls
c
c$OMP PARALLEL DO PRIVATE(j,i)
        do k=1,48
          nflxdi(k)=0
          do j=1,jj
            do i=1,ii
              tflxdi(i,j,k)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      endif
c
      if (disflx) then
c
        write (fname,'(2a,i4.4,a)')
     .    runid(1:runid_len),'_sflxdi_',nyear-1,'.nc'
        if (mnproc.eq.1) then
          write (lp,'(2a)')
     .      'Writing diagnostic salt flux to ',trim(fname)
        endif
        call ncfopn(fname,'w','c',1,iotype)
c
#ifdef COMPDIA
        call ncdimc('pcomp',ip,0)
#else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
#endif 
        call ncdims('week',48)
c
        call ncputr('time',time)
        call ncputi('nflxdi',nflxdi)
c
#ifdef COMPDIA
        call nccomp('sflxdi','pcomp week',sflxdi,ip,1.,0.,8)
#else
        call ncwrtr('sflxdi','x y week',sflxdi,ip,1,1.,0.,8)
#endif 
c       
        call ncfcls
c
c$OMP PARALLEL DO PRIVATE(j,i)
        do k=1,48
          nflxdi(k)=0
          do j=1,jj
            do i=1,ii
              sflxdi(i,j,k)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      endif
c
      return
      end
