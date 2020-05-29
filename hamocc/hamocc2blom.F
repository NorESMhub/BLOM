c Copyright (C) 2020  J. Schwinger, M. Bentsen
c
c This file is part of BLOM/iHAMOCC.
c
c BLOM is free software: you can redistribute it and/or modify it under the
c terms of the GNU Lesser General Public License as published by the Free 
c Software Foundation, either version 3 of the License, or (at your option) 
c any later version. 
c
c BLOM is distributed in the hope that it will be useful, but WITHOUT ANY 
c WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
c FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
c more details. 
c
c You should have received a copy of the GNU Lesser General Public License 
c along with BLOM. If not, see https://www.gnu.org/licenses/.


      subroutine hamocc2blom(m,n,mm,nn)
c***********************************************************************
c
c**** *SUBROUTINE hamocc2blom* - Interface between BLOM and HAMOCC.
c
c     J. Schwinger      *GFI, UiB        2014-05-21 initial version
c      - 
c
c     Modified
c     --------
c     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
c     - added sediment bypass preprocessor option
c
c     M. Bentsen,       *NORCE, Bergen*          2020-05-03
c     - changed ocean model from MICOM to BLOM
c  
c     Purpose
c     -------
c      Pass flux and tracer fields back from HAMOCC to BLOM.
c      The local HAMOCC arrays are copied back in the appropriate
c      time-level of the tracer field. Note that also sediment fields 
c      are copied back, since a two time-level copy of sediment fields
c      is kept outside HAMOCC. For the sediment fields the same time-
c      smothing as for the tracer field (i.e. analog to tmsmt2.F) is 
c      performed to avoid a seperation of the two time levels.
c
c***********************************************************************
c
      use mod_xc
      use mo_carbch,     only: ocetra,atm
      use mo_param1_bgc, only: ks,nsedtra,npowtra,natm
      use mo_sedmnt,     only: sedlay,powtra,sedhpl,burial
      use mo_intfcblom,  only: bgc_rho,sedlay2,powtra2,burial2,atm2
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_forc.h"
c
      integer       :: m,n,mm,nn
c
      integer       :: i,j,k,l,nns,mms,kn,km
c
c
c --- ------------------------------------------------------------------
c --- pass tracer fields to ocean model; convert kmol/m^3 -> mol/kg
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
      do k=1,kk
      kn=k+nn
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        trc(i,j,kn,itrbgc:itrbgc+ntrbgc-1)=ocetra(i,j,k,:)/bgc_rho(i,j,k)
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c
c
c --- ------------------------------------------------------------------
c --- apply time smoothing for sediment fields and pass them back
c --- ------------------------------------------------------------------
c
#ifndef sedbypass
      nns=(n-1)*ks
      mms=(m-1)*ks

c$OMP PARALLEL DO PRIVATE(km,kn)
      do j=1,jj
      do k=1,ks
      km=k+mms
      kn=k+nns
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))          ! time smoothing (analog to tmsmt2.F)
        sedlay2(i,j,km,:) = wts1*sedlay2(i,j,km,:)   ! mid timelevel
     .                     +wts2*sedlay2(i,j,kn,:)   ! old timelevel
     .                     +wts2*sedlay(i,j,k,:)     ! new timelevel
        powtra2(i,j,km,:) = wts1*powtra2(i,j,km,:)
     .                     +wts2*powtra2(i,j,kn,:)
     .                     +wts2*powtra(i,j,k,:)
        burial2(i,j,m,:)  = wts1*burial2(i,j,m,:)
     .                     +wts2*burial2(i,j,n,:)
     .                     +wts2*burial(i,j,:)
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
c    
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
      do k=1,ks
      kn=k+nns
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        sedlay2(i,j,kn,:) = sedlay(i,j,k,:)  ! new time level replaces old time level here
        powtra2(i,j,kn,:) = powtra(i,j,k,:)
        burial2(i,j,n,:)  = burial(i,j,:) 
      enddo
      enddo
      enddo
      enddo
c$OMP END PARALLEL DO
#endif
c
c --- ------------------------------------------------------------------
c --- apply time smoothing for atmosphere fields if required
c --- ------------------------------------------------------------------
c
#if defined(BOXATM)
c$OMP PARALLEL DO
      do j=1,jj
      do i=1,ii                            ! time smoothing (analog to tmsmt2.F)
        atm2(i,j,m,:) = wts1*atm2(i,j,m,:) ! mid timelevel
     .                 +wts2*atm2(i,j,n,:) ! old timelevel
     .                 +wts2*atm(i,j,:)    ! new timelevel
      enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
      do i=1,ii
        atm2(i,j,n,:) = atm(i,j,:)  ! new time level replaces old time level here
      enddo
      enddo
c$OMP END PARALLEL DO
#endif

      return
      end
