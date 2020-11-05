! ------------------------------------------------------------------------------
! Copyright (C) 2015-2019 Mats Bentsen
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

      subroutine cmnfld_nslope(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Estimate slope of local neutral surface.
c --- ------------------------------------------------------------------
c
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real rho0,pm,rho_x,phi_x,bfsqm,rho_y,phi_y
      integer i,j,k,l,km,kintr,kmax,knnsl
c
      real p_alpha,rho
      external p_alpha,rho
c
c --- ------------------------------------------------------------------
c --- Compute geopotential at layer interfaces.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
      do j=-1,jj+2
        do k=kk,1,-1
          km=k+mm
          do l=1,isp(j)
          do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
            if (dp(i,j,km).lt.epsil) then
              phi(i,j,k)=phi(i,j,k+1)
            else
              phi(i,j,k)=phi(i,j,k+1)
     .                  -p_alpha(p(i,j,k+1),p(i,j,k),
     .                           temp(i,j,km),saln(i,j,km))
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute slope vector of local neutral surfaces and also slope
c --- vector times Brunt-Vaisala frequency (optionally used in the
c --- computation of eddy growth rate). The latter is not computed when
c --- the gradient of the geopotential is expected to be influenced by
c --- the gradient of the bathymetry and in this case values are
c --- extrapolated from above.
c --- ------------------------------------------------------------------
c
      rho0=1./alpha0
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  l,i,k,kmax,km,kintr,knnsl,pm,rho_x,phi_x,bfsqm)
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
c
c --- --- Set the x-component of the slope vector to zero initially.
          do k=1,kk
            nslpx(i,j,k)=0.
            nnslpx(i,j,k)=0.
          enddo
c
          if     (kfpla(i-1,j,m).le.kk.or.kfpla(i,j,m).le.kk) then
c
c --- ----- Index of last layer containing mass at either of the scalar
c --- ----- points adjacent to the velocity point.
            kmax=1
            do k=3,kk
              km=k+mm
              if (dp(i-1,j,km).gt.epsil.or.dp(i,j,km).gt.epsil) kmax=k
            enddo
c
c --- ----- The first interior interface where the x-component of the
c --- ----- slope vector is estimated is at index kintr+1.
            kintr=max(kfpla(i-1,j,m),kfpla(i,j,m))
c
c --- ----- Index of last interface where slope vector times
c --- ----- Brunt-Vaisala frequency is computed.
            knnsl=2
c
c --- ----- Compute the x-component of the slope vector at the mixed
c --- ----- layer base.
            pm=.5*(p(i-1,j,3)+p(i,j,3))
            rho_x=rho(pm,temp(i  ,j,2+mm),saln(i  ,j,2+mm))
     .           -rho(pm,temp(i-1,j,2+mm),saln(i-1,j,2+mm))
            phi_x=phi(i,j,3)-phi(i-1,j,3)
            bfsqm=.5*(bfsqf(i-1,j,3)+bfsqf(i,j,3))
            nslpx(i,j,3)=(g*rho_x/(rho0*bfsqm)+phi_x/g)*scuxi(i,j)
            if (phi(i  ,j,3).gt.phi(i-1,j,kk+1).and.
     .          phi(i-1,j,3).gt.phi(i  ,j,kk+1)) then
              nnslpx(i,j,3)=sqrt(bfsqm)*nslpx(i,j,3)
              knnsl=3
            endif
c
c --- ----- Compute the x-component of the slope vector at interior
c --- ----- interfaces.
            do k=kintr+1,kmax
              km=k+mm
              pm=.5*(p(i-1,j,k)+p(i,j,k))
              rho_x=.5*(rho(pm,temp(i  ,j,km-1),saln(i  ,j,km-1))
     .                 -rho(pm,temp(i-1,j,km-1),saln(i-1,j,km-1))
     .                 +rho(pm,temp(i  ,j,km  ),saln(i  ,j,km  ))
     .                 -rho(pm,temp(i-1,j,km  ),saln(i-1,j,km  )))
              phi_x=phi(i,j,k)-phi(i-1,j,k)
              bfsqm=.5*(bfsqf(i-1,j,k)+bfsqf(i,j,k))
              nslpx(i,j,k)=(g*rho_x/(rho0*bfsqm)+phi_x/g)*scuxi(i,j)
              if (phi(i  ,j,k).gt.phi(i-1,j,kk+1).and.
     .            phi(i-1,j,k).gt.phi(i  ,j,kk+1)) then
                nnslpx(i,j,k)=sqrt(bfsqm)*nslpx(i,j,k)
                knnsl=k
              endif
            enddo
            do k=knnsl+1,kmax
              nnslpx(i,j,k)=nnslpx(i,j,knnsl)
            enddo
            if (kintr.lt.kmax) then
              do k=4,kintr
                nslpx(i,j,k)=nslpx(i,j,kintr+1)
                nnslpx(i,j,k)=nnslpx(i,j,kintr+1)
              enddo
            else
              do k=4,kmax
                nslpx(i,j,k)=nslpx(i,j,3)
                nnslpx(i,j,k)=nnslpx(i,j,3)
              enddo
            endif
c
          endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  l,i,k,kmax,km,kintr,knnsl,pm,rho_y,phi_y,bfsqm)
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
c
c --- --- Set the y-component of the slope vector to zero initially.
          do k=1,kk
            nslpy(i,j,k)=0.
            nnslpy(i,j,k)=0.
          enddo
c
          if     (kfpla(i,j-1,m).le.kk.or.kfpla(i,j,m).le.kk) then
c
c --- ----- Index of last layer containing mass at either of the scalar
c --- ----- points adjacent to the velocity point.
            kmax=1
            do k=3,kk
              km=k+mm
              if (dp(i,j-1,km).gt.epsil.or.dp(i,j,km).gt.epsil) kmax=k
            enddo
c
c --- ----- The first interior interface where the y-component of the
c --- ----- slope vector is estimated is at index kintr+1.
            kintr=max(kfpla(i,j-1,m),kfpla(i,j,m))
c
c --- ----- Index of last interface where slope vector times
c --- ----- Brunt-Vaisala frequency is computed.
            knnsl=2
c
c --- ----- Compute the y-component of the slope vector at the mixed
c --- ----- layer base.
            pm=.5*(p(i,j-1,3)+p(i,j,3))
            rho_y=rho(pm,temp(i,j  ,2+mm),saln(i,j  ,2+mm))
     .           -rho(pm,temp(i,j-1,2+mm),saln(i,j-1,2+mm))
            phi_y=phi(i,j,3)-phi(i,j-1,3)
            bfsqm=.5*(bfsqf(i,j-1,3)+bfsqf(i,j,3))
            nslpy(i,j,3)=(g*rho_y/(rho0*bfsqm)+phi_y/g)*scvyi(i,j)
            if (phi(i,j  ,3).gt.phi(i,j-1,kk+1).and.
     .          phi(i,j-1,3).gt.phi(i,j  ,kk+1)) then
              nnslpy(i,j,3)=sqrt(bfsqm)*nslpy(i,j,3)
              knnsl=3
            endif
c
c --- ----- Compute the y-component of the slope vector at interior
c --- ----- interfaces.
            do k=kintr+1,kmax
              km=k+mm
              pm=.5*(p(i,j-1,k)+p(i,j,k))
              rho_y=.5*(rho(pm,temp(i,j  ,km-1),saln(i,j  ,km-1))
     .                 -rho(pm,temp(i,j-1,km-1),saln(i,j-1,km-1))
     .                 +rho(pm,temp(i,j  ,km  ),saln(i,j  ,km  ))
     .                 -rho(pm,temp(i,j-1,km  ),saln(i,j-1,km  )))
              phi_y=phi(i,j,k)-phi(i,j-1,k)
              bfsqm=.5*(bfsqf(i,j-1,k)+bfsqf(i,j,k))
              nslpy(i,j,k)=(g*rho_y/(rho0*bfsqm)+phi_y/g)*scvyi(i,j)
              if (phi(i,j  ,k).gt.phi(i,j-1,kk+1).and.
     .            phi(i,j-1,k).gt.phi(i,j  ,kk+1)) then
                nnslpy(i,j,k)=sqrt(bfsqm)*nslpy(i,j,k)
                knnsl=k
              endif
            enddo
            do k=knnsl+1,kmax
              nnslpy(i,j,k)=nnslpy(i,j,knnsl)
            enddo
            if (kintr.lt.kmax) then
              do k=4,kintr
                nslpy(i,j,k)=nslpy(i,j,kintr+1)
                nnslpy(i,j,k)=nnslpy(i,j,kintr+1)
              enddo
            else
              do k=4,kmax
                nslpy(i,j,k)=nslpy(i,j,3)
                nnslpy(i,j,k)=nnslpy(i,j,3)
              enddo
            endif
c
          endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'cmnfld_nslope:'
        endif
        call chksummsk(nslpx,iu,kk,'nslpx')
        call chksummsk(nslpy,iv,kk,'nslpy')
        call chksummsk(nnslpx,iu,kk,'nnslpx')
        call chksummsk(nnslpy,iv,kk,'nnslpy')
      endif
c
      return
      end
