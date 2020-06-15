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

      subroutine cmnfld_bfsqf(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Compute buoyancy frequency squared (BFSQ) on layer interfaces and
c --- representative of the layer itself. Also compute a filtered BFSQ
c --- on interfaces.
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
c --- Parameters:
c ---   sls0   - minimum smoothing length scale in the computation of
c ---            filtered BFSQ [g/cm/s**2].
c ---   slsmfq - factor to be multiplied with the mixed layer depth to
c ---            find the smoothing length scale at the base of the
c ---            mixed layer in the computation of filtered BFSQ [].
c ---   slsels - factor to be multiplied with the mixed layer depth to
c ---            find the e-folding length scale of the smoothing length
c ---            scale in the computation of filtered BFSQ [].
c ---   bfsqmn - minimum value of BFSQ used in the computation of
c ---            neutral slope [s-2].
c
      real sls0,slsmfq,slsels,bfsqmn
      parameter (sls0=10.*98060.,slsmfq=2.,slsels=2.,bfsqmn=1.e-7)
c
      real, dimension(kdm) :: delp,bfsq,sls2,atd,btd,ctd,rtd,gam
      real pml,q,pup,tup,sup,plo,tlo,slo,bei
      integer i,j,k,l,km,kfpl
c
      real rho
      external rho
c
c --- ------------------------------------------------------------------
c --- The BFSQ is estimated locally at layer interfaces. The filtered
c --- BFSQ is smoothed in the vertical direction by solving a diffusion
c --- equation. At the mixed layer base the diffusion length scale is
c --- set to the maximum of sls0 and mixed layer depth (MLD) times
c --- slsmfq. Below the mixed layer, the diffusion length scale
c --- approaches sls0 with an e-folding length scale of MLD times
c --- slsels.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  l,i,kfpl,k,pml,delp,bfsq,q,sls2,pup,tup,sup,km,plo,tlo,slo,ctd,
c$OMP+  btd,rtd,atd,bei,gam)
      do j=-1,jj+2
        do l=1,isp(j)
        do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
c
c --- --- Compute BFSQ in the mixed layer.
          bfsqi(i,j,1)=
     .      .5*g*g*(rho(p(i,j,2),temp(i,j,2+mm),saln(i,j,2+mm))
     .             -rho(p(i,j,2),temp(i,j,1+mm),saln(i,j,1+mm)))
     .      /(dp(i,j,1+mm)+dp(i,j,2+mm))
          bfsqi(i,j,2)=bfsqi(i,j,1)
          bfsql(i,j,1)=bfsqi(i,j,1)
          bfsql(i,j,2)=bfsqi(i,j,1)
c
          kfpl=kfpla(i,j,m)
c
          if (kfpl.gt.kk) then
c
c --- ----- If the mixed layer extends to the bottom, propagate the
c --- ----- interface and layer BFSQ of the mixed layer downwards while
c --- ----- the filtered BFSQ is set to a minimum value.
            do k=3,kk
              bfsqi(i,j,k)=bfsqi(i,j,1)
              bfsql(i,j,k)=bfsqi(i,j,1)
            enddo
            do k=1,kk
              bfsqf(i,j,k)=bfsqmn
            enddo
c
          else
c
c --- ----- At layer interfaces, compute BFSQ and length scale for the
c --- ----- subsequent smoothing.
            pml=max(.5*(p(i,j,3)+p(i,j,1)),
     .              .5*(3.*p(i,j,3)-p(i,j,kfpl+1)))
            delp(kfpl-1)=pml-p(i,j,1)
            bfsqi(i,j,kfpl-1)=bfsqi(i,j,2)
            bfsq(kfpl-1)=bfsqmn
            q=max(sls0,delp(kfpl-1)*slsmfq)
            sls2(kfpl-1)=q*q
            pup=pml
            tup=temp(i,j,2+mm)
            sup=saln(i,j,2+mm)
            do k=kfpl,kk
              km=k+mm
              if (p(i,j,kk+1)-p(i,j,k).lt.epsil) then
                delp(k)=onemm
                bfsqi(i,j,k)=bfsqi(i,j,k-1)
                bfsq(k)=bfsqmn
                q=exp(-(p(i,j,kk+1)-pml)/(slsels*delp(kfpl-1)))
                q=max(sls0,delp(kfpl-1)*slsmfq*q+sls0*(1.-q))
                sls2(k)=q*q
              else
                if (p(i,j,kk+1)-p(i,j,k+1).lt.epsil) then
                  plo=p(i,j,kk+1)
                else
                  plo=.5*(p(i,j,k)+p(i,j,k+1))
                endif
                tlo=temp(i,j,km)
                slo=saln(i,j,km)
                delp(k)=max(onemm,plo-pup)
                bfsqi(i,j,k)=g*g*(rho(p(i,j,k),tlo,slo)
     .                           -rho(p(i,j,k),tup,sup))/delp(k)
                bfsq(k)=max(bfsqmn,bfsqi(i,j,k))
                bfsqi(i,j,k)=bfsqi(i,j,k)*delp(k)/max(onem,delp(k))
                if (p(i,j,kk+1)-p(i,j,k).lt.onem) then
                  bfsqi(i,j,k)=bfsqi(i,j,k-1)
                endif
                q=exp(-(p(i,j,k)-pml)/(slsels*delp(kfpl-1)))
                q=max(sls0,delp(kfpl-1)*slsmfq*q+sls0*(1.-q))
                sls2(k)=q*q
                pup=plo
                tup=tlo
                sup=slo
              endif
            enddo
c
c --- ----- Compute the layer BFSQ as the arithmetic mean of the layer
c --- ----- interface BFSQ.
            do k=kfpl,kk-1
              bfsql(i,j,k)=.5*(bfsqi(i,j,k)+bfsqi(i,j,k+1))
            enddo
            bfsql(i,j,kk)=bfsqi(i,j,kk)
            do k=3,kfpl-1
              bfsqi(i,j,k)=bfsqi(i,j,kfpl)
              bfsql(i,j,k)=bfsql(i,j,kfpl)
            enddo
c
c --- ----- For the filtered BFSQ, compute the coefficients for the
c --- ----- tridiagonal set of equations arising from the implicit
c --- ----- backward discretization.
            k=kfpl-1
            ctd(k)=-2.*sls2(k  )/(delp(k)*(delp(k  )+delp(k+1)))
            btd(k)=1.-ctd(k)
            rtd(k)=bfsq(k)
            do k=kfpl,kk-1
              atd(k)=-2.*sls2(k-1)/(delp(k)*(delp(k-1)+delp(k  )))
              ctd(k)=-2.*sls2(k  )/(delp(k)*(delp(k  )+delp(k+1)))
              btd(k)=1.-atd(k)-ctd(k)
              rtd(k)=bfsq(k)
            enddo
            k=kk
            atd(k)=-2.*sls2(k-1)/(delp(k)*(delp(k-1)+delp(k  )))
            btd(k)=1.-atd(k)
            rtd(k)=bfsq(k)
c
c --- ----- Solve the tridiagonal set of equations.
            bei=1./btd(kfpl-1)
            bfsqf(i,j,kfpl-1)=rtd(kfpl-1)*bei
            do k=kfpl,kk
              gam(k)=ctd(k-1)*bei
              bei=1./(btd(k)-atd(k)*gam(k))
              bfsqf(i,j,k)=(rtd(k)-atd(k)*bfsqf(i,j,k-1))*bei
            enddo
            do k=kk-1,kfpl-1,-1
              bfsqf(i,j,k)=bfsqf(i,j,k)-gam(k+1)*bfsqf(i,j,k+1)
            enddo
            do k=1,kfpl-2
              bfsqf(i,j,k)=bfsqf(i,j,kfpl-1)
            enddo
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
          write (lp,*) 'cmnfld_bfsqf:'
        endif
        call chksummsk(bfsqi,ip,kk,'bfsqi')
        call chksummsk(bfsql,ip,kk,'bfsql')
        call chksummsk(bfsqf,ip,kk,'bfsqf')
      endif
c
      return
      end
