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

      subroutine eddtra_gm(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Estimate eddy-induced transport following the Gent-McWilliams
c --- parameterization.
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
c ---   ffac   - fraction of the mass of a grid cell a mass flux is
c ---            allowed to deplete [].
c ---   fface  - (1-epsilon)*ffac [].
c
      real ffac,fface
      parameter (ffac=.0625,fface=.99*ffac)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ptu,ptv
      real, dimension(kdm+1) :: upsilon,mfl
      real, dimension(kdm) :: dlm,dlp
      real rho0,q,et2mf,kappa,fhi,flo
      integer i,j,k,l,km,kn,kintr,kmax,kmin,niter,kdir
      logical changed
c
      real rho
      external rho
c
      rho0=1./alpha0
c
      call xctilr(difint, 1,kk, 2,2, halo_ps)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
c
c --- ------------------------------------------------------------------
c --- Compute top pressure at velocity points.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
          ptu(i,j)=max(p(i-1,j,1),p(i,j,1))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
          ptv(i,j)=max(p(i,j-1,1),p(i,j,1))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute u-component of eddy-induced mass fluxes.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  l,i,k,km,et2mf,kmax,kn,kintr,kappa,upsilon,
c$OMP+  kmin,mfl,dlm,dlp,fhi,flo,changed,niter,kdir,q)
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
c
c --- ------------------------------------------------------------------
c --- --- Set eddy-induced mass fluxes to zero initially.
c --- ------------------------------------------------------------------
c
          do k=1,kk
            km=k+mm
            umfltd(i,j,km)=0.
          enddo
c
c --- --- Eddy transport to mass flux conversion factor.
          et2mf=-g*rho0*delt1*scuy(i,j)
c
c --- --- Index of last layer containing mass at either of the scalar
c --- --- points adjacent to the velocity point.
          kmax=1
          do k=3,kk
            kn=k+nn
            if (dp(i-1,j,kn).gt.epsil.or.dp(i,j,kn).gt.epsil) kmax=k
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Proceed with mass flux computation if at least one of the
c --- --- adjacent scalar points to the velocity point has a mass
c --- --- containing interior layer. Mass fluxes will be assigned at
c --- --- layer interface corresponding to the eddy induced transport.
c --- --- The final layer mass flux will be the lower minus the upper
c --- --- interface flux.  The mass fluxes are limited to keep
c --- --- interfaces within the water column. There are 3 cases to
c --- --- consider:
c --- ---   Case 1: The mixed layer extends to the bottom at both
c --- ---           adjacent scalar points to the velocity point
c --- ---   Case 2: The mixed layer extends to the bottom at scalar
c --- ---           point (i  ,j).
c --- ---   Case 3: The mixed layer extends to the bottom at scalar
c --- ---           point (i-1,j).
c --- ---   Case 4: The mixed layer does not reach the bottom at neither
c --- ---           of the scalar points adjacent to the velocity point.
c --- ------------------------------------------------------------------
c
          if     (kfpla(i-1,j,n).gt.kk.and.kfpla(i,j,n).gt.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Case 1:
c --- ------------------------------------------------------------------
c
c --- ----- Keep the initial zero mass fluxes for this column.
            cycle
c
          elseif (kfpla(i-1,j,n).le.kk.and.kfpla(i,j,n).gt.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Case 2:
c --- ------------------------------------------------------------------
c
c --- ----- Find the index of the first layer at (i-1,j) that is
c --- ----- hydrostatically stable at the mixed layer base at (i  ,j).
            km=2+nn
            kintr=kfpla(i-1,j,n)
            kn=kintr+nn
            do while
     .        (rho(p(i  ,j,3),temp(i-1,j,kn),saln(i-1,j,kn)).lt.
     .         rho(p(i  ,j,3),temp(i  ,j,km),saln(i  ,j,km)).or.
     .         dp(i-1,j,kn).lt.epsil)
              kintr=kintr+1
              if (kintr.eq.kmax+1) exit
              kn=kintr+nn
            enddo
c
c --- ----- If a physical layer cannot be found, keep the initial zero
c --- ----- mass fluxes for this column.
            if (kintr.eq.kmax+1) cycle
c
c --- ----- Compute the eddy induced transport (upsilon) at the mixed
c --- ----- layer base.
            kappa=.5*(difint(i-1,j,2)+difint(i,j,2))
            upsilon(3)=-kappa*nslpx(i,j,3)
c
c --- ----- If the eddy-induced transport at the base of the mixed layer
c --- ----- would cause a negative mass flux below the mixed layer, keep
c --- ----- the initial zero mass fluxes for this column.
            if (upsilon(3).le.0.) cycle
c
c --- ----- Assign interface mass fluxes.
            kmin=kintr-1
            mfl(kmin)=0.
            mfl(kintr)=et2mf*upsilon(3)
            do k=kintr+1,kmax+1
              mfl(k)=0.
            enddo
c
          elseif (kfpla(i-1,j,n).gt.kk.and.kfpla(i,j,n).le.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Case 3:
c --- ------------------------------------------------------------------
c
c --- ----- Find the index of the first layer at (i  ,j) that is
c --- ----- hydrostatically stable at the mixed layer base at (i-1,j).
            km=2+nn
            kintr=kfpla(i  ,j,n)
            kn=kintr+nn
            do while
     .        (rho(p(i-1,j,3),temp(i  ,j,kn),saln(i  ,j,kn)).lt.
     .         rho(p(i-1,j,3),temp(i-1,j,km),saln(i-1,j,km)).or.
     .         dp(i  ,j,kn).lt.epsil)
              kintr=kintr+1
              if (kintr.eq.kmax+1) exit
              kn=kintr+nn
            enddo
c
c --- ----- If a physical layer cannot be found, keep the initial zero
c --- ----- mass fluxes for this column.
            if (kintr.eq.kmax+1) cycle
c
c --- ----- Compute the eddy induced transport (upsilon) at the mixed
c --- ----- layer base.
            kappa=.5*(difint(i-1,j,2)+difint(i,j,2))
            upsilon(3)=-kappa*nslpx(i,j,3)
c
c --- ----- If the eddy-induced transport at the base of the mixed layer
c --- ----- would cause a positive mass flux below the mixed layer, keep
c --- ----- the initial zero mass fluxes for this column.
            if (upsilon(3).ge.0.) cycle
c
c --- ----- Assign interface mass fluxes.
            kmin=kintr-1
            mfl(kmin)=0.
            mfl(kintr)=et2mf*upsilon(3)
            do k=kintr+1,kmax+1
              mfl(k)=0.
            enddo
c
          else
c
c --- ------------------------------------------------------------------
c --- ----- Case 4:
c --- ------------------------------------------------------------------
c
c --- ----- The first interior interface where the eddy induced
c --- ----- transport is estimated is at index kintr+1.
            kintr=max(kfpla(i-1,j,n),kfpla(i,j,n))
c
c --- ----- Compute the eddy induced transport (upsilon) at the mixed
c --- ----- layer base.
            kappa=.5*(difint(i-1,j,2)+difint(i,j,2))
            upsilon(3)=-kappa*nslpx(i,j,3)
c
c --- ----- Compute the eddy induced transport at interior interfaces.
            do k=kintr+1,kmax
              kn=k+nn
              kappa=.25*(difint(i-1,j,k-1)+difint(i  ,j,k-1)
     .                  +difint(i-1,j,k  )+difint(i  ,j,k  ))
              upsilon(k)=-kappa*nslpx(i,j,k)
            enddo
            upsilon(kmax+1)=0.
c
c --- ----- If the layer kintr-1 is a physical layer at either of the
c --- ----- adjacent scalar points to the velocity point, then apply
c --- ----- an upper interface mass flux corresponding to the eddy
c --- ----- induced transport at the mixed layer base and a lower
c --- ----- interface mass flux corresponding to the eddy induced
c --- ----- transport at the kintr+1 interface if this would lead to a
c --- ----- hydrostatically stable layer arrangement.
            km=2+nn
            kn=kintr-1+nn
            if ((kfpla(i-1,j,n).lt.kintr.and.
     .           upsilon(3)-upsilon(kintr+1).gt.0..and.
     .           rho(p(i  ,j,3),temp(i-1,j,kn),saln(i-1,j,kn)).gt.
     .           rho(p(i  ,j,3),temp(i  ,j,km),saln(i  ,j,km))).or.
     .          (kfpla(i  ,j,n).lt.kintr.and.
     .           upsilon(3)-upsilon(kintr+1).lt.0..and.
     .           rho(p(i-1,j,3),temp(i  ,j,kn),saln(i  ,j,kn)).gt.
     .           rho(p(i-1,j,3),temp(i-1,j,km),saln(i-1,j,km)))) then
              kintr=kintr-1
              upsilon(kintr+1)=upsilon(kintr+2)
            endif
c
c --- ----- Assign interface mass fluxes.
            kmin=kintr-1
            mfl(kmin)=0.
            mfl(kintr)=et2mf*upsilon(3)
            do k=kintr+1,kmax
              mfl(k)=et2mf*upsilon(k)
            enddo
            mfl(kmax+1)=0.
c
          endif
c
c --- ------------------------------------------------------------------
c --- --- Ensure that mass fluxes do not create negative layer
c --- --- thicknesses.
c --- ------------------------------------------------------------------
c
c --- --- Compute the layer thicknesses available to be depleted by mass
c --- --- fluxes at the scalar points adjacent to the velocity point.
c --- --- These bounded layer thicknesses are consistent with the
c --- --- transport algorithm.
          dlm(kmin)=max(0.,min(p(i-1,j,3),pbu(i,j,n))
     .                    -max(p(i-1,j,1),ptu(i,j)))
          dlp(kmin)=max(0.,min(p(i  ,j,3),pbu(i,j,n))
     .                    -max(p(i  ,j,1),ptu(i,j)))
          do k=kintr,kmax
            dlm(k)=max(0.,min(p(i-1,j,k+1),pbu(i,j,n))
     .                   -max(p(i-1,j,k  ),ptu(i,j)))
            dlp(k)=max(0.,min(p(i  ,j,k+1),pbu(i,j,n))
     .                   -max(p(i  ,j,k  ),ptu(i,j)))
          enddo
c
c --- --- If excessive depletion of layers occur beneath the mixed layer
c --- --- base, try to adjust interface fluxes other that the mixed
c --- --- layer base interface flux.
          fhi= fface
     .         *max(0.,min((p(i-1,j,3)-ptu(i,j))*scp2(i-1,j),
     .                     (pbu(i,j,n)-p(i  ,j,kintr))*scp2(i  ,j)))
          flo=-fface
     .         *max(0.,min((p(i  ,j,3)-ptu(i,j))*scp2(i  ,j),
     .                     (pbu(i,j,n)-p(i-1,j,kintr))*scp2(i-1,j)))
          mfl(kmin+1)=min(fhi,max(flo,mfl(kmin+1)))
          do k=kmin+1,kmax-1
            if     (mfl(k+1)-mfl(k).gt.
     .              ffac*max(epsil,dlm(k))*scp2(i-1,j)) then
              mfl(k+1)=mfl(k  )+fface*dlm(k)*scp2(i-1,j)
            elseif (mfl(k+1)-mfl(k).lt.
     .             -ffac*max(epsil,dlp(k))*scp2(i  ,j)) then
              mfl(k+1)=mfl(k  )-fface*dlp(k)*scp2(i  ,j)
            else
              exit
            endif
          enddo
c
c --- --- Apply an iterative procedure for flux limiting by alternate
c --- --- upward and downward propagation through the layers.
c
          changed=.true.
          niter=0
          kdir=1
c
          do while (changed)
c
            niter=niter+1
            if (niter.eq.1000) then
              k=kmin
              write (lp,*)
              write (lp,'(i3,3e16.8)') 1,mfl(k+1),mfl(k),
     . (mfl(k+1)-mfl(k))
     . /(max(onemm,dpu(i,j,1+nn)+dpu(i,j,2+nn))*delt1*scuy(i,j))
              do k=kintr,kmax
                kn=k+nn
                write (lp,'(i3,3e16.8)') k,mfl(k+1),mfl(k),
     . (mfl(k+1)-mfl(k))/(max(onemm,dpu(i,j,kn))*delt1*scuy(i,j))
              enddo
              write (lp,*) 'no convergence u',i+i0,j+j0
              call xchalt('(eddtra_gm)')
                     stop '(eddtra_gm)'
            endif
c
            changed=.false.
            kdir=-kdir
c
            do k=((1-kdir)*kmax+(1+kdir)*kmin)/2,
     .           ((1-kdir)*kmin+(1+kdir)*kmax)/2,kdir
c
c --- ------- Proceed with flux limiting of this layer if the mass flux
c --- ------- difference between lower and upper interface is beyond the
c --- ------- floating point accuracy limitation.
              if (abs(mfl(k+1)-mfl(k)).gt.
     .            1.e-14*max(epsil*scu2(i,j),abs(mfl(k+1)+mfl(k)))) then
c
                if     (mfl(k+1)-mfl(k).gt.
     .                  ffac*max(epsil,dlm(k))*scp2(i-1,j)) then
c
c --- ----------- In this case, the mass fluxes are removing too much
c --- ----------- mass from the grid cell at (i-1,j,k). Limit the
c --- ----------- dominating interface flux.
                  q=fface*dlm(k)*scp2(i-1,j)
                  if (mfl(k+1).gt.-mfl(k)) then
                    if (mfl(k  ).gt.-.5*q) then
                      mfl(k+1)=mfl(k  )+q
                    else
                      mfl(k+1)= .5*q
                      mfl(k  )=-mfl(k+1)
                    endif
                  else
                    if (mfl(k+1).lt. .5*q) then
                      mfl(k  )=mfl(k+1)-q
                    else
                      mfl(k  )=-.5*q
                      mfl(k+1)=-mfl(k  )
                    endif
                  endif
                  changed=.true.
                elseif (mfl(k+1)-mfl(k).lt.
     .                 -ffac*max(epsil,dlp(k))*scp2(i  ,j)) then
c
c --- ----------- In this case, the mass fluxes are removing too much
c --- ----------- mass from the grid cell at (i  ,j,k). Limit the
c --- ----------- dominating interface flux.
                  q=fface*dlp(k)*scp2(i  ,j)
                  if (mfl(k+1).lt.-mfl(k)) then
                    if (mfl(k  ).lt. .5*q) then
                      mfl(k+1)=mfl(k  )-q
                    else
                      mfl(k+1)=-.5*q
                      mfl(k  )=-mfl(k+1)
                    endif
                  else
                    if (mfl(k+1).gt.-.5*q) then
                      mfl(k  )=mfl(k+1)+q
                    else
                      mfl(k  )= .5*q
                      mfl(k+1)=-mfl(k  )
                    endif
                  endif
                  changed=.true.
                endif
              endif
c
            enddo
c
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Compute the final mass fluxes.
c --- ------------------------------------------------------------------
c
          k=kmin
          if (abs(mfl(k+1)-mfl(k)).gt.
     .        1.e-14*max(epsil*scu2(i,j),abs(mfl(k+1)+mfl(k)))) then
            umfltd(i,j,2+mm)=mfl(k+1)-mfl(k)
            umfltd(i,j,1+mm)=umfltd(i,j,2+mm)
     .                       *dpu(i,j,1+nn)/(dpu(i,j,1+nn)
     .                                      +dpu(i,j,2+nn))
            umfltd(i,j,2+mm)=umfltd(i,j,2+mm)-umfltd(i,j,1+mm)
          else
            umfltd(i,j,1+mm)=0.
            umfltd(i,j,2+mm)=0.
          endif
          do k=kintr,kmax
            km=k+mm
            if (abs(mfl(k+1)-mfl(k)).gt.
     .          1.e-14*max(epsil*scu2(i,j),abs(mfl(k+1)+mfl(k)))) then
              umfltd(i,j,km)=mfl(k+1)-mfl(k)
            else
              umfltd(i,j,km)=0.
            endif
            if (umfltd(i,j,km).gt.
     .          ffac*max(epsil,dlm(k))*scp2(i-1,j)) then
              write (lp,*) 'eddtra_gm u >',i+i0,j+j0,k,umfltd(i,j,km),
     .                     ffac*max(epsil,dlm(k))*scp2(i-1,j)
              call xchalt('(eddtra_gm)')
                     stop '(eddtra_gm)'
            endif
            if (umfltd(i,j,km).lt.
     .         -ffac*max(epsil,dlp(k))*scp2(i  ,j)) then
              write (lp,*) 'eddtra_gm u <',i+i0,j+j0,k,umfltd(i,j,km),
     .                    -ffac*max(epsil,dlp(k))*scp2(i  ,j)
              call xchalt('(eddtra_gm)')
                     stop '(eddtra_gm)'
            endif
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute v-component of eddy-induced mass fluxes.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  l,i,k,km,et2mf,kmax,kn,kintr,kappa,upsilon,
c$OMP+  kmin,mfl,dlm,dlp,fhi,flo,changed,niter,kdir,q)
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
c
c --- ------------------------------------------------------------------
c --- --- Set eddy-induced mass fluxes to zero initially.
c --- ------------------------------------------------------------------
c
          do k=1,kk
            km=k+mm
            vmfltd(i,j,km)=0.
          enddo
c
c --- --- Eddy transport to mass flux conversion factor.
          et2mf=-g*rho0*delt1*scvx(i,j)
c
c --- --- Index of last layer containing mass at either of the scalar
c --- --- points adjacent to the velocity point.
          kmax=1
          do k=3,kk
            kn=k+nn
            if (dp(i,j-1,kn).gt.epsil.or.dp(i,j,kn).gt.epsil) kmax=k
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Proceed with mass flux computation if at least one of the
c --- --- adjacent scalar points to the velocity point has a mass
c --- --- containing interior layer. Mass fluxes will be assigned at
c --- --- layer interface corresponding to the eddy induced transport.
c --- --- The final layer mass flux will be the lower minus the upper
c --- --- interface flux.  The mass fluxes are limited to keep
c --- --- interfaces within the water column. There are 3 cases to
c --- --- consider:
c --- ---   Case 1: The mixed layer extends to the bottom at both
c --- ---           adjacent scalar points to the velocity point
c --- ---   Case 2: The mixed layer extends to the bottom at scalar
c --- ---           point (i,j  ).
c --- ---   Case 3: The mixed layer extends to the bottom at scalar
c --- ---           point (i,j-1).
c --- ---   Case 4: The mixed layer does not reach the bottom at neither
c --- ---           of the scalar points adjacent to the velocity point.
c --- ------------------------------------------------------------------
c
          if     (kfpla(i,j-1,n).gt.kk.and.kfpla(i,j,n).gt.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Case 1:
c --- ------------------------------------------------------------------
c
c --- ----- Keep the initial zero mass fluxes for this column.
            cycle
c
          elseif (kfpla(i,j-1,n).le.kk.and.kfpla(i,j,n).gt.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Case 2:
c --- ------------------------------------------------------------------
c
c --- ----- Find the index of the first layer at (i,j-1) that is
c --- ----- hydrostatically stable at the mixed layer base at (i,j  ).
            km=2+nn
            kintr=kfpla(i,j-1,n)
            kn=kintr+nn
            do while
     .        (rho(p(i,j  ,3),temp(i,j-1,kn),saln(i,j-1,kn)).lt.
     .         rho(p(i,j  ,3),temp(i,j  ,km),saln(i,j  ,km)).or.
     .         dp(i,j-1,kn).lt.epsil)
              kintr=kintr+1
              if (kintr.eq.kmax+1) exit
              kn=kintr+nn
            enddo
c
c --- ----- If a physical layer cannot be found, keep the initial zero
c --- ----- mass fluxes for this column.
            if (kintr.eq.kmax+1) cycle
c
c --- ----- Compute the eddy induced transport (upsilon) at the mixed
c --- ----- layer base.
            kappa=.5*(difint(i,j-1,2)+difint(i,j,2))
            upsilon(3)=-kappa*nslpy(i,j,3)
c
c --- ----- If the eddy-induced transport at the base of the mixed layer
c --- ----- would cause a negative mass flux below the mixed layer, keep
c --- ----- the initial zero mass fluxes for this column.
            if (upsilon(3).le.0.) cycle
c
c --- ----- Assign interface mass fluxes.
            kmin=kintr-1
            mfl(kmin)=0.
            mfl(kintr)=et2mf*upsilon(3)
            do k=kintr+1,kmax+1
              mfl(k)=0.
            enddo
c
          elseif (kfpla(i,j-1,n).gt.kk.and.kfpla(i,j,n).le.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Case 3:
c --- ------------------------------------------------------------------
c
c --- ----- Find the index of the first layer at (i,j  ) that is
c --- ----- hydrostatically stable at the mixed layer base at (i,j-1).
            km=2+nn
            kintr=kfpla(i,j  ,n)
            kn=kintr+nn
            do while
     .        (rho(p(i,j-1,3),temp(i,j  ,kn),saln(i,j  ,kn)).lt.
     .         rho(p(i,j-1,3),temp(i,j-1,km),saln(i,j-1,km)).or.
     .         dp(i,j  ,kn).lt.epsil)
              kintr=kintr+1
              if (kintr.eq.kmax+1) exit
              kn=kintr+nn
            enddo
c
c --- ----- If a physical layer cannot be found, keep the initial zero
c --- ----- mass fluxes for this column.
            if (kintr.eq.kmax+1) cycle
c
c --- ----- Compute the eddy induced transport (upsilon) at the mixed
c --- ----- layer base.
            kappa=.5*(difint(i,j-1,2)+difint(i,j,2))
            upsilon(3)=-kappa*nslpy(i,j,3)
c
c --- ----- If the eddy-induced transport at the base of the mixed layer
c --- ----- would cause a positive mass flux below the mixed layer, keep
c --- ----- the initial zero mass fluxes for this column.
            if (upsilon(3).ge.0.) cycle
c
c --- ----- Assign interface mass fluxes.
            kmin=kintr-1
            mfl(kmin)=0.
            mfl(kintr)=et2mf*upsilon(3)
            do k=kintr+1,kmax+1
              mfl(k)=0.
            enddo
c
          else
c
c --- ------------------------------------------------------------------
c --- ----- Case 4:
c --- ------------------------------------------------------------------
c
c --- ----- The first interior interface where the eddy induced
c --- ----- transport is estimated is at index kintr+1.
            kintr=max(kfpla(i,j-1,n),kfpla(i,j,n))
c
c --- ----- Compute the eddy induced transport (upsilon) at the mixed
c --- ----- layer base.
            kappa=.5*(difint(i,j-1,2)+difint(i,j,2))
            upsilon(3)=-kappa*nslpy(i,j,3)
c
c --- ----- Compute the eddy induced transport at interior interfaces.
            do k=kintr+1,kmax
              kn=k+nn
              kappa=.25*(difint(i,j-1,k-1)+difint(i,j  ,k-1)
     .                  +difint(i,j-1,k  )+difint(i,j  ,k  ))
              upsilon(k)=-kappa*nslpy(i,j,k)
            enddo
            upsilon(kmax+1)=0.
c
c --- ----- If the layer kintr-1 is a physical layer at either of the
c --- ----- adjacent scalar points to the velocity point, then apply
c --- ----- an upper interface mass flux corresponding to the eddy
c --- ----- induced transport at the mixed layer base and a lower
c --- ----- interface mass flux corresponding to the eddy induced
c --- ----- transport at the kintr+1 interface if this would lead to a
c --- ----- hydrostatically stable layer arrangement.
            km=2+nn
            kn=kintr-1+nn
            if ((kfpla(i,j-1,n).lt.kintr.and.
     .           upsilon(3)-upsilon(kintr+1).gt.0..and.
     .           rho(p(i,j  ,3),temp(i,j-1,kn),saln(i,j-1,kn)).gt.
     .           rho(p(i,j  ,3),temp(i,j  ,km),saln(i,j  ,km))).or.
     .          (kfpla(i,j  ,n).lt.kintr.and.
     .           upsilon(3)-upsilon(kintr+1).lt.0..and.
     .           rho(p(i,j-1,3),temp(i,j  ,kn),saln(i,j  ,kn)).gt.
     .           rho(p(i,j-1,3),temp(i,j-1,km),saln(i,j-1,km)))) then
              kintr=kintr-1
              upsilon(kintr+1)=upsilon(kintr+2)
            endif
c
c --- ----- Assign interface mass fluxes.
            kmin=kintr-1
            mfl(kmin)=0.
            mfl(kintr)=et2mf*upsilon(3)
            do k=kintr+1,kmax
              mfl(k)=et2mf*upsilon(k)
            enddo
            mfl(kmax+1)=0.
c
          endif
c
c --- ------------------------------------------------------------------
c --- --- Ensure that mass fluxes do not create negative layer
c --- --- thicknesses.
c --- ------------------------------------------------------------------
c
c --- --- Compute the layer thicknesses available to be depleted by mass
c --- --- fluxes at the scalar points adjacent to the velocity point.
c --- --- These bounded layer thicknesses are consistent with the
c --- --- transport algorithm.
          dlm(kmin)=max(0.,min(p(i,j-1,3),pbv(i,j,n))
     .                    -max(p(i,j-1,1),ptv(i,j)))
          dlp(kmin)=max(0.,min(p(i,j  ,3),pbv(i,j,n))
     .                    -max(p(i,j  ,1),ptv(i,j)))
          do k=kintr,kmax
            dlm(k)=max(0.,min(p(i,j-1,k+1),pbv(i,j,n))
     .                   -max(p(i,j-1,k  ),ptv(i,j)))
            dlp(k)=max(0.,min(p(i,j  ,k+1),pbv(i,j,n))
     .                   -max(p(i,j  ,k  ),ptv(i,j)))
          enddo
c
c --- --- If excessive depletion of layers occur beneath the mixed layer
c --- --- base, try to adjust interface fluxes other that the mixed
c --- --- layer base interface flux.
          fhi= fface
     .         *max(0.,min((p(i,j-1,3)-ptv(i,j))*scp2(i,j-1),
     .                     (pbv(i,j,n)-p(i,j  ,kintr))*scp2(i,j  )))
          flo=-fface
     .         *max(0.,min((p(i,j  ,3)-ptv(i,j))*scp2(i,j  ),
     .                     (pbv(i,j,n)-p(i,j-1,kintr))*scp2(i,j-1)))
          mfl(kmin+1)=min(fhi,max(flo,mfl(kmin+1)))
          do k=kmin+1,kmax-1
            if     (mfl(k+1)-mfl(k).gt.
     .              ffac*max(epsil,dlm(k))*scp2(i,j-1)) then
              mfl(k+1)=mfl(k  )+fface*dlm(k)*scp2(i,j-1)
            elseif (mfl(k+1)-mfl(k).lt.
     .             -ffac*max(epsil,dlp(k))*scp2(i,j  )) then
              mfl(k+1)=mfl(k  )-fface*dlp(k)*scp2(i,j  )
            else
              exit
            endif
          enddo
c
c --- --- Apply an iterative procedure for flux limiting by alternate
c --- --- upward and downward propagation through the layers.
c
          changed=.true.
          niter=0
          kdir=1
c
          do while (changed)
c
            niter=niter+1
            if (niter.eq.1000) then
              k=kmin
              write (lp,*)
              write (lp,'(i3,3e16.8)') 1,mfl(k+1),mfl(k),
     . (mfl(k+1)-mfl(k))
     . /(max(onemm,dpv(i,j,1+nn)+dpv(i,j,2+nn))*delt1*scvx(i,j))
              do k=kintr,kmax
                kn=k+nn
                write (lp,'(i3,3e16.8)') k,mfl(k+1),mfl(k),
     . (mfl(k+1)-mfl(k))/(max(onemm,dpv(i,j,kn))*delt1*scvx(i,j))
              enddo
              write (lp,*) 'no convergence v',i+i0,j+j0
              call xchalt('(eddtra_gm)')
                     stop '(eddtra_gm)'
            endif
c
            changed=.false.
            kdir=-kdir
c
            do k=((1-kdir)*kmax+(1+kdir)*kmin)/2,
     .           ((1-kdir)*kmin+(1+kdir)*kmax)/2,kdir
c
c --- ------- Proceed with flux limiting of this layer if the mass flux
c --- ------- difference between lower and upper interface is beyond the
c --- ------- floating point accuracy limitation.
              if (abs(mfl(k+1)-mfl(k)).gt.
     .            1.e-14*max(epsil*scv2(i,j),abs(mfl(k+1)+mfl(k)))) then
c
                if     (mfl(k+1)-mfl(k).gt.
     .                  ffac*max(epsil,dlm(k))*scp2(i,j-1)) then
c
c --- ----------- In this case, the mass fluxes are removing too much
c --- ----------- mass from the grid cell at (i,j-1,k). Limit the
c --- ----------- dominating interface flux.
                  q=fface*dlm(k)*scp2(i,j-1)
                  if (mfl(k+1).gt.-mfl(k)) then
                    if (mfl(k  ).gt.-.5*q) then
                      mfl(k+1)=mfl(k  )+q
                    else
                      mfl(k+1)= .5*q
                      mfl(k  )=-mfl(k+1)
                    endif
                  else
                    if (mfl(k+1).lt. .5*q) then
                      mfl(k  )=mfl(k+1)-q
                    else
                      mfl(k  )=-.5*q
                      mfl(k+1)=-mfl(k  )
                    endif
                  endif
                  changed=.true.
                elseif (mfl(k+1)-mfl(k).lt.
     .                 -ffac*max(epsil,dlp(k))*scp2(i,j  )) then
c
c --- ----------- In this case, the mass fluxes are removing too much
c --- ----------- mass from the grid cell at (i,j  ,k). Limit the
c --- ----------- dominating interface flux.
                  q=fface*dlp(k)*scp2(i,j  )
                  if (mfl(k+1).lt.-mfl(k)) then
                    if (mfl(k  ).lt. .5*q) then
                      mfl(k+1)=mfl(k  )-q
                    else
                      mfl(k+1)=-.5*q
                      mfl(k  )=-mfl(k+1)
                    endif
                  else
                    if (mfl(k+1).gt.-.5*q) then
                      mfl(k  )=mfl(k+1)+q
                    else
                      mfl(k  )= .5*q
                      mfl(k+1)=-mfl(k  )
                    endif
                  endif
                  changed=.true.
                endif
              endif
c
            enddo
c
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Compute the final mass fluxes.
c --- ------------------------------------------------------------------
c
          k=kmin
          if (abs(mfl(k+1)-mfl(k)).gt.
     .        1.e-14*max(epsil*scv2(i,j),abs(mfl(k+1)+mfl(k)))) then
            vmfltd(i,j,2+mm)=mfl(k+1)-mfl(k)
            vmfltd(i,j,1+mm)=vmfltd(i,j,2+mm)
     .                       *dpv(i,j,1+nn)/(dpv(i,j,1+nn)
     .                                      +dpv(i,j,2+nn))
            vmfltd(i,j,2+mm)=vmfltd(i,j,2+mm)-vmfltd(i,j,1+mm)
          else
            vmfltd(i,j,1+mm)=0.
            vmfltd(i,j,2+mm)=0.
          endif
          do k=kintr,kmax
            km=k+mm
            if (abs(mfl(k+1)-mfl(k)).gt.
     .          1.e-14*max(epsil*scv2(i,j),abs(mfl(k+1)+mfl(k)))) then
              vmfltd(i,j,km)=mfl(k+1)-mfl(k)
            else
              vmfltd(i,j,km)=0.
            endif
            if (vmfltd(i,j,km).gt.
     .          ffac*max(epsil,dlm(k))*scp2(i,j-1)) then
              write (lp,*) 'eddtra_gm v >',i+i0,j+j0,k,vmfltd(i,j,km),
     .                     ffac*max(epsil,dlm(k))*scp2(i,j-1)
              call xchalt('(eddtra_gm)')
                     stop '(eddtra_gm)'
            endif
            if (vmfltd(i,j,km).lt.
     .         -ffac*max(epsil,dlp(k))*scp2(i,j  )) then
              write (lp,*) 'eddtra_gm v <',i+i0,j+j0,k,vmfltd(i,j,km),
     .                    -ffac*max(epsil,dlp(k))*scp2(i,j  )
              call xchalt('(eddtra_gm)')
                     stop '(eddtra_gm)'
            endif
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
