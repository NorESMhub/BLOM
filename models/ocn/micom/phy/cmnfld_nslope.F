      subroutine cmnfld_nslope(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Estimate slope of local neutral surface.
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
      real rho0,pm,rho_x,phi_x,bfsqm,rho_y,phi_y
      integer i,j,k,l,kn,kintr,kmax
c
      real p_alpha,rho
      external p_alpha,rho
c
c --- ------------------------------------------------------------------
c --- Compute geopotential at layer interfaces.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=-1,jj+2
        do k=kk,1,-1
          kn=k+nn
          do l=1,isp(j)
          do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
            if (dp(i,j,kn).lt.epsil) then
              phi(i,j,k)=phi(i,j,k+1)
            else
              phi(i,j,k)=phi(i,j,k+1)
     .                  -p_alpha(p(i,j,k+1),p(i,j,k),
     .                           temp(i,j,kn),saln(i,j,kn))
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute slope vector of local neutral surfaces.
c --- ------------------------------------------------------------------
c
      rho0=1./alpha0
c
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
c --- --- Index of last layer containing mass at either of the scalar
c --- --- points adjacent to the velocity point.
          kmax=1
          do k=3,kk
            kn=k+nn
            if (dp(i-1,j,kn).gt.epsil.or.dp(i,j,kn).gt.epsil) kmax=k
          enddo
c
          if     (kfpla(i-1,j,n).le.kk.or.kfpla(i,j,n).le.kk) then
c
c --- ----- The first interior interface where the x-component of the
c --- ----- slope vector is estimated is at index kintr+1.
            kintr=max(kfpla(i-1,j,n),kfpla(i,j,n))
c
c --- ----- Compute the x-component of the slope vector at the mixed
c --- ----- layer base.
            pm=.5*(p(i-1,j,3)+p(i,j,3))
            rho_x=rho(pm,temp(i  ,j,2+nn),saln(i  ,j,2+nn))
     .           -rho(pm,temp(i-1,j,2+nn),saln(i-1,j,2+nn))
            phi_x=phi(i,j,3)-phi(i-1,j,3)
            bfsqm=.5*(bfsqf(i-1,j,3)+bfsqf(i,j,3))
            nslpx(i,j,3)=(g*rho_x/(rho0*bfsqm)+phi_x/g)*scuxi(i,j)
            nnslpx(i,j,3)=sqrt(bfsqm)*nslpx(i,j,3)
c
c --- ----- Compute the x-component of the slope vector at interior
c --- ----- interfaces.
            do k=kintr+1,kmax
              kn=k+nn
              pm=.5*(p(i-1,j,k)+p(i,j,k))
              rho_x=.5*(rho(pm,temp(i  ,j,kn-1),saln(i  ,j,kn-1))
     .                 -rho(pm,temp(i-1,j,kn-1),saln(i-1,j,kn-1))
     .                 +rho(pm,temp(i  ,j,kn  ),saln(i  ,j,kn  ))
     .                 -rho(pm,temp(i-1,j,kn  ),saln(i-1,j,kn  )))
              phi_x=phi(i,j,k)-phi(i-1,j,k)
              bfsqm=.5*(bfsqf(i-1,j,k)+bfsqf(i,j,k))
              nslpx(i,j,k)=(g*rho_x/(rho0*bfsqm)+phi_x/g)*scuxi(i,j)
              nnslpx(i,j,k)=sqrt(bfsqm)*nslpx(i,j,k)
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
c
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
c --- --- Index of last layer containing mass at either of the scalar
c --- --- points adjacent to the velocity point.
          kmax=1
          do k=3,kk
            kn=k+nn
            if (dp(i,j-1,kn).gt.epsil.or.dp(i,j,kn).gt.epsil) kmax=k
          enddo
c
          if     (kfpla(i,j-1,n).le.kk.or.kfpla(i,j,n).le.kk) then
c
c --- ----- The first interior interface where the y-component of the
c --- ----- slope vector is estimated is at index kintr+1.
            kintr=max(kfpla(i,j-1,n),kfpla(i,j,n))
c
c --- ----- Compute the y-component of the slope vector at the mixed
c --- ----- layer base.
            pm=.5*(p(i,j-1,3)+p(i,j,3))
            rho_y=rho(pm,temp(i,j  ,2+nn),saln(i,j  ,2+nn))
     .           -rho(pm,temp(i,j-1,2+nn),saln(i,j-1,2+nn))
            phi_y=phi(i,j,3)-phi(i,j-1,3)
            bfsqm=.5*(bfsqf(i,j-1,3)+bfsqf(i,j,3))
            nslpy(i,j,3)=(g*rho_y/(rho0*bfsqm)+phi_y/g)*scvyi(i,j)
            nnslpy(i,j,3)=sqrt(bfsqm)*nslpy(i,j,3)
c
c --- ----- Compute the y-component of the slope vector at interior
c --- ----- interfaces.
            do k=kintr+1,kmax
              kn=k+nn
              pm=.5*(p(i,j-1,k)+p(i,j,k))
              rho_y=.5*(rho(pm,temp(i,j  ,kn-1),saln(i,j  ,kn-1))
     .                 -rho(pm,temp(i,j-1,kn-1),saln(i,j-1,kn-1))
     .                 +rho(pm,temp(i,j  ,kn  ),saln(i,j  ,kn  ))
     .                 -rho(pm,temp(i,j-1,kn  ),saln(i,j-1,kn  )))
              phi_y=phi(i,j,k)-phi(i,j-1,k)
              bfsqm=.5*(bfsqf(i,j-1,k)+bfsqf(i,j,k))
              nslpy(i,j,k)=(g*rho_y/(rho0*bfsqm)+phi_y/g)*scvyi(i,j)
              nnslpy(i,j,k)=sqrt(bfsqm)*nslpy(i,j,k)
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
c
      return
      end
