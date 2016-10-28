      subroutine cmnfld_bfsqf(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Compute filtered buoyancy frequency squared.
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
c ---            buoyancy frequency squared [g/cm/s**2].
c ---   slsmfq - factor to be multiplied with the mixed layer depth to
c ---            find the smoothing length scale at the base of the
c ---            mixed layer in the computation of buoyancy frequency
c ---            squared [].
c ---   slsels - factor to be multiplied with the mixed layer depth to
c ---            find the e-folding length scale of the smoothing length
c ---            scale in the computation of buoyancy frequency squared
c ---            [].
c ---   bfsqmn - minimum value of buoyancy frequency squared used in the
c ---            computation of neutral slope [s-2].
c
      real sls0,slsmfq,slsels,bfsqmn
      parameter (sls0=10.*98060.,slsmfq=2.,slsels=2.,bfsqmn=1.e-7)
c
      real, dimension(kdm) :: delp,bfsq,sls2,atd,btd,ctd,rtd,gam
      real rho0,pml,q,pup,tup,sup,plo,tlo,slo,bei,pm,bfsqm
      integer i,j,k,l,kn,kfpl,kintr,kmax
c
      real rho
      external rho
c
c --- ------------------------------------------------------------------
c --- The buoyancy frequency squared is esitmated locally at each layer
c --- interface and smoothed in the vertical direction by solving a
c --- diffusion equation. At the mixed layer base the diffusion length
c --- scale is set to the maximum of sls0 and mixed layer depth (MLD)
c --- times slsmfq. Below the mixed layer, the diffusion length scale
c --- approaches sls0 with an e-folding length scale of MLD times
c --- slsels.
c --- ------------------------------------------------------------------
c
      do j=-1,jj+2
        do l=1,isp(j)
        do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
c
          kfpl=kfpla(i,j,n)
c
          if (kfpl.gt.kk) then
c
c --- ----- If the mixed layer extends to the bottom, set the buoyancy
c --- ----- frequency to a minimum value.
            do k=1,kk
              bfsqf(i,j,k)=bfsqmn
            enddo
c
          else
c
c --- ----- At layer interfaces, compute buoyancy frequency squared and
c --- ----- length scale for the subsequent smoothing.
            pml=max(.5*(p(i,j,3)+p(i,j,1)),
     .              .5*(3.*p(i,j,3)-p(i,j,kfpl+1)))
            delp(kfpl-1)=pml-p(i,j,1)
            bfsq(kfpl-1)=bfsqmn
            q=max(sls0,delp(kfpl-1)*slsmfq)
            sls2(kfpl-1)=q*q
            pup=pml
            tup=temp(i,j,2+nn)
            sup=saln(i,j,2+nn)
            do k=kfpl,kk
              kn=k+nn
              if (p(i,j,kk+1)-p(i,j,k).lt.epsil) then
                delp(k)=onemm
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
                tlo=temp(i,j,kn)
                slo=saln(i,j,kn)
                delp(k)=max(onemm,plo-pup)
                bfsq(k)=max(bfsqmn,
     .                      g*g*max(0.,rho(p(i,j,k),tlo,slo)
     .                                -rho(p(i,j,k),tup,sup))/delp(k))
                q=exp(-(p(i,j,k)-pml)/(slsels*delp(kfpl-1)))
                q=max(sls0,delp(kfpl-1)*slsmfq*q+sls0*(1.-q))
                sls2(k)=q*q
                pup=plo
                tup=tlo
                sup=slo
              endif
            enddo
c
c --- ----- Compute the coefficients for the tridiagonal set of
c --- ----- equations arising from the implicit backward discretization.
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
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'cmnfld_bfsqf:'
        endif
        call chksummsk(bfsqf,ip,kk,'bfsqf')
      endif
c
      return
      end
