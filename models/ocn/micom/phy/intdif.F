      subroutine intdif(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate diffusive fluxes due to interface diffusion
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
      integer i,j,k,l,kn
      real flxhi,flxlo,q
c
      call xctilr(difint, 1,kk, 2,2, halo_ps)
c
c$OMP PARALLEL DO
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
          uflxdf(i,j,1)=0.
          uflxdf(i,j,2)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
          vflxdf(i,j,1)=0.
          vflxdf(i,j,2)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=3,kk
        kn=k+nn
c
c$OMP PARALLEL DO PRIVATE(flxhi,flxlo,q)
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            flxhi= .125*min(dp(i-1,j,kn-1)*scp2(i-1,j),
     .                      dp(i  ,j,kn  )*scp2(i  ,j))
            flxlo=-.125*min(dp(i  ,j,kn-1)*scp2(i  ,j),
     .                      dp(i-1,j,kn  )*scp2(i-1,j))
            q=.25*(difint(i-1,j,k-1)+difint(i,j,k-1)
     .            +difint(i-1,j,k  )+difint(i,j,k  ))
            q=min(flxhi,max(flxlo,
     .            delt1*q*(p(i-1,j,k)-p(i,j,k))*scuy(i,j)*scuxi(i,j)))
            uflxdf(i,j,k-1)=uflxdf(i,j,k-1)+q
            uflxdf(i,j,k  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(flxhi,flxlo,q)
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            flxhi= .125*min(dp(i,j-1,kn-1)*scp2(i,j-1),
     .                      dp(i,j  ,kn  )*scp2(i,j  ))
            flxlo=-.125*min(dp(i,j  ,kn-1)*scp2(i,j  ),
     .                      dp(i,j-1,kn  )*scp2(i,j-1))
            q=.25*(difint(i,j-1,k-1)+difint(i,j,k-1)
     .            +difint(i,j-1,k  )+difint(i,j,k  ))
            q=min(flxhi,max(flxlo,
     .            delt1*q*(p(i,j-1,k)-p(i,j,k))*scvx(i,j)*scvyi(i,j)))
            vflxdf(i,j,k-1)=vflxdf(i,j,k-1)+q
            vflxdf(i,j,k  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'intdif:'
        endif
        call chksummsk(uflxdf,iu,kk,'uflxdf')
        call chksummsk(vflxdf,iv,kk,'vflxdf')
      endif
c
      return
      end
