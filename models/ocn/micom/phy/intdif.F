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
      integer i,j,k,l,km,kn
      real flxhi,flxlo,q
c
      call xctilr(difint, 1,kk, 2,2, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1m), 1,kk, 1,1, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1m), 1,kk, 1,1, halo_ps)
c
c$OMP PARALLEL DO
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
          umfltd(i,j,1+mm)=0.
          umfltd(i,j,2+mm)=0.
          umfltd(i,j,3+mm)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
          vmfltd(i,j,1+mm)=0.
          vmfltd(i,j,2+mm)=0.
          vmfltd(i,j,3+mm)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=4,kk
        km=k+mm
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
            umfltd(i,j,km-1)=umfltd(i,j,km-1)+q
            umfltd(i,j,km  )=-q
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
            vmfltd(i,j,km-1)=vmfltd(i,j,km-1)+q
            vmfltd(i,j,km  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c --- diagnose diffusive flux components of heat and salt
c
      do k=1,kk
        km=k+mm
c
c$OMP PARALLEL DO
        do j=1,jj
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
c$OMP END PARALLEL DO
c
      enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'intdif:'
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
