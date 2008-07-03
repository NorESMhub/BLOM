      subroutine icedif(h,scuy,scvx,scuxi,scvyi,scp2i,dt,cdiff,cdfac)
c
c --- ------------------------------------------------------------------
c --- Laplacian diffusion of sea ice scalar field with enhanced
c --- diffusivity in channels of one grid point width
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  h,scuy,scvx,scuxi,scvyi,scp2i
      real dt,cdiff,cdfac
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fu,fv
      real q
      integer i,j,l
c
c --- compute fluxes at grid cell boundaries
c$OMP PARALLEL DO
      do j=1,jj+1
        do i=1,ii+1
          fu(i,j)=0.
          fv(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(h, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO PRIVATE(q)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          q=max(0.,1.-iu(i,j-1)-iu(i,j+1))
          fu(i,j)=cdiff*(cdfac*q+1.-q)*(h(i,j)-h(i-1,j))
     .            *scuy(i,j)*scuxi(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(q)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          q=max(0.,1.-iv(i-1,j)-iv(i+1,j))
          fv(i,j)=cdiff*(cdfac*q+1.-q)*(h(i,j)-h(i,j-1))
     .            *scvx(i,j)*scvyi(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- update scalar field (note that diffusivity is in m^2/s while scale
c --- factors are in cm)
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          h(i,j)=h(i,j)+(fu(i+1,j)-fu(i,j)+fv(i,j+1)-fv(i,j))
     .                  *dt*scp2i(i,j)*1.e4
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
