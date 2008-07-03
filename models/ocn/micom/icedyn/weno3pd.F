      subroutine weno3pd(g,sao,u,v,scuy,scvx,scp2i,scp2,dt)
c
c --- ------------------------------------------------------------------
c --- By a weighted essentially non-oscillatory scheme with up to 3th
c --- order accuracy, obtain the spatial advective operator of a
c --- 2-dimensional field defined at the scalar points of a C-grid. The
c --- fluxes are limited to make the scheme positive definite.
c --- Advective velocities in the i- and j-direction are defined at u-
c --- and v-points, respectively.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  g,sao,u,v,scuy,scvx,scp2i,scp2
      real dt
c
      real cq00,cq01,cq10,cq11,ca0,ca1,eps
      parameter (cq00=-1./2.,cq01= 3./2.,
     .           cq10= 1./2.,cq11= 1./2.,
     .           ca0=1./3.,ca1=2./3.,
     .           eps=1.e-12)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ful,fuh,fvl,fvh,gt
      real q0,q1,a0,a1,q
      integer i,j,l,im1,im2,ip1,jm1,jm2,jp1
c
c --- Compute grid cell boundary fluxes. Split in a low order flux
c --- (donor cell) and a high order correction flux.
c
c$OMP PARALLEL DO
      do j=0,jj+2
        do i=0,ii+2
          ful(i,j)=0.
          fuh(i,j)=0.
          fvl(i,j)=0.
          fvh(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(g, 1,1, 3,3, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(im1,im2,q0,q1,a0,a1,ip1)
      do j=0,jj+1
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
          im1=i-1
c
          if (u(i,j).gt.0.) then
            im2=im1-iu(im1,j)
c
            q0=cq00*g(im2,j)+cq01*g(im1,j)
            q1=cq10*g(im1,j)+cq11*g(i  ,j)
c
            a0=ca0
            a1=ca1*(abs(g(im2,j)-g(im1,j))+eps)
     .            /(abs(g(im1,j)-g(i  ,j))+eps)
c
            ful(i,j)=u(i,j)*g(im1,j)*scuy(i,j)
c
          else
            ip1=i+iu(i+1,j)
c
            q0=cq11*g(im1,j)+cq10*g(i  ,j)
            q1=cq01*g(i  ,j)+cq00*g(ip1,j)
c
            a0=ca1
            a1=ca0*(abs(g(im1,j)-g(i  ,j))+eps)
     .            /(abs(g(i  ,j)-g(ip1,j))+eps)
c
            ful(i,j)=u(i,j)*g(i  ,j)*scuy(i,j)
c
          endif
c
          fuh(i,j)=u(i,j)*(a0*q0+a1*q1)/(a0+a1)*scuy(i,j)-ful(i,j)
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(jm1,q0,q1,a0,a1,jm2,jp1)
      do j=0,jj+2
        jm1=j-1
        do l=1,isv(j)
        do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
c
          if (v(i,j).gt.0.) then
            jm2=jm1-iv(i,jm1)
c
            q0=cq00*g(i,jm2)+cq01*g(i,jm1)
            q1=cq10*g(i,jm1)+cq11*g(i,j  )
c
            a0=ca0
            a1=ca1*(abs(g(i,jm2)-g(i,jm1))+eps)
     .            /(abs(g(i,jm1)-g(i,j  ))+eps)
c
            fvl(i,j)=v(i,j)*g(i,jm1)*scvx(i,j)
c
          else
            jp1=j+iv(i,j+1)
c
            q0=cq11*g(i,jm1)+cq10*g(i,j  )
            q1=cq01*g(i,j  )+cq00*g(i,jp1)
c
            a0=ca1
            a1=ca0*(abs(g(i,jm1)-g(i,j  ))+eps)
     .            /(abs(g(i,j  )-g(i,jp1))+eps)
c
            fvl(i,j)=v(i,j)*g(i,j  )*scvx(i,j)
c
          endif
c
          fvh(i,j)=v(i,j)*(a0*q0+a1*q1)/(a0+a1)*scvx(i,j)-fvl(i,j)
c        
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Update field with low order fluxes.
c$OMP PARALLEL DO
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          gt(i,j)=g(i,j)-dt*(ful(i+1,j)-ful(i,j)
     .                      +fvl(i,j+1)-fvl(i,j))*scp2i(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Obtain fluxes with limited high order correction fluxes.
      q=.25/dt
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          fuh(i,j)=ful(i,j)+max(-q*gt(i  ,j)*scp2(i  ,j),
     .                      min( q*gt(i-1,j)*scp2(i-1,j),fuh(i,j)))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          fvh(i,j)=fvl(i,j)+max(-q*gt(i,j  )*scp2(i,j  ),
     .                      min( q*gt(i,j-1)*scp2(i,j-1),fvh(i,j)))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute the spatial advective operator.
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          sao(i,j)=-(fuh(i+1,j)-fuh(i,j)+fvh(i,j+1)-fvh(i,j))*scp2i(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
