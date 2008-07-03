      subroutine iceadv(h,u,v,scuy,scvx,scp2i,scp2,dt)
c
c --- ------------------------------------------------------------------
c --- Advection is done with flux limited 3rd order WENO in space and
c --- 2nd order Runge-Kutta in time
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  h,u,v,scuy,scvx,scp2i,scp2
      real dt
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sao,hp
      real dtm
      integer i,j,l
c
c --- Use a modified time step since velocities are in m/s while scale
c --- factors are in cm
      dtm=dt*1.e2
c
c --- Prediction step
      call weno3pd(h,sao,u,v,scuy,scvx,scp2i,scp2,dtm)
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          hp(i,j)=h(i,j)+dtm*sao(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Correction step
      call weno3pd(hp,sao,u,v,scuy,scvx,scp2i,scp2,dtm)
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          h(i,j)=.5*(h(i,j)+hp(i,j)+dtm*sao(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
