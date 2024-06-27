! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2024 Mats Bentsen, Lars Inge Enstad, Mariana Vertenstein
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

module mod_barotp

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to time
  ! integration of the barotropic equations.
  ! ------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm, itdm, nreg
  use mod_types,     only: r8
  use mod_constants, only: spval, onem
  use mod_time,      only: lstep, dlt
  use mod_xc,        only: xctilr, xcstop, ii, jj, kk, isq, ifq, ilq, &
                           ii, jj, isp, ifp, ilp, jsp, jfp, jlp, lp, &
                           isu, ifu, ilu, ilv, isv, ifv, i0, ip, iu, iv, jpr, &
                           iu, iv, mnproc, nproc, nbdy, &
                           halo_ps, halo_us, halo_vs, halo_uv, halo_vv, halo_qs
  use mod_grid,      only: scuy, scvx, scp2i, scuxi, scuyi, scvxi, scvyi, &
                           corioq
  use mod_state,     only: u, v, ubflxs, vbflxs, ub, vb, pb, pbu, pbv, &
                           ubflxs_p, vbflxs_p, pb_p, pbu_p, pbv_p, &
                           ubcors_p, vbcors_p
  use mod_pgforc,    only: pgfxm, pgfym, xixp, xixm, xiyp, xiym, &
                           pgfxm_o, pgfym_o, &
                           xixp_o, xixm_o, xiyp_o, xiym_o
  use mod_momtum,    only: mommth
  use mod_tmsmt,     only: wbaro
  use mod_utility,   only: utotn, vtotn, umax, vmax
  use mod_checksum,  only: csdiag, chksummsk

  implicit none
  private

  public :: cwbdts, cwbdls, ubflx, vbflx, pb_mn, ubflx_mn, vbflx_mn, &
            pvtrop, inivar_barotp, barotp

  ! Variables to be set in namelist:
  real(r8) :: cwbdts   ! Coastal wave breaking damping resiprocal time scale [s-1]
  real(r8) :: cwbdls    ! Coastal wave breaking damping length scale [m].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  ubflx    ! u-component of barotropic mass flux [g cm s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  vbflx    ! v-component of barotropic mass flux [g cm s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  pb_mn    ! Bottom pressure [g cm-1 s-2].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  ubflx_mn ! u-component of barotropic mass flux [g cm s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  vbflx_mn ! v-component of barotropic mass flux [g cm s-3].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  pvtrop   ! Potential vorticity of barotropic flow [cm s g-1].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::   pvtrop_o  ! Potential vorticity of barotropic flow at old time
                                                                      ! level [cm s g-1].

contains

  subroutine inivar_barotp()

    ! ------------------------------------------------------------------
    ! Initialize arrays.
    ! ------------------------------------------------------------------

    ! Local variables
    integer :: i,j,k,l

    !$omp parallel do private(i)
    do j = 1-nbdy,jj+nbdy
      do i = 1-nbdy,ii+nbdy
        do k = 1,2
          ubflx(i,j,k) = spval
          vbflx(i,j,k) = spval
          pb_mn(i,j,k) = spval
          ubflx_mn(i,j,k) = spval
          vbflx_mn(i,j,k) = spval
          pvtrop(i,j,k) = spval
        end do
        pvtrop_o(i,j) = spval
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,k)
    do j = 1,jj+1
      do l = 1,isq(j)
        do i = max(1,ifq(j,l)),min(ii+1,ilq(j,l))
          do k = 1,2
            pb_mn(i  ,j  ,k) = 0.
            pb_mn(i-1,j  ,k) = 0.
            pb_mn(i  ,j-1,k) = 0.
            pb_mn(i-1,j-1,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(pb_mn,1,   2, nbdy,nbdy, halo_ps)

    ! initialize  ubflx,ubflx_mn  at points located upstream and
    ! downstream (in i direction) of p points.

    !$omp parallel do private(l,i,k)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
          do k = 1,2
            ubflx(i,j,k) = 0.
            ubflx_mn(i,j,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(ubflx,  1,   2, nbdy,nbdy, halo_us)
    call xctilr(ubflx_mn, 1,   2, nbdy,nbdy, halo_us)

    ! initialize  vbflx,vbflx_mn  at points located upstream and
    ! downstream (in j direction) of p points.

    !$omp parallel do private(l,j,k)
    do i = 1,ii
      do l = 1,jsp(i)
        do j = max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
          do k = 1,2
            vbflx(i,j,k) = 0.
            vbflx_mn(i,j,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(vbflx,  1,   2, nbdy,nbdy, halo_vs)
    call xctilr(vbflx_mn, 1,   2, nbdy,nbdy, halo_vs)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'inivar_barotp:'
      end if
      ! call chksummsk(pb_mn,ip,2,'pb')
      ! call chksummsk(ubflx,iu,2,'ubflx')
      ! call chksummsk(vbflx,iv,2,'vbflx')
      ! call chksummsk(ubflx_mn,iu,2,'ubflx')
      ! call chksummsk(vbflx_mn,iv,2,'vbflx')
    end if

  end subroutine inivar_barotp

  ! ------------------------------------------------------------------

  subroutine barotp(m,n,mm,nn,k1m,k1n)

    ! Local variables
    integer :: m,n,mm,nn,k1m,k1n

    ! Local variables
    ! Note - the following equal declaration automatically makes them a 'save' variable
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: pb_t    =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: ubflx_t =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: vbflx_t =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: umaxb   =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: uminb   =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: vmaxb   =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: vminb   =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: uglue   =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: vglue   =spval
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: ubflxs_t
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: vbflxs_t
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: ubcors_t
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)   :: vbcors_t
    real    :: q,woa,wob,wna,wnb,wo,wm,wn,utndcy,vtndcy
    integer :: i,j,k,l,kn,nb,lll0,lll,ml,nl,ll

    ! ------------------------------------------------------------------
    ! initialize barotropic velocity sums, determine maximum allowable
    ! barotropic velocities, and determine coefficients for coastal wave
    ! breaking parameterization
    ! ------------------------------------------------------------------

    !$omp parallel do private(l,i,k,kn)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          umaxb(i,j) = 0.
          uminb(i,j) = 0.
          uglue(i,j) = cwbdts*exp(1.-pbu(i,j,m)/(cwbdls*onem))
        end do
      end do
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            umaxb(i,j) = max(umaxb(i,j),u(i,j,kn))
            uminb(i,j) = min(uminb(i,j),u(i,j,kn))
          end do
        end do
      end do
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          umaxb(i,j) = (umax(i,j)-umaxb(i,j))*pbu(i,j,m)*scuy(i,j)
          uminb(i,j) = (umax(i,j)+uminb(i,j))*pbu(i,j,m)*scuy(i,j)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          vmaxb(i,j) = 0.
          vminb(i,j) = 0.
          vglue(i,j) = cwbdts*exp(1.-pbv(i,j,m)/(cwbdls*onem))
        end do
      end do
      do k = 1,kk
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vmaxb(i,j) = max(vmaxb(i,j),v(i,j,kn))
            vminb(i,j) = min(vminb(i,j),v(i,j,kn))
          end do
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          vmaxb(i,j) = (vmax(i,j)-vmaxb(i,j))*pbv(i,j,m)*scvx(i,j)
          vminb(i,j) = (vmax(i,j)+vminb(i,j))*pbv(i,j,m)*scvx(i,j)
        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------------
    ! potential vorticity of barotropic flow
    ! ------------------------------------------------------------------

    !$omp parallel do private(i)
    do j = -2,jj+3
      do i = 0,ii+1
        pvtrop_o(i,j) = pvtrop(i,j,n)
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i,q)
    do j = 0,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          q = 2./(pb_p(i,j)+pb_p(i-1,j))
          pvtrop(i,j  ,n) = corioq(i,j  )*q
          pvtrop(i,j+1,n) = corioq(i,j+1)*q
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i,q)
    do j = 1,jj
      do l = 1,isv(j)
        do i = max(0,ifv(j,l)),min(ii,ilv(j,l))
          q = 2./(pb_p(i,j)+pb_p(i,j-1))
          pvtrop(i  ,j,n) = corioq(i  ,j)*q
          pvtrop(i+1,j,n) = corioq(i+1,j)*q
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isq(j)
        do i = max(1,ifq(j,l)),min(ii,ilq(j,l))
          pvtrop(i,j,n) = corioq(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  ) &
               +pb_p(i,j-1)+pb_p(i-1,j-1))
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(uglue, 1,1, 1,2, halo_us)
    call xctilr(utotn, 1,1, 1,2, halo_uv)
    call xctilr(umaxb, 1,1, 1,2, halo_us)
    call xctilr(uminb, 1,1, 1,2, halo_us)
    call xctilr(vglue, 1,1, 1,2, halo_vs)
    call xctilr(vtotn, 1,1, 1,2, halo_vv)
    call xctilr(vmaxb, 1,1, 1,2, halo_vs)
    call xctilr(vminb, 1,1, 1,2, halo_vs)
    call xctilr(pvtrop(1-nbdy,1-nbdy,n), 1,1, 1,3, halo_qs)
    call xctilr(pgfxm(1-nbdy,1-nbdy,n), 1,1, 1,2, halo_uv)
    call xctilr(xixp(1-nbdy,1-nbdy,n), 1,1, 1,2, halo_us)
    call xctilr(xixm(1-nbdy,1-nbdy,n), 1,1, 1,2, halo_us)
    call xctilr(pgfym(1-nbdy,1-nbdy,n), 1,1, 1,2, halo_vv)
    call xctilr(xiyp(1-nbdy,1-nbdy,n), 1,1, 1,2, halo_vs)
    call xctilr(xiym(1-nbdy,1-nbdy,n), 1,1, 1,2, halo_vs)

    ! with arctic patch, switch umaxb and uminb, vmaxb and vminb, xixp
    ! and xixm, and xiyp and xiym in the halo region adjacent to the
    ! arctic grid intersection
    if (nreg == 2.and.nproc == jpr) then
      do j = jj,jj+2
        do i = 0,ii+1
          q = umaxb(i,j)
          umaxb(i,j) = uminb(i,j)
          uminb(i,j) = q
          q = xixp(i,j,n)
          xixp(i,j,n) = xixm(i,j,n)
          xixm(i,j,n) = q
        end do
      end do
      do i = max(0,itdm/2-i0+1),ii+1
        q = vmaxb(i,jj)
        vmaxb(i,jj) = vminb(i,jj)
        vminb(i,jj) = q
        q = xiyp(i,jj,n)
        xiyp(i,jj,n) = xiym(i,jj,n)
        xiym(i,jj,n) = q
      end do
      do j = jj+1,jj+2
        do i = 0,ii+1
          q = vmaxb(i,j)
          vmaxb(i,j) = vminb(i,j)
          vminb(i,j) = q
          q = xiyp(i,j,n)
          xiyp(i,j,n) = xiym(i,j,n)
          xiym(i,j,n) = q
        end do
      end do
    end if

    ! ------------------------------------------------------------------
    ! advance barotropic equations from baroclinic time level -m- to
    ! level -n- then advance barotropic equations another baroclinic
    ! time level so that the average barotropic transport for a
    ! leap-frog baroclinic step can be predicted
    ! ------------------------------------------------------------------

    do nb = 1,5

      if     (nb == 1) then
        lll0 = 1
        ml = 1
        nl = 2
        woa = -1./lstep
        wob = .5+(lll0-.5)/lstep
        wna = 0.
        wnb = 0.
        !$omp parallel do private(i)
        do j = 1,jj
          do i = 1,ii
            pb_t(i,j,ml) = pb_mn(i,j,ml)
            pb_t(i,j,nl) = pb_mn(i,j,nl)
            ubflx_t(i,j,ml) = ubflx_mn(i,j,ml)
            ubflx_t(i,j,nl) = ubflx_mn(i,j,nl)
            vbflx_t(i,j,ml) = vbflx_mn(i,j,ml)
            vbflx_t(i,j,nl) = vbflx_mn(i,j,nl)
          end do
        end do
        !$omp end parallel do
      else if (nb == 2) then
        woa = 0.
        wob = 0.
        wna = 1./lstep
        wnb = -(lll0-.5)/lstep
      else if (nb == 4) then
        wna = 0.
        wnb = 1.
      end if

      !$omp parallel do private(l,i)
      do j = -1,jj+2
        do l = 1,isu(j)
          do i = max(0,ifu(j,l)),min(ii+1,ilu(j,l))
            ubflxs_t(i,j) = 0.
            ubcors_t(i,j) = 0.
          end do
        end do
      end do
      !$omp end parallel do

      !$omp parallel do private(l,i)
      do j = 0,jj+2
        do l = 1,isv(j)
          do i = max(0,ifv(j,l)),min(ii,ilv(j,l))
            vbflxs_t(i,j) = 0.
            vbcors_t(i,j) = 0.
          end do
        end do
      end do
      !$omp end parallel do

      ! explicit time integration of barotropic flow (forward-backward
      ! scheme) in order to combine forward-backward scheme with
      ! leapfrog treatment of coriolis term, v-eqn must be solved before
      ! u-eqn every other time step

      do lll = lll0,lll0+lstep/2-1

        if (mod(lll,2) == 1) then

          wo = woa*lll+wob
          wn = wna*lll+wnb
          wm = 1.-wo-wn

          call xctilr(pb_t, 1,2, 2,2, halo_ps)
          call xctilr(ubflx_t, 1,2, 2,2, halo_uv)
          call xctilr(vbflx_t, 1,2, 2,3, halo_vv)

          ! continuity equation

          !$omp parallel do private(l,i)
          do j = -1,jj+2
            do l = 1,isp(j)
              do i = max(-1,ifp(j,l)),min(ii+1,ilp(j,l))
                pb_t(i,j,nl) = (1.-wbaro)*pb_t(i,j,ml)+wbaro*pb_t(i,j,nl) &
                              -(1.+wbaro)*dlt*(ubflx_t(i+1,j,ml)-ubflx_t(i,j,ml) &
                                              +vbflx_t(i,j+1,ml)-vbflx_t(i,j,ml)) &
                               *scp2i(i,j)
              end do
            end do
          end do
          !$omp end parallel do

          ! u momentum equation

          if (mommth == 'enscon') then

            ! Sadourny (1975) enstrophy conserving scheme

            !$omp parallel do private(l,i,q,utndcy)
            do j = -1,jj+2
              do l = 1,isu(j)
                do i = max(0,ifu(j,l)),min(ii+1,ilu(j,l))

                  ubflxs_t(i,j) = ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl) &
                                 +(1.+wbaro)*ubflx_t(i,j,ml)

                  q= (vbflx_t(i  ,j  ,ml)*scvxi(i  ,j  ) &
                     +vbflx_t(i  ,j+1,ml)*scvxi(i  ,j+1) &
                     +vbflx_t(i-1,j  ,ml)*scvxi(i-1,j  ) &
                     +vbflx_t(i-1,j+1,ml)*scvxi(i-1,j+1)) &
                     *(wo*(pvtrop_o(i,j)+pvtrop_o(i,j+1)) &
                     +wm*(pvtrop(i,j,m)+pvtrop(i,j+1,m)) &
                     +wn*(pvtrop(i,j,n)+pvtrop(i,j+1,n)))*.125

                  ubcors_t(i,j) = ubcors_t(i,j)+q

                  utndcy = q+ &
                       (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl) &
                           -xixm_o(i,j)*pb_t(i-1,j,nl))) &
                       +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl) &
                           -xixm(i,j,m)*pb_t(i-1,j,nl))) &
                       +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl) &
                           -xixm(i,j,n)*pb_t(i-1,j,nl)))) &
                       *scuxi(i,j)

                  ubflx_t(i,j,nl)= &
                       (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl) &
                       +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j) &
                                        *min(pb_t(i-1,j,nl),pb_t(i,j,nl)) &
                                        -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl) = max(-uminb(i,j),min(umaxb(i,j), &
                                        ubflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else if (mommth == 'enecon'.or.mommth == 'enedis') then

            ! Sadourny (1975) energy conserving scheme

            !$omp parallel do private(l,i,q,utndcy)
            do j = -1,jj+2
              do l = 1,isu(j)
                do i = max(0,ifu(j,l)),min(ii+1,ilu(j,l))

                  ubflxs_t(i,j) = ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl) &
                                +(1.+wbaro)*ubflx_t(i,j,ml)

                  q= .25*( (vbflx_t(i  ,j  ,ml)*scvxi(i  ,j  ) &
                           +vbflx_t(i-1,j  ,ml)*scvxi(i-1,j  )) &
                          *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m) &
                           +wn*pvtrop(i,j,n)) &
                          +(vbflx_t(i  ,j+1,ml)*scvxi(i  ,j+1) &
                          + vbflx_t(i-1,j+1,ml)*scvxi(i-1,j+1)) &
                          *(wo*pvtrop_o(i,j+1)+wm*pvtrop(i,j+1,m) &
                           +wn*pvtrop(i,j+1,n)))

                  ubcors_t(i,j) = ubcors_t(i,j)+q

                  utndcy = q+ &
                       (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl) &
                                         -xixm_o(i,j)*pb_t(i-1,j,nl))) &
                       +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl) &
                                         -xixm(i,j,m)*pb_t(i-1,j,nl))) &
                       +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl) &
                                         -xixm(i,j,n)*pb_t(i-1,j,nl)))) &
                       *scuxi(i,j)

                  ubflx_t(i,j,nl)= &
                       (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl) &
                      +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j) &
                                       *min(pb_t(i-1,j,nl),pb_t(i,j,nl)) &
                                       -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl) = max(-uminb(i,j),min(umaxb(i,j), &
                                        ubflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else
            if (mnproc == 1) then
              write (lp,'(3a)') ' mommth = ',trim(mommth),' is unsupported!'
            end if
            call xcstop('(barotp)')
            stop '(barotp)'
          end if

          ! v momentum equation

          if (mommth == 'enscon') then

            ! Sadourny (1975) enstrophy conserving scheme

            !$omp parallel do private(l,i,q,vtndcy)
            do j = 0,jj+2
              do l = 1,isv(j)
                do i = max(0,ifv(j,l)),min(ii,ilv(j,l))

                  vbflxs_t(i,j) = vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl) &
                       +(1.+wbaro)*vbflx_t(i,j,ml)

                  q = -(ubflx_t(i  ,j  ,nl)*scuyi(i  ,j  ) &
                       +ubflx_t(i+1,j  ,nl)*scuyi(i+1,j  ) &
                       +ubflx_t(i  ,j-1,nl)*scuyi(i  ,j-1) &
                       +ubflx_t(i+1,j-1,nl)*scuyi(i+1,j-1)) &
                       *(wo*(pvtrop_o(i,j)+pvtrop_o(i+1,j)) &
                        +wm*(pvtrop(i,j,m)+pvtrop(i+1,j,m)) &
                        +wn*(pvtrop(i,j,n)+pvtrop(i+1,j,n)))*.125

                  vbcors_t(i,j) = vbcors_t(i,j)+q

                  vtndcy = q+ &
                       (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl) &
                                         -xiym_o(i,j)*pb_t(i,j-1,nl))) &
                       +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,m)*pb_t(i,j-1,nl))) &
                       +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,n)*pb_t(i,j-1,nl)))) &
                       *scvyi(i,j)

                  vbflx_t(i,j,nl)= &
                       (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl) &
                      +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j) &
                                       *min(pb_t(i,j-1,nl),pb_t(i,j,nl)) &
                                       -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl) = max(-vminb(i,j),min(vmaxb(i,j), &
                                         vbflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else if (mommth == 'enecon'.or.mommth == 'enedis') then

            ! Sadourny (1975) energy conserving scheme

            !$omp parallel do private(l,i,q,vtndcy)
            do j = 0,jj+2
              do l = 1,isv(j)
                do i = max(0,ifv(j,l)),min(ii,ilv(j,l))

                  vbflxs_t(i,j) = vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl) &
                                 +(1.+wbaro)*vbflx_t(i,j,ml)

                  q = -.25*( (ubflx_t(i  ,j  ,nl)*scuyi(i  ,j  ) &
                             +ubflx_t(i  ,j-1,nl)*scuyi(i  ,j-1)) &
                            *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m) &
                             +wn*pvtrop(i,j,n)) &
                             +(ubflx_t(i+1,j  ,nl)*scuyi(i+1,j  ) &
                              +ubflx_t(i+1,j-1,nl)*scuyi(i+1,j-1)) &
                             *(wo*pvtrop_o(i+1,j)+wm*pvtrop(i+1,j,m) &
                              +wn*pvtrop(i+1,j,n)))

                  vbcors_t(i,j) = vbcors_t(i,j)+q

                  vtndcy = q+ &
                       (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl) &
                                         -xiym_o(i,j)*pb_t(i,j-1,nl))) &
                       +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,m)*pb_t(i,j-1,nl))) &
                       +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,n)*pb_t(i,j-1,nl)))) &
                       *scvyi(i,j)

                  vbflx_t(i,j,nl)= &
                       (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl) &
                      +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j) &
                                       *min(pb_t(i,j-1,nl),pb_t(i,j,nl)) &
                                       -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl) = max(-vminb(i,j),min(vmaxb(i,j), &
                                        vbflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else
            if (mnproc == 1) then
              write (lp,'(3a)') ' mommth = ',trim(mommth), &
                   ' is unsupported!'
            end if
            call xcstop('(barotp)')
            stop '(barotp)'
          end if


          ll = ml
          ml = nl
          nl = ll

        else

          wo = woa*lll+wob
          wn = wna*lll+wnb
          wm = 1.-wo-wn

          ! continuity equation

          !$omp parallel do private(l,i)
          do j = 0,jj+1
            do l = 1,isp(j)
              do i = max(0,ifp(j,l)),min(ii,ilp(j,l))
                pb_t(i,j,nl) = (1.-wbaro)*pb_t(i,j,ml)+wbaro*pb_t(i,j,nl) &
                              -(1.+wbaro)*dlt*(ubflx_t(i+1,j,ml)-ubflx_t(i,j,ml) &
                              +vbflx_t(i,j+1,ml)-vbflx_t(i,j,ml)) &
                              *scp2i(i,j)
              end do
            end do
          end do
          !$omp end parallel do

          ! v momentum equation

          if     (mommth == 'enscon') then

            ! Sadourny (1975) enstrophy conserving scheme

            !$omp parallel do private(l,i,q,vtndcy)
            do j = 1,jj+1
              do l = 1,isv(j)
                do i = max(0,ifv(j,l)),min(ii,ilv(j,l))

                  vbflxs_t(i,j) = vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl) &
                                 +(1.+wbaro)*vbflx_t(i,j,ml)

                  q = -(ubflx_t(i  ,j  ,ml)*scuyi(i  ,j  ) &
                       +ubflx_t(i+1,j  ,ml)*scuyi(i+1,j  ) &
                       +ubflx_t(i  ,j-1,ml)*scuyi(i  ,j-1) &
                       +ubflx_t(i+1,j-1,ml)*scuyi(i+1,j-1)) &
                       *(wo*(pvtrop_o(i,j)+pvtrop_o(i+1,j)) &
                        +wm*(pvtrop(i,j,m)+pvtrop(i+1,j,m)) &
                        +wn*(pvtrop(i,j,n)+pvtrop(i+1,j,n)))*.125

                  vbcors_t(i,j) = vbcors_t(i,j)+q

                  vtndcy = q+ &
                       (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl) &
                                         -xiym_o(i,j)*pb_t(i,j-1,nl))) &
                       +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,m)*pb_t(i,j-1,nl))) &
                       +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,n)*pb_t(i,j-1,nl)))) &
                       *scvyi(i,j)

                  vbflx_t(i,j,nl)= &
                       (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl) &
                       +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j) &
                                        *min(pb_t(i,j-1,nl),pb_t(i,j,nl)) &
                                        -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl) = max(-vminb(i,j),min(vmaxb(i,j), &
                                        vbflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else if (mommth == 'enecon'.or.mommth == 'enedis') then

            ! Sadourny (1975) energy conserving scheme

            !$omp parallel do private(l,i,q,vtndcy)
            do j = 1,jj+1
              do l = 1,isv(j)
                do i = max(0,ifv(j,l)),min(ii,ilv(j,l))

                  vbflxs_t(i,j) = vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl) &
                                 +(1.+wbaro)*vbflx_t(i,j,ml)

                  q = -.25*( (ubflx_t(i  ,j  ,ml)*scuyi(i  ,j  ) &
                             +ubflx_t(i  ,j-1,ml)*scuyi(i  ,j-1)) &
                             *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m) &
                              +wn*pvtrop(i,j,n)) &
                             +(ubflx_t(i+1,j  ,ml)*scuyi(i+1,j  ) &
                              +ubflx_t(i+1,j-1,ml)*scuyi(i+1,j-1)) &
                             *(wo*pvtrop_o(i+1,j)+wm*pvtrop(i+1,j,m) &
                              +wn*pvtrop(i+1,j,n)))

                  vbcors_t(i,j) = vbcors_t(i,j)+q

                  vtndcy = q+ &
                       (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl) &
                                         -xiym_o(i,j)*pb_t(i,j-1,nl))) &
                       +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,m)*pb_t(i,j-1,nl))) &
                       +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl) &
                                         -xiym(i,j,n)*pb_t(i,j-1,nl)))) &
                       *scvyi(i,j)

                  vbflx_t(i,j,nl)= &
                       (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl) &
                      +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j) &
                                       *min(pb_t(i,j-1,nl),pb_t(i,j,nl)) &
                                       -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl) = max(-vminb(i,j),min(vmaxb(i,j), &
                                        vbflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else
            if (mnproc == 1) then
              write (lp,'(3a)') ' mommth = ',trim(mommth),' is unsupported!'
            end if
            call xcstop('(barotp)')
            stop '(barotp)'
          end if

          ! u momentum equation

          if (mommth == 'enscon') then

            ! Sadourny (1975) enstrophy conserving scheme

            !$omp parallel do private(l,i,q,utndcy)
            do j = 1,jj
              do l = 1,isu(j)
                do i = max(1,ifu(j,l)),min(ii,ilu(j,l))

                  ubflxs_t(i,j) = ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl) &
                       +(1.+wbaro)*ubflx_t(i,j,ml)

                  q= (vbflx_t(i  ,j  ,nl)*scvxi(i  ,j  ) &
                     +vbflx_t(i  ,j+1,nl)*scvxi(i  ,j+1) &
                     +vbflx_t(i-1,j  ,nl)*scvxi(i-1,j  ) &
                     +vbflx_t(i-1,j+1,nl)*scvxi(i-1,j+1)) &
                     *(wo*(pvtrop_o(i,j)+pvtrop_o(i,j+1)) &
                      +wm*(pvtrop(i,j,m)+pvtrop(i,j+1,m)) &
                      +wn*(pvtrop(i,j,n)+pvtrop(i,j+1,n)))*.125

                  ubcors_t(i,j) = ubcors_t(i,j)+q

                  utndcy = q+ &
                       (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl) &
                                         -xixm_o(i,j)*pb_t(i-1,j,nl))) &
                       +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl) &
                                         -xixm(i,j,m)*pb_t(i-1,j,nl))) &
                       +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl) &
                                         -xixm(i,j,n)*pb_t(i-1,j,nl)))) &
                       *scuxi(i,j)

                  ubflx_t(i,j,nl)= &
                       (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl) &
                      +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j) &
                                       *min(pb_t(i-1,j,nl),pb_t(i,j,nl)) &
                                       -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl) = max(-uminb(i,j),min(umaxb(i,j), &
                                        ubflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else if (mommth == 'enecon'.or.mommth == 'enedis') then

            ! Sadourny (1975) energy conserving scheme

            !$omp parallel do private(l,i,q,utndcy)
            do j = 1,jj
              do l = 1,isu(j)
                do i = max(1,ifu(j,l)),min(ii,ilu(j,l))

                  ubflxs_t(i,j) = ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl) &
                                 +(1.+wbaro)*ubflx_t(i,j,ml)

                  q= .25*( (vbflx_t(i  ,j  ,nl)*scvxi(i  ,j  ) &
                           +vbflx_t(i-1,j  ,nl)*scvxi(i-1,j  )) &
                           *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m) &
                           +wn*pvtrop(i,j,n)) &
                           +(vbflx_t(i  ,j+1,nl)*scvxi(i  ,j+1) &
                            +vbflx_t(i-1,j+1,nl)*scvxi(i-1,j+1)) &
                           *(wo*pvtrop_o(i,j+1)+wm*pvtrop(i,j+1,m) &
                            +wn*pvtrop(i,j+1,n)))

                  ubcors_t(i,j) = ubcors_t(i,j)+q

                  utndcy = q+ &
                       (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl) &
                                         -xixm_o(i,j)*pb_t(i-1,j,nl))) &
                       +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl) &
                                         -xixm(i,j,m)*pb_t(i-1,j,nl))) &
                       +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl) &
                                         -xixm(i,j,n)*pb_t(i-1,j,nl)))) &
                       *scuxi(i,j)

                  ubflx_t(i,j,nl)= &
                       (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl) &
                      +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j) &
                       *min(pb_t(i-1,j,nl),pb_t(i,j,nl)) &
                           -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl) = max(-uminb(i,j),min(umaxb(i,j), &
                                        ubflx_t(i,j,nl)))

                end do
              end do
            end do
            !$omp end parallel do

          else
            if (mnproc == 1) then
              write (lp,'(3a)') ' mommth = ',trim(mommth),' is unsupported!'
            end if
            call xcstop('(barotp)')
            stop '(barotp)'
          end if

          ll = ml
          ml = nl
          nl = ll

        end if

      end do ! lll

      lll0 = lll0+lstep/2

      if     (nb == 1) then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb(i,j,m) = pb_t(i,j,ml)
            end do
          end do
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              pbu(i,j,m) = min(pb_t(i,j,ml),pb_t(i-1,j,ml))
              ubflx(i,j,m) = ubflx_t(i,j,ml)
              ub(i,j,m) = ubflx(i,j,m)/(pbu(i,j,m)*scuy(i,j))
              ubflxs(i,j,n) = ubflxs(i,j,n)+ubflxs_t(i,j)
              ubflxs(i,j,m) = ubflxs(i,j,3)+ubflxs_t(i,j)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              pbv(i,j,m) = min(pb_t(i,j,ml),pb_t(i,j-1,ml))
              vbflx(i,j,m) = vbflx_t(i,j,ml)
              vb(i,j,m) = vbflx(i,j,m)/(pbv(i,j,m)*scvx(i,j))
              vbflxs(i,j,n) = vbflxs(i,j,n)+vbflxs_t(i,j)
              vbflxs(i,j,m) = vbflxs(i,j,3)+vbflxs_t(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (nb == 2) then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb_mn(i,j,ml) = pb_t(i,j,ml)
              pb_mn(i,j,nl) = pb_t(i,j,nl)
            end do
          end do
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              ubflx_mn(i,j,ml) = ubflx_t(i,j,ml)
              ubflx_mn(i,j,nl) = ubflx_t(i,j,nl)
              ubflxs(i,j,m) = ubflxs(i,j,m)+ubflxs_t(i,j)
              ubflxs(i,j,3) = ubflxs_t(i,j)
              ubflxs_p(i,j,n) = ubflxs_t(i,j)
              ubcors_p(i,j) = ubcors_t(i,j)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vbflx_mn(i,j,ml) = vbflx_t(i,j,ml)
              vbflx_mn(i,j,nl) = vbflx_t(i,j,nl)
              vbflxs(i,j,m) = vbflxs(i,j,m)+vbflxs_t(i,j)
              vbflxs(i,j,3) = vbflxs_t(i,j)
              vbflxs_p(i,j,n) = vbflxs_t(i,j)
              vbcors_p(i,j) = vbcors_t(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (nb == 3) then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb(i,j,n) = pb_t(i,j,ml)
            end do
          end do
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              pbu(i,j,n) = min(pb_t(i,j,ml),pb_t(i-1,j,ml))
              ubflx(i,j,n) = ubflx_t(i,j,ml)
              ub(i,j,n) = ubflx(i,j,n)/(pbu(i,j,n)*scuy(i,j))
              ubflxs_p(i,j,m) = ubflxs(i,j,m)+ubflxs_t(i,j)
              ubflxs_p(i,j,n) = ubflxs_p(i,j,n)+ubflxs_t(i,j)
              ubcors_p(i,j) = ubcors_p(i,j)+ubcors_t(i,j)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              pbv(i,j,n) = min(pb_t(i,j,ml),pb_t(i,j-1,ml))
              vbflx(i,j,n) = vbflx_t(i,j,ml)
              vb(i,j,n) = vbflx(i,j,n)/(pbv(i,j,n)*scvx(i,j))
              vbflxs_p(i,j,m) = vbflxs(i,j,m)+vbflxs_t(i,j)
              vbflxs_p(i,j,n) = vbflxs_p(i,j,n)+vbflxs_t(i,j)
              vbcors_p(i,j) = vbcors_p(i,j)+vbcors_t(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (nb == 4) then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              ubflxs_p(i,j,n) = ubflxs_p(i,j,n)+ubflxs_t(i,j)
              ubcors_p(i,j) = ubcors_p(i,j)+ubcors_t(i,j)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vbflxs_p(i,j,n) = vbflxs_p(i,j,n)+vbflxs_t(i,j)
              vbcors_p(i,j) = vbcors_p(i,j)+vbcors_t(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (nb == 5) then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb_p(i,j) = pb_t(i,j,ml)
            end do
          end do
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              pbu_p(i,j) = min(pb_t(i,j,ml),pb_t(i-1,j,ml))
              ubflxs_p(i,j,n) = ubflxs_p(i,j,n)+ubflxs_t(i,j)
              ubcors_p(i,j) = ubcors_p(i,j)+ubcors_t(i,j)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              pbv_p(i,j) = min(pb_t(i,j,ml),pb_t(i,j-1,ml))
              vbflxs_p(i,j,n) = vbflxs_p(i,j,n)+vbflxs_t(i,j)
              vbcors_p(i,j) = vbcors_p(i,j)+vbcors_t(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      end if

    end do ! nb

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'barotp:'
      end if
      call chksummsk(pb,ip,2,'pb')
      call chksummsk(pbu,iu,2,'pbu')
      call chksummsk(ubflx,iu,2,'ubflx')
      call chksummsk(ub,iu,2,'ub')
      call chksummsk(ubflxs,iu,3,'ubflxs')
      call chksummsk(pbv,iv,2,'pbv')
      call chksummsk(vbflx,iv,2,'vbflx')
      call chksummsk(vb,iv,2,'vb')
      call chksummsk(vbflxs,iv,3,'vbflxs')
      call chksummsk(pb_p,ip,1,'pb_p')
      call chksummsk(pbu_p,iu,1,'pbu_p')
      call chksummsk(ubflxs_p,iu,2,'ubflxs_p')
      call chksummsk(ubcors_p,iu,1,'ubcors_p')
      call chksummsk(pbv_p,iv,1,'pbv_p')
      call chksummsk(vbflxs_p,iv,2,'vbflxs_p')
      call chksummsk(vbcors_p,iv,1,'vbcors_p')
    end if

  end subroutine barotp

end module mod_barotp
