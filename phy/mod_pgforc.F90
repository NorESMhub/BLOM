! ------------------------------------------------------------------------------
! Copyright (C) 2005-2022 Mats Bentsen, Mehmet Ilicak

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_pgforc

  ! --- ------------------------------------------------------------------
  ! --- This module contains variables and procedures related to time
  ! --- integration of the baroclinic momentum equation.
  ! --- ------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm
  use mod_xc,        only: xctilr, nbdy, ii, jj, kk, &
                           isp, ifp, ilp, isu, ifu, ilu, isv, ifv, ilv, &
                           halo_ps, mnproc, lp, ip, iu, iv
  use mod_types,     only: r8
  use mod_constants, only: g, epsilp, spval
  use mod_state,     only: dp, dpu, dpv, temp, saln, p, pu, pv, phi, &
                           pb_p, pbu_p, pbv_p, sealv
  use mod_eos,       only: delphi
  use mod_checksum,  only: csdiag, chksummsk

  implicit none
  private

  ! --- Constants.
  real(r8) :: &
       wpgf = .25_r8 ! Weight for time averaging of pressure gradient
                     ! force [].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) :: &
       pgfx  ! x-component of baroclinic pressure gradient
             ! force [cm2 s-2]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) :: &
       pgfy  ! y-component of baroclinic pressure gradient
             ! force [cm2 s-2].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       pgfxo         ! 'pgfx' at old time level [cm2 s-2]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       pgfyo         ! 'pgfy' at old time level [cm2 s-2].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       pgfxm  ! x-component of barotropic pressure gradient
              ! force, not dependent on bottom pressure
              ! [cm2 s-2]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       pgfym  ! y-component of barotropic pressure gradient
              ! force, not dependent on bottom pressure
              ! [cm2 s-2]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       xixp   ! Dependeny of x-component of barotropic pressure
              ! gradient force on bottom pressure at (i, j)
              ! [g cm s-4]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       xixm   ! Dependeny of x-component of barotropic pressure
              ! gradient force on bottom pressure at (i - 1, j)
              ! [g cm s-4]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       xiyp   ! Dependeny of y-component of barotropic pressure
              ! gradient force on bottom pressure at (i, j)
              ! [g cm s-4]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       xiym   ! Dependeny of y-component of barotropic pressure
              ! gradient force on bottom pressure at (i, j - 1)
              ! [g cm s-4].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       pgfxm_o       ! 'pgfxm' at old time level [cm2 s-2]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       pgfym_o       ! 'pgfym' at old time level [cm2 s-2]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       xixp_o        ! 'xixp' at old time level [g cm s-4]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       xixm_o        ! 'xixm' at old time level [g cm s-4]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       xiyp_o        ! 'xiyp' at old time level [g cm s-4]. &
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       xiym_o        ! 'xiym' at old time level [g cm s-4].

  ! Public variables
  public :: wpgf, pgfx, pgfy, pgfxo, pgfyo, pgfxm, pgfym
  public :: xixp, xixm, xiyp, xiym, pgfxm_o, pgfym_o
  public :: xixp_o, xixm_o, xiyp_o, xiym_o

  ! Public routines
  public :: inivar_pgforc, pgforc

contains

  ! --- ------------------------------------------------------------------

  subroutine inivar_pgforc()

    ! --- ------------------------------------------------------------------
    ! --- Initialize arrays.
    ! --- ------------------------------------------------------------------

    ! Local variables
    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,k)
    do j = 1-nbdy,jj+nbdy
      do k = 1,2*kk
        do i = 1-nbdy,ii+nbdy
          pgfx(i,j,k) = spval
          pgfy(i,j,k) = spval
        end do
      end do
      do k = 1,kk
        do i = 1-nbdy,ii+nbdy
          pgfxo(i,j,k) = spval
          pgfyo(i,j,k) = spval
        end do
      end do
      do k = 1,2
        do i = 1-nbdy,ii+nbdy
          pgfxm(i,j,k) = spval
          pgfym(i,j,k) = spval
          xixp(i,j,k) = spval
          xixm(i,j,k) = spval
          xiyp(i,j,k) = spval
          xiym(i,j,k) = spval
        end do
      end do
      do i = 1-nbdy,ii+nbdy
        pgfxm_o(i,j) = spval
        pgfym_o(i,j) = spval
        xixp_o(i,j) = spval
        xixm_o(i,j) = spval
        xiyp_o(i,j) = spval
        xiym_o(i,j) = spval
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine inivar_pgforc

  ! --- ------------------------------------------------------------------

  subroutine pgforc(m,n,mm,nn,k1m,k1n)

    ! --- ------------------------------------------------------------------
    ! --- compute the pressure gradient force
    ! --- ------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: phip
    real :: dphi,alpl,alpu,prs,dphip,dphim,alplp,alpup,alplm,alpum,cp,cm, &
         phi_p,phi_m,q
    integer :: kup(idm),kum(idm),kvp(idm),kvm(idm)
    integer :: i,j,k,l,kn

    ! --- compute new -dpu,dpv- field.

    !$OMP PARALLEL DO PRIVATE(k,kn,l,i)
    do j = -2,jj+2
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(k,kn,l,i,q)
    do j = -1,jj+2
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q = min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)= &
                 .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                 +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
            pu(i,j,k+1) = pu(i,j,k)+dpu(i,j,kn)
          end do
        end do
        do l = 1,isv(j)
          do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q = min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)= &
                 .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                 +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
            pv(i,j,k+1) = pv(i,j,k)+dpv(i,j,kn)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(l,i,k,kn,dphi,alpu,alpl)
    do j = 0,jj
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii,ilp(j,l))
          phip(i,j,kk+1) = 0.
        end do
      end do
      do k = kk,1,-1
        kn = k+nn
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii,ilp(j,l))
            if (dp(i,j,kn) < epsilp) then
              phi (i,j,k) = phi (i,j,k+1)
              phip(i,j,k) = phip(i,j,k+1)
            else
              call delphi(p(i,j,k),p(i,j,k+1),temp(i,j,kn),saln(i,j,kn), &
                   dphi,alpu,alpl)
              phi (i,j,k) = phi (i,j,k+1)-dphi
              phip(i,j,k) = phip(i,j,k+1)+p(i,j,k+1)*alpl-p(i,j,k)*alpu
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !$omp parallel do private(l,i)
    do j = -1,jj+2
      do l = 1,isu(j)
        do i = max(0,ifu(j,l)),min(ii+1,ilu(j,l))
          xixp_o(i,j) = xixp(i,j,n)
          xixm_o(i,j) = xixm(i,j,n)
          pgfxm_o(i,j) = pgfxm(i,j,n)
        end do
      end do
      do l = 1,isv(j)
        do i = max(0,ifv(j,l)),min(ii+1,ilv(j,l))
          xiyp_o(i,j) = xiyp(i,j,n)
          xiym_o(i,j) = xiym(i,j,n)
          pgfym_o(i,j) = pgfym(i,j,n)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private( &
    !$omp l,i,kup,kum,kvp,kvm,k,kn,prs,dphip,alpup,alplp,dphim,alpum,alplm, &
    !$omp cp,cm,q,phi_p,phi_m)
    do j = 1,jj

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          kup(i) = kk
          kum(i) = kk
          xixp(i,j,n) = 0.
          xixm(i,j,n) = 0.
          pgfxm(i,j,n) = 0.
        end do
      end do

      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          kvp(i) = kk
          kvm(i) = kk
          xiyp(i,j,n) = 0.
          xiym(i,j,n) = 0.
          pgfym(i,j,n) = 0.
        end do
      end do

      do k = kk,1,-1
        kn = k+nn

        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))

            prs = pu(i,j,k+1)-.5*dpu(i,j,kn)

            do while (p(i  ,j,kup(i)) > prs)
              kup(i) = kup(i)-1
            end do

            do while (p(i-1,j,kum(i)) > prs)
              kum(i) = kum(i)-1
            end do

            call delphi(prs,p(i  ,j,kup(i)+1), &
                 temp(i  ,j,kup(i)+nn),saln(i  ,j,kup(i)+nn), &
                 dphip,alpup,alplp)

            call delphi(prs,p(i-1,j,kum(i)+1), &
                 temp(i-1,j,kum(i)+nn),saln(i-1,j,kum(i)+nn), &
                 dphim,alpum,alplm)

            cp = .25*(p(i  ,j,k+1)+p(i  ,j,k))
            cm = .25*(p(i-1,j,k+1)+p(i-1,j,k))
            q = prs/(cp+cm)
            !           if (i.eq.itest.and.j.eq.jtest) write (lp,*) 'u',k,q
            cp = q*cp
            cm = q*cm

            phi_p = phi(i  ,j,kup(i)+1)-dphip
            xixp(i,j,n) = xixp(i,j,n) &
                 +(phip(i  ,j,kup(i)+1) &
                 +p(i  ,j,kup(i)+1)*alplp-cp*(alpup-alpum)) &
                 *dpu(i,j,kn)

            phi_m = phi(i-1,j,kum(i)+1)-dphim
            xixm(i,j,n) = xixm(i,j,n) &
                 +(phip(i-1,j,kum(i)+1) &
                 +p(i-1,j,kum(i)+1)*alplm-cm*(alpum-alpup)) &
                 *dpu(i,j,kn)

            pgfxo(i,j,k) = pgfx(i,j,kn)
            pgfx(i,j,kn) = -(phi_p-phi_m)
            pgfxm(i,j,n) = pgfxm(i,j,n)+pgfx(i,j,kn)*dpu(i,j,kn)

          end do
        end do

        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))

            prs = pv(i,j,k+1)-.5*dpv(i,j,kn)

            do while (p(i,j  ,kvp(i)) > prs)
              kvp(i) = kvp(i)-1
            end do

            do while (p(i,j-1,kvm(i)) > prs)
              kvm(i) = kvm(i)-1
            end do

            call delphi(prs,p(i,j  ,kvp(i)+1), &
                 temp(i,j  ,kvp(i)+nn),saln(i,j  ,kvp(i)+nn), &
                 dphip,alpup,alplp)

            call delphi(prs,p(i,j-1,kvm(i)+1), &
                 temp(i,j-1,kvm(i)+nn),saln(i,j-1,kvm(i)+nn), &
                 dphim,alpum,alplm)

            cp = .25*(p(i,j  ,k+1)+p(i,j  ,k))
            cm = .25*(p(i,j-1,k+1)+p(i,j-1,k))
            q = prs/(cp+cm)
            !           if (i.eq.itest.and.j.eq.jtest) write (lp,*) 'v',k,q
            cp = q*cp
            cm = q*cm

            phi_p = phi(i,j  ,kvp(i)+1)-dphip
            xiyp(i,j,n) = xiyp(i,j,n) &
                 +(phip(i,j  ,kvp(i)+1) &
                 +p(i,j  ,kvp(i)+1)*alplp-cp*(alpup-alpum)) &
                 *dpv(i,j,kn)

            phi_m = phi(i,j-1,kvm(i)+1)-dphim
            xiym(i,j,n) = xiym(i,j,n) &
                 +(phip(i,j-1,kvm(i)+1) &
                 +p(i,j-1,kvm(i)+1)*alplm-cm*(alpum-alpup)) &
                 *dpv(i,j,kn)

            pgfyo(i,j,k) = pgfy(i,j,kn)
            pgfy(i,j,kn) = -(phi_p-phi_m)
            pgfym(i,j,n) = pgfym(i,j,n)+pgfy(i,j,kn)*dpv(i,j,kn)

          end do
        end do

      end do

    end do
    !$omp end parallel do

    call xctilr(pb_p, 1,1, 1,1, halo_ps)

    !$omp parallel do private(l,i,q,k,kn)
    do j = 1,jj

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          q = 1./pbu_p(i,j)
          pgfxm(i,j,n) = pgfxm(i,j,n)*q
          xixp(i,j,n) = xixp(i,j,n)*q
          xixm(i,j,n) = xixm(i,j,n)*q
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          q = 1./pbv_p(i,j)
          pgfym(i,j,n) = pgfym(i,j,n)*q
          xiyp(i,j,n) = xiyp(i,j,n)*q
          xiym(i,j,n) = xiym(i,j,n)*q
        end do
      end do

      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            pgfx(i,j,kn) = pgfx(i,j,kn)-pgfxm(i,j,n)
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            pgfy(i,j,kn) = pgfy(i,j,kn)-pgfym(i,j,n)
          end do
        end do
      end do

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          pgfxm(i,j,n) = pgfxm(i,j,n)+xixp(i,j,n)-xixm(i,j,n)
          xixp(i,j,n) = xixp(i,j,n)/pb_p(i  ,j)
          xixm(i,j,n) = xixm(i,j,n)/pb_p(i-1,j)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          pgfym(i,j,n) = pgfym(i,j,n)+xiyp(i,j,n)-xiym(i,j,n)
          xiyp(i,j,n) = xiyp(i,j,n)/pb_p(i,j  )
          xiym(i,j,n) = xiym(i,j,n)/pb_p(i,j-1)
        end do
      end do

      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          sealv(i,j) = phi(i,j,1)/g
        end do
      end do

    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'pgforc:'
      end if
      call chksummsk(phi,ip,kk+1,'phi')
      call chksummsk(pgfx,iu,2*kk,'pgfx')
      call chksummsk(pgfy,iv,2*kk,'pgfy')
      call chksummsk(pgfxm,iu,2,'pgfxm')
      call chksummsk(pgfym,iv,2,'pgfym')
      call chksummsk(xixp,iu,2,'xixp')
      call chksummsk(xixm,iu,2,'xixm')
      call chksummsk(xiyp,iv,2,'xiyp')
      call chksummsk(xiym,iv,2,'xiym')
    end if

  end subroutine pgforc

end module mod_pgforc
