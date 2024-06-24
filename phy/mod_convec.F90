! ------------------------------------------------------------------------------
! Copyright (C) 2009-2024 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_convec

  use dimensions,    only: idm, jdm, kdm
  use mod_constants, only: epsilp
  use mod_xc,        only: xctilr, ii, jj, kk, isp, ifp, ilp, &
                           isu, ifu, ilu, ifv, ilv, ip, iu, &
                           lp, i0, j0, isv, iv, mnproc, nbdy, halo_ps
  use mod_vcoord,    only: sigmar
  use mod_eos,       only: rho, sig, sofsig
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                           p, pu, pv, kfpla
  use mod_checksum,  only: csdiag, chksummsk
  use mod_tracers,   only: ntr, trc
  use mod_ifdefs,    only: use_TRC

  implicit none
  private

  public :: convec

contains

  subroutine convec(m,n,mm,nn,k1m,k1n)

    ! ------------------------------------------------------------------
    ! Remove static instabilitites between the mixed layer and interior
    ! layers
    ! ------------------------------------------------------------------

    ! Arguments
    integer :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(kdm+1) :: po,pn
    real, dimension(kdm) :: ttem,ssal,delp,dens,densr,uo,un
    real    :: tdps,sdps,dps,ttmp,stmp,dtmp,q,udpn
    integer :: i,j,k,l,kn,kfpl,kfplo,kmix,ko
    logical :: done
    integer :: niter
    real, dimension(ntr,kdm) :: ttrc
    real, dimension(ntr) :: trdps
    integer :: nt

    !$omp parallel do private( &
    !$omp l,i,k,kn,ttem,ssal,delp,dens,densr,dps,kfpl,kfplo,tdps,sdps,q, &
    !$omp ttmp,stmp,dtmp,done,niter,kmix, &
    !$omp nt,ttrc,trdps)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! Copy variables into 1d arrays
          do k = 1,kk
            kn = k+nn
            ttem(k) = temp(i,j,kn)
            ssal(k) = saln(i,j,kn)
            delp(k) = dp(i,j,kn)
            dens(k) = sigma(i,j,kn)
            densr(k) = sigmar(i,j,k)
            if (use_TRC) then
              do nt = 1,ntr
                ttrc(nt,k) = trc(i,j,kn,nt)
              end do
            end if
          end do

          !---------------------------------------------------------------
          ! Define first physical interior layer
          !---------------------------------------------------------------

          k = 3
          dps = 0.
          do while (delp(k) < epsilp)
            dps = dps+delp(k)
            delp(k) = 0.
            k = k+1
            if (k > kk) exit
          end do
          if (k > kk) then
            delp(2) = delp(2)+dps
          else
            delp(k) = delp(k)+dps
          end if
          kfpl = k
          kfplo = kfpla(i,j,n)
          if (kfpl < kfplo) then
            if (kfplo <= kk) then
              tdps = 0.
              sdps = 0.
              dps = 0.
              if (use_TRC) then
                do nt = 1,ntr
                  trdps(nt) = 0.
                end do
              end if
              do k = kfpl,kfplo
                tdps = tdps+ttem(k)*delp(k)
                sdps = sdps+ssal(k)*delp(k)
                dps = dps+delp(k)
                if (use_TRC) then
                  do nt = 1,ntr
                    trdps(nt) = trdps(nt)+ttrc(nt,k)*delp(k)
                  end do
                end if
              end do
              q = 1./dps
              ttmp = tdps*q
              stmp = sdps*q
              dtmp = sig(ttmp,stmp)
              if (dtmp > densr(kfplo)) then
                do k = kfpl,kfplo-1
                  delp(k) = 0.
                end do
                kfpl = kfplo
                ttem(kfpl) = ttmp
                ssal(kfpl) = stmp
                dens(kfpl) = dtmp
                delp(kfpl) = dps
                if (use_TRC) then
                  do nt = 1,ntr
                    ttrc(nt,kfpl) = trdps(nt)*q
                  end do
                end if
              end if
            else
              tdps = 0.
              sdps = 0.
              dps = 0.
              if (use_TRC) then
                do nt = 1,ntr
                  trdps(nt) = 0.
                end do
              end if
              do k = kfpl,kk
                tdps = tdps+ttem(k)*delp(k)
                sdps = sdps+ssal(k)*delp(k)
                dps = dps+delp(k)
                if (use_TRC) then
                  do nt = 1,ntr
                    trdps(nt) = trdps(nt)+ttrc(nt,k)*delp(k)
                  end do
                end if
                delp(k) = 0.
              end do
              q = 1./dps
              ttmp = tdps*q
              stmp = sdps*q
              dtmp = sig(ttmp,stmp)
              kfpl = kk
              do while (dtmp < densr(kfpl))
                if (kfpl == 3) exit
                kfpl = kfpl-1
              end do
              ttem(kfpl) = ttmp
              ssal(kfpl) = stmp
              dens(kfpl) = dtmp
              delp(kfpl) = dps
              if (use_TRC) then
                do nt = 1,ntr
                  ttrc(nt,kfpl) = trdps(nt)*q
                end do
              end if
            end if
          end if

          if (kfpl <= kk) then

            !---------------------------------------------------------------
            ! Remove static instabilities
            !---------------------------------------------------------------

            done = .false.

            niter = 0
            do while (.not.done)
              niter = niter+1
              if (niter == 100) then
                write (lp,*) 'blom: convec: no convergence!',i+i0,j+j0
                exit
              end if

              done = .true.

              ! Remove instabilities between the lower mixed layer and
              ! interior layers by considering the potential density jump
              ! across the mixed layer base with reference pressure at the
              ! interface
              tdps = ttem(2)*delp(2)
              sdps = ssal(2)*delp(2)
              dps = delp(2)
              if (use_TRC) then
                do nt = 1,ntr
                  trdps(nt) = ttrc(nt,2)*delp(2)
                end do
              end if
              ttmp = ttem(2)
              stmp = ssal(2)
              k = kfpl
              do while (rho(dps,ttmp,stmp) > &
                        rho(dps,ttem(k),ssal(k)).or. &
                        delp(k) < epsilp)
                tdps = tdps+ttem(k)*delp(k)
                sdps = sdps+ssal(k)*delp(k)
                dps = dps+delp(k)
                q = 1./dps
                ttmp = tdps*q
                stmp = sdps*q
                if (use_TRC) then
                  do nt = 1,ntr
                    trdps(nt) = trdps(nt)+ttrc(nt,k)*delp(k)
                  end do
                end if
                k = k+1
                if (k > kk) exit
              end do
              kmix = k-1
              if (kmix >= kfpl) then
                ttem(2) = ttmp
                ssal(2) = stmp
                dens(2) = sig(ttem(2),ssal(2))
                if (use_TRC) then
                  do nt = 1,ntr
                    ttrc(nt,2) = trdps(nt)*q
                  end do
                end if
                dps = 0.
                do k = kfpl,kmix
                  dps = dps+delp(k)
                  delp(k) = 0.
                end do
                k = kmix
                do while (dens(2) < densr(k))
                  if (k == 3) exit
                  k = k-1
                end do
                kfpl = k
                ttem(kfpl) = ttem(2)
                ssal(kfpl) = ssal(2)
                dens(kfpl) = dens(2)
                delp(kfpl) = dps
                if (use_TRC) then
                  do nt = 1,ntr
                    ttrc(nt,kfpl) = ttrc(nt,2)
                  end do
                end if
                do k = kfpl+1,kmix
                  ttem(k) = ttem(2)
                  dens(k) = densr(k)
                  ssal(k) = sofsig(dens(k),ttem(k))
                end do
              end if

            end do

          end if

          kfpla(i,j,n) = kfpl

          ! Copy 1d arrays to 3d arrays
          do k = 1,kk
            kn = k+nn
            temp(i,j,kn) = ttem(k)
            saln(i,j,kn) = ssal(k)
            sigma(i,j,kn) = dens(k)
            dp(i,j,kn) = delp(k)
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
            if (use_TRC) then
              do nt = 1,ntr
                trc(i,j,kn,nt) = ttrc(nt,k)
              end do
            end if
          end do

        end do
      end do
    end do
    !$omp end parallel do

    !---------------------------------------------------------------
    ! Redistribute momentum
    !---------------------------------------------------------------

    call xctilr(p, 1,kk+1, 1,1, halo_ps)

    !$omp parallel do private(l,i,k,kn,uo,po,pn,ko,un,udpn)
    do j = 1,jj

      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k = 1,kk
            kn = k+nn
            uo(k) = u(i,j,kn)
          end do
          po(1) = 0.
          pn(1) = 0.
          do k = 2,kk+1
            po(k) = pu(i,j,k)
            pn(k) = .5*(min(pu(i,j,kk+1),p(i  ,j,k)) &
                       +min(pu(i,j,kk+1),p(i-1,j,k)))
          end do

          ko = 1
          do kn = 1,kk
            if (pn(kn+1)-pn(kn) == 0.) then
              un(kn) = 0.
            else
              udpn = 0.
              do while (pn(kn+1) > po(ko+1))
                udpn = udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko = ko+1
              end do
              un(kn) = (udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn)))) &
                      /(pn(kn+1)-pn(kn))
            end if
          end do
          do k = 1,kk
            kn = k+nn
            u(i,j,kn) = un(k)
          end do

        end do
      end do

      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          do k = 1,kk
            kn = k+nn
            uo(k) = v(i,j,kn)
          end do
          po(1) = 0.
          pn(1) = 0.
          do k = 2,kk+1
            po(k) = pv(i,j,k)
            pn(k) = .5*(min(pv(i,j,kk+1),p(i,j  ,k)) &
                       +min(pv(i,j,kk+1),p(i,j-1,k)))
          end do

          ko = 1
          do kn = 1,kk
            if (pn(kn+1)-pn(kn) == 0.) then
              un(kn) = 0.
            else
              udpn = 0.
              do while (pn(kn+1) > po(ko+1))
                udpn = udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko = ko+1
              end do
              un(kn) = (udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn)))) &
                      /(pn(kn+1)-pn(kn))
            end if
          end do
          do k = 1,kk
            kn = k+nn
            v(i,j,kn) = un(k)
          end do

        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,kn,l,i,q)
    do j = 1,jj
      do k = 1,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            q = min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)= &
                 .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                    +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            q = min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)= &
                 .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                    +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          end do
        end do
      end do
    end do
    !$omp end parallel do
    !     do j=1,jj
    !       do l=1,isu(j)
    !       do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
    !         q=0.
    !         do k=1,kk
    !           kn=k+nn
    !           q=q+u(i,j,kn)*dpu(i,j,kn)
    !         enddo
    !         if (abs(q).gt.1.e-4) then
    !           write (lp,*) 'convec: u imbalance:',q,i,j
    !         endif
    !       enddo
    !       enddo
    !       do l=1,isv(j)
    !       do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
    !         q=0.
    !         do k=1,kk
    !           kn=k+nn
    !           q=q+v(i,j,kn)*dpv(i,j,kn)
    !         enddo
    !         if (abs(q).gt.1.e-4) then
    !           write (lp,*) 'convec: v imbalance:',q,i,j
    !         endif
    !       enddo
    !       enddo
    !     enddo

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'convec:'
      end if
      call chksummsk(dp,ip,2*kk,'dp')
      call chksummsk(temp,ip,2*kk,'temp')
      call chksummsk(saln,ip,2*kk,'saln')
      call chksummsk(sigma,ip,2*kk,'sigma')
      call chksummsk(u,iu,2*kk,'u')
      call chksummsk(v,iv,2*kk,'v')
      if (use_TRC) then
        do nt = 1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        end do
      end if
    end if

  end subroutine convec

end module mod_convec
