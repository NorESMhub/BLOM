! ------------------------------------------------------------------------------
! Copyright (C) 2005-2025 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_tmsmt

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to time
  ! smoothing. Note that time smoothing of baroclinic velocity is
  ! carried out in the baroclinic momentum equation solver (routine
  ! momtum of mod_momtum.F90).
  ! ------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm
  use mod_xc,        only: xctilr, nbdy, ii, jj, kk, isu, ifu, ilu, &
                           isp, ifp, ilp, jsv, jfv, jsp, jfp, &
                           jlv, jlp, halo_ps, halo_us, halo_vs, mnproc, &
                           lp, iu, iv, ip, isv, ifv, isv, ilv, nbdy
  use mod_types,     only: r8
  use mod_constants, only: epsilp, spval
  use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml
  use mod_state,     only: dp, dpu, dpv, temp, saln, p, pb
  use mod_checksum,  only: csdiag, chksum
  use mod_tracers,   only: ntr, trc, trcold
  use mod_ifdefs,    only: use_TRC

  implicit none
  private

  ! Weights for time smoothing.
  real(r8) ::   &
       wuv1 = .75_r8,   &
       wuv2 = .125_r8,  &
       wts1 = .875_r8,  &
       wts2 = .0625_r8, &
       wbaro = .125_r8

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) :: &
       dpold   ! Layer pressure thickness at old time level
               ! [kg m-1 s-2].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       dpuold  ! Layer pressure thickness at u-point at old time
               ! level [kg m-1 s-2].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       dpvold  ! Layer pressure thickness at v-point at old time
               ! level [kg m-1 s-2].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       told    ! Potential temperature at old time level
               ! [deg C].
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       sold    ! Salinity at old time level [g kg-1].

  ! Public module variables
  public :: wuv1, wuv2, wts1, wts2, wbaro, dpold, dpuold, dpvold

  ! Public routines
  public :: inivar_tmsmt, initms, tmsmt1, tmsmt2

contains

  subroutine inivar_tmsmt()

    ! ------------------------------------------------------------------
    ! Initialize arrays.
    ! ------------------------------------------------------------------

    integer :: i,j,k,l

    dpold(:,:,:) = spval
    dpuold(:,:,:) = spval
    dpvold(:,:,:) = spval
    told(:,:,:) = spval
    sold(:,:,:) = spval

    ! initialize  dpuold  upstream and downstream of p-points as well as
    ! at lateral neighbors of interior u-points.

    !$omp parallel do private(l,i,k)
    do j = 0,jj+1
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k = 1,kk
            dpuold(i,j-1,k) = 0.
            dpuold(i,j+1,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,k)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
          do k = 1,kk
            dpuold(i,j,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! initialize  dpvold  upstream and downstream of p-points as well as
    ! at lateral neighbors of interior v-points.

    !$omp parallel do private(l,j,k)
    do i = 0,ii+1
      do l = 1,jsv(i)
        do j = max(1,jfv(i,l)),min(jj,jlv(i,l))
          do k = 1,kk
            dpvold(i-1,j,k) = 0.
            dpvold(i+1,j,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,j,k)
    do i = 1,ii
      do l = 1,jsp(i)
        do j = max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
          do k = 1,kk
            dpvold(i,j,k) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(dpuold, 1,  kk, nbdy,nbdy, halo_us)
    call xctilr(dpvold, 1,  kk, nbdy,nbdy, halo_vs)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'inivar_tmsmt:'
      end if
      call chksum(dpuold, kk, halo_us, 'dpuold')
      call chksum(dpvold, kk, halo_vs, 'dpvold')
    end if

  end subroutine inivar_tmsmt

  ! --- ------------------------------------------------------------------

  subroutine initms(mm)

    ! save old layer thickness, temperature and salinity for time
    ! smoothing

    ! Arguments
    integer, intent(in)  :: mm

    ! Local variables
    integer :: i,j,k,l,km,nt
    character(len = 2) cnt

    !$omp parallel do private(k,km,l,i,nt)
    do j = 1,jj
      do k = 1,kk
        km = k+mm
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            told(i,j,k) = temp(i,j,km)
            sold(i,j,k) = saln(i,j,km)
            if (use_TRC) then
              do nt = 1,ntr
                trcold(i,j,k,nt) = trc(i,j,km,nt)
              end do
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'initms:'
      end if
      call chksum(dpold, 2*kk, halo_ps, 'dpold')
      call chksum(told , kk  , halo_ps, 'told' )
      call chksum(sold , kk  , halo_ps, 'sold' )
      do nt = 1,ntr
        write(cnt, '(i2.2)') nt
        call chksum(trcold(1-nbdy,1-nbdy,1,nt), kk, halo_ps, 'trcold'//cnt)
      end do
    end if

  end subroutine initms

  ! --- ------------------------------------------------------------------

  subroutine tmsmt1(nn)

    ! save old layer thickness at velocity points for time smoothing in
    ! momentum equation.

    ! Arguments
    integer, intent(in) :: nn

    ! Local variables
    integer :: i,j,k,l,kn,nt
    character(len = 2) cnt

    !$omp parallel do private(k,kn,l,i,nt)
    do j = 1,jj
      do k = 1,kk
        kn = k+nn
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpold(i,j,kn)=dp(i,j,kn)
            told(i,j,k)=temp(i,j,kn)
            sold(i,j,k)=saln(i,j,kn)
            if (use_TRC) then
              do nt=1,ntr
                trcold(i,j,k,nt)=trc(i,j,kn,nt)
              enddo
            end if
          enddo
        enddo
      enddo
    end do
    !$omp end parallel do
    if (vcoord_tag == vcoord_isopyc_bulkml) then
      !$omp parallel do private(k,kn,l,i)
      do j = 1,jj
        do k = 1,kk
          kn = k+nn
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              dpuold(i,j,k) = dpu(i,j,kn)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              dpvold(i,j,k) = dpv(i,j,kn)
            end do
          end do
        end do
      end do
      !$omp end parallel do
    endif

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'tmsmt1:'
      end if
      call chksum(dpold , 2*kk, halo_ps, 'dpold' )
      call chksum(told  , kk  , halo_ps, 'told'  )
      call chksum(sold  , kk  , halo_ps, 'sold'  )
      do nt = 1,ntr
        write(cnt, '(i2.2)') nt
        call chksum(trcold(1-nbdy,1-nbdy,1,nt), kk, halo_ps, 'trcold'//cnt)
      end do
      if (vcoord_tag == vcoord_isopyc_bulkml) then
        call chksum(dpuold, kk  , halo_us, 'dpuold')
        call chksum(dpvold, kk  , halo_vs, 'dpvold')
      end if
    end if

  end subroutine tmsmt1

  ! --- ------------------------------------------------------------------

  subroutine tmsmt2(m,mm,nn,k1m)

    ! time smoothing of layer thickness, temperature and salinity

    ! Arguments
    integer, intent(in) :: m,mm,nn,k1m

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy) :: pbfaco,pbfacn
    integer :: i,j,k,l,kn,km,kp
    real :: pold,pmid,pnew,q
    integer :: nt
    character(len = 2) cnt

    !$omp parallel do private( &
    !$omp l,i,pbfaco,pbfacn,k,kn,km,kp,pold,pmid,pnew,q,nt)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i) = 0.
          pbfacn(i) = 0.
        end do
      end do
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            pbfaco(i) = pbfaco(i)+dpold(i,j,kn)
            pbfacn(i) = pbfacn(i)+dp(i,j,kn)
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i) = pb(i,j,m)/pbfaco(i)
          pbfacn(i) = pb(i,j,m)/pbfacn(i)
        end do
      end do
      do k = 1,kk
        km = k+mm
        kn = k+nn
        kp = min(k+1,kk)
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            pold = max(0.,dpold(i,j,kn)*pbfaco(i))
            pmid = max(0.,dp(i,j,km))
            pnew = max(0.,dp(i,j,kn)*pbfacn(i))
            dp(i,j,km) = wts1*pmid+wts2*(pold+pnew)
            pold = pold+epsilp
            pmid = pmid+epsilp
            pnew = pnew+epsilp
            temp(i,j,km) = (wts1*pmid*temp(i,j,km) &
                           +wts2*(pold*told(i,j,k)+pnew*temp(i,j,kn))) &
                           /(dp(i,j,km)+epsilp)
            saln(i,j,km) = (wts1*pmid*saln(i,j,km) &
                           +wts2*(pold*sold(i,j,k)+pnew*saln(i,j,kn))) &
                           /(dp(i,j,km)+epsilp)
            if (use_TRC) then
              do nt = 1,ntr
                trc(i,j,km,nt) = (wts1*pmid*trc(i,j,km,nt) &
                                 +wts2*(pold*trcold(i,j,k,nt) &
                                 +pnew*trc(i,j,kn,nt))) &
                                 /(dp(i,j,km)+epsilp)
              end do
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

    call xctilr(dp(1-nbdy,1-nbdy,k1m), 1,kk, 3,3, halo_ps)

    !$omp parallel do private(k,km,l,i)
    do j = -2,jj+2
      do k = 1,kk
        km = k+mm
        do l = 1,isp(j)
          do i = max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,km)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (vcoord_tag == vcoord_isopyc_bulkml) then

      !$omp parallel do private(k,km,l,i,q)
      do j = -1,jj+2
        do k = 1,kk
          km = k+mm
          do l = 1,isu(j)
            do i = max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
              q = min(p(i,j,kk+1),p(i-1,j,kk+1))
              dpu(i,j,km)= &
                   .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                      +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
            end do
          end do
          do l = 1,isv(j)
            do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              q = min(p(i,j,kk+1),p(i,j-1,kk+1))
              dpv(i,j,km)= &
                   .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                      +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
            end do
          end do
        end do
      end do
      !$omp end parallel do

    end if

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'tmsmt2:'
      end if
      call chksum(dp  , 2*kk, halo_ps, 'dp'  )
      call chksum(temp, 2*kk, halo_ps, 'temp')
      call chksum(saln, 2*kk, halo_ps, 'saln')
      call chksum(dpu , 2*kk, halo_us, 'dpu' )
      call chksum(dpv , 2*kk, halo_vs, 'dpv' )
      do nt = 1,ntr
        write(cnt, '(i2.2)') nt
        call chksum(trc(1-nbdy,1-nbdy,1,nt), 2*kk, halo_ps, 'trc'//cnt)
      end do
    end if

  end subroutine tmsmt2

end module mod_tmsmt
