! ------------------------------------------------------------------------------
! Copyright (C) 2005-2020 Mats Bentsen, Mehmet Ilicak
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

      module mod_pbcor
c
c --- ------------------------------------------------------------------
c --- This module contains variables and procedures related to the
c --- application of corrected layer-wise mass fluxes so that the
c --- vertical sum of layer thicknesses match the bottom pressure found
c --- by the integration of the barotropic equations.
c --- ------------------------------------------------------------------
c
      use mod_types, only: r8
      use mod_constants, only: epsilp, P_mks2cgs
      use mod_time, only: dlt
      use mod_xc
      use mod_grid, only: scp2i
      use mod_eos, only: sig
      use mod_state, only: dp, temp, saln, sigma, uflx, vflx,
     .                     utflx, vtflx, usflx, vsflx, p,
     .                     ubflxs, vbflxs, pb, ubflxs_p, vbflxs_p, pb_p
      use mod_utility, only: utotm, vtotm, utotn, vtotn, uflux, vflux,
     .                       uflux2, vflux2, uflux3, vflux3
      use mod_checksum, only: csdiag, chksummsk
#ifdef TRC
      use mod_tracers, only: ntr, trc, uflxtr, vflxtr
#endif
c
      implicit none
c
      private
c
c --- Variables to be set in namelist:
      character(len = 80) ::
     .  bmcmth ! Baroclinic mass flux correction method. Valid methods:
               ! 'uc' (upstream column), 'dluc' (depth limited upstream
               ! column).
c
c --- Parameters:
      real(r8), parameter ::
     .  dpeps1 = 1.e-5_r8*P_mks2cgs, ! Small layer pressure thickness
                                     ! [g cm-1 s-2].
     .  dpeps2 = 1.e-7_r8*P_mks2cgs  ! Small layer pressure thickness
                                     ! [g cm-1 s-2].
c
      public :: bmcmth, pbcor1, pbcor2
c
      contains
c
c --- ------------------------------------------------------------------
c
      subroutine pbcor1(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Correct the layer thicknesses to match the predictive bottom
c --- pressure from the barotropic solution.
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu_t,pbv_t
      real, dimension(1-nbdy:idm+nbdy) :: pbfac
      real dpo,dpni
      integer i,j,k,l,kn,km
#ifdef TRC
      integer nt
#endif
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=0,jj+1
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,k,km)
      do j=1,jj
        if     (bmcmth.eq.'uc') then
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j)=dlt*ubflxs_p(i,j,m)
          enddo
          enddo
        elseif (bmcmth.eq.'dluc') then
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j)=dlt*ubflxs_p(i,j,m)
            pbu_t(i,j)=min(p(i,j,kk+1),p(i-1,j,kk+1))
          enddo
          enddo
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' bmcmth=',trim(bmcmth),' is unsupported!'
          endif
          call xcstop('(pbcor1)')
                 stop '(pbcor1)'
        endif
        do k=1,kk
          km=k+mm
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j)=utotm(i,j)-uflx(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,k,km)
      do j=1,jj+1
        if     (bmcmth.eq.'uc') then
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j)=dlt*vbflxs_p(i,j,m)
          enddo
          enddo
        elseif (bmcmth.eq.'dluc') then
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j)=dlt*vbflxs_p(i,j,m)
            pbv_t(i,j)=min(p(i,j,kk+1),p(i,j-1,kk+1))
          enddo
          enddo
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' bmcmth=',trim(bmcmth),' is unsupported!'
          endif
          call xcstop('(pbcor1)')
                 stop '(pbcor1)'
        endif
        do k=1,kk
          km=k+mm
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j)=vtotm(i,j)-vflx(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        km=k+mm
        kn=k+nn
c
        if     (bmcmth.eq.'uc') then
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotm(i,j).gt.0.) then
                uflux(i,j)=utotm(i,j)*dp(i-1,j,kn)/p(i-1,j,kk+1)
                uflux2(i,j)=uflux(i,j)*saln(i-1,j,kn)
                uflux3(i,j)=uflux(i,j)*temp(i-1,j,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,kn,nt)
                enddo
#endif
              else
                uflux(i,j)=utotm(i,j)*dp(i  ,j,kn)/p(i  ,j,kk+1)
                uflux2(i,j)=uflux(i,j)*saln(i  ,j,kn)
                uflux3(i,j)=uflux(i,j)*temp(i  ,j,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,kn,nt)
                enddo
#endif
              endif
              uflx(i,j,km)=uflx(i,j,km)+uflux(i,j)
              usflx(i,j,km)=usflx(i,j,km)+uflux2(i,j)
              utflx(i,j,km)=utflx(i,j,km)+uflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj+1
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotm(i,j).gt.0.) then
                vflux(i,j)=vtotm(i,j)*dp(i,j-1,kn)/p(i,j-1,kk+1)
                vflux2(i,j)=vflux(i,j)*saln(i,j-1,kn)
                vflux3(i,j)=vflux(i,j)*temp(i,j-1,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,kn,nt)
                enddo
#endif
              else
                vflux(i,j)=vtotm(i,j)*dp(i,j  ,kn)/p(i,j  ,kk+1)
                vflux2(i,j)=vflux(i,j)*saln(i,j  ,kn)
                vflux3(i,j)=vflux(i,j)*temp(i,j  ,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,kn,nt)
                enddo
#endif
              endif
              vflx(i,j,km)=vflx(i,j,km)+vflux(i,j)
              vsflx(i,j,km)=vsflx(i,j,km)+vflux2(i,j)
              vtflx(i,j,km)=vtflx(i,j,km)+vflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (bmcmth.eq.'dluc') then
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotm(i,j).gt.0.) then
                uflux(i,j)=
     .            utotm(i,j)
     .            *max(0.,min(pbu_t(i,j),p(i-1,j,k+1))-p(i-1,j,k))
     .            /pbu_t(i,j)
                uflux2(i,j)=uflux(i,j)*saln(i-1,j,kn)
                uflux3(i,j)=uflux(i,j)*temp(i-1,j,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,kn,nt)
                enddo
#endif
              else
                uflux(i,j)=
     .            utotm(i,j)
     .            *max(0.,min(pbu_t(i,j),p(i  ,j,k+1))-p(i  ,j,k))
     .            /pbu_t(i,j)
                uflux2(i,j)=uflux(i,j)*saln(i  ,j,kn)
                uflux3(i,j)=uflux(i,j)*temp(i  ,j,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,kn,nt)
                enddo
#endif
              endif
              uflx(i,j,km)=uflx(i,j,km)+uflux(i,j)
              usflx(i,j,km)=usflx(i,j,km)+uflux2(i,j)
              utflx(i,j,km)=utflx(i,j,km)+uflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj+1
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotm(i,j).gt.0.) then
                vflux(i,j)=
     .            vtotm(i,j)
     .            *max(0.,min(pbv_t(i,j),p(i,j-1,k+1))-p(i,j-1,k))
     .            /pbv_t(i,j)
                vflux2(i,j)=vflux(i,j)*saln(i,j-1,kn)
                vflux3(i,j)=vflux(i,j)*temp(i,j-1,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,kn,nt)
                enddo
#endif
              else
                vflux(i,j)=
     .            vtotm(i,j)
     .            *max(0.,min(pbv_t(i,j),p(i,j  ,k+1))-p(i,j  ,k))
     .            /pbv_t(i,j)
                vflux2(i,j)=vflux(i,j)*saln(i,j  ,kn)
                vflux3(i,j)=vflux(i,j)*temp(i,j  ,kn)
#ifdef TRC
                do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
                  if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,kn,nt)
                enddo
#endif
              endif
              vflx(i,j,km)=vflx(i,j,km)+vflux(i,j)
              vsflx(i,j,km)=vsflx(i,j,km)+vflux2(i,j)
              vtflx(i,j,km)=vtflx(i,j,km)+vflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' bmcmth=',trim(bmcmth),' is unsupported!'
          endif
          call xcstop('(pbcor1)')
                 stop '(pbcor1)'
        endif
c
c$OMP PARALLEL DO PRIVATE(l,i,dpo,dpni
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpo=dp(i,j,kn)
            dp(i,j,kn)=max(0.,dpo-(uflux(i+1,j)-uflux(i,j)
     .                            +vflux(i,j+1)-vflux(i,j))*scp2i(i,j))
            dpo=dpo+dpeps1
            dpni=1./(dp(i,j,kn)+dpeps1)
            saln(i,j,kn)=(dpo*saln(i,j,kn)
     .                   -(uflux2(i+1,j)-uflux2(i,j)
     .                    +vflux2(i,j+1)-vflux2(i,j))*scp2i(i,j))*dpni
            temp(i,j,kn)=(dpo*temp(i,j,kn)
     .                   -(uflux3(i+1,j)-uflux3(i,j)
     .                    +vflux3(i,j+1)-vflux3(i,j))*scp2i(i,j))*dpni
#ifdef TRC
            do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
              if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
              trc(i,j,kn,nt)=(dpo*trc(i,j,kn,nt)
     .          -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .           +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
            enddo
#endif
            if (dp(i,j,kn).lt.dpeps2) dp(i,j,kn)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i,pbfac)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfac(i)=pb_p(i,j)/p(i,j,kk+1)
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'pbcor1:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(uflx,iu,2*kk,'uflx')
        call chksummsk(vflx,iv,2*kk,'vflx')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      end subroutine pbcor1
c
c --- ------------------------------------------------------------------
c
      subroutine pbcor2(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Correct the layer thicknesses to better match the corrected bottom
c --- pressure from the barotropic solution.
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu_t,pbv_t
      real, dimension(1-nbdy:idm+nbdy) :: pbfac
      real dpo,dpni
      integer i,j,k,l,kn,km
#ifdef TRC
      integer nt
#endif
c
      call xctilr(ubflxs(1-nbdy,1-nbdy,n), 1,1, 1,1, halo_uv)
      call xctilr(vbflxs(1-nbdy,1-nbdy,n), 1,1, 1,1, halo_vv)
#ifdef TRC
      do nt=1,ntr
        call xctilr(trc(1-nbdy,1-nbdy,k1m,nt), 1,kk, 1,1, halo_ps)
      enddo
#endif
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
      do j=0,jj+1
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,km)=max(0.,dp(i,j,km))+epsilp
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,k,kn)
      do j=1,jj
        if     (bmcmth.eq.'uc') then
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j)=dlt*ubflxs(i,j,n)
          enddo
          enddo
        elseif (bmcmth.eq.'dluc') then
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j)=dlt*ubflxs(i,j,n)
            pbu_t(i,j)=min(p(i,j,kk+1),p(i-1,j,kk+1))
          enddo
          enddo
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' bmcmth=',trim(bmcmth),' is unsupported!'
          endif
          call xcstop('(pbcor2)')
                 stop '(pbcor2)'
        endif
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j)=utotn(i,j)-uflx(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,k,kn)
      do j=1,jj+1
        if     (bmcmth.eq.'uc') then
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j)=dlt*vbflxs(i,j,n)
          enddo
          enddo
        elseif (bmcmth.eq.'dluc') then
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j)=dlt*vbflxs(i,j,n)
            pbv_t(i,j)=min(p(i,j,kk+1),p(i,j-1,kk+1))
          enddo
          enddo
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' bmcmth=',trim(bmcmth),' is unsupported!'
          endif
          call xcstop('(pbcor2)')
                 stop '(pbcor2)'
        endif
        do k=1,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j)=vtotn(i,j)-vflx(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
        km=k+mm
c
        if     (bmcmth.eq.'uc') then
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotn(i,j).gt.0.) then
                uflux(i,j)=utotn(i,j)*dp(i-1,j,km)/p(i-1,j,kk+1)
                uflux2(i,j)=uflux(i,j)*saln(i-1,j,km)
                uflux3(i,j)=uflux(i,j)*temp(i-1,j,km)
#ifdef TRC
                do nt=1,ntr
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,km,nt)
                enddo
#endif
              else
                uflux(i,j)=utotn(i,j)*dp(i  ,j,km)/p(i  ,j,kk+1)
                uflux2(i,j)=uflux(i,j)*saln(i  ,j,km)
                uflux3(i,j)=uflux(i,j)*temp(i  ,j,km)
#ifdef TRC
                do nt=1,ntr
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,km,nt)
                enddo
#endif
              endif
              uflx(i,j,kn)=uflx(i,j,kn)+uflux(i,j)
              usflx(i,j,kn)=usflx(i,j,kn)+uflux2(i,j)
              utflx(i,j,kn)=utflx(i,j,kn)+uflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj+1
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotn(i,j).gt.0.) then
                vflux(i,j)=vtotn(i,j)*dp(i,j-1,km)/p(i,j-1,kk+1)
                vflux2(i,j)=vflux(i,j)*saln(i,j-1,km)
                vflux3(i,j)=vflux(i,j)*temp(i,j-1,km)
#ifdef TRC
                do nt=1,ntr
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,km,nt)
                enddo
#endif
              else
                vflux(i,j)=vtotn(i,j)*dp(i,j  ,km)/p(i,j  ,kk+1)
                vflux2(i,j)=vflux(i,j)*saln(i,j  ,km)
                vflux3(i,j)=vflux(i,j)*temp(i,j  ,km)
#ifdef TRC
                do nt=1,ntr
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,km,nt)
                enddo
#endif
              endif
              vflx(i,j,kn)=vflx(i,j,kn)+vflux(i,j)
              vsflx(i,j,kn)=vsflx(i,j,kn)+vflux2(i,j)
              vtflx(i,j,kn)=vtflx(i,j,kn)+vflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (bmcmth.eq.'dluc') then
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              if (utotn(i,j).gt.0.) then
                uflux(i,j)=
     .            utotn(i,j)
     .            *max(0.,min(pbu_t(i,j),p(i-1,j,k+1))-p(i-1,j,k))
     .            /pbu_t(i,j)
                uflux2(i,j)=uflux(i,j)*saln(i-1,j,km)
                uflux3(i,j)=uflux(i,j)*temp(i-1,j,km)
#ifdef TRC
                do nt=1,ntr
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,km,nt)
                enddo
#endif
              else
                uflux(i,j)=
     .            utotn(i,j)
     .            *max(0.,min(pbu_t(i,j),p(i  ,j,k+1))-p(i  ,j,k))
     .            /pbu_t(i,j)
                uflux2(i,j)=uflux(i,j)*saln(i  ,j,km)
                uflux3(i,j)=uflux(i,j)*temp(i  ,j,km)
#ifdef TRC
                do nt=1,ntr
                  uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,km,nt)
                enddo
#endif
              endif
              uflx(i,j,kn)=uflx(i,j,kn)+uflux(i,j)
              usflx(i,j,kn)=usflx(i,j,kn)+uflux2(i,j)
              utflx(i,j,kn)=utflx(i,j,kn)+uflux3(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
          do j=1,jj+1
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (vtotn(i,j).gt.0.) then
                vflux(i,j)=
     .            vtotn(i,j)
     .            *max(0.,min(pbv_t(i,j),p(i,j-1,k+1))-p(i,j-1,k))
     .            /pbv_t(i,j)
                vflux2(i,j)=vflux(i,j)*saln(i,j-1,km)
                vflux3(i,j)=vflux(i,j)*temp(i,j-1,km)
#ifdef TRC
                do nt=1,ntr
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,km,nt)
                enddo
#endif
              else
                vflux(i,j)=
     .            vtotn(i,j)
     .            *max(0.,min(pbv_t(i,j),p(i,j  ,k+1))-p(i,j  ,k))
     .            /pbv_t(i,j)
                vflux2(i,j)=vflux(i,j)*saln(i,j  ,km)
                vflux3(i,j)=vflux(i,j)*temp(i,j  ,km)
#ifdef TRC
                do nt=1,ntr
                  vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,km,nt)
                enddo
#endif
              endif
              vflx(i,j,kn)=vflx(i,j,kn)+vflux(i,j)
              vsflx(i,j,kn)=vsflx(i,j,kn)+vflux2(i,j)
              vtflx(i,j,kn)=vtflx(i,j,kn)+vflux3(i,j)
            enddo
            enddo
          enddo
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' bmcmth=',trim(bmcmth),' is unsupported!'
          endif
          call xcstop('(pbcor2)')
                 stop '(pbcor2)'
        endif
c
c$OMP PARALLEL DO PRIVATE(l,i,dpo,dpni
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpo=dp(i,j,km)
            dp(i,j,km)=dpo-scp2i(i,j)*(uflux(i+1,j)-uflux(i,j)
     .                                +vflux(i,j+1)-vflux(i,j))
            dpni=1./dp(i,j,km)
            saln(i,j,km)=(dpo*saln(i,j,km)
     .                   -scp2i(i,j)*(uflux2(i+1,j)-uflux2(i,j)
     .                               +vflux2(i,j+1)-vflux2(i,j)))*dpni
            temp(i,j,km)=(dpo*temp(i,j,km)
     .                   -scp2i(i,j)*(uflux3(i+1,j)-uflux3(i,j)
     .                               +vflux3(i,j+1)-vflux3(i,j)))*dpni
#ifdef TRC
            do nt=1,ntr
              trc(i,j,km,nt)=(dpo*trc(i,j,km,nt)
     .          -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .           +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
            enddo
#endif
            sigma(i,j,km)=sig(temp(i,j,km),saln(i,j,km))
            dp(i,j,km)=dp(i,j,km)-epsilp
            if (dp(i,j,km).lt.dpeps2) dp(i,j,km)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i,pbfac)
      do j=1,jj
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfac(i)=pb(i,j,m)/p(i,j,kk+1)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,km)=dp(i,j,km)*pbfac(i)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'pbcor2:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(p,ip,kk+1,'p')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(uflx,iu,2*kk,'uflx')
        call chksummsk(vflx,iv,2*kk,'vflx')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      end subroutine pbcor2
c
      end module mod_pbcor
