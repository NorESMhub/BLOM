! ------------------------------------------------------------------------------
! Copyright (C) 2005-2018 Mats Bentsen, Mehmet Ilicak
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

      subroutine pbtcor(m,n,mm,nn,k1m,k1n)
c
c --- correct the layer thicknesses to better match the bottom pressure
c --- from the barotropic solution
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu_t,pbv_t
      real, dimension(1-nbdy:idm+nbdy) :: pbfac
      real dpo,dpni
      integer i,j,k,l,kn,km
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real sig
      external sig
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
            dp(i,j,km)=max(0.,dp(i,j,km))+epsil
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
          call xcstop('(pbtcor)')
                 stop '(pbtcor)'
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
          call xcstop('(pbtcor)')
                 stop '(pbtcor)'
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
          call xcstop('(pbtcor)')
                 stop '(pbtcor)'
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
            dp(i,j,km)=dp(i,j,km)-epsil
            if (dp(i,j,km).lt.1.e-6) dp(i,j,km)=0.
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
          write (lp,*) 'pbtcor:'
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
      return
      end
