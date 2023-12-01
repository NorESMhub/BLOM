! ------------------------------------------------------------------------------
! Copyright (C) 2005-2022 Mats Bentsen, Mehmet Ilicak
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
c
c --- ------------------------------------------------------------------
c --- This module contains variables and procedures related to time
c --- smoothing. Note that time smoothing of baroclinic velocity is
c --- carried out in the baroclinic momentum equation solver (routine
c --- momtum of mod_momtum.F).
c --- ------------------------------------------------------------------
c
      use mod_types, only: r8
      use mod_constants, only: epsilp, spval
      use mod_xc
      use mod_vcoord, only: vcoord_type_tag, isopyc_bulkml
      use mod_state, only: dp, dpu, dpv, temp, saln, p, pb
      use mod_checksum, only: csdiag, chksummsk
#ifdef TRC
      use mod_tracers, only: ntr, trc, trcold
#endif
c
      implicit none
c
      private
c
c --- Weights for time smoothing.
      real(r8) ::
     .  wuv1 = .75_r8, 
     .  wuv2 = .125_r8, 
     .  wts1 = .875_r8, 
     .  wts2 = .0625_r8,
     .  wbaro = .125_r8
c
      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .  dpold          ! Layer pressure thickness at old time level
                       ! [g cm-1 s-2].
c
      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  dpuold,        ! Layer pressure thickness at u-point at old time
                       ! level [g cm-1 s-2].
     .  dpvold,        ! Layer pressure thickness at v-point at old time
                       ! level [g cm-1 s-2].
     .  told,          ! Potential temperature at old time level
                       ! [deg C].
     .  sold           ! Salinity at old time level [g kg-1].
c
      public :: wuv1, wuv2, wts1, wts2, wbaro, dpold, dpuold, dpvold,
     .          inivar_tmsmt, initms, tmsmt1, tmsmt2
c
      contains
c
c --- ------------------------------------------------------------------
c
      subroutine inivar_tmsmt
c
c --- ------------------------------------------------------------------
c --- Initialize arrays.
c --- ------------------------------------------------------------------
c
      integer :: i,j,k,l
c
c$OMP PARALLEL DO PRIVATE(i,k)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          do k=1,kk
            dpold(i,j,k   )=spval
            dpold(i,j,k+kk)=spval
            dpuold(i,j,k)=spval
            dpvold(i,j,k)=spval
            told(i,j,k)=spval
            sold(i,j,k)=spval
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- initialize  dpuold  upstream and downstream of p-points as well as
c --- at lateral neighbors of interior u-points.
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=0,jj+1
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k=1,kk
            dpuold(i,j-1,k)=0.
            dpuold(i,j+1,k)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
          do k=1,kk
            dpuold(i,j,k)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- initialize  dpvold  upstream and downstream of p-points as well as
c --- at lateral neighbors of interior v-points.
c
c$OMP PARALLEL DO PRIVATE(l,j,k)
      do i=0,ii+1
        do l=1,jsv(i)
        do j=max(1,jfv(i,l)),min(jj,jlv(i,l))
          do k=1,kk
            dpvold(i-1,j,k)=0.
            dpvold(i+1,j,k)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,j,k)
      do i=1,ii
        do l=1,jsp(i)
        do j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
          do k=1,kk
            dpvold(i,j,k)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(dpuold, 1,  kk, nbdy,nbdy, halo_us)
      call xctilr(dpvold, 1,  kk, nbdy,nbdy, halo_vs)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'inivar_tmsmt:'
        endif
        call chksummsk(dpuold,iu,kk,'dpuold')
        call chksummsk(dpvold,iv,kk,'dpvold')
      endif
c
      end subroutine inivar_tmsmt
c
c --- ------------------------------------------------------------------
c
      subroutine initms(m,n,mm,nn,k1m,k1n)
c
c --- save old layer thickness, temperature and salinity for time
c --- smoothing
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
#ifdef TRC
      integer nt
#endif
c
      integer i,j,k,l,km
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1,jj
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpold(i,j,km)=dp(i,j,km)
            told(i,j,k)=temp(i,j,km)
            sold(i,j,k)=saln(i,j,km)
#ifdef TRC
            do nt=1,ntr
              trcold(i,j,k,nt)=trc(i,j,km,nt)
            enddo
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'initms:'
        endif
        call chksummsk(dpold,ip,2*kk,'dpold')
        call chksummsk(told,ip,kk,'told')
        call chksummsk(sold,ip,kk,'sold')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trcold(1-nbdy,1-nbdy,1,nt),ip,kk,'trcold')
        enddo
#endif
      endif
c
      end subroutine initms
c
      subroutine tmsmt1(m,n,mm,nn,k1m,k1n)
c
c --- save old layer thickness at velocity points for time smoothing in
c --- momentum equation.
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
      integer i,j,k,l,kn
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            dpuold(i,j,k)=dpu(i,j,kn)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            dpvold(i,j,k)=dpv(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'tmsmt1:'
        endif
        call chksummsk(dpuold,iu,kk,'dpuold')
        call chksummsk(dpvold,iv,kk,'dpvold')
      endif
c
      end subroutine tmsmt1
c
c --- ------------------------------------------------------------------
c
      subroutine tmsmt2(m,n,mm,nn,k1m,k1n)
c
c --- time smoothing of layer thickness, temperature and salinity
c
      use mod_constants, only: epsilp
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy) :: pbfaco,pbfacn
      integer i,j,k,l,kn,km,kp
      real pold,pmid,pnew,q
#ifdef TRC
      integer nt
#endif
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,pbfaco,pbfacn,k,kn,km,kp,pold,pmid,pnew,q
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i)=0.
          pbfacn(i)=0.
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pbfaco(i)=pbfaco(i)+dpold(i,j,kn)
            pbfacn(i)=pbfacn(i)+dp(i,j,kn)
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i)=pb(i,j,m)/pbfaco(i)
          pbfacn(i)=pb(i,j,m)/pbfacn(i)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          kn=k+nn
          kp=min(k+1,kk)
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pold=max(0.,dpold(i,j,kn)*pbfaco(i))
            pmid=max(0.,dp(i,j,km))
            pnew=max(0.,dp(i,j,kn)*pbfacn(i))
            dp(i,j,km)=wts1*pmid+wts2*(pold+pnew)
            dpold(i,j,km)=dp(i,j,km)
            pold=pold+epsilp
            pmid=pmid+epsilp
            pnew=pnew+epsilp
            temp(i,j,km)=(wts1*pmid*temp(i,j,km)
     .                   +wts2*(pold*told(i,j,k)+pnew*temp(i,j,kn)))
     .                   /(dp(i,j,km)+epsilp)
            told(i,j,k)=temp(i,j,km)
            saln(i,j,km)=(wts1*pmid*saln(i,j,km)
     .                   +wts2*(pold*sold(i,j,k)+pnew*saln(i,j,kn)))
     .                   /(dp(i,j,km)+epsilp)
            sold(i,j,k)=saln(i,j,km)
#ifdef TRC
            do nt=1,ntr
              trc(i,j,km,nt)=(wts1*pmid*trc(i,j,km,nt)
     .                       +wts2*(pold*trcold(i,j,k,nt)
     .                             +pnew*trc(i,j,kn,nt)))
     .                       /(dp(i,j,km)+epsilp)
              trcold(i,j,k,nt)=trc(i,j,km,nt)
            enddo
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (vcoord_type_tag == isopyc_bulkml) then
c
        call xctilr(dp(1-nbdy,1-nbdy,k1m), 1,kk, 3,3, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=-2,jj+2
          do k=1,kk
            do l=1,isp(j)
            do i=max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i,q)
        do j=-1,jj+2
          do k=1,kk
            km=k+mm
            do l=1,isu(j)
            do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
              q=min(p(i,j,kk+1),p(i-1,j,kk+1))
              dpu(i,j,km)=
     .          .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k)))
     .             +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
            enddo
            enddo
            do l=1,isv(j)
            do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              q=min(p(i,j,kk+1),p(i,j-1,kk+1))
              dpv(i,j,km)=
     .          .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k)))
     .             +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      else
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'tmsmt2:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(dpu,iu,2*kk,'dpu')
        call chksummsk(dpv,iv,2*kk,'dpv')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      end subroutine tmsmt2
c
      end module mod_tmsmt
