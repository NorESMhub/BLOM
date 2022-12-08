! ------------------------------------------------------------------------------
! Copyright (C) 2007-2020 Mats Bentsen, Mehmet Ilicak
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

      module mod_advect
c
c --- ------------------------------------------------------------------
c --- This module contains variables and procedures related to advection
c --- of layer pressure thickness and tracers by calling incremental
c --- remapping routines.
c --- ------------------------------------------------------------------
c
      use mod_types, only: r8
      use mod_constants, only: onemm
      use mod_time, only: delt1, dlt
      use mod_xc
      use mod_grid, only: scuy, scvx, scp2i, scp2
      use mod_state, only: u, v, dp, dpu, dpv, temp, saln, sigma,
     .                     uflx, vflx, utflx, vtflx, usflx, vsflx,
     .                     p, pbu, pbv, ubflxs_p, vbflxs_p
      use mod_diffusion, only: umfltd, vmfltd
      use mod_remap, only: remap_eitvel, remap_eitflx
      use mod_utility, only: utotm, vtotm, umax, vmax
      use mod_checksum, only: csdiag, chksummsk
#ifdef TRC
      use mod_tracers, only: ntr, itrtke, itrgls, trc, uflxtr, vflxtr
#endif
c
      implicit none
c
      private
c
c --- Variables to be set in namelist:
      character(len = 80) ::
     .  rmpmth ! Method of applying eddy-induced transport in the remap
               ! transport algorithm. Valid methods: 'eitvel', 'eitflx'.
c
      public :: rmpmth, advect
c
      contains
c
c --- ------------------------------------------------------------------
c
      subroutine advect(m,n,mm,nn,k1m,k1n)
c
      integer m,n,mm,nn,k1m,k1n
c
      integer i,j,k,l,km,kn,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
#ifdef TRC
      integer nt
#endif
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  l,i,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne)
      do j=-1,jj+2
        do l=1,isp(j)
        do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
          iw=i-iu(i  ,j)
          ie=i+iu(i+1,j)
          js=j-iv(i,j  )
          jn=j+iv(i,j+1)
          isw=i*(1-ip(iw,js))+iw*ip(iw,js)
          jsw=j*(1-ip(iw,js))+js*ip(iw,js)
          ise=i*(1-ip(ie,js))+ie*ip(ie,js)
          jse=j*(1-ip(ie,js))+js*ip(ie,js)
          inw=i*(1-ip(iw,jn))+iw*ip(iw,jn)
          jnw=j*(1-ip(iw,jn))+jn*ip(iw,jn)
          ine=i*(1-ip(ie,jn))+ie*ip(ie,jn)
          jne=j*(1-ip(ie,jn))+jn*ip(ie,jn)
          pbmin(i,j)=
     .      min(p(isw,jsw,kk+1),p(i  ,js ,kk+1),p(ise,jse,kk+1),
     .          p(iw ,j  ,kk+1),p(i  ,j  ,kk+1),p(ie ,j  ,kk+1),
     .          p(inw,jnw,kk+1),p(i  ,jn ,kk+1),p(ine,jne,kk+1))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
#ifdef TRC
      do nt=1,ntr
#  if defined(TKE) && !defined(TKEADV)
        if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
        call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1,kk, 3,3, halo_ps)
      enddo
#endif
c
      if     (rmpmth.eq.'eitvel') then
c$OMP PARALLEL DO PRIVATE(km,kn,j,l,i)
c$OMP+FIRSTPRIVATE(utotm,vtotm)
        do k=1,kk
          km=k+mm
          kn=k+nn
c
c --- --- advective and diffusive velocity at mid time level
c
          do j=-1,jj+2
            do l=1,isu(j)
            do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
              utotm(i,j)=u(i,j,km)
     .          +(ubflxs_p(i,j,m)*dlt/pbu(i,j,m)
     .           +umfltd(i,j,km)/max(onemm,dpu(i,j,kn)))
     .           /(delt1*scuy(i,j))
              utotm(i,j)=max(-umax(i,j),min(umax(i,j),utotm(i,j)))
            enddo
            enddo
          enddo
          do j=0,jj+2
            do l=1,isv(j)
            do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              vtotm(i,j)=v(i,j,km)
     .          +(vbflxs_p(i,j,m)*dlt/pbv(i,j,m)
     .           +vmfltd(i,j,km)/max(onemm,dpv(i,j,kn)))
     .           /(delt1*scvx(i,j))
              vtotm(i,j)=max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
            enddo
            enddo
          enddo
c
          call remap_eitvel(scuy,scvx,scp2i,scp2,pbmin,
     .                      pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n),
     .                      p(1-nbdy,1-nbdy,k+1),utotm,vtotm,delt1,1,
     .                      dp(1-nbdy,1-nbdy,kn),
     .                      temp(1-nbdy,1-nbdy,kn),
     .                      saln(1-nbdy,1-nbdy,kn),
     .                      uflx(1-nbdy,1-nbdy,km),
     .                      vflx(1-nbdy,1-nbdy,km),
     .                      utflx(1-nbdy,1-nbdy,km),
     .                      vtflx(1-nbdy,1-nbdy,km),
     .                      usflx(1-nbdy,1-nbdy,km),
     .                      vsflx(1-nbdy,1-nbdy,km)
#ifdef TRC
     .                     ,kn,trc
#endif
     .                     )
c
        enddo
c$OMP END PARALLEL DO
      elseif (rmpmth.eq.'eitflx') then
c$OMP PARALLEL DO PRIVATE(km,kn,j,l,i)
c$OMP+FIRSTPRIVATE(utotm,vtotm)
        do k=1,kk
          km=k+mm
          kn=k+nn
c
c --- --- advective velocity at mid time level
c
          do j=-1,jj+2
            do l=1,isu(j)
            do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
              utotm(i,j)=u(i,j,km)
     .                  +dlt*ubflxs_p(i,j,m)
     .                   /(delt1*pbu(i,j,m)*scuy(i,j))
              utotm(i,j)=max(-umax(i,j),min(umax(i,j),utotm(i,j)))
            enddo
            enddo
          enddo
          do j=0,jj+2
            do l=1,isv(j)
            do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
              vtotm(i,j)=v(i,j,km)
     .                  +dlt*vbflxs_p(i,j,m)
     .                   /(delt1*pbv(i,j,m)*scvx(i,j))
              vtotm(i,j)=max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
            enddo
            enddo
          enddo
c
          call remap_eitflx(scuy,scvx,scp2i,scp2,pbmin,
     .                      pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n),
     .                      p(1-nbdy,1-nbdy,k+1),utotm,vtotm,
     .                      umfltd(1-nbdy,1-nbdy,km),
     .                      vmfltd(1-nbdy,1-nbdy,km),
     .                      delt1,1,
     .                      dp(1-nbdy,1-nbdy,kn),
     .                      temp(1-nbdy,1-nbdy,kn),
     .                      saln(1-nbdy,1-nbdy,kn),
     .                      uflx(1-nbdy,1-nbdy,km),
     .                      vflx(1-nbdy,1-nbdy,km),
     .                      utflx(1-nbdy,1-nbdy,km),
     .                      vtflx(1-nbdy,1-nbdy,km),
     .                      usflx(1-nbdy,1-nbdy,km),
     .                      vsflx(1-nbdy,1-nbdy,km)
#ifdef TRC
     .                     ,kn,trc
#endif
     .                     )
c
        enddo
c$OMP END PARALLEL DO
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' rmpmth=',trim(rmpmth),' is unsupported!'
        endif
        call xcstop('(advect)')
               stop '(advect)'
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'advect:'
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
      end subroutine advect
c
      end module mod_advect
