! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2020 Mats Bentsen, Lars Inge Enstad
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

      subroutine barotp(m,n,mm,nn,k1m,k1n)
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  pb_t=spval,ubflx_t=spval,vbflx_t=spval
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ubflxs_t,vbflxs_t,ubcors_t,vbcors_t
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  umaxb=spval,uminb=spval,vmaxb=spval,vminb=spval,
     .  uglue=spval,vglue=spval
      real q,woa,wob,wna,wnb,wo,wm,wn,utndcy,vtndcy
      integer i,j,k,l,kn,nb,lll0,lll,ml,nl,ll
c
c --- ------------------------------------------------------------------
c --- initialize barotropic velocity sums, determine maximum allowable
c --- barotropic velocities, and determine coefficients for coastal wave
c --- breaking parameterization
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,k,kn)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          umaxb(i,j)=0.
          uminb(i,j)=0.
          uglue(i,j)=cwbdts*exp(1.-pbu(i,j,m)/(cwbdls*onem))
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            umaxb(i,j)=max(umaxb(i,j),u(i,j,kn))
            uminb(i,j)=min(uminb(i,j),u(i,j,kn))
          enddo
          enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          umaxb(i,j)=(umax(i,j)-umaxb(i,j))*pbu(i,j,m)*scuy(i,j)
          uminb(i,j)=(umax(i,j)+uminb(i,j))*pbu(i,j,m)*scuy(i,j)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vmaxb(i,j)=0.
          vminb(i,j)=0.
          vglue(i,j)=cwbdts*exp(1.-pbv(i,j,m)/(cwbdls*onem))
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vmaxb(i,j)=max(vmaxb(i,j),v(i,j,kn))
            vminb(i,j)=min(vminb(i,j),v(i,j,kn))
          enddo
          enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vmaxb(i,j)=(vmax(i,j)-vmaxb(i,j))*pbv(i,j,m)*scvx(i,j)
          vminb(i,j)=(vmax(i,j)+vminb(i,j))*pbv(i,j,m)*scvx(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- potential vorticity of barotropic flow
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=-2,jj+3
        do i=0,ii+1
          pvtrop_o(i,j)=pvtrop(i,j,n)
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q)
      do j=0,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          q=2./(pb_p(i,j)+pb_p(i-1,j))
          pvtrop(i,j  ,n)=corioq(i,j  )*q
          pvtrop(i,j+1,n)=corioq(i,j+1)*q
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q)
      do j=1,jj
        do l=1,isv(j)
        do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
          q=2./(pb_p(i,j)+pb_p(i,j-1))
          pvtrop(i  ,j,n)=corioq(i  ,j)*q
          pvtrop(i+1,j,n)=corioq(i+1,j)*q
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isq(j)
        do i=max(1,ifq(j,l)),min(ii,ilq(j,l))
          pvtrop(i,j,n)=corioq(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  )
     .                                 +pb_p(i,j-1)+pb_p(i-1,j-1))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
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
c
c --- with arctic patch, switch umaxb and uminb, vmaxb and vminb, xixp
c --- and xixm, and xiyp and xiym in the halo region adjacent to the
c --- arctic grid intersection
      if (nreg.eq.2.and.nproc.eq.jpr) then
        do j=jj,jj+2
          do i=0,ii+1
            q=umaxb(i,j)
            umaxb(i,j)=uminb(i,j)
            uminb(i,j)=q
            q=xixp(i,j,n)
            xixp(i,j,n)=xixm(i,j,n)
            xixm(i,j,n)=q
          enddo
        enddo
        do i=max(0,itdm/2-i0+1),ii+1
          q=vmaxb(i,jj)
          vmaxb(i,jj)=vminb(i,jj)
          vminb(i,jj)=q
          q=xiyp(i,jj,n)
          xiyp(i,jj,n)=xiym(i,jj,n)
          xiym(i,jj,n)=q
        enddo
        do j=jj+1,jj+2
          do i=0,ii+1
            q=vmaxb(i,j)
            vmaxb(i,j)=vminb(i,j)
            vminb(i,j)=q
            q=xiyp(i,j,n)
            xiyp(i,j,n)=xiym(i,j,n)
            xiym(i,j,n)=q
          enddo
        enddo
      endif
c
c --- ------------------------------------------------------------------
c --- advance barotropic equations from baroclinic time level -m- to
c --- level -n- then advance barotropic equations another baroclinic
c --- time level so that the average barotropic transport for a
c --- leap-frog baroclinic step can be predicted
c --- ------------------------------------------------------------------
c
      do nb=1,5
c
        if     (nb.eq.1) then
          lll0=1
          ml=1
          nl=2
          woa=-1./lstep
          wob=.5+(lll0-.5)/lstep
          wna=0.
          wnb=0.
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,jj
            do i=1,ii
              pb_t(i,j,ml)=pb_mn(i,j,ml)
              pb_t(i,j,nl)=pb_mn(i,j,nl)
              ubflx_t(i,j,ml)=ubflx_mn(i,j,ml)
              ubflx_t(i,j,nl)=ubflx_mn(i,j,nl)
              vbflx_t(i,j,ml)=vbflx_mn(i,j,ml)
              vbflx_t(i,j,nl)=vbflx_mn(i,j,nl)
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (nb.eq.2) then
          woa=0.
          wob=0.
          wna=1./lstep
          wnb=-(lll0-.5)/lstep
        elseif (nb.eq.4) then
          wna=0.
          wnb=1.
        endif
c
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
            ubflxs_t(i,j)=0.
            ubcors_t(i,j)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
            vbflxs_t(i,j)=0.
            vbcors_t(i,j)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - explicit time integration of barotropic flow (forward-backward
c --- - scheme) in order to combine forward-backward scheme with
c --- - leapfrog treatment of coriolis term, v-eqn must be solved before
c --- - u-eqn every other time step
c
        do lll=lll0,lll0+lstep/2-1
c
          if (mod(lll,2).eq.1) then
c
            wo=woa*lll+wob
            wn=wna*lll+wnb
            wm=1.-wo-wn
c
            call xctilr(pb_t, 1,2, 2,2, halo_ps)
            call xctilr(ubflx_t, 1,2, 2,2, halo_uv)
            call xctilr(vbflx_t, 1,2, 2,3, halo_vv)
c
c --- ----- continuity equation
c
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=-1,jj+2
              do l=1,isp(j)
              do i=max(-1,ifp(j,l)),min(ii+1,ilp(j,l))
                pb_t(i,j,nl)=(1.-wbaro)*pb_t(i,j,ml)+wbaro*pb_t(i,j,nl)
     .           -(1.+wbaro)*dlt*(ubflx_t(i+1,j,ml)-ubflx_t(i,j,ml)
     .                           +vbflx_t(i,j+1,ml)-vbflx_t(i,j,ml))
     .            *scp2i(i,j)
              enddo
              enddo
            enddo
c$OMP END PARALLEL DO
c
c --- ----- u momentum equation
c
            if     (mommth.eq.'enscon') then
c
c --- ------- Sadourny (1975) enstrophy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,utndcy)
              do j=-1,jj+2
                do l=1,isu(j)
                do i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
c
                  ubflxs_t(i,j)=ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl)
     .                         +(1.+wbaro)*ubflx_t(i,j,ml)
c
                  q= (vbflx_t(i  ,j  ,ml)*scvxi(i  ,j  )
     .               +vbflx_t(i  ,j+1,ml)*scvxi(i  ,j+1)
     .               +vbflx_t(i-1,j  ,ml)*scvxi(i-1,j  )
     .               +vbflx_t(i-1,j+1,ml)*scvxi(i-1,j+1))
     .               *(wo*(pvtrop_o(i,j)+pvtrop_o(i,j+1))
     .                +wm*(pvtrop(i,j,m)+pvtrop(i,j+1,m))
     .                +wn*(pvtrop(i,j,n)+pvtrop(i,j+1,n)))*.125
c
                  ubcors_t(i,j)=ubcors_t(i,j)+q
c
                  utndcy=q+
     .             (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl)
     .                               -xixm_o(i,j)*pb_t(i-1,j,nl)))
     .             +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,m)*pb_t(i-1,j,nl)))
     .             +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,n)*pb_t(i-1,j,nl))))
     .             *scuxi(i,j)
c
                  ubflx_t(i,j,nl)=
     .              (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j)
     .                              *min(pb_t(i-1,j,nl),pb_t(i,j,nl))
     .                             -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl)=max(-uminb(i,j),min(umaxb(i,j),
     .                                ubflx_t(i,j,nl)))
c
                enddo
                enddo
              enddo
c$OMP END PARALLEL DO
c
            elseif (mommth.eq.'enecon'.or.mommth.eq.'enedis') then
c
c --- ------- Sadourny (1975) energy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,utndcy)
              do j=-1,jj+2
                do l=1,isu(j)
                do i=max(0,ifu(j,l)),min(ii+1,ilu(j,l))
c
                  ubflxs_t(i,j)=ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl)
     .                         +(1.+wbaro)*ubflx_t(i,j,ml)
c
                  q= .25*( (vbflx_t(i  ,j  ,ml)*scvxi(i  ,j  )
     .                     +vbflx_t(i-1,j  ,ml)*scvxi(i-1,j  ))
     .                    *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m)
     .                     +wn*pvtrop(i,j,n))
     .                    +(vbflx_t(i  ,j+1,ml)*scvxi(i  ,j+1)
     .                     +vbflx_t(i-1,j+1,ml)*scvxi(i-1,j+1))
     .                    *(wo*pvtrop_o(i,j+1)+wm*pvtrop(i,j+1,m)
     .                     +wn*pvtrop(i,j+1,n)))
c
                  ubcors_t(i,j)=ubcors_t(i,j)+q
c
                  utndcy=q+
     .             (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl)
     .                               -xixm_o(i,j)*pb_t(i-1,j,nl)))
     .             +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,m)*pb_t(i-1,j,nl)))
     .             +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,n)*pb_t(i-1,j,nl))))
     .             *scuxi(i,j)
c
                  ubflx_t(i,j,nl)=
     .              (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j)
     .                              *min(pb_t(i-1,j,nl),pb_t(i,j,nl))
     .                             -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl)=max(-uminb(i,j),min(umaxb(i,j),
     .                                ubflx_t(i,j,nl)))
c
                enddo
                enddo
              enddo
c$OMP END PARALLEL DO
c
            else
              if (mnproc.eq.1) then
                write (lp,'(3a)') ' mommth=',trim(mommth),
     .                            ' is unsupported!'
              endif
              call xcstop('(barotp)')
                     stop '(barotp)'
            endif
c
c --- ----- v momentum equation
c
            if     (mommth.eq.'enscon') then
c
c --- ------- Sadourny (1975) enstrophy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,vtndcy)
            do j=0,jj+2
              do l=1,isv(j)
              do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
c
                vbflxs_t(i,j)=vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl)
     .                       +(1.+wbaro)*vbflx_t(i,j,ml)
c
                q=-(ubflx_t(i  ,j  ,nl)*scuyi(i  ,j  )
     .             +ubflx_t(i+1,j  ,nl)*scuyi(i+1,j  )
     .             +ubflx_t(i  ,j-1,nl)*scuyi(i  ,j-1)
     .             +ubflx_t(i+1,j-1,nl)*scuyi(i+1,j-1))
     .            *(wo*(pvtrop_o(i,j)+pvtrop_o(i+1,j))
     .             +wm*(pvtrop(i,j,m)+pvtrop(i+1,j,m))
     .             +wn*(pvtrop(i,j,n)+pvtrop(i+1,j,n)))*.125
c
                vbcors_t(i,j)=vbcors_t(i,j)+q
c
                vtndcy=q+
     .           (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl)
     .                             -xiym_o(i,j)*pb_t(i,j-1,nl)))
     .           +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl)
     .                             -xiym(i,j,m)*pb_t(i,j-1,nl)))
     .           +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl)
     .                             -xiym(i,j,n)*pb_t(i,j-1,nl))))
     .           *scvyi(i,j)
c
                vbflx_t(i,j,nl)=
     .            (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl)
     .           +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j)
     .                            *min(pb_t(i,j-1,nl),pb_t(i,j,nl))
     .                           -vglue(i,j)*vbflx_t(i,j,ml))
                vbflx_t(i,j,nl)=max(-vminb(i,j),min(vmaxb(i,j),
     .                              vbflx_t(i,j,nl)))
c
              enddo 
              enddo 
            enddo 
c$OMP END PARALLEL DO
c
            elseif (mommth.eq.'enecon'.or.mommth.eq.'enedis') then
c
c --- ------- Sadourny (1975) energy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,vtndcy)
              do j=0,jj+2
                do l=1,isv(j)
                do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
c
                  vbflxs_t(i,j)=vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl)
     .                         +(1.+wbaro)*vbflx_t(i,j,ml)
c
                  q=-.25*( (ubflx_t(i  ,j  ,nl)*scuyi(i  ,j  )
     .                     +ubflx_t(i  ,j-1,nl)*scuyi(i  ,j-1))
     .                    *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m)
     .                     +wn*pvtrop(i,j,n))
     .                    +(ubflx_t(i+1,j  ,nl)*scuyi(i+1,j  )
     .                     +ubflx_t(i+1,j-1,nl)*scuyi(i+1,j-1))
     .                    *(wo*pvtrop_o(i+1,j)+wm*pvtrop(i+1,j,m)
     .                     +wn*pvtrop(i+1,j,n)))
c
                  vbcors_t(i,j)=vbcors_t(i,j)+q
c
                  vtndcy=q+
     .             (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl)
     .                               -xiym_o(i,j)*pb_t(i,j-1,nl)))
     .             +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl)
     .                               -xiym(i,j,m)*pb_t(i,j-1,nl)))
     .             +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl)
     .                               -xiym(i,j,n)*pb_t(i,j-1,nl))))
     .             *scvyi(i,j)
c
                  vbflx_t(i,j,nl)=
     .              (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j)
     .                              *min(pb_t(i,j-1,nl),pb_t(i,j,nl))
     .                             -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl)=max(-vminb(i,j),min(vmaxb(i,j),
     .                                vbflx_t(i,j,nl)))
c
                enddo 
                enddo 
              enddo 
c$OMP END PARALLEL DO
c
            else
              if (mnproc.eq.1) then
                write (lp,'(3a)') ' mommth=',trim(mommth),
     .                            ' is unsupported!'
              endif
              call xcstop('(barotp)')
                     stop '(barotp)'
            endif
c
c
            ll=ml
            ml=nl
            nl=ll
c
          else
c
            wo=woa*lll+wob
            wn=wna*lll+wnb
            wm=1.-wo-wn
c
c --- ----- continuity equation
c
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=0,jj+1
              do l=1,isp(j)
              do i=max(0,ifp(j,l)),min(ii,ilp(j,l))
                pb_t(i,j,nl)=(1.-wbaro)*pb_t(i,j,ml)+wbaro*pb_t(i,j,nl)
     .           -(1.+wbaro)*dlt*(ubflx_t(i+1,j,ml)-ubflx_t(i,j,ml)
     .                           +vbflx_t(i,j+1,ml)-vbflx_t(i,j,ml))
     .            *scp2i(i,j)
              enddo
              enddo
            enddo
c$OMP END PARALLEL DO
c
c --- ----- v momentum equation
c
            if     (mommth.eq.'enscon') then
c
c --- ------- Sadourny (1975) enstrophy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,vtndcy)
              do j=1,jj+1
                do l=1,isv(j)
                do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
c
                  vbflxs_t(i,j)=vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl)
     .                         +(1.+wbaro)*vbflx_t(i,j,ml)
c
                  q=-(ubflx_t(i  ,j  ,ml)*scuyi(i  ,j  )
     .               +ubflx_t(i+1,j  ,ml)*scuyi(i+1,j  )
     .               +ubflx_t(i  ,j-1,ml)*scuyi(i  ,j-1)
     .               +ubflx_t(i+1,j-1,ml)*scuyi(i+1,j-1))
     .              *(wo*(pvtrop_o(i,j)+pvtrop_o(i+1,j))
     .               +wm*(pvtrop(i,j,m)+pvtrop(i+1,j,m))
     .               +wn*(pvtrop(i,j,n)+pvtrop(i+1,j,n)))*.125
c
                  vbcors_t(i,j)=vbcors_t(i,j)+q
c
                  vtndcy=q+
     .             (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl)
     .                               -xiym_o(i,j)*pb_t(i,j-1,nl)))
     .             +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl)
     .                               -xiym(i,j,m)*pb_t(i,j-1,nl)))
     .             +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl)
     .                               -xiym(i,j,n)*pb_t(i,j-1,nl))))
     .             *scvyi(i,j)
c
                  vbflx_t(i,j,nl)=
     .              (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j)
     .                              *min(pb_t(i,j-1,nl),pb_t(i,j,nl))
     .                             -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl)=max(-vminb(i,j),min(vmaxb(i,j),
     .                                vbflx_t(i,j,nl)))
c
                enddo 
                enddo 
              enddo 
c$OMP END PARALLEL DO
c
            elseif (mommth.eq.'enecon'.or.mommth.eq.'enedis') then
c
c --- ------- Sadourny (1975) energy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,vtndcy)
              do j=1,jj+1
                do l=1,isv(j)
                do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
c
                  vbflxs_t(i,j)=vbflxs_t(i,j)-wbaro*vbflx_t(i,j,nl)
     .                         +(1.+wbaro)*vbflx_t(i,j,ml)
c
                  q=-.25*( (ubflx_t(i  ,j  ,ml)*scuyi(i  ,j  )
     .                     +ubflx_t(i  ,j-1,ml)*scuyi(i  ,j-1))
     .                    *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m)
     .                     +wn*pvtrop(i,j,n))
     .                    +(ubflx_t(i+1,j  ,ml)*scuyi(i+1,j  )
     .                     +ubflx_t(i+1,j-1,ml)*scuyi(i+1,j-1))
     .                    *(wo*pvtrop_o(i+1,j)+wm*pvtrop(i+1,j,m)
     .                     +wn*pvtrop(i+1,j,n)))
c
                  vbcors_t(i,j)=vbcors_t(i,j)+q
c
                  vtndcy=q+
     .             (wo*(pgfym_o(i,j)-(xiyp_o(i,j)*pb_t(i,j  ,nl)
     .                               -xiym_o(i,j)*pb_t(i,j-1,nl)))
     .             +wm*(pgfym(i,j,m)-(xiyp(i,j,m)*pb_t(i,j  ,nl)
     .                               -xiym(i,j,m)*pb_t(i,j-1,nl)))
     .             +wn*(pgfym(i,j,n)-(xiyp(i,j,n)*pb_t(i,j  ,nl)
     .                               -xiym(i,j,n)*pb_t(i,j-1,nl))))
     .             *scvyi(i,j)
c
                  vbflx_t(i,j,nl)=
     .              (1.-wbaro)*vbflx_t(i,j,ml)+wbaro*vbflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((vtndcy+vtotn(i,j))*scvx(i,j)
     .                              *min(pb_t(i,j-1,nl),pb_t(i,j,nl))
     .                             -vglue(i,j)*vbflx_t(i,j,ml))
                  vbflx_t(i,j,nl)=max(-vminb(i,j),min(vmaxb(i,j),
     .                                vbflx_t(i,j,nl)))
c
                enddo 
                enddo 
              enddo 
c$OMP END PARALLEL DO
c
            else
              if (mnproc.eq.1) then
                write (lp,'(3a)') ' mommth=',trim(mommth),
     .                            ' is unsupported!'
              endif
              call xcstop('(barotp)')
                     stop '(barotp)'
            endif
c
c --- ----- u momentum equation
c
            if     (mommth.eq.'enscon') then
c
c --- ------- Sadourny (1975) enstrophy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,utndcy)
              do j=1,jj
                do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c
                  ubflxs_t(i,j)=ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl)
     .                         +(1.+wbaro)*ubflx_t(i,j,ml)
c
                  q= (vbflx_t(i  ,j  ,nl)*scvxi(i  ,j  )
     .               +vbflx_t(i  ,j+1,nl)*scvxi(i  ,j+1)
     .               +vbflx_t(i-1,j  ,nl)*scvxi(i-1,j  )
     .               +vbflx_t(i-1,j+1,nl)*scvxi(i-1,j+1))
     .               *(wo*(pvtrop_o(i,j)+pvtrop_o(i,j+1))
     .                +wm*(pvtrop(i,j,m)+pvtrop(i,j+1,m))
     .                +wn*(pvtrop(i,j,n)+pvtrop(i,j+1,n)))*.125
c
                  ubcors_t(i,j)=ubcors_t(i,j)+q
c
                  utndcy=q+
     .             (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl)
     .                               -xixm_o(i,j)*pb_t(i-1,j,nl)))
     .             +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,m)*pb_t(i-1,j,nl)))
     .             +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,n)*pb_t(i-1,j,nl))))
     .             *scuxi(i,j)
c
                  ubflx_t(i,j,nl)=
     .              (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j)
     .                              *min(pb_t(i-1,j,nl),pb_t(i,j,nl))
     .                             -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl)=max(-uminb(i,j),min(umaxb(i,j),
     .                                ubflx_t(i,j,nl)))
c
                enddo
                enddo
              enddo
c$OMP END PARALLEL DO
c
            elseif (mommth.eq.'enecon'.or.mommth.eq.'enedis') then
c
c --- ------- Sadourny (1975) energy conserving scheme
c
c$OMP PARALLEL DO PRIVATE(l,i,q,utndcy)
              do j=1,jj
                do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c
                  ubflxs_t(i,j)=ubflxs_t(i,j)-wbaro*ubflx_t(i,j,nl)
     .                         +(1.+wbaro)*ubflx_t(i,j,ml)
c
                  q= .25*( (vbflx_t(i  ,j  ,nl)*scvxi(i  ,j  )
     .                     +vbflx_t(i-1,j  ,nl)*scvxi(i-1,j  ))
     .                    *(wo*pvtrop_o(i,j)+wm*pvtrop(i,j,m)
     .                     +wn*pvtrop(i,j,n))
     .                    +(vbflx_t(i  ,j+1,nl)*scvxi(i  ,j+1)
     .                     +vbflx_t(i-1,j+1,nl)*scvxi(i-1,j+1))
     .                    *(wo*pvtrop_o(i,j+1)+wm*pvtrop(i,j+1,m)
     .                     +wn*pvtrop(i,j+1,n)))
c
                  ubcors_t(i,j)=ubcors_t(i,j)+q
c
                  utndcy=q+
     .             (wo*(pgfxm_o(i,j)-(xixp_o(i,j)*pb_t(i  ,j,nl)
     .                               -xixm_o(i,j)*pb_t(i-1,j,nl)))
     .             +wm*(pgfxm(i,j,m)-(xixp(i,j,m)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,m)*pb_t(i-1,j,nl)))
     .             +wn*(pgfxm(i,j,n)-(xixp(i,j,n)*pb_t(i  ,j,nl)
     .                               -xixm(i,j,n)*pb_t(i-1,j,nl))))
     .             *scuxi(i,j)
c
                  ubflx_t(i,j,nl)=
     .              (1.-wbaro)*ubflx_t(i,j,ml)+wbaro*ubflx_t(i,j,nl)
     .             +(1.+wbaro)*dlt*((utndcy+utotn(i,j))*scuy(i,j)
     .                              *min(pb_t(i-1,j,nl),pb_t(i,j,nl))
     .                             -uglue(i,j)*ubflx_t(i,j,ml))
                  ubflx_t(i,j,nl)=max(-uminb(i,j),min(umaxb(i,j),
     .                                ubflx_t(i,j,nl)))
c
                enddo
                enddo
              enddo
c$OMP END PARALLEL DO
c
            else
              if (mnproc.eq.1) then
                write (lp,'(3a)') ' mommth=',trim(mommth),
     .                            ' is unsupported!'
              endif
              call xcstop('(barotp)')
                     stop '(barotp)'
            endif
c
            ll=ml
            ml=nl
            nl=ll
c
          endif
c
        enddo ! lll
c
        lll0=lll0+lstep/2
c
        if     (nb.eq.1) then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb(i,j,m)=pb_t(i,j,ml)
            enddo
            enddo
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              pbu(i,j,m)=min(pb_t(i,j,ml),pb_t(i-1,j,ml))
              ubflx(i,j,m)=ubflx_t(i,j,ml)
              ub(i,j,m)=ubflx(i,j,m)/(pbu(i,j,m)*scuy(i,j))
              ubflxs(i,j,n)=ubflxs(i,j,n)+ubflxs_t(i,j)
              ubflxs(i,j,m)=ubflxs(i,j,3)+ubflxs_t(i,j)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              pbv(i,j,m)=min(pb_t(i,j,ml),pb_t(i,j-1,ml))
              vbflx(i,j,m)=vbflx_t(i,j,ml)
              vb(i,j,m)=vbflx(i,j,m)/(pbv(i,j,m)*scvx(i,j))
              vbflxs(i,j,n)=vbflxs(i,j,n)+vbflxs_t(i,j)
              vbflxs(i,j,m)=vbflxs(i,j,3)+vbflxs_t(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (nb.eq.2) then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb_mn(i,j,ml)=pb_t(i,j,ml)
              pb_mn(i,j,nl)=pb_t(i,j,nl)
            enddo
            enddo
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              ubflx_mn(i,j,ml)=ubflx_t(i,j,ml)
              ubflx_mn(i,j,nl)=ubflx_t(i,j,nl)
              ubflxs(i,j,m)=ubflxs(i,j,m)+ubflxs_t(i,j)
              ubflxs(i,j,3)=ubflxs_t(i,j)
              ubflxs_p(i,j,n)=ubflxs_t(i,j)
              ubcors_p(i,j)=ubcors_t(i,j)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vbflx_mn(i,j,ml)=vbflx_t(i,j,ml)
              vbflx_mn(i,j,nl)=vbflx_t(i,j,nl)
              vbflxs(i,j,m)=vbflxs(i,j,m)+vbflxs_t(i,j)
              vbflxs(i,j,3)=vbflxs_t(i,j)
              vbflxs_p(i,j,n)=vbflxs_t(i,j)
              vbcors_p(i,j)=vbcors_t(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (nb.eq.3) then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb(i,j,n)=pb_t(i,j,ml)
            enddo
            enddo
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              pbu(i,j,n)=min(pb_t(i,j,ml),pb_t(i-1,j,ml))
              ubflx(i,j,n)=ubflx_t(i,j,ml)
              ub(i,j,n)=ubflx(i,j,n)/(pbu(i,j,n)*scuy(i,j))
              ubflxs_p(i,j,m)=ubflxs(i,j,m)+ubflxs_t(i,j)
              ubflxs_p(i,j,n)=ubflxs_p(i,j,n)+ubflxs_t(i,j)
              ubcors_p(i,j)=ubcors_p(i,j)+ubcors_t(i,j)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              pbv(i,j,n)=min(pb_t(i,j,ml),pb_t(i,j-1,ml))
              vbflx(i,j,n)=vbflx_t(i,j,ml)
              vb(i,j,n)=vbflx(i,j,n)/(pbv(i,j,n)*scvx(i,j))
              vbflxs_p(i,j,m)=vbflxs(i,j,m)+vbflxs_t(i,j)
              vbflxs_p(i,j,n)=vbflxs_p(i,j,n)+vbflxs_t(i,j)
              vbcors_p(i,j)=vbcors_p(i,j)+vbcors_t(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (nb.eq.4) then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              ubflxs_p(i,j,n)=ubflxs_p(i,j,n)+ubflxs_t(i,j)
              ubcors_p(i,j)=ubcors_p(i,j)+ubcors_t(i,j)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vbflxs_p(i,j,n)=vbflxs_p(i,j,n)+vbflxs_t(i,j)
              vbcors_p(i,j)=vbcors_p(i,j)+vbcors_t(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (nb.eq.5) then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              pb_p(i,j)=pb_t(i,j,ml)
            enddo
            enddo
            do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              pbu_p(i,j)=min(pb_t(i,j,ml),pb_t(i-1,j,ml))
              ubflxs_p(i,j,n)=ubflxs_p(i,j,n)+ubflxs_t(i,j)
              ubcors_p(i,j)=ubcors_p(i,j)+ubcors_t(i,j)
            enddo
            enddo
            do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              pbv_p(i,j)=min(pb_t(i,j,ml),pb_t(i,j-1,ml))
              vbflxs_p(i,j,n)=vbflxs_p(i,j,n)+vbflxs_t(i,j)
              vbcors_p(i,j)=vbcors_p(i,j)+vbcors_t(i,j)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
c
      enddo ! nb
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'barotp:'
        endif
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
      endif
c
      return
      end
