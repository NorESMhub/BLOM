c routine tested with TRC off, and no copyout from device to host, just interested in the timing of loop construct.
#undef TRC
      
!------------------------------------------------------------------------------
! Copyright (C) 2006-2018 Mats Bentsen
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


      module mod_diffus
      use mod_xc

      real delt1
      logical csdiag,cnsvdi,edsprs
      
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .     u,v,                 ! velocity components
     .     dp,                  ! layer thickness
     .     temp,                ! temperature
     .     saln,                ! salinity
     .     sigma                !,         ! potential density

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .     difiso               !,        ! isopycnal diffusivity

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .     utflld,vtflld,       ! horizontal heat fluxes due to lateral diffusion
     .     usfltd,vsfltd,       ! horizontal salt fluxes due to thickness diffusion
     .     usflld,vsflld        ! horizontal salt fluxes due to lateral diffusion

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .     utflx,vtflx,         ! horizontal heat fluxes
     .     usflx,vsflx          ! horizontal salt fluxes

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .     scux,scuy,           ! mesh size at u-points in x,y direction
     .     scvx,scvy,           ! mesh size at v-points in x,y direction
     .     scq2,scp2,           ! grid box size at q- and p-points
     .     scuxi,scvyi          !,   ! inverses of scux,scvy

c openacc2.6:      
c If the declare directive appears in a global context, then the data in var-list is statically allocated
c in device memory and the structured reference counter is set to one.
c (is a module considered global context?)
      
cacc declare device_resident(      
c$acc declare create(
c$acc.     delt1,csdiag,cnsvdi,edsprs,
c$acc.     u,v,                 ! velocity components
c$acc.     dp,                  ! layer thickness
c$acc.     temp,                ! temperature
c$acc.     saln,                ! salinity
c$acc.     sigma,               !,         ! potential density
c$acc.     difiso,              !,        ! isopycnal diffusivity
c$acc.     utflld,vtflld,       ! horizontal heat fluxes due to lateral diffusion
c$acc.     usfltd,vsfltd,       ! horizontal salt fluxes due to thickness diffusion
c$acc.     usflld,vsflld,       ! horizontal salt fluxes due to lateral diffusion
c$acc.     utflx,vtflx,         ! horizontal heat fluxes
c$acc.     usflx,vsflx,          ! horizontal salt fluxes
c$acc.     scux,scuy,           ! mesh size at u-points in x,y direction
c$acc.     scvx,scvy,           ! mesh size at v-points in x,y direction
c$acc.     scq2,scp2,           ! grid box size at q- and p-points
c$acc.     scuxi,scvyi)         !,   ! inverses of scux,scvy
      
      end module

      
      subroutine update_mod_diffus
c      
c copy common block variables needed in diffus.F into module mod_diffus,
c and initiate an update of the variables also on the gpu(device). when
c we have everything in modules we only need the "update host/device" part.
c

C  USE <MODULE NAME>, <LOCAL ALIAS> => <VARIABLE NAME IN MODULE>
      
      use mod_diffus,
     . md_delt1 => delt1, md_csdiag=>csdiag,md_cnsvdi=>cnsvdi,md_edsprs=>edsprs,     
     . md_u       => u      ,
     . md_v       => v      ,
     . md_dp      => dp     ,
     . md_temp    => temp   ,
     . md_saln    => saln   ,
     . md_sigma   => sigma  ,
     . md_difiso  => difiso ,
     . md_utflld  => utflld ,
     . md_vtflld  => vtflld ,
     . md_usfltd  => usfltd ,
     . md_vsfltd  => vsfltd ,
     . md_usflld  => usflld ,
     . md_vsflld  => vsflld ,
     . md_utflx   => utflx  ,
     . md_vtflx   => vtflx  ,
     . md_usflx   => usflx  ,
     . md_vsflx   => vsflx  ,
     . md_scux    => scux   ,
     . md_scuy    => scuy   ,
     . md_scvx    => scvx   ,
     . md_scvy    => scvy   ,
     . md_scq2    => scq2   ,
     . md_scp2    => scp2   ,
     . md_scuxi   => scuxi  ,
     . md_scvyi   => scvyi  
      
      implicit none
      
#include "common_blocks.h"

c update host module variables
      
      md_u       = u      
      md_v       = v      
      md_dp      = dp     
      md_temp    = temp   
      md_saln    = saln   
      md_sigma   = sigma  
      md_difiso  = difiso 
      md_utflld  = utflld 
      md_vtflld  = vtflld 
      md_usfltd  = usfltd 
      md_vsfltd  = vsfltd 
      md_usflld  = usflld 
      md_vsflld  = vsflld 
      md_utflx   = utflx  
      md_vtflx   = vtflx  
      md_usflx   = usflx  
      md_vsflx   = vsflx  
      md_scux    = scux   
      md_scuy    = scuy   
      md_scvx    = scvx   
      md_scvy    = scvy   
      md_scq2    = scq2   
      md_scp2    = scp2   
      md_scuxi   = scuxi  
      md_scvyi   = scvyi  

c then the device module variables
      
c$acc update device(md_u,
c$acc. md_v,     
c$acc. md_dp,    
c$acc. md_temp,  
c$acc. md_saln,  
c$acc. md_sigma, 
c$acc. md_difiso,
c$acc. md_utflld,
c$acc. md_vtflld,
c$acc. md_usfltd,
c$acc. md_vsfltd,
c$acc. md_usflld,
c$acc. md_vsflld,
c$acc. md_utflx, 
c$acc. md_vtflx, 
c$acc. md_usflx, 
c$acc. md_vsflx, 
c$acc. md_scux , 
c$acc. md_scuy,  
c$acc. md_scvx,  
c$acc. md_scvy,  
c$acc. md_scq2,  
c$acc. md_scp2,  
c$acc. md_scuxi, 
c$acc. md_scvyi)
      
      end subroutine
    
  
      
      
      subroutine diffus(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- diffusion of tracers
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_diffus
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
      integer i,j,k,l,kn,km
      real q
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real sig
      external sig
!$acc routine (sig)
      
c
      real dpeps
      parameter (dpeps=1.e-4)
c
      call xctilr(dp  (1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
#ifdef TRC
      do nt=1,ntr
#  if defined(TKE) && !defined(TKEIDF)
        if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
        call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1,kk, 2,2, halo_ps)
      enddo
#endif
      call xctilr(difiso, 1,kk, 2,2, halo_ps)
c
      do k=1,kk
        kn=k+nn
        km=k+mm
c
c$OMP PARALLEL DO PRIVATE(l,i,q
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )       
c        do j=0,jj+1
c          do l=1,isu(j)
c          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
        
c tiles gave ca the same perf as collapse.        
cacc parallel loop private(i,j,q) tile(4,4)
c$acc parallel loop private(i,j,q) collapse(2)
        do j=0,jj+1
           do i=0,ii+2
              if(iu(i,j) > 0)then
             
                 q=delt1*.5*(difiso(i-1,j,k)+difiso(i,j,k))
     .                *scuy(i,j)*scuxi(i,j)
     .                *max(min(dp(i-1,j,kn),dp(i,j,kn)),dpeps)
                 usflld(i,j,km)=q*(saln(i-1,j,kn)-saln(i,j,kn))
                 utflld(i,j,km)=q*(temp(i-1,j,kn)-temp(i,j,kn))
#ifdef TRC
                 do nt=1,ntr
#  if defined(TKE) && !defined(TKEIDF)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                    uflxtr(nt,i,j)=q*(trc(i-1,j,kn,nt)-trc(i,j,kn,nt))
                 enddo
#endif
                 usflx(i,j,km)=usflx(i,j,km)+usflld(i,j,km)
                 utflx(i,j,km)=utflx(i,j,km)+utflld(i,j,km)

              end if            ! iu
            
           enddo
        enddo
c$acc end parallel loop        
c     enddo
c$OMP END PARALLEL DO




        
c$OMP PARALLEL DO PRIVATE(l,i,q
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
c        do j=0,jj+2
c          do l=1,isv(j)
c          do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
cacc parallel loop private(i,j,q) tile(4,4)
c$acc parallel loop private(i,j,q) collapse(2)
        do j=0,jj+2
           do i=0,ii+1
              if(iv(i,j) > 0)then
                 q=delt1*.5*(difiso(i,j-1,k)+difiso(i,j,k))
     .                *scvx(i,j)*scvyi(i,j)
     .                *max(min(dp(i,j-1,kn),dp(i,j,kn)),dpeps)
                 vsflld(i,j,km)=q*(saln(i,j-1,kn)-saln(i,j,kn))
                 vtflld(i,j,km)=q*(temp(i,j-1,kn)-temp(i,j,kn))
#ifdef TRC
                 do nt=1,ntr
#  if defined(TKE) && !defined(TKEIDF)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                    vflxtr(nt,i,j)=q*(trc(i,j-1,kn,nt)-trc(i,j,kn,nt))
                 enddo
#endif
                 vsflx(i,j,km)=vsflx(i,j,km)+vsflld(i,j,km)
                 vtflx(i,j,km)=vtflx(i,j,km)+vtflld(i,j,km)
              end if
           enddo
        enddo
c$acc end parallel loop        
        
c     enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
c        do j=0,jj+1
c          do l=1,isp(j)
c          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
cacc parallel loop private(i,j,q) tile(4,4)
c$acc parallel loop private(i,j,q) collapse(2)
        do j=0,jj+2
           do i=0,ii+1
              if(ip(i,j) > 0)then
                 
                 q=1./(scp2(i,j)*max(dp(i,j,kn),dpeps))
                 saln(i,j,kn)=saln(i,j,kn)
     .                -q*(usflld(i+1,j,km)-usflld(i,j,km)
     .                +vsflld(i,j+1,km)-vsflld(i,j,km))
                 temp(i,j,kn)=temp(i,j,kn)
     .                -q*(utflld(i+1,j,km)-utflld(i,j,km)
     .                +vtflld(i,j+1,km)-vtflld(i,j,km))
#ifdef TRC
                 do nt=1,ntr
#  if defined(TKE) && !defined(TKEIDF)
                    if (nt.eq.itrtke.or.nt.eq.itrgls) cycle
#  endif
                    trc(i,j,kn,nt)=trc(i,j,kn,nt)
     .                   -q*(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .                   +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))
                 enddo
#endif
                 sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
              end if
           enddo
        enddo
c$acc end parallel loop        
c     enddo
c$OMP END PARALLEL DO
c
      enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'diffus:'
        endif
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      return
      end
