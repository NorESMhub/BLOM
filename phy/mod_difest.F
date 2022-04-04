! ------------------------------------------------------------------------------
! Copyright (C) 2009-2022 Mats Bentsen, Mehmet Ilicak
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
! ----------------------------------------------------------------------------

      module mod_difest
c
      use mod_types, only: r8
      use mod_constants, only: g, alpha0, pi, epsil, spval, onem, onecm
      use mod_time, only: delt1
      use mod_xc
      use mod_vcoord, only: vcoord_type_tag, isopyc_bulkml,
     .                      cntiso_hybrid, sigmar
      use mod_grid, only: scpx, scpy, scp2,
     .                    plat, coriop, betafp, cosang, sinang
      use mod_eos, only: rho
      use mod_state, only:  u, v, dp, dpu, dpv, temp, saln, sigma, p,
     .                      pbu, pbv, ubflxs_p, vbflxs_p, kfpla
      use mod_diffusion, only: egc, eggam, eglsmn, egmndf, egmxdf,
     .                         egidfq, ri0, bdmc1, bdmc2, tkepf, bdmtyp,
     .                         edsprs, edritp_opt, edritp_shear,
     .                         edritp_large_scale, edwmth_opt,
     .                         edwmth_smooth, edwmth_step,
     .                         difint, difiso, difdia, difmxp, difwgt,
     .                         Kvisc_m, Kdiff_t, Kdiff_s, 
     .                         t_ns_nonloc, s_nonloc
      use mod_cmnfld, only: bfsqi, nnslpx, nnslpy, mlts
      use mod_forcing, only: ustar, ustarb, ustar3, buoyfl, t_sw_nonloc
     .                       , surflx, sswflx, salflx
      use mod_tidaldissip, only: twedon
      use mod_niw, only: niwgf, niwbf, niwlf, idkedt, niw_ke_tendency
      use mod_seaice, only: ficem
      use mod_utility, only: util1
      use mod_checksum, only: csdiag, chksummsk
      use CVMix_kpp, only : CVMix_coeffs_kpp
      use CVMix_kpp, only : CVMix_kpp_compute_turbulent_scales
      use CVMix_kpp, only : CVMix_kpp_compute_bulk_Richardson
      use CVMix_kpp, only : CVMix_kpp_compute_OBL_depth
      use CVMix_kpp, only : CVmix_kpp_compute_unresolved_shear
      use CVMix_kpp, only : CVMix_kpp_compute_kOBL_depth
      use CVMix_shear, only : CVMix_init_shear, CVMix_coeffs_shear
      use CVMix_background, only : CVMix_init_bkgnd, CVMix_coeffs_bkgnd
      use CVMix_convection, only : CVMix_init_conv, CVMix_coeffs_conv
      use CVMix_tidal, only : CVMix_init_tidal, CVMix_coeffs_tidal
      use CVMix_tidal, only : CVMix_compute_Simmons_invariant,
     .                        CVMix_tidal_params_type
      use CVMix_kinds_and_types,  only : CVMix_global_params_type
      use CVMix_kpp,  only : CVMix_kpp_params_type
      use CVMix_kpp,  only : CVMix_put_kpp
      use CVMix_kpp,  only : CVMix_init_kpp
      use CVMix_put_get,  only : CVMix_put
#if defined(TRC) && defined(TKE)
      use mod_tracers, only: itrtke, itrgls, trc
      use mod_tke, only: gls_cmu0, Pr_t, tke_min, gls_psi_min, gls_p,
     .                   gls_m, gls_n, gls_c1, gls_c2, gls_c3plus,
     .                   gls_c3minus, gls_Gh0, gls_Ghmin, gls_Ghcri,
     .                   Ls_unlmt_min, Prod, Buoy, Shear2, L_scale,
     .                   gls_s0, gls_s1, gls_s2, gls_s4, gls_s5, gls_s6,
     .                   gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5,
     .                   sqrt2, cmu_fac1, cmu_fac2, cmu_fac3, tke_exp1,
     .                   gls_exp1, gls_fac6
#endif
c
      implicit none
c
      private
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     .  rig
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  du2l,drhol,up,vp
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  OBLdepth
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  mskv,msku
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  kmax,kfil,kOBL
c
      type(CVMix_tidal_params_type)   :: CVMix_tidal_params
      type(CVMix_global_params_type)  :: CVMix_glb_params   !<CVMix-specific for Prandtl number only
      type(CVMix_kpp_params_type)     :: KPP_params
c      type(CVMix_kpp_params_type), pointer :: CVmix_kpp_params_in
c      type(CVMix_kpp_params_type) ::  CVmix_kpp_params_in
c
c --- parameters:
c ---   iidtyp - type of interface and isopycnal diffusivities. If
c ---            iidtyp=1 the diffusivities are diffusive velocities
c ---            multiplied by the local horizontal grid scale, if
c ---            iidtyp=2 the diffusivities are parameterized according
c ---            to Eden and Greatbatch (2008).
c ---   bdmldp - If bdmldp=1, make the background mixing latitude
c ---            dependent according to Gregg et al. (2003).
c ---   tdmflg - If tdmflg=1, apply tidally driven diapycnal mixing.
c ---   iwdflg - If iwdflg=1, reduce background diapycnal diffusivity
c ---            due to internal wave damping under sea-ice.
c ---   dpbmin - smallest layer thickness allowed in evaluating
c ---            local gradient richardson number [g/cm/s**2].
c ---   drhomn - minimum density difference in evaluations the
c ---            Brunt-Vaisala frequency and the local gradient
c ---            Richardson number [g/cm*3].
c ---   thkdff - diffusive velocity for thickness diffusion [cm/s].
c ---   temdff - diffusive velocity for tracer isopycnal diffusion
c ---            [cm/s].
c ---   nu0    - diapycnal diffusivity when range of isopycnic physical
c ---            layers is restricted [cm**2/s].
c ---   nus0   - maximum shear driven diapycnal diffusivity
c ---            [cm**2/s].
c ---   nug0   - maximum gravity current diapycnal diffusivity
c ---            [cm**2/s].
c ---   drho0  - critical local interface density difference [g/cm**3]
c ---   nuls0  - maximum diapycnal diffusivity applied when local
c ---            stability is weak [cm**2/s].
c ---   iwdfac - internal wave dissipation factor under sea ice [].
c ---   dmxeff - diapycnal mixing efficiency [].
c ---   tdmq   - tidal dissipation efficiency [].
c ---   tdmls0 - tidal driven mixing length scale below critical
c ---            latitude [g/cm/s**2].
c ---   tdmls1 - tidal driven mixing length scale above critical
c ---            latitude [g/cm/s**2].
c ---   tdclat - critical latitude for tide M2 propagation [].
c ---   tddlat - latitudinal transition zone for different tidal driven
c ---            mixing length scales near the critical latitude.
c ---   tkepls - length scale of surface TKE penetration beneath the
c ---            mixed layer [g/cm/s**2]
c ---   niwls  - near-inertial waves driven mixing length scale
c ---            beneath the mixed layer [g/cm/s**2].
c ---   cori30 - coriolis parameter at 30N [1/s].
c ---   bvf0   - reference stratification in the parameterization of
c ---            latitude dependent background diapycnal mixing [1/s].
c ---   nubmin - minimum background diapycnal diffusivity [cm**2/s].
c ---   dpgc   - thickness of region near the bottom where the maximum
c ---            diffusivity is increased due to gravity current mixing
c ---            processes [g/cm/s**2].
c ---   dpgrav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate upper ocean Eady growth rate
c ---            [g/cm/s**2].
c ---   dpdiav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate lateral diffusivities in the
c ---            non-isopycnic layers [g/cm/s**2].
c ---   dpddav - thickness of region below the non-isopycnic surface
c ---            layers used to estimate diapycnal diffusivities in the
c ---            non-isopycnic layers [g/cm/s**2].
c ---   dpnbav - thickness of region near the bottom used to estimate
c ---            bottom Brunt-Vaisala frequency [g/cm/s**2].
      integer iidtyp,bdmldp,tdmflg,iwdflg
      real dptmin,dpbmin,drhomn,thkdff,temdff,nu0,nus0,nug0,drho0,nuls0,
     .     iwdfac,dmxeff,tdmq,tdmls0,tdmls1,tdclat,tddlat,tkepls,niwls,
     .     cori30,bvf0,nubmin,dpgc,dpgrav,dpdiav,dpddav,dpnbav,ustmin,
     .     kappa,bfeps,sleps,zetas,as,cs,minOBLdepth
      parameter (iidtyp=2,bdmldp=1,tdmflg=1,iwdflg=1,dptmin=98060.,
     .           dpbmin=980.6,drhomn=6.e-6,thkdff=.5,temdff=.35,nu0=.1,
     .           nus0=50.,nug0=2500.,drho0=6.e-6,nuls0=500.,iwdfac=.06,
     .           dmxeff=.2,tdmq=1./3.,tdmls0=500.*98060.,
     .           tdmls1=100.*98060.,tdclat=74.5,tddlat=3.,
     .           tkepls=20.*98060.,niwls=300.*98060.,cori30=7.2722e-5,
     .           bvf0=5.24e-3,nubmin=.01,dpgc=300.*98060.,
     .           dpgrav=100.*98060.,dpdiav=100.*98060.,
     .           dpddav=10.*98060.,dpnbav=250.*98060.,ustmin=.1,
     .           kappa=.4,bfeps=1.e-12,sleps=.1,zetas=-1.,as=-28.86,
     .           cs=98.96,minOBLdepth=1.0)
c
      public :: OBLdepth, inivar_difest, init_difest, difest_isobml,
     .          difest_lateral_hybrid, difest_vertical_hybrid, kOBL
c
      contains
c
      subroutine inivar_difest
c
c --- ------------------------------------------------------------------
c --- Initialize arrays.
c --- ------------------------------------------------------------------
c
      integer :: i,j,k
c
c$OMP PARALLEL DO PRIVATE(i,k)
      do j=1-nbdy,jj+nbdy
        do k=1,kk+1
          do i=1-nbdy,ii+nbdy
            rig(i,j,k)=spval
          enddo
        enddo
        do k=1,kk
          do i=1-nbdy,ii+nbdy
            du2l(i,j,k)=spval
            drhol(i,j,k)=spval
            up(i,j,k)=spval
            vp(i,j,k)=spval
          enddo
        enddo
        do i=1-nbdy,ii+nbdy
          OBLdepth(i,j)=spval
          kOBL(i,j)=2
       enddo
      enddo
c$OMP END PARALLEL DO
c
      end subroutine inivar_difest
c
      subroutine init_difest
c
c --- ------------------------------------------------------------------
c --- Initialize CVmix variables.
c --- ------------------------------------------------------------------
c
c -- ------- Background diapycnal mixing.
c         The Bryan-Lewis parameterization is based on the following:
c         \begin{eqnarray*}
c         \kappa_{BL} &=& \textrm{bl1} +
c                   \frac{\textrm{bl2}}{\pi}\tan^{-1}\bigg(
c                        \textrm{bl3}(|z|-\textrm{bl4})\bigg)\\
c         \nu_{BL} &=& \textrm{Pr}\cdot\kappa_{BL}
c         \end{eqnarray*}
c --- ------ Diapycnal mixing when local stability is weak
c --- ------ convection routine based on N2 not rho
c --- ------ if lBruntVaisala is TRUE, otherwise based on rho
c --- ------ convert nuls0 to m2/s
          call CVMix_init_conv(convect_diff=20.0*nuls0*1e-4,
     .                         convect_visc=20.0*nuls0*1e-4,
     .                         lBruntVaisala=.true.,
     .                         BVsqr_convect=0.0)
          call CVMix_put(CVMix_glb_params,'max_nlev',kk)
          call CVMix_put(CVMix_glb_params,'Prandtl',1.0)
          call CVMix_put(CVMix_glb_params,'FreshWaterDensity',1000.0)
          call CVMix_put(CVMix_glb_params,'SaltWaterDensity',1025.0)
          call cvmix_init_shear(mix_scheme='KPP',
     .                           KPP_nu_zero=nus0*1e-4,
     .                           KPP_Ri_zero=ri0,
     .                           KPP_exp=3.0)
         !  CVmix_kpp_params_in => CVmix_kpp_params_user
          call CVMix_init_kpp(Ri_crit=0.3,
     .          minOBLdepth=minOBLdepth,
     .          minVtsqr=1e-10,
     .          vonKarman=0.4,
     .          surf_layer_ext=0.1,
     .          interp_type='quadratic',
     .          interp_type2='LMD94',
     .          lEkman=.false.,
     .          lMonOb=.false.,
     .          MatchTechnique='SimpleShapes',
     .          lenhanced_diff=.true.,
     .          lnonzero_surf_nonlocal=.false.  ,
     .          lnoDGat1=.true.                 ,
     .          CVMix_kpp_params_user=KPP_params )
c         call CVMix_init_kpp(Ri_crit=0.3,
c    .          minOBLdepth=minOBLdepth,
c    .          minVtsqr=1e-10,
c    .          vonKarman=0.4,
c    .          surf_layer_ext=0.1,
c    .          interp_type='quadratic',
c    .          interp_type2='LMD94',
c    .          lEkman=.false.,
c    .          lMonOb=.false.,
c    .          MatchTechnique='MatchGradient',
c    .          lenhanced_diff=.true.,
c    .          lnonzero_surf_nonlocal=.false.  ,
c    .          lnoDGat1=.false.                 ,
c    .          CVMix_kpp_params_user=KPP_params )
c         call CVMix_init_kpp(Ri_crit=0.3,
c    .          minOBLdepth=minOBLdepth,
c    .          minVtsqr=1e-10,
c    .          vonKarman=0.4,
c    .          surf_layer_ext=0.1,
c    .          interp_type='quadratic',
c    .          interp_type2='LMD94',
c    .          lEkman=.false.,
c    .          lMonOb=.false.,
c    .          MatchTechnique='ParabolicNonLocal',
c    .          lenhanced_diff=.true.,
c    .          lnonzero_surf_nonlocal=.true.  ,
c    .          lnoDGat1=.true.                 ,
c    .          CVMix_kpp_params_user=KPP_params )
c
c
      end subroutine init_difest
c
      subroutine difest_common_iso(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Obtain common fields for the estimation of lateral and vertical
c --- diffusivities diapycnal diffusivities when vcoord_type_tag ==
c --- isopyc_bulkml.
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dv2
      real, dimension(1-nbdy:idm+nbdy,kdm) :: du2
      real, dimension(1-nbdy:idm+nbdy) :: tup
      integer, dimension(1-nbdy:idm+nbdy) :: kfpl,klpl
      integer i,j,k,l,kn
      real q
c
c --- Locate the range of layers to be considered in the computation of
c --- diffusivities.
      do j=0,jj+1
        do i=0,ii+1
          kmax(i,j)=0
        enddo
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmax(i,j)=1
          do k=3,kk
            kn=k+nn
            if (dp(i,j,kn).gt.dpbmin) kmax(i,j)=k
          enddo
          if (kfpla(i,j,n).ge.kmax(i,j)) then
            kfil(i,j)=kfpla(i,j,n)+1
          else
            if (sigma(i,j,kfpla(i,j,n)+nn).lt.
     .          .5*(sigmar(i,j,kfpla(i,j,n)  )
     .             +sigmar(i,j,kfpla(i,j,n)+1))) then
              kfil(i,j)=kfpla(i,j,n)+1
            else
              kfil(i,j)=kfpla(i,j,n)+2
            endif
          endif
        enddo
        enddo
      enddo
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kfil(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kfil(i,j)=nint(util1(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute squared vertical velocity difference of v-component
c$OMP PARALLEL DO PRIVATE(l,i,kfpl,klpl,k,kn,q,tup)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          kfpl(i)=kk+1
          klpl(i)=1
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (dpv(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,4,-1
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k.ge.max(kfil(i,j-1),kfil(i,j)).and.
     .          dpv(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
        do k=1,kk
          kn=k+nn
          do i=1,ii
            dv2(i,j,k)=0.
            mskv(i,j,k)=0
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k.ge.kfpl(i).and.k.le.klpl(i).and.
     .          klpl(i)-kfpl(i).ge.1) then
              if     (k.eq.kfpl(i)) then
                q=v(i,j,kn+1)-v(i,j,kn)
                q=q*q
                dv2(i,j,k)=q
                tup(i)=q
              elseif (k.lt.klpl(i)) then
                q=v(i,j,kn+1)-v(i,j,kn)
                q=q*q
                dv2(i,j,k)=.5*(tup(i)+q)
                tup(i)=q
              else
                dv2(i,j,k)=tup(i)
              endif
              mskv(i,j,k)=1
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,kfpl,klpl,k,kn,du2,q,tup)
      do j=1,jj
c
c ----- Compute squared vertical velocity difference of u-component
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          kfpl(i)=kk+1
          klpl(i)=1
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (dpu(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=kk,4,-1
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k.ge.min(kfil(i-1,j),kfil(i,j)).and.
     .          dpu(i,j,kn).gt.dptmin) kfpl(i)=k
          enddo
          enddo
        enddo
        do k=1,kk
          kn=k+nn
          do i=1,ii+1
            du2(i,k)=0.
            msku(i,j,k)=0
          enddo
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k.ge.kfpl(i).and.k.le.klpl(i).and.
     .          klpl(i)-kfpl(i).ge.1) then
              if     (k.eq.kfpl(i)) then
                q=u(i,j,kn+1)-u(i,j,kn)
                q=q*q
                du2(i,k)=q
                tup(i)=q
              elseif (k.lt.klpl(i)) then
                q=u(i,j,kn+1)-u(i,j,kn)
                q=q*q
                du2(i,k)=.5*(tup(i)+q)
                tup(i)=q
              else
                du2(i,k)=tup(i)
              endif
              msku(i,j,k)=1
            endif
          enddo
          enddo
        enddo
c
c --- - Centered at layers, compute vertical in-situ density difference,
c --- - squared vertical velocity difference and local gradient
c --- - Richardson number.
        do k=4,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .          kmax(i,j)-kfil(i,j).ge.1) then
c
c --- ------- Vertical in-situ density difference.
              if     (k.eq.kfil(i,j)) then
                q=max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1))
     .                  -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drhol(i,j,k)=q
                tup(i)=q
              elseif (k.lt.kmax(i,j)) then
                q=max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1))
     .                  -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drhol(i,j,k)=2.*tup(i)*q/max(1.e-14,tup(i)+q)
                tup(i)=q
              else
                drhol(i,j,k)=tup(i)
              endif
c
c --- ------- Vertical squared velocity difference.
              du2l(i,j,k)=(msku(i  ,j,k)*du2(i  ,k)
     .                    +msku(i+1,j,k)*du2(i+1,k))
     .                    /max(1,msku(i,j,k)+msku(i+1,j,k))
     .                   +(mskv(i,j  ,k)*dv2(i,j  ,k)
     .                    +mskv(i,j+1,k)*dv2(i,j+1,k))
     .                    /max(1,mskv(i,j,k)+mskv(i,j+1,k))
c
c --- ------- Local gradient Richardson number.
              rig(i,j,k)=alpha0*alpha0
     .                   *max(drhomn,drhol(i,j,k))*dp(i,j,kn)
     .                   /max(1.e-9,du2l(i,j,k))
c
            endif
          enddo
          enddo
        enddo
c
      enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_common_iso:'
        endif
        call chksummsk(drhol,ip,kk,'drhol')
        call chksummsk(du2l,ip,kk,'du2l')
        call chksummsk(rig,ip,kk,'rig')
      endif
c
      end subroutine difest_common_iso
c
      subroutine difest_common_hyb(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Obtain common fields for the estimation of lateral and vertical
c --- diffusivities diapycnal diffusivities when vcoord_type_tag ==
c --- isopyc_bulkml.
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dv2
      real, dimension(1-nbdy:idm+nbdy,kdm) :: du2
      integer, dimension(1-nbdy:idm+nbdy) :: klpl
      integer i,j,k,l,kn
      real q,dz
c
c --- Compute squared vertical velocity difference of v-component
c$OMP PARALLEL DO PRIVATE(l,i,klpl,k,kn,q)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          klpl(i)=1
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (dpv(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=2,kk
          do i=1,ii
            dv2(i,j,k)=0.
            mskv(i,j,k)=0
          enddo
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (klpl(i).ge.2.and.k.le.klpl(i)) then
              q=v(i,j,kn)-v(i,j,kn-1)
              dv2(i,j,k)=q*q
              mskv(i,j,k)=1
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,klpl,k,kn,du2,q,dz)
      do j=1,jj
c
c ----- Compute squared vertical velocity difference of u-component
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          klpl(i)=1
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (dpu(i,j,kn).gt.dpbmin) klpl(i)=k
          enddo
          enddo
        enddo
        do k=2,kk
          do i=1,ii+1
            du2(i,k)=0.
            msku(i,j,k)=0
          enddo
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (klpl(i).ge.2.and.k.le.klpl(i)) then
              q=u(i,j,kn)-u(i,j,kn-1)
              du2(i,k)=q*q
              msku(i,j,k)=1
            endif
          enddo
          enddo
        enddo
c
c --- - Compute local gradient Richardson number at interfaces.
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          rig(i,j,1)=0.
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (msku(i,j,k)+msku(i+1,j,k)
     .         +mskv(i,j,k)+mskv(i,j+1,k).gt.0) then
              q=(msku(i,j,k)*du2(i,k)  +msku(i+1,j,k)*du2(i+1,k))
     .          /max(1,msku(i,j,k)+msku(i+1,j,k))
     .         +(mskv(i,j,k)*dv2(i,j,k)+mskv(i,j+1,k)*dv2(i,j+1,k))
     .          /max(1,mskv(i,j,k)+mskv(i,j+1,k))
              dz=.5*(dp(i,j,kn-1)+dp(i,j,kn))*alpha0/g
              rig(i,j,k)=max(0.,bfsqi(i,j,k)*dz*dz)/max(1.e-9,q)
            else
              rig(i,j,k)=rig(i,j,k-1)
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          rig(i,j,1)=rig(i,j,2)
          rig(i,j,kk+1)=rig(i,j,kk)
        enddo
        enddo
c
c --- - Compute velocity components at p-points.
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            up(i,j,k)=(msku(i,j,k)*u(i,j,kn)+msku(i+1,j,k)*u(i+1,j,kn))
     .                /max(1,msku(i,j,k)+msku(i+1,j,k))
            vp(i,j,k)=(mskv(i,j,k)*v(i,j,kn)+mskv(i,j+1,k)*v(i,j+1,kn))
     .                /max(1,mskv(i,j,k)+mskv(i,j+1,k))
          enddo
          enddo
        enddo
c
      enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_common_hyb:'
        endif
        call chksummsk(rig,ip,kk+1,'rig')
        call chksummsk(up,ip,kk,'up')
        call chksummsk(vp,ip,kk,'vp')
      endif
c
      end subroutine difest_common_hyb
c
      subroutine difest_isobml(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate diffusivities for eddy-induced transport, layer-wise
c --- diffusion and vertical diffusion
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      integer i,j,k,l,kn
c
c --- ------------------------------------------------------------------
c --- update halos of various fields
c --- ------------------------------------------------------------------
c
      call xctilr(u, 1,2*kk, 2,2, halo_uv)
      call xctilr(v, 1,2*kk, 2,2, halo_vv)
      call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
c
c --- ------------------------------------------------------------------
c --- Update layer interface pressure.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=-2,jj+3
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(-2,ifp(j,l)),min(ii+3,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Estimate friction velocity cubed.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          ustar3(i,j)=ustar(i,j)**3
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Estimate energy input by near-inertial waves.
      call niw_ke_tendency(m,n,mm,nn,k1m,k1n)
c
c --- Obtain common fields for the estimation of lateral and vertical
c --- diffusivities diapycnal diffusivities.
      call difest_common_iso(m,n,mm,nn,k1m,k1n)
c
c --- Estimate vertical diffusivity.
      call difest_vertical_iso(m,n,mm,nn,k1m,k1n)
c
c --- Estimate diffusivities for eddy-induced transport and layer-wise
c --- diffusion.
      call difest_lateral_iso(m,n,mm,nn,k1m,k1n)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_isobml:'
        endif
        call chksummsk(ustar3,ip,1,'ustar3')
      endif
c
      end subroutine difest_isobml
c
      subroutine difest_lateral_hybrid(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate diffusivities for eddy-induced transport, layer-wise
c --- diffusion and vertical diffusion
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      integer i,j,k,l,kn
c
c --- ------------------------------------------------------------------
c --- update halos of various fields
c --- ------------------------------------------------------------------
c
      call xctilr(u, 1,2*kk, 2,2, halo_uv)
      call xctilr(v, 1,2*kk, 2,2, halo_vv)
      call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
c
c --- ------------------------------------------------------------------
c --- Update layer interface pressure.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=-2,jj+3
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(-2,ifp(j,l)),min(ii+3,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Estimate friction velocity cubed.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          ustar3(i,j)=ustar(i,j)**3
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Obtain common fields for the estimation of lateral and vertical
c --- diffusivities diapycnal diffusivities.
      call difest_common_hyb(m,n,mm,nn,k1m,k1n)
c
c --- Estimate diffusivities for eddy-induced transport and layer-wise
c --- diffusion.
      call difest_lateral_hyb(m,n,mm,nn,k1m,k1n)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_lateral_hybrid:'
        endif
        call chksummsk(ustar3,ip,1,'ustar3')
      endif
c
      end subroutine difest_lateral_hybrid
c
      subroutine difest_vertical_hybrid(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate diffusivities for eddy-induced transport, layer-wise
c --- diffusion and vertical diffusion
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      integer i,j,k,l,kn
c
c --- ------------------------------------------------------------------
c --- update halos of various fields
c --- ------------------------------------------------------------------
c
      call xctilr(u(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
      call xctilr(v(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)
c
c --- ------------------------------------------------------------------
c --- Update layer interface pressure.
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=-2,jj+3
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(-2,ifp(j,l)),min(ii+3,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Obtain common fields for the estimation of lateral and vertical
c --- diffusivities diapycnal diffusivities.
      call difest_common_hyb(m,n,mm,nn,k1m,k1n)
c
c --- Estimate vertical diffusivities..
      call difest_vertical_hyb(m,n,mm,nn,k1m,k1n)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_vertical_hybrid:'
        endif
      endif
c
      end subroutine difest_vertical_hybrid
c
      subroutine difest_vertical_hyb(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate layer diapycnal, diffusivities for hybrid
c --- coordinates
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(kdm+1) :: rig_i
      integer i,j,k,l,kn
      real q
c
      type(CVMix_tidal_params_type)   :: CVMix_tidal_params
      real, dimension(kdm+1) :: depth_int
      real, dimension(kdm+1) :: Kv_col, Kd_col ! background visc/diff
      real, dimension(kdm+1) :: Kv_shr, Kd_shr ! shear driven visc/diff
      real, dimension(kdm+1) :: Kv_conv, Kd_conv ! convection visc/diff
      real, dimension(kdm+1) :: vert_dep       ! vertical deposition
      real, dimension(kdm+1) :: Kv_tidal, Kd_tidal ! tidal viscosity,diffusivity
      real, dimension(kdm+1) :: Kv_kpp, Kt_kpp, Ks_kpp ! vertical viscosity,diffusivity temp/salt
      real, dimension(kdm+1) :: iFaceHeight    ! Height of interfaces [m]
      real, dimension(kdm+1) :: bvfsq_i, bvf_i ! N2, N at interfaces
      real, dimension(kdm) :: cellHeight       ! Height of cell centers [m]
      real, dimension(kdm) :: rho_zeros, rho_lwr  ! dummy vars for convection
      real, dimension(kdm) :: rho_1d           ! 1D density at the layer center
      real, dimension(kdm) :: Ws_1d            ! Profile of vertical velocity scale for scalars [m s-1]
      real, dimension(kdm) :: surfBuoyFlux2
      real, dimension(kdm) :: BulkRi_1d        ! Bulk Richardson number for each layer
      real, dimension(kdm) :: deltaU2          ! square of delta U (shear) in denominator of Bulk Ri [m2 s-2]
      real, dimension(kdm) :: VT2              ! unresolved shear used for  Bulk Ri
      real, dimension(kdm) :: deltaRho         ! delta Rho [g/cm3] in numerator of Bulk Ri number
      real, dimension(kdm+1,2) :: nonLocalTrans  ! Non-local transport for heat/salt at interfaces [nondim]
      real :: surf_layer_ext, surfFricVel
      real :: surfBuoyFlux
      real :: delH, bvfbot, dps
      real :: dh, hcorr
      real :: Uk, Vk
      real :: surfU, surfV
      real :: surfHu, surfHv
      real :: surfTemp, surfSalt, surfRho
      real :: surfHtemp, surfHsalt
      real :: SLdepth_0d, hTot
      real :: Simmons_coeff, zBottomMinusOffset
      real :: bl1, bl2, bl3, bl4
      integer ki, ksfc, ktmp, kn1
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  hOBL
c
      surf_layer_ext = 0.1
      bl1 = 8e-5
      bl2 = 1.05e-4
      bl3 = 4.5e-3
      bl4 = 2500.0
c
c single column diffusivity
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c -- ------- CVMix variables computed below
c
             surfBuoyFlux2 = 0.0
             surfBuoyFlux = 0.0
c             Ws_1d = 0.0
             bvfbot    = 0.
             dps       = 0.
             depth_int = 0.
             hcorr     = 0.
             Kv_col    = 0.
             Kd_col    = 0.
             vert_dep  = 0.
             Kv_tidal  = 0.
             Kd_tidal  = 0.
             Kv_conv   = 0.
             Kd_conv   = 0.
             Kv_shr    = 0.
             Kd_shr    = 0.
             iFaceHeight = 0.
             cellHeight = 0.
             bvfsq_i   = 0.
             rho_lwr(:)= drho0
             rho_zeros(:)= 0.
             rho_1d    = 0.
             nonLocalTrans(:,:) = 0.0
             rig_i = 1.e8 !Initialize w/ large Richardson value
             Kv_kpp = 0.0
             Kt_kpp = 0.0
             Ks_kpp = 0.0
             do k=1,kk+1
                 Kv_kpp(k) = Kvisc_m(i,j,k)*1e-4
                 Kt_kpp(k) = Kdiff_t(i,j,k)*1e-4
                 Ks_kpp(k) = Kdiff_s(i,j,k)*1e-4
             enddo
             depth_int(1) = p(i,j,1)/onem
             iFaceHeight(1) = -depth_int(1)
             ! convert cm/s to m/s
             surfFricVel = ustar(i,j) * 1e-2
             ! convert cm2/s3 to m2/s3
             surfBuoyFlux = - buoyfl(i,j,1) * 1e-4
             do k=1,kk
                kn = k + nn
                kn1 = max(nn+1,kn-1)
                ! Old method to compute interface location, thicknesses
c                depth_int(k+1) = p(i,j,k+1)/onem
c                iFaceHeight(k+1) = -depth_int(k+1)
c                cellHeight(k) = 0.5*(iFaceHeight(k+1) +
c     .                                  iFaceHeight(k))
                ! New method to compute interface location, thicknesses
                dh = dp(i,j,kn)/onem
                dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
                hcorr = min(dh - 1e-10, 0.) ! If inflating then hcorr<0
                dh = max(dh, 1e-10) ! Limit increment  dh>=min_thicknes
                cellHeight(k) = iFaceHeight(k) - 0.5 * dh
                iFaceHeight(k+1) = iFaceHeight(k) - dh
                depth_int(k+1) = -iFaceHeight(k+1)

                ! compute rho_1d at the interfaces
                rho_1d(k)=rho(p(i,j,k),temp(i,j,kn),saln(i,j,kn))

                ! find ksfc for cell where "surface layer" sits
                SLdepth_0d = surf_layer_ext*
     .           max(max(-cellHeight(k),-iFaceHeight(2)),
     .                minOBLdepth)
                ksfc = k
                do ki = 1,k
                   if (-1.0*iFaceHeight(ki+1) >= SLdepth_0d) then
                      ksfc = ki
                      exit
                   endif
                enddo
                surfHu    = 0.0
                surfHv    = 0.0
                surfHtemp = 0.0
                surfHsalt = 0.0
                hTot      = 0.0
                do ki = 1,ksfc
                   ktmp = ki+nn
                   ! SLdepth_0d can be between cell interfaces
                   delH = min( max(0.0, SLdepth_0d - hTot),
     .                             dp(i,j,ktmp)/onem )
                   ! surface layer thickness
                   hTot = hTot + delH
                   ! surface averaged fields
                   surfHtemp = surfHtemp + temp(i,j,ktmp)*delH
                   surfHsalt = surfHsalt + saln(i,j,ktmp)*delH
                   surfHu = surfHu+up(i,j,ki)*delH
                   surfHv = surfHv+vp(i,j,ki)*delH
                enddo
                surfTemp = surfHtemp / hTot
                surfSalt = surfHsalt / hTot
                surfU    = surfHu / hTot
                surfV    = surfHv / hTot
                surfRho  = rho(p(i,j,k),surfTemp,surfSalt)
                if (p(i,j,kk+1)-p(i,j,k) < epsil) then
                    deltaRho(k) = deltaRho(k-1)
                else
                    deltaRho(k) = rho_1d(k) - surfRho
                endif
                ! vertical shear between present layer and
                ! surface layer averaged surfU,surfV.
                ! C-grid average to get Uk and Vk on T-points.
                Uk = up(i,j,k) - surfU
                Vk = vp(i,j,k) - surfV
                deltaU2(k) = (Uk**2 + Vk**2)

                ! XXX: Temporary de-scaling of N2_int(i,:) into a
                ! temporary variable
                bvfsq_i(k) = bfsqi(i,j,k)
                bvf_i(k) = sqrt( max( bvfsq_i(k), 0.) )
c --- ------- Accumulate Brunt-Vaisala frequency in a region near the
c --- ------- bottom
                q=max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpnbav,p(i,j,k)))
                if (q.gt.0.) then
                   bvfbot=bvfbot+bvf_i(k)*q
                   dps=dps+q
                endif

c --- ------- Local gradient Richardson number
                rig_i(k)=rig(i,j,k)

                surfBuoyFlux2(k) = ( buoyfl(i,j,k+1) 
     .                             - buoyfl(i,j,1  )) * 1e-4
c
             enddo  ! k
             if(dps.gt.0.) bvfbot=bvfbot/dps
             ! convert cm2/s2 to m2/s2
             deltaU2 = deltaU2*1e-4

             ! bottom values for the Ri, N2, and N
             rig_i(kk+1) = rig_i(kk)
             bvfsq_i(kk+1) = bfsqi(i,j,kk+1)
             bvf_i(kk+1) = sqrt( max( bvfsq_i(kk+1), 0.) )

c -- ------- Background diapycnal mixing.
             if (bdmtyp.eq.1) then
c zw interface depths relative to the surface in m, must be positive.
               call CVMix_init_bkgnd(max_nlev=kk, zw = depth_int(:),
     .                   bl1 = bl1, bl2 = bl2, bl3 = bl3, bl4 = bl4,
     .                   prandtl = CVMix_glb_params%Prandtl)
               call CVMix_coeffs_bkgnd(Mdiff_out=Kv_col,
     .                          Tdiff_out=Kd_col, nlev=kk, max_nlev=kk)
             elseif (bdmtyp.eq.2) then
c --- --------- Type 2: Background diffusivity is a constant
             ! convert cm2/s2 to m2/s2
                Kv_col(:) = bdmc2*1e-4
                Kd_col(:) = bdmc2*1e-4
             else
                Kv_col(:) = 0.
                Kd_col(:) = 0.
             endif
             if (iwdflg.eq.1) then
               Kv_col=Kv_col*(1.+(iwdfac-1.)*ficem(i,j))
               Kd_col=Kd_col*(1.+(iwdfac-1.)*ficem(i,j))
             endif

c --- ------ Tidally driven diapycnal mixing
c
             if (tdmflg.eq.1) then
               call CVMix_init_tidal(
     .                CVmix_tidal_params_user=CVMix_tidal_params,
     .                mix_scheme='Simmons',
     .                efficiency=dmxeff, local_mixing_frac=tdmq)

               call CVMix_compute_Simmons_invariant(nlev=kk,
     .                energy_flux=twedon(i,j)*bvfbot*1e-3,
     .                rho=CVMix_glb_params%FreshWaterDensity,
     .                SimmonsCoeff = Simmons_coeff, VertDep = vert_dep,
     .                zw = iFaceHeight, zt = cellHeight,
     .                CVmix_tidal_params_user=CVMix_tidal_params)


               call CVMix_coeffs_tidal(Mdiff_out=Kv_tidal,
     .                Tdiff_out=Kd_tidal, Nsqr = bvfsq_i,
     .                OceanDepth = -iFaceHeight(kk+1),
     .                SimmonsCoeff = Simmons_coeff,
     .                vert_dep = vert_dep,
     .                nlev=kk, max_nlev=kk,
     .                cvmix_params = CVMix_glb_params,
     .                CVmix_tidal_params_user=CVMix_tidal_params)
             else
               Kd_tidal=0.
             endif

! Call to CVMix wrapper for computing interior mixing coefficients.
             call  CVMix_coeffs_shear(Mdiff_out=Kv_shr(:),
     .                                Tdiff_out=Kd_shr(:),
     .                                RICH=rig_i(:),
     .                                nlev=kk,
     .                                max_nlev=kk)


c --- ------ turbulent velocity scales w_s and w_m computed at the cell
c --- ------ centers.
             call CVMix_kpp_compute_turbulent_scales(
     .              surf_layer_ext, ! (in)  Normalized surface layer Cdepth; sigma = CS%surf_layer_ext
     .              -cellHeight,       ! (in)  Assume here that OBL depth [m] = -cellHeight(k)
     .              surfBuoyFlux2,     ! (in)  Buoyancy flux at surface [m2 s-3]
     .              surfFricVel,       ! (in)  Turbulent friction  velocity at surface [m s-1]
     .              w_s=Ws_1d,         ! (out) Turbulent velocity scale profile [m s-1]
     .              CVMix_kpp_params_user=KPP_params)

             ! Compute unresolved shear for CVMix
             VT2(:) = CVmix_kpp_compute_unresolved_shear(
     .                 zt_cntr=cellHeight,   ! Depth ofcell center [m]
     .                 ws_cntr=Ws_1d,        ! Turbulent velocity scale profile, at centers [m s-1]
     .                 N_iface=bvf_i,        ! Buoyancy frequency at the interface [s-1]
     .                 CVMix_kpp_params_user=KPP_params)

             ! Calculate Bulk Richardson number from eq (21) of LMD94
             BulkRi_1d = CVmix_kpp_compute_bulk_Richardson(
     .              zt_cntr = cellHeight,                    ! Depth of cell  center [m]
     .              delta_buoy_cntr=g*alpha0*deltaRho*1e-2,  ! Bulk buoyancy difference, Br-B(z) [m s-2]
     .              delta_Vsqr_cntr=deltaU2,                 ! Square of resolved velocity difference [m2 s-2]
     .              Vt_sqr_cntr=VT2(:),                      ! Unresolved shear [m2 s-2]
     .              ws_cntr=Ws_1d,                           ! Turbulent velocity scale profile [m s-1]
     .              N_iface=bvf_i,                           ! Buoyancy frequency at the interface [s-1]
     .              CVMix_kpp_params_user=KPP_params ) ! KPP parameters

             ! Compute OBL depth for KPP
             call CVMix_kpp_compute_OBL_depth(
     .              BulkRi_1d,              ! (in) Bulk Richardson number
     .              iFaceHeight,            ! (in)  Height of interfaces [m]
     .              OBLdepth(i,j),          ! (out) OBL depth [m]
     .              hOBL(i,j),              ! (out) level (+fraction) of OBL extent
     .              zt_cntr = cellHeight,   ! Depth of cell  center [m]
     .              surf_fric=surfFricVel,  ! (in)  Turbulent friction  velocity at surface [m s-1]
     .              surf_buoy=surfBuoyFlux, ! (in) Buoyancy flux at surface [m2 s-3]
     .              Coriolis=coriop(i,j),   ! (in) Coriolis parameter [s-1]
     .              CVMix_kpp_params_user=KPP_params ) ! KPP parameters

             ! Avoid KPP reaching bottom
             zBottomMinusOffset = iFaceHeight(kk+1)
     .                          + min(1.0,-0.1*iFaceHeight(kk+1))
             OBLdepth(i,j) = min(OBLdepth(i,j), -zBottomMinusOffset)
             ! no shallower than top layer
             OBLdepth(i,j) = max(OBLdepth(i,j), -iFaceHeight(2))
             ! no deeper than bottom
             OBLdepth(i,j) = min(OBLdepth(i,j), -iFaceHeight(kk+1))
             ! gets index of the level and interface above hbl
             hOBL(i,j) = CVMix_kpp_compute_kOBL_depth(iFaceHeight,
     .                               cellHeight,OBLdepth(i,j))

             ! gets index of the level and interface above hbl
             kOBL(i,j) = CVMix_kpp_compute_kOBL_depth(iFaceHeight,
     .                               cellHeight,OBLdepth(i,j))
c --- ------ Diapycnal mixing when local stability is weak
c --- ------ convection routine based on N2 not rho
c --- ------ make sure it is in metrics if stability depends on rho
             call CVMix_coeffs_conv(Mdiff_out=Kv_conv,
     .                        Tdiff_out=Kd_conv, Nsqr = bvfsq_i,
     .                        dens=rho_zeros,dens_lwr=rho_lwr,
     .                        nlev=kk, max_nlev=kk,
     .                        OBL_ind=kOBL(i,j))
             ! Do not apply mixing due to convection within the boundary layer
             do k = 1,kOBL(i,j)
                Kv_conv(k) = 0.0
                Kd_conv(k) = 0.0
             enddo

             ! total diffusivities without KPP
             Kv_kpp(:) = Kv_col(:)+Kv_conv(:)+Kv_shr(:)
             Kt_kpp(:) = Kd_col(:)+Kd_conv(:)+Kd_shr(:)+Kd_tidal(:)
             Ks_kpp(:) = Kd_col(:)+Kd_conv(:)+Kd_shr(:)+Kd_tidal(:)

             ! Buoyancy flux acting on the OBL
             surfBuoyFlux = ( buoyfl(i,j,kOBL+1)
     .                      - buoyfl(i,j,1  )) * 1e-4

             ! Compute KPP using CVMix
             call CVMix_coeffs_kpp(Kv_kpp(:),          ! (inout) Total viscosity [m2 s-1]
     .              Kt_kpp(:),                         ! (inout) Total temp diffusivity [m2 s-1]
     .              Ks_kpp(:),                         ! (inout) Total salt  diffusivity [m2 s-1]
     .              iFaceHeight,                       ! (in)  Height of interfaces [m]
     .              cellHeight,                        ! (in)   Height of level centers [m]
     .              Kv_kpp(:),                         ! (in) Original viscosity [m2 s-1]
     .              Kt_kpp(:),                         ! (in) Original temp diffusivity [m2 s-1]
     .              Ks_kpp(:),                         ! (in) Original salt  diffusivity [m2 s-1]
     .              OBLdepth(i,j),                     ! (in) OBL depth [m]
     .              hOBL(i,j),                         ! (in) level (+fraction) of OBL extent
     .              nonLocalTrans(:,1),                ! (out) Non-local heat transport [nondim]
     .              nonLocalTrans(:,2),                ! (out) Non-local salt transport [nondim]
     .              surfFricVel,                       ! (in)  Turbulent friction  velocity at surface [m s-1]
     .              surfBuoyFlux,                      ! (in) Buoyancy flux at surface [m2 s-3]
     .              kk,                                ! (in) Number of levels to compute coeffs for
     .              kk,                                ! (in) Number of levels in array shape
     .              CVMix_kpp_params_user=KPP_params ) ! KPP parameters

             ! gets index of the level and interface above hbl
             kOBL(i,j) = CVMix_kpp_compute_kOBL_depth(iFaceHeight,
     .                               cellHeight,OBLdepth(i,j))

c ---- ccc -------
             ! convert m2/s to cm2/s
             Kv_kpp = Kv_kpp*1e4
             Kt_kpp = Kt_kpp*1e4
             Ks_kpp = Ks_kpp*1e4
             Kv_kpp=max(nubmin,Kv_kpp)
             Kt_kpp=max(nubmin,Kt_kpp)
             Ks_kpp=max(nubmin,Ks_kpp)
             Kvisc_m(i,j,:) = Kv_kpp(:)
             Kdiff_t(i,j,:) = Kt_kpp(:)
             Kdiff_s(i,j,:) = Ks_kpp(:)
             t_ns_nonloc(i,j,:) = nonLocalTrans(:,1)
             s_nonloc(i,j,:) = nonLocalTrans(:,2)
             do k = 1, kk+1
                t_sw_nonloc(i,j,k) = max(t_sw_nonloc(i,j,k),
     .                                   nonLocalTrans(k,1))
             enddo

        enddo
        enddo
c end of single column
c
      enddo ! j-index
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_vertical_hyb:'
        endif
        call chksummsk(Kvisc_m,ip,kk+1,'Kvisc_m')
        call chksummsk(Kdiff_t,ip,kk+1,'Kdiff_t')
        call chksummsk(Kdiff_s,ip,kk+1,'Kdiff_s')
      endif
c
      end subroutine difest_vertical_hyb
c
      subroutine difest_lateral_hyb(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate layer interface, isopycnal, diffusivities for hybrid
c --- coordinates
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
c
      real, dimension(1-nbdy:idm+nbdy,kdm) :: egr
      real, dimension(1-nbdy:idm+nbdy) ::
     .  tup,pup,sup,cr,bcrrd,afeql,dps,egrs,egrup,dfints,udps,vdps,
     .  umlzon,urmse,cpse
      integer i,j,k,l,kn
      real q,plo,tlo,slo,rhisc,els,egrlo,umnsc,esfac
c
c --- Locate the range of layers to be considered in the computation of
c --- diffusivities.
      do j=0,jj+1
        do i=0,ii+1
          kmax(i,j)=0
        enddo
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmax(i,j)=1
          do k=2,kk
            kn=k+nn
            if (dp(i,j,kn).gt.dpbmin) kmax(i,j)=k
          enddo
        enddo
        enddo
      enddo
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          kfil(i,j)=kk+1
          do k=kk,2,-1
            if (p(i,j,k).gt.mlts(i,j)*onecm) kfil(i,j)=k
          enddo
        enddo
        enddo
      enddo
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,k,kn,q,tup,pup,sup,cr,plo,tlo,slo,bcrrd,
c$OMP+ afeql,dps,egrs,egr,egrup,egrlo,dfints,
c$OMP+ rhisc,els,udps,vdps,umlzon,urmse,cpse,umnsc,esfac)
      do j=1,jj
c
c ----- Compute the first baroclinic rossby radius of deformation using
c ----- the WKB approximation by Chelton at al. (1998).
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pup(i)=p(i,j,1)
          kn=1+nn
          tup(i)=temp(i,j,kn)
          sup(i)=saln(i,j,kn)
          cr(i)=0.
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (p(i,j,kk+1)-p(i,j,k+1).lt.epsil) then
              plo=p(i,j,kk+1)
            else
              plo=.5*(p(i,j,k)+p(i,j,k+1))
            endif
            tlo=temp(i,j,kn)
            slo=saln(i,j,kn)
            cr(i)=cr(i)
     .           +sqrt(max(0.,(rho(p(i,j,k),tlo,slo)
     .                        -rho(p(i,j,k),tup(i),sup(i)))
     .                        *(plo-pup(i))))
            pup(i)=plo
            tup(i)=tlo
            sup(i)=slo
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cr(i)=alpha0*cr(i)/pi
          bcrrd(i)=
     .      sqrt(cr(i)*cr(i)
     .           /max(coriop(i,j)*coriop(i,j)+2.*betafp(i,j)*cr(i),
     .                1.e-24))
          afeql(i)=max(abs(coriop(i,j)),sqrt(2.*betafp(i,j)*cr(i)))
        enddo
        enddo
c
c --- - Compute diffusivity weigth to reduce eddy diffusivity when the
c --- - Rossby radius is resolved by the grid.
        if     (edwmth_opt.eq.edwmth_smooth) then
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j)
     .                         +scpy(i,j)*scpy(i,j)))
            difwgt(i,j)=1./(1.+.25*q**4)
          enddo
          enddo
        elseif (edwmth_opt.eq.edwmth_step) then
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j)
     .                         +scpy(i,j)*scpy(i,j)))
            if (q.le.2.) then
              difwgt(i,j)=1.
            else
              difwgt(i,j)=0.
            endif
          enddo
          enddo
        endif
c
c --- ------------------------------------------------------------------
c --- - Compute layer interface and isopycnal diffusivities
c --- ------------------------------------------------------------------
c
        if (iidtyp.eq.1) then
c
c --- --- Type 1: Diffusivities are diffusive velocities multiplied by
c --- --- the local horizontal grid scale.
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=sqrt(scp2(i,j))
            difint(i,j,1)=thkdff*q
            difiso(i,j,1)=temdff*q
          enddo
          enddo
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              difint(i,j,k)=difint(i,j,1)
              difiso(i,j,k)=difiso(i,j,1)
            enddo
            enddo
          enddo
c
        else
c
c --- --- Type 2: Diffusivities are parameterized according to Eden and
c --- --- Greatbatch (2008).
c
c --- --- Eady growth rate.
          if (edsprs) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              egrs(i)=0.
              dps(i)=0.
            enddo
            enddo
          endif
          if     (edritp_opt.eq.edritp_shear) then
            do k=2,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .              kmax(i,j)-kfil(i,j).ge.1) then
                  egr(i,k)=afeql(i)
     .                     /sqrt(.5*(rig(i,j,k)+rig(i,j,k+1))+eggam)
                  if (edsprs) then
                    q=max(0.,min(p(i,j,kfil(i,j))+dpgrav,
     .                           p(i,j,k+1))-p(i,j,k))
                    dps(i)=dps(i)+q
                    egrs(i)=egrs(i)+egr(i,k)*q
                  endif
                endif
              enddo
              enddo
            enddo
          elseif (edritp_opt.eq.edritp_large_scale) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (kmax(i,j)-kfil(i,j).ge.1) then
                k=kfil(i,j)
                if     (kmax(i-1,j).ge.k.and.kmax(i+1,j).ge.k) then
                  q=.25*(nnslpx(i,j,k)+nnslpx(i+1,j,k))**2
                elseif (kmax(i-1,j).ge.k) then
                  q=nnslpx(i,j,k)**2
                elseif (kmax(i+1,j).ge.k) then
                  q=nnslpx(i+1,j,k)**2
                else
                  q=0.
                endif
                if     (kmax(i,j-1).ge.k.and.kmax(i,j+1).ge.k) then
                  q=q+.25*(nnslpy(i,j,k)+nnslpy(i,j+1,k))**2
                elseif (kmax(i,j-1).ge.k) then
                  q=q+nnslpy(i,j,k)**2
                elseif (kmax(i,j+1).ge.k) then
                  q=q+nnslpy(i,j+1,k)**2
                endif
                egrup(i)=sqrt(q)
              endif
            enddo
            enddo
            do k=2,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (kmax(i,j)-kfil(i,j).ge.1) then
                  if     (k.ge.kfil(i,j).and.k.lt.kmax(i,j)) then
                    if     (kmax(i-1,j).gt.k.and.kmax(i+1,j).gt.k) then
                      q=.25*(nnslpx(i,j,k+1)+nnslpx(i+1,j,k+1))**2
                    elseif (kmax(i-1,j).gt.k) then
                      q=nnslpx(i,j,k+1)**2
                    elseif (kmax(i+1,j).gt.k) then
                      q=nnslpx(i+1,j,k+1)**2
                    else
                      q=0.
                    endif
                    if     (kmax(i,j-1).gt.k.and.kmax(i,j+1).gt.k) then
                      q=q+.25*(nnslpy(i,j,k+1)+nnslpy(i,j+1,k+1))**2
                    elseif (kmax(i,j-1).gt.k) then
                      q=q+nnslpy(i,j,k+1)**2
                    elseif (kmax(i,j+1).gt.k) then
                      q=q+nnslpy(i,j+1,k+1)**2
                    endif
                    egrlo=sqrt(q)
                    egr(i,k)=.5*(egrup(i)+egrlo)
                    egrup(i)=egrlo
                    if (edsprs) then
                      q=max(0.,min(p(i,j,kfil(i,j))+dpgrav,
     .                             p(i,j,k+1))-p(i,j,k))
                      dps(i)=dps(i)+q
                      egrs(i)=egrs(i)+egr(i,k)*q
                    endif
                  elseif (k.eq.kmax(i,j)) then
                    egr(i,k)=egr(i,k-1)
                  endif
                endif
              enddo
              enddo
            enddo
          endif
          if (edsprs) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (dps(i).gt.0.) then
                egrs(i)=egrs(i)/dps(i)
              else
                egrs(i)=0.
              endif
            enddo
            enddo
          endif
c
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            difint(i,j,1)=egmndf
            dfints(i)=0.
            dps(i)=0.
          enddo
          enddo
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
c --- --------- Rhines scale.
                rhisc=egr(i,k)/max(1.e-24,betafp(i,j))
c
c --- --------- Eddy length scale.
                els=max(eglsmn,min(bcrrd(i),rhisc))
c
c --- --------- Temporary layer interface diffusivity.
                difint(i,j,k)=egc*egr(i,k)*els*els
c
c --- --------- Accumulate diffusivities in a region below the first
c --- --------- physical layer.
                q=max(0.,min(p(i,j,kfil(i,j))+dpdiav,
     .                       p(i,j,k+1))-p(i,j,k))
                dps(i)=dps(i)+q
                dfints(i)=dfints(i)+difint(i,j,k)*q
c
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
c
c --- --- Apply eddy diffusivity limiting, suppression when the Rossby
c --- --- radius is resolved by the grid, and suppression away from
c --- --- steering levels if requested.
c
c --- --- Eddy diffusivity modification of surface non-isopycnic
c --- --- layers.
c
          if (edsprs) then
c
c --- ----- Zonal mixed layer velocity.
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              udps(i)=0.
              vdps(i)=0.
            enddo
            enddo
            do k=1,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                q=max(0.,min(p(i,j,k+1),OBLdepth(i,j)*onem)-p(i,j,k))
                udps(i)=udps(i)+up(i,j,k)*q
                vdps(i)=vdps(i)+vp(i,j,k)*q
              enddo
              enddo
            enddo
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              umlzon(i)=(udps(i)*cosang(i,j)-vdps(i)*sinang(i,j))
     .                  /(OBLdepth(i,j)*onem)
            enddo
            enddo
          endif
c
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
            if (edsprs) then
c
c --- ------- RMS eddy velocity estimated from K = Gamma*u_rms*L, where
c --- ------- a mixing efficiency of Gamma = 0.35 is used (Klocker and
c --- ------- Abernathey, 2014).
              rhisc=egrs(i)/max(1.e-24,betafp(i,j))
              els=max(eglsmn,min(bcrrd(i),rhisc))
              urmse(i)=2.86*egc*egrs(i)*els
c
c --- ------- Zonal eddy phase speed minus zonal barotropic velocity
c --- ------- with a lower bound of -20 cm s-1.
              cpse(i)=max(-20.,-betafp(i,j)*bcrrd(i)**2)
c
            endif
c
            if (dps(i).gt.0.) then
c
              if (edsprs) then
c
c --- --------- Zonal mixed layer velocity minus eddy phase speed. Note
c --- --------- that only the baroclinic component is used since the
c --- --------- barotropic velocity is subtracted from the estimate of
c --- --------- eddy phase speed.
                umnsc=umlzon(i)-cpse(i)
c
c --- --------- Eddy mixing suppresion factor where lower bounds of
c --- --------- zonal velocity minus eddy phase speed and absolute value
c --- --------- of RMS eddy velocity is set to -20 cm s-1 and 5 cm s-1,
c --- --------- respectively.
                esfac=1./(1.+4.*(umnsc/max(5.,abs(urmse(i))))**2)
c
              else
                esfac=1.
              endif
c
              dfints(i)=dfints(i)/dps(i)
              dfints(i)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,dfints(i)*difwgt(i,j)*esfac))
            else
              dfints(i)=egmndf
            endif
          enddo
          enddo
c
c --- --- Eddy diffusivity modification of isopycnic layers.
          do k=2,kk
            kn=k+nn
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
                if (edsprs) then
c
c --- ----------- Zonal velocity minus eddy phase speed.
                  umnsc=up(i,j,k)*cosang(i,j)-vp(i,j,k)*sinang(i,j)
     .                 -cpse(i)
c
c --- ----------- Eddy mixing suppresion factor.
                  esfac=1./(1.+4.*(umnsc/max(5.,abs(urmse(i))))**2)
c
                else
                  esfac=1.
                endif
c
                difint(i,j,k)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,difint(i,j,k)*difwgt(i,j)*esfac))
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
c
c --- --- Set isopycnal tracer diffusivity proportional to the layer
c --- --- interface diffusivity by the factor EGIDFQ.
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.lt.kfil(i,j)) then
                difint(i,j,k)=dfints(i)
              endif
              difiso(i,j,k)=difint(i,j,k)*egidfq
            enddo
            enddo
          enddo
c
        endif
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_lateral_hyb:'
        endif
        call chksummsk(difint,ip,kk,'difint')
        call chksummsk(difiso,ip,kk,'difiso')
      endif
c
      end subroutine difest_lateral_hyb
c
      subroutine difest_lateral_iso(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate layer interface, isopycnal, diffusivities for isopycnal
c --- coordinates
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy,kdm) :: egr
      real, dimension(1-nbdy:idm+nbdy) ::
     .  tup,pup,sup,cr,bcrrd,afeql,dps,egrs,egrup,dfints,urmse,cpse
      integer i,j,k,l,kn
      real q,plo,tlo,slo,rhisc,els,egrlo,umnsc,esfac
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,k,kn,q,tup,pup,sup,cr,plo,tlo,slo,bcrrd,
c$OMP+ afeql,dps,egrs,egr,egrup,egrlo,dfints,
c$OMP+ rhisc,els,urmse,cpse,umnsc,esfac)
      do j=1,jj
c
c ----- Compute the first baroclinic rossby radius of deformation using
c ----- the WKB approximation by Chelton at al. (1998).
c ----- !!! Could include top layer in computation !!!
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pup(i)=.5*(3.*p(i,j,3)-p(i,j,min(kk,kfpla(i,j,n))+1))
          kn=2+nn
          tup(i)=temp(i,j,kn)
          sup(i)=saln(i,j,kn)
          cr(i)=0.
        enddo
        enddo
        do k=3,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfpla(i,j,n)) then
              if (p(i,j,kk+1)-p(i,j,k+1).lt.epsil) then
                plo=p(i,j,kk+1)
              else
                plo=.5*(p(i,j,k)+p(i,j,k+1))
              endif
              tlo=temp(i,j,kn)
              slo=saln(i,j,kn)
              cr(i)=cr(i)
     .             +sqrt(max(0.,(rho(p(i,j,k),tlo,slo)
     .                          -rho(p(i,j,k),tup(i),sup(i)))
     .                          *(plo-pup(i))))
              pup(i)=plo
              tup(i)=tlo
              sup(i)=slo
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cr(i)=alpha0*cr(i)/pi
          bcrrd(i)=
     .      sqrt(cr(i)*cr(i)
     .           /max(coriop(i,j)*coriop(i,j)+2.*betafp(i,j)*cr(i),
     .                1.e-24))
          afeql(i)=max(abs(coriop(i,j)),sqrt(2.*betafp(i,j)*cr(i)))
        enddo
        enddo
c
c --- - Compute diffusivity weigth to reduce eddy diffusivity when the
c --- - Rossby radius is resolved by the grid.
        if     (edwmth_opt.eq.edwmth_smooth) then
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j)
     .                         +scpy(i,j)*scpy(i,j)))
            difwgt(i,j)=1./(1.+.25*q**4)
          enddo
          enddo
        elseif (edwmth_opt.eq.edwmth_step) then
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j)
     .                         +scpy(i,j)*scpy(i,j)))
            if (q.le.2.) then
              difwgt(i,j)=1.
            else
              difwgt(i,j)=0.
            endif
          enddo
          enddo
        endif
c
c --- ------------------------------------------------------------------
c --- - Compute layer interface and isopycnal diffusivities
c --- ------------------------------------------------------------------
c
        if (iidtyp.eq.1) then
c
c --- --- Type 1: Diffusivities are diffusive velocities multiplied by
c --- --- the local horizontal grid scale.
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=sqrt(scp2(i,j))
            difint(i,j,1)=thkdff*q
            difiso(i,j,1)=temdff*q
          enddo
          enddo
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              difint(i,j,k)=difint(i,j,1)
              difiso(i,j,k)=difiso(i,j,1)
            enddo
            enddo
          enddo
c
        else
c
c --- --- Type 2: Diffusivities are parameterized according to Eden and
c --- --- Greatbatch (2008).
c
c --- --- Eady growth rate.
          if (edsprs) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              egrs(i)=0.
              dps(i)=0.
            enddo
            enddo
          endif
          if     (edritp_opt.eq.edritp_shear) then
            do k=2,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .              kmax(i,j)-kfil(i,j).ge.1) then
                  egr(i,k)=afeql(i)/sqrt(rig(i,j,k)+eggam)
                  if (edsprs) then
                    q=max(0.,min(p(i,j,kfil(i,j))+dpgrav,
     .                           p(i,j,k+1))-p(i,j,k))
                    dps(i)=dps(i)+q
                    egrs(i)=egrs(i)+egr(i,k)*q
                  endif
                endif
              enddo
              enddo
            enddo
          elseif (edritp_opt.eq.edritp_large_scale) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (kmax(i,j)-kfil(i,j).ge.1) then
                k=kfil(i,j)
                if     (kmax(i-1,j).ge.k.and.kmax(i+1,j).ge.k) then
                  q=.25*(nnslpx(i,j,k)+nnslpx(i+1,j,k))**2
                elseif (kmax(i-1,j).ge.k) then
                  q=nnslpx(i,j,k)**2
                elseif (kmax(i+1,j).ge.k) then
                  q=nnslpx(i+1,j,k)**2
                else
                  q=0.
                endif
                if     (kmax(i,j-1).ge.k.and.kmax(i,j+1).ge.k) then
                  q=q+.25*(nnslpy(i,j,k)+nnslpy(i,j+1,k))**2
                elseif (kmax(i,j-1).ge.k) then
                  q=q+nnslpy(i,j,k)**2
                elseif (kmax(i,j+1).ge.k) then
                  q=q+nnslpy(i,j+1,k)**2
                endif
                egrup(i)=sqrt(q)
              endif
            enddo
            enddo
            do k=2,kk
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (kmax(i,j)-kfil(i,j).ge.1) then
                  if     (k.ge.kfil(i,j).and.k.lt.kmax(i,j)) then
                    if     (kmax(i-1,j).gt.k.and.kmax(i+1,j).gt.k) then
                      q=.25*(nnslpx(i,j,k+1)+nnslpx(i+1,j,k+1))**2
                    elseif (kmax(i-1,j).gt.k) then
                      q=nnslpx(i,j,k+1)**2
                    elseif (kmax(i+1,j).gt.k) then
                      q=nnslpx(i+1,j,k+1)**2
                    else
                      q=0.
                    endif
                    if     (kmax(i,j-1).gt.k.and.kmax(i,j+1).gt.k) then
                      q=q+.25*(nnslpy(i,j,k+1)+nnslpy(i,j+1,k+1))**2
                    elseif (kmax(i,j-1).gt.k) then
                      q=q+nnslpy(i,j,k+1)**2
                    elseif (kmax(i,j+1).gt.k) then
                      q=q+nnslpy(i,j+1,k+1)**2
                    endif
                    egrlo=sqrt(q)
                    egr(i,k)=.5*(egrup(i)+egrlo)
                    egrup(i)=egrlo
                    if (edsprs) then
                      q=max(0.,min(p(i,j,kfil(i,j))+dpgrav,
     .                             p(i,j,k+1))-p(i,j,k))
                      dps(i)=dps(i)+q
                      egrs(i)=egrs(i)+egr(i,k)*q
                    endif
                  elseif (k.eq.kmax(i,j)) then
                    egr(i,k)=egr(i,k-1)
                  endif
                endif
              enddo
              enddo
            enddo
          endif
          if (edsprs) then
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (dps(i).gt.0.) then
                egrs(i)=egrs(i)/dps(i)
              else
                egrs(i)=0.
              endif
            enddo
            enddo
          endif
c
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            difint(i,j,1)=egmndf
            dfints(i)=0.
            dps(i)=0.
          enddo
          enddo
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
c --- --------- Rhines scale.
                rhisc=egr(i,k)/max(1.e-24,betafp(i,j))
c
c --- --------- Eddy length scale.
                els=max(eglsmn,min(bcrrd(i),rhisc))
c
c --- --------- Temporary layer interface diffusivity.
                difint(i,j,k)=egc*egr(i,k)*els*els
c
c --- --------- Accumulate diffusivities in a region below the first
c --- --------- physical layer.
                q=max(0.,min(p(i,j,kfil(i,j))+dpdiav,
     .                       p(i,j,k+1))-p(i,j,k))
                dps(i)=dps(i)+q
                dfints(i)=dfints(i)+difint(i,j,k)*q
c
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
c
c --- --- Apply eddy diffusivity limiting, suppression when the Rossby
c --- --- radius is resolved by the grid, and suppression away from
c --- --- steering levels if requested.
c
c --- --- Eddy diffusivity modification of surface non-isopycnic
c --- --- layers.
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
            if (edsprs) then
c
c --- ------- RMS eddy velocity estimated from K = Gamma*u_rms*L, where
c --- ------- a mixing efficiency of Gamma = 0.35 is used (Klocker and
c --- ------- Abernathey, 2014).
              rhisc=egrs(i)/max(1.e-24,betafp(i,j))
              els=max(eglsmn,min(bcrrd(i),rhisc))
              urmse(i)=2.86*egc*egrs(i)*els
c
c --- ------- Zonal eddy phase speed minus zonal barotropic velocity
c --- ------- with a lower bound of -20 cm s-1.
              cpse(i)=max(-20.,-betafp(i,j)*bcrrd(i)**2)
c
            endif
c
            if (dps(i).gt.0.) then
c
              if (edsprs) then
c
c --- --------- Zonal mixed layer velocity minus eddy phase speed. Note
c --- --------- that only the baroclinic component is used since the
c --- --------- barotropic velocity is subtracted from the estimate of
c --- --------- eddy phase speed.
                if     (ip(i-1,j)+ip(i+1,j).eq.2) then
                  q=.5*((u(i  ,j,1+nn)*dpu(i  ,j,1+nn)
     .                  +u(i  ,j,2+nn)*dpu(i  ,j,2+nn))
     .                  /(dpu(i  ,j,1+nn)+dpu(i  ,j,2+nn))
     .                 +(u(i+1,j,1+nn)*dpu(i+1,j,1+nn)
     .                  +u(i+1,j,2+nn)*dpu(i+1,j,2+nn))
     .                  /(dpu(i+1,j,1+nn)+dpu(i+1,j,2+nn)))
                elseif (ip(i-1,j).eq.1) then
                  q=(u(i  ,j,1+nn)*dpu(i  ,j,1+nn)
     .              +u(i  ,j,2+nn)*dpu(i  ,j,2+nn))
     .              /(dpu(i  ,j,1+nn)+dpu(i  ,j,2+nn))
                elseif (ip(i+1,j).eq.1) then
                  q=(u(i+1,j,1+nn)*dpu(i+1,j,1+nn)
     .              +u(i+1,j,2+nn)*dpu(i+1,j,2+nn))
     .              /(dpu(i+1,j,1+nn)+dpu(i+1,j,2+nn))
                else
                  q=0.
                endif
                umnsc=q*cosang(i,j)
                if     (ip(i,j-1)+ip(i,j+1).eq.2) then
                  q=.5*((v(i,j  ,1+nn)*dpv(i,j  ,1+nn)
     .                  +v(i,j  ,2+nn)*dpv(i,j  ,2+nn))
     .                  /(dpv(i,j  ,1+nn)+dpv(i,j  ,2+nn))
     .                 +(v(i,j+1,1+nn)*dpv(i,j+1,1+nn)
     .                  +v(i,j+1,2+nn)*dpv(i,j+1,2+nn))
     .                  /(dpv(i,j+1,1+nn)+dpv(i,j+1,2+nn)))
                elseif (ip(i,j-1).eq.1) then
                  q=(v(i,j  ,1+nn)*dpv(i,j  ,1+nn)
     .              +v(i,j  ,2+nn)*dpv(i,j  ,2+nn))
     .              /(dpv(i,j  ,1+nn)+dpv(i,j  ,2+nn))
                elseif (ip(i,j+1).eq.1) then
                  q=(v(i,j+1,1+nn)*dpv(i,j+1,1+nn)
     .              +v(i,j+1,2+nn)*dpv(i,j+1,2+nn))
     .              /(dpv(i,j+1,1+nn)+dpv(i,j+1,2+nn))
                else
                  q=0.
                endif
                umnsc=umnsc-q*sinang(i,j)-cpse(i)
c
c --- --------- Eddy mixing suppresion factor where lower bounds of
c --- --------- zonal velocity minus eddy phase speed and absolute value
c --- --------- of RMS eddy velocity is set to -20 cm s-1 and 5 cm s-1,
c --- --------- respectively.
                esfac=1./(1.+4.*(umnsc/max(5.,abs(urmse(i))))**2)
c
              else
                esfac=1.
              endif
c
              dfints(i)=dfints(i)/dps(i)
              dfints(i)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,dfints(i)*difwgt(i,j)*esfac))
            else
              dfints(i)=egmndf
            endif
          enddo
          enddo
c
c --- --- Eddy diffusivity modification of isopycnic layers.
          do k=2,kk
            kn=k+nn
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .            kmax(i,j)-kfil(i,j).ge.1) then
c
                if (edsprs) then
c
c --- ----------- Zonal velocity minus eddy phase speed.
                  umnsc=
     .              (msku(i,j,k)*u(i,j,kn)+msku(i+1,j,k)*u(i+1,j,kn))
     .              /max(1,msku(i,j,k)+msku(i+1,j,k))*cosang(i,j)
     .             -(mskv(i,j,k)*v(i,j,kn)+mskv(i,j+1,k)*v(i,j+1,kn))
     .              /max(1,mskv(i,j,k)+mskv(i,j+1,k))*sinang(i,j)
     .             -cpse(i)
c
c --- ----------- Eddy mixing suppresion factor.
                  esfac=1./(1.+4.*(umnsc/max(5.,abs(urmse(i))))**2)
c
                else
                  esfac=1.
                endif
c
                difint(i,j,k)=
     .            min(difmxp(i,j),egmxdf,
     .                max(egmndf,difint(i,j,k)*difwgt(i,j)*esfac))
              else
                difint(i,j,k)=difint(i,j,k-1)
              endif
            enddo
            enddo
          enddo
c
c --- --- Set isopycnal tracer diffusivity proportional to the layer
c --- --- interface diffusivity by the factor EGIDFQ.
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k.lt.kfil(i,j)) then
                difint(i,j,k)=dfints(i)
              endif
              difiso(i,j,k)=difint(i,j,k)*egidfq
            enddo
            enddo
          enddo
c
        endif
      enddo
c$OMP END PARALLEL DO
c
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_lateral_iso:'
        endif
        call chksummsk(difint,ip,kk,'difint')
        call chksummsk(difiso,ip,kk,'difiso')
      endif
c
      end subroutine difest_lateral_iso
c
      subroutine difest_vertical_iso(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- estimate diapycnal diffusivities for isopycnal model
c --- ------------------------------------------------------------------
c
      integer m,n,mm,nn,k1m,k1n
c
      real, dimension(1-nbdy:idm+nbdy,kdm) :: bvfsq,bvf
      real, dimension(1-nbdy:idm+nbdy) :: bvfbot,dps,dfddsu,dfddsl
      integer i,j,k,l,kn
      real q,nus,nub,nut,nuls,vsf,nusm,ust,mols,h,sg,zeta,phis,ws
c
#if defined(TRC) && defined(TKE)
      real gls_c3,tke_prod,tke_buoy,tke_epsilon,Ls_unlmt,Ls_lmt,tke_Q,
     .     Gm,Gh,Sm,Sh,cff,ql
# ifdef GLS
      real gls_prod,gls_buoy,gls_diss,gls_Q
# endif
#endif
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,k,kn,q,bvfbot,dps,bvfsq,bvf,dfddsu,dfddsl,nub,nus,ust,vsf,
c$OMP+ nut,nuls,nusm,mols,h,sg,zeta,phis,ws
#if defined(TRC) && defined(TKE)
c$OMP+ ,gls_c3,tke_epsilon,tke_prod,tke_buoy,tke_Q,Ls_unlmt,Ls_lmt,Gh,
c$OMP+ Gm,cff,Sm,Sh,ql
# ifdef GLS
c$OMP+ ,gls_prod,gls_buoy,gls_diss,gls_Q
# endif
#endif
c$OMP+ )
      do j=1,jj
c
c ----- Compute Brunt-Vaisala frequency.
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          bvfbot(i)=0.
          dps(i)=0.
        enddo
        enddo
        do k=4,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .          kmax(i,j)-kfil(i,j).ge.1) then
c
c --- ------- Brunt-Vaisala frequency squared
              bvfsq(i,k)=g*g*max(drhomn,drhol(i,j,k))
     .                   /max(epsil,dp(i,j,kn))
c
c --- ------- Brunt-Vaisala frequency
              bvf(i,k)=sqrt(bvfsq(i,k))
c
#if defined(TRC) && defined(TKE)
              if (dp(i,j,kn).gt.dpbmin) then
                Buoy(i,j,k)=-difdia(i,j,k)*bvfsq(i,k)
                h=max(onem,dp(i,j,kn))*alpha0/g
c               h=max(onem*1e-8,dp(i,j,kn))*alpha0/g
c               h=max(onemm,dp(i,j,kn))*alpha0/g
                Shear2(i,j,k)=max(1.e-9,du2l(i,j,k))/(h*h)
                Prod(i,j,k)=difdia(i,j,k)*Pr_t*Shear2(i,j,k)
              else
                Buoy(i,j,k)=0.
                Shear2(i,j,k)=1.e-9
                Prod(i,j,k)=0.
              endif
#endif
c
c --- ------- Accumulate Brunt-Vaisala frequency in a region near the
c --- ------- bottom
              q=max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpnbav,p(i,j,k)))
              if (q.gt.0.) then
                bvfbot(i)=bvfbot(i)+bvf(i,k)*q
                dps(i)=dps(i)+q
              endif
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i).gt.0.) then
            bvfbot(i)=bvfbot(i)/dps(i)
          endif
        enddo
        enddo
c
c --- ------------------------------------------------------------------
c --- - Compute diapycnal diffusivity.
c --- ------------------------------------------------------------------
c
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          difdia(i,j,1)=nu0
          dfddsu(i)=0.
          dfddsl(i)=0.
          dps(i)=0.
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.ge.kfil(i,j).and.k.le.kmax(i,j).and.
     .          kmax(i,j)-kfil(i,j).ge.1) then
c
c --- ------- Background diapycnal mixing.
              if     (bdmtyp.eq.1) then
c
c --- --------- Type 1: Background diffusivity is a constant divided by
c --- --------- Brunt-Vaisala frequency.
                nub=bdmc1/bvf(i,k)
              elseif (bdmtyp.eq.2) then
c
c --- --------- Type 2: Background diffusivity is a constant
                nub=bdmc2
              else
                nub=0.
              endif
              if (iwdflg.eq.1) then
                nub=nub*(1.+(iwdfac-1.)*ficem(i,j))
              endif
c
c --- ------- Latitude dependency of background diapycnal mixing
              if (bdmldp.eq.1) then
                q=max(1.e-9,abs(coriop(i,j)))
                nub=nub*q/cori30*log(2.*bvf0/q)/log(2.*bvf0/cori30)
              endif
c
              nub=max(nubmin,nub)
c
#if !defined(TRC) || !defined(TKE)
c --- ------- Shear driven diapycnal mixing.
              if (rig(i,j,k).lt.ri0) then
c
c --- --------- Maximum diffusivity is increased near the bottom to
c --- --------- provide additional mixing of gravity currents.
                q=(p(i,j,kk+1)-p(i,j,k)+.5*dp(i,j,kn))
     .            /min(dpgc,.5*p(i,j,kk+1))
                q=max(0.,1.-q*q)
                q=q*q*q
                nus=q*nug0+(1.-q)*nus0
c
c --- --------- Parameterization of diffusivity as a function of local
c --- --------- gradient richardson number.
                q=rig(i,j,k)/ri0
                q=max(0.,1.-q*q)
                nus=nus*q*q*q
              else
                nus=0.
              endif
#else
              if (bvfsq(i,k).gt.0.) then  ! stable stratification
                gls_c3=gls_c3minus
              else                        ! unstable stratification
                gls_c3=gls_c3plus
              endif
#  ifndef GLS
              trc(i,j,kn,itrgls)=max((gls_c1*Prod(i,j,k)
     .                               +gls_c3*Buoy(i,j,k))/gls_c2,
     .                               gls_psi_min)
#  endif
              tke_epsilon=cmu_fac2*trc(i,j,kn,itrtke)**(1.5+gls_m/gls_n)
     .                    *trc(i,j,kn,itrgls)**(-1./gls_n)
              tke_prod=Prod(i,j,k)
              tke_buoy=Buoy(i,j,k)
              tke_Q=tke_epsilon/trc(i,j,kn,itrtke)
#  ifdef GLS
              gls_prod=(trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke))
     .                 *gls_c1*Prod(i,j,k)
              gls_buoy=(trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke))
     .                 *gls_c3*Buoy(i,j,k)
              gls_diss=(trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke))
     .                 *gls_c2*tke_epsilon
              gls_Q=gls_diss/trc(i,j,kn,itrgls)
#  endif
#  ifdef GLS
              if (gls_prod+gls_buoy.ge.0.) then
                trc(i,j,kn,itrgls)=
     .            (trc(i,j,kn,itrgls)+delt1*(gls_prod+gls_buoy))
     .            /(1.+delt1*gls_Q)
              else
                trc(i,j,kn,itrgls)=
     .            (trc(i,j,kn,itrgls)+delt1*gls_prod)
     .            /(1.+delt1*(gls_Q-(gls_buoy/trc(i,j,kn,itrgls))))
              endif
              trc(i,j,kn,itrgls)=max(trc(i,j,kn,itrgls),gls_psi_min)
              q=.56**(.5*gls_n)*gls_cmu0**gls_p
     .          *trc(i,j,kn,itrtke)**(gls_m+.5*gls_n)
     .          *bvf(i,k)**(-gls_n)
              if (gls_n.lt.0.) then
                trc(i,j,kn,itrgls)=max(trc(i,j,kn,itrgls),q)
              else
                trc(i,j,kn,itrgls)=min(trc(i,j,kn,itrgls),q)
              endif
#  endif
c
              tke_epsilon=cmu_fac2*trc(i,j,kn,itrtke)**(1.5+gls_m/gls_n)
     .                    *trc(i,j,kn,itrgls)**(-1./gls_n)
              tke_Q=tke_epsilon/trc(i,j,kn,itrtke)
c
              if (tke_prod+tke_buoy.ge.0.) then
                trc(i,j,kn,itrtke)=
     .            (trc(i,j,kn,itrtke)+delt1*(tke_prod+tke_buoy))
     .            /(1.+delt1*tke_Q)
              else
                trc(i,j,kn,itrtke)=
     .            (trc(i,j,kn,itrtke)+delt1*tke_prod)
     .            /(1.+delt1*(tke_Q-(tke_buoy/trc(i,j,kn,itrtke))))
                trc(i,j,kn,itrtke)=max(trc(i,j,kn,itrtke),tke_min)
              endif
c
c --- ------- Penetration of surface TKE below mixed layer.
              if (tkepf.gt.0.) then
                if (dp(i,j,kn).lt.epsil) then
                  q=exp(-p(i,j,k)/tkepls)
                else
                  q=tkepls*(exp(-p(i,j,k  )/tkepls)
     .                     -exp(-p(i,j,k+1)/tkepls))/dp(i,j,kn)
                endif
                trc(i,j,kn,itrtke)=trc(i,j,kn,itrtke)
     .                            +67.83*tkepf*q*ustar(i,j)**2
              endif
c
c --- ------- Set TKE and GLS to prescribed minimum values in surface
c --- ------- mixed layers and thin layers
              if (dp(i,j,kn).lt.epsil) then
                trc(i,j,kn,itrtke)=tke_min
                trc(i,j,kn,itrgls)=gls_psi_min
              endif
              trc(i,j,1+nn,itrtke)=tke_min
              trc(i,j,2+nn,itrtke)=tke_min
              trc(i,j,1+nn,itrgls)=gls_psi_min
              trc(i,j,2+nn,itrgls)=gls_psi_min
c
c --- ------- Bottom Boundary Conditions
              if (k.eq.kmax(i,j)) then
                ust=max(ustarb(i,j),ustmin)
                trc(i,j,kn,itrtke)=max(tke_min,(ust/gls_cmu0)**2)
#  ifdef GLS
                trc(i,j,kn,itrgls)=max(gls_psi_min,
     .                                 (gls_cmu0**(gls_p-2.*gls_m))
     .                                 *(ust**(2.*gls_m))
     .                                 *(kappa*1.e2)**gls_n)
#  endif
              endif
c
              Ls_unlmt=max(Ls_unlmt_min,
     .                     cmu_fac1*trc(i,j,kn,itrgls)**(gls_exp1)
     .                     *trc(i,j,kn,itrtke)**(-tke_exp1))

              if (bvfsq(i,k).gt.0.) then  ! stable stratification
c               Ls_lmt=min(Ls_unlmt,
c    .                     sqrt(.56*trc(i,j,kn,itrtke)
c    .                          /max(bvfsq(i,k),1.e-10)))

                Ls_lmt=min(Ls_unlmt,trc(i,j,kn,itrtke)**(-gls_m/gls_n)
     .                    *trc(i,j,kn,itrgls)**gls_n)
c               Ls_lmt=Ls_unlmt
              else                        ! unstable stratification
                Ls_lmt=Ls_unlmt
              endif
c
c --- ------- Compute nondimensional stability functions for tracers
c --- ------- (Sh) and momentum (Sm). Canuto-A
              Gh=min(gls_Gh0,-bvfsq(i,k)*Ls_lmt*Ls_lmt
     .                        /(2.*trc(i,j,kn,itrtke)))
              Gh=min(Gh,(Gh-(Gh-gls_Ghcri)**2)
     .                  /(Gh+gls_Gh0-2.*gls_Ghcri))
              Gh=max(Gh,gls_Ghmin)
              Gh=min(Gh,gls_Gh0)
c
c --- ------- Compute shear number.
              Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh**2))
     .           /(gls_b2-gls_b4*gls_fac6*Gh)
              Gm=min(Gm,Shear2(i,j,k)*Ls_lmt*Ls_lmt
     .                  /(2.*trc(i,j,kn,itrtke)))
c
c --- ------- Compute stability functions
              cff=gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm
     .           +gls_b3*gls_fac6**2*Gh**2-gls_b4*gls_fac6**2*Gh*Gm
     .           +gls_b5*gls_fac6**2*Gm*Gm
              Sm=(gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/cff
              Sh=(gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/cff
              Sm=max(Sm,0.)
              Sh=max(Sh,0.)
c
c --- ------- Relate Canuto stability to BLOM notation
              Sm=Sm*cmu_fac3/gls_cmu0**3
              Sh=Sh*cmu_fac3/gls_cmu0**3
c
              ql=sqrt2*(Ls_lmt)
     .           *sqrt(trc(i,j,kn,itrtke))
c             ql=sqrt2*.5*(Ls_lmt+L_scale(i,j,k))
c    .           *sqrt(trc(i,j,kn,itrtke))
c
c             nus=Sh*ql
c             nus=min(0.1*ql,4.05*nug0)
              nus=min(Sh*ql,4.05*nug0)
c             nus=Sh*(trc(i,j,k,itrtke)*trc(i,j,k,itrtke))
c    .            /trc(i,j,k,itrgls)
              L_scale(i,j,k)=max(Ls_lmt,Ls_unlmt_min)
#  ifdef GLS
c
c --- ------- Recompute gls based on limited length scale
              trc(i,j,kn,itrgls)=
     .          max(gls_cmu0**gls_p*trc(i,j,kn,itrtke)**gls_m
     .              *L_scale(i,j,k)**gls_n,gls_psi_min)
#  endif
#endif
c
c --- ------- Tidally driven diapycnal mixing
              if (tdmflg.eq.1) then
                q=.5*(tanh(4.*(abs(plat(i,j))-tdclat)/tddlat-2.)+1.)
                q=(1.-q)*tdmls0+q*tdmls1
                if (dp(i,j,kn).lt.epsil) then
                  vsf=exp(p(i,j,k)/q)/(q*(exp(p(i,j,kk+1)/q)-1.))
                else
                  vsf=(exp(p(i,j,k+1)/q)-exp(p(i,j,k)/q))
     .                /(dp(i,j,kn)*(exp(p(i,j,kk+1)/q)-1.))
                endif
                nut=g*tdmq*dmxeff*twedon(i,j)*bvfbot(i)*vsf/bvfsq(i,k)
              else
                nut=0.
              endif
c
c --- ------- Diapycnal mixing when local stability is weak
              if (drhol(i,j,k).lt.drho0) then
                q=drhol(i,j,k)/drho0
                q=max(0.,1.-q*q)
                nuls=nuls0*q*q*q
              else
               nuls=0.
              endif
c
c --- ------- Total diapycnal diffusivity.
              difdia(i,j,k)=nub+nus+nut+nuls
c
c --- ------- Accumulate diffusivities in a region below the first
c --- ------- physical layer
              q=max(0.,min(p(i,j,kfil(i,j))+dpddav,p(i,j,k+1))-p(i,j,k))
              dps(i)=dps(i)+q
              dfddsu(i)=dfddsu(i)+nub*q
              dfddsl(i)=dfddsl(i)+difdia(i,j,k)*q
c
            else
              difdia(i,j,k)=difdia(i,j,k-1)
#if defined(TRC) && defined(TKE)
c             trc(i,j,kn,itrtke)=tke_min
c             L_scale(i,j,k)=Ls_unlmt_min
              trc(i,j,kn,itrtke)=trc(i,j,kn-1,itrtke)
              L_scale(i,j,k)=L_scale(i,j,k-1)
#  ifdef GLS
c             trc(i,j,kn,itrgls)=gls_psi_min
              trc(i,j,kn,itrgls)=trc(i,j,kn-1,itrgls)
#  endif
#endif
            endif
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i).gt.0.) then
            dfddsu(i)=dfddsu(i)/dps(i)
            dfddsl(i)=dfddsl(i)/dps(i)
          else
            dfddsu(i)=nu0
            dfddsl(i)=nu0
          endif
        enddo
        enddo
        do k=2,kk-1
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.lt.kfil(i,j)) then
              if (k.gt.2.and.kfil(i,j).le.kk.and.
     .            p(i,j,min(kk,kfil(i,j)))-p(i,j,3).gt.epsil) then
                q=.5*(p(i,j,k+1)+p(i,j,k))
                difdia(i,j,k)=((q-p(i,j,3))*dfddsl(i)
     .                        +(p(i,j,kfil(i,j))-q)*dfddsu(i))
     .                        /(p(i,j,kfil(i,j))-p(i,j,3))
              else
                difdia(i,j,k)=dfddsu(i)
              endif
            endif
          enddo
          enddo
        enddo
c
c --- - Diapycnal diffusivity beneath mixed layer by dissipation of
c --- - energy originating from near-inertial waves.
        do k=2,kk-1
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k.le.kmax(i,j).and.kmax(i,j)-kfil(i,j).ge.1) then
              q=niwls
              if (k.eq.2.or.dp(i,j,kn).lt.epsil) then
                vsf=exp((p(i,j,3)-p(i,j,k+1))/q)
     .              /(q*(1.-exp((p(i,j,3)-p(i,j,kk+1))/q)))
              else
                vsf=(exp((p(i,j,3)-p(i,j,k  ))/q)
     .              -exp((p(i,j,3)-p(i,j,k+1))/q))
     .              /(dp(i,j,kn)*(1.-exp((p(i,j,3)-p(i,j,kk+1))/q)))
              endif
              nusm=g*niwgf*(1.-niwbf)*niwlf*dmxeff*idkedt(i,j)*vsf
     .             /(alpha0*bvfsq(i,max(k,kfil(i,j))))
              difdia(i,j,k)=difdia(i,j,k)+nusm
            endif
          enddo
          enddo
        enddo
c
c --- - Diffusivity at the lower interface of the top layer
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Lower bounded friction velocity
          ust=max(ustmin,ustar(i,j))
c
c --- --- Monin-Obukhov length scale
          mols=ust**3/(kappa*sign(max(abs(buoyfl(i,j,1)),bfeps),
     .                           -buoyfl(i,j,1)))
c
c --- --- Mixed layer thickness
          h=(p(i,j,3)-p(i,j,1))/onecm
c
c --- --- Dimensionless vertical coordinate in the boundary layer
          sg=(p(i,j,2)-p(i,j,1))/(p(i,j,3)-p(i,j,1))
c
c --- --- Velocity scale
          if (mols.lt.0.) then
            zeta=min(sleps,sg)*h/mols
            if (zeta.gt.zetas) then
              phis=(1.-16.*zeta)**(-1./2.)
            else
              phis=(as-cs*zeta)**(-1./3.)
            endif
          else
            zeta=sg*h/mols
            phis=1.+5.*zeta
          endif
          ws=kappa*ust/phis
c
          difdia(i,j,1)=h*ws*sg*(1.-sg)**2
        enddo
        enddo
c
      enddo
c$OMP END PARALLEL DO
c
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'difest_vertical_iso:'
        endif
        call chksummsk(idkedt,ip,1,'idkedt')
        call chksummsk(difdia,ip,kk,'difdia')
#if defined(TRC) && defined(TKE)
        call chksummsk(trc(1-nbdy,1-nbdy,1,itrtke),ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(trc(1-nbdy,1-nbdy,1,itrgls),ip,2*kk,'gls_psi')
#  endif
#endif
      endif
c
      end subroutine difest_vertical_iso
c
      end module mod_difest
