! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2020 Mats Bentsen, Lars Inge Enstad, Ingo Bethke,
!                         Mehmet Ilicak, Alok Kumar Gupta
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

c --- ------------------------------------------------------------------
c --- common blocks related to the dynamic/thermodynamic core
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .  u,v,           ! velocity components
     .  dp,            ! layer thickness
     .  dpold,         ! layer thickness at old time level
     .  dpu,dpv,       ! layer thickness at u- and v-points
     .  temp,          ! temperature
     .  saln,          ! salinity
     .  sigma,         ! potential density
     .  absvor,        ! absolute vorticity
     .  dpvor          ! layer pressure thickness used in vorticity computation
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     .  p,             ! interface pressure
     .  pu,pv,         ! interface pressure at u- and v-points
     .  phi            ! interface geopotential
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  sigmar,        ! reference potential density
     .  temmin,        ! minimum temperature allowed in an isopycnic layer
     .  dpuold,dpvold, ! layer thickness at u- and v-points at old time level
     .  told,          ! temperature at old time level
     .  sold           ! salinity at old time level
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  corioq,        ! coriolis parameter at q-point
     .  coriop,        ! coriolis parameter at p-point
     .  betafp,        ! latitudinal variation of the coriolis param. at p-point
     .  potvor         ! potential vorticity
c
      common /blom1/ u,v,dp,dpold,dpu,dpv,temp,saln,sigma,absvor,dpvor,
     .               p,pu,pv,phi,sigmar,temmin,dpuold,dpvold,told,sold,
     .               corioq,coriop,betafp,potvor
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .  uflx,vflx,     ! horizontal mass fluxes
     .  utflx,vtflx,   ! horizontal heat fluxes
     .  usflx,vsflx    ! horizontal salt fluxes
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .  umfltd,vmfltd, ! horizontal mass fluxes due to thickness diffusion
     .  utfltd,vtfltd, ! horizontal heat fluxes due to thickness diffusion
     .  utflld,vtflld, ! horizontal heat fluxes due to lateral diffusion
     .  usfltd,vsfltd, ! horizontal salt fluxes due to thickness diffusion
     .  usflld,vsflld  ! horizontal salt fluxes due to lateral diffusion
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3) ::
     .  ubflxs,vbflxs  ! barotropic mass flux sums
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  pb,            ! bottom pressure
     .  ubflx,vbflx,   ! barotropic mass fluxes
     .  pb_mn,         ! bottom pressure
     .  ubflx_mn,vbflx_mn, ! barotropic mass fluxes
     .  pbu,pbv,       ! bottom pressure at velocity points
     .  ub,vb,         ! barotropic velocity components
     .  ubflxs_p,vbflxs_p, ! predicted barotropic mass flux sums
     .  pvtrop         ! potential vorticity of barotropic flow
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  pb_p,              ! predicted bottom pressure
     .  pbu_p,pbv_p,       ! predicted bottom pressure at velocity points
     .  pvtrop_o,      ! potential vorticity of barotropic flow at old time lev.
     .  ubcors_p,vbcors_p, ! predicted sums of barotropic coriolis terms
     .  defor1,defor2, ! deformation components
     .  utotm,vtotm,   ! total (barotropic+baroclinic) ...
     .  utotn,vtotn,   ! ... velocities at 2 time levels
     .  uflux,vflux,   ! horizontal mass fluxes
     .  uflux2,vflux2, ! more mass fluxes
     .  uflux3,vflux3  ! more mass fluxes
c
      common /blom2/ uflx,vflx,utflx,vtflx,usflx,vsflx,umfltd,vmfltd,
     .               utfltd,vtfltd,utflld,vtflld,usfltd,vsfltd,usflld,
     .               vsflld,ubflxs,vbflxs,pb,ubflx,vbflx,pb_mn,
     .               ubflx_mn,vbflx_mn,pbu,pbv,ub,vb,ubflxs_p,vbflxs_p,
     .               pvtrop,pb_p,pbu_p,pbv_p,pvtrop_o,ubcors_p,
     .               vbcors_p,defor1,defor2,utotm,vtotm,utotn,vtotn,
     .               uflux,vflux,uflux2,vflux2,uflux3,vflux3
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) ::
     .  pgfx,pgfy      ! horizontal pressure gradient force
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  pgfxo,pgfyo    ! horizontal pressure gradient force at old time level
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  pgfxm,pgfym,   ! PGF terms not dependent on barotropic bottom pressure
     .  xixp,xixm,     ! PGF terms dependent on barotropic bottom pressure
     .  xiyp,xiym      ! PGF terms dependent on barotropic bottom pressure
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  pgfxm_o,pgfym_o,   ! PGF terms not dep. on barotr. bot. pres. at old ti.
     .  xixp_o,xixm_o,     ! PGF terms dep. on barotr. bot. pres. at old tim. l.
     .  xiyp_o,xiym_o,     ! PGF terms dep. on barotr. bot. pres. at old tim. l.
     .  util1,util2,   ! arrays for temporary storage
     .  util3,util4,   ! arrays for temporary storage
     .  scqx,scqy,     ! mesh size at q-points in x,y direction
     .  scpx,scpy,     ! mesh size at p-points in x,y direction
     .  scux,scuy,     ! mesh size at u-points in x,y direction
     .  scvx,scvy,     ! mesh size at v-points in x,y direction
     .  scq2,scp2,     ! grid box size at q- and p-points
     .  scu2,scv2,     ! grid box size at u- and v-points
     .  scq2i,scp2i,   ! inverses of scq2,scp2
     .  scuxi,scvyi,   ! inverses of scux,scvy
     .  scuyi,scvxi,   ! inverses of scuy,scvx
     .  umax,vmax,     ! maximum allowable velocities
     .  depths         ! water depth
c
      common /blom3/ pgfx,pgfy,pgfxo,pgfyo,pgfxm,pgfym,xixp,xixm,
     .               xiyp,xiym,pgfxm_o,pgfym_o,xixp_o,xixm_o,
     .               xiyp_o,xiym_o,util1,util2,util3,util4,scqx,scqy,
     .               scpx,scpy,scux,scuy,scvx,scvy,scq2,scp2,scu2,scv2,
     .               scq2i,scp2i,scuxi,scvyi,scuyi,scvxi,umax,vmax,
     .               depths
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  bfsqi,         ! interface buoyancy frequency squared
     .  bfsql,         ! layer buoyancy frequency squared
     .  bfsqf,         ! filtered interface buoyancy frequency squared
     .  nslpx,nslpy,   ! local neutral slope
     .  nnslpx,nnslpy, ! local neutral slope times buoyancy frequency
     .  difint,        ! layer interface diffusivity
     .  difiso,        ! isopycnal diffusivity
     .  difdia         ! diapycnal diffusivity
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::
     .  uml,vml        ! total mixed layer 1 velocity
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  umlres,vmlres  ! mixed layer velocity reservoar for temporal smoothing
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  uja,ujb,       ! velocities at lateral ...
     .  via,vib,       ! ... neighbor points
     .  difmxp,        ! maximum lateral diffusivity at p-points
     .  difmxq,        ! maximum lateral diffusivity at q-points
     .  difwgt,        ! eddy diffusivity weight
     .  sealv,         ! sea surface height
     .  surflx,        ! surface thermal energy flux
     .  surrlx,        ! surface relaxation thermal energy flux
     .  sswflx,        ! surface solar energy flux
     .  salflx,        ! surface salinity flux
     .  brnflx,        ! surface brine flux
     .  salrlx,        ! surface relaxation salinity flux
     .  taux,tauy,     ! surface stress components
     .  ustar,         ! surface friction velocity
     .  ustarb,        ! bottom friction velocity
     .  idkedt,        ! vertically integrated inertial kinetic energy tendency
     .  ustar3,        ! friction velocity cubed
     .  buoyfl,        ! surface buoyancy flux
     .  mtkeus,        ! mixed layer TKE tendency related to friction velocity
     .  mtkeni,        ! mixed layer TKE tendency related to near inertial mot.
     .  mtkebf,        ! mixed layer TKE tendency related to buoyancy forcing
     .  mtkers,        ! mixed layer TKE tendency related to eddy restratific.
     .  mtkepe,        ! mixed layer TKE tendency related to pot. energy change
     .  mtkeke,        ! mixed layer TKE tendency related to kin. energy change
     .  twedon,        ! tidal wave energy diffipation over buoyancy frequency
     .  pbrnda         ! brine plume pressure depth
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  kfpla          ! index of first physical layer
c
      common /blom4/ bfsqi,bfsql,bfsqf,nslpx,nslpy,nnslpx,nnslpy,
     .               difint,difiso,difdia,uml,vml,umlres,vmlres,
     .               uja,ujb,via,vib,difmxp,difmxq,difwgt,sealv,
     .               surflx,surrlx,sswflx,salflx,brnflx,salrlx,
     .               taux,tauy,ustar,ustarb,idkedt,ustar3,buoyfl,
     .               mtkeus,mtkeni,mtkebf,mtkers,mtkepe,mtkeke,
     .               twedon,pbrnda,kfpla
!!!!!!!!!!!!$acc declare device_resident(/blom4/)
c
      real time,delt1,dlt,area,avgbot
      integer nstep,nstep0,nstep1,nstep2,lstep
c
      common /varbls/ time,delt1,dlt,area,avgbot,
     .                nstep,nstep0,nstep1,nstep2,lstep
!!!!!!!!!!!!!!!!!!!$acc declare device_resident(/varbls/)
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'mdv2hi' = Laplacian diffusion velocity (cm/s) for momentum dissipation
c --- 'mdv2lo' = same as mdv2hi but used when Rossby radius is resolved
c --- 'mdv4hi' = Biharmonic diffusion velocity (cm/s) for momentum dissipation
c --- 'mdv4lo' = same as mdv4hi but used when Rossby radius is resolved
c --- 'mdc2hi' = Laplacian diffusivity (cm**2/s) for momentum dissipation
c --- 'mdc2lo' = same as mdc2hi but used when Rossby radius is resolved
c --- 'vsc2hi' = parameter used in deformation-dependent Laplacian viscosity
c --- 'vsc2lo' = same as vsc2hi but used when Rossby radius is resolved
c --- 'vsc4hi' = parameter used in deformation-dependent Biharmonic viscosity
c --- 'vsc4lo' = same as vsc4hi but used when Rossby radius is resolved
c --- slip = +1  for free-slip boundary cond., slip = -1  for non-slip cond.
c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
c --- 'cb'     = coefficient of quadratic bottom friction
c --- 'cwbdts' = coastal wave breaking damping resiprocal time scale (1/s)
c --- 'cwbdls' = coastal wave breaking damping length scale (m)
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'wpgf'   = weight for time averaging of pressure gradient force
c --- 'mltmin' = minimum mixed-layer thickness (m)
c --- 'thktop' = thickness of top layer (m)
c --- 'thkbot' = thickness of bottom boundary layer (pressure units)
c --- 'egc'    = the parameter c in the Eden and Greatbatch (2008)
c ---            parameterization
c --- 'eggam'  = the parameter gamma in the Eden and Greatbatch (2008)
c ---            parameterization []
c --- 'eglsmn' = minimum eddy length scale in Eden and Greatbatch (2008)
c ---            parameterization [cm]
c --- 'egmndf' = minimum diffusivity in the Eden and Greatbatch (2008)
c ---            parameterization [cm**2/s]
c --- 'egmxdf' = maximum diffusivity in the Eden and Greatbatch (2008)
c ---            parameterization [cm**2/s]
c --- 'egidfq' = factor relating the isopycnal diffusivity to the layer
c ---            interface diffusivity in the Eden and Greatbatch (2008)
c ---            parameterization. egidfq=difint/difiso
c --- 'ri0'    = critical gradient richardson number for shear driven
c ---            vertical mixing
c --- 'rm0'    = efficiency factor of wind TKE generation in the
c ---            Oberhuber (1993) TKE closure
c --- 'rm5'    = efficiency factor of TKE generation by momentum
c ---            entrainment in the Oberhuber (1993) TKE closure
c --- 'ce'     = efficiency factor for the restratification by mixed
c ---            layer eddies (Fox-Kemper et al., 2008)
c --- 'bdmtyp' = type of background diapycnal mixing
c --- 'bdmc1'  = background diapycnal diffusivity times buoyancy
c ---            frequency [cm**2/s**2]
c --- 'bdmc2'  = background diapycnal diffusivity [cm**2/s]
c --- 'tkepf'  = fraction of surface TKE that penetrates beneath mixed
c ---            layer
c --- 'niwgf'  = global factor applied to the energy input by
c ---            near-inertial motions
c --- 'niwbf'  = fraction of near-inertial energy dissipated in the
c ---            boundary layer
c --- 'niwlf'  = fraction of near-inertial energy dissipated locally
c ---            beneath the boundary layer
c --- 'csdiag' = if set to .true., then output check sums
c --- 'cnsvdi' = if set to .true., then output conservation diagnostics
c --- 'expcnf' = experiment configuration
c --- 'mommth' = momentum equation discretization method
c --- 'eitmth' = eddy-induced transport parameterization method
c --- 'edritp' = type of Richardson number used in eddy diffusivity
c ---            computation
c --- 'bmcmth' = baroclinic mass flux correction method
c --- 'rmpmth' = method of applying eddy-induced transport in the remap
c ---            transport algorithm
c --- 'edwmth' = method to estimate eddy diffusivity weight as a
c ---            function of the ration of Rossby radius of deformation
c ---            to the horizontal grid spacing
c --- 'mlrttp' = type of mixed layer restratification time scale
c --- 'edsprs' = if set to .true,, apply eddy mixing suppression away
c ---            from steering level
c --- 'grfile' = name of file containing grid specification
c --- 'icfile' = name of file containing initial conditions
c --- 'tdfile' = name of file containing tidal wave energy dissipation
c ---            divided by by bottom buoyancy frequency
c
      real baclin,batrop,mdv2hi,mdv2lo,mdv4hi,mdv4lo,mdc2hi,mdc2lo,
     .     vsc2hi,vsc2lo,vsc4hi,vsc4lo,slip,cbar,cb,cwbdts,cwbdls,
     .     wuv1,wuv2,wts1,wts2,wbaro,wpgf,mltmin,thktop,thkbot,egc,
     .     eggam,eglsmn,egmndf,egmxdf,egidfq,ri0,rm0,rm5,ce,
     .     bdmc1,bdmc2,tkepf,niwgf,niwbf,niwlf
      integer bdmtyp
      logical csdiag,cnsvdi,edsprs
      character*80 expcnf,mommth,eitmth,edritp,bmcmth,rmpmth,edwmth,
     .             mlrttp
      character*256 grfile,icfile,tdfile
c
      common /parms1/ baclin,batrop,mdv2hi,mdv2lo,mdv4hi,mdv4lo,
     .                mdc2hi,mdc2lo,vsc2hi,vsc2lo,vsc4hi,vsc4lo,slip,
     .                cbar,cb,cwbdts,cwbdls,wuv1,wuv2,wts1,wts2,wbaro,
     .                wpgf,mltmin,thktop,thkbot,egc,eggam,eglsmn,egmndf,
     .                egmxdf,egidfq,ri0,rm0,rm5,ce,bdmc1,bdmc2,tkepf,
     .                niwgf,niwbf,niwlf,bdmtyp,csdiag,cnsvdi,edsprs,
     .                expcnf,mommth,eitmth,edritp,bmcmth,rmpmth,edwmth,
     .                mlrttp,grfile,icfile,tdfile
!!!!!!!!!!!!!!!$acc declare device_resident(/parms1/)
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'g'      = gravity acceleration
c --- 'rearth' = radius of the earth
c --- 'spcifh' = specific heat of sea water (j/g/deg)
c --- 't0deg'  = zero degree celsius in kelvin (K)
c --- 'alpha0' = reference value of specific volume (cm**3/g)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c
      real tenm,onem,tencm,onecm,onemm,g,rearth,spcifh,t0deg,alpha0,
     .     epsil,radian,pi
c
      common /consts/ tenm,onem,tencm,onecm,onemm,g,rearth,spcifh,t0deg,
     .                alpha0,epsil,radian,pi
!!!!!!!!!!!!!!!$acc declare device_resident(/consts/)
c
c --- grid point where detailed diagnostics are desired:
c
      integer itest,jtest,ptest
c
      common /testpt/ itest,jtest,ptest
!!!!!!!!!!!!!!$acc declare device_resident(/testpt/)
c
      character*256 runid
      integer runid_len
c
      common /iovars/ runid,runid_len
!!!!!!!!!!$acc declare device_resident(/iovars/)
