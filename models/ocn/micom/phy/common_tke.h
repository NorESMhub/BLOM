c --- ------------------------------------------------------------------
c --- common blocks related to second order turbulence closure
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) :: 
     .  tke,           ! turbulent kinetic energy
     .  gls_psi        ! generic length scale 
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: 
     .  tke_old,       ! turbulent kinetic energy at old time level
     .  gls_psi_old,   ! generic length scale at old time level
     .  Prod,          ! shear production
     .  Buoy,          ! buoyancy production
     .  Shear2,        ! square of the shear frequency
     .  L_scale        ! dissipative length scale
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  uflxtke,vflxtke,         ! horizontal fluxes of tke
     .  uflxgls_psi,vflxgls_psi, ! horizontal fluxes of generic length scale
     .  tkeflx,        ! surface flux of tke
     .  gls_psiflx     ! surface flux of generic length scale
c
      common /tke1/ tke,gls_psi,tke_old,gls_psi_old,
     .              Prod,Buoy,Shear2,L_scale,
     .              uflxtke,vflxtke,uflxgls_psi,vflxgls_psi,
     .              tkeflx,gls_psiflx
c
c --- various coefficients
      real gls_s0,gls_s1,gls_s2,gls_s3,gls_s4,gls_s5,gls_s6,
     .     gls_b0,gls_b1,gls_b2,gls_b3,gls_b4,gls_b5,
     .     sqrt2,cmu_fac1,cmu_fac2,cmu_fac3,tke_exp1,gls_exp1,gls_fac6
c
      common /tke2/ gls_s0,gls_s1,gls_s2,gls_s3,gls_s4,gls_s5,gls_s6, 
     .              gls_b0,gls_b1,gls_b2,gls_b3,gls_b4,gls_b5,
     .              sqrt2,cmu_fac1,cmu_fac2,cmu_fac3,tke_exp1,
     .              gls_exp1,gls_fac6
