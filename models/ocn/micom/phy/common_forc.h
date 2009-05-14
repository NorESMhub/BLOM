c --- ------------------------------------------------------------------
c --- common blocks related to the application forcing fields that is
c --- shared among the various sources of forcing
c --- ------------------------------------------------------------------
c
c --- fields for diagnosed relaxation fluxes
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,48) ::
     .  tflxap,        ! heat flux to be applied
     .  sflxap,        ! fsalt flux to be applied
     .  tflxdi,        ! diagnosed heat flux
     .  sflxdi         ! diagnosed salt flux
c
c --- monthly climatological fields used in the computation of
c --- climatological fluxes and relaxation of sea surface temperature
c --- and salinity
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) ::
     .  sstclm,        ! surface temperature
     .  ricclm,        ! ice concentration
     .  sssclm         ! sea surface salinity
c
c --- accumulation counter for diagnosed relaxation fluxes
      integer, dimension(48) :: nflxdi
c
c --- e-folding relaxation time scales
      real trxday,srxday
c
c --- interpolation parameters for monthly climatological fields
      real x
      integer l1,l2,l3,l4,l5
c
c --- flags concerning diagnosed heat and salt fluxes
      logical aptflx,apsflx,ditflx,disflx
c
      common /frc1/ tflxap,sflxap,tflxdi,sflxdi,sstclm,ricclm,sssclm,
     .              nflxdi,trxday,srxday,x,l1,l2,l3,l4,l5,aptflx,
     .              apsflx,ditflx,disflx
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
c
c --- fields to be specified in getflux or import_mct
     .  swa,           ! solar heat flux
     .  nsf,           ! non-solar heat flux
     .  dfl,           ! derivative of non-solar heat flux by surface temp.
     .  lip,           ! liquid water flux
     .  sop,           ! solid precipitation
     .  eva,           ! evaporation
     .  ztx,           ! u component of wind stress
     .  mty,           ! v component of wind stress
     .  rnf,           ! runoff
     .  ustarw,        ! friction velocity for open water
     .  tsi,           ! averaged snow/ice surface temperature
     .  slp,           ! sea level pressure
     .  abswnd,        ! wind speed at measurement height -zu-
     .  albw,          ! daily mean open water albedo
     .  frzpot,        ! freezing potential
     .  mltpot,        ! melting potential
     .  sfl,           ! salt flux
c
c --- fields to be specified in thermf
     .  tsi_tda,       ! accumulated snow/ice surface temperature
     .  tml_tda,       ! accumulated mixed layer temperature
     .  sml_tda,       ! accumulated mixed layer salinity
     .  alb_tda,       ! accumulated albedo
     .  fice_tda,      ! accumulated sea ice concentration
     .  ssu_tda,       ! accumulated sea surface velocity (u-comp.)
     .  ssv_tda,       ! accumulated sea surface velocity (v-comp.)
c
c --- albedo
     .  alb,         
c
     .  rnfres,        ! runoff reservoar
     .  rnfflx         ! runoff freshwater flux taken out of the reservoar
c
c --- accumulation number
      integer ntda
c
c --- time steps between surface forcing updates
      integer nfrco
c
      common /frc2/ swa,nsf,dfl,lip,sop,eva,ztx,mty,rnf,ustarw,tsi,slp,
     .              abswnd,albw,frzpot,mltpot,sfl,tsi_tda,tml_tda,
     .              sml_tda,alb_tda,fice_tda,ssu_tda,ssv_tda,alb,rnfres,
     .              rnfflx,ntda,nfrco
c
c --- constants set in 'frcdat'
      real albw_d,rhowat,t0deg
      integer nrfets
c
      common /frcpar/ albw_d,rhowat,t0deg,nrfets
