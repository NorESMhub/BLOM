!-----------------------------------------------------------------------
! Grid size related parameters
!   dx:         grid size in x-direction (before use of metric!)[m]
!   dy:         grid size in y-direction [m]
!   dxsq:       dx*dx [m**2]
!   dysq:       dy*dy [m**2]
!   sx2:        0.5/dxsq [1/m**2]
!   sy2:        0.5/dysq [1/m**2]
!   sxy:        0.25/(dx*dy) [1/m**2]
!-----------------------------------------------------------------------
      common /drv/ dx,dy,dxsq,dysq,sx2,sy2,sxy
      real dx,dy,dxsq,dysq,sx2,sy2,sxy
!-----------------------------------------------------------------------
!  dt:          length of time step [s]
!-----------------------------------------------------------------------
      common /step/ dt
      real dt
!-----------------------------------------------------------------------
!  pm,pn:       metric coefficients for x- and y-direction
!-----------------------------------------------------------------------
      common /coord/ pm(0:l,0:m),pn(0:l,0:m)
      real pm,pn
!-----------------------------------------------------------------------
!  f      Coriolis parameter [1/s]
!-----------------------------------------------------------------------
      common /coriolis/ f(0:L,0:M)
      real f
!-----------------------------------------------------------------------
!  h:           mean ice thickness
!  A:           ice concentration
!-----------------------------------------------------------------------
      common /thck/ h(0:L,0:M),A(0:L,0:M)
      real h,A
!-----------------------------------------------------------------------
! Ice drift velocity in x and y direction (in model coordinate system!)
!-----------------------------------------------------------------------
      common /vel/ u(l,0:m,3),v(l,0:m,3)
      real u,v
!-----------------------------------------------------------------------
! Sea ice related parameters
!   rhoice:     density of sea ice [kg/m**3]
!   rhowat:     density of sea water [kg/m**3]
!   cdwat:      oceanic drag coefficient []
!   angwin:     abl turning angle [deg]
!   angwat:     obl turning angle [deg]
!   sinwin:     sin of abl turning angle []
!   coswin:     cos of abl turning angle []
!   sinwat:     sin of obl turning angle []
!   coswat:     cos of obl turning angle []
!   g:          gravity constant [m/s**2]
!-----------------------------------------------------------------------
      common /parm/ sinwin(L,0:M),sinwat(L,0:M),
     &              rhoice,rhowat,cdwat,angwin,angwat,
     &              coswin,coswat,g
      real sinwin,sinwat,rhoice,rhowat,cdwat,angwin,angwat,
     &     coswin,coswat,g
!-----------------------------------------------------------------------
! Forcing variables
!   tauxw:      x-component of wind stress [N/m**2]
!   tauyw:      y-component of wind stress [N/m**2]
!   uwat:       x-component of current velocity [m/s]
!   vwat:       y-component of current velocity [m/s]
!   gradhx:     x-component of ocean surface tilt
!   gradhy:     y-component of ocean surface tilt
!-----------------------------------------------------------------------
      common /frc/ tauxw(L,0:M),tauyw(L,0:M),
     &             uwat2l(L,0:M,2),vwat2l(L,0:M,2),
     &             uwat(L,0:M),vwat(L,0:M),
     &             gradhx2l(L,0:M,2),gradhy2l(L,0:M,2),
     &             gradhx(L,0:M),gradhy(L,0:M)
      real tauxw,tauyw,uwat2l,vwat2l,uwat,vwat,gradhx2l,gradhy2l,
     &     gradhx,gradhy
!-----------------------------------------------------------------------
! Viscosity parameters
!   Pstar:      empirical ice strength parameter [N/m**2]
!   Cstar:      empirical ice strength constant []
!   eccen:      ratio of compressive to shear strength []
!   gmin:       maximum viscous creep rate [1/s]
!   ecm2:       1/eccen**2
!-----------------------------------------------------------------------
      common /viscp/ Pstar(0:L,0:M),Cstar,eccen,gmin,ecm2
      real Pstar,Cstar,eccen,gmin,ecm2
!-----------------------------------------------------------------------
! Viscosities and delta term
!   zeta:       bulk  viscosity
!   eta:        shear viscosity
!-----------------------------------------------------------------------
      common /viscos/ zeta(0:L,0:M),eta(0:L,0:M)
      real zeta,eta
!-----------------------------------------------------------------------
! Ice strength [N/m]
!-----------------------------------------------------------------------
      common /press/ P(0:L,0:M)
      real P
!-----------------------------------------------------------------------
! Terms of momentum equation (subroutine  relax)
!-----------------------------------------------------------------------
      common /rel/ ru(l,0:m),rv(l,0:m),den(l,0:m),amas(l,0:m),
     &             bu(l,0:m),bv(l,0:m),fx(l,0:m),fy(l,0:m),asy(l,0:m)
      real ru,rv,den,amas,bu,bv,fx,fy,asy
!-----------------------------------------------------------------------
! Parameters for overrelaxation routine
!   vrmax:      cut off velocity difference between iteration steps[m/s]
!   wt:         relaxation factor []
!   vrwt:       threshold to stop overrelaxation
!   amasmin:    minimum ice mass for relax
!   mmax:       maximum iteration steps []
!   mwt:        iteration steps with overrelaxation []
!-----------------------------------------------------------------------
      common /relaxp/ vrmax,wt,vrwt,amasmin,mmax,mwt
      real vrmax,wt,vrwt,amasmin
      integer mmax,mwt
!-----------------------------------------------------------------------
! Information about the number of iterations done by RELAX.
! For computation diagnostic purposes only.
!   mrelax:     number of iterations
!   numax:      number of cut-offs of large velocities
!   numaxi:     number of cut-offs of large velocities, ice cells only
!   nosolution: flag:   0=solution found, 1=not found
!-----------------------------------------------------------------------
      common /relaxi/ mrelax,numax,numaxi,nosolution
      integer mrelax,numax,numaxi,nosolution
!-----------------------------------------------------------------------
! Constant fields used by  relax (Markus Harder, 1994)
! These fields hold metric coefficients and grid cell dimensions
!   pmpnv:      PM * PN at vector grid points (dimensionless cell area)
!-----------------------------------------------------------------------
      common /relaxc/ pmpnv(L,0:M),
     &                sx2p(L,0:M),sx2m(L,0:M),sy2p(L,0:M),sy2m(L,0:M)
      real pmpnv,sx2p,sx2m,sy2p,sy2m
!-----------------------------------------------------------------------
! Masks
!   vm:         vector mask
!   hm:         scalar mask
!   om:         scalar mask without outflow points
!-----------------------------------------------------------------------
      common /mask/ vm(L,0:M),hm(0:L,0:M),om(0:L,0:M)
      integer vm,hm,om
!-----------------------------------------------------------------------
! Run control
!-----------------------------------------------------------------------
      common /run/ lold,lnew
      integer lold,lnew
