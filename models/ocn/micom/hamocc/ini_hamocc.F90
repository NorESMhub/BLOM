      SUBROUTINE INI_HAMOCC(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke&
     &            ,kpbe,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu,ptiestw    &
     &            ,kplyear,kplmonth,kplday,kpldtoce,pmonts             &
     &            ,pgila,pgiph,omask,dummy_tr,ndtr                     &
     &            ,rid,rid_len,path,path_len,path2,path2_len)      
!      SUBROUTINE INI_HAMOCC(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke&
!     &           ,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu,ptiestw          &
!     &           ,kplyear,kplmonth,kplday,kpldtoce,pyears,pmonts       &
!     &           ,pgila,pgiph)

!****************************************************************
!
!**** *INI_BGC* - initialize marine bio-geo-chemistry module.
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - initialize sediment layering
!     - initialize bgc variables
!     - calculation of chemical constants
!     - read restart fields
!     - calculate budgets
!
!     Method
!     -------
!     - Noch zu tun : Biharmonic tracer diffusion (mit AULAPTS > ALMZER)
!                     Convective adjustment
!                     sea ice growth influence on tracer concentration
!
!**   Interface.
!     ----------
!
!     *CALL*       *INI_BGC(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke
!                          ,pddpo,ptho,psao,pdlxp,pdlyp,ptiestu
!                          ,kplyear,kplmonth,kplday,kpldtoce)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpaufr*   - 1/0 for read / do not read restart file
!     *INTEGER* *kpicycli* - flag for cyclicity.
!     *REAL*    *pdt*      - ocean model time step [sec].
!     *INTEGER* *kpndtrun* - total no. of time steps of run.
!     *INTEGER* *kpie*     - zonal dimension of model grid.
!     *INTEGER* *kpje*     - meridional dimension of model grid.
!     *INTEGER* *kpke*     - vertical dimension of model grid.
!     *REAL*    *pddpo*    - size of scalar grid cell (3rd REAL) [m].
!     *REAL*    *ptho*     - potential temperature [deg C].
!     *REAL*    *psao*     - salinity [psu.].
!     *REAL*    *pdlxp*    - size of scalar grid cell (zonal) [m].
!     *REAL*    *pdlyp*    - size of scalar grid cell (meridional) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu*  - depth of level [m].
!     *INTEGER* *kplyear*  - year  in ocean restart date
!     *INTEGER* *kplmonth* - month in ocean restart date
!     *INTEGER* *kplday*   - day   in ocean restart date
!     *INTEGER* *kpldtoce* - step  in ocean restart date
!
!     Externals
!     ---------
!     READ_NAMELIST CALL INI_TIMESER_BGC BODENSED BELEG_BGC AUFR_BGC 
!     CHEMIN CHCK_BGC PRINT_BGC
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc, only: mnproc,lp,nfu
#ifdef PDYNAMIC_BGC
      use mo_dynamic
#endif /* PDYNAMIC_BGC */ 
 
      implicit none
      INTEGER :: kpie,kpje,kpke,kpbe,pmonts,ndtr
!      INTEGER :: pyears,pmonts,kpie,kpje,kpke
      INTEGER :: kplyear,kplmonth,kplday,kpldtoce
      INTEGER :: kpaufr,kpicycli,kpndtrun,k,l
      INTEGER :: ii,jj,kk,i,j
      
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: psao (kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
!      REAL :: pgila(kpie,kpje)
!      REAL :: pgiph(kpie,kpje)
      REAL :: pgila(kpie*2,kpje*2)
      REAL :: pgiph(kpie*2,kpje*2)
      REAL :: ptiestu(kpie,kpje,kpke+1),ptiestw(kpie,kpje,kpke+1)
      REAL :: omask(kpie,kpje)
      REAL :: dummy_tr(1-kpbe:kpie+kpbe,1-kpbe:kpje+kpbe,kpke,ndtr)
!      REAL :: zo(kpie,kpje),sicsno(kpie,kpje),sictho(kpie,kpje)
      REAL :: pdt
      character*(*) rid,path,path2
      integer rid_len,path_len,path2_len

! Define io units

      io_stdo_bgc = lp      !  standard out.
      io_stdi_bgc = 5       !  standard in.
      io_rsti_bgc = nfu     !  restart in. 
      io_rsto_bgc = nfu     !  restart out.

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 'HAMOCC initialisation'
      write(io_stdo_bgc,*) 'restart',kpaufr,kpicycli,kpndtrun
      write(io_stdo_bgc,*) 'dims',kpie,kpje,kpke,pmonts
      write(io_stdo_bgc,*) 'time',kplyear,kplmonth,kplday,kpldtoce
      write(io_stdo_bgc,*) 'time step',pdt
      write(io_stdo_bgc,*) 'saln,temp',psao(40,40,1),ptho(40,40,1)
      write(io_stdo_bgc,*) 'dp,dx,dy',pddpo(40,40,1),pdlxp(40,40),pdlyp(40,40)
      write(io_stdo_bgc,*) 'pu,pw',ptiestu(40,40,1),ptiestw(40,40,2)
      write(io_stdo_bgc,*) 'lat,lon',pgiph(40,40),pgila(40,40),omask(40,40)
      write(io_stdo_bgc,*) 'dummy_tr in',dummy_tr(40,40,1,3)
      endif
!                    
! Set control constants ( mo_control_bgc )
!
      dtbgc = pdt                   !  time step length [sec].
      ndtdaybgc=NINT(86400./dtbgc)  !  time steps per day [No].
      dtb=1./ndtdaybgc              !  time step length [days].
      
      icyclibgc = kpicycli
      ndtrunbgc = kpndtrun

!
! Initialize some namelist parameters
!
      isac = 1
      mean_2D_day = 0
      mean_3D_month = 0

!
! Initialize time step counter of run.
!
      ldtrunbgc = 0

      CALL ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

!                        
! Allocate memory : biology
!
      CALL ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

!                        
! Allocate memory : sediment
!
      CALL ALLOC_MEM_SEDMNT(kpie,kpje)

!                        
! Allocate memory : inorganic carbon cycle
!
      CALL ALLOC_MEM_CARBCH(kpie,kpje,kpke)

#ifdef PDYNAMIC_BGC
      CALL ALLOC_MEM_DYNAMIC(kpie,kpje,kpke)
#endif /* PDYNAMIC_BGC */

!                        
! Initialize sediment layering
!
      CALL BODENSED(kpie,kpje,kpke,pddpo)

!                        
! Initialize sediment and ocean tracer.
! 
!     ocetra(:,:,:,:)=dummy_tr(1:kpie,1:kpje,:,:)
      ocetra(:,:,:,:)=0.
 
      CALL BELEG_BGC(kpie,kpje,kpke,psao,ptho,pddpo,ptiestu,    &
     &               omask,pgila,pgiph,path,path_len)

!                        
! Initialize chemical constants: this routine needs 3D temp./sal. 
! input from the ocean model. For C-HOPE this means, that it must not be 
! called before SBR AUFR.
! 
! If kpaufr.eq.1 the initial chem. constants are read from the restart file.
! If kpaufr.eq.0 they have to be calculated now.
! Parameter "-13" is initializing all month.

      CALL CHEMCON(-13,kpie,kpje,kpke,psao,ptho,                &
     &     pddpo,pdlxp,pdlyp,ptiestu,kplmonth,omask)
     
!                        
! Read restart fields from restart file
!
      IF(kpaufr.eq.1) THEN
         CALL AUFR_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmonth,   &
     &                 kplday,kpldtoce,omask,rid,rid_len,       &
     &                 path2,path2_len)
      ENDIF

! aufsetz! (for initialization of 14C)
!      call c14_correction(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,psao,  &
!     &                    ptho,omask)

!#ifdef DIFFAT
! correction of alkalinity during spin-up of pco2 
!
!      CALL SPINUP_BGC(kpie,kpje,kpke,omask,pdlxp,pdlyp)
!#endif

!
! Global inventory of all tracers
!
      
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,1)

!      CALL CHCK_BGC(io_stdo_bgc,kpicycli,                         &
!     &'Check values of ocean tracer at exit from SBR INI_BGC :',  &
!     & kpie,kpje,kpke,pddpo) 

!
! Fill dummy_tr with HAMOCC initialised ocetra to pass to MICOM       
!
      dummy_tr(1:kpie,1:kpje,:,:)=ocetra(:,:,:,:)

      DO j=1,kpje
      DO i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
      DO k=1,kpke

           bgcm3d(i,j,k,jphyto)  = 			&
     &      bgcm3d(i,j,k,jphyto)  + ocetra(i,j,k,iphy)
            bgcm3d(i,j,k,jgrazer) = 			&
     &      bgcm3d(i,j,k,jgrazer) + ocetra(i,j,k,izoo)
            bgcm3d(i,j,k,jphosph) = 			&
     &      bgcm3d(i,j,k,jphosph) + ocetra(i,j,k,iphosph)
            bgcm3d(i,j,k,joxygen) = 			&
     &      bgcm3d(i,j,k,joxygen) + ocetra(i,j,k,ioxygen)
            bgcm3d(i,j,k,jiron)   = 			&
     &      bgcm3d(i,j,k,jiron)   + ocetra(i,j,k,iiron) 
            bgcm3d(i,j,k,jano3)   = 			&
     &      bgcm3d(i,j,k,jano3)   + ocetra(i,j,k,iano3)
            bgcm3d(i,j,k,jalkali) = 			&
     &      bgcm3d(i,j,k,jalkali) + ocetra(i,j,k,ialkali) 
            bgcm3d(i,j,k,jsilica) = 			&
     &      bgcm3d(i,j,k,jsilica) + ocetra(i,j,k,isilica)
            bgcm3d(i,j,k,jdic)    = 			&
     &      bgcm3d(i,j,k,jdic)    + ocetra(i,j,k,isco212)
            bgcm3d(i,j,k,jdoc)    = 			&
     &      bgcm3d(i,j,k,jdoc)    + ocetra(i,j,k,idoc)  
!            bgcm3d(i,j,k,jdms)    = 			&
!     &      bgcm3d(i,j,k,jdms)    + ocetra(i,j,k,idms) 
            bgcm3d(i,j,k,jpoc)    = 			&
     &      bgcm3d(i,j,k,jpoc)    + ocetra(i,j,k,idet)
            bgcm3d(i,j,k,jcalc)    = 			&
     &      bgcm3d(i,j,k,jcalc)    + ocetra(i,j,k,icalc) 
            bgcm3d(i,j,k,jopal)    = 			&
     &      bgcm3d(i,j,k,jopal)    + ocetra(i,j,k,iopal) 
#ifdef __c_isotopes
            bgcm3d(i,j,k,jdic13)  =                     &
     &      bgcm3d(i,j,k,jdic13)  + ocetra(i,j,k,isco213)
            bgcm3d(i,j,k,jdic14)  =                     &
     &      bgcm3d(i,j,k,jdic14)  + ocetra(i,j,k,isco214)
#endif
#ifdef AGG
            bgcm3d(i,j,k,jnos)    = 			&
	    bgcm3d(i,j,k,jnos)    + ocetra(i,j,k,inos)  
#endif   

      ENDDO
         ENDIF
      ENDDO
      ENDDO
!
!      call ncwrt3_bgc(bgcm3d)
!
      DO l=1,nbgcm3d
      DO j=1,kpje
      DO i=1,kpie
      DO k=1,kpke
            bgcm3d(i,j,k,l)  = 0.
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
