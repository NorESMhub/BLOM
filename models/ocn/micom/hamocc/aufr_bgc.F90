      SUBROUTINE AUFR_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmon     &
     &                    ,kplday,kpldtoce,omask,rid,rid_len      &
     &                    ,rstfnm_ext,path,path_len)

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/aufr_bgc.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

!****************************************************************
!
!**** *AUFR_BGC* - reads marine bgc restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - extra SBR for reading bgc data from the restart file.
!     S.Legutke,        *MPI-MaD, HH*    15.08.01
!     - netCDF version (with cond.comp. PNETCDF)
!     - no use of chemc values from netCDF restart
!
!     Patrick Wetzel,    *MPI-Met, HH*    16.04.02
!     - read chemcm(i,j,7,12) from netCDF restart
!
!     Purpose
!     -------
!     Read restart data to continue an interrupted integration.
!
!     Method
!     -------
!     The bgc data are read from an extra file, other than the ocean data.
!     The time stamp of the bgc restart file (idate) is specified from the
!     ocean time stamp through the SBR parameter list of AUFW_BGC. The only 
!     time control variable proper to the bgc is the time step number 
!     (idate(5)). It can differ from that of the ocean (idate(4)) by the 
!     difference of the offsets of restart files.
!
!**   Interface.
!     ----------
!
!     *CALL*       *AUFR_BGC(kpie,kpje,kpke,pddpo
!                            ,kplyear,kplmon,kplday,kpldtoce)*
!
!     *COMMON*     *MO_PARAM1_BGC* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!     *INTEGER* *kplyear*   - year  in ocean restart date
!     *INTEGER* *kplmon*  - month in ocean restart date
!     *INTEGER* *kplday*    - day   in ocean restart date
!     *INTEGER* *kpldtoce*  - step  in ocean restart date
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik IK introduced ocetra array elements for phyto, grazer, poc (=det), calciu 
!ik array indices are: iphy, izoo, idet, icalc
!iktodo IK introduced new variable opal (index iopal)
!ik nocetra is the number of all BGC element (size of ocetra(,,,l))
!ik nosedi is the number of all elements interacting with the sediment

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc

      implicit none
      
      INTEGER  kpie,kpje,kpke
      INTEGER  kplyear,kplmon,kplday,kpldtoce
      REAL pddpo(kpie,kpje,kpke)
      REAL chemcm_t(kpie,kpje,8)
      REAL omask(kpie,kpje)
      INTEGER  i,j,k,l,kmon

      INTEGER idate(5)

      INTEGER :: restyear            !  year of restart file
      INTEGER :: restmonth           !  month of restart file
      INTEGER :: restday             !  day of restart file
      INTEGER :: restdtoce           !  time step number from bgc ocean file

      character rstfnm*80
      character*(*) rid,rstfnm_ext,path
      integer rid_len,path_len

#ifdef PNETCDF                
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat,ncvarid

!
! Open netCDF data file
!
      IF(mnproc==1) THEN

#ifdef CCSMCOUPLED
      rstfnm=rid(1:rid_len)//'.micom.rbgc'//rstfnm_ext
#else
      rstfnm=rid(1:rid_len)//'_rest_b'//rstfnm_ext
#endif

        ncstat = NF_OPEN(path(1:path_len)//rstfnm,NF_NOWRITE, ncid)
        IF ( ncstat .NE. NF_NOERR ) THEN
             CALL xchalt('(AUFR: Problem with netCDF1)')
                    stop '(AUFR: Problem with netCDF1)'
        ENDIF

!
! Read restart data : date
!

        ncstat = NF_GET_ATT_INT(ncid, NF_GLOBAL,'date', idate)
        IF ( ncstat .NE. NF_NOERR ) THEN
          CALL xchalt('(AUFR: Problem reading date of restart file)')
                 stop '(AUFR: Problem reading date of restart file)'
        ENDIF
        restyear  = idate(1)
        restmonth = idate(2)
        restday   = idate(3)
        restdtoce = idate(4)
        ldtbgc = idate(5)
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'Date of bgc restart file : '
        WRITE(io_stdo_bgc,*) ' year  = ',restyear
        WRITE(io_stdo_bgc,*) ' month = ',restmonth
        WRITE(io_stdo_bgc,*) ' day   = ',restday
        WRITE(io_stdo_bgc,*) ' dtoce = ',restdtoce
        WRITE(io_stdo_bgc,*) ' dtbgc = ',ldtbgc
        WRITE(io_stdo_bgc,*) ' '
      ENDIF

!
! Compare with date read from ocean restart file
!
! As the ocean is already in its first step, its counter has 
! gone up one step already for the year and month. The ocean day
! counter is still at its restart date. Therefore:

!      restmonth = restmonth + 1
!
!      IF (restmonth .GT. 12) THEN
!         restmonth=1
!         restyear=restyear+1
!      ENDIF

! Memorize ocean start date :   
      bgcstartyear  = kplyear
      bgcstartmonth = kplmon
      bgcstartday   = kplday

      IF ( kplyear  .NE. restyear  ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart years in oce/bgc are not the same : '  &
     &   ,kplyear,'/',restyear,' !!!'
         ENDIF
      ENDIF

      IF ( kplmon .NE. restmonth ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart months in oce/bgc are not the same : '   &
     &   ,kplmon,'/',restmonth,' !!!'
         ENDIF
!         STOP 'Stop : restart months in oce/bgc are not the same.'
      ENDIF

      IF ( kplday   .NE. restday   ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart days in oce/bgc are not the same : '   &
     &   ,kplday,'/',restday,' !!!'
         ENDIF
!         STOP 'Stop : restart days in oce/bgc are not the same.'
      ENDIF 

!      IF ( kpldtoce .NE. ldtbgc   ) THEN
!         WRITE(io_stdo_bgc,*)                                       & 
!     &   'WARNING: restart step no.  in oce/bgc are not the same : '&
!     &   ,kpldtoce,'/',ldtbgc,' !!!'
!      ENDIF

!
! Read restart data : ocean aquateous tracer
!                
      CALL read_netcdf_var(ncid,'sco212',ocetra(1,1,1,isco212),kpke,0)
!      call chksumbgc(ocetra(1,1,1,isco212),kpke,'sco212')
! comment out for first run (no 13c in restart file)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sco213',ocetra(1,1,1,isco213),kpke,0)
!      call chksumbgc(ocetra(1,1,1,isco213),kpke,'sco213')
      CALL read_netcdf_var(ncid,'sco214',ocetra(1,1,1,isco214),kpke,0)
!      call chksumbgc(ocetra(1,1,1,isco214),kpke,'sco214')
#endif
      CALL read_netcdf_var(ncid,'alkali',ocetra(1,1,1,ialkali),kpke,0)
!      call chksumbgc(ocetra(1,1,1,ialkali),kpke,'alkali')
      CALL read_netcdf_var(ncid,'phosph',ocetra(1,1,1,iphosph),kpke,0)
      CALL read_netcdf_var(ncid,'oxygen',ocetra(1,1,1,ioxygen),kpke,0)
      CALL read_netcdf_var(ncid,'gasnit',ocetra(1,1,1,igasnit),kpke,0)
      CALL read_netcdf_var(ncid,'ano3',ocetra(1,1,1,iano3),kpke,0)
      CALL read_netcdf_var(ncid,'silica',ocetra(1,1,1,isilica),kpke,0)
      CALL read_netcdf_var(ncid,'doc',ocetra(1,1,1,idoc),kpke,0)
      CALL read_netcdf_var(ncid,'phyto',ocetra(1,1,1,iphy),kpke,0)
      CALL read_netcdf_var(ncid,'grazer',ocetra(1,1,1,izoo),kpke,0)
      CALL read_netcdf_var(ncid,'poc',ocetra(1,1,1,idet),kpke,0)
! comment out for first run (no 13c in restart file)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'poc13',ocetra(1,1,1,idet13),kpke,0)
      CALL read_netcdf_var(ncid,'poc14',ocetra(1,1,1,idet14),kpke,0)
#endif
      CALL read_netcdf_var(ncid,'calciu',ocetra(1,1,1,icalc),kpke,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'calciu13',ocetra(1,1,1,icalc13),kpke,0)
      CALL read_netcdf_var(ncid,'calciu14',ocetra(1,1,1,icalc14),kpke,0)
#endif
      CALL read_netcdf_var(ncid,'opal',ocetra(1,1,1,iopal),kpke,0)
      CALL read_netcdf_var(ncid,'n2o',ocetra(1,1,1,ian2o),kpke,0)
      CALL read_netcdf_var(ncid,'dms',ocetra(1,1,1,idms),kpke,0)
      CALL read_netcdf_var(ncid,'fdust',ocetra(1,1,1,ifdust),kpke,0)
      CALL read_netcdf_var(ncid,'iron',ocetra(1,1,1,iiron),kpke,0)

#ifdef AGG
      CALL read_netcdf_var(ncid,'snos',ocetra(1,1,1,inos),kpke,0)
      CALL read_netcdf_var(ncid,'adust',ocetra(1,1,1,iadust),kpke,0)
#endif /*AGG*/

#ifdef ANTC14
      CALL read_netcdf_var(ncid,'antc14',ocetra(1,1,1,iantc14),kpke,0)
#endif
#ifdef PCFC
      CALL read_netcdf_var(ncid,'cfc11',ocetra(1,1,1,icfc11),kpke,0)
      CALL read_netcdf_var(ncid,'cfc12',ocetra(1,1,1,icfc12),kpke,0)
#endif

!
!Check aquateous restart data for topography
! 

!      DO i    =1,kpie
!      DO j    =1,kpje
!      DO k    =1,kpke      
!      DO l    =1,nocetra 
!	IF (omask(i,j) .le. 0.5 ) THEN
!	  IF ( ocetra(i,j,k,l) .NE. rmasko ) THEN
!             WRITE(io_stdo_bgc,*) 'ocetra not properly masked at :'  &
!     &                             ,i,j,k,l,ocetra(i,j,k,l)
!	  ENDIF
!	ELSE
!	  IF ( ocetra(i,j,k,l) .EQ. rmasko ) THEN
!             WRITE(io_stdo_bgc,*) 'land mask values at wet points :'  &
!     &                             ,i,j,k,l,ocetra(i,j,k,l)
!             call STOP_ALL('Stop : Restart file with different topography:     &
!     &	     land mask values at wet points')
!
!	  ENDIF	 
!	ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO

!
! Read restart data : other fields
!
     CALL read_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'akw3',akw3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'akb3',akb3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'ak13',ak13(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'ak23',ak23(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'aksp',aksp(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'satn2o',satn2o(1,1),1,0)

!
! Read restart data : chemical constants
!
      DO kmon =1,12
      CALL read_netcdf_var(ncid,'chemcm',chemcm_t,8,kmon)
      DO i    =1,kpie
      DO j    =1,kpje
	 IF (omask(i,j) .lt. 0.5 ) THEN
	    DO k=1,8
	       chemcm(i,j,k,kmon)  =   rmasko
	    ENDDO
	 ELSE
	   chemcm(i,j,1,kmon) = chemcm_t(i,j,1)
	   chemcm(i,j,2,kmon) = chemcm_t(i,j,2)
	   chemcm(i,j,3,kmon) = chemcm_t(i,j,3)
	   chemcm(i,j,4,kmon) = chemcm_t(i,j,4)
	   chemcm(i,j,5,kmon) = chemcm_t(i,j,5)
	   chemcm(i,j,6,kmon) = chemcm_t(i,j,6)
	   chemcm(i,j,7,kmon) = chemcm_t(i,j,7)
	   chemcm(i,j,8,kmon) = chemcm_t(i,j,8)
	 ENDIF
      ENDDO
      ENDDO
      ENDDO
!      call chksumbgc(chemcm,8*12,'chemcm')

!
! Read restart data : sediment variables.
! js: reading of terrigenous sediment was missing until 02.08.2005
      CALL read_netcdf_var(ncid,'ssso12',sedlay(1,1,1,issso12),ks,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'ssso13',sedlay(1,1,1,issso13),ks,0)
      CALL read_netcdf_var(ncid,'ssso14',sedlay(1,1,1,issso14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'sssc12',sedlay(1,1,1,isssc12),ks,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sssc13',sedlay(1,1,1,isssc13),ks,0)
      CALL read_netcdf_var(ncid,'sssc14',sedlay(1,1,1,isssc14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'ssssil',sedlay(1,1,1,issssil),ks,0)
      CALL read_netcdf_var(ncid,'ssster',sedlay(1,1,1,issster),ks,0)
      CALL read_netcdf_var(ncid,'bur_o12',burial(1,1,issso12),1,0)
      CALL read_netcdf_var(ncid,'bur_c12',burial(1,1,isssc12),1,0)
      CALL read_netcdf_var(ncid,'bur_sil',burial(1,1,issssil),1,0)
      CALL read_netcdf_var(ncid,'bur_clay',burial(1,1,issster),1,0)
      CALL read_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0)
      CALL read_netcdf_var(ncid,'powaic',powtra(1,1,1,ipowaic),ks,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'powc13',powtra(1,1,1,ipowc13),ks,0)
      CALL read_netcdf_var(ncid,'powc14',powtra(1,1,1,ipowc14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'powaal',powtra(1,1,1,ipowaal),ks,0)
      CALL read_netcdf_var(ncid,'powaph',powtra(1,1,1,ipowaph),ks,0)
      CALL read_netcdf_var(ncid,'powaox',powtra(1,1,1,ipowaox),ks,0)
      CALL read_netcdf_var(ncid,'pown2',powtra(1,1,1,ipown2),ks,0)
      CALL read_netcdf_var(ncid,'powno3',powtra(1,1,1,ipowno3),ks,0)
      CALL read_netcdf_var(ncid,'powasi',powtra(1,1,1,ipowasi),ks,0)

#ifdef DIFFAT 
!
! Read restart data : co2 diffusion
!
      CALL read_netcdf_var(ncid,'atmco2',atm(1,1,iatmco2),1,0)
      CALL read_netcdf_var(ncid,'atmo2',atm(1,1,iatmo2),1,0)
      CALL read_netcdf_var(ncid,'atmn2',atm(1,1,iatmn2),1,0)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'atmc13',atm(1,1,iatmc13),1,0)
      CALL read_netcdf_var(ncid,'atmc14',atm(1,1,iatmc14),1,0)
#endif

#endif

      if(mnproc==1) ncstat = NF_CLOSE(ncid)

#else

     ! Attention - this doesn't work any more for MPI
     ! One big FORTRAN READ like here is a very, very bad idea
     ! for parallel programs !!!!!

#ifndef NOMPI
     call stop_all("can't read in bgc restart with MPI and no NETCDF")
! would need to read these in on p_io and scatter if we were actually i
! going to do this - c.f. read_netcdf_var.F90

#else

      OPEN(io_rsti_bgc,FILE='restart_bgc',STATUS='UNKNOWN'             &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      READ(io_rsti_bgc)                                                &
     &           (((ocetra(i,j,k,iphosph),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,isilica),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,ioxygen),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,iphy   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,izoo   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,idet   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,idoc   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,icalc  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,iano3  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,igasnit),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,iopal  ),i=1,kpie),j=1,kpje),k=1,kpke)


      READ(io_rsti_bgc) chemcm,hi,co3,aksp                             &
     &          ,(((ocetra(i,j,k,isco212),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,ialkali),i=1,kpie),j=1,kpje),k=1,kpke)

      READ(io_rsti_bgc) sedlay,sedhpl

      READ(io_rsti_bgc) 
     &           (((powtra(i,j,k,ipowaic),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowaal),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowaph),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowaox),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowasi),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowno3),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipown2) ,i=1,kpie),j=1,kpje),k=1,ks)

      CLOSE (io_rsti_bgc)
#endif

#endif

!     
!  Masking aqueous sea water tracer.
!
!      DO l=1,nocetra
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(pddpo(i,j,k) .LT. 0.5) THEN
!         ocetra(i,j,k,l)=rmasko
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
! only for the first run!!      
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!        ocetra(i,j,k,isco214)=ocetra(i,j,k,isco212)*0.75
!      ENDDO
!      ENDDO
!      ENDDO      
!
!  Masking other sea water tracer.
!
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(pddpo(i,j,k) .LT. 0.5) THEN
!ik      phyto(i,j,k)=rmasko
!ik      grazer(i,j,k)=rmasko
!ik      poc(i,j,k)=rmasko
!         hi(i,j,k)=rmasko
!ik      calciu(i,j,k)=rmasko
!         co3(i,j,k)=rmasko
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO

!
!  Masking sediment pore water tracer.
!
!      DO  k=1,ks
!      DO  j=1,kpje
!      DO  i=1,kpie 
!      IF(kbo(i,j) .LE. 1) THEN
!         powtra(i,j,k,ipowaic)=rmasks
!         powtra(i,j,k,ipowaal)=rmasks
!         powtra(i,j,k,ipowaph)=rmasks
!         powtra(i,j,k,ipowaox)=rmasks
!         powtra(i,j,k,ipown2)=rmasks
!         powtra(i,j,k,ipowno3)=rmasks
!         powtra(i,j,k,ipowasi)=rmasks
!         sedlay(i,j,k,issso12)=rmasks
!         sedlay(i,j,k,isssc12)=rmasks
!         sedlay(i,j,k,issssil)=rmasks
!         burial(i,j,issso12)=rmasks
!         burial(i,j,isssc12)=rmasks
!         burial(i,j,issssil)=rmasks
!         sedhpl(i,j,k)=rmasks
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO

!     
!  Restrict to positive values (until error is found only !!!!)
!

!      DO l=1,nocetra
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(omask(i,j) .GT. 0.5) THEN
!	 ocetra(i,j,k,l)=MAX(ocetra(i,j,k,l),0.)
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO

!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(omask(i,j) .GT. 0.5) THEN
!ik         phyto (i,j,k)=MAX(phyto (i,j,k),0.)
!ik         grazer(i,j,k)=MAX(grazer(i,j,k),0.)
!ik         poc   (i,j,k)=MAX(poc   (i,j,k),0.)
!        hi    (i,j,k)=MAX(hi    (i,j,k),1.e-12)
!ik         calciu(i,j,k)=MAX(calciu(i,j,k),0.)
!ka multiply DIC to fix atmospheric CO2
!         ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)*0.995
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO


!      DO  k=1,ks
!      DO  j=1,kpje
!      DO  i=1,kpie 
!      IF(omask(i,j) .GT. 0.5) THEN
!         powtra(i,j,k,ipowaic)=MAX(powtra(i,j,k,ipowaic),0.)
!#ifdef __c_isotopes
!         powtra(i,j,k,ipowc13)=MAX(powtra(i,j,k,ipowc13),0.)
!         powtra(i,j,k,ipowc14)=MAX(powtra(i,j,k,ipowc14),0.)
!#endif
!         powtra(i,j,k,ipowaal)=MAX(powtra(i,j,k,ipowaal),0.)
!         powtra(i,j,k,ipowaph)=MAX(powtra(i,j,k,ipowaph),0.)
!         powtra(i,j,k,ipowaox)=MAX(powtra(i,j,k,ipowaox),0.)
!         powtra(i,j,k,ipown2) =MAX(powtra(i,j,k,ipown2) ,0.)
!         powtra(i,j,k,ipowno3)=MAX(powtra(i,j,k,ipowno3),0.)
!         powtra(i,j,k,ipowasi)=MAX(powtra(i,j,k,ipowasi),0.)
!         sedlay(i,j,k,issso12)=MAX(sedlay(i,j,k,issso12),0.)
!         sedlay(i,j,k,isssc12)=MAX(sedlay(i,j,k,isssc12),0.)
!#ifdef __c_isotopes
!         sedlay(i,j,k,isssc13)=MAX(sedlay(i,j,k,isssc13),0.)
!         sedlay(i,j,k,isssc14)=MAX(sedlay(i,j,k,isssc14),0.)
!#endif
!         sedlay(i,j,k,issssil)=MAX(sedlay(i,j,k,issssil),0.)
!         sedhpl(i,j,k)        =MAX(sedhpl(i,j,k)    ,1.e-12)
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO

      RETURN
      END
