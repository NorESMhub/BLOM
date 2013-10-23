      SUBROUTINE AUFR_BGC(kpie,kpje,kpke,                       &
     &                    kplyear,kplmon,kplday,kpldtoce,       &
     &                    rstfnm_ocn,path)

!****************************************************************
!
!**** *AUFR_BGC* - reads marine bgc restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*     10.04.01
!     - extra SBR for reading bgc data from the restart file.
!     S.Legutke,        *MPI-MaD, HH*     15.08.01
!     - netCDF version (with cond.comp. PNETCDF)
!     - no use of chemc values from netCDF restart
!
!     Patrick Wetzel,   *MPI-Met, HH*     16.04.02
!     - read chemcm(i,j,7,12) from netCDF restart
!
!     J.Schwinger,      *GFI, Bergen*     2013-10-21
!     - removed reading of chemcm and ak* fields
!     - code cleanup, remoded preprocessor option "PNETCDF"
!       and "NOMPI"
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
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*       - 1st dimension of model grid.
!     *INTEGER* *kpje*       - 2nd dimension of model grid.
!     *INTEGER* *kpke*       - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmon*     - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *INTEGER* *kpldtoce*   - step  in ocean restart date
!     *CHAR*    *rstfnm_ocn* - restart file name-informations
!     *CHAR*    *path*       - path to restart files
!
!
!**************************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc, only: mnproc,xchalt

      implicit none
      
      INTEGER          :: kpie,kpje,kpke
      INTEGER          :: kplyear,kplmon,kplday,kpldtoce
      character(len=*) :: rstfnm_ocn,path

      ! Local variables
      INTEGER   :: idate(5)
      INTEGER   :: i,j,k,l
      INTEGER   :: restyear            !  year of restart file
      INTEGER   :: restmonth           !  month of restart file
      INTEGER   :: restday             !  day of restart file
      INTEGER   :: restdtoce           !  time step number from bgc ocean file
      character :: rstfnm*80

      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat,ncvarid

!
! Open netCDF data file
!
      IF(mnproc==1) THEN

#ifdef CCSMCOUPLED
        i=1
        do while (rstfnm_ocn(i:i+8).ne.'.micom.r.')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (io_stdo_bgc,*)                                    &
     &        'Could not generate restart file name!'
            call xchalt('(aufr_bgc)')
            stop '(aufr_bgc)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'.micom.rbgc.'//rstfnm_ocn(i+9:)
#else
        i=1
        do while (rstfnm_ocn(i:i+8).ne.'_restphy_')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (io_stdo_bgc,*)                                    &
     &        'Could not generate restart file name!'
            call xchalt('(aufr_bgc)')
            stop '(aufr_bgc)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'_rest_b_'//rstfnm_ocn(i+9:)
#endif

        ncstat = NF_OPEN(trim(path)//rstfnm,NF_NOWRITE, ncid)
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


!
! Read restart data : ocean aquateous tracer
!                
      CALL read_netcdf_var(ncid,'sco212',ocetra(1,1,1,isco212),kpke,0)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sco213',ocetra(1,1,1,isco213),kpke,0)
      CALL read_netcdf_var(ncid,'sco214',ocetra(1,1,1,isco214),kpke,0)
#endif
      CALL read_netcdf_var(ncid,'alkali',ocetra(1,1,1,ialkali),kpke,0)
      CALL read_netcdf_var(ncid,'phosph',ocetra(1,1,1,iphosph),kpke,0)
      CALL read_netcdf_var(ncid,'oxygen',ocetra(1,1,1,ioxygen),kpke,0)
      CALL read_netcdf_var(ncid,'gasnit',ocetra(1,1,1,igasnit),kpke,0)
      CALL read_netcdf_var(ncid,'ano3',ocetra(1,1,1,iano3),kpke,0)
      CALL read_netcdf_var(ncid,'silica',ocetra(1,1,1,isilica),kpke,0)
      CALL read_netcdf_var(ncid,'doc',ocetra(1,1,1,idoc),kpke,0)
      CALL read_netcdf_var(ncid,'phyto',ocetra(1,1,1,iphy),kpke,0)
      CALL read_netcdf_var(ncid,'grazer',ocetra(1,1,1,izoo),kpke,0)
      CALL read_netcdf_var(ncid,'poc',ocetra(1,1,1,idet),kpke,0)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'poc13',ocetra(1,1,1,idet13),kpke,0)
      CALL read_netcdf_var(ncid,'poc14',ocetra(1,1,1,idet14),kpke,0)
#endif
      CALL read_netcdf_var(ncid,'calciu',ocetra(1,1,1,icalc),kpke,0)
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
#ifdef CFC
      CALL read_netcdf_var(ncid,'cfc11',ocetra(1,1,1,icfc11),kpke,0)
      CALL read_netcdf_var(ncid,'cfc12',ocetra(1,1,1,icfc12),kpke,0)
      CALL read_netcdf_var(ncid,'sf6',ocetra(1,1,1,isf6),kpke,0)
#endif

!
! Read restart data : other fields
!
     CALL read_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'satn2o',satn2o(1,1),1,0)

!
! Read restart data : sediment variables.
      CALL read_netcdf_var(ncid,'ssso12',sedlay(1,1,1,issso12),ks,0)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'ssso13',sedlay(1,1,1,issso13),ks,0)
      CALL read_netcdf_var(ncid,'ssso14',sedlay(1,1,1,issso14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'sssc12',sedlay(1,1,1,isssc12),ks,0)
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


      RETURN
      END
