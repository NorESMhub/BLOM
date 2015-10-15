      SUBROUTINE AUFR_BGC(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,     &
     &                    trc,sedlay2,powtra2,burial2,          &
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
!     J.Schwinger,      *GFI, Bergen*     2014-05-21
!     - adapted code for writing of two time level tracer 
!       and sediment fields
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
!     *INTEGER* *ntr*        - number of tracers in tracer field
!     *INTEGER* *ntrbgc*     - number of biogechemical tracers in tracer field
!     *INTEGER* *itrbgc*     - start index for biogeochemical tracers in tracer field
!     *REAL*    *trc*        - initial/restart tracer field to be passed to the 
!                              ocean model [mol/kg]
!     *REAL*    *sedlay2*    - initial/restart sediment (two time levels) field
!     *REAL*    *powtra2*    - initial/restart pore water tracer (two time levels) field
!     *REAL*    *sedhpl2*    - initial/restart pore water ph (two time levels) field
!     *REAL*    *burial2*    - initial/restart sediment burial (two time levels) field
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
      USE mo_control_bgc
      use mo_param1_bgc 
      USE mo_sedmnt, only:sedhpl
      use mod_xc
      use mod_dia, only : iotype
      use netcdf
      implicit none
#ifdef PNETCDF
#include <pnetcdf.inc>
#endif
#include <mpif.h>      
      INTEGER          :: kpie,kpje,kpke,ntr,ntrbgc,itrbgc
      REAL             :: trc(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,2*kpke,ntr)
      REAL             :: sedlay2(kpie,kpje,2*ks,nsedtra)
      REAL             :: powtra2(kpie,kpje,2*ks,npowtra)
      REAL             :: burial2(kpie,kpje,2,   nsedtra)
      INTEGER          :: kplyear,kplmon,kplday,kpldtoce
      character(len=*) :: rstfnm_ocn,path

      ! Local variables
      REAL      :: locetra(kpie,kpje,2*kpke,nocetra) ! local array for reading 
      INTEGER   :: restyear                          !  year of restart file
      INTEGER   :: restmonth                         !  month of restart file
      INTEGER   :: restday                           !  day of restart file
      INTEGER   :: restdtoce                         !  time step number from bgc ocean file
      INTEGER   :: idate(5),i,j
      character :: rstfnm*80

      INTEGER ncid,ncstat,ncvarid
      integer*4 ,save :: info=MPI_INFO_NULL
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),                      &
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
      character(len=3) :: stripestr
      character(len=9) :: stripestr2
      integer ierr,testio
      locetra(:,:,:,:) = 0.0
!
! Open netCDF data file
!
      testio=0
      IF(mnproc==1 .AND. IOTYPE==0) THEN

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

        ncstat = NF90_OPEN(trim(path)//rstfnm,NF90_NOWRITE, ncid)
        IF ( ncstat .NE. NF90_NOERR ) THEN
             CALL xchalt('(AUFR: Problem with netCDF1)')
                    stop '(AUFR: Problem with netCDF1)'
        ENDIF

!
! Read restart data : date
!

        ncstat = NF90_GET_ATT(ncid, NF90_GLOBAL,'date', idate)
        IF ( ncstat .NE. NF90_NOERR ) THEN
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
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
        testio=1
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
        write(stripestr,('(i3)')) 16
        write(stripestr2,('(i9)')) 1024*1024
        call mpi_info_create(info,ierr)
        call mpi_info_set(info,'romio_ds_read','disable',ierr)
        call mpi_info_set(info,'romio_ds_write','disable',ierr)
        call mpi_info_set(info,"striping_factor",stripestr,ierr)
        call mpi_info_set(info,"striping_unit",stripestr2,ierr)

        ncstat = NFMPI_OPEN(mpicomm,trim(path)//rstfnm,NF_NOWRITE,   &
     &    INFO, ncid)
        IF ( ncstat .NE. NF_NOERR ) THEN
             CALL xchalt('(AUFR: Problem with netCDF1)')
                    stop '(AUFR: Problem with netCDF1)'
        ENDIF

!
! Read restart data : date
!

        ncstat = NFMPI_GET_ATT_INT(ncid, NF_GLOBAL,'date', idate)
        IF ( ncstat .NE. NF_NOERR ) THEN
          CALL xchalt('(AUFR: Problem reading date of restart file)')
                 stop '(AUFR: Problem reading date of restart file)'
        ENDIF
        restyear  = idate(1)
        restmonth = idate(2)
        restday   = idate(3)
        restdtoce = idate(4)
        ldtbgc = idate(5)
        IF(mnproc==1) THEN
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'Date of bgc restart file : '
        WRITE(io_stdo_bgc,*) ' year  = ',restyear
        WRITE(io_stdo_bgc,*) ' month = ',restmonth
        WRITE(io_stdo_bgc,*) ' day   = ',restday
        WRITE(io_stdo_bgc,*) ' dtoce = ',restdtoce
        WRITE(io_stdo_bgc,*) ' dtbgc = ',ldtbgc
        WRITE(io_stdo_bgc,*) ' '
        ENDIF
#endif
      if(testio .eq. 0) then
      CALL xchalt('(AUFR: Problem with namelist iotype)')
                    stop '(AUFR: Problem with namelist iotype)'
      endif
      ENDIF

!
! Compare with date read from ocean restart file

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
      CALL read_netcdf_var(ncid,'sco212',locetra(1,1,1,isco212),2*kpke,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sco213',locetra(1,1,1,isco213),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'sco214',locetra(1,1,1,isco214),2*kpke,0,iotype)
#endif
      CALL read_netcdf_var(ncid,'alkali',locetra(1,1,1,ialkali),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'phosph',locetra(1,1,1,iphosph),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'oxygen',locetra(1,1,1,ioxygen),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'gasnit',locetra(1,1,1,igasnit),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'ano3',locetra(1,1,1,iano3),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'silica',locetra(1,1,1,isilica),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'doc',locetra(1,1,1,idoc),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'phyto',locetra(1,1,1,iphy),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'grazer',locetra(1,1,1,izoo),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'poc',locetra(1,1,1,idet),2*kpke,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'poc13',locetra(1,1,1,idet13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'poc14',locetra(1,1,1,idet14),2*kpke,0,iotype)
#endif
      CALL read_netcdf_var(ncid,'calciu',locetra(1,1,1,icalc),2*kpke,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'calciu13',locetra(1,1,1,icalc13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'calciu14',locetra(1,1,1,icalc14),2*kpke,0,iotype)
#endif
      CALL read_netcdf_var(ncid,'opal',locetra(1,1,1,iopal),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'n2o',locetra(1,1,1,ian2o),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'dms',locetra(1,1,1,idms),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'fdust',locetra(1,1,1,ifdust),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'iron',locetra(1,1,1,iiron),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefo2',locetra(1,1,1,iprefo2),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefpo4',locetra(1,1,1,iprefpo4),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefalk',locetra(1,1,1,iprefalk),2*kpke,0,iotype)

#ifdef AGG
      CALL read_netcdf_var(ncid,'snos',locetra(1,1,1,inos),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'adust',locetra(1,1,1,iadust),2*kpke,0,iotype)
#endif /*AGG*/

#ifdef ANTC14
      CALL read_netcdf_var(ncid,'antc14',locetra(1,1,1,iantc14),2*kpke,0,iotype)
#endif
#ifdef CFC
#ifdef RESTART_CFC
      CALL read_netcdf_var(ncid,'cfc11',locetra(1,1,1,icfc11),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'cfc12',locetra(1,1,1,icfc12),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'sf6',locetra(1,1,1,isf6),2*kpke,0,iotype)
#endif
#endif
!
! Read restart data : diagnostic ocean fields
!
     CALL read_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0,iotype)
     CALL read_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0,iotype)
     CALL read_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0,iotype)
     CALL read_netcdf_var(ncid,'satn2o',satn2o(1,1),1,0,iotype)

!
! Read restart data : sediment variables.
!
      CALL read_netcdf_var(ncid,'ssso12',sedlay2(1,1,1,issso12),2*ks,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'ssso13',sedlay2(1,1,1,issso13),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'ssso14',sedlay2(1,1,1,issso14),2*ks,0,iotype)
#endif
      CALL read_netcdf_var(ncid,'sssc12',sedlay2(1,1,1,isssc12),2*ks,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sssc13',sedlay2(1,1,1,isssc13),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'sssc14',sedlay2(1,1,1,isssc14),2*ks,0,iotype)
#endif
      CALL read_netcdf_var(ncid,'ssssil',sedlay2(1,1,1,issssil),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'ssster',sedlay2(1,1,1,issster),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'bur_o12',burial2(1,1,1,issso12),2,0,iotype)
      CALL read_netcdf_var(ncid,'bur_c12',burial2(1,1,1,isssc12),2,0,iotype)
      CALL read_netcdf_var(ncid,'bur_sil',burial2(1,1,1,issssil),2,0,iotype)
      CALL read_netcdf_var(ncid,'bur_clay',burial2(1,1,1,issster),2,0,iotype)
      CALL read_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaic',powtra2(1,1,1,ipowaic),2*ks,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'powc13',powtra2(1,1,1,ipowc13),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powc14',powtra2(1,1,1,ipowc14),2*ks,0,iotype)
#endif
      CALL read_netcdf_var(ncid,'powaal',powtra2(1,1,1,ipowaal),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaph',powtra2(1,1,1,ipowaph),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaox',powtra2(1,1,1,ipowaox),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'pown2',powtra2(1,1,1,ipown2),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powno3',powtra2(1,1,1,ipowno3),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powasi',powtra2(1,1,1,ipowasi),2*ks,0,iotype)

#ifdef DIFFAT 
!
! Read restart data : co2 diffusion
!
      CALL read_netcdf_var(ncid,'atmco2',atm(1,1,iatmco2),1,0,iotype)
      CALL read_netcdf_var(ncid,'atmo2',atm(1,1,iatmo2),1,0,iotype)
      CALL read_netcdf_var(ncid,'atmn2',atm(1,1,iatmn2),1,0,iotype)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'atmc13',atm(1,1,iatmc13),1,0,iotype)
      CALL read_netcdf_var(ncid,'atmc14',atm(1,1,iatmc14),1,0,iotype)
#endif

#endif
      IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_CLOSE(ncid)
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_CLOSE(ncid)
#endif
      ENDIF


! return tracer fields to ocean model (both timelevels); No unit
! conversion here, since tracers in the restart file are in 
! MICOM units (mol/kg) 
!--------------------------------------------------------------------
!
      trc(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)=locetra(:,:,:,:)


      RETURN
      END
