! Copyright (C) 2002  Ernst Maier-Reimer, S. Legutke, P. Wetzel
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, A. Moree
!                     M. Bentsen
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.


      SUBROUTINE AUFR_BGC(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,trc,               &
                          kplyear,kplmon,kplday,omask,rstfnm_ocn)
!******************************************************************************
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
!     J.Schwinger,      *GFI, Bergen*     2014-05-21
!     - adapted code for writing of two time level tracer 
!       and sediment fields
!
!     A.Moree,          *GFI, Bergen*   2018-04-12
!     - new version of carbon isotope code
!
!     J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
!     - added preformed and saturated DIC tracers
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - added cappability to restart c-isotopes from scratch (from
!       observed d13C and d14C). This is used if c-isotope fields are
!       not found in the restart file.
!     - consistently organised restart of CFC and natural tracers
!       from scratch, i.e. for the case that CFC and natural tracers are 
!       not found in the restart file.
!     - removed satn2o which is not needed to restart the model
!     - added sediment bypass preprocessor option
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-08-23
!     - added reading of atmosphere field for BOXATM
!
!     M. Bentsen,       *NORCE, Bergen*          2020-05-03
!     - changed ocean model from MICOM to BLOM
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
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmon*     - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *REAL*    *omask*      - land/ocean mask
!     *CHAR*    *rstfnm_ocn* - restart file name-informations
!
!
!**************************************************************************
      use mod_xc
      use netcdf
      USE mo_carbch
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc
      use mo_vgrid,     only: kbo
      USE mo_sedmnt,    only: sedhpl
      use mod_dia,      only : iotype
      use mo_intfcblom, only: sedlay2,powtra2,burial2,atm2
      implicit none
#ifdef PNETCDF
#include <pnetcdf.inc>
#endif
#include <mpif.h>      
      INTEGER          :: kpie,kpje,kpke,ntr,ntrbgc,itrbgc
      REAL             :: trc(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,2*kpke,ntr)
      REAL             :: omask(kpie,kpje)    
      INTEGER          :: kplyear,kplmon,kplday
      character(len=*) :: rstfnm_ocn

      ! Local variables
      REAL      :: locetra(kpie,kpje,2*kpke,nocetra) ! local array for reading 
      INTEGER   :: restyear                          !  year of restart file
      INTEGER   :: restmonth                         !  month of restart file
      INTEGER   :: restday                           !  day of restart file
      INTEGER   :: restdtoce                         !  time step number from bgc ocean file
      INTEGER   :: idate(5),i,j,k
      character :: rstfnm*256
      logical   :: lread_cfc,lread_nat,lread_iso,lread_atm
#ifdef cisonew
      REAL :: rco213,rco214,alpha14,beta13,beta14,d13C_atm,d14cat
#endif
      INTEGER ncid,ncstat,ncvarid
#ifdef PNETCDF
      integer*4 ,save :: info=MPI_INFO_NULL
      integer         :: mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),                      &
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
      character(len=3) :: stripestr
      character(len=9) :: stripestr2
      integer :: ierr,testio

      locetra(:,:,:,:) = 0.0
!
! Open netCDF data file
!
      testio=0
      IF(mnproc==1 .AND. IOTYPE==0) THEN

        i=1
        do while (rstfnm_ocn(i:i+7).ne.'.blom.r.' .AND.              &
     &            rstfnm_ocn(i:i+8).ne.'.micom.r.')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (io_stdo_bgc,*)                                    &
     &        'Could not generate restart file name!'
            call xchalt('(aufr_bgc)')
            stop '(aufr_bgc)'
          endif
        enddo
        if (rstfnm_ocn(i:i+7).eq.'.blom.r.') then
          rstfnm=rstfnm_ocn(1:i-1)//'.blom.rbgc.'//rstfnm_ocn(i+8:)
        else
          rstfnm=rstfnm_ocn(1:i-1)//'.micom.rbgc.'//rstfnm_ocn(i+9:)
        endif

        ncstat = NF90_OPEN(rstfnm,NF90_NOWRITE, ncid)
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
        do while (rstfnm_ocn(i:i+7).ne.'.blom.r.' .AND.              &
     &            rstfnm_ocn(i:i+8).ne.'.micom.r.')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (io_stdo_bgc,*)                                    &
     &        'Could not generate restart file name!'
            call xchalt('(aufr_bgc)')
            stop '(aufr_bgc)'
          endif
        enddo
        if (rstfnm_ocn(i:i+7).eq.'.blom.r.') then
          rstfnm=rstfnm_ocn(1:i-1)//'.blom.rbgc.'//rstfnm_ocn(i+8:)
        else
          rstfnm=rstfnm_ocn(1:i-1)//'.micom.rbgc.'//rstfnm_ocn(i+9:)
        endif
        write(stripestr,('(i3)')) 16
        write(stripestr2,('(i9)')) 1024*1024
        call mpi_info_create(info,ierr)
        call mpi_info_set(info,'romio_ds_read','disable',ierr)
        call mpi_info_set(info,'romio_ds_write','disable',ierr)
        call mpi_info_set(info,"striping_factor",stripestr,ierr)
        call mpi_info_set(info,"striping_unit",stripestr2,ierr)

        ncstat = NFMPI_OPEN(mpicomm,rstfnm,NF_NOWRITE,INFO, ncid)
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
!
      IF (mnproc.eq.1) THEN

         IF ( kplyear  .NE. restyear  ) THEN
            WRITE(io_stdo_bgc,*)                                     &
     &      'WARNING: restart years in oce/bgc are not the same : '  &
     &      ,kplyear,'/',restyear,' !!!'
         ENDIF

         IF ( kplmon .NE. restmonth ) THEN
            WRITE(io_stdo_bgc,*)                                     &
     &      'WARNING: restart months in oce/bgc are not the same : ' &
     &      ,kplmon,'/',restmonth,' !!!'
         ENDIF

         IF ( kplday   .NE. restday   ) THEN
            WRITE(io_stdo_bgc,*)                                     &
     &      'WARNING: restart days in oce/bgc are not the same : '   &
     &      ,kplday,'/',restday,' !!!'
         ENDIF

      ENDIF 

! Find out whether to restart CFCs
#ifdef CFC
      lread_cfc=.true.
      IF(IOTYPE==0) THEN
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'cfc11',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_cfc=.false.
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'cfc11',ncvarid)
        if(ncstat.ne.nf_noerr) lread_cfc=.false.
#endif
      ENDIF
      IF(mnproc==1 .and. .not. lread_cfc) THEN
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'AUFR_BGC info: CFC tracers not in restart file, '
        WRITE(io_stdo_bgc,*) ' CFCs initialised to zero.'
      ENDIF
#endif

! Find out whether to restart natural tracers
#ifdef natDIC
      lread_nat=.true.
      IF(IOTYPE==0) THEN
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'natsco212',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_nat=.false.
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'natsco212',ncvarid)
        if(ncstat.ne.nf_noerr) lread_nat=.false.
#endif
      ENDIF
      IF(mnproc==1 .and. .not. lread_nat) THEN
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'AUFR_BGC info: natural tracers not in restart file. '
        WRITE(io_stdo_bgc,*) ' Initialising natural tracers with their non-natural '
        WRITE(io_stdo_bgc,*) ' counterpart.'
      ENDIF
#endif

! Find out whether to restart marine carbon isotopes
#ifdef cisonew
      lread_iso=.true.
      IF(IOTYPE==0) THEN
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'sco213',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_iso=.false.
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'sco213',ncvarid)
        if(ncstat.ne.nf_noerr) lread_iso=.false.
#endif
      ENDIF
      IF(mnproc==1 .and. .not. lread_iso) THEN
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'AUFR_BGC info: carbon isotopes not in restart file. '
        WRITE(io_stdo_bgc,*) ' Initialising carbon isotopes from scratch '
      ENDIF
#endif

! Find out whether to restart atmosphere
#if defined(BOXATM)
      lread_atm=.true.
      IF(IOTYPE==0) THEN
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'atmco2',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_atm=.false.
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'atmco2',ncvarid)
        if(ncstat.ne.nf_noerr) lread_atm=.false.
#endif
      ENDIF
      IF(mnproc==1 .and. .not. lread_atm) THEN
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'AUFR_BGC info: atmosphere fields not in restart file. '
        WRITE(io_stdo_bgc,*) ' Initialising atmosphere from scratch '
      ENDIF
#endif

!
! Read restart data : ocean aquateous tracer
!                
      CALL read_netcdf_var(ncid,'sco212',locetra(1,1,1,isco212),2*kpke,0,iotype)
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
      CALL read_netcdf_var(ncid,'calciu',locetra(1,1,1,icalc),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'opal',locetra(1,1,1,iopal),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'n2o',locetra(1,1,1,ian2o),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'dms',locetra(1,1,1,idms),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'fdust',locetra(1,1,1,ifdust),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'iron',locetra(1,1,1,iiron),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefo2',locetra(1,1,1,iprefo2),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefpo4',locetra(1,1,1,iprefpo4),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefalk',locetra(1,1,1,iprefalk),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'prefdic',locetra(1,1,1,iprefdic),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'dicsat',locetra(1,1,1,idicsat),2*kpke,0,iotype)

#ifdef cisonew
      IF(lread_iso) THEN
      CALL read_netcdf_var(ncid,'sco213',locetra(1,1,1,isco213),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'sco214',locetra(1,1,1,isco214),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'doc13',locetra(1,1,1,idoc13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'doc14',locetra(1,1,1,idoc14),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'phyto13',locetra(1,1,1,iphy13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'phyto14',locetra(1,1,1,iphy14),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'grazer13',locetra(1,1,1,izoo13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'grazer14',locetra(1,1,1,izoo14),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'poc13',locetra(1,1,1,idet13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'poc14',locetra(1,1,1,idet14),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'calciu13',locetra(1,1,1,icalc13),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'calciu14',locetra(1,1,1,icalc14),2*kpke,0,iotype)
      ENDIF
#endif
#ifdef AGG
      CALL read_netcdf_var(ncid,'snos',locetra(1,1,1,inos),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'adust',locetra(1,1,1,iadust),2*kpke,0,iotype)
#endif /*AGG*/
#ifdef CFC
      IF(lread_cfc) THEN
      CALL read_netcdf_var(ncid,'cfc11',locetra(1,1,1,icfc11),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'cfc12',locetra(1,1,1,icfc12),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'sf6',locetra(1,1,1,isf6),2*kpke,0,iotype)
      ENDIF
#endif
#ifdef natDIC
      IF(lread_nat) THEN
      CALL read_netcdf_var(ncid,'natsco212',locetra(1,1,1,inatsco212),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'natalkali',locetra(1,1,1,inatalkali),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'natcalciu',locetra(1,1,1,inatcalc),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'nathi',nathi(1,1,1),kpke,0,iotype)
      ELSE
      CALL read_netcdf_var(ncid,'sco212',locetra(1,1,1,inatsco212),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'alkali',locetra(1,1,1,inatalkali),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'calciu',locetra(1,1,1,inatcalc),2*kpke,0,iotype)
      CALL read_netcdf_var(ncid,'hi',nathi(1,1,1),kpke,0,iotype)
      ENDIF
#endif
!
! Read restart data : diagnostic ocean fields (needed for bit to bit reproducability)
!
      CALL read_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0,iotype)
      CALL read_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0,iotype)
      CALL read_netcdf_var(ncid,'co2star',co2star(1,1,1),kpke,0,iotype)
      CALL read_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0,iotype)

!
! Read restart data : sediment variables.
!
#ifndef sedbypass
      CALL read_netcdf_var(ncid,'ssso12',sedlay2(1,1,1,issso12),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'sssc12',sedlay2(1,1,1,isssc12),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'ssssil',sedlay2(1,1,1,issssil),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'ssster',sedlay2(1,1,1,issster),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'bur_o12',burial2(1,1,1,issso12),2,0,iotype)
      CALL read_netcdf_var(ncid,'bur_c12',burial2(1,1,1,isssc12),2,0,iotype)
      CALL read_netcdf_var(ncid,'bur_sil',burial2(1,1,1,issssil),2,0,iotype)
      CALL read_netcdf_var(ncid,'bur_clay',burial2(1,1,1,issster),2,0,iotype)
      CALL read_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaic',powtra2(1,1,1,ipowaic),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaal',powtra2(1,1,1,ipowaal),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaph',powtra2(1,1,1,ipowaph),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powaox',powtra2(1,1,1,ipowaox),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'pown2',powtra2(1,1,1,ipown2),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powno3',powtra2(1,1,1,ipowno3),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powasi',powtra2(1,1,1,ipowasi),2*ks,0,iotype)
#ifdef cisonew
      IF(lread_iso) THEN
      ! Burial fields for c-isotopes still missing
      CALL read_netcdf_var(ncid,'ssso13',sedlay2(1,1,1,issso13),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'ssso14',sedlay2(1,1,1,issso14),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'sssc13',sedlay2(1,1,1,isssc13),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'sssc14',sedlay2(1,1,1,isssc14),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powc13',powtra2(1,1,1,ipowc13),2*ks,0,iotype)
      CALL read_netcdf_var(ncid,'powc14',powtra2(1,1,1,ipowc14),2*ks,0,iotype)
      ENDIF
#endif
#endif

!
! Read restart data: atmosphere
!
#if defined(BOXATM)
      IF(lread_atm) THEN
      CALL read_netcdf_var(ncid,'atmco2',atm2(1,1,1,iatmco2),2,0,iotype)
      CALL read_netcdf_var(ncid,'atmo2',atm2(1,1,1,iatmo2),2,0,iotype)
      CALL read_netcdf_var(ncid,'atmn2',atm2(1,1,1,iatmn2),2,0,iotype)
#ifdef cisonew
      IF(lread_iso) THEN
      CALL read_netcdf_var(ncid,'atmc13',atm2(1,1,1,iatmc13),2,0,iotype)
      CALL read_netcdf_var(ncid,'atmc14',atm2(1,1,1,iatmc14),2,0,iotype)
      ELSE
      ! If atm isotopes are not in restart but boxatm is on, calculate initial value using atmco2
      ! that is just read in from restart files. Normalize atmc14 using beleg c14fac.
      DO j=1,kpje
      DO i=1,kpie
        beta13 = (prei13/1000.)+1.
        alpha14 = 2.*(prei13+25.)
        d14cat  = (prei14+alpha14)/(1.-alpha14/1000.)
        atm(i,j,iatmc13) = beta13*re1312*atm2(i,j,1,iatmco2)/(1.+beta13*re1312)
        atm(i,j,iatmc14) = ((d14cat/1000.)+1.)*re14to*atm2(i,j,1,iatmco2)/c14fac
      ENDDO
      ENDDO
      ! Copy the isotope atmosphere fields into both timelevels of atm2.
      atm2(:,:,1,iatmc13) = atm(:,:,iatmc13)
      atm2(:,:,2,iatmc13) = atm(:,:,iatmc13)
      atm2(:,:,1,iatmc14) = atm(:,:,iatmc14)
      atm2(:,:,2,iatmc14) = atm(:,:,iatmc14)
      ENDIF
#endif
#ifdef natDIC
      CALL read_netcdf_var(ncid,'atmnco2',atm2(1,1,1,iatmnco2),2,0,iotype)
#endif
      ELSE
      ! If atmosphere field is not in restart, copy the atmosphere field
      ! (initialised in beleg.F90) into both timelevels of atm2.
      atm2(:,:,1,:) = atm(:,:,:)
      atm2(:,:,2,:) = atm(:,:,:)
      ENDIF
#endif

      IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_CLOSE(ncid)
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_CLOSE(ncid)
#endif
      ENDIF


#ifdef cisonew
      IF(.NOT. lread_iso) THEN
      ! If carbon isotope fields are not read from restart file, copy the d13C
      ! d14C fields (initialised in beleg.F90) into both timelevels of locetra.
      locetra(:,:,1:kpke,       isco213)=ocetra(:,:,:,isco213)
      locetra(:,:,kpke+1:2*kpke,isco213)=ocetra(:,:,:,isco213)
      locetra(:,:,1:kpke,       isco214)=ocetra(:,:,:,isco214)
      locetra(:,:,kpke+1:2*kpke,isco214)=ocetra(:,:,:,isco214)
      ! Initialise 13C and 14C fields in the same way as in beleg.F90
      DO k=1,2*kpke
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j) .GT. 0.5) THEN
        ! 13C is read in as delta13C, convert to 13C using model restart total C
        beta13=locetra(i,j,k,isco213)/1000.+1.
        locetra(i,j,k,isco213)=locetra(i,j,k,isco212)*beta13*re1312/(1.+beta13*re1312)

        ! 14C is read in as delta14C, convert to 14C using model restart total C, 
        ! normalize 14C by c14fac to prevent numerical errors
        beta14=locetra(i,j,k,isco214)/1000.+1.
        locetra(i,j,k,isco214)=locetra(i,j,k,isco212)*beta14*re14to/c14fac
        
        ! Initialise the remaining 13C and 14C fields, using the restart isco212 field
        rco213=locetra(i,j,k,isco213)/(locetra(i,j,k,isco212)+safediv)
        rco214=locetra(i,j,k,isco214)/(locetra(i,j,k,isco212)+safediv)
        locetra(i,j,k,idoc13)=locetra(i,j,k,idoc)*rco213*bifr13
        locetra(i,j,k,idoc14)=locetra(i,j,k,idoc)*rco214*bifr14
        locetra(i,j,k,iphy13)=locetra(i,j,k,iphy)*rco213*bifr13
        locetra(i,j,k,iphy14)=locetra(i,j,k,iphy)*rco214*bifr14
        locetra(i,j,k,izoo13)=locetra(i,j,k,izoo)*rco213*bifr13
        locetra(i,j,k,izoo14)=locetra(i,j,k,izoo)*rco214*bifr14
        locetra(i,j,k,idet13)=locetra(i,j,k,idet)*rco213*bifr13
        locetra(i,j,k,idet14)=locetra(i,j,k,idet)*rco214*bifr14
        locetra(i,j,k,icalc13)=locetra(i,j,k,icalc)*rco213
        locetra(i,j,k,icalc14)=locetra(i,j,k,icalc)*rco214

        ENDIF
      ENDDO
      ENDDO
      ENDDO
#ifndef sedbypass
      ! Burial fields for c-isotopes still missing
      DO  k=1,2*ks
      DO  j=1,kpje
      DO  i=1,kpie 
        IF(omask(i,j) .GT. 0.5) THEN
        rco213=ocetra(i,j,kbo(i,j),isco213)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
        rco214=ocetra(i,j,kbo(i,j),isco214)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
        powtra2(i,j,k,ipowc13)=powtra2(i,j,k,ipowc)*rco213*bifr13
        powtra2(i,j,k,ipowc14)=powtra2(i,j,k,ipowc)*rco214*bifr14
        sedlay2(i,j,k,issso13)=sedlay2(i,j,k,issso)*rco213*bifr13
        sedlay2(i,j,k,issso14)=sedlay2(i,j,k,issso)*rco214*bifr14
        sedlay2(i,j,k,isssc13)=sedlay2(i,j,k,isssc)*rco213
        sedlay2(i,j,k,isssc14)=sedlay2(i,j,k,isssc)*rco214
        ENDIF
      ENDDO
      ENDDO
      ENDDO
#endif
      ENDIF ! .NOT. lread_iso
#endif

! return tracer fields to ocean model (both timelevels); No unit
! conversion here, since tracers in the restart file are in 
! BLOM units (mol/kg) 
!--------------------------------------------------------------------
!
      trc(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)=locetra(:,:,:,:)


      RETURN
      END
