! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
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


      SUBROUTINE AUFW_BGC(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,trc,              &
                          kplyear,kplmon,kplday,kpldtoce,omask,rstfnm_ocn)
!******************************************************************************
!
!**** *AUFW_BGC* - write marine bgc restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - extra SBR for writing bgc data to the restart file.
!     S.Legutke,        *MPI-MaD, HH*    15.08.01
!     - netCDF version (cond.comp. PNETCDF)
!     - chemcm is multiplied with layer-dependent constant in order
!       to be displayable by ncview. It is not read in AUFR_BGC!
!
!     J.Schwinger,      *GFI, Bergen*     2013-10-21
!     - tracer field is passed from ocean model for writing now
!     - removed writing of chemcm and ak* fields
!     - code cleanup, removed preprocessor option "PNETCDF"
!
!     J.Schwinger,      *GFI, Bergen*     2014-05-21
!     - adapted code for writing of two time level tracer and 
!       sediment fields
!
!     A.Moree,          *GFI, Bergen*   2018-04-12
!     - new version of carbon isotope code
!
!     J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
!     - added preformed and saturated DIC tracers
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - removed satn2o which is not needed to restart the model
!     - added sediment bypass preprocessor option
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-08-23
!     - added writing of atmosphere field for BOXATM
!
!     M. Bentsen,       *NORCE, Bergen*          2020-05-03
!     - changed ocean model from MICOM to BLOM
!
!     Purpose
!     -------
!     Write restart data for continuation of interrupted integration.
!
!     Method
!     -------
!     The bgc data are written to an extra file, other than the ocean data.
!     The time stamp of the bgc restart file (idate) is taken from the
!     ocean time stamp through the SBR parameter list. The only time
!     control variable proper to the bgc is the time step number (idate(5)). 
!     It can differ from that of the ocean (idate(4)) by the difference
!     of the offsets of restart files.
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
!     *REAL*    *trc*        - initial/restart tracer field to be passed from the 
!                              ocean model [mol/kg]
!     *REAL*    *sedlay2*    - initial/restart sediment (two time levels) field
!     *REAL*    *powtra2*    - initial/restart pore water tracer (two time levels) field
!     *REAL*    *burial2*    - initial/restart sediment burial (two time levels) field
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmon*     - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *INTEGER* *kpldtoce*   - step  in ocean restart date
!     *REAL*    *omask*      - land/ocean mask
!     *CHAR*    *rstfnm_ocn* - restart file name-informations
!
!**************************************************************************
      USE netcdf
      USE mo_carbch
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc 
      USE mo_sedmnt,    only: sedhpl
      use mod_xc,       only: nbdy,itdm,jtdm,mnproc,xchalt
      use mo_intfcblom, only: sedlay2,powtra2,burial2,atm2
      use mod_dia

      implicit none

      INTEGER,          intent(in) :: kpie,kpje,kpke,ntr,ntrbgc,itrbgc
      REAL,             intent(in) :: trc(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,2*kpke,ntr)
      REAL,             intent(in) :: omask(kpie,kpje)    
      INTEGER,          intent(in) :: kplyear,kplmon,kplday,kpldtoce
      character(len=*), intent(in) :: rstfnm_ocn

      ! Local variables
      INTEGER             :: i,j
      REAL                :: locetra(kpie,kpje,2*kpke,nocetra)
      CHARACTER(LEN=256)  :: rstfnm

      ! Variables for netcdf
      INTEGER             :: ncid,ncvarid,ncstat,ncoldmod,ncdimst(4)
      INTEGER             :: nclatid,nclonid,nclevid,nclev2id,ncksid,ncks2id,nctlvl2id
      INTEGER             :: idate(5),ierr,testio
      REAL                :: rmissing
      character(len=3)    :: stripestr
      character(len=9)    :: stripestr2

#ifdef PNETCDF
      integer*4 ,save     :: info=MPI_INFO_NULL
      integer             :: mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),                          &
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif

! pass tracer fields in from ocean model, note that both timelevels 
! are passed into the local array locetra; No unit conversion here, 
! tracers in the restart file are written in mol/kg 
!--------------------------------------------------------------------
!
      testio=0
      locetra(:,:,:,:)=trc(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)

      idate(1) = kplyear
      idate(2) = kplmon
      idate(3) = kplday
      idate(4) = kpldtoce
      idate(5) = ldtbgc
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Writing restart file at date :'              &
     &,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)
      ENDIF

      rmissing = rmasko

!
! Open netCDF data file
!
      IF(mnproc==1 .AND. IOTYPE==0) THEN

      i=1
      do while (rstfnm_ocn(i:i+7).ne.'.blom.r.')
        i=i+1
        if (i+7.gt.len(rstfnm_ocn)) then
          write (io_stdo_bgc,*)                                      &
     &      'Could not generate restart file name!'
          call xchalt('(aufw_bgc)')
          stop '(aufw_bgc)'
        endif
      enddo
      rstfnm=rstfnm_ocn(1:i-1)//'.blom.rbgc.'//rstfnm_ocn(i+8:)

      write(io_stdo_bgc,*) 'BGC RESTART   ',rstfnm
      ncstat = NF90_CREATE(rstfnm,NF90_64BIT_OFFSET,ncid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF1)')
               stop '(AUFW: Problem with netCDF1)'
      ENDIF
      ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      testio=1
      i=1
      do while (rstfnm_ocn(i:i+7).ne.'.blom.r.')
        i=i+1
        if (i+7.gt.len(rstfnm_ocn)) then
          write (io_stdo_bgc,*)                                      &
     &      'Could not generate restart file name!'
          call xchalt('(aufw_bgc)')
          stop '(aufw_bgc)'
        endif
      enddo
      rstfnm=rstfnm_ocn(1:i-1)//'.blom.rbgc.'//rstfnm_ocn(i+8:)

      IF(mnproc==1) write(io_stdo_bgc,*) 'BGC RESTART   ',rstfnm
      write(stripestr,('(i3)')) 16
      write(stripestr2,('(i9)')) 1024*1024
      call mpi_info_create(info,ierr)
      call mpi_info_set(info,'romio_ds_read','disable',ierr)
      call mpi_info_set(info,'romio_ds_write','disable',ierr)
      call mpi_info_set(info,"striping_factor",stripestr,ierr)
      call mpi_info_set(info,"striping_unit",stripestr2,ierr)
      ncstat = NFMPI_CREATE(mpicomm,rstfnm,                         &
     &       IOR(nf_clobber,nf_64bit_offset),info,ncid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF1)')
               stop '(AUFW: Problem with netCDF1)'
      ENDIF
#endif

      if(testio .eq. 0) then
      CALL xchalt('(AUFW: Problem with namelist iotype)')
                    stop '(AUFW: Problem with namelist iotype)'
      endif

      ENDIF
!
! Define dimension
! ----------------------------------------------------------------------    
!
      IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_DEF_DIM(ncid, 'lon', itdm, nclonid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF2)')
               stop '(AUFW: Problem with netCDF2)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'lat', jtdm, nclatid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF3)')
               stop '(AUFW: Problem with netCDF3)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF4)')
               stop '(AUFW: Problem with netCDF4)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'depth2', 2*kpke, nclev2id)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF4)')
               stop '(AUFW: Problem with netCDF4)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'nks', ks, ncksid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF7)')
               stop '(AUFW: Problem with netCDF7)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'nks2', 2*ks, ncks2id)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF7)')
               stop '(AUFW: Problem with netCDF7)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'tlvl2', 2, nctlvl2id)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF7)')
               stop '(AUFW: Problem with netCDF7)'
      ENDIF
      ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      clen=itdm
      ncstat = NFMPI_DEF_DIM(ncid, 'lon', clen, nclonid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF2)')
               stop '(AUFW: Problem with PnetCDF2)'
      ENDIF

      clen=jtdm
      ncstat = NFMPI_DEF_DIM(ncid, 'lat', clen, nclatid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF3)')
               stop '(AUFW: Problem with PnetCDF3)'
      ENDIF

      clen=kpke
      ncstat = NFMPI_DEF_DIM(ncid, 'depth', clen, nclevid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF4)')
               stop '(AUFW: Problem with PnetCDF4)'
      ENDIF

      clen=2*kpke
      ncstat = NFMPI_DEF_DIM(ncid, 'depth2', clen, nclev2id)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF4)')
               stop '(AUFW: Problem with PnetCDF4)'
      ENDIF

      clen=ks
      ncstat = NFMPI_DEF_DIM(ncid, 'nks', clen, ncksid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF5)')
               stop '(AUFW: Problem with PnetCDF5)'
      ENDIF

      clen=2*ks
      ncstat = NFMPI_DEF_DIM(ncid, 'nks2', clen, ncks2id)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF6)')
               stop '(AUFW: Problem with PnetCDF6)'
      ENDIF

      clen=2
      ncstat = NFMPI_DEF_DIM(ncid, 'tlvl2', clen, nctlvl2id)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF7)')
               stop '(AUFW: Problem with PnetCDF7)'
      ENDIF
#endif
      ENDIF

!
! Define global attributes
! ----------------------------------------------------------------------    
!
      IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'title'               &
     &, 'Restart data for marine bgc modules') 
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'history'             &
     &, 'Restart data for marine bgc modules')
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'conventions'         &
     &,'COARDS')
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'source'              &
     &, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF10)')
               stop '(AUFW: Problem with netCDF10)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'date', idate)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF11)')
               stop '(AUFW: Problem with netCDF11)'
      ENDIF

!PNETCDF
      ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      clen=len('Restart data for marine bgc modules')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'             &
     &, clen,'Restart data for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF9)')
               stop '(AUFW: Problem with PnetCDF9)'
      ENDIF
      clen=len('Restart data for marine bgc modules')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'history'           &
     &, clen,'Restart data for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF9)')
               stop '(AUFW: Problem with PnetCDF9)'
      ENDIF
      clen=6
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'conventions'       &
     &,clen, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF9)')
               stop '(AUFW: Problem with PnetCDF9)'
      ENDIF
      clen=len('Marine bgc model output HOPC68/grob')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'            &
     &,clen, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF10)')
               stop '(AUFW: Problem with PnetCDF10)'
      ENDIF
      clen=5
      ncstat = NFMPI_PUT_ATT_INT(ncid, NF_GLOBAL, 'date',              &
     &                           nf_int, clen, idate)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF11)')
               stop '(AUFW: Problem with netCDF11)'

      ENDIF
#endif
      ENDIF
! 
! Define variables : advected ocean tracer
! ----------------------------------------------------------------------    
!
      IF((mnproc==1 .AND. IOTYPE==0) .OR. IOTYPE==1) THEN
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = nclev2id
        ncdimst(4) = 0
      ENDIF

      CALL NETCDF_DEF_VARDB(ncid,6,'sco212',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',13, 'Dissolved CO2',rmissing,10,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'alkali',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',10,'Alkalinity',rmissing,11,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'phosph',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',19,'Dissolved phosphate',rmissing,12,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'oxygen',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',16,'Dissolved oxygen',                             &
          rmissing,13,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'gasnit',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',21,'Gaseous nitrogen (N2)',                        &
          rmissing,14,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'ano3',3,ncdimst,ncvarid,            &
     &    6,'mol/kg',17,'Dissolved nitrate',                            &
          rmissing,15,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'silica',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',22,'Silicid acid (Si(OH)4)',                       &
          rmissing,16,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'doc',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',24,'Dissolved organic carbon',                     &
     &    rmissing,17,io_stdo_bgc) 

      CALL NETCDF_DEF_VARDB(ncid,3,'poc',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',25,'Particulate organic carbon',                   &
     &    rmissing,18,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'phyto',3,ncdimst,ncvarid,           &
     &    7,'molP/kg',27,'Phytoplankton concentration',                 &
     &    rmissing,19,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'grazer',3,ncdimst,ncvarid,          &
     &    7,'molP/kg',25,'Zooplankton concentration',                   &
     &    rmissing,20,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'calciu',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',17,'Calcium carbonate',                            &
     &    rmissing,21,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'opal',3,ncdimst,ncvarid,            &
     &    6,'mol/kg',15,'Biogenic silica',                              &
     &    rmissing,22,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'n2o',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',12,'laughing gas',                                 &
     &    rmissing,23,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'dms',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',15 ,'DiMethylSulfide',                             &
     &    rmissing,24,io_stdo_bgc)            

      CALL NETCDF_DEF_VARDB(ncid,5,'fdust',3,ncdimst,ncvarid,           &
     &    5,'kg/kg',19,'Non-aggregated dust',                           &
     &    rmissing,25,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'iron',3,ncdimst,ncvarid,            &
     &    6,'mol/kg',14,'Dissolved iron',                               &
     &    rmissing,26,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'prefo2',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',16,'Preformed oxygen',                             &
          rmissing,27,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'prefpo4',3,ncdimst,ncvarid,         &
     &    6,'mol/kg',19,'Preformed phosphate',                          &
          rmissing,28,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'prefalk',3,ncdimst,ncvarid,         &
     &    6,'mol/kg',20,'Preformed alkalinity',                         &
          rmissing,29,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'prefdic',3,ncdimst,ncvarid,         &
     &    6,'mol/kg',13,'Preformed dic',                                &
          rmissing,30,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'dicsat',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',13,'Saturated dic',                                &
          rmissing,31,io_stdo_bgc)

#ifdef cisonew
      CALL NETCDF_DEF_VARDB(ncid,6,'sco213',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',15, 'Dissolved CO213',rmissing,32,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sco214',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',15, 'Dissolved CO214',rmissing,33,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'doc13',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',24,'Dissolved organic carb13',                     &
     &    rmissing,34,io_stdo_bgc) 

      CALL NETCDF_DEF_VARDB(ncid,5,'doc14',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',24,'Dissolved organic carb14',                     &
     &    rmissing,35,io_stdo_bgc) 

      CALL NETCDF_DEF_VARDB(ncid,5,'poc13',3,ncdimst,ncvarid,           &
     &    7,'molC/kg',28,'Particulate organic carbon13',                &
     &    rmissing,36,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'poc14',3,ncdimst,ncvarid,           &
     &    7,'molC/kg',28,'Particulate organic carbon14',                &
     &    rmissing,37,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'phyto13',3,ncdimst,ncvarid,         &
     &    7,'molP/kg',27,'Phytoplankton concentr. 13c',                 &
     &    rmissing,38,io_stdo_bgc)

       CALL NETCDF_DEF_VARDB(ncid,7,'phyto14',3,ncdimst,ncvarid,        &
     &    7,'molP/kg',27,'Phytoplankton concentr. 14c',                 &
     &    rmissing,39,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'grazer13',3,ncdimst,ncvarid,        &
     &    7,'molP/kg',25,'Zooplankton concentr. 13c',                   &
     &    rmissing,40,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'grazer14',3,ncdimst,ncvarid,        &
     &    7,'molP/kg',25,'Zooplankton concentr. 14c',                   &
     &    rmissing,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'calciu13',3,ncdimst,ncvarid,        &
     &    7,'molC/kg',19,'Calcium carbonate13',                         &
     &    rmissing,42,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'calciu14',3,ncdimst,ncvarid,        &
     &    7,'molC/kg',19,'Calcium carbonate14',                         &
     &    rmissing,43,io_stdo_bgc)
#endif
#ifdef AGG
      CALL NETCDF_DEF_VARDB(ncid,4,'snos',3,ncdimst,ncvarid,            &
     &    3,'1/g',38,'marine snow aggregates per g sea water',          &
     &    rmissing,44,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'adust',3,ncdimst,ncvarid,           &
     &    4,'g/kg',15,'Aggregated dust',                                &
     &    rmissing,45,io_stdo_bgc)
#endif
#ifdef CFC
      CALL NETCDF_DEF_VARDB(ncid,5,'cfc11',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',5,'CFC11',                                         &
     &    rmissing,47,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'cfc12',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',5,'CFC12',                                         &
     &    rmissing,48,io_stdo_bgc)     

      CALL NETCDF_DEF_VARDB(ncid,3,'sf6',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',4,'SF-6',                                          &
     &    rmissing,49,io_stdo_bgc)     
#endif
#ifdef natDIC
      CALL NETCDF_DEF_VARDB(ncid,9,'natsco212',3,ncdimst,ncvarid,       &
     &   6,'mol/kg',21, 'Natural dissolved CO2',rmissing,50,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,9,'natalkali',3,ncdimst,ncvarid,       &
     &    6,'mol/kg',18,'Natural alkalinity',rmissing,51,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,9,'natcalciu',3,ncdimst,ncvarid,       &
     &    6,'mol/kg',25,'Natural calcium carbonate',                    &
     &    rmissing,52,io_stdo_bgc)
#endif

!
! Define variables : diagnostic ocean fields
! ----------------------------------------------------------------------    
!
      IF((mnproc==1 .AND. IOTYPE==0) .OR. IOTYPE==1) THEN
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = nclevid
        ncdimst(4) = 0
      ENDIF

      CALL NETCDF_DEF_VARDB(ncid,2,'hi',3,ncdimst,ncvarid,              &
     &    6,'mol/kg',26,'Hydrogen ion concentration',                   &
     &    rmissing,60,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'co3',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',25,'Dissolved carbonate (CO3)',                    &
     &    rmissing,61,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'co2star',3,ncdimst,ncvarid,         &
     &    6,'mol/kg',20,'Dissolved CO2 (CO2*)',                         &
     &    rmissing,62,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'satoxy',3,ncdimst,ncvarid,          &
     &    9,'xxxxxxxxx',9 ,'xxxxxxxxx',  &
     &    rmissing,63,io_stdo_bgc)

#ifdef natDIC
      CALL NETCDF_DEF_VARDB(ncid,5,'nathi',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',34,'Natural hydrogen ion concentration',           &
     &    rmissing,64,io_stdo_bgc)
#endif


!
! Define variables : sediment
! ----------------------------------------------------------------------    
!
#ifndef sedbypass
      IF((mnproc==1 .AND. IOTYPE==0) .OR. IOTYPE==1) THEN
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = ncks2id
        ncdimst(4) = 0
      ENDIF

      CALL NETCDF_DEF_VARDB(ncid,6,'ssso12',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',35,'Sediment accumulated organic carbon',       &
     &    rmissing,70,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc12',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',38,'Sediment accumulated calcium carbonate',    &
     &    rmissing,71,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'ssssil',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',25,'Sediment accumulated opal',                 &
     &    rmissing,72,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'ssster',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',25,'Sediment accumulated clay',                 &
     &    rmissing,73,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaic',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',23,'Sediment pore water CO2',                   &
     &    rmissing,74,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaal',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',30,'Sediment pore water alkalinity',            &
     &    rmissing,75,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaph',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',29,'Sediment pore water phosphate',             &
     &    rmissing,76,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaox',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',26,'Sediment pore water oxygen',                &
     &    rmissing,77,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'pown2',3,ncdimst,ncvarid,           &
     &    9,'kmol/m**3',36,'Sediment pore water gaseous nitrogen',      &
     &    rmissing,78,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powno3',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',33,'Sediment pore water nitrate (NO3)',         &
     &    rmissing,79,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powasi',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',42,'Sediment pore water silicid acid (Si(OH)4)',&
     &    rmissing,80,io_stdo_bgc)

#ifdef cisonew
      CALL NETCDF_DEF_VARDB(ncid,6,'ssso13',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',37,'Sediment accumulated organic carbon13',     &
     &    rmasks,81,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'ssso14',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',37,'Sediment accumulated organic carbon14',     &
     &    rmasks,82,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc13',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',40,'Sediment accumulated calcium carbonate13',  &
     &    rmasks,83,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc14',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',40,'Sediment accumulated calcium carbonate14',  &
     &    rmasks,84,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powc13',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',25,'Sediment pore water DIC13',                 &
     &    rmasks,85,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powc14',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**3',25,'Sediment pore water DIC14',                 &
     &    rmasks,86,io_stdo_bgc)
#endif

      IF((mnproc==1 .AND. IOTYPE==0) .OR. IOTYPE==1) THEN
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = ncksid
        ncdimst(4) = 0
      ENDIF

      CALL NETCDF_DEF_VARDB(ncid,6,'sedhpl',3,ncdimst,ncvarid,          &
     &    9,'kmol/m**2',34,'Sediment accumulated hydrogen ions',        &
     &    rmissing,87,io_stdo_bgc)

!
! Define variables : sediment burial
! ----------------------------------------------------------------------    
!
      IF((mnproc==1 .AND. IOTYPE==0) .OR. IOTYPE==1) THEN
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = nctlvl2id
        ncdimst(4) = 0
      ENDIF

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_o12',3,ncdimst,ncvarid,         &
     &    9,'kmol/m**2',30,'Burial layer of organic carbon',            &
     &    rmissing,90,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_c12',3,ncdimst,ncvarid,         &
     &    9,'kmol/m**2',33,'Burial layer of calcium carbonate',         &
     &    rmissing,91,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_sil',3,ncdimst,ncvarid,         &
     &    9,'kmol/m**2',20,'Burial layer of opal',                      &
     &    rmissing,92,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'bur_clay',3,ncdimst,ncvarid,        &
     &    9,'kmol/m**2',20,'Burial layer of clay',                      &
     &    rmissing,93,io_stdo_bgc)

#endif /* sedbypass */

!
! Define variables: atmosphere
! ----------------------------------------------------------------------    
!
#if defined(BOXATM)
    IF((mnproc==1 .AND. IOTYPE==0) .OR. IOTYPE==1) THEN
      ncdimst(1) = nclonid
      ncdimst(2) = nclatid
      ncdimst(3) = nctlvl2id
      ncdimst(4) = 0
    ENDIF
     
      CALL NETCDF_DEF_VARDB(ncid,6,'atmco2',3,ncdimst,ncvarid,          &
     &    3,'ppm',15,'atmospheric CO2',                                 &
     &    rmissing,101,io_stdo_bgc)         
     
      CALL NETCDF_DEF_VARDB(ncid,5,'atmo2',3,ncdimst,ncvarid,           &
     &    3,'ppm',14,'atmospheric O2',                                  &
     &    rmissing,102,io_stdo_bgc)    
     
      CALL NETCDF_DEF_VARDB(ncid,5,'atmn2',3,ncdimst,ncvarid,           &
     &    3,'ppm',14,'atmospheric N2',                                  &
     &    rmissing,103,io_stdo_bgc)    

#ifdef cisonew
      CALL NETCDF_DEF_VARDB(ncid,6,'atmc13',3,ncdimst,ncvarid,          &
     &    3,'ppm',17,'atmospheric 13CO2',                               &
     &    rmissing,104,io_stdo_bgc)      
      CALL NETCDF_DEF_VARDB(ncid,6,'atmc14',3,ncdimst,ncvarid,          &
     &    3,'ppm',17,'atmospheric 14CO2',                               &
     &    rmissing,105,io_stdo_bgc)      
#endif
#ifdef natDIC
      CALL NETCDF_DEF_VARDB(ncid,7,'atmnco2',3,ncdimst,ncvarid,         &
     &    3,'ppm',23,'natural atmospheric CO2',                         &
     &    rmissing,106,io_stdo_bgc)      
#endif
#endif
      IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_ENDDEF(ncid)
      
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF00)')
               stop '(AUFW: Problem with netCDF00)'
      ENDIF

!
! Set fill mode
! ----------------------------------------------------------------------    
!
      ncstat = NF90_SET_FILL(ncid,NF90_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF97)')
               stop '(AUFW: Problem with netCDF97)'
      ENDIF



      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_ENDDEF(ncid)

      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF00)')
               stop '(AUFW: Problem with PnetCDF00)'
      ENDIF
#endif
      ENDIF

!
! Write restart data : ocean aquateous tracer
!--------------------------------------------------------------------
!
      CALL write_netcdf_var(ncid,'sco212',locetra(1,1,1,isco212),2*kpke,0)
      CALL write_netcdf_var(ncid,'alkali',locetra(1,1,1,ialkali),2*kpke,0)
      CALL write_netcdf_var(ncid,'phosph',locetra(1,1,1,iphosph),2*kpke,0)
      CALL write_netcdf_var(ncid,'oxygen',locetra(1,1,1,ioxygen),2*kpke,0)
      CALL write_netcdf_var(ncid,'gasnit',locetra(1,1,1,igasnit),2*kpke,0)
      CALL write_netcdf_var(ncid,'ano3',locetra(1,1,1,iano3),2*kpke,0)
      CALL write_netcdf_var(ncid,'silica',locetra(1,1,1,isilica),2*kpke,0)
      CALL write_netcdf_var(ncid,'doc',locetra(1,1,1,idoc),2*kpke,0)
      CALL write_netcdf_var(ncid,'poc',locetra(1,1,1,idet),2*kpke,0)
      CALL write_netcdf_var(ncid,'phyto',locetra(1,1,1,iphy),2*kpke,0)   
      CALL write_netcdf_var(ncid,'grazer',locetra(1,1,1,izoo),2*kpke,0)
      CALL write_netcdf_var(ncid,'calciu',locetra(1,1,1,icalc),2*kpke,0)
      CALL write_netcdf_var(ncid,'opal',locetra(1,1,1,iopal),2*kpke,0)
      CALL write_netcdf_var(ncid,'n2o',locetra(1,1,1,ian2o),2*kpke,0)
      CALL write_netcdf_var(ncid,'dms',locetra(1,1,1,idms),2*kpke,0)
      CALL write_netcdf_var(ncid,'fdust',locetra(1,1,1,ifdust),2*kpke,0)
      CALL write_netcdf_var(ncid,'iron',locetra(1,1,1,iiron),2*kpke,0)
      CALL write_netcdf_var(ncid,'prefo2',locetra(1,1,1,iprefo2),2*kpke,0)
      CALL write_netcdf_var(ncid,'prefpo4',locetra(1,1,1,iprefpo4),2*kpke,0)
      CALL write_netcdf_var(ncid,'prefalk',locetra(1,1,1,iprefalk),2*kpke,0)
      CALL write_netcdf_var(ncid,'prefdic',locetra(1,1,1,iprefdic),2*kpke,0)
      CALL write_netcdf_var(ncid,'dicsat',locetra(1,1,1,idicsat),2*kpke,0)
#ifdef cisonew
      CALL write_netcdf_var(ncid,'sco213',locetra(1,1,1,isco213),2*kpke,0)
      CALL write_netcdf_var(ncid,'sco214',locetra(1,1,1,isco214),2*kpke,0)
      CALL write_netcdf_var(ncid,'doc13',locetra(1,1,1,idoc13),2*kpke,0)
      CALL write_netcdf_var(ncid,'doc14',locetra(1,1,1,idoc14),2*kpke,0)
      CALL write_netcdf_var(ncid,'poc13',locetra(1,1,1,idet13),2*kpke,0)
      CALL write_netcdf_var(ncid,'poc14',locetra(1,1,1,idet14),2*kpke,0)
      CALL write_netcdf_var(ncid,'phyto13',locetra(1,1,1,iphy13),2*kpke,0)
      CALL write_netcdf_var(ncid,'phyto14',locetra(1,1,1,iphy14),2*kpke,0)
      CALL write_netcdf_var(ncid,'grazer13',locetra(1,1,1,izoo13),2*kpke,0)
      CALL write_netcdf_var(ncid,'grazer14',locetra(1,1,1,izoo14),2*kpke,0)
      CALL write_netcdf_var(ncid,'calciu13',locetra(1,1,1,icalc13),2*kpke,0)
      CALL write_netcdf_var(ncid,'calciu14',locetra(1,1,1,icalc14),2*kpke,0)
#endif
#ifdef AGG
      CALL write_netcdf_var(ncid,'snos',locetra(1,1,1,inos),2*kpke,0)
      CALL write_netcdf_var(ncid,'adust',locetra(1,1,1,iadust),2*kpke,0)
#endif /*AGG*/
#ifdef CFC
      CALL write_netcdf_var(ncid,'cfc11',locetra(1,1,1,icfc11),2*kpke,0)
      CALL write_netcdf_var(ncid,'cfc12',locetra(1,1,1,icfc12),2*kpke,0)
      CALL write_netcdf_var(ncid,'sf6',locetra(1,1,1,isf6),2*kpke,0)
#endif
#ifdef natDIC
      CALL write_netcdf_var(ncid,'natsco212',locetra(1,1,1,inatsco212),2*kpke,0)
      CALL write_netcdf_var(ncid,'natalkali',locetra(1,1,1,inatalkali),2*kpke,0)
      CALL write_netcdf_var(ncid,'natcalciu',locetra(1,1,1,inatcalc),2*kpke,0)
#endif
!
! Write restart data : diagtnostic ocean fields
!
      CALL write_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
      CALL write_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
      CALL write_netcdf_var(ncid,'co2star',co2star(1,1,1),kpke,0)
      CALL write_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0)
#ifdef natDIC
      CALL write_netcdf_var(ncid,'nathi',nathi(1,1,1),kpke,0)
#endif
!
! Write restart data : sediment variables.
!
#ifndef sedbypass
      CALL write_netcdf_var(ncid,'ssso12',sedlay2(1,1,1,issso12),2*ks,0)
      CALL write_netcdf_var(ncid,'sssc12',sedlay2(1,1,1,isssc12),2*ks,0)
      CALL write_netcdf_var(ncid,'ssssil',sedlay2(1,1,1,issssil),2*ks,0)
      CALL write_netcdf_var(ncid,'ssster',sedlay2(1,1,1,issster),2*ks,0)
      CALL write_netcdf_var(ncid,'bur_o12',burial2(1,1,1,issso12),2,0)
      CALL write_netcdf_var(ncid,'bur_c12',burial2(1,1,1,isssc12),2,0)
      CALL write_netcdf_var(ncid,'bur_sil',burial2(1,1,1,issssil),2,0)
      CALL write_netcdf_var(ncid,'bur_clay',burial2(1,1,1,issster),2,0)
      CALL write_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0)
      CALL write_netcdf_var(ncid,'powaic',powtra2(1,1,1,ipowaic),2*ks,0)
      CALL write_netcdf_var(ncid,'powaal',powtra2(1,1,1,ipowaal),2*ks,0)
      CALL write_netcdf_var(ncid,'powaph',powtra2(1,1,1,ipowaph),2*ks,0)
      CALL write_netcdf_var(ncid,'powaox',powtra2(1,1,1,ipowaox),2*ks,0)
      CALL write_netcdf_var(ncid,'pown2',powtra2(1,1,1,ipown2),2*ks,0)
      CALL write_netcdf_var(ncid,'powno3',powtra2(1,1,1,ipowno3),2*ks,0)
      CALL write_netcdf_var(ncid,'powasi',powtra2(1,1,1,ipowasi),2*ks,0)
#ifdef cisonew
      CALL write_netcdf_var(ncid,'ssso13',sedlay2(1,1,1,issso13),2*ks,0)
      CALL write_netcdf_var(ncid,'ssso14',sedlay2(1,1,1,issso14),2*ks,0)
      CALL write_netcdf_var(ncid,'sssc13',sedlay2(1,1,1,isssc13),2*ks,0)
      CALL write_netcdf_var(ncid,'sssc14',sedlay2(1,1,1,isssc14),2*ks,0)
      CALL write_netcdf_var(ncid,'powc13',powtra2(1,1,1,ipowc13),2*ks,0)
      CALL write_netcdf_var(ncid,'powc14',powtra2(1,1,1,ipowc14),2*ks,0)
#endif
#endif
!
! Write restart data: atmosphere.
!
#if defined(BOXATM)
      CALL write_netcdf_var(ncid,'atmco2',atm2(1,1,1,iatmco2),2,0)
      CALL write_netcdf_var(ncid,'atmo2',atm2(1,1,1,iatmo2),2,0)
      CALL write_netcdf_var(ncid,'atmn2',atm2(1,1,1,iatmn2),2,0)
#ifdef cisonew
      CALL write_netcdf_var(ncid,'atmc13',atm2(1,1,1,iatmc13),2,0)
      CALL write_netcdf_var(ncid,'atmc14',atm2(1,1,1,iatmc14),2,0)
#endif
#ifdef natDIC
      CALL write_netcdf_var(ncid,'atmnco2',atm2(1,1,1,iatmnco2),2,0)
#endif
#endif


      IF(mnproc==1 .AND. IOTYPE==0) THEN
        ncstat = NF90_CLOSE(ncid)
        IF ( ncstat .NE. NF90_NOERR ) THEN
          call xchalt('(AUFW: netCDF200)')
                 stop '(AUFW: netCDF200)'
        ENDIF
      ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
        ncstat = NFMPI_CLOSE(ncid)
        IF ( ncstat .NE. NF_NOERR ) THEN
          call xchalt('(AUFW: PnetCDF200)')
                 stop '(AUFW: PnetCDF200)'
        ENDIF
#endif
      ENDIF

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'End of AUFW_BGC'
      WRITE(io_stdo_bgc,*) '***************'
      ENDIF

      RETURN
      END
