! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, A. Moree
!                     M. Bentsen, P.-G. Chiu
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

module mo_aufw_bgc

  implicit none
  private

  public :: aufw_bgc

contains

  subroutine aufw_bgc(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,trc,                                        &
       &              kplyear,kplmon,kplday,kpldtoce,omask,rstfnm)

    !***********************************************************************************************
    ! Write marine bgc restart data.
    !
    ! Write restart data for continuation of interrupted integration.
    ! The bgc data are written to an extra file, other than the ocean data.
    ! The time stamp of the bgc restart file (idate) is taken from the
    ! ocean time stamp through the SBR parameter list. The only time
    ! control variable proper to the bgc is the time step number (idate(5)).
    ! It can differ from that of the ocean (idate(4)) by the difference
    ! of the offsets of restart files.
    !
    ! Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    ! Modified
    ! S.Legutke,        *MPI-MaD, HH*    10.04.01
    ! - extra SBR for writing bgc data to the restart file.
    ! S.Legutke,        *MPI-MaD, HH*    15.08.01
    ! - netCDF version (cond.comp. PNETCDF)
    ! - chemcm is multiplied with layer-dependent constant in order
    !   to be displayable by ncview. It is not read in AUFR_BGC!
    ! J.Schwinger,      *GFI, Bergen*     2013-10-21
    ! - tracer field is passed from ocean model for writing now
    ! - removed writing of chemcm and ak* fields
    ! - code cleanup, removed preprocessor option "PNETCDF"
    ! J.Schwinger,      *GFI, Bergen*     2014-05-21
    ! - adapted code for writing of two time level tracer and
    !   sediment fields
    ! A.Moree,          *GFI, Bergen*   2018-04-12
    ! - new version of carbon isotope code
    ! J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
    ! - added preformed and saturated DIC tracers
    ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    ! - removed satn2o which is not needed to restart the model
    ! - added sediment bypass preprocessor option
    ! J.Schwinger,      *Uni Research, Bergen*   2018-08-23
    ! - added writing of atmosphere field for BOXATM
    ! M. Bentsen,       *NORCE, Bergen*          2020-05-03
    ! - changed ocean model from MICOM to BLOM
    !***********************************************************************************************

    use netcdf,         only: nf90_64bit_offset,nf90_global,nf90_noerr,nf90_nofill,nf90_def_dim,   &
                              nf90_enddef,nf90_close,nf90_create,nf90_put_att,nf90_set_fill
    use mod_xc,         only: nbdy,itdm,jtdm,mnproc,iqr,jqr,xchalt
    use mod_dia,        only: iotype
    use mo_carbch,      only: co2star,co3,hi,satoxy,nathi
    use mo_control_bgc, only: io_stdo_bgc,ldtbgc,rmasks,rmasko,use_cisonew,use_AGG,use_BOXATM,     &
                              use_BROMO,use_CFC,use_natDIC,use_sedbypass,use_extNcycle,            &
                              use_pref_tracers,use_shelfsea_res_time,use_sediment_quality,         &
                              use_river2omip,use_DOMclasses
    use mo_sedmnt,      only: sedhpl
    use mo_intfcblom,   only: sedlay2,powtra2,burial2,atm2,prorca_mavg2
    use mo_param1_bgc,  only: ialkali, ian2o,iano3,icalc,idet,idicsat,idms,idoc,ifdust,igasnit,    &
                              iiron,iopal,ioxygen,iphosph,iphy,iprefalk,iprefdic,iprefo2,iprefpo4, &
                              isco212,isilica,izoo,ks,nocetra,iadust, inos,iatmco2,iatmn2,iatmo2,  &
                              ibromo,icfc11,icfc12,isf6,icalc13,icalc14,idet13,idet14,idoc13,      &
                              idoc14,iphy13,iphy14,isco213,isco214,izoo13,izoo14,issso13,issso14,  &
                              isssc13,isssc14,ipowc13,ipowc14,iatmnco2,iatmc13,iatmc14,inatalkali, &
                              inatcalc,inatsco212,ipowaal,ipowaic,ipowaox,ipowaph,ipowasi,ipown2,  &
                              ipowno3,isssc12,issso12,issssil,issster,iprefsilica,ianh4,iano2,     &
                              ipownh4,ipown2o,ipowno2,ishelfage,issso12_age,itdoc_lc,itdoc_hc,     &
                              itdoc_lc13,itdoc_hc13,itdoc_lc14,itdoc_hc14,idocsl,idocsr,idocr,     &
                              iprefdoc,iprefdocsl,iprefdocsr,iprefdocr
    use mo_netcdf_bgcrw,only: write_netcdf_var,netcdf_def_vardb
#ifdef PNETCDF
    use mod_xc,         only: mpicomm
#endif

    ! Arguments
    integer,          intent(in) :: kpie             ! 1st dimension of model grid.
    integer,          intent(in) :: kpje             ! 2nd dimension of model grid.
    integer,          intent(in) :: kpke             ! 3rd (vertical) dimension of model grid.
    integer,          intent(in) :: ntr              ! number of tracers in tracer field
    integer,          intent(in) :: ntrbgc           ! number of bgc tracers in tracer field
    integer,          intent(in) :: itrbgc           ! start index for bgc tracers in tracer field
    integer,          intent(in) :: kplyear          ! year  in ocean restart date
    integer,          intent(in) :: kplmon           ! month in ocean restart date
    integer,          intent(in) :: kplday           ! day   in ocean restart date
    integer,          intent(in) :: kpldtoce         ! step  in ocean restart date
    real,             intent(in) :: omask(kpie,kpje) ! land/ocean mask
    character(len=*), intent(in) :: rstfnm           ! restart file name-informations
    ! initial/restart tracer field to be passed to the ocean model [mol/kg]
    real,             intent(in) :: trc(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,2*kpke,ntr)

    ! Local variables
    integer             :: i,j
    real                :: locetra(kpie,kpje,2*kpke,nocetra)
    integer             :: errstat

    ! Variables for netcdf
    integer             :: ncid,ncvarid,ncstat,ncoldmod,ncdimst(4)
    integer             :: nclatid,nclonid,nclevid,nclev2id,ncksid,ncks2id,nctlvl2id
    integer             :: idate(5),ierr,testio
    real                :: rmissing
    character(len=3)    :: stripestr
    character(len=9)    :: stripestr2
#ifdef PNETCDF
#   include <pnetcdf.inc>
#   include <mpif.h>
    integer(kind=MPI_OFFSET_KIND) :: clen
    integer*4 ,save               :: info=MPI_INFO_NULL
#endif

    ! pass tracer fields in from ocean model, note that both timelevels
    ! are passed into the local array locetra; No unit conversion here,
    ! tracers in the restart file are written in mol/kg
    !--------------------------------------------------------------------
    locetra(:,:,:,:) = trc(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)


    idate(1) = kplyear
    idate(2) = kplmon
    idate(3) = kplday
    idate(4) = kpldtoce
    idate(5) = ldtbgc
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ' '
      write(io_stdo_bgc,*) 'Writing restart file at date : YY=',idate(1), &
           &                                             ' MM=',idate(2), &
           &                                            ' day=',idate(3)
      write(io_stdo_bgc,*) 'Ocean model step number is ',       idate(4)
      write(io_stdo_bgc,*) 'Bgc   model step number is ',       idate(5)
    endif

    testio=0
    rmissing = rmasko
    !
    ! Open netCDF data file
    !
    if(mnproc==1 .and. IOTYPE==0) then
      write(io_stdo_bgc,*) 'BGC RESTART   ',rstfnm
      ncstat = NF90_CREATE(rstfnm,NF90_64BIT_OFFSET,ncid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF1)')
        stop        '(AUFW_BGC: Problem with netCDF1)'
      endif
    else if (IOTYPE==1) then
#ifdef PNETCDF
      testio=1
      if(mnproc==1) write(io_stdo_bgc,*) 'BGC RESTART   ',rstfnm
      write(stripestr,('(i3)')) 16
      write(stripestr2,('(i9)')) 1024*1024
      call mpi_info_create(info,ierr)
      call mpi_info_set(info,'romio_ds_read','disable',ierr)
      call mpi_info_set(info,'romio_ds_write','disable',ierr)
      call mpi_info_set(info,"striping_factor",stripestr,ierr)
      call mpi_info_set(info,"striping_unit",stripestr2,ierr)
      ncstat = NFMPI_CREATE(mpicomm,rstfnm, &
           &       IOR(nf_clobber,nf_64bit_offset),info,ncid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF1)')
        stop        '(AUFW_BGC: Problem with PnetCDF1)'
      endif
#endif
      if(testio .eq. 0) then
        call xchalt('(AUFW_BGC: Problem with namelist iotype)')
        stop        '(AUFW_BGC: Problem with namelist iotype)'
      endif

    endif
    !
    ! Define dimension
    ! ----------------------------------------------------------------------
    !
    if(mnproc==1 .and. IOTYPE==0) then

      ncstat = NF90_DEF_DIM(ncid, 'lon', itdm, nclonid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF2)')
        stop        '(AUFW_BGC: Problem with netCDF2)'
      endif

      ncstat = NF90_DEF_DIM(ncid, 'lat', jtdm, nclatid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF3)')
        stop        '(AUFW_BGC: Problem with netCDF3)'
      endif

      ncstat = NF90_DEF_DIM(ncid, 'depth', kpke, nclevid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF4)')
        stop        '(AUFW_BGC: Problem with netCDF4)'
      endif

      ncstat = NF90_DEF_DIM(ncid, 'depth2', 2*kpke, nclev2id)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF5)')
        stop        '(AUFW_BGC: Problem with netCDF5)'
      endif

      ncstat = NF90_DEF_DIM(ncid, 'nks', ks, ncksid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF6)')
        stop        '(AUFW_BGC: Problem with netCDF6)'
      endif

      ncstat = NF90_DEF_DIM(ncid, 'nks2', 2*ks, ncks2id)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF7)')
        stop        '(AUFW_BGC: Problem with netCDF7)'
      endif

      ncstat = NF90_DEF_DIM(ncid, 'tlvl2', 2, nctlvl2id)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF8)')
        stop        '(AUFW_BGC: Problem with netCDF8)'
      endif

    else if (IOTYPE==1) then
#ifdef PNETCDF
      clen=itdm
      ncstat = NFMPI_DEF_DIM(ncid, 'lon', clen, nclonid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF2)')
        stop        '(AUFW_BGC: Problem with PnetCDF2)'
      endif

      clen=jtdm
      ncstat = NFMPI_DEF_DIM(ncid, 'lat', clen, nclatid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF3)')
        stop        '(AUFW_BGC: Problem with PnetCDF3)'
      endif

      clen=kpke
      ncstat = NFMPI_DEF_DIM(ncid, 'depth', clen, nclevid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF4)')
        stop        '(AUFW_BGC: Problem with PnetCDF4)'
      endif

      clen=2*kpke
      ncstat = NFMPI_DEF_DIM(ncid, 'depth2', clen, nclev2id)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF5)')
        stop        '(AUFW_BGC: Problem with PnetCDF5)'
      endif

      clen=ks
      ncstat = NFMPI_DEF_DIM(ncid, 'nks', clen, ncksid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF6)')
        stop        '(AUFW_BGC: Problem with PnetCDF6)'
      endif

      clen=2*ks
      ncstat = NFMPI_DEF_DIM(ncid, 'nks2', clen, ncks2id)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF7)')
        stop        '(AUFW_BGC: Problem with PnetCDF7)'
      endif

      clen=2
      ncstat = NFMPI_DEF_DIM(ncid, 'tlvl2', clen, nctlvl2id)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF8)')
        stop        '(AUFW_BGC: Problem with PnetCDF8)'
      endif
#endif
    endif  !mnproc==1 .and. IOTYPE==0

    !
    ! Define global attributes
    ! ----------------------------------------------------------------------
    !
    if (mnproc==1 .and. IOTYPE==0) then

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'title'          &
           &, 'Restart data for marine bgc modules')
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF9)')
        stop        '(AUFW_BGC: Problem with netCDF9)'
      endif

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'history'        &
           &, 'Restart data for marine bgc modules')
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF10)')
        stop        '(AUFW_BGC: Problem with netCDF10)'
      endif

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'conventions'    &
           &,'COARDS')
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF11)')
        stop        '(AUFW_BGC: Problem with netCDF11)'
      endif

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'source'         &
           &, 'Marine bgc model output HOPC68/grob')
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF12)')
        stop        '(AUFW_BGC: Problem with netCDF12)'
      endif

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'date', idate)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF13)')
        stop        '(AUFW_BGC: Problem with netCDF13)'
      endif

    else if (IOTYPE==1) then
#ifdef PNETCDF
      clen=len('Restart data for marine bgc modules')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'      &
           &, clen,'Restart data for marine bgc modules')
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF9)')
        stop        '(AUFW_BGC: Problem with PnetCDF9)'
      endif

      clen=len('Restart data for marine bgc modules')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'history'    &
           &, clen,'Restart data for marine bgc modules')
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF10)')
        stop        '(AUFW_BGC: Problem with PnetCDF10)'
      endif

      clen=6
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'conventions'&
           &,clen, 'COARDS')
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF11)')
        stop        '(AUFW_BGC: Problem with PnetCDF11)'
      endif

      clen=len('Marine bgc model output HOPC68/grob')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'     &
           &,clen, 'Marine bgc model output HOPC68/grob')
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF12)')
        stop        '(AUFW_BGC: Problem with PnetCDF12)'
      endif

      clen=5
      ncstat = NFMPI_PUT_ATT_INT(ncid, NF_GLOBAL, 'date',      &
           &                           nf_int, clen, idate)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF13)')
        stop        '(AUFW_BGC: Problem with netCDF13)'

      endif
#endif
    endif ! IOTYPE == 1
    !
    ! Define variables : advected ocean tracer
    ! ----------------------------------------------------------------------
    !
    if((mnproc==1 .and. IOTYPE==0) .OR. IOTYPE==1) then
      ncdimst(1) = nclonid
      ncdimst(2) = nclatid
      ncdimst(3) = nclev2id
      ncdimst(4) = 0
    endif

    call NETCDF_DEF_VARDB(ncid,6,'sco212',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',13, 'Dissolved CO2',rmissing,10,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'alkali',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',10,'Alkalinity',rmissing,11,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'phosph',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',19,'Dissolved phosphate',rmissing,12,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'oxygen',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',16,'Dissolved oxygen',rmissing,13,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'gasnit',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',21,'Gaseous nitrogen (N2)',rmissing,14,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,4,'ano3',3,ncdimst,ncvarid,                                         &
         &    6,'mol/kg',17,'Dissolved nitrate',rmissing,15,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'silica',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',22,'Silicid acid (Si(OH)4)',rmissing,16,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,3,'doc',3,ncdimst,ncvarid,                                          &
         &    6,'mol/kg',24,'Dissolved organic carbon',rmissing,17,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,3,'poc',3,ncdimst,ncvarid,                                          &
         &    6,'mol/kg',26,'Particulate organic carbon',rmissing,18,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,5,'phyto',3,ncdimst,ncvarid,                                        &
         &    7,'molP/kg',27,'Phytoplankton concentration',rmissing,19,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'grazer',3,ncdimst,ncvarid,                                       &
         &    7,'molP/kg',25,'Zooplankton concentration',rmissing,20,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'calciu',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',17,'Calcium carbonate',rmissing,21,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,4,'opal',3,ncdimst,ncvarid,                                         &
         &    6,'mol/kg',15,'Biogenic silica',rmissing,22,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,3,'n2o',3,ncdimst,ncvarid,                                          &
         &    6,'mol/kg',12,'laughing gas',rmissing,23,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,3,'dms',3,ncdimst,ncvarid,                                          &
         &    6,'mol/kg',15 ,'DiMethylSulfide',rmissing,24,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,5,'fdust',3,ncdimst,ncvarid,                                        &
         &    5,'kg/kg',19,'Non-aggregated dust',rmissing,25,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,4,'iron',3,ncdimst,ncvarid,                                         &
         &    6,'mol/kg',14,'Dissolved iron',rmissing,26,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'dicsat',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',13,'Saturated dic',rmissing,31,io_stdo_bgc)

    if (use_pref_tracers) then
      call NETCDF_DEF_VARDB(ncid,6,'prefo2',3,ncdimst,ncvarid,                                     &
           &    6,'mol/kg',16,'Preformed oxygen',rmissing,27,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'prefpo4',3,ncdimst,ncvarid,                                    &
           &    6,'mol/kg',19,'Preformed phosphate',rmissing,28,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,10,'prefsilica',3,ncdimst,ncvarid,                                &
           &    6,'mol/kg',16,'Preformed silica',rmissing,28,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'prefalk',3,ncdimst,ncvarid,                                    &
           &    6,'mol/kg',20,'Preformed alkalinity',rmissing,29,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'prefdic',3,ncdimst,ncvarid,                                    &
           &    6,'mol/kg',13,'Preformed dic',rmissing,30,io_stdo_bgc)
    endif
    if (use_cisonew) then
      call NETCDF_DEF_VARDB(ncid,6,'sco213',3,ncdimst,ncvarid,                                     &
           &    6,'mol/kg',15, 'Dissolved CO213',rmissing,32,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'sco214',3,ncdimst,ncvarid,                                     &
           &    6,'mol/kg',15, 'Dissolved CO214',rmissing,33,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'doc13',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',24,'Dissolved organic carb13',rmissing,34,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'doc14',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',24,'Dissolved organic carb14',rmissing,35,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'poc13',3,ncdimst,ncvarid,                                      &
           &    7,'molC/kg',28,'Particulate organic carbon13',rmissing,36,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'poc14',3,ncdimst,ncvarid,                                      &
           &    7,'molC/kg',28,'Particulate organic carbon14',rmissing,37,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'phyto13',3,ncdimst,ncvarid,                                    &
           &    7,'molP/kg',27,'Phytoplankton concentr. 13c',rmissing,38,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'phyto14',3,ncdimst,ncvarid,                                    &
           &    7,'molP/kg',27,'Phytoplankton concentr. 14c',rmissing,39,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,8,'grazer13',3,ncdimst,ncvarid,                                   &
           &    7,'molP/kg',25,'Zooplankton concentr. 13c',rmissing,40,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,8,'grazer14',3,ncdimst,ncvarid,                                   &
           &    7,'molP/kg',25,'Zooplankton concentr. 14c',rmissing,41,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,8,'calciu13',3,ncdimst,ncvarid,                                   &
           &    7,'molC/kg',19,'Calcium carbonate13',rmissing,42,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,8,'calciu14',3,ncdimst,ncvarid,                                   &
           &    7,'molC/kg',19,'Calcium carbonate14',rmissing,43,io_stdo_bgc)
    endif
    if (use_AGG) then
      call NETCDF_DEF_VARDB(ncid,4,'snos',3,ncdimst,ncvarid,                                       &
           &    3,'1/g',38,'marine snow aggregates per g sea water',rmissing,44,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'adust',3,ncdimst,ncvarid,                                      &
           &    4,'g/kg',15,'Aggregated dust',rmissing,45,io_stdo_bgc)
    endif
    if (use_CFC) then
      call NETCDF_DEF_VARDB(ncid,5,'cfc11',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',5,'CFC11',rmissing,47,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'cfc12',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',5,'CFC12',rmissing,48,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,3,'sf6',3,ncdimst,ncvarid,                                        &
           &    6,'mol/kg',4,'SF-6',rmissing,49,io_stdo_bgc)
    endif
    if (use_natDIC) then
      call NETCDF_DEF_VARDB(ncid,9,'natsco212',3,ncdimst,ncvarid,                                  &
           &   6,'mol/kg',21, 'Natural dissolved CO2',rmissing,50,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,9,'natalkali',3,ncdimst,ncvarid,                                  &
           &    6,'mol/kg',18,'Natural alkalinity',rmissing,51,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,9,'natcalciu',3,ncdimst,ncvarid,                                  &
           &    6,'mol/kg',25,'Natural calcium carbonate',rmissing,52,io_stdo_bgc)
    endif
    if (use_BROMO) then
      call NETCDF_DEF_VARDB(ncid,5,'bromo',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',9,'Bromoform',rmissing,47,io_stdo_bgc)
    endif
    if (use_extNcycle) then
      call NETCDF_DEF_VARDB(ncid,4,'anh4',3,ncdimst,ncvarid,                                       &
           &    6,'mol/kg',18,'Dissolved ammonium',rmissing,54,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,4,'ano2',3,ncdimst,ncvarid,                                       &
           &    6,'mol/kg',17,'Dissolved nitrite',rmissing,55,io_stdo_bgc)
    endif
    if (use_shelfsea_res_time) then
      call NETCDF_DEF_VARDB(ncid,8,'shelfage',3,ncdimst,ncvarid,                                   &
           &    1,'d',25,'Shelfwater residence time',rmissing,56,io_stdo_bgc)
    endif
    if (use_river2omip) then
      call NETCDF_DEF_VARDB(ncid,7,'tdoc_lc',3,ncdimst,ncvarid,                                    &
           &    6,'mol/kg',52,'Terrestrial dissolved organic carbon (low C content)',              &
           &    rmissing,57,io_stdo_bgc)
      call NETCDF_DEF_VARDB(ncid,7,'tdoc_hc',3,ncdimst,ncvarid,                                    &
           &    6,'mol/kg',53,'Terrestrial dissolved organic carbon (high C content)',             &
           &    rmissing,58,io_stdo_bgc)
      if (use_cisonew) then
        call NETCDF_DEF_VARDB(ncid,9,'tdoc_lc13',3,ncdimst,ncvarid,                                &
             &    6,'mol/kg',54,'Terrestrial dissolved organic carbon13 (low C content)',          &
             &    rmissing,59,io_stdo_bgc)
        call NETCDF_DEF_VARDB(ncid,9,'tdoc_hc13',3,ncdimst,ncvarid,                                &
             &    6,'mol/kg',55,'Terrestrial dissolved organic carbon13 (high C content)',         &
             &    rmissing,60,io_stdo_bgc)
        call NETCDF_DEF_VARDB(ncid,9,'tdoc_lc14',3,ncdimst,ncvarid,                                &
             &    6,'mol/kg',54,'Terrestrial dissolved organic carbon14 (low C content)',          &
             &    rmissing,61,io_stdo_bgc)
        call NETCDF_DEF_VARDB(ncid,9,'tdoc_hc14',3,ncdimst,ncvarid,                                &
             &    6,'mol/kg',55,'Terrestrial dissolved organic carbon14 (high C content)',         &
             &    rmissing,62,io_stdo_bgc)
      endif
    endif
    if (use_DOMclasses) then
      call NETCDF_DEF_VARDB(ncid,5,'docsl',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',36,'Semi labile dissolved organic carbon',rmissing,56,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'docsr',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',40,'Semi refractory dissolved organic carbon',rmissing,57,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,4,'docr',3,ncdimst,ncvarid,                                       &
           &    6,'mol/kg',35,'Refractory dissolved organic carbon',rmissing,58,io_stdo_bgc)
    endif
    if (use_DOMclasses .and. use_pref_tracers) then
      call NETCDF_DEF_VARDB(ncid,7,'prefdoc',3,ncdimst,ncvarid,                                    &
        &    6,'mol/kg',20,'Preformed labile DOC',rmissing,59,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,9,'prefdocsl',3,ncdimst,ncvarid,                                  &
        &    6,'mol/kg',25,'Preformed semi-labile DOC',rmissing,60,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,9,'prefdocsr',3,ncdimst,ncvarid,                                  &
        &    6,'mol/kg',29,'Preformed semi-refractory DOC',rmissing,61,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,8,'prefdocr',3,ncdimst,ncvarid,                                   &
        &    6,'mol/kg',23,'Preformed refactory DOC',rmissing,62,io_stdo_bgc)
    endif

    !
    ! Define variables : diagnostic ocean fields
    ! ----------------------------------------------------------------------
    !
    if((mnproc==1 .and. IOTYPE==0) .OR. IOTYPE==1) then
      ncdimst(1) = nclonid
      ncdimst(2) = nclatid
      ncdimst(3) = nclevid
      ncdimst(4) = 0
    endif

    call NETCDF_DEF_VARDB(ncid,2,'hi',3,ncdimst,ncvarid,                                           &
         &    6,'mol/kg',26,'Hydrogen ion concentration',rmissing,63,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,3,'co3',3,ncdimst,ncvarid,                                          &
         &    6,'mol/kg',25,'Dissolved carbonate (CO3)',rmissing,64,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,7,'co2star',3,ncdimst,ncvarid,                                      &
         &    6,'mol/kg',20,'Dissolved CO2 (CO2*)',rmissing,65,io_stdo_bgc)

    call NETCDF_DEF_VARDB(ncid,6,'satoxy',3,ncdimst,ncvarid,                                       &
         &    6,'mol/kg',16 ,'Saturated oxygen',rmissing,66,io_stdo_bgc)

    if (use_natDIC) then
      call NETCDF_DEF_VARDB(ncid,5,'nathi',3,ncdimst,ncvarid,                                      &
           &    6,'mol/kg',34,'Natural hydrogen ion concentration',rmissing,67,io_stdo_bgc)
    endif
    !
    ! Define variables : sediment
    ! ----------------------------------------------------------------------
    !
    if (.not. use_sedbypass) then

      if((mnproc==1 .and. IOTYPE==0) .OR. IOTYPE==1) then
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = ncks2id
        ncdimst(4) = 0
      endif

      call NETCDF_DEF_VARDB(ncid,6,'ssso12',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',35,'Sediment accumulated organic carbon',rmissing,70,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'sssc12',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',38,'Sediment accumulated calcium carbonate',rmissing,71,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'ssssil',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',25,'Sediment accumulated opal',rmissing,72,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'ssster',3,ncdimst,ncvarid,                                     &
           &    7,'kg/m**3',25,'Sediment accumulated clay',rmissing,73,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'powaic',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',23,'Sediment pore water CO2',rmissing,74,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'powaal',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',30,'Sediment pore water alkalinity',rmissing,75,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'powaph',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',29,'Sediment pore water phosphate',rmissing,76,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'powaox',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',26,'Sediment pore water oxygen',rmissing,77,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'pown2',3,ncdimst,ncvarid,                                      &
           &    9,'kmol/m**3',36,'Sediment pore water gaseous nitrogen',rmissing,78,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'powno3',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',33,'Sediment pore water nitrate (NO3)',rmissing,79,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,6,'powasi',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**3',42,'Sediment pore water silicid acid (Si(OH)4)',                     &
           &    rmissing,80,io_stdo_bgc)

      if (use_cisonew) then
        call NETCDF_DEF_VARDB(ncid,6,'ssso13',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',37,'Sediment accumulated organic carbon13',rmissing,81,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'ssso14',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',37,'Sediment accumulated organic carbon14',rmissing,82,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'sssc13',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',40,'Sediment accumulated calcium carbonate13',                     &
             &    rmissing,83,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'sssc14',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',40,'Sediment accumulated calcium carbonate14',                     &
             &    rmissing,84,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'powc13',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',25,'Sediment pore water DIC13', rmissing,85,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'powc14',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',25,'Sediment pore water DIC14',rmissing,86,io_stdo_bgc)
      endif

      if (use_extNcycle) then
        call NETCDF_DEF_VARDB(ncid,6,'pownh4',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',34,'Sediment pore water ammonium (NH4)',rmissing,79,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'pown2o',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',39,'Sediment pore water nitrous oxide (N2O)',rmissing,79,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,6,'powno2',3,ncdimst,ncvarid,                                   &
             &    9,'kmol/m**3',33,'Sediment pore water nitrite (NO2)',rmissing,79,io_stdo_bgc)
      endif
      if (use_sediment_quality) then
        call NETCDF_DEF_VARDB(ncid,10,'ssso12_age',3,ncdimst,ncvarid,                              &
           &    2,'yr',39,'Sediment accumulated organic carbon age',rmissing,70,io_stdo_bgc)
        call NETCDF_DEF_VARDB(ncid,11,'prorca_mavg',3,ncdimst,ncvarid,                             &
           &    11,'kmol/m**2/d',51,'Moving average of organic carbon sedimentation flux',rmissing,&
           &    70,io_stdo_bgc)
      endif

      if((mnproc==1 .and. IOTYPE==0) .OR. IOTYPE==1) then
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = ncksid
        ncdimst(4) = 0
      endif

      call NETCDF_DEF_VARDB(ncid,6,'sedhpl',3,ncdimst,ncvarid,                                     &
           &    9,'kmol/m**2',34,'Sediment accumulated hydrogen ions',rmissing,87,io_stdo_bgc)
      !
      ! Define variables : sediment burial
      ! ----------------------------------------------------------------------
      !
      if((mnproc==1 .and. IOTYPE==0) .OR. IOTYPE==1) then
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = nctlvl2id
        ncdimst(4) = 0
      endif

      call NETCDF_DEF_VARDB(ncid,7,'bur_o12',3,ncdimst,ncvarid,                                    &
           &    9,'kmol/m**2',30,'Burial layer of organic carbon',rmissing,90,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'bur_c12',3,ncdimst,ncvarid,                                    &
           &    9,'kmol/m**2',33,'Burial layer of calcium carbonate',rmissing,91,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,7,'bur_sil',3,ncdimst,ncvarid,                                    &
           &    9,'kmol/m**2',20,'Burial layer of opal',rmissing,92,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,8,'bur_clay',3,ncdimst,ncvarid,                                   &
           &    7,'kg/m**2',20,'Burial layer of clay',rmissing,93,io_stdo_bgc)

      if (use_cisonew) then
        call NETCDF_DEF_VARDB(ncid,8,'bur_o13',3,ncdimst,ncvarid,                                  &
             &    9,'kmol/m**2',27,'Burial layer of organic 13C',rmissing,94,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,8,'bur_o14',3,ncdimst,ncvarid,                                  &
             &    9,'kmol/m**2',27,'Burial layer of organic 14C',rmissing,95,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,8,'bur_c13',3,ncdimst,ncvarid,                                  &
             &    9,'kmol/m**2',23,'Burial layer of Ca13CO3',rmissing,96,io_stdo_bgc)

        call NETCDF_DEF_VARDB(ncid,8,'bur_c14',3,ncdimst,ncvarid,                                  &
             &    9,'kmol/m**2',23,'Burial layer of Ca14CO3',rmissing,97,io_stdo_bgc)
      endif

      if (use_sediment_quality) then
        call NETCDF_DEF_VARDB(ncid,11,'bur_o12_age',3,ncdimst,ncvarid,                             &
             &    2,'yr',34,'Burial layer of organic carbon age',rmissing,97,io_stdo_bgc)
      endif

    endif ! not sedbypass
    !
    ! Define variables: atmosphere
    ! ----------------------------------------------------------------------
    !
    if (use_BOXATM) then

      if((mnproc==1 .and. IOTYPE==0) .OR. IOTYPE==1) then
        ncdimst(1) = nclonid
        ncdimst(2) = nclatid
        ncdimst(3) = nctlvl2id
        ncdimst(4) = 0
      endif

      call NETCDF_DEF_VARDB(ncid,6,'atmco2',3,ncdimst,ncvarid,                                     &
           &    3,'ppm',15,'atmospheric CO2',rmissing,101,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'atmo2',3,ncdimst,ncvarid,                                      &
           &    3,'ppm',14,'atmospheric O2',rmissing,102,io_stdo_bgc)

      call NETCDF_DEF_VARDB(ncid,5,'atmn2',3,ncdimst,ncvarid,                                      &
           &    3,'ppm',14,'atmospheric N2',rmissing,103,io_stdo_bgc)

      if (use_cisonew) then
        call NETCDF_DEF_VARDB(ncid,6,'atmc13',3,ncdimst,ncvarid,                                   &
             &    3,'ppm',17,'atmospheric 13CO2',rmissing,104,io_stdo_bgc)
        call NETCDF_DEF_VARDB(ncid,6,'atmc14',3,ncdimst,ncvarid,                                   &
             &    3,'ppm',17,'atmospheric 14CO2',rmissing,105,io_stdo_bgc)
      endif
      if (use_natDIC) then
        call NETCDF_DEF_VARDB(ncid,7,'atmnco2',3,ncdimst,ncvarid,                                  &
             &    3,'ppm',23,'natural atmospheric CO2',rmissing,106,io_stdo_bgc)
      endif
    endif ! if (use_BOXATM)

    if (mnproc==1 .and. IOTYPE==0) then

      ncstat = NF90_ENDDEF(ncid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF00)')
        stop        '(AUFW_BGC: Problem with netCDF00)'
      endif
      !
      ! Set fill mode
      ! ----------------------------------------------------------------------
      !
      ncstat = NF90_SET_FILL(ncid,NF90_NOFILL, ncoldmod)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with netCDF97)')
        stop        '(AUFW_BGC: Problem with netCDF97)'
      endif

    else if (IOTYPE==1) then

#ifdef PNETCDF
      ncstat = NFMPI_ENDDEF(ncid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: Problem with PnetCDF00)')
        stop        '(AUFW_BGC: Problem with PnetCDF00)'
      endif
#endif

    endif
    !
    ! Write restart data : ocean aquateous tracer
    !--------------------------------------------------------------------
    !
    call write_netcdf_var(ncid,'sco212',locetra(1,1,1,isco212),2*kpke,0)
    call write_netcdf_var(ncid,'alkali',locetra(1,1,1,ialkali),2*kpke,0)
    call write_netcdf_var(ncid,'phosph',locetra(1,1,1,iphosph),2*kpke,0)
    call write_netcdf_var(ncid,'oxygen',locetra(1,1,1,ioxygen),2*kpke,0)
    call write_netcdf_var(ncid,'gasnit',locetra(1,1,1,igasnit),2*kpke,0)
    call write_netcdf_var(ncid,'ano3',locetra(1,1,1,iano3),2*kpke,0)
    call write_netcdf_var(ncid,'silica',locetra(1,1,1,isilica),2*kpke,0)
    call write_netcdf_var(ncid,'doc',locetra(1,1,1,idoc),2*kpke,0)
    call write_netcdf_var(ncid,'poc',locetra(1,1,1,idet),2*kpke,0)
    call write_netcdf_var(ncid,'phyto',locetra(1,1,1,iphy),2*kpke,0)
    call write_netcdf_var(ncid,'grazer',locetra(1,1,1,izoo),2*kpke,0)
    call write_netcdf_var(ncid,'calciu',locetra(1,1,1,icalc),2*kpke,0)
    call write_netcdf_var(ncid,'opal',locetra(1,1,1,iopal),2*kpke,0)
    call write_netcdf_var(ncid,'n2o',locetra(1,1,1,ian2o),2*kpke,0)
    call write_netcdf_var(ncid,'dms',locetra(1,1,1,idms),2*kpke,0)
    call write_netcdf_var(ncid,'fdust',locetra(1,1,1,ifdust),2*kpke,0)
    call write_netcdf_var(ncid,'iron',locetra(1,1,1,iiron),2*kpke,0)
    call write_netcdf_var(ncid,'dicsat',locetra(1,1,1,idicsat),2*kpke,0)
    if (use_pref_tracers) then
      call write_netcdf_var(ncid,'prefo2',locetra(1,1,1,iprefo2),2*kpke,0)
      call write_netcdf_var(ncid,'prefpo4',locetra(1,1,1,iprefpo4),2*kpke,0)
      call write_netcdf_var(ncid,'prefsilica',locetra(1,1,1,iprefsilica),2*kpke,0)
      call write_netcdf_var(ncid,'prefalk',locetra(1,1,1,iprefalk),2*kpke,0)
      call write_netcdf_var(ncid,'prefdic',locetra(1,1,1,iprefdic),2*kpke,0)
    endif
    if (use_shelfsea_res_time) then
      call write_netcdf_var(ncid,'shelfage',locetra(1,1,1,ishelfage),2*kpke,0)
    endif
    if (use_cisonew) then
      call write_netcdf_var(ncid,'sco213'   ,locetra(1,1,1,isco213) ,2*kpke,0)
      call write_netcdf_var(ncid,'sco214'   ,locetra(1,1,1,isco214) ,2*kpke,0)
      call write_netcdf_var(ncid,'doc13'    ,locetra(1,1,1,idoc13)  ,2*kpke,0)
      call write_netcdf_var(ncid,'doc14'    ,locetra(1,1,1,idoc14)  ,2*kpke,0)
      call write_netcdf_var(ncid,'poc13'    ,locetra(1,1,1,idet13)  ,2*kpke,0)
      call write_netcdf_var(ncid,'poc14'    ,locetra(1,1,1,idet14)  ,2*kpke,0)
      call write_netcdf_var(ncid,'phyto13'  ,locetra(1,1,1,iphy13)  ,2*kpke,0)
      call write_netcdf_var(ncid,'phyto14'  ,locetra(1,1,1,iphy14)  ,2*kpke,0)
      call write_netcdf_var(ncid,'grazer13' ,locetra(1,1,1,izoo13)  ,2*kpke,0)
      call write_netcdf_var(ncid,'grazer14' ,locetra(1,1,1,izoo14)  ,2*kpke,0)
      call write_netcdf_var(ncid,'calciu13' ,locetra(1,1,1,icalc13) ,2*kpke,0)
      call write_netcdf_var(ncid,'calciu14' ,locetra(1,1,1,icalc14) ,2*kpke,0)
    endif
    if (use_AGG) then
      call write_netcdf_var(ncid,'snos',locetra(1,1,1,inos),2*kpke,0)
      call write_netcdf_var(ncid,'adust',locetra(1,1,1,iadust),2*kpke,0)
    endif
    if (use_CFC) then
      call write_netcdf_var(ncid,'cfc11',locetra(1,1,1,icfc11),2*kpke,0)
      call write_netcdf_var(ncid,'cfc12',locetra(1,1,1,icfc12),2*kpke,0)
      call write_netcdf_var(ncid,'sf6',locetra(1,1,1,isf6),2*kpke,0)
    endif
    if (use_natDIC) then
      call write_netcdf_var(ncid,'natsco212',locetra(1,1,1,inatsco212),2*kpke,0)
      call write_netcdf_var(ncid,'natalkali',locetra(1,1,1,inatalkali),2*kpke,0)
      call write_netcdf_var(ncid,'natcalciu',locetra(1,1,1,inatcalc),2*kpke,0)
    endif
    if (use_BROMO) then
      call write_netcdf_var(ncid,'bromo',locetra(1,1,1,ibromo),2*kpke,0)
    endif
    if (use_extNcycle) then
      call write_netcdf_var(ncid,'anh4',locetra(1,1,1,ianh4),2*kpke,0)
      call write_netcdf_var(ncid,'ano2',locetra(1,1,1,iano2),2*kpke,0)
    endif
    if (use_DOMclasses) then
      call write_netcdf_var(ncid,'docsl',locetra(1,1,1,idocsl),2*kpke,0)
      call write_netcdf_var(ncid,'docsr',locetra(1,1,1,idocsr),2*kpke,0)
      call write_netcdf_var(ncid,'docr' ,locetra(1,1,1,idocr),2*kpke,0)
    endif
    if (use_DOMclasses .and. use_pref_tracers) then
      call write_netcdf_var(ncid,'prefdoc',locetra(1,1,1,iprefdoc),2*kpke,0)
      call write_netcdf_var(ncid,'prefdocsl',locetra(1,1,1,iprefdocsl),2*kpke,0)
      call write_netcdf_var(ncid,'prefdocsr',locetra(1,1,1,iprefdocsr),2*kpke,0)
      call write_netcdf_var(ncid,'prefdocr',locetra(1,1,1,iprefdocr),2*kpke,0)
    endif

    !
    ! Write restart data : diagnostic ocean fields
    !
    call write_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
    call write_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
    call write_netcdf_var(ncid,'co2star',co2star(1,1,1),kpke,0)
    call write_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0)
    if (use_natDIC) then
      call write_netcdf_var(ncid,'nathi',nathi(1,1,1),kpke,0)
    endif
    !
    ! Write restart data : sediment variables.
    !
    if (.not. use_sedbypass) then
      call write_netcdf_var(ncid,'ssso12',sedlay2(1,1,1,issso12),2*ks,0)
      call write_netcdf_var(ncid,'sssc12',sedlay2(1,1,1,isssc12),2*ks,0)
      call write_netcdf_var(ncid,'ssssil',sedlay2(1,1,1,issssil),2*ks,0)
      call write_netcdf_var(ncid,'ssster',sedlay2(1,1,1,issster),2*ks,0)
      call write_netcdf_var(ncid,'bur_o12',burial2(1,1,1,issso12),2,0)
      call write_netcdf_var(ncid,'bur_c12',burial2(1,1,1,isssc12),2,0)
      call write_netcdf_var(ncid,'bur_sil',burial2(1,1,1,issssil),2,0)
      call write_netcdf_var(ncid,'bur_clay',burial2(1,1,1,issster),2,0)
      call write_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0)
      call write_netcdf_var(ncid,'powaic',powtra2(1,1,1,ipowaic),2*ks,0)
      call write_netcdf_var(ncid,'powaal',powtra2(1,1,1,ipowaal),2*ks,0)
      call write_netcdf_var(ncid,'powaph',powtra2(1,1,1,ipowaph),2*ks,0)
      call write_netcdf_var(ncid,'powaox',powtra2(1,1,1,ipowaox),2*ks,0)
      call write_netcdf_var(ncid,'pown2',powtra2(1,1,1,ipown2),2*ks,0)
      call write_netcdf_var(ncid,'powno3',powtra2(1,1,1,ipowno3),2*ks,0)
      call write_netcdf_var(ncid,'powasi',powtra2(1,1,1,ipowasi),2*ks,0)
      if (use_cisonew) then
        call write_netcdf_var(ncid,'ssso13',sedlay2(1,1,1,issso13),2*ks,0)
        call write_netcdf_var(ncid,'ssso14',sedlay2(1,1,1,issso14),2*ks,0)
        call write_netcdf_var(ncid,'sssc13',sedlay2(1,1,1,isssc13),2*ks,0)
        call write_netcdf_var(ncid,'sssc14',sedlay2(1,1,1,isssc14),2*ks,0)
        call write_netcdf_var(ncid,'bur_o13',burial2(1,1,1,issso13),2,0)
        call write_netcdf_var(ncid,'bur_o14',burial2(1,1,1,issso14),2,0)
        call write_netcdf_var(ncid,'bur_c13',burial2(1,1,1,isssc13),2,0)
        call write_netcdf_var(ncid,'bur_c14',burial2(1,1,1,isssc14),2,0)
        call write_netcdf_var(ncid,'powc13',powtra2(1,1,1,ipowc13),2*ks,0)
        call write_netcdf_var(ncid,'powc14',powtra2(1,1,1,ipowc14),2*ks,0)
      endif
      if (use_extNcycle) then
        call write_netcdf_var(ncid,'pownh4',powtra2(1,1,1,ipownh4),2*ks,0)
        call write_netcdf_var(ncid,'pown2o',powtra2(1,1,1,ipown2o),2*ks,0)
        call write_netcdf_var(ncid,'powno2',powtra2(1,1,1,ipowno2),2*ks,0)
      endif
      if (use_sediment_quality) then
        call write_netcdf_var(ncid,'ssso12_age',sedlay2(1,1,1,issso12_age),2*ks,0)
        call write_netcdf_var(ncid,'bur_o12_age',burial2(1,1,1,issso12_age),2,0)
        call write_netcdf_var(ncid,'prorca_mavg',prorca_mavg2(1,1,1),2,0)
      endif
    endif
    !
    ! Write restart data: atmosphere.
    !
    if (use_BOXATM) then
      call write_netcdf_var(ncid,'atmco2',atm2(1,1,1,iatmco2),2,0)
      call write_netcdf_var(ncid,'atmo2',atm2(1,1,1,iatmo2),2,0)
      call write_netcdf_var(ncid,'atmn2',atm2(1,1,1,iatmn2),2,0)
      if (use_cisonew) then
        call write_netcdf_var(ncid,'atmc13',atm2(1,1,1,iatmc13),2,0)
        call write_netcdf_var(ncid,'atmc14',atm2(1,1,1,iatmc14),2,0)
      endif
      if (use_natDIC) then
        call write_netcdf_var(ncid,'atmnco2',atm2(1,1,1,iatmnco2),2,0)
      endif
    endif

    if(mnproc==1 .and. IOTYPE==0) then
      ncstat = NF90_CLOSE(ncid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFW_BGC: netCDF200)')
        stop        '(AUFW_BGC: netCDF200)'
      endif
    else if(IOTYPE==1) then
#ifdef PNETCDF
      ncstat = NFMPI_CLOSE(ncid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFW_BGC: PnetCDF200)')
        stop        '(AUFW_BGC: PnetCDF200)'
      endif
#endif
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 'End of AUFW_BGC'
      write(io_stdo_bgc,*) '***************'
    endif

  end subroutine aufw_bgc

end module mo_aufw_bgc
