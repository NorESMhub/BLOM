! Copyright (C) 2002  Ernst Maier-Reimer, S. Legutke, P. Wetzel
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

module mo_aufr_bgc

  implicit none
  private

  public :: aufr_bgc

CONTAINS

  subroutine aufr_bgc(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,trc,kplyear,kplmon,kplday,omask,rstfnm)

    !***********************************************************************************************
    ! Reads marine bgc restart data.
    !
    ! Read restart data to continue an interrupted integration.
    ! The bgc data are read from an extra file, other than the ocean data.
    ! The time stamp of the bgc restart file (idate) is specified from the
    ! ocean time stamp through the SBR parameter list of AUFW_BGC. The only
    ! time control variable proper to the bgc is the time step number
    ! (idate(5)). It can differ from that of the ocean (idate(4)) by the
    ! difference of the offsets of restart files.
    !
    ! Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    ! Modified:
    ! S.Legutke,        *MPI-MaD, HH*     10.04.01
    ! - extra SBR for reading bgc data from the restart file.
    ! S.Legutke,        *MPI-MaD, HH*     15.08.01
    ! - netCDF version (with cond.comp. PNETCDF)
    ! - no use of chemc values from netCDF restart
    ! Patrick Wetzel,   *MPI-Met, HH*     16.04.02
    ! - read chemcm(i,j,7,12) from netCDF restart
    ! J.Schwinger,      *GFI, Bergen*     2013-10-21
    ! - removed reading of chemcm and ak* fields
    ! - code cleanup, remoded preprocessor option "PNETCDF"
    !   and "NOMPI"
    ! J.Schwinger,      *GFI, Bergen*     2014-05-21
    ! - adapted code for writing of two time level tracer
    !   and sediment fields
    ! A.Moree,          *GFI, Bergen*   2018-04-12
    ! - new version of carbon isotope code
    ! J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
    ! - added preformed and saturated DIC tracers
    ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    ! - added cappability to restart c-isotopes from scratch (from
    !   observed d13C and d14C). This is used if c-isotope fields are
    !   not found in the restart file.
    ! - consistently organised restart of CFC and natural tracers
    !   from scratch, i.e. for the case that CFC and natural tracers are
    !   not found in the restart file.
    ! - removed satn2o which is not needed to restart the model
    ! - added sediment bypass preprocessor option
    ! J.Schwinger,      *Uni Research, Bergen*   2018-08-23
    ! - added reading of atmosphere field for BOXATM
    ! M. Bentsen,       *NORCE, Bergen*          2020-05-03
    ! - changed ocean model from MICOM to BLOM
    !***********************************************************************************************

    use netcdf,             only: nf90_global,nf90_noerr,nf90_nowrite,nf90_close,                  &
                                  nf90_open,nf90_get_att,nf90_inq_varid
    use mod_xc,             only: nbdy,mnproc,iqr,jqr,xcbcst,xchalt
    use mod_dia,            only: iotype
    use mo_carbch,          only: co2star,co3,hi,satoxy,ocetra,atm,nathi
    use mo_control_bgc,     only: io_stdo_bgc,ldtbgc,use_cisonew,use_AGG,                          &
                                  use_BOXATM,use_BROMO,use_CFC,use_natDIC,use_sedbypass,           &
                                  use_extNcycle,use_pref_tracers,use_shelfsea_res_time,            &
                                  use_sediment_quality
    use mo_param1_bgc,      only: ialkali,ian2o,iano3,icalc,idet,idicsat,                          &
                                  idms,idoc,ifdust,igasnit,iiron,iopal,ioxygen,iphosph,iphy,       &
                                  iprefalk,iprefdic,iprefo2,iprefpo4,iprefsilica,ishelfage,        &
                                  isco212,isilica,izoo,nocetra,                                    &
                                  iadust,inos,iatmco2,iatmn2,iatmo2,ibromo,icfc11,icfc12,isf6,     &
                                  icalc13,icalc14,idet13,idet14,idoc13,idoc14,iphy13,iphy14,       &
                                  isco213,isco214,izoo13,izoo14,safediv,                           &
                                  issso13,issso14,isssc13,isssc14,ipowc13,ipowc14,                 &
                                  iatmc13,iatmc14,iatmnco2,inatalkali,inatcalc,inatsco212,         &
                                  ipowaal,ipowaic,ipowaox,ipowaph,ipowasi,ipown2,ipowno3,          &
                                  isssc12,issso12,issssil,issster,ks,ianh4,iano2,ipownh4,ipown2o,  &
                                  ipowno2,issso12_age
    use mo_vgrid,           only: kbo
    use mo_sedmnt,          only: sedhpl,prorca_mavg,burial,sedlay
    use mo_intfcblom,       only: sedlay2,powtra2,burial2,atm2,prorca_mavg2
    use mo_param_bgc,       only: bifr13_ini,bifr14_ini,c14fac,re1312,re14to,prei13,prei14
    use mo_netcdf_bgcrw,    only: read_netcdf_var
#ifdef PNETCDF
    use mod_xc,             only: mpicomm
#endif

    ! Arguments
    integer,          intent(in)    :: kpie                                              ! 1st dimension of model grid.
    integer,          intent(in)    :: kpje                                              ! 2nd dimension of model grid.
    integer,          intent(in)    :: kpke                                              ! 3rd (vertical) dimension of model grid.
    integer,          intent(in)    :: ntr                                               ! number of tracers in tracer field
    integer,          intent(in)    :: ntrbgc                                            ! number of biogechemical tracers in tracer field
    integer,          intent(in)    :: itrbgc                                            ! start index for biogeochemical tracers in tracer field
    real,             intent(inout) :: trc(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,2*kpke,ntr) ! initial/restart tracer field to be passed to the ocean model [mol/kg]
    integer,          intent(in)    :: kplyear                                           ! year  in ocean restart date
    integer,          intent(in)    :: kplmon                                            ! month in ocean restart date
    integer,          intent(in)    :: kplday                                            ! day   in ocean restart date
    real,             intent(in)    :: omask(kpie,kpje)                                  ! land/ocean mask
    character(len=*), intent(in)    :: rstfnm                                            ! restart file name-informations

    ! Local variables
    real, allocatable :: locetra(:,:,:,:)          ! local array for reading
    integer   :: errstat
    integer   :: restyear                          !  year of restart file
    integer   :: restmonth                         !  month of restart file
    integer   :: restday                           !  day of restart file
    integer   :: restdtoce                         !  time step number from bgc ocean file
    integer   :: idate(5),i,j,k
    logical   :: lread_cfc,lread_nat,lread_iso,lread_atm,lread_bro,lread_extn,lread_prefsi,        &
              &  lread_sedqual
    logical   :: lread_shelfage
    real      :: rco213,rco214,alpha14,beta13,beta14,d13C_atm,d14cat
    integer   :: ncid,ncstat,ncvarid

#ifdef PNETCDF
#   include <pnetcdf.inc>
#   include <mpif.h>
    integer*4, save  :: info=MPI_INFO_NULL
#endif
    character(len=3) :: stripestr
    character(len=9) :: stripestr2
    integer          :: ierr,testio
    integer          :: leninrstfn
    !
    ! Allocate and initialize local array for reading (locetra)
    !
    allocate(locetra(kpie,kpje,2*kpke,nocetra),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory for locetra allocation'
    locetra(:,:,:,:) = 0.0
    !
    ! Open netCDF data file
    !
    testio=0
    if(mnproc==1 .and. IOTYPE==0) then

      ncstat = NF90_OPEN(rstfnm,NF90_NOWRITE, ncid)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFR_BGC: Problem with netCDF1)')
        stop        '(AUFR_BGC: Problem with netCDF1)'
      endif
      !
      ! Read restart data : date
      !
      ncstat = NF90_GET_ATT(ncid, NF90_GLOBAL,'date', idate)
      if ( ncstat  /=  NF90_NOERR ) then
        call xchalt('(AUFR_BGC: Problem reading date of restart file)')
        stop        '(AUFR_BGC: Problem reading date of restart file)'
      endif
      restyear  = idate(1)
      restmonth = idate(2)
      restday   = idate(3)
      restdtoce = idate(4)
      ldtbgc = idate(5)
      write(io_stdo_bgc,*) ' '
      write(io_stdo_bgc,*) 'Date of bgc restart file : '
      write(io_stdo_bgc,*) ' year  = ',restyear
      write(io_stdo_bgc,*) ' month = ',restmonth
      write(io_stdo_bgc,*) ' day   = ',restday
      write(io_stdo_bgc,*) ' dtoce = ',restdtoce
      write(io_stdo_bgc,*) ' dtbgc = ',ldtbgc
      write(io_stdo_bgc,*) ' '

    else if(IOTYPE==1) then

#ifdef PNETCDF
      testio=1
      write(stripestr,('(i3)')) 16
      write(stripestr2,('(i9)')) 1024*1024
      call mpi_info_create(info,ierr)
      call mpi_info_set(info,'romio_ds_read','disable',ierr)
      call mpi_info_set(info,'romio_ds_write','disable',ierr)
      call mpi_info_set(info,"striping_factor",stripestr,ierr)
      call mpi_info_set(info,"striping_unit",stripestr2,ierr)

      ncstat = NFMPI_OPEN(mpicomm,rstfnm,NF_NOWRITE,INFO, ncid)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFR_BGC: Problem with netCDF1)')
        stop        '(AUFR_BGC: Problem with netCDF1)'
      endif
      !
      ! Read restart data : date
      !
      ncstat = NFMPI_GET_ATT_INT(ncid, NF_GLOBAL,'date', idate)
      if ( ncstat  /=  NF_NOERR ) then
        call xchalt('(AUFR_BGC: Problem reading date of restart file)')
        stop        '(AUFR_BGC: Problem reading date of restart file)'
      endif
      restyear  = idate(1)
      restmonth = idate(2)
      restday   = idate(3)
      restdtoce = idate(4)
      ldtbgc = idate(5)
      if(mnproc==1) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'Date of bgc restart file : '
        write(io_stdo_bgc,*) ' year  = ',restyear
        write(io_stdo_bgc,*) ' month = ',restmonth
        write(io_stdo_bgc,*) ' day   = ',restday
        write(io_stdo_bgc,*) ' dtoce = ',restdtoce
        write(io_stdo_bgc,*) ' dtbgc = ',ldtbgc
        write(io_stdo_bgc,*) ' '
      endif
#endif
      if(testio .eq. 0) then
        call xchalt('(AUFR_BGC: Problem with namelist iotype)')
        stop        '(AUFR_BGC: Problem with namelist iotype)'
      endif

    endif ! mnproc==1 .and. IOTYPE==0

    !
    ! Compare with date read from ocean restart file
    !
    if (mnproc.eq.1) then

      if ( kplyear  /=  restyear  ) write(io_stdo_bgc,*)                                           &
           & 'WARNING: restart years in oce/bgc are not the same : ', kplyear,'/',restyear,' !!!'

      if ( kplmon   /=  restmonth ) write(io_stdo_bgc,*)                                           &
           & 'WARNING: restart months in oce/bgc are not the same : ',kplmon,'/',restmonth,' !!!'

      if ( kplday   /=  restday   ) write(io_stdo_bgc,*)                                           &
           & 'WARNING: restart days in oce/bgc are not the same : ',  kplday,'/',restday,' !!!'

    endif

    ! Find out whether to restart CFCs
    if (use_CFC) then
      lread_cfc=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'cfc11',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_cfc=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'cfc11',ncvarid)
        if(ncstat.ne.nf_noerr) lread_cfc=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_cfc) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: CFC tracers not in restart file, '
        write(io_stdo_bgc,*) ' CFCs initialised to zero.'
      endif
    endif

    ! Find out whether to restart natural tracers
    if (use_natDIC) then
      lread_nat=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'natsco212',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_nat=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'natsco212',ncvarid)
        if(ncstat.ne.nf_noerr) lread_nat=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_nat) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: natural tracers not in restart file. '
        write(io_stdo_bgc,*) ' Initialising natural tracers with their non-natural '
        write(io_stdo_bgc,*) ' counterpart.'
      endif
    endif

    ! Find out whether to restart marine carbon isotopes
    if (use_cisonew) then
      lread_iso=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'sco213',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_iso=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'sco213',ncvarid)
        if(ncstat.ne.nf_noerr) lread_iso=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_iso) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: carbon isotopes not in restart file. '
        write(io_stdo_bgc,*) ' Initialising carbon isotopes from scratch '
      endif
    endif

    ! Find out whether to restart Bromoform
    if (use_BROMO) then
      lread_bro=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'bromo',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_bro=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'bromo',ncvarid)
        if(ncstat.ne.nf_noerr) lread_bro=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_bro) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: Bromoform tracer not in restart file, '
        write(io_stdo_bgc,*) ' Initialised to 0.01 pmol L-1 (Stemmler et al., 2015).'
      endif
    endif

! Find out whether to restart extended nitrogen cycle
    if (use_extNcycle) then
      lread_extn=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'anh4',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_extn=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'anh4',ncvarid)
        if(ncstat.ne.nf_noerr) lread_extn=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_extn) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: extended nitrogen cycle tracer not in restart file '
        write(io_stdo_bgc,*) 'Initialising extended nitrogen cycle from scratch'
      endif
    endif

    ! Find out whether to restart atmosphere
    if (use_BOXATM) then
      lread_atm=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'atmco2',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_atm=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'atmco2',ncvarid)
        if(ncstat.ne.nf_noerr) lread_atm=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_atm) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: atmosphere fields not in restart file. '
        write(io_stdo_bgc,*) ' Initialising atmosphere from scratch '
      endif
    endif

      lread_prefsi=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'prefsilica',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_prefsi=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'prefsilica',ncvarid)
        if(ncstat.ne.nf_noerr) lread_prefsi=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_prefsi) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: preformed silica not in restart file '
        write(io_stdo_bgc,*) 'Initialising preformed tracer from scratch'
      endif

    if (use_shelfsea_res_time) then
      lread_shelfage=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'shelfage',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_shelfage=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'shelfage',ncvarid)
        if(ncstat.ne.nf_noerr) lread_shelfage=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_shelfage) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: shelfage not in restart file '
        write(io_stdo_bgc,*) 'Initialising shelfage from scratch'
      endif
    endif

    if (use_sediment_quality) then
      lread_sedqual=.true.
      if(IOTYPE==0) then
        if(mnproc==1) ncstat=nf90_inq_varid(ncid,'ssso12_age',ncvarid)
        call xcbcst(ncstat)
        if(ncstat.ne.nf90_noerr) lread_sedqual=.false.
      else if(IOTYPE==1) then
#ifdef PNETCDF
        ncstat=nfmpi_inq_varid(ncid,'ssso12_age',ncvarid)
        if(ncstat.ne.nf_noerr) lread_sedqual=.false.
#endif
      endif
      if(mnproc==1 .and. .not. lread_sedqual) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'AUFR_BGC info: ssso12_age not in restart file '
        write(io_stdo_bgc,*) 'Initialising ssso12_age from scratch'
      endif
    endif

    !
    ! Read restart data : ocean aquateous tracer
    !
    call read_netcdf_var(ncid,'sco212',locetra(1,1,1,isco212),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'alkali',locetra(1,1,1,ialkali),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'phosph',locetra(1,1,1,iphosph),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'oxygen',locetra(1,1,1,ioxygen),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'gasnit',locetra(1,1,1,igasnit),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'ano3',locetra(1,1,1,iano3),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'silica',locetra(1,1,1,isilica),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'doc',locetra(1,1,1,idoc),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'phyto',locetra(1,1,1,iphy),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'grazer',locetra(1,1,1,izoo),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'poc',locetra(1,1,1,idet),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'calciu',locetra(1,1,1,icalc),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'opal',locetra(1,1,1,iopal),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'n2o',locetra(1,1,1,ian2o),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'dms',locetra(1,1,1,idms),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'fdust',locetra(1,1,1,ifdust),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'iron',locetra(1,1,1,iiron),2*kpke,0,iotype)
    call read_netcdf_var(ncid,'dicsat',locetra(1,1,1,idicsat),2*kpke,0,iotype)
    if (use_pref_tracers) then
      call read_netcdf_var(ncid,'prefo2',locetra(1,1,1,iprefo2),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'prefpo4',locetra(1,1,1,iprefpo4),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'prefalk',locetra(1,1,1,iprefalk),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'prefdic',locetra(1,1,1,iprefdic),2*kpke,0,iotype)
      if(lread_prefsi) then
        call read_netcdf_var(ncid,'prefsilica',locetra(1,1,1,iprefsilica),2*kpke,0,iotype)
      endif
    endif
    if (use_shelfsea_res_time .and. lread_shelfage) then
      call read_netcdf_var(ncid,'shelfage',locetra(1,1,1,ishelfage),2*kpke,0,iotype)
    endif
    if (use_cisonew .and. lread_iso) then
      call read_netcdf_var(ncid,'sco213',locetra(1,1,1,isco213),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'sco214',locetra(1,1,1,isco214),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'doc13',locetra(1,1,1,idoc13),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'doc14',locetra(1,1,1,idoc14),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'phyto13',locetra(1,1,1,iphy13),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'phyto14',locetra(1,1,1,iphy14),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'grazer13',locetra(1,1,1,izoo13),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'grazer14',locetra(1,1,1,izoo14),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'poc13',locetra(1,1,1,idet13),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'poc14',locetra(1,1,1,idet14),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'calciu13',locetra(1,1,1,icalc13),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'calciu14',locetra(1,1,1,icalc14),2*kpke,0,iotype)
    endif
    if (use_AGG)then
      call read_netcdf_var(ncid,'snos',locetra(1,1,1,inos),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'adust',locetra(1,1,1,iadust),2*kpke,0,iotype)
    endif
    if (use_CFC .and. lread_cfc) then
      call read_netcdf_var(ncid,'cfc11',locetra(1,1,1,icfc11),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'cfc12',locetra(1,1,1,icfc12),2*kpke,0,iotype)
      call read_netcdf_var(ncid,'sf6',locetra(1,1,1,isf6),2*kpke,0,iotype)
    endif
    if (use_natDIC) then
      if (lread_nat) then
        call read_netcdf_var(ncid,'natsco212',locetra(1,1,1,inatsco212),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'natalkali',locetra(1,1,1,inatalkali),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'natcalciu',locetra(1,1,1,inatcalc),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'nathi',nathi(1,1,1),kpke,0,iotype)
      else
        call read_netcdf_var(ncid,'sco212',locetra(1,1,1,inatsco212),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'alkali',locetra(1,1,1,inatalkali),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'calciu',locetra(1,1,1,inatcalc),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'hi',nathi(1,1,1),kpke,0,iotype)
      endif
    endif
    if (use_BROMO .and. lread_bro) then
      call read_netcdf_var(ncid,'bromo',locetra(1,1,1,ibromo),2*kpke,0,iotype)
    endif
    if (use_extNcycle) then
      if(lread_extn) then
        call read_netcdf_var(ncid,'anh4',locetra(1,1,1,ianh4),2*kpke,0,iotype)
        call read_netcdf_var(ncid,'ano2',locetra(1,1,1,iano2),2*kpke,0,iotype)
      endif
    endif

    !
    ! Read restart data : diagnostic ocean fields (needed for bit to bit reproducability)
    !
    call read_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0,iotype)
    call read_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0,iotype)
    call read_netcdf_var(ncid,'co2star',co2star(1,1,1),kpke,0,iotype)
    call read_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0,iotype)
    !
    ! Read restart data : sediment variables.
    !
    if (.not. use_sedbypass) then
      call read_netcdf_var(ncid,'ssso12',sedlay2(1,1,1,issso12),2*ks,0,iotype)
      call read_netcdf_var(ncid,'sssc12',sedlay2(1,1,1,isssc12),2*ks,0,iotype)
      call read_netcdf_var(ncid,'ssssil',sedlay2(1,1,1,issssil),2*ks,0,iotype)
      call read_netcdf_var(ncid,'ssster',sedlay2(1,1,1,issster),2*ks,0,iotype)
      call read_netcdf_var(ncid,'bur_o12',burial2(1,1,1,issso12),2,0,iotype)
      call read_netcdf_var(ncid,'bur_c12',burial2(1,1,1,isssc12),2,0,iotype)
      call read_netcdf_var(ncid,'bur_sil',burial2(1,1,1,issssil),2,0,iotype)
      call read_netcdf_var(ncid,'bur_clay',burial2(1,1,1,issster),2,0,iotype)
      call read_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0,iotype)
      call read_netcdf_var(ncid,'powaic',powtra2(1,1,1,ipowaic),2*ks,0,iotype)
      call read_netcdf_var(ncid,'powaal',powtra2(1,1,1,ipowaal),2*ks,0,iotype)
      call read_netcdf_var(ncid,'powaph',powtra2(1,1,1,ipowaph),2*ks,0,iotype)
      call read_netcdf_var(ncid,'powaox',powtra2(1,1,1,ipowaox),2*ks,0,iotype)
      call read_netcdf_var(ncid,'pown2',powtra2(1,1,1,ipown2),2*ks,0,iotype)
      call read_netcdf_var(ncid,'powno3',powtra2(1,1,1,ipowno3),2*ks,0,iotype)
      call read_netcdf_var(ncid,'powasi',powtra2(1,1,1,ipowasi),2*ks,0,iotype)
      if (use_cisonew .and. lread_iso) then
        call read_netcdf_var(ncid,'ssso13',sedlay2(1,1,1,issso13),2*ks,0,iotype)
        call read_netcdf_var(ncid,'ssso14',sedlay2(1,1,1,issso14),2*ks,0,iotype)
        call read_netcdf_var(ncid,'sssc13',sedlay2(1,1,1,isssc13),2*ks,0,iotype)
        call read_netcdf_var(ncid,'sssc14',sedlay2(1,1,1,isssc14),2*ks,0,iotype)
        call read_netcdf_var(ncid,'bur_o13',burial2(1,1,1,issso13),2,0,iotype)
        call read_netcdf_var(ncid,'bur_o14',burial2(1,1,1,issso14),2,0,iotype)
        call read_netcdf_var(ncid,'bur_c13',burial2(1,1,1,isssc13),2,0,iotype)
        call read_netcdf_var(ncid,'bur_c14',burial2(1,1,1,isssc14),2,0,iotype)
        call read_netcdf_var(ncid,'powc13',powtra2(1,1,1,ipowc13),2*ks,0,iotype)
        call read_netcdf_var(ncid,'powc14',powtra2(1,1,1,ipowc14),2*ks,0,iotype)
      endif
      if (use_sediment_quality .and. lread_sedqual) then
        call read_netcdf_var(ncid,'ssso12_age',sedlay2(1,1,1,issso12_age),2*ks,0,iotype)
        call read_netcdf_var(ncid,'bur_o12_age',burial2(1,1,1,issso12_age),2,0,iotype)
        call read_netcdf_var(ncid,'prorca_mavg',prorca_mavg2(1,1,1),2,0,iotype)
      endif
      if (use_extNcycle) then
        if(lread_extn) then
          call read_netcdf_var(ncid,'pownh4',powtra2(1,1,1,ipownh4),2*ks,0,iotype)
          call read_netcdf_var(ncid,'pown2o',powtra2(1,1,1,ipown2o),2*ks,0,iotype)
          call read_netcdf_var(ncid,'powno2',powtra2(1,1,1,ipowno2),2*ks,0,iotype)
        endif
      endif
    endif
    !
    ! Read restart data: atmosphere
    !
    if (use_BOXATM) then
      if(lread_atm) then
        call read_netcdf_var(ncid,'atmco2',atm2(1,1,1,iatmco2),2,0,iotype)
        call read_netcdf_var(ncid,'atmo2',atm2(1,1,1,iatmo2),2,0,iotype)
        call read_netcdf_var(ncid,'atmn2',atm2(1,1,1,iatmn2),2,0,iotype)
        if (use_cisonew) then
          if(lread_iso) then
            call read_netcdf_var(ncid,'atmc13',atm2(1,1,1,iatmc13),2,0,iotype)
            call read_netcdf_var(ncid,'atmc14',atm2(1,1,1,iatmc14),2,0,iotype)
          else
            ! If atm isotopes are not in restart but boxatm is on, calculate initial value using atmco2
            ! that is just read in from restart files. Normalize atmc14 using beleg c14fac.
            do j=1,kpje
              do i=1,kpie
                beta13 = (prei13/1000.)+1.
                alpha14 = 2.*(prei13+25.)
                d14cat  = (prei14+alpha14)/(1.-alpha14/1000.)
                atm(i,j,iatmc13) = beta13*re1312*atm2(i,j,1,iatmco2)/(1.+beta13*re1312)
                atm(i,j,iatmc14) = ((d14cat/1000.)+1.)*re14to*atm2(i,j,1,iatmco2)/c14fac
              enddo
            enddo
            ! Copy the isotope atmosphere fields into both timelevels of atm2.
            atm2(:,:,1,iatmc13) = atm(:,:,iatmc13)
            atm2(:,:,2,iatmc13) = atm(:,:,iatmc13)
            atm2(:,:,1,iatmc14) = atm(:,:,iatmc14)
            atm2(:,:,2,iatmc14) = atm(:,:,iatmc14)
          endif
        endif
        if (use_natDIC) then
          call read_netcdf_var(ncid,'atmnco2',atm2(1,1,1,iatmnco2),2,0,iotype)
        endif
      else
        ! If atmosphere field is not in restart, copy the atmosphere field
        ! (initialised in beleg.F90) into both timelevels of atm2.
        atm2(:,:,1,:) = atm(:,:,:)
        atm2(:,:,2,:) = atm(:,:,:)
      endif
    endif

    if(mnproc==1 .and. IOTYPE==0) then
      ncstat = NF90_CLOSE(ncid)
    else if(IOTYPE==1) then
#ifdef PNETCDF
      ncstat = NFMPI_CLOSE(ncid)
#endif
    endif

    if (use_cisonew .and. .not. lread_iso) then
      ! If carbon isotope fields are not read from restart file, copy the d13C
      ! d14C fields (initialised in beleg.F90) into both timelevels of locetra.
      locetra(:,:,1:kpke,       isco213)=ocetra(:,:,:,isco213)
      locetra(:,:,kpke+1:2*kpke,isco213)=ocetra(:,:,:,isco213)
      locetra(:,:,1:kpke,       isco214)=ocetra(:,:,:,isco214)
      locetra(:,:,kpke+1:2*kpke,isco214)=ocetra(:,:,:,isco214)
      ! Initialise 13C and 14C fields in the same way as in beleg.F90
      do k=1,2*kpke
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j)  >  0.5) then
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
              locetra(i,j,k,idoc13)=locetra(i,j,k,idoc)*rco213*bifr13_ini
              locetra(i,j,k,idoc14)=locetra(i,j,k,idoc)*rco214*bifr14_ini
              locetra(i,j,k,iphy13)=locetra(i,j,k,iphy)*rco213*bifr13_ini
              locetra(i,j,k,iphy14)=locetra(i,j,k,iphy)*rco214*bifr14_ini
              locetra(i,j,k,izoo13)=locetra(i,j,k,izoo)*rco213*bifr13_ini
              locetra(i,j,k,izoo14)=locetra(i,j,k,izoo)*rco214*bifr14_ini
              locetra(i,j,k,idet13)=locetra(i,j,k,idet)*rco213*bifr13_ini
              locetra(i,j,k,idet14)=locetra(i,j,k,idet)*rco214*bifr14_ini
              locetra(i,j,k,icalc13)=locetra(i,j,k,icalc)*rco213
              locetra(i,j,k,icalc14)=locetra(i,j,k,icalc)*rco214
            endif
          enddo
        enddo
      enddo

      if (.not. use_sedbypass) then
        do  k=1,2*ks
          do  j=1,kpje
            do  i=1,kpie
              if(omask(i,j)  >  0.5) then
                rco213=locetra(i,j,kbo(i,j),isco213)/(locetra(i,j,kbo(i,j),isco212)+safediv)
                rco214=locetra(i,j,kbo(i,j),isco214)/(locetra(i,j,kbo(i,j),isco212)+safediv)
                powtra2(i,j,k,ipowc13)=powtra2(i,j,k,ipowaic)*rco213
                powtra2(i,j,k,ipowc14)=powtra2(i,j,k,ipowaic)*rco214
                sedlay2(i,j,k,issso13)=sedlay2(i,j,k,issso12)*rco213*bifr13_ini
                sedlay2(i,j,k,issso14)=sedlay2(i,j,k,issso12)*rco214*bifr14_ini
                sedlay2(i,j,k,isssc13)=sedlay2(i,j,k,isssc12)*rco213
                sedlay2(i,j,k,isssc14)=sedlay2(i,j,k,isssc12)*rco214
              endif
            enddo
          enddo
        enddo

        do  k=1,2
          do  j=1,kpje
            do  i=1,kpie
              if(omask(i,j)  >  0.5) then
                rco213=locetra(i,j,kbo(i,j),isco213)/(locetra(i,j,kbo(i,j),isco212)+safediv)
                rco214=locetra(i,j,kbo(i,j),isco214)/(locetra(i,j,kbo(i,j),isco212)+safediv)
                burial2(i,j,k,issso13)=burial2(i,j,k,issso12)*rco213*bifr13_ini
                burial2(i,j,k,issso14)=burial2(i,j,k,issso12)*rco214*bifr14_ini
                burial2(i,j,k,isssc13)=burial2(i,j,k,isssc12)*rco213
                burial2(i,j,k,isssc14)=burial2(i,j,k,isssc12)*rco214
              endif
            enddo
          enddo
        enddo

      endif  ! .not. use_sedbypass
    endif ! use_cisonew .and. .not. lread_iso

    if (.not. use_sedbypass) then
      if (use_sediment_quality .and. .not. lread_sedqual) then
        ! if (hybrid) restart, but age and mavgs should be started with filled values, if provided
        do j=1,kpje
          do i=1,kpie
            if (omask(i,j) > 0.) then
              burial2(i,j,1,issso12_age)      = burial(i,j,issso12_age)
              burial2(i,j,2,issso12_age)      = burial(i,j,issso12_age)
              prorca_mavg2(i,j,1)             = prorca_mavg(i,j)
              prorca_mavg2(i,j,2)             = prorca_mavg(i,j)
              do k=1,ks
                sedlay2(i,j,k,issso12_age)    = sedlay(i,j,k,issso12_age)
                sedlay2(i,j,ks+k,issso12_age) = sedlay(i,j,k,issso12_age)
              enddo
            endif
          enddo
        enddo
      endif
    endif

    ! return tracer fields to ocean model (both timelevels); No unit
    ! conversion here, since tracers in the restart file are in
    ! BLOM units (mol/kg)
    !--------------------------------------------------------------------
    trc(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)=locetra(:,:,:,:)
    deallocate(locetra)

  end subroutine aufr_bgc

end module mo_aufr_bgc
