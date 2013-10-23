      SUBROUTINE AUFW_BGC(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,             &
     &                    trc,pddpo,pglon,pglat,                        &
     &                    kplyear,kplmon,kplday,kpldtoce,omask,         &
     &                    rstfnm_ocn,path)

!****************************************************************
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
!     - code cleanup, remoded preprocessor option "PNETCDF"
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
!     *REAL*    *trc*        - initial/restart tracer field to be passed to the 
!                              ocean model [mol/kg]
!     *REAL*    *pddpo*      - size of grid cell (3rd dimension) [m].
!     *REAL*    *pglon*      - geographical longitude of grid points [degree E].
!     *REAL*    *pglat*      - geographical latitude  of grid points [degree N].

!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmon*     - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *INTEGER* *kpldtoce*   - step  in ocean restart date
!     *REAL*    *omask*      - land/ocean mask
!     *CHAR*    *rstfnm_ocn* - restart file name-informations
!     *CHAR*    *path*       - path to restart files
!
!**************************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc, only: nbdy,itdm,jtdm,mnproc,xchalt

      implicit none

      INTEGER           :: kpie,kpje,kpke,ntr,ntrbgc,itrbgc
      REAL              :: trc(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,kpke,ntr)
      REAL              :: pddpo(kpie,kpje,kpke)    
      REAL              :: pglon(kpie,kpje)
      REAL              :: pglat(kpie,kpje)
      REAL              :: omask(kpie,kpje)    
      INTEGER           :: kplyear,kplmon,kplday,kpldtoce
      character(len=*)  :: rstfnm_ocn,path

      INTEGER           :: i,j,k,l,jj,ii,kk,kt
      CHARACTER(LEN=80) :: err_text,rstfnm

      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncvarid,ncstat,ncoldmod,ncdims(4)                  &
     &       ,nclatid,nclonid,nclevid,nclev1id                        &
     &       ,nctraid,ncksid,ncsedid                                  &
     &       ,nstart2(2),ncount2(2),nstride2(2),idate(5)
      REAL zfield(kpie,kpje)
      REAL rmissing


! pass tracer fields in from ocean model; No unit conversion here, 
! tracers in the restart file are written in mol/kg 
!--------------------------------------------------------------------
!
      ocetra(:,:,:,:)=trc(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)



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

#ifdef DIFFAT
!
!  Masking co2.
!      
      DO  j=1,kpje
      DO  i=1,kpie 
      IF(omask(i,j) .LT. 0.5) THEN
      suppco2(i,j)=rmissing
      ENDIF
      ENDDO
      ENDDO      
#endif      

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
            call xchalt('(aufw_bgc)')
            stop '(aufw_bgc)'
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
            call xchalt('(aufw_bgc)')
            stop '(aufw_bgc)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'_rest_b_'//rstfnm_ocn(i+9:)
#endif

      write(io_stdo_bgc,*) 'BGC RESTART   ',rstfnm
      ncstat = NF_CREATE(trim(path)//rstfnm,NF_CLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF1)')
               stop '(AUFW: Problem with netCDF1)'
      ENDIF

!
! Define dimension
! ----------------------------------------------------------------------    
!
      ncstat = NF_DEF_DIM(ncid, 'lon', itdm, nclonid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF2)')
               stop '(AUFW: Problem with netCDF2)'
      ENDIF

      ncstat = NF_DEF_DIM(ncid, 'lat', jtdm, nclatid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF3)')
               stop '(AUFW: Problem with netCDF3)'
      ENDIF

      ncstat = NF_DEF_DIM(ncid, 'depth', kpke, nclevid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF4)')
               stop '(AUFW: Problem with netCDF4)'
      ENDIF

      ncstat = NF_DEF_DIM(ncid, 'ntra', nocetra, nctraid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF5)')
               stop '(AUFW: Problem with netCDF5)'
      ENDIF

      ncstat = NF_DEF_DIM(ncid, 'nks', ks, ncksid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF7)')
               stop '(AUFW: Problem with netCDF7)'
      ENDIF

      ncstat = NF_DEF_DIM(ncid, 'nsed',nsedtra, ncsedid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF8)')
               stop '(AUFW: Problem with netCDF8)'
      ENDIF
     
      ncstat = NF_DEF_DIM(ncid, 'lev1', 1, nclev1id)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF8b)')
               stop '(AUFW: Problem with netCDF8b)'
      ENDIF

!
! Define global attributes
! ----------------------------------------------------------------------    
!
      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'               &
     &,35, 'Restart data for marine bgc modules') 
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'history'             &
     &,35, 'Restart data for marine bgc modules')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'conventions'         &
     &,6, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'              &
     &,24, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF10)')
               stop '(AUFW: Problem with netCDF10)'
      ENDIF

      ncstat = NF_PUT_ATT_INT(ncid, NF_GLOBAL, 'date', NF_INT, 5, idate)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF11)')
               stop '(AUFW: Problem with netCDF11)'
      ENDIF

!
! Define variables : grid
! ----------------------------------------------------------------------    
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid

      ncstat = NF_DEF_VAR(ncid,'scal_lon',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF13)')
               stop '(AUFW: Problem with netCDF13)'
      ENDIF
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree E')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF14)')
               stop '(AUFW: Problem with netCDF14)'
      ENDIF
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,34, '2-d longitude of scalar grid cells')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF15)')
               stop '(AUFW: Problem with netCDF15)'
      ENDIF

      ncstat = NF_DEF_VAR(ncid,'scal_lat',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF16)')
               stop '(AUFW: Problem with netCDF16)'
      ENDIF
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree N')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF17)')
               stop '(AUFW: Problem with netCDF17)'
      ENDIF
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,33, '2-d latitude of scalar grid cells')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF18)')
               stop '(AUFW: Problem with netCDF18)'
      ENDIF

      ncstat = NF_DEF_VAR(ncid,'scal_wdep',NF_DOUBLE,2,ncdims,ncvarid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF16a)')
               stop '(AUFW: Problem with netCDF16a)'
      ENDIF
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF17a)')
               stop '(AUFW: Problem with netCDF17a)'
      ENDIF
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,37, '2-d water depth at scalar grid points')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF18a)')
               stop '(AUFW: Problem with netCDF18a)'
      ENDIF

!
! Define variables : advected ocean tracer
! ----------------------------------------------------------------------    
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid

      CALL NETCDF_DEF_VARDB(ncid,6,'sco212',3,ncdims,ncvarid,           &
     &    6,'mol/kg',13, 'Dissolved CO2',rmissing,22,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'sco213',3,ncdims,ncvarid,           &
     &    6,'mol/kg',15, 'Dissolved CO213',rmissing,22,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sco214',3,ncdims,ncvarid,           &
     &    6,'mol/kg',15, 'Dissolved CO214',rmissing,22,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'alkali',3,ncdims,ncvarid,           &
     &    6,'mol/kg',10,'Alkalinity',rmissing,25,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'phosph',3,ncdims,ncvarid,           &
     &    6,'mol/kg',19,'Dissolved phosphate',rmissing,28,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'oxygen',3,ncdims,ncvarid,           &
     &    6,'mol/kg',16,'Dissolved oxygen',                             &
          rmissing,31,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'gasnit',3,ncdims,ncvarid,           &
     &    6,'mol/kg',21,'Gaseous nitrogen (N2)',                        &
          rmissing,34,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'ano3',3,ncdims,ncvarid,             &
     &    6,'mol/kg',17,'Dissolved nitrate',                            &
          rmissing,34,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'silica',3,ncdims,ncvarid,           &
     &    6,'mol/kg',22,'Silicid acid (Si(OH)4)',                       &
          rmissing,40,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'doc',3,ncdims,ncvarid,              &
     &    6,'mol/kg',24,'Dissolved organic carbon',                     &
     &    rmissing,40,io_stdo_bgc) 

      CALL NETCDF_DEF_VARDB(ncid,3,'poc',3,ncdims,ncvarid,              &
     &    6,'mol/kg',25,'Particulate organic carbon',                   &
     &    rmissing,46,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,5,'poc13',3,ncdims,ncvarid,            &
     &    7,'molC/kg',28,'Particulate organic carbon13',                &
     &    rmissing,46,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'poc14',3,ncdims,ncvarid,            &
     &    7,'molC/kg',28,'Particulate organic carbon14',                &
     &    rmissing,46,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,2,'hi',3,ncdims,ncvarid,               &
     &    6,'mol/kg',26,'Hydrogen ion concentration',                   &
     &    rmissing,46,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'co3',3,ncdims,ncvarid,              &
     &    6,'mol/kg',25,'Dissolved carbonate (CO3)',                    &
     &    rmissing,52,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'phyto',3,ncdims,ncvarid,            &
     &    7,'molP/kg',27,'Phytoplankton concentration',                 &
     &    rmissing,28,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'grazer',3,ncdims,ncvarid,           &
     &    7,'molP/kg',25,'Zooplankton concentration',                   &
     &    rmissing,29,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'calciu',3,ncdims,ncvarid,           &
     &    6,'mol/kg',17,'Calcium carbonate',                            &
     &    rmissing,30,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,8,'calciu13',3,ncdims,ncvarid,         &
     &    7,'molC/kg',19,'Calcium carbonate13',                         &
     &    rmissing,30,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'calciu14',3,ncdims,ncvarid,         &
     &    7,'molC/kg',19,'Calcium carbonate14',                         &
     &    rmissing,30,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,4,'opal',3,ncdims,ncvarid,             &
     &    6,'mol/kg',15,'Biogenic silica',                              &
     &    rmissing,31,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'n2o',3,ncdims,ncvarid,              &
     &    6,'mol/kg',12,'laughing gas',                                 &
     &    rmissing,32,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'dms',3,ncdims,ncvarid,              &
     &    6,'mol/kg',15 ,'DiMethylSulfide',                             &
     &    rmissing,33,io_stdo_bgc)            

      CALL NETCDF_DEF_VARDB(ncid,5,'fdust',3,ncdims,ncvarid,            &
     &    5,'kg/kg',19,'Non-aggregated dust',                           &
     &    rmissing,34,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'iron',3,ncdims,ncvarid,             &
     &    6,'mol/kg',14,'Dissolved iron',                               &
     &    rmissing,35,io_stdo_bgc)

#ifdef AGG
      CALL NETCDF_DEF_VARDB(ncid,4,'snos',3,ncdims,ncvarid,             &
     &    3,'1/g',38,'marine snow aggregates per g sea water',          &
     &    rmissing,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'adust',3,ncdims,ncvarid,            &
     &    4,'g/kg',15,'Aggregated dust',                                &
     &    rmissing,42,io_stdo_bgc)
#endif /*AGG*/   

#ifdef ANTC14
      CALL NETCDF_DEF_VARDB(ncid,6,'antc14',3,ncdims,ncvarid,           &
     &    6,'mol/kg',17,'anthropogenic C14',                            &
     &    rmissing,41,io_stdo_bgc)
#endif
#ifdef CFC
      CALL NETCDF_DEF_VARDB(ncid,5,'cfc11',3,ncdims,ncvarid,            &
     &    6,'mol/kg',5,'CFC11',                                         &
     &    rmissing,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'cfc12',3,ncdims,ncvarid,            &
     &    6,'mol/kg',5,'CFC12',                                         &
     &    rmissing,41,io_stdo_bgc)     

      CALL NETCDF_DEF_VARDB(ncid,3,'sf6',3,ncdims,ncvarid,              &
     &    6,'mol/kg',4,'SF-6',                                          &
     &    rmissing,41,io_stdo_bgc)     
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'satoxy',3,ncdims,ncvarid,           &
     &    9,'xxxxxxxxx',9 ,'xxxxxxxxx',  &
     &    rmissing,64,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'satn2o',2,ncdims,ncvarid,           &
     &    9,'xxxxxxxxx',9 ,'xxxxxxxxx',  &
     &    rmissing,64,io_stdo_bgc)

!
! Define variables : sediment
! ----------------------------------------------------------------------    
!
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = ncksid
      ncdims(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,6,'ssso12',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',35,'Sediment accumulated organic carbon',      &
     &    rmissing,69,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'ssso13',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',37,'Sediment accumulated organic carbon13',    &
     &    rmasks,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'ssso14',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',37,'Sediment accumulated organic carbon14',    &
     &    rmasks,69,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc12',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',38,'Sediment accumulated calcium carbonate',   &
     &    rmissing,69,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'sssc13',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',40,'Sediment accumulated calcium carbonate13', &
     &    rmasks,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sssc14',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',40,'Sediment accumulated calcium carbonate14', &
     &    rmasks,69,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'ssssil',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',25,'Sediment accumulated opal',                &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'ssster',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',25,'Sediment accumulated clay',                &
     &    rmissing,69,io_stdo_bgc)

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = 0
      ncdims(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_o12',2,ncdims,ncvarid,         &
     &    9,'kmol/m**2',30,'Burial layer of organic carbon',           &
     &    rmissing,70,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_c12',2,ncdims,ncvarid,         &
     &    9,'kmol/m**2',33,'Burial layer of calcium carbonate',        &
     &    rmissing,71,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'bur_sil',2,ncdims,ncvarid,         &
     &    9,'kmol/m**2',20,'Burial layer of opal',                     &
     &    rmissing,72,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'bur_clay',2,ncdims,ncvarid,        &
     &    9,'kmol/m**2',20,'Burial layer of clay',                     &
     &    rmissing,72,io_stdo_bgc)

      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = ncksid

      CALL NETCDF_DEF_VARDB(ncid,6,'sedhpl',3,ncdims,ncvarid,          &
     &    9,'kmol/m**2',34,'Sediment accumulated hydrogen ions',       &
     &    rmissing,72,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaic',3,ncdims,ncvarid,          &
     &    9,'kmol/m**3',23,'Sediment pore water CO2',                  &
     &    rmissing,75,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'powc13',3,ncdims,ncvarid,          &
     &    9,'kmol/m**3',25,'Sediment pore water DIC13',                &
     &    rmasks,75,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powc14',3,ncdims,ncvarid,          &
     &    9,'kmol/m**3',25,'Sediment pore water DIC13',                &
     &    rmasks,75,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'powaal',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',30,'Sediment pore water alkalinity',        &
     &    rmissing,78,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaph',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',29,'Sediment pore water phosphate',         &
     &    rmissing,78,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powaox',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',26,'Sediment pore water oxygen',            &
     &    rmissing,84,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'pown2',3,ncdims,ncvarid,        &
     &    9,'kmol/m**3',36,'Sediment pore water gaseous nitrogen',  &
     &    rmissing,87,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powno3',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',33,'Sediment pore water nitrate (NO3)',     &
     &    rmissing,90,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'powasi',3,ncdims,ncvarid,       &
     &    9,'kmol/m**3',42,'Sediment pore water silicid acid (Si(OH)4)', &
     &    rmissing,91,io_stdo_bgc)

#ifdef DIFFAT
!
! Define variables : co2 diffusion
! ----------------------------------------------------------------------    
!
      CALL NETCDF_DEF_VARDB(ncid,7,'suppco2',2,ncdims,ncvarid,      &
     &    4,'ppmv',42,'pCO2 from total dissolved inorganic carbon', &
     &    rmissing,92,io_stdo_bgc)
     
      CALL NETCDF_DEF_VARDB(ncid,6,'atmco2',2,ncdims,ncvarid,       &
     &    3,'ppm',15,'atmospheric CO2',                             &
     &    rmissing,93,io_stdo_bgc)         
     
      CALL NETCDF_DEF_VARDB(ncid,5,'atmo2',2,ncdims,ncvarid,        &
     &    3,'ppm',14,'atmospheric O2',                              &
     &    rmissing,94,io_stdo_bgc)    
     
      CALL NETCDF_DEF_VARDB(ncid,5,'atmn2',2,ncdims,ncvarid,        &
     &    3,'ppm',14,'atmospheric N2',                              &
     &    rmissing,95,io_stdo_bgc)    

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'atmc13',2,ncdims,ncvarid,       &
     &    3,'ppm',15,'atmospheric CO2',                             &
     &    rmissing,93,io_stdo_bgc)      
      CALL NETCDF_DEF_VARDB(ncid,6,'atmc14',2,ncdims,ncvarid,       &
     &    3,'ppm',15,'atmospheric CO2',                             &
     &    rmissing,93,io_stdo_bgc)      
#endif

#endif

      ncstat = NF_ENDDEF(ncid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF00)')
               stop '(AUFW: Problem with netCDF00)'
      ENDIF


!
! Set fill mode
! ----------------------------------------------------------------------    
!
      ncstat = NF_SET_FILL(ncid,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF97)')
               stop '(AUFW: Problem with netCDF97)'
      ENDIF

      ENDIF ! mnproc == 1

!
! Write grid describing data
! ----------------------------------------------------------------------    
!

      CALL write_netcdf_var(ncid,'scal_lon', pglon(1,1),1,0)
      CALL write_netcdf_var(ncid,'scal_lat', pglat(1,1),1,0)

      zfield(:,:) = 0.0
      DO k=1,kpke
         DO j=1,kpje
            DO i=1,kpie
               zfield(i,j) = zfield(i,j) + pddpo(i,j,k)
            ENDDO
         ENDDO
      ENDDO
             
      CALL write_netcdf_var(ncid,'scal_wdep',zfield(1,1),1,0)

!
! Write restart data : ocean aquateous tracer
!--------------------------------------------------------------------
!
      CALL write_netcdf_var(ncid,'sco212',ocetra(1,1,1,isco212),kpke,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'sco213',ocetra(1,1,1,isco213),kpke,0)
      CALL write_netcdf_var(ncid,'sco214',ocetra(1,1,1,isco214),kpke,0)
#endif
      CALL write_netcdf_var(ncid,'alkali',ocetra(1,1,1,ialkali),kpke,0)
      CALL write_netcdf_var(ncid,'phosph',ocetra(1,1,1,iphosph),kpke,0)
      CALL write_netcdf_var(ncid,'oxygen',ocetra(1,1,1,ioxygen),kpke,0)
      CALL write_netcdf_var(ncid,'gasnit',ocetra(1,1,1,igasnit),kpke,0)
      CALL write_netcdf_var(ncid,'ano3',ocetra(1,1,1,iano3),kpke,0)
      CALL write_netcdf_var(ncid,'silica',ocetra(1,1,1,isilica),kpke,0)
      CALL write_netcdf_var(ncid,'doc',ocetra(1,1,1,idoc),kpke,0)
      CALL write_netcdf_var(ncid,'poc',ocetra(1,1,1,idet),kpke,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'poc13',ocetra(1,1,1,idet13),kpke,0)
      CALL write_netcdf_var(ncid,'poc14',ocetra(1,1,1,idet14),kpke,0)
#endif
      CALL write_netcdf_var(ncid,'phyto',ocetra(1,1,1,iphy),kpke,0)
      CALL write_netcdf_var(ncid,'grazer',ocetra(1,1,1,izoo),kpke,0)
      CALL write_netcdf_var(ncid,'calciu',ocetra(1,1,1,icalc),kpke,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'calciu13',ocetra(1,1,1,icalc13),kpke,0)
      CALL write_netcdf_var(ncid,'calciu14',ocetra(1,1,1,icalc14),kpke,0)
#endif
      CALL write_netcdf_var(ncid,'opal',ocetra(1,1,1,iopal),kpke,0)
      CALL write_netcdf_var(ncid,'n2o',ocetra(1,1,1,ian2o),kpke,0)
      CALL write_netcdf_var(ncid,'dms',ocetra(1,1,1,idms),kpke,0)
      CALL write_netcdf_var(ncid,'fdust',ocetra(1,1,1,ifdust),kpke,0)
      CALL write_netcdf_var(ncid,'iron',ocetra(1,1,1,iiron),kpke,0)
#ifdef AGG
      CALL write_netcdf_var(ncid,'snos',ocetra(1,1,1,inos),kpke,0)
      CALL write_netcdf_var(ncid,'adust',ocetra(1,1,1,iadust),kpke,0)
#endif /*AGG*/
#ifdef ANTC14
      CALL write_netcdf_var(ncid,'antc14',ocetra(1,1,1,iantc14),kpke,0)
#endif
#ifdef CFC
      CALL write_netcdf_var(ncid,'cfc11',ocetra(1,1,1,icfc11),kpke,0)
      CALL write_netcdf_var(ncid,'cfc12',ocetra(1,1,1,icfc12),kpke,0)
      CALL write_netcdf_var(ncid,'sf6',ocetra(1,1,1,isf6),kpke,0)
#endif


!
! Write restart data : diagtnostic ocean tracer
!
      CALL write_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
      CALL write_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
!
! Write restart data : other fields
!
      CALL write_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0)
      CALL write_netcdf_var(ncid,'satn2o',satn2o(1,1),1,0)

!
! Write restart data : sediment variables.
!
      CALL write_netcdf_var(ncid,'ssso12',sedlay(1,1,1,issso12),ks,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'ssso13',sedlay(1,1,1,issso13),ks,0)
      CALL write_netcdf_var(ncid,'ssso14',sedlay(1,1,1,issso14),ks,0)
#endif
      CALL write_netcdf_var(ncid,'sssc12',sedlay(1,1,1,isssc12),ks,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'sssc13',sedlay(1,1,1,isssc13),ks,0)
      CALL write_netcdf_var(ncid,'sssc14',sedlay(1,1,1,isssc14),ks,0)
#endif
      CALL write_netcdf_var(ncid,'ssssil',sedlay(1,1,1,issssil),ks,0)
      CALL write_netcdf_var(ncid,'ssster',sedlay(1,1,1,issster),ks,0)
      CALL write_netcdf_var(ncid,'bur_o12',burial(1,1,issso12),1,0)
      CALL write_netcdf_var(ncid,'bur_c12',burial(1,1,isssc12),1,0)
      CALL write_netcdf_var(ncid,'bur_sil',burial(1,1,issssil),1,0)
      CALL write_netcdf_var(ncid,'bur_clay',burial(1,1,issster),1,0)
      CALL write_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0)
      CALL write_netcdf_var(ncid,'powaic',powtra(1,1,1,ipowaic),ks,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'powc13',powtra(1,1,1,ipowc13),ks,0)
      CALL write_netcdf_var(ncid,'powc14',powtra(1,1,1,ipowc14),ks,0)
#endif
      CALL write_netcdf_var(ncid,'powaal',powtra(1,1,1,ipowaal),ks,0)
      CALL write_netcdf_var(ncid,'powaph',powtra(1,1,1,ipowaph),ks,0)
      CALL write_netcdf_var(ncid,'powaox',powtra(1,1,1,ipowaox),ks,0)
      CALL write_netcdf_var(ncid,'pown2',powtra(1,1,1,ipown2),ks,0)
      CALL write_netcdf_var(ncid,'powno3',powtra(1,1,1,ipowno3),ks,0)
      CALL write_netcdf_var(ncid,'powasi',powtra(1,1,1,ipowasi),ks,0)
#ifdef DIFFAT
      CALL write_netcdf_var(ncid,'suppco2',suppco2(1,1),1,0)
      CALL write_netcdf_var(ncid,'atmco2',atm(1,1,iatmco2),1,0)
      CALL write_netcdf_var(ncid,'atmo2',atm(1,1,iatmo2),1,0)
      CALL write_netcdf_var(ncid,'atmn2',atm(1,1,iatmn2),1,0)
#ifdef __c_isotopes
      CALL write_netcdf_var(ncid,'atmc13',atm(1,1,iatmc13),1,0)
      CALL write_netcdf_var(ncid,'atmc14',atm(1,1,iatmc14),1,0)
#endif
#endif


      IF(mnproc==1) THEN
        ncstat = NF_CLOSE(ncid)
        IF ( ncstat .NE. NF_NOERR ) THEN
          call xchalt('(AUFW: netCDF200)')
                 stop '(AUFW: netCDF200)'
        ENDIF
      ENDIF


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'End of AUFW_BGC'
      WRITE(io_stdo_bgc,*) '***************'
      ENDIF

      RETURN
      END
