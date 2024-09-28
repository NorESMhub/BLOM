module ocn_stream_rivin

#ifdef HAMOCC
  !*************************************************************************************************
  ! Routine for reading riverine nutrient data
  !
  ! Riverine carbon and nutrient input is activated through a logical switch 'do_rivinpt' read
  ! from HAMOCC's bgcnml namelist. When coupled to NorESM, this is achieved by setting
  ! BLOM_RIVER_NUTRIENTS to TRUE in env_run.xml.
  !
  ! The model attempts to read nutrient fluxes from a NetCDF file
  ! derived from the GNEWS 2000 data base, which is specified through the
  ! namelist. The nutrient fluxes are then interpolated to the model ocean grid.
  !
  ! The nutrient discharge is distributed on the ocean grid in manner that is
  ! consistent with how model distributes its freshwater runoff.
  !
  ! Since only alkalinity is available from measurements, DIC is updated using
  ! the assumtions that a_t=a_c+a_n and DIC=a_c (a_t: total alkalinity,
  ! a_c: carbonate alkalinity, a_n: contribution of nutrients to a_t).
  !*************************************************************************************************

   use ESMF              , only: ESMF_Finalize, ESMF_LogFoundError
   use ESMF              , only: ESMF_Mesh, ESMF_MeshCreate, ESMF_MeshGet, ESMF_FieldCreate
   use ESMF              , only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
   use ESMF              , only: ESMF_FieldSMMStore, ESMF_FieldRedistStore
   use ESMF              , only: ESMF_FieldRegrid, ESMF_RouteHandle
   use ESMF              , only: ESMF_DistGrid, ESMF_DistGridGet
   use ESMF              , only: ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
   use ESMF              , only: ESMF_FILEFORMAT_ESMFMESH, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
   use ESMF              , only: ESMF_REGION_TOTAL, ESMF_TERMORDER_SRCSEQ
   use pio               , only: pio_openfile, pio_closefile, pio_nowrite, file_desc_t
   use shr_pio_mod       , only: shr_pio_getiosys, shr_pio_getiotype
   use shr_kind_mod      , only: r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_log_mod       , only: errMsg => shr_log_errMsg
   use shr_sys_mod       , only: shr_sys_abort
   use shr_nl_mod        , only: shr_nl_find_group_name
   use nuopc_shr_methods , only: chkerr
   use mod_blom_pio      , only: pio_iosystem, pio_iotype, blom_pio_getdata
   use mod_forcing       , only: rivflx_forcing
   use mod_output_forcing, only: output_forcing
   use mo_param1_bgc     , only: nriv,irdin,irdip,irsi,iralk,iriron,irdoc,irdet
   use mo_intfcblom      , only : omask
   use mod_xc

   implicit none
   private

   public :: ocn_stream_rivin_init      ! position datasets for dynamic rivin

   character(len=*), parameter :: sourcefile = &
        __FILE__

!==============================================================================
contains
!==============================================================================

   subroutine ocn_stream_rivin_init(mesh_model, rc)

      ! input/output variables
      type(ESMF_Mesh) , intent(in)  :: mesh_model
      integer         , intent(out) :: rc

      ! local variables
      integer                :: nu_nml                ! unit for namelist file
      integer                :: nml_error             ! namelist i/o error flag
      integer                :: ierr                  ! error status
      type(file_desc_t)      :: pioid
      character(len=CL)      :: filein
      character(len=CL)      :: stream_rivin_mesh_filename
      character(len=CL)      :: stream_rivin_data_filename
      character(len=CL)      :: stream_rivin_map_filename
      character(len=CL)      :: stream_rivin_varname
      integer                :: errstat
      type(ESMF_MESH)        :: mesh_file
      type(ESMF_DistGrid)    :: distgrid
      type(ESMF_Field)       :: field_file
      type(ESMF_Field)       :: field_model
      type(ESMF_RouteHandle) :: rhandle_file2model
      integer                :: srcMaskValue = 0          ! ignore source points where the mesh mask is 0
      integer                :: dstMaskValue = -987987    ! don't ingore any destination points
      integer                :: srcTermProcessing_Value = 0
      real(r8), pointer      :: dataptr_file(:)
      real(r8), pointer      :: dataptr_model(:)
      real(r8), allocatable  :: data_file(:)
      integer                :: nf,i,j,n
      integer                :: lsize
      integer                :: jjcpl
      integer                :: index
      character(len=CL)      :: output_filename
      character(len=3)       :: fldnames(7)
      integer                :: indices(7)
      integer                :: rcode
      logical                :: checkflag = .false.
      real(r8), parameter    :: mval = -1.e12_r8
      real(r8), parameter    :: fval = -1.e13_r8
      character(*), parameter :: subName = "('ocn_stream_rivin_init')"
      !-----------------------------------------------------------------------

      namelist /stream_rivin/           &
           stream_rivin_data_filename,  &
           stream_rivin_mesh_filename,  &
           stream_rivin_map_filename

      rc = ESMF_SUCCESS

      fldnames = (/'DIN', 'DIP', 'DSi', 'DIC',   'Fe ', 'DOC', 'DET'/)
      indices  = (/irdin, irdip,  irsi,  iralk, iriron, irdoc, irdet/)

      ! Default values for namelist
      stream_rivin_data_filename(:) = ' '
      stream_rivin_mesh_filename = ' '
      stream_rivin_map_filename = ' '

      ! Read stream_rivin namelist
      if (mnproc == 1) then
         filein = "ocn_in"
         open( newunit=nu_nml, file=trim(filein), status='old', iostat=nml_error )
         if (nml_error /= 0) then
            call shr_sys_abort(subName//': ERROR opening '//trim(filein)//errMsg(sourcefile, __LINE__))
         end if
         call shr_nl_find_group_name(nu_nml, 'stream_rivin', status=nml_error)
         if (nml_error == 0) then
            read(nu_nml, nml=stream_rivin, iostat=nml_error)
            if (nml_error /= 0) then
               call shr_sys_abort(' ERROR reading stream_rivin namelist: '//errMsg(sourcefile, __LINE__))
            end if
         else
            call shr_sys_abort(' ERROR finding stream_rivin namelist: '//errMsg(sourcefile, __LINE__))
         end if
         close(nu_nml)
      endif

      call xcbcst(stream_rivin_data_filename)
      call xcbcst(stream_rivin_mesh_filename)
      call xcbcst(stream_rivin_map_filename)

      ! Write out info
      if (mnproc == 1) then
         write(lp,'(a)'   ) ' '
         write(lp,'(a)'   )  'stream rivin settings:'
         write(lp,'(a,a)' )  '  stream_rivin_data_filename = ',trim(stream_rivin_data_filename)
         write(lp,'(a,a)' )  '  stream_rivin_mesh_filename = ',trim(stream_rivin_mesh_filename)
         write(lp,'(a,a)' )  '  stream_rivin_map_filename  = ',trim(stream_rivin_map_filename)
         write(lp,'(a)'   )  ' '
      endif

      ! allocate field to hold rivin fields
      allocate (rivflx_forcing(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nriv),stat=errstat)
      if(errstat /= 0) stop 'not enough memory rivflx allocation'
      rivflx_forcing(:,:,:) = 0.0

      ! read riverine nutrient fluxes from file
      if (mnproc.eq.1) then
         write(lp,*) ''
         write(lp,'(a)') 'ini_read_rivin: read riverine nutrients from ',trim(stream_rivin_data_filename)
      endif

      ! Read in the mesh
      mesh_file = ESMF_MeshCreate(filename=stream_rivin_mesh_filename, fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
      if (ChkErr(rc, __LINE__, sourcefile)) return

      ! Creat fields used in remapping
      field_file  = ESMF_FieldCreate(mesh_file, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
      if (ChkErr(rc,__LINE__,sourcefile)) return
      field_model = ESMF_FieldCreate(mesh_model, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
      if (ChkErr(rc,__LINE__,sourcefile)) return

      ! Compute a route handle between the model mesh and data mesh
      if (stream_rivin_map_filename == 'unset') then
         ! create a redist route handle
         call ESMF_FieldRedistStore(field_file, field_model, routehandle=rhandle_file2model, &
              ignoreUnmatchedIndices = .true., rc=rc)
         if (chkerr(rc,__LINE__,sourcefile)) return
      else
         ! create a route handle from the input mapping file
         call ESMF_FieldSMMStore(field_file, field_model, stream_rivin_map_filename, routehandle=rhandle_file2model, &
              ignoreUnmatchedIndices=.true., srcTermProcessing=srcTermProcessing_Value, rc=rc)
         if (chkerr(rc,__LINE__,sourcefile)) return
      end if

      ! Determine local size of input data from file
      call ESMF_MeshGet(mesh_file, elementdistGrid=distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,sourcefile)) return
      call ESMF_DistGridGet(distgrid, localDe=0, elementCount=lsize, rc=rc)
      if (ChkErr(rc,__LINE__,sourcefile)) return

      ! Allocate data_file
      allocate(data_file(lsize), stat=errstat)
      if (errstat /= 0) stop 'not enough memory for data_fileallocation'

      ! Open input data file
      rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(stream_rivin_data_filename), pio_nowrite)

      do nf = 1,size(fldnames)

         ! Read the data (fill in the data_file array)
         call blom_pio_getdata(pioid, fldnames(nf), mesh_file, data_file, rc=rc)
         if (chkerr(rc,__LINE__,sourcefile)) return

         ! Get pointer to data_file array and fill it in with data_file valus
         call ESMF_FieldGet(field_file, farrayptr=dataptr_file, rc=rc)
         if (chkerr(rc,__LINE__,sourcefile)) return
         dataptr_file(:) = data_file(:)

         ! Regrid the data
         call ESMF_FieldRegrid(field_file, field_model, routehandle=rhandle_file2model, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
         if (chkerr(rc,__LINE__,sourcefile)) return

         ! Get pointer to mapped field
         call ESMF_FieldGet(field_model, farrayptr=dataptr_model, rc=rc)
         if (chkerr(rc,__LINE__,sourcefile)) return

         ! Fill in the blom data structures
         if (nreg == 2 .and. nproc == jpr) then
            jjcpl = jj - 1
         else
            jjcpl = jj
         endif

         index = indices(nf)
         do j = 1, jjcpl
            do i = 1, ii
               if (ip(i,j) == 0) then
                  rivflx_forcing(i,j,index) = mval
               elseif (cplmsk(i,j) == 0) then
                  rivflx_forcing(i,j,index) = fval
               else
                  n = (j - 1)*ii + i
                  rivflx_forcing(i,j,index) = dataptr_model(n)
               end if
               ! set flux to zero over land
               if (omask(i,j) < 0.5) then
                  rivflx_forcing(i,j,index) = 0.0
               end if
            end do
         end do

         output_filename = 'frivin_stream_'//trim(fldnames(nf))//'.nc'
         call output_forcing(trim(output_filename), trim(fldnames(nf)), rivflx_forcing(:,:,index))
      end do

      call pio_closefile(pioid)

   end subroutine ocn_stream_rivin_init

#endif
end module ocn_stream_rivin
