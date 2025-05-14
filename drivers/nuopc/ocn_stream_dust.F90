module ocn_stream_dust
#ifdef HAMOCC
   !-----------------------------------------------------------------------
   ! Contains methods for reading in dust forcing file
   !-----------------------------------------------------------------------
   !
   use ESMF              , only : ESMF_Clock, ESMF_Mesh, ESMF_Time, ESMF_ClockGet, ESMF_TimeGet
   use ESMF              , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
   use ESMF              , only : ESMF_Finalize, ESMF_LogFoundError
   use nuopc_shr_methods , only : chkerr
   use dshr_strdata_mod  , only : shr_strdata_type
   use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use shr_sys_mod       , only : shr_sys_abort
   use mod_forcing       , only : dust_stream
   use mo_intfcblom      , only : omask
   use mod_config        , only : inst_suffix
   use mod_xc

   implicit none
   private

   public :: ocn_stream_dust_init      ! position datasets for dynamic dust
   public :: ocn_stream_dust_interp    ! interpolates between two years of dust file data

   type(shr_strdata_type) :: sdat_dust ! input data stream
   character(len=CS)      :: stream_dust_varname='dustdep'

   integer, parameter :: nfiles_max = 100  ! maximum number of stream files
   character(len=*), parameter :: sourcefile = &
        __FILE__

!==============================================================================
contains
!==============================================================================

   subroutine ocn_stream_dust_init(model_mesh, model_clock, rc)
      !
      ! Initialize data stream information.

      ! Uses
      use shr_nl_mod       , only: shr_nl_find_group_name
      use dshr_strdata_mod , only: shr_strdata_init_from_inline

      ! input/output variables
      type(ESMF_CLock), intent(in)  :: model_clock
      type(ESMF_Mesh) , intent(in)  :: model_mesh
      integer         , intent(out) :: rc

      ! local variables
      integer                        :: nu_nml                ! unit for namelist file
      integer                        :: nml_error             ! namelist i/o error flag
      integer                        :: ierr                  ! error status
      character(len=CL)              :: filein                ! ocn namelist file
      integer                        :: stream_dust_year_first ! first year in stream to use
      integer                        :: stream_dust_year_last  ! last year in stream to use
      integer                        :: stream_dust_year_align ! align stream_year_firstdust with
      character(len=CL)              :: stream_dust_mesh_filename
      character(len=CL)              :: stream_dust_data_filename(nfiles_max)
      character(len=CL), allocatable :: stream_filenames(:)
      character(len=4)               :: tmpchar
      integer                        :: nfiles
      integer                        :: nf
      integer                        :: errstat
      character(*), parameter :: subName = "('ocn_stream_dust_init')"
      !-----------------------------------------------------------------------

      namelist /stream_dust/           &
           stream_dust_data_filename,  &
           stream_dust_mesh_filename,  &
           stream_dust_varname,        &
           stream_dust_year_first,     &
           stream_dust_year_last,      &
           stream_dust_year_align

      rc = ESMF_SUCCESS

      ! Default values for namelist
      stream_dust_data_filename(:) = ' '
      stream_dust_mesh_filename = ' '
      stream_dust_year_first    = 1 ! first year in stream to use
      stream_dust_year_last     = 1 ! last  year in stream to use
      stream_dust_year_align    = 1 ! align stream_dust_year_first with this model year

      ! Read stream_dust namelist
      if (mnproc == 1) then
         filein = trim("ocn_in"//inst_suffix)
         open( newunit=nu_nml, file=trim(filein), status='old', iostat=nml_error )
         if (nml_error /= 0) then
            call shr_sys_abort(subName//': ERROR opening '//trim(filein)//errMsg(sourcefile, __LINE__))
         end if
         call shr_nl_find_group_name(nu_nml, 'stream_dust', status=nml_error)
         if (nml_error == 0) then
            read(nu_nml, nml=stream_dust, iostat=nml_error)
            if (nml_error /= 0) then
               call shr_sys_abort(' ERROR reading stream_dust namelist: '//errMsg(sourcefile, __LINE__))
            end if
         else
            call shr_sys_abort(' ERROR finding stream_dust namelist: '//errMsg(sourcefile, __LINE__))
         end if
         close(nu_nml)
      endif

      do nf = 1,nfiles_max
         call xcbcst(stream_dust_data_filename(nf))
      end do
      call xcbcst(stream_dust_mesh_filename)
      call xcbcst(stream_dust_varname)
      call xcbcst(stream_dust_year_first)
      call xcbcst(stream_dust_year_last)
      call xcbcst(stream_dust_year_align)

      ! Determine the actual number and filenames of stream files that will be used
      nfiles = 0
      do nf = 1,nfiles_max
         if (stream_dust_data_filename(nf) == '') then
            exit
         else
            nfiles = nfiles + 1
         end if
      end do
      allocate(stream_filenames(nfiles), stat=errstat)
      if (errstat /= 0) then
         write(lp,*) 'Failed to allocate stream_filenames'
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if
      end if
      do nf = 1,nfiles
         stream_filenames(nf) = trim(stream_dust_data_filename(nf))
      end do

      ! Write out info
      if (mnproc == 1) then
         write(lp,'(a)'   ) ' '
         write(lp,'(a,i8)')  'stream dust settings:'
         do nf = 1,nfiles
            write(tmpchar,'(i0)') nf
            write(lp,'(a,a)' )  '  stream_dust_data_filename('//trim(tmpchar)//') = ',trim(stream_filenames(nf))
         end do
         write(lp,'(a,a)' )  '  stream_dust_mesh_filename = ',trim(stream_dust_mesh_filename)
         write(lp,'(a,a,a)') '  stream_dust_varname       = ',trim(stream_dust_varname)
         write(lp,'(a,i8)')  '  stream_dust_year_first    = ',stream_dust_year_first
         write(lp,'(a,i8)')  '  stream_dust_year_last     = ',stream_dust_year_last
         write(lp,'(a,i8)')  '  stream_dust_year_align    = ',stream_dust_year_align
         write(lp,'(a)'   )  ' '
      endif

      ! Initialize the cdeps data type sdat_dust
      ! tintalgo can be one of [lower,upper,nearest,linear,coszen]
      call shr_strdata_init_from_inline(sdat_dust,                 &
           my_task             = mnproc,                           &
           logunit             = lp,                               &
           compname            = 'OCN',                            &
           model_clock         = model_clock,                      &
           model_mesh          = model_mesh,                       &
           stream_meshfile     = trim(stream_dust_mesh_filename),  &
           stream_filenames    = stream_filenames,                 &
           stream_yearFirst    = stream_dust_year_first,           &
           stream_yearLast     = stream_dust_year_last,            &
           stream_yearAlign    = stream_dust_year_align,           &
           stream_fldlistFile  = (/stream_dust_varname/),          &
           stream_fldlistModel = (/stream_dust_varname/),          &
           stream_lev_dimname  = 'null',                           &
           stream_mapalgo      = 'bilinear',                       &
           stream_offset       = 0,                                &
           stream_taxmode      = 'cycle',                          &
           stream_dtlimit      = 1.0e30_r8,                        &
           stream_tintalgo     = 'linear',                         &
           stream_name         = 'Fe Dust Deposition',             &
           rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! allocate field to hold dust fields
      dust_stream(:,:) = 0.0

   end subroutine ocn_stream_dust_init

   !================================================================
   subroutine ocn_stream_dust_interp(model_clock, rc)

      use dshr_strdata_mod   , only : shr_strdata_advance
      use dshr_methods_mod   , only : dshr_fldbun_getfldptr
      use mod_checksum       , only : csdiag, chksummsk

      ! input/output variables
      type(ESMF_Clock), intent(in)  :: model_clock
      integer         , intent(out) :: rc

      ! local variables
      type(ESMF_Time)     :: date
      integer             :: i,j,n
      integer             :: jjcpl
      integer             :: year    ! year (0, ...) for nstep+1
      integer             :: mon     ! month (1, ..., 12) for nstep+1
      integer             :: day     ! day of month (1, ..., 31) for nstep+1
      integer             :: sec     ! seconds into current date for nstep+1
      integer             :: mcdate  ! Current model date (yyyymmdd)
      real(r8), pointer   :: dataptr1(:)
      real(r8), parameter :: mval = -1.e12_r8
      real(r8), parameter :: fval = -1.e13_r8
      integer             :: seconds_per_month = 2592000 ! assuming 30 days per month
      !-----------------------------------------------------------------------

      ! Advance sdat stream
      call ESMF_ClockGet( model_clock, currTime=date, rc=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
      call ESMF_TimeGet(date, yy=year, mm=mon, dd=day, s=sec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
      mcdate = year*10000 + mon*100 + day

      call shr_strdata_advance(sdat_dust, ymd=mcdate, tod=sec, logunit=lp, istr='dust', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Get pointer for stream data that is time and spatially interpolated to model time and grid
      call dshr_fldbun_getFldPtr(sdat_dust%pstrm(1)%fldbun_model, stream_dust_varname, fldptr1=dataptr1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Set the j-extent of the local ocean domain to be exchanged. Needed because of
      ! duplication of the last global domain row when using a tripolar grid.
      ! dimensions.F is created at build time by the bash script blom_dimensions and contains

      if (nreg == 2 .and. nproc == jpr) then
         jjcpl = jj - 1
      else
         jjcpl = jj
      endif

      do j = 1, jjcpl
         do i = 1, ii
            if (ip(i,j) == 0) then
               dust_stream(i,j) = mval
            elseif (cplmsk(i,j) == 0) then
               dust_stream(i,j) = fval
            else
               n = (j - 1)*ii + i
               ! need dust_stream in units of kg/m2/month not kg/m2/s which is in the file
               dust_stream(i,j) = dataptr1(n)*seconds_per_month
            end if
            ! set flux to zero over land
            if (omask(i,j) < 0.5) then
               dust_stream(i,j) = 0.0
            end if
         end do
      end do

      if (csdiag) then
         if (mnproc == 1) then
            write(lp,*) 'ocn_stream_dust_interp:'
         end if
         call chksummsk(dust_stream,ip,1,'dust_stream')
      end if

   end subroutine ocn_stream_dust_interp

#endif
end module ocn_stream_dust
