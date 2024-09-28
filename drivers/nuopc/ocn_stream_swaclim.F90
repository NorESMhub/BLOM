module ocn_stream_swaclim
#ifdef HAMOCC

   !-----------------------------------------------------------------------
   ! Contains methods for reading in swaclim forcing file (HAMOCC only)
   !-----------------------------------------------------------------------

   use ESMF              , only : ESMF_Clock, ESMF_Mesh, ESMF_Time, ESMF_ClockGet, ESMF_TimeGet
   use ESMF              , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
   use ESMF              , only : ESMF_Finalize, ESMF_LogFoundError
   use nuopc_shr_methods , only : chkerr
   use dshr_strdata_mod  , only : shr_strdata_type
   use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use shr_sys_mod       , only : shr_sys_abort
   use mo_intfcblom      , only : omask
   use mod_forcing       , only : swaclim_forcing
   use mod_xc

   implicit none
   private

   public :: ocn_stream_swaclim_init      ! position datasets for dynamic swaclim
   public :: ocn_stream_swaclim_interp    ! interpolates between two years of swaclim file data

   type(shr_strdata_type) :: sdat_swaclim ! input data stream
   character(len=CS)      :: stream_swaclim_varname='swaclim'

   integer, parameter :: nfiles_max = 100  ! maximum number of stream files
   character(len=*), parameter :: sourcefile = &
        __FILE__

!==============================================================================
contains
!==============================================================================

   subroutine ocn_stream_swaclim_init(model_mesh, model_clock, rc)
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
      integer                        :: nu_nml                    ! unit for namelist file
      integer                        :: nml_error                 ! namelist i/o error flag
      integer                        :: ierr                      ! error status
      character(len=CL)              :: filein                    ! ocn namelist file
      integer                        :: stream_swaclim_year_first ! first year in stream to use
      integer                        :: stream_swaclim_year_last  ! last year in stream to use
      integer                        :: stream_swaclim_year_align ! align stream_year_first swaclim with first model year
      character(len=CL)              :: stream_swaclim_mesh_filename
      character(len=CL)              :: stream_swaclim_data_filename(nfiles_max)
      character(len=CL), allocatable :: stream_filenames(:)
      character(len=4)               :: tmpchar
      integer                        :: nfiles
      integer                        :: nf
      integer                        :: errstat
      character(*), parameter :: subName = "('ocn_stream_swaclim_init')"
      !-----------------------------------------------------------------------

      namelist /stream_swaclim/           &
           stream_swaclim_data_filename,  &
           stream_swaclim_mesh_filename,  &
           stream_swaclim_varname,        &
           stream_swaclim_year_first,     &
           stream_swaclim_year_last,      &
           stream_swaclim_year_align

      rc = ESMF_SUCCESS

      ! Default values for namelist
      stream_swaclim_data_filename(:) = ' '
      stream_swaclim_mesh_filename = ' '
      stream_swaclim_varname = ' '
      stream_swaclim_year_first    = 1 ! first year in stream to use
      stream_swaclim_year_last     = 1 ! last  year in stream to use
      stream_swaclim_year_align    = 1 ! align stream_swaclim_year_first with this model year

      ! Read stream_swaclim namelist
      if (mnproc == 1) then
         filein = "ocn_in"
         open( newunit=nu_nml, file=trim(filein), status='old', iostat=nml_error )
         if (nml_error /= 0) then
            call shr_sys_abort(subName//': ERROR opening '//trim(filein)//errMsg(sourcefile, __LINE__))
         end if
         call shr_nl_find_group_name(nu_nml, 'stream_swaclim', status=nml_error)
         if (nml_error == 0) then
            read(nu_nml, nml=stream_swaclim, iostat=nml_error)
            if (nml_error /= 0) then
               call shr_sys_abort(' ERROR reading stream_swaclim namelist: '//errMsg(sourcefile, __LINE__))
            end if
         else
            call shr_sys_abort(' ERROR finding stream_swaclim namelist: '//errMsg(sourcefile, __LINE__))
         end if
         close(nu_nml)
      endif

      do nf = 1,nfiles_max
         call xcbcst(stream_swaclim_data_filename(nf))
      end do
      call xcbcst(stream_swaclim_mesh_filename)
      call xcbcst(stream_swaclim_varname)
      call xcbcst(stream_swaclim_year_first)
      call xcbcst(stream_swaclim_year_last)
      call xcbcst(stream_swaclim_year_align)

      ! Determine the actual number and filenames of stream files that will be used
      nfiles = 0
      do nf = 1,nfiles_max
         if (stream_swaclim_data_filename(nf) == '') then
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
         stream_filenames(nf) = trim(stream_swaclim_data_filename(nf))
      end do

      ! Write out info
      if (mnproc == 1) then
         write(lp,'(a)'   ) ' '
         write(lp,'(a,i8)')  'stream swaclim settings:'
         do nf = 1,nfiles
            write(tmpchar,'(i0)') nf
            write(lp,'(a,a)' )  '  stream_swaclim_data_filename('//trim(tmpchar)//') = ',trim(stream_filenames(nf))
         end do
         write(lp,'(a,a)' )  '  stream_swaclim_mesh_filename = ',trim(stream_swaclim_mesh_filename)
         write(lp,'(a,a,a)') '  stream_swaclim_varnames      = ',trim(stream_swaclim_varname)
         write(lp,'(a,i8)')  '  stream_swaclim_year_first    = ',stream_swaclim_year_first
         write(lp,'(a,i8)')  '  stream_swaclim_year_last     = ',stream_swaclim_year_last
         write(lp,'(a,i8)')  '  stream_swaclim_year_align    = ',stream_swaclim_year_align
         write(lp,'(a)'   )  ' '
      endif

      ! Initialize the cdeps data type sdat_swaclim
      call shr_strdata_init_from_inline(sdat_swaclim,                &
           my_task             = mnproc,                             &
           logunit             = lp,                                 &
           compname            = 'OCN',                              &
           model_clock         = model_clock,                        &
           model_mesh          = model_mesh,                         &
           stream_meshfile     = trim(stream_swaclim_mesh_filename), &
           stream_filenames    = stream_filenames,                   &
           stream_yearFirst    = stream_swaclim_year_first,          &
           stream_yearLast     = stream_swaclim_year_last,           &
           stream_yearAlign    = stream_swaclim_year_align,          &
           stream_fldlistFile  = (/stream_swaclim_varname/),         &
           stream_fldlistModel = (/stream_swaclim_varname/),         &
           stream_lev_dimname  = 'null',                             &
           stream_mapalgo      = 'bilinear',                         &
           stream_offset       = 0,                                  &
           stream_taxmode      = 'cycle',                            &
           stream_dtlimit      = 1.0e30_r8,                          &
           stream_tintalgo     = 'linear',                           &
           stream_name         = 'Climatological SWA',               &
           rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

   end subroutine ocn_stream_swaclim_init

   !================================================================
   subroutine ocn_stream_swaclim_interp(model_clock, rc)

      use dshr_strdata_mod , only : shr_strdata_advance
      use dshr_methods_mod , only : dshr_fldbun_getfldptr
      use mod_checksum     , only : csdiag, chksummsk

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

      call shr_strdata_advance(sdat_swaclim, ymd=mcdate, tod=sec, logunit=lp, istr='swaclim', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Get pointer for stream data that is time and spatially interpolated to model time and grid
      call dshr_fldbun_getFldPtr(sdat_swaclim%pstrm(1)%fldbun_model, stream_swaclim_varname, fldptr1=dataptr1, rc=rc)
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
               swaclim_forcing(i,j) = mval
            elseif (cplmsk(i,j) == 0) then
               swaclim_forcing(i,j) = fval
            else
               n = (j - 1)*ii + i
               swaclim_forcing(i,j) = dataptr1(n)
            end if
         end do
      end do

      if (csdiag) then
         if (mnproc == 1) then
            write(lp,*) 'ocn_forcing_swaclim_interp:'
         end if
         call chksummsk(swaclim_forcing(1-nbdy,1-nbdy),ip,1,'swaclim_forcing')
      end if

      ! set flux to zero over land
      do j=1,jdm
         do i=1,idm
            if (omask(i,j) < 0.5) then
               swaclim_forcing(i,j) = 0.0
            end if
         enddo
      enddo

   end subroutine ocn_stream_swaclim_interp

#endif
end module ocn_stream_swaclim
