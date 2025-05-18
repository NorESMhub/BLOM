module mod_io_input

   use ESMF              , only : ESMF_Clock, ESMF_Mesh, ESMF_Time, ESMF_ClockGet, ESMF_TimeGet
   use ESMF              , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
   use ESMF              , only : ESMF_Finalize, ESMF_LogFoundError
   use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleGet
   use ESMF              , only : ESMF_Field, ESMF_RouteHandle, ESMF_DistGrid
   use ESMF              , only : ESMF_FieldStatus_Flag, ESMF_TYPEKIND_R8, ESMF_FieldCreate
   use ESMF              , only : ESMF_FieldRegridStore, ESMF_FieldRegrid, ESMF_FieldGet
   use ESMF              , only : ESMF_MeshGet, ESMF_DistGridGet, ESMF_FieldStatus_Flag, ESMF_FIELDSTATUS_COMPLETE
   use ESMF              , only : ESMF_MESHLOC_ELEMENT, ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_DSTAREA
   use ESMF              , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_TERMORDER_SRCSEQ, ESMF_REGION_TOTAL
   use ESMF              , only : ESMF_REGRIDMETHOD_BILINEAR,ESMF_POLEMETHOD_ALLAVG
   use ESMF              , only : ESMF_MAXSTR, ESMF_ITEMORDER_ADDORDER, ESMF_KIND_R8, ESMF_REGION_EMPTY
   use ESMF              , only : ESMF_FieldDestroy
   use ESMF              , only : ESMF_DYNAMICMASK, ESMF_DynamicMaskSetR8R8R8, ESMF_DYNAMICMASKELEMENTR8R8R8
   use ESMF              , only : operator(/=), operator(==)
   use nuopc_shr_methods , only : chkerr
   use shr_kind_mod      , only : r8 => shr_kind_r8, r4 => shr_kind_r4, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_log_mod       , only : shr_log_error
   use shr_const_mod     , only : r8fill => SHR_CONST_SPVAL
   use shr_infnan_mod    , only : shr_infnan_isnan
   use shr_pio_mod       , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
   use pio               , only : file_desc_t, iosystem_desc_t, io_desc_t, var_desc_t
   use pio               , only : pio_openfile, pio_closefile, pio_nowrite
   use pio               , only : pio_seterrorhandling, pio_initdecomp, pio_freedecomp
   use pio               , only : pio_inquire, pio_inq_varid, pio_inq_varndims, pio_inq_vardimid
   use pio               , only : pio_inq_dimlen, pio_inq_vartype, pio_inq_dimname, pio_inq_dimid
   use pio               , only : pio_double, pio_real, pio_int, pio_offset_kind, pio_get_var
   use pio               , only : pio_read_darray, pio_setframe, pio_fill_double, pio_get_att, pio_inq_att
   use pio               , only : PIO_BCAST_ERROR, PIO_RETURN_ERROR, PIO_NOERR, PIO_INTERNAL_ERROR
   use mod_xc            , only : mnproc, lp

   implicit none
   private

   public :: read_map_input_data
   public :: field_getfldptr

   private :: set_iodesc

   type(ESMF_DynamicMask)         :: dynamicOcnMask
   type(iosystem_desc_t), pointer :: pio_subsystem => null()     ! pio info
   integer                        :: io_type                     ! pio info
   integer                        :: io_format                   ! pio info

   character(len=*), parameter :: u_FILE_u = &
      __FILE__

!===============================================================================
contains
!===============================================================================

   subroutine read_map_input_data(mesh_input, filename, fldlist, nlev, regrid_method, fldbun_blom, rc)

      ! Read and map the input data

      ! input/output variables
      type(ESMF_Mesh)             , intent(in)    :: mesh_input
      character(len=*)            , intent(in)    :: filename
      character(len=*)            , intent(in)    :: fldlist(:)
      integer                     , intent(in)    :: nlev
      character(len=*)            , intent(in)    :: regrid_method ! conserve or bilinear for now
      type(ESMF_FieldBundle)      , intent(inout) :: fldbun_blom   ! input data interpolated to model grid
      integer                     , intent(out)   :: rc

      ! local variables
      real(r8)                :: dynamicSrcMaskValue
      type(ESMF_Field)        :: field_data
      type(ESMF_Field)        :: field_dst
      type(file_desc_t)       :: pioid
      type(var_desc_t)        :: varid
      integer                 :: nt
      integer                 :: nf
      integer                 :: rCode
      real(r4)                :: fillvalue_r4
      real(r8)                :: fillvalue_r8
      logical                 :: handlefill = .false.
      integer                 :: old_error_handle
      real(r8), pointer       :: dataptr(:)
      real(r8), pointer       :: dataptr1d(:)        ! field bundle data
      real(r8), pointer       :: dataptr2d(:,:)      ! field bundle data
      real(r4), allocatable   :: data_real1d(:)      ! input data
      real(r4), allocatable   :: data_real2d(:,:)    ! input data
      real(r8), allocatable   :: data_dbl1d(:)       ! input data
      real(r8), allocatable   :: data_dbl2d(:,:)     ! input data
      integer                 :: lsize, n
      integer                 :: spatialDim, numOwnedElements
      integer                 :: pio_iovartype
      integer                 :: i, lev
      logical                 :: checkflag = .false.
      type(io_desc_t)         :: pio_iodesc             ! stream pio descriptor
      type(ESMF_RouteHandle)  :: routehandle
      !integer                 :: srcMaskValue = 0
      integer                 :: srcMaskValue = -987987
      integer                 :: dstMaskValue = -987987 ! spval for RH mask values
      integer                 :: srcTermProcessing_Value = 0
      character(*), parameter :: subname = '(read_map_input_data) '
      character(*), parameter :: F00   = "('(read_map_input_data) ',8a)"
      character(*), parameter :: F02   = "('(read_map_input_data) ',2a,i8)"
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! nullify local pointers
      nullify(dataptr)
      nullify(dataptr1d)
      nullify(dataptr2d)

      ! Initialize PIO
      ! TODO - print io_type and io_format to make sure they make sense (numbers between 1 and 4)
      pio_subsystem => shr_pio_getiosys('OCN')
      io_type       =  shr_pio_getiotype('OCN')
      io_format     =  shr_pio_getioformat('OCN')

      ! Open file
      if (mnproc == 1) then
         write(lp,F00) 'opening   : ',trim(filename)
      endif
      rcode = pio_openfile(pio_subsystem, pioid, io_type, trim(filename), pio_nowrite)

      ! ******************************************************************************
      ! Determine the pio io descriptor
      ! ******************************************************************************

      if (mnproc == 1) then
         write(lp,F00) 'setting pio descriptor : ',trim(filename)
      end if
      call set_iodesc(mesh_input, fldlist(1), nlev, pioid, pio_iodesc, rc=rc)

      ! ******************************************************************************
      ! Create a data field that can be used for reading in all input data
      ! ******************************************************************************
      if (nlev > 1) then
         field_data = ESMF_FieldCreate(mesh_input, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
              ungriddedLbound=(/1/), ungriddedUbound=(/nlev/), gridToFieldMap=(/2/), rc=rc)
         call field_getfldptr(field_data, fldptr2=dataptr2d, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
      else
         field_data = ESMF_FieldCreate(mesh_input, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
         call field_getfldptr(field_data, fldptr1=dataptr1d, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
      end if

      ! ******************************************************************************
      ! Read in the  data for field names in fldlist
      ! Note that the dimensions in the field (nlev,i) whereas
      ! the dimension of the input data are (i,nlev)
      ! ******************************************************************************

      if (mnproc == 1) then
         write(lp,F02) 'reading file '//trim(filename)
      endif

      do nf = 1,size(fldlist)

         ! determine type of the variable
         rcode = pio_inq_varid(pioid, trim(fldlist(nf)), varid)
         rcode = pio_inq_vartype(pioid, varid, pio_iovartype)

         ! allocate memory for input read
         if (nlev > 1) then
            lsize = size(dataptr2d, dim=2)
            if (pio_iovartype == PIO_REAL .and. .not. allocated(data_real2d)) then
               allocate(data_real2d(lsize, nlev))
            else if (pio_iovartype == PIO_DOUBLE .and. .not. allocated(data_dbl2d)) then
               allocate(data_dbl2d(lsize, nlev))
            endif
         else
            lsize = size(dataptr1d)
            if (pio_iovartype == PIO_REAL .and. .not. allocated(data_real1d)) then
               allocate(data_real1d(lsize))
            else if (pio_iovartype == PIO_DOUBLE .and. .not. allocated(data_dbl1d)) then
               allocate(data_dbl1d(lsize))
            endif
         end if

         ! determine if will handle fill
         handlefill = .false.
         call PIO_seterrorhandling(pioid, PIO_BCAST_ERROR, old_error_handle)
         if (pio_iovartype == PIO_REAL) then
            rcode = pio_get_att(pioid, varid, "_FillValue", fillvalue_r4)
         else if (pio_iovartype == PIO_DOUBLE) then
            rcode = pio_get_att(pioid, varid, "_FillValue", fillvalue_r8)
         endif
         if(rcode == PIO_NOERR) handlefill=.true.
         call PIO_seterrorhandling(pioid, old_error_handle)

         if (mnproc == 1) then
            write(lp,F02)' reading '// trim(fldlist(nf))//' at time index 1 '
         end if

         ! Set only 1 time index
         nt = 1

         ! read the data
         call pio_setframe(pioid, varid, int(nt,kind=Pio_Offset_Kind))
         if (pio_iovartype == PIO_REAL) then
            ! -----------------------------
            ! pio_iovartype is PIO_REAL
            ! -----------------------------
            if (nlev > 1) then
               call pio_read_darray(pioid, varid, pio_iodesc, data_real2d, rcode)
               if ( rcode /= PIO_NOERR ) then
                  rc = rcode
                  call shr_log_error(' ERROR: reading in variable: '// trim(fldlist(nf)), rc=rc)
                  return
               end if
               if (handlefill) then
                  do lev = 1,nlev
                     do n = 1,size(dataptr2d, dim=2)
                        if (.not. shr_infnan_isnan(data_real2d(n,lev)) .and. data_real2d(n,lev) .ne. fillvalue_r4) then
                           dataptr2d(lev,n) = real(data_real2d(n,lev), kind=r8) ! Note the order of indices
                        else
                           dataptr2d(lev,n) = r8fill
                        endif
                     enddo
                  end do
               else
                  do lev = 1,nlev
                     do n = 1,size(dataptr2d, dim=2)
                        dataptr2d(lev,n) = real(data_real2d(n,lev), kind=r8)
                     end do
                  end do
               end if
            else ! nlev == 1
               call pio_read_darray(pioid, varid, pio_iodesc, data_real1d, rcode)
               if ( rcode /= PIO_NOERR ) then
                  rc = rcode
                  call shr_log_error(' ERROR: reading in variable: '// trim(fldlist(nf)), rc=rc)
                  return
               end if
               if (handlefill) then
                  do n=1,size(dataptr1d)
                     if(.not. shr_infnan_isnan(data_real1d(n)) .and. data_real1d(n) .ne. fillvalue_r4) then
                        dataptr1d(n) = real(data_real1d(n), kind=r8)
                     else
                        dataptr1d(n) = r8fill
                     endif
                  enddo
               else
                  dataptr1d(:) = real(data_real1d(:),kind=r8)
               endif
            end if

         else if (pio_iovartype == PIO_DOUBLE) then
            ! -----------------------------
            ! pio_iovartype is PIO_DOUBLE
            ! -----------------------------
            if (nlev > 1) then
               call pio_read_darray(pioid, varid, pio_iodesc, data_dbl2d, rcode)
               if ( rcode /= PIO_NOERR ) then
                  rc = rcode
                  call shr_log_error(' ERROR: reading in 2d double variable: '// trim(fldlist(nf)), rc=rc)
                  return
               end if
               if (handlefill) then
                  do lev = 1,nlev
                     do n = 1,size(dataptr2d, dim=2)
                        if (.not. shr_infnan_isnan(data_dbl2d(n,lev)) .and. data_dbl2d(n,lev) .ne. fillvalue_r8) then
                           dataptr2d(lev,n) = data_dbl2d(n,lev)
                        else
                           dataptr2d(lev,n) = r8fill
                        endif
                     enddo
                  end do
               else
                  do lev = 1,nlev
                     do n = 1,size(dataptr2d, dim=2)
                        dataptr2d(lev,n) = data_dbl2d(n,lev)
                     end do
                  end do
               end if
            else ! nlev == 1
               call pio_read_darray(pioid, varid, pio_iodesc, data_dbl1d, rcode)
               if ( rcode /= PIO_NOERR ) then
                  rc = rcode
                  call shr_log_error(' ERROR: reading in variable: '// trim(fldlist(nf)), rc=rc)
                  return
               end if
               if (handlefill) then
                  do n = 1,size(dataptr1d)
                     if (.not. shr_infnan_isnan(data_dbl1d(n)) .and. data_dbl1d(n) .ne. fillvalue_r8) then
                        dataptr1d(n) = data_dbl1d(n)
                     else
                        dataptr1d(n) = r8fill
                     end if
                  enddo
               else
                  do n = 1,size(dataptr1d)
                     dataptr1d(n) = data_dbl1d(n)
                  end do
               endif
            end if

         else
            ! -----------------------------
            ! pio_iovartype is not supported
            ! -----------------------------
            call shr_log_error(subName//"ERROR: only double and real types are supported for read", rc=rc)
            return

         end if

         ! ******************************************************************************
         ! Regrid the input data to the blom grid
         ! ******************************************************************************

         ! Get the field from the input field bundle
         call fldbun_getfieldN(fldbun_blom, nf, field_dst, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         ! create route handle to map input data to blom mesh
         if (nf == 1) then
            ! Create a dynamic mask object
            ! The dynamic mask object further holds a pointer to the routine that will be called in order to
            ! handle dynamically masked elements - in this case its DynOcnMaskProc (see below)
            dynamicSrcMaskValue = 1.e30_r8
            call ESMF_DynamicMaskSetR8R8R8(dynamicOcnMask, dynamicMaskRoutine=DynOcnMaskProc, &
                 dynamicSrcMaskValue=dynamicSrcMaskValue,  handleAllElements=.true., rc=rc)
            if (chkerr(rc,__LINE__,u_FILE_u)) return

            if (regrid_method == 'conserve') then
               call ESMF_FieldRegridStore(field_data, field_dst, routehandle=routehandle, &
                    srcMaskValues=(/srcMaskValue/), &
                    dstMaskValues=(/dstMaskValue/), &
                    regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
                    normType=ESMF_NORMTYPE_DSTAREA, &
                    srcTermProcessing=srcTermProcessing_Value, &
                    ignoreDegenerate=.true., &
                    unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
               if (chkerr(rc,__LINE__,u_FILE_u)) return
            else if (regrid_method == 'bilinear') then
               call ESMF_FieldRegridStore(field_data, field_dst, routehandle=routehandle, &
                    srcMaskValues=(/srcMaskValue/), &
                    dstMaskValues=(/dstMaskValue/), &
                    regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
                    polemethod=ESMF_POLEMETHOD_ALLAVG, &
                    srcTermProcessing=srcTermProcessing_Value, &
                    ignoreDegenerate=.true., &
                    unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
               if (chkerr(rc,__LINE__,u_FILE_u)) return
            else
               call shr_log_error(' ERROR: regrid_method must be conserve or bilinear - not '//trim(regrid_method), rc=rc)
               return
            end if
         end if

         call ESMF_FieldRegrid(field_data, field_dst, routehandle=routehandle, dynamicMask=dynamicOcnMask, &
              zeroregion=ESMF_REGION_EMPTY, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
      enddo

      ! Close the file
      call pio_closefile(pioid)

      if (pio_iovartype == PIO_REAL) then
         if (allocated(data_real1d)) deallocate(data_real1d)
         if (allocated(data_real2d)) deallocate(data_real2d)
      else if (pio_iovartype == PIO_DOUBLE) then
         if (allocated(data_dbl1d)) deallocate(data_dbl1d)
         if (allocated(data_dbl2d)) deallocate(data_dbl2d)
      endif

      call ESMF_FieldDestroy(field_data)

   end subroutine read_map_input_data

   !===============================================================================
   subroutine set_iodesc(mesh_data, fldname, nlev, pioid, pio_iodesc, rc)

      ! Set pio_iodesc

      ! input/output variables
      type(ESMF_Mesh)   , intent(in)    :: mesh_data
      character(len=*)  , intent(in)    :: fldname
      integer           , intent(in)    :: nlev
      type(file_desc_t) , intent(inout) :: pioid
      type(io_desc_t)   , intent(inout) :: pio_iodesc  ! pio descriptor
      integer           , intent(out)   :: rc

      ! local variables
      integer                 :: gsize2d
      integer                 :: pio_iovartype
      integer                 :: n, m, cnt
      type(var_desc_t)        :: varid
      integer                 :: ndims
      character(CS)           :: dimname
      integer, allocatable    :: dimids(:)
      integer, allocatable    :: dimlens(:)
      type(ESMF_DistGrid)     :: distGrid
      integer                 :: lsize
      integer, pointer        :: compdof(:)
      integer, pointer        :: compdof3d(:)
      integer                 :: rCode ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
      character(*), parameter :: subname = '(set_iodesc) '
      character(*), parameter :: F00  = "('(set_iodesc) ',a,i8,2x,i8,2x,a)"
      character(*), parameter :: F01  = "('(set_iodesc) ',a,i8,2x,i8,2x,a)"
      character(*), parameter :: F02  = "('(set_iodesc) ',a,i8,2x,i8,2x,i8,2x,a)"
      character(*), parameter :: F03  = "('(set_iodesc) ',a,i8,2x,a)"
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! nullify local pointers
      nullify(compdof)
      nullify(compdof3d)

      ! query the variable fldname in the input dataset for its dimensions
      rcode = pio_inq_varid(pioid, trim(fldname), varid)
      rcode = pio_inq_varndims(pioid, varid, ndims)

      ! allocate memory for dimids and dimlens
      allocate(dimids(ndims))
      allocate(dimlens(ndims))
      rcode = pio_inq_vardimid(pioid, varid, dimids(1:ndims))
      do n = 1, ndims
         rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
      end do

      ! determine compdof for input data
      call ESMF_MeshGet(mesh_data, elementdistGrid=distGrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(compdof(lsize))
      call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (nlev > 1) then
         allocate(compdof3d(nlev*lsize))
         ! Assume that first 2 dimensions correspond to the compdof
         gsize2d = dimlens(1)*dimlens(2)
         cnt = 0
         do n = 1,nlev
            do m = 1,size(compdof)
               cnt = cnt + 1
               compdof3d(cnt) = (n-1)*gsize2d + compdof(m)
            enddo
         enddo
      end if

      ! determine type of the variable
      rcode = pio_inq_vartype(pioid, varid, pio_iovartype)

      ! determine io descriptor
      if (ndims == 2) then
         rcode = pio_inq_dimname(pioid, dimids(ndims), dimname)
         if (trim(dimname) == 'time') then
            if (mnproc == 1) then
               write(lp,F03) 'setting iodesc for : '//trim(fldname)// &
                    ' with dimlens(1) = ',dimlens(1),' and the variable has a time dimension '
            end if
            call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1)/), compdof, pio_iodesc)
         else
            if (mnproc == 1) then
               write(lp,F00) 'setting iodesc for : '//trim(fldname)// &
                    ' with dimlens(1), dimlens(2) = ',dimlens(1),dimlens(2),&
                    ' variable has no time dimension '
            end if
            call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
         end if

      else if (ndims == 3) then
         rcode = pio_inq_dimname(pioid, dimids(ndims), dimname)
         if (nlev > 1) then
            if (mnproc == 1) then
               write(lp,F01) 'setting iodesc for : '//trim(fldname)// &
                    ' with dimlens(1), dimlens(2), dimlens(3) = ',dimlens(1),dimlens(2), dimlens(3), &
                    ' variable has no time dimension '//trim(dimname)
            end if
            call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
         else if (trim(dimname) == 'time' .or. trim(dimname) == 'nt') then
            if (mnproc == 1) then
               write(lp,F01) 'setting iodesc for : '//trim(fldname)// &
                    ' with dimlens(1), dimlens(2) = ',dimlens(1),dimlens(2),&
                    ' variable as time dimension '//trim(dimname)
            end if
            call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
         end if

      else if (ndims == 4) then
         rcode = pio_inq_dimname(pioid, dimids(ndims), dimname)
         if (nlev > 1 .and. (trim(dimname) == 'time' .or. trim(dimname) == 'nt')) then
            if (mnproc == 1) then
               write(lp,F02) 'setting iodesc for : '//trim(fldname)// &
                    ' with dimlens(1), dimlens(2),dimlens(3) = ',dimlens(1),dimlens(2),dimlens(3),&
                    ' variable has time dimension '
            end if
            call pio_initdecomp(pio_subsystem, pio_iovartype, (/dimlens(1),dimlens(2),dimlens(3)/), compdof3d, pio_iodesc)
         else
            write(6,*)'ERROR: dimlens= ',dimlens
            call shr_log_error(trim(subname)//' dimlens = 4 assumes a time dimension', rc=rc)
            return
         end if

      else
         write(6,*)'ERROR: dimlens= ',dimlens
         call shr_log_error(trim(subname)//' only ndims of 2 and 3 and 4 are currently supported', rc=rc)
         return
      end if

      ! deallocate memory
      deallocate(compdof)
      if (associated(compdof3d)) deallocate(compdof3d)
      deallocate(dimids)
      deallocate(dimlens)

   end subroutine set_iodesc

   !===============================================================================
   subroutine field_getfldptr(field, fldptr1, fldptr2, rc)

      ! ----------------------------------------------
      ! for a field, determine rank and return fldptr1 or fldptr2
      ! abort is true by default and will abort if fldptr is not yet allocated in field
      ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_Field)  , intent(in)              :: field
      real(r8), pointer , intent(inout), optional :: fldptr1(:)
      real(r8), pointer , intent(inout), optional :: fldptr2(:,:)
      integer           , intent(out)             :: rc

      ! local variables
      type(ESMF_FieldStatus_Flag) :: status
      integer                     :: ungriddedUBound(1)
      character(len=CS)           :: name
      character(len=*), parameter :: subname='(field_getfldptr)'
      ! ----------------------------------------------
      rc = ESMF_SUCCESS

      call ESMF_FieldGet(field, status=status, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (status /= ESMF_FIELDSTATUS_COMPLETE) then
         call ESMF_FieldGet(field, name=name, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         call shr_log_error(trim(subname)//": field "//trim(name)//" has no data not allocated ", rc=rc)
         return
      else
         call ESMF_FieldGet(field, ungriddedUBound=ungriddedUBound, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         if (ungriddedUBound(1) > 0) then
            if (.not.present(fldptr2)) then
               call shr_log_error(trim(subname)//": ERROR missing rank=2 array for "//trim(name), &
                    line=__LINE__, file=u_FILE_u, rc=rc)
               return
            endif
            call ESMF_FieldGet(field, farrayptr=fldptr2, rc=rc)
            if (chkerr(rc,__LINE__,u_FILE_u)) return
         else
            if (.not.present(fldptr1)) then
               call shr_log_error(trim(subname)//": ERROR missing rank=1 array for "//trim(name), &
                    line=__LINE__, file=u_FILE_u, rc=rc)
               return
            endif
            call ESMF_FieldGet(field, farrayptr=fldptr1, rc=rc)
            if (chkerr(rc,__LINE__,u_FILE_u)) return
         end if
      endif  ! status

   end subroutine field_getfldptr

   !===============================================================================
   subroutine fldbun_getFieldN(FB, fieldnum, field, rc)

      ! ----------------------------------------------
      ! Get field with number fieldnum in input field bundle FB
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_FieldBundle), intent(in)    :: FB
      integer               , intent(in)    :: fieldnum
      type(ESMF_Field)      , intent(inout) :: field
      integer               , intent(out)   :: rc

      ! local variables
      character(len=ESMF_MAXSTR) :: name
      character(len=*),parameter :: subname='(dshr_fldbun_getFieldN)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      call fldbun_getNameN(FB, fieldnum, name, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      call ESMF_FieldBundleGet(FB, fieldName=name, field=field, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

   end subroutine fldbun_getFieldN

   !===============================================================================
   subroutine fldbun_getNameN(FB, fieldnum, fieldname, rc)

      ! ----------------------------------------------
      ! Get name of field number fieldnum in input field bundle FB
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_FieldBundle), intent(in)    :: FB
      integer               , intent(in)    :: fieldnum
      character(len=*)      , intent(out)   :: fieldname
      integer               , intent(out)   :: rc

      ! local variables
      integer                         :: fieldCount
      character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
      character(len=*),parameter      :: subname='(getNameN)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      fieldname = ' '
      call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (fieldnum > fieldCount) then
         call shr_log_error(trim(subname)//": ERROR fieldnum > fieldCount ", rc=rc)
         return
      endif

      allocate(lfieldnamelist(fieldCount))
      call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      fieldname = lfieldnamelist(fieldnum)
      deallocate(lfieldnamelist)

   end subroutine fldbun_getNameN

   subroutine dynOcnMaskProc(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

      use ESMF, only : ESMF_RC_ARG_BAD

      ! input/output arguments
      type(ESMF_DynamicMaskElementR8R8R8) , pointer :: dynamicMaskList(:)
      real(ESMF_KIND_R8), intent(in), optional  :: dynamicSrcMaskValue
      real(ESMF_KIND_R8), intent(in), optional  :: dynamicDstMaskValue
      integer           , intent(out)           :: rc

      ! local variables
      integer  :: no, ni
      real(ESMF_KIND_R8)  :: renorm
      !---------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! Below - ONLY if you do NOT have the source masked out then do
      ! the regridding (which is done explicitly here)

      if (associated(dynamicMaskList)) then
         do no = 1, size(dynamicMaskList)
            dynamicMaskList(no)%dstElement = 0.0_r8 ! set to zero
            renorm = 0.d0 ! reset
            do ni = 1, size(dynamicMaskList(no)%factor)
               ! Need to multiply by .90 to handle averaging of input fields before remapping is called
               if ( dynamicMaskList(no)%srcElement(ni) < dynamicSrcMaskValue*.90) then
                  dynamicMaskList(no)%dstElement = dynamicMaskList(no)%dstElement + &
                       (dynamicMaskList(no)%factor(ni) * dynamicMaskList(no)%srcElement(ni))
                  renorm = renorm + dynamicMaskList(no)%factor(ni)
               endif
            enddo
            if (renorm > 0.d0) then
               dynamicMaskList(no)%dstElement = dynamicMaskList(no)%dstElement / renorm
            else if (present(dynamicSrcMaskValue)) then
               dynamicMaskList(no)%dstElement = dynamicSrcMaskValue
            else
               rc = ESMF_RC_ARG_BAD  ! error detected
               return
            endif
         enddo
      endif

   end subroutine DynOcnMaskProc

end module mod_io_input
