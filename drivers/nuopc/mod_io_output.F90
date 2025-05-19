module mod_io_output

   use shr_kind_mod      , only : r8 => shr_kind_r8, r4 => shr_kind_r4, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_log_mod       , only : shr_log_error
   use shr_const_mod     , only : SHR_CONST_SPVAL
   use shr_pio_mod       , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
   use nuopc_shr_methods , only : chkerr
   use ESMF              , only : operator(==)
   use ESMF              , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
   use ESMF              , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundle, ESMF_Mesh, ESMF_DistGrid
   use ESMF              , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
   use ESMF              , only : ESMF_Field, ESMF_FieldGet, ESMF_AttributeGet
   use ESMF              , only : ESMF_MAXSTR
   use pio               , only : file_desc_t, iosystem_desc_t, io_desc_t, var_desc_t
   use pio               , only : pio_openfile, pio_closefile, pio_nowrite, pio_enddef, pio_createfile
   use pio               , only : pio_seterrorhandling, pio_initdecomp, pio_freedecomp
   use pio               , only : pio_inquire, pio_inq_varid, pio_inq_varndims, pio_inq_vardimid
   use pio               , only : pio_inq_dimlen, pio_inq_vartype, pio_inq_dimname, pio_inq_dimid
   use pio               , only : pio_put_att, pio_redef, pio_get_att
   use pio               , only : pio_double, pio_real, pio_int, pio_offset_kind, pio_get_var
   use pio               , only : pio_read_darray, pio_setframe, pio_fill_double, pio_get_att, pio_inq_att
   use pio               , only : PIO_BCAST_ERROR, PIO_RETURN_ERROR, PIO_NOERR, PIO_INTERNAL_ERROR
   use pio               , only : var_desc_t, io_desc_t, pio_offset_kind
   use pio               , only : pio_def_dim, pio_inq_dimid, pio_real, pio_def_var, pio_put_att, pio_double
   use pio               , only : pio_inq_varid, pio_setframe, pio_write_darray, pio_initdecomp, pio_freedecomp
   use pio               , only : pio_syncfile
   use pio               , only : PIO_NOCLOBBER, PIO_IOTYPE_NETCDF, PIO_IOTYPE_PNETCDF, PIO_CLOBBER
   use mod_xc            , only : mnproc, lp

   implicit none
   private

   public :: io_write

   type(iosystem_desc_t), pointer :: pio_subsystem => null()     ! pio info
   integer                        :: io_type                     ! pio info
   integer                        :: io_format                   ! pio info

   character(len=*), parameter :: u_FILE_u = &
      __FILE__

!===============================================================================
contains
!===============================================================================

   subroutine io_write(filename, fldbun, use_float, rc)

      use dimensions, only : itdm, jtdm

      !---------------
      ! Write fldbun to netcdf file
      !---------------

      ! input/output variables
      character(len=*),            intent(in) :: filename
      type(ESMF_FieldBundle)     , intent(in) :: fldbun    ! data to be written
      logical,          optional , intent(in) :: use_float ! write output as float rather than double
      integer                    , intent(out):: rc

      ! local variables
      real(r8)                      :: fillvalue = SHR_CONST_SPVAL ! TODO - make this an input argument
      integer                       :: nx        ! global longitude size
      integer                       :: ny        ! global latitude size
      integer                       :: nz        ! number of vertical levels
      integer                       :: rcode
      integer                       :: nmode
      type(file_desc_t)             :: io_file
      type(ESMF_Field)              :: field
      type(ESMF_Mesh)               :: mesh
      type(ESMF_Distgrid)           :: distgrid
      integer                       :: nf,ng,lsize
      integer                       :: k,m,n,lev
      integer                       :: ndims, nelements
      integer    ,target            :: dimid2(2)
      integer    ,target            :: dimid3(3)
      integer    ,pointer           :: dimid(:)
      type(var_desc_t)              :: varid
      integer, pointer              :: dof(:)
      integer, pointer              :: dof3d(:)
      type(io_desc_t)               :: iodesc
      type(io_desc_t)               :: iodesc3d
      logical                       :: create_iodesc
      logical                       :: create_iodesc3d
      character(CL)                 :: itemc            ! string converted to char
      logical                       :: luse_float
      real(r8), pointer             :: fldptr1(:)
      real(r8), pointer             :: fldptr2(:,:)
      real(r4), allocatable         :: data_real2d(:,:)    ! input data
      real(r8), allocatable         :: data_dbl2d(:,:)     ! input data
      character(CS)                 :: cnumber
      character(CL)                 :: tmpstr
      type(ESMF_Field)              :: lfield
      integer                       :: rank
      integer                       :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fields
      integer                       :: gsize2d
      integer                       :: cnt
      character(CL), allocatable    :: fieldNameList(:)
      character(*), parameter :: F01  = "('(io_write_fldbun) ',a,i8,2x,i8)"
      character(*), parameter :: subName = '(io_write_fldbun) '
      !-------------------------------------------------------------------------------

      rc = ESMF_Success

      luse_float = .false.
      if (present(use_float)) luse_float = use_float

      ! Error check
      if (.not. ESMF_FieldBundleIsCreated(fldbun, rc=rc)) then
         call ESMF_LogWrite(trim(subname)//" fldbun not created", ESMF_LOGMSG_INFO)
         return
      endif

      ! Initialize PIO
      pio_subsystem => shr_pio_getiosys('OCN')
      io_type       =  shr_pio_getiotype('OCN')
      io_format     =  shr_pio_getioformat('OCN')

      nmode = pio_clobber
      ! only applies to classic NETCDF files.
      if (io_type == PIO_IOTYPE_NETCDF .or. io_type == PIO_IOTYPE_PNETCDF) then
         nmode = ior(nmode,io_format)
      endif
      rcode = pio_createfile(pio_subsystem, io_file, io_type, trim(filename), nmode)
      if (mnproc == 1) then
         write(lp,'(a)') trim(subname) //' creating file '// trim(filename)
         write(lp,'(a)') 'opening   : '//trim(filename)
      endif

      ! Get number of fields
      call ESMF_FieldBundleGet(fldbun, fieldCount=nf, rc=rc)
      write(tmpstr,*) subname//' field count = ', nf
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
      if (nf < 1) then
         rc = ESMF_Success
         return
      endif
      allocate(fieldNameList(nf))
      call ESMF_FieldBundleGet(fldbun, fieldNameList=fieldNameList, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Get field bundle mesh from first field
      call fldbun_getFieldN(fldbun, 1, field, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(field, mesh=mesh, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Get mesh distgrid and number of elements
      call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_MeshGet(mesh, spatialDim=ndims, numOwnedElements=nelements, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      write(tmpstr,*) subname, 'ndims, nelements = ', ndims, nelements
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

      ! Define coordinate

      ! *** ASSUME that all variables in the field bundle have the same number of dimensions ***
      ! Determine rank of first field in input fldbun
      call ESMF_FieldBundleGet(fldbun, fieldnamelist(1),  field=lfield, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(lfield, rank=rank, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Determine global nx, ny
      nx = itdm
      ny = jtdm

      if (rank == 2) then
         rcode = pio_def_dim(io_file, 'lon', nx, dimid3(1))
         rcode = pio_def_dim(io_file, 'lat', ny, dimid3(2))
         call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         write(cnumber,'(i0)') ungriddedUbound(1)
         if (mnproc == 1) then
            write(lp, '(a)') trim(subname)//':'//trim(fieldnamelist(1))// &
                 ' has an griddedUBound of  '//trim(cnumber)
         end if
         rcode = pio_def_dim(io_file, 'depth', ungriddedUbound(1), dimid3(3))
         dimid => dimid3
         nz = ungriddedUBound(1)
      else
         rcode = pio_def_dim(io_file, 'lon', nx, dimid2(1))
         rcode = pio_def_dim(io_file, 'lat', ny, dimid2(2))
         dimid => dimid2
         nz = 1
      end if
      ndims = size(dimid)

      ! -------------------------
      ! Write header
      ! -------------------------

      ! Define variable and attribute on output file
      do k = 1,size(fieldNameList)
         itemc = trim(fieldNameList(k))
         if (luse_float) then
            rcode = pio_def_var(io_file, trim(itemc), PIO_REAL, dimid, varid)
            rcode = pio_put_att(io_file, varid, "_FillValue", real(fillvalue, r4))
         else
            rcode = pio_def_var(io_file, trim(itemc), PIO_DOUBLE, dimid, varid)
            rcode = pio_put_att(io_file, varid, "_FillValue", fillvalue)
         end if
      end do

      ! -------------------------------
      ! End define mode
      ! -------------------------------

      rcode = pio_enddef(io_file)

      ! -------------------------------
      ! Create iodesc and iodesc3d
      ! -------------------------------

      ! Determine dof
      nullify(dof)
      nullify(dof3d)

      ! Use distgrid extracted from field 1 above
      call ESMF_DistGridGet(distgrid, localDE=0, elementCount=lsize, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(dof(lsize))
      call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
      if (nz > 1) then
         allocate(dof3d(lsize*nz))
         ! Assume that first 2 dimensions correspond to the dof
         gsize2d = nx*ny
         cnt = 0
         do n = 1,nz
            do m = 1,size(dof)
               cnt = cnt + 1
               dof3d(cnt) = (n-1)*gsize2d + dof(m)
            enddo
         enddo
      end if

      ! Create iodesc and iodesc3d
      create_iodesc = .true.
      create_iodesc3d = .true.
      do k = 1,size(fieldNameList)
         itemc = trim(fieldNameList(k))

         ! Determine rank of field with name itemc
         call ESMF_FieldBundleGet(fldbun, itemc,  field=lfield, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         call ESMF_FieldGet(lfield, rank=rank, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         if (create_iodesc) then
            ! Always create 2d iodesc
            if (mnproc == 1) then
               write(lp,F01) 'setting iodesc for : '//trim(itemc)//' with dims(1:2) = ',nx,ny
            end if
            if (luse_float) then
               call pio_initdecomp(pio_subsystem, pio_real, (/nx,ny/), dof, iodesc)
            else
               call pio_initdecomp(pio_subsystem, pio_double, (/nx,ny/), dof, iodesc)
            end if
            create_iodesc = .false.
         end if
         if (rank < 2) then
            create_iodesc3d = .false.
         else if (rank == 2) then
            if (create_iodesc3d) then
               if (luse_float) then
                  call pio_initdecomp(pio_subsystem, pio_real, (/nx,ny,nz/), dof3d, iodesc3d)
               else
                  call pio_initdecomp(pio_subsystem, pio_double, (/nx,ny,nz/), dof3d, iodesc3d)
               end if
               create_iodesc3d = .false.
            end if
         end if
      end do

      deallocate(dof)
      deallocate(dof3d)

      ! -------------------------------
      ! Write data
      ! -------------------------------

      do k = 1,size(fieldNameList)
         ! Determine field name
         itemc = trim(fieldNameList(k))

         ! Get rank and fldptr1 or fldptr2 depending on rank
         call fldbun_getFldPtr(fldbun, itemc, fldptr1=fldptr1, fldptr2=fldptr2, rank=rank, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         if (rank == 2) then

            rcode = pio_inq_varid(io_file, trim(itemc), varid)
            ! Need to swap depth to be second dimension
            if (luse_float) then
               allocate(data_real2d(lsize, nz))
               do lev = 1,nz
                  do n = 1,size(fldptr2, dim=2)
                     data_real2d(n,lev) = fldptr2(lev,n)
                  end do
               end do
               call pio_write_darray(io_file, varid, iodesc3d, real(data_real2d,r4), rcode, fillval=real(fillvalue,r4))
               deallocate(data_real2d)
            else
               allocate(data_dbl2d(lsize, nz))
               do lev = 1,nz
                  do n = 1,size(fldptr2, dim=2)
                     data_dbl2d(n,lev) = fldptr2(lev,n)
                  end do
               end do
               call pio_write_darray(io_file, varid, iodesc3d, data_dbl2d, rcode, fillval=fillvalue)
               deallocate(data_dbl2d)
            end if

         else if (rank == 1) then

            rcode = pio_inq_varid(io_file, trim(itemc), varid)
            if (luse_float) then
               call pio_write_darray(io_file, varid, iodesc, real(fldptr1,r4), rcode, fillval=real(fillvalue,r4))
            else
               call pio_write_darray(io_file, varid, iodesc, fldptr1, rcode, fillval=fillvalue)
            end if
         end if  ! end if rank is 2 or 1

      end do  ! end loop over fields in fldbun

      !call pio_syncfile(io_file)
      call pio_freedecomp(io_file, iodesc)
      call pio_freedecomp(io_file, iodesc3d)

      call pio_closefile(io_file)

   end subroutine io_write

   !===============================================================================
   subroutine fldbun_getFieldN(fldbun, fieldnum, field, rc)

      ! ----------------------------------------------
      ! Get field with number fieldnum in input field bundle fldbun
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_FieldBundle), intent(in)    :: fldbun
      integer               , intent(in)    :: fieldnum
      type(ESMF_Field)      , intent(inout) :: field
      integer               , intent(out)   :: rc

      ! local variables
      character(len=ESMF_MAXSTR) :: name
      character(len=*),parameter :: subname='(fldbun_getFieldN)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      call fldbun_getNameN(fldbun, fieldnum, name, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleGet(fldbun, fieldName=name, field=field, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

   end subroutine fldbun_getFieldN

   !===============================================================================
   subroutine fldbun_getNameN(fldbun, fieldnum, fieldname, rc)

      ! ----------------------------------------------
      ! Get name of field number fieldnum in input field bundle fldbun
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_FieldBundle), intent(in)    :: fldbun
      integer               , intent(in)    :: fieldnum
      character(len=*)      , intent(out)   :: fieldname
      integer               , intent(out)   :: rc

      ! local variables
      integer                         :: fieldCount
      character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
      character(len=*),parameter      :: subname='(fldbun_getNameN)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      fieldname = ' '
      call ESMF_FieldBundleGet(fldbun, fieldCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (fieldnum > fieldCount) then
         call shr_log_error(trim(subname)//": ERROR fieldnum > fieldCount ", rc=rc)
         return
      endif
      allocate(lfieldnamelist(fieldCount))
      call ESMF_FieldBundleGet(fldbun, fieldNameList=lfieldnamelist, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      fieldname = lfieldnamelist(fieldnum)
      deallocate(lfieldnamelist)

   end subroutine fldbun_getNameN

   !===============================================================================
   subroutine fldbun_GetFldPtr(fldbun, fldname, fldptr1, fldptr2, rank, field, rc)

      ! ----------------------------------------------
      ! Get pointer to a field bundle field
      ! ----------------------------------------------

      type(ESMF_FieldBundle) , intent(in)              :: fldbun
      character(len=*)       , intent(in)              :: fldname
      real(R8), pointer      , intent(inout), optional :: fldptr1(:)
      real(R8), pointer      , intent(inout), optional :: fldptr2(:,:)
      integer                , intent(out),   optional :: rank
      integer                , intent(out),   optional :: rc
      type(ESMF_Field)       , intent(out),   optional :: field

      ! local variables
      type(ESMF_Field) :: lfield
      integer          :: lrank
      character(len=*), parameter :: subname='(fldbun_GetFldPtr)'
      ! ----------------------------------------------

      if (.not.present(rc)) then
         call shr_log_error(trim(subname)//": ERROR rc not present "//trim(fldname), &
              line=__LINE__, file=u_FILE_u, rc=rc)
         return
      endif

      rc = ESMF_SUCCESS

      if (.not. fldbun_FldChk(fldbun, trim(fldname), rc=rc)) then
         call shr_log_error(trim(subname)//": ERROR field "//trim(fldname)//" not in fldbun ", &
              line=__LINE__, file=u_FILE_u, rc=rc)
         return
      endif

      call ESMF_FieldBundleGet(fldbun, fieldName=trim(fldname), field=lfield, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call field_GetFldPtr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (present(rank)) then
         rank = lrank
      endif
      if (present(field)) then
         field = lfield
      endif

   end subroutine fldbun_GetFldPtr

   !===============================================================================
   subroutine field_GetFldPtr(field, fldptr1, fldptr2, rank, abort, rc)

      ! ----------------------------------------------
      ! for a field, determine rank and return fldptr1 or fldptr2
      ! abort is true by default and will abort if fldptr is not yet allocated in field
      ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_Field)  , intent(in)              :: field
      real(R8), pointer , intent(inout), optional :: fldptr1(:)
      real(R8), pointer , intent(inout), optional :: fldptr2(:,:)
      integer           , intent(out)  , optional :: rank
      logical           , intent(in)   , optional :: abort
      integer           , intent(out)  , optional :: rc

      ! local variables
      type(ESMF_Mesh)          :: lmesh
      integer                  :: lrank, nnodes, nelements
      logical                  :: labort
      character(len=*), parameter :: subname='(field_GetFldPtr)'
      ! ----------------------------------------------

      if (.not.present(rc)) then
         call shr_log_error(trim(subname)//": ERROR rc not present ", &
              line=__LINE__, file=u_FILE_u, rc=rc)
         return
      endif

      rc = ESMF_SUCCESS

      labort = .true.
      if (present(abort)) then
         labort = abort
      endif
      lrank = -99

      call ESMF_FieldGet(field, rank=lrank, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (lrank == 1) then
         if (.not.present(fldptr1)) then
            call shr_log_error(trim(subname)//": ERROR missing rank=1 array ", &
                 line=__LINE__, file=u_FILE_u, rc=rc)
            return
         endif
         call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
      elseif (lrank == 2) then
         if (.not.present(fldptr2)) then
            call shr_log_error(trim(subname)//": ERROR missing rank=2 array ", &
                 line=__LINE__, file=u_FILE_u, rc=rc)
            return
         endif
         call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
      else
         call shr_log_error(trim(subname)//": ERROR in rank ", &
              line=__LINE__, file=u_FILE_u, rc=rc)
         return
      endif

      if (present(rank)) then
         rank = lrank
      endif
   end subroutine field_GetFldPtr

   !===============================================================================
   logical function fldbun_FldChk(fldbun, fldname, rc)

      ! ----------------------------------------------
      ! Determine if field with fldname is in input field bundle
      ! ----------------------------------------------

      ! input/output variables
      type(ESMF_FieldBundle), intent(in)  :: fldbun
      character(len=*)      , intent(in)  :: fldname
      integer               , intent(out) :: rc

      ! local variables
      logical :: ispresent
      character(len=*), parameter :: subname='(fldbun_FldChk)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! If field bundle is not created then set return to .false.
      if (.not. ESMF_FieldBundleIsCreated(fldbun)) then
         fldbun_FldChk = .false.
         return
      end if

      ! If field bundle is created determine if fldname is present in field bundle
      fldbun_FldChk = .false.

      call ESMF_FieldBundleGet(fldbun, fieldName=trim(fldname), isPresent=isPresent, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) then
         call shr_log_error(string=trim(subname)//" Error checking field: "//trim(fldname), &
              line=__LINE__,file=u_FILE_u, rc=rc)
         return
      endif
      if (isPresent) then
         fldbun_FldChk = .true.
      endif

   end function fldbun_FldChk

end module mod_io_output
