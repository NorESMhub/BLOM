module mod_blom_pio

   use ESMF
   use pio
   use mod_types         , only : r8, r4
   use mod_xc            , only : lp, mnproc
   use shr_pio_mod       , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
   use shr_sys_mod       , only : shr_sys_abort
   use nuopc_shr_methods , only : chkerr

   implicit none
   private

   public  :: blom_pio_init
   public  :: blom_pio_getdata
   private :: blom_pio_iodesc

  integer                        , public :: pio_iotype
  integer                        , public :: pio_ioformat
  type(iosystem_desc_t), pointer , public :: pio_iosystem

   character(len=*), parameter :: u_FILE_u = &
      __FILE__

#include <mpif.h>

!================================================================================
contains
!================================================================================

   subroutine blom_pio_init()
      pio_iosystem => shr_pio_getiosys('OCN')
      pio_iotype   =  shr_pio_getiotype('OCN')
      pio_ioformat =  shr_pio_getioformat('OCN')
   end subroutine blom_pio_init

   !===============================================================
   subroutine blom_pio_getdata(pioid, varname, mesh_i, data_i, rc)

      ! input/output variables
      type(file_desc_t), intent(inout) :: pioid
      character(len=*) , intent(in)    :: varname   ! field name in rawdata file
      type(ESMF_Mesh)  , intent(in)    :: mesh_i
      real(r8)         , intent(inout) :: data_i(:) ! input raw data
      integer          , intent(out)   :: rc

      ! local variables
      type(io_desc_t)       :: pio_iodesc
      type(var_desc_t)      :: pio_varid
      integer               :: pio_vartype
      real(r4), allocatable :: data_real(:)
      integer               :: lsize
      integer               :: rcode
      character(len=*), parameter :: subname = 'blom_pio_getdata'
      !-------------------------------------------------

      rc = ESMF_SUCCESS

      ! Create io descriptor for input raw data
      ! This will query the raw data file for the dimensions of the variable varname and
      ! create iodesc for either single or multi level input data
      call blom_pio_iodesc(mesh_i, trim(varname), pioid, pio_varid, pio_vartype, pio_iodesc, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Read the input data
      if (pio_vartype == PIO_REAL) then
         lsize = size(data_i)
         allocate(data_real(lsize))
         call pio_read_darray(pioid, pio_varid, pio_iodesc, data_real, rcode)
         data_i(:) = real(data_real(:), kind=r8)
         deallocate(data_real)
      else if (pio_vartype == PIO_DOUBLE) then
         call pio_read_darray(pioid, pio_varid, pio_iodesc, data_i, rcode)
      else
         call shr_sys_abort(subName//" ERROR: supported variable type not found for "//trim(varname))
      end if

      call pio_freedecomp(pioid, pio_iodesc)

   end subroutine blom_pio_getdata

   !===============================================================
   subroutine blom_pio_iodesc( mesh, varname, pioid, pio_varid,  pio_vartype, pio_iodesc, rc)

      ! Determine pio io descriptor for variable on rawdata file

      ! input/output variables
      type(ESMF_Mesh)   , intent(in)    :: mesh
      character(len=*)  , intent(in)    :: varname
      type(file_desc_t) , intent(inout) :: pioid
      type(var_desc_t)  , intent(out)   :: pio_varid
      integer           , intent(out)   :: pio_vartype
      type(io_desc_t)   , intent(inout) :: pio_iodesc
      integer           , intent(out)   :: rc

      ! local variables
      type(ESMF_DistGrid)     :: distGrid
      integer                 :: n, ndims
      integer, allocatable    :: compdof(:)
      integer, allocatable    :: dimids(:)
      integer, allocatable    :: dimlens(:)
      integer                 :: lsize
      integer                 :: rCode ! pio return code (only used when pio error handling is PIO_BCAST_ERROR)
      character(*), parameter :: subname = '(pio_iodesc) '
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(compdof(lsize))
      call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=compdof, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! get pio variable id, type and number of dimensions
      rcode = pio_inq_varid(pioid, trim(varname), pio_varid)
      rcode = pio_inq_vartype(pioid, pio_varid, pio_vartype)
      rcode = pio_inq_varndims(pioid, pio_varid, ndims)

      ! get variable dimension sizes
      allocate(dimids(ndims))
      allocate(dimlens(ndims))
      rcode = pio_inq_vardimid(pioid, pio_varid, dimids(1:ndims))
      do n = 1, ndims
         rcode = pio_inq_dimlen(pioid, dimids(n), dimlens(n))
      end do

      ! determine io descriptor for this variable
      call pio_initdecomp(pio_iosystem, pio_vartype, (/dimlens(1),dimlens(2)/), compdof, pio_iodesc)
      if (mnproc == 1) then
         write(lp,'(a,i8,i8)') ' set iodesc for rawdata: '//trim(varname)//' with dim(1),dim(2) = ',&
              dimlens(1),dimlens(2)
      end if

      ! deallocate memory
      deallocate(compdof)

   end subroutine blom_pio_iodesc

end module mod_blom_pio
