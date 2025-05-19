module ocn_map_woa18

   use ESMF              , only : ESMF_Mesh, ESMF_MeshCreate, ESMF_FILEFORMAT_ESMFMESH
   use ESMF              , only : ESMF_Field, ESMF_FieldCreate
   use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleCreate
   use ESMF              , only : ESMF_FieldBundleAdd, ESMF_FieldBundleGet
   use ESMF              , only : ESMF_SUCCESS, ESMF_LogFoundError
   use ESMF              , only : ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
   use nuopc_shr_methods , only : chkerr
   use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use mod_io_input      , only : read_map_input_data, field_getfldptr
   use mod_io_output     , only : io_write
   use mod_xc

   implicit none
   private

   public :: map_woa18

   real(r8), public, allocatable :: woa18_t_depth(:,:,:)
   real(r8), public, allocatable :: woa18_s_depth(:,:,:)
   real(r8), public, allocatable :: depth_bnds(:,:)
   integer , public              :: nlev

   character(len=*), parameter :: u_FILE_u = &
      __FILE__

!===============================================================================
contains
!===============================================================================

   subroutine map_woa18(mesh_blom, rc)

      ! input/out variables
      type(ESMF_Mesh)     , intent(in)  :: mesh_blom
      integer             , intent(out) :: rc

      ! local variables:
      character(len=CL)      :: mesh_input_file
      character(len=CL)      :: filename_t
      character(len=CL)      :: filename_s
      character(len=4)       :: fldlist_input_t(1)
      character(len=4)       :: fldlist_input_s(1)
      type(ESMF_Mesh)        :: mesh_input
      type(ESMF_Field)       :: field_blom
      type(ESMF_FieldBundle) :: fldbun_blom
      integer                :: nf,n,i,j,ko,l
      integer                :: jjcpl
      real(r8), pointer      :: dataptr2d(:,:)
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! ---------------------------
      ! Filenames and field names (TODO: this should be moved to namelist)
      ! ---------------------------

      mesh_input_file = '/cluster/projects/nn9560k/matsbn/WOA_mesh/WOA_1.00_degree_ESMFmesh_20250506_cdf5.nc'
      filename_t = '/cluster/work/users/matsbn/WOA18/woa18_decav_t13_01.nc'
      fldlist_input_t(1) = 't_an'
      filename_s = '/cluster/work/users/matsbn/WOA18/woa18_decav_s13_01.nc'
      fldlist_input_s(1) = 's_an'
      nlev = 102

      ! ---------------------------
      ! Allocate module arrays
      ! ---------------------------

      allocate(woa18_t_depth(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlev))
      allocate(woa18_s_depth(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlev))
      allocate(depth_bnds(2,nlev))

      ! ---------------------------
      ! Create input data mesh
      ! ---------------------------
      mesh_input = ESMF_MeshCreate(filename=trim(mesh_input_file), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! ---------------------------
      ! Create a field bundle on the blom mesh
      ! ---------------------------

      fldbun_blom = ESMF_FieldBundleCreate(name='fldbun_blom', rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      field_blom = ESMF_FieldCreate(mesh_blom, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, name='field_blom',  &
           ungriddedLbound=(/1/), ungriddedUbound=(/nlev/), gridToFieldMap=(/2/), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleAdd(fldbun_blom, (/field_blom/), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! ---------------------------
      ! Read and map temperature - the output will be fldbun_blom which is on the blom mesh
      ! ---------------------------

      ! Read and map the data using bilinear interpolation - and also get depth_bnds
      call read_map_input_data(mesh_input, filename_t, fldlist_input_t, nlev, 'bilinear', &
           fldbun_blom, depth_bnds=depth_bnds, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Plot mapped fldbun temperature
      call io_write(filename="woa18_t_an.nc", fldbun=fldbun_blom, use_float=.false., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Extract the data from the field bundle
      call ESMF_FieldBundleGet(fldbun_blom, fieldName='field_blom', field=field_blom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call field_getfldptr(field_blom, fldptr2=dataptr2d, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Set the j-extent of the local ocean domain to be exchanged. Needed
      ! because of duplication of the last global domain row when using a
      ! tripolar grid.
      if (nreg == 2 .and. nproc == jpr) then
         jjcpl = jj - 1
      else
         jjcpl = jj
      endif

      ! now set woa18_t_depth
      do j = 1, jjcpl
         do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
               n = (j - 1)*ii + i
               do ko = 1,nlev
                  woa18_t_depth(i,j,ko) = dataptr2d(ko,n)
               end do
            end do
         end do
      end do

      ! ---------------------------
      ! Read and map salinity - the output will again be in fldbun_blom on the blom mesh
      ! ---------------------------

      ! Read and map the data using bilinear interpolation
      call read_map_input_data(mesh_input, filename_s, fldlist_input_s, nlev, 'bilinear', &
           fldbun_blom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Extract the data from the field bundle
      call ESMF_FieldBundleGet(fldbun_blom, fieldName='field_blom', field=field_blom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call field_getfldptr(field_blom, fldptr2=dataptr2d, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Plot mapped fldbun salinity
      call io_write(filename="woa18_s_an.nc", fldbun=fldbun_blom, use_float=.false., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! now set woa18_s_depth
      do j = 1, jjcpl
         do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
               n = (j - 1)*ii + i
               do ko = 1,nlev
                  woa18_s_depth(i,j,ko) = dataptr2d(ko,n)
               end do
            end do
         end do
      end do

   end subroutine map_woa18

end module ocn_map_woa18
