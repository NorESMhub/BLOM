module ocn_map_woa

   use ESMF              , only : ESMF_Mesh, ESMF_MeshCreate, ESMF_FILEFORMAT_ESMFMESH
   use ESMF              , only : ESMF_Field, ESMF_FieldCreate
   use ESMF              , only : ESMF_FieldBundle, ESMF_FieldBundleCreate
   use ESMF              , only : ESMF_FieldBundleAdd, ESMF_FieldBundleGet
   use ESMF              , only : ESMF_SUCCESS, ESMF_LogFoundError
   use ESMF              , only : ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
   use nuopc_shr_methods , only : chkerr
   use shr_kind_mod      , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_const_mod     , only : SHR_CONST_SPVAL
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use mod_io_input      , only : read_map_input_data, field_getfldptr
   use mod_io_output     , only : io_write
   use mod_cesm          , only : woa_nuopc_icfile_mesh, woa_nuopc_icfile_data
   use shr_pio_mod       , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
   use pio               , only : file_desc_t, iosystem_desc_t
   use pio               , only : pio_openfile, pio_closefile, pio_nowrite
   use pio               , only : pio_inq_dimid, pio_inq_dimlen
   use mod_inicon        , only : t_woa, s_woa, depth_bnds_woa, depth_woa
   use mod_inicon        , only : t_woa_fval, s_woa_fval, kdm_woa
   use mod_cesm          , only : woa_nuopc_icfile_mesh, woa_nuopc_icfile_data
   use mod_xc

   implicit none
   private

   public :: map_woa

   character(len=*), parameter :: u_FILE_u = &
      __FILE__

!===============================================================================
contains
!===============================================================================

   subroutine map_woa(mesh_blom, rc)

      ! input/out variables
      type(ESMF_Mesh)     , intent(in)  :: mesh_blom
      integer             , intent(out) :: rc

      ! local variables:
      character(len=CL)      :: mesh_input_file
      character(len=CL)      :: filename_ts
      character(len=4)       :: fldlist_input(2)
      type(ESMF_Mesh)        :: mesh_input
      type(ESMF_Field)       :: field_blom
      type(ESMF_FieldBundle) :: fldbun_blom
      integer                :: nf,n,i,j,ko,l
      integer                :: jjcpl
      real(r8), pointer      :: dataptr2d(:,:)
      type(file_desc_t)      :: pioid
      integer                :: dimid
      integer                :: rcode
      integer                :: io_type                         ! pio info
      integer                :: io_format                       ! pio info
      type(iosystem_desc_t), pointer :: pio_subsystem => null() ! pio info
      character(*), parameter :: subname = '(map_woa) '
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! ---------------------------
      ! Filenames and field names
      ! ---------------------------

      mesh_input_file = trim(woa_nuopc_icfile_mesh)
      filename_ts = trim(woa_nuopc_icfile_data)
      fldlist_input(1) = 't_an'
      fldlist_input(2) = 's_an'

      ! ---------------------------
      ! Determine vertical dimension
      ! ---------------------------

      pio_subsystem => shr_pio_getiosys('OCN')
      io_type       =  shr_pio_getiotype('OCN')
      io_format     =  shr_pio_getioformat('OCN')
      if (mnproc == 1) then
         write(lp,'(a)') trim(subname) // ' determining vertical dimension of WOA climatology from '//trim(filename_ts)
      end if
      rcode = pio_openfile(pio_subsystem, pioid, io_type, trim(filename_ts), pio_nowrite)
      rcode = pio_inq_dimid(pioid, 'depth', dimid)
      rcode = pio_inq_dimlen(pioid, dimid, kdm_woa)
      call pio_closefile(pioid)

      ! ---------------------------
      ! Set fill values
      ! ---------------------------

      t_woa_fval = SHR_CONST_SPVAL
      s_woa_fval = SHR_CONST_SPVAL

      ! ---------------------------
      ! Allocate module arrays in mod_inicon.F90
      ! ---------------------------

      allocate(t_woa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm_woa), &
               s_woa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm_woa), &
               depth_woa(kdm_woa), &
               depth_bnds_woa(2,kdm_woa), &
               stat = rcode)
      if (rcode /= 0) then
         write(lp,*) 'Failed to allocate WOA arrays!'
         call xchalt('(ocn_map_woa)')
         stop '(ocn_map_woa)'
      endif

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
           ungriddedLbound=(/1/), ungriddedUbound=(/kdm_woa/), gridToFieldMap=(/2/), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleAdd(fldbun_blom, (/field_blom/), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! ---------------------------
      ! Read and map temperature - the output will be fldbun_blom which is on the blom mesh
      ! ---------------------------

      ! Read and map the data using bilinear interpolation - and also get depth_bnds and depths
      if (mnproc == 1) then
         write(lp,*)
         write(lp,'(a)') trim(subname) // ' calling read_map_inputdata for '//trim(fldlist_input(1))
      end if
      call read_map_input_data(mesh_input, filename_ts, (/fldlist_input(1)/), kdm_woa, 'conserve', &
           fldbun_blom, depth=depth_woa, depth_bnds=depth_bnds_woa, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Plot mapped fldbun temperature
      if (mnproc == 1) then
         write(lp,*)
         write(lp,'(a)') trim(subname) //' plotting mapped data for '//trim(fldlist_input(1))
      end if
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

      ! now set t_woa
      do j = 1, jjcpl
         do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
               n = (j - 1)*ii + i
               do ko = 1,kdm_woa
                  t_woa(i,j,ko) = dataptr2d(ko,n)
               end do
            end do
         end do
      end do

      ! ---------------------------
      ! Read and map salinity - the output will again be in fldbun_blom on the blom mesh
      ! ---------------------------

      ! Read and map the data using bilinear interpolation
      if (mnproc == 1) then
         write(lp,*)
         write(lp,'(a)') trim(subname) // ' calling read_map_inputdata for '//trim(fldlist_input(2))
      end if
      call read_map_input_data(mesh_input, filename_ts, (/fldlist_input(2)/), kdm_woa, 'conserve', &
           fldbun_blom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Plot mapped fldbun salinity
      if (mnproc == 1) then
         write(lp,*)
         write(lp,'(a)') trim(subname) // ' plotting mapped data for '//trim(fldlist_input(2))
      end if
      call io_write(filename="woa18_s_an.nc", fldbun=fldbun_blom, use_float=.false., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Extract the data from the field bundle
      call ESMF_FieldBundleGet(fldbun_blom, fieldName='field_blom', field=field_blom, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call field_getfldptr(field_blom, fldptr2=dataptr2d, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! now set s_woa
      do j = 1, jjcpl
         do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
               n = (j - 1)*ii + i
               do ko = 1,kdm_woa
                  s_woa(i,j,ko) = dataptr2d(ko,n)
               end do
            end do
         end do
      end do

   end subroutine map_woa

end module ocn_map_woa
