module ocn_map_river_nutrients

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
   use mo_read_rivin     , only : rivflx(:,:,:)

   use mod_xc

   implicit none
   private

   public :: map_river_nutrients

   character(len=*), parameter :: u_FILE_u = &
        __FILE__

!===============================================================================
contains
!===============================================================================

   subroutine map_river_nutrients(mesh_blom, rc)

      ! input/out variables
      type(ESMF_Mesh)     , intent(in)  :: mesh_blom
      integer             , intent(out) :: rc

      ! local variables:
      character(len=CL)      :: mesh_input_file
      character(len=CL)      :: filename
      character(len=3)       :: fldlist_input(7)
      character(len=CL)      :: mapfile
      type(ESMF_Mesh)        :: mesh_input
      type(ESMF_Field)       :: field_blom
      type(ESMF_FieldBundle) :: fldbun_blom
      integer                :: nf,n,i,j,ko,l
      integer                :: jjcpl
      real(r8), pointer      :: dataptr1d(:)
      integer                :: nlev
      integer                :: nfld
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      filename = '/cluster/shared/noresm/inputdata/ocn/blom/bndcon/river_nutrients_GNEWS2000c00_r05_20250220.nc'
      mesh_input_file = '/cluster/shared/noresm/inputdata/share/meshes/r05_nomask_c110308_ESMFmesh.nc'
      fldlist_input = (/'DIN','DIP','DSi','DIC','Fe ','DOC','DET'/)
      mapfile = '/cluster/shared/noresm/inputdata/cpl/cpl6/map_r05_to_tnx1v4_e1000r300_170609.nc'
      nlev = 1

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
      do nfld = 1,size(fldlist_input)
         field_blom = ESMF_FieldCreate(mesh_blom, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, name=fldlist_input(nfld), rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         call ESMF_FieldBundleAdd(fldbun_blom, (/field_blom/), rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
      end do

      ! ---------------------------
      ! Read and map riverin nutrients
      ! ---------------------------

      call read_map_input_data(mesh_input, filename, fldlist_input, nlev, 'mapfile', &
           fldbun_blom, mapfile=trim(mapfile), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Plot mapped fldbun fields
      call io_write(filename="river_nutrients.nc", fldbun=fldbun_blom, use_float=.false., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Allocate field to hold river fluxes
      if (mnproc == 1) then
         write(io_stdo_bgc,*)'Memory allocation for variable rivflx ...'
         write(io_stdo_bgc,*)'First dimension    : ',idm
         write(io_stdo_bgc,*)'Second dimension   : ',kdm
         write(io_stdo_bgc,*)'Third  dimension   : ',size(inputdata_list)
      endif

      allocate (rivflx(idm,kdm,nriv),stat=errstat)
      if(errstat./=0) stop 'not enough memory rivflx'
      rivflx(:,:,:) = 0.0

      ! Loop over the input fields
      do nfld = 1,size(fldlist_input)
         ! Extract the field from fldbun_blom
         call ESMF_FieldBundleGet(fldbun_blom, fieldName=trim(fldlist_input(nfld)), field=field_blom, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         call field_getfldptr(field_blom, fldptr1=dataptr1d, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         fldlist_input = (/'DIN','DIP','DSi','DIC','Fe ','DOC','DET'/)

         if (trim(fldlist_input(nfld)) == 'DIN') then
            index = 1
         else if (trim(fldlist_input(nfld)) == 'DIP') then
            index = 2
         else if (trim(fldlist_input(nfld)) == 'DSi') then
            index = 3
         else if (trim(fldlist_input(nfld)) == 'DIC') then
            index = 4
         else if (trim(fldlist_input(nfld)) == 'Fe') then
            index = 5
         else if (trim(fldlist_input(nfld)) == 'DOC') then
            index = 6
         else if (trim(fldlist_input(nfld)) == 'DOT') then
            index = 7
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
               n = (j - 1)*ii + i
               if (omask(i,j) > 0.5) then
                  rivflx(i,j,index)  = dataptr1d(n)
               end if
            end do
         end do
      end do

   end subroutine map_river_nutrients

end module ocn_map_river_nutrients
