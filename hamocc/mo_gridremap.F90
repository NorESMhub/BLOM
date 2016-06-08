module mo_gridremap

   use types

   implicit none

   integer, dimension(:), allocatable :: &
      o_weight_index, a_weight_index
   ! Weights for conservative remapping from ocean to atmospheric grid
   real (r8), dimension(:), allocatable :: o_weight_cnsrv
   ! Weights for scalar remapping from ocean to atmospheric grid
   real (r8), dimension(:), allocatable :: o_weight_scalr
   ! Weights for conservative remapping from atmospheric to ocean grid
   real (r8), dimension(:), allocatable :: a_weight_cnsrv
   ! Weights for scalar remapping from atmospheric to ocean grid
   real (r8), dimension(:), allocatable :: a_weight_scalr

   ! Number of polygons for ocean, atmosphere, and intersecting ocean and
   ! atmosphere polygons
   integer :: nop, nap, noap

   real (r8), parameter :: &
      area_eps = 1.e-8_r8 ! Small area in radians squared
   
contains

   subroutine remap_init(file_name, nop_req, nap_req)
      ! Initialize remapping of fields between an ocean and atmosphere grid

      use netcdf
      use ncutils

      implicit none

      ! Input variables:
      !    file_name: name of netcdf file containing remapping areas and
      !               addresses
      !    nop_req:   the required number of ocean polygons
      !    nap_req:   the required number of atmosphere polygons
      character (len = *), intent(in) :: file_name
      integer, intent(in) :: nop_req, nap_req

      ! Local variables
      real (r8), dimension(:), allocatable :: &
         o_poly_area, a_poly_area, oa_poly_area, o_oa_poly_area, a_oa_poly_area
      integer :: ncid, ioap, iop, iap
      logical :: lerr

      call ncopen(file_name, ncid)

      ! Check number of polygons in the two grids
      call ncgetdimlen(ncid, 'num_a_poly', nop)
      call ncgetdimlen(ncid, 'num_b_poly', nap)
      lerr = .false.
      if (nop_req /= nop) then
         write (*,*) &
            ' Inconsistent number of ocean polygons: required: ', nop_req, &
            ', in file:',nop
         lerr = .true.
      endif
      if (nap_req /= nap) then
         write (*,*) &
           ' Inconsistent number of atmosphere polygons: required: ', nap_req, &
           ', in file:',nap
         lerr = .true.
      endif
      if (lerr) then
         stop
      endif

      ! Get number of intersection grid polygons
      call ncgetdimlen(ncid, 'num_ab_poly', noap)

      ! Allocate some arrays
      allocate(o_poly_area(nop))
      allocate(a_poly_area(nap))
      allocate(oa_poly_area(noap))
      allocate(o_weight_index(noap))
      allocate(a_weight_index(noap))
      allocate(o_oa_poly_area(nop))
      allocate(a_oa_poly_area(nap))
      allocate(o_weight_cnsrv(noap))
      allocate(o_weight_scalr(noap))
      allocate(a_weight_cnsrv(noap))
      allocate(a_weight_scalr(noap))

      ! Read polygon areas and indexes
      call ncgetvar(ncid, 'a_poly_area', o_poly_area)
      call ncgetvar(ncid, 'b_poly_area', a_poly_area)
      call ncgetvar(ncid, 'ab_poly_area', oa_poly_area)
      call ncgetvar(ncid, 'a_poly_index', o_weight_index)
      call ncgetvar(ncid, 'b_poly_index', a_weight_index)

      call ncclose(ncid)

      ! Compute area of grid cells overlapped by the other grid
      o_oa_poly_area(:) = 0._r8
      a_oa_poly_area(:) = 0._r8
      do ioap = 1, noap
         iop = o_weight_index(ioap)
         iap = a_weight_index(ioap)
         o_oa_poly_area(iop) = o_oa_poly_area(iop) + oa_poly_area(ioap)
         a_oa_poly_area(iap) = a_oa_poly_area(iap) + oa_poly_area(ioap)
      enddo

!$OMP PARALLEL DO
      do iop = 1, nop
         if (o_poly_area(iop) - o_oa_poly_area(iop) < area_eps) then
            o_poly_area(iop) = o_oa_poly_area(iop)
         endif
      enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      do iap = 1, nap
         if (a_poly_area(iap) - a_oa_poly_area(iap) < area_eps) then
            a_poly_area(iap) = a_oa_poly_area(iap)
         endif
      enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      do ioap = 1, noap
         iop = o_weight_index(ioap)
         iap = a_weight_index(ioap)
         o_weight_cnsrv(ioap) = oa_poly_area(ioap)/o_poly_area(iop)
         o_weight_scalr(ioap) = oa_poly_area(ioap) &
                                /max(area_eps,a_oa_poly_area(iap))
         a_weight_cnsrv(ioap) = oa_poly_area(ioap)/a_poly_area(iap)
         a_weight_scalr(ioap) = oa_poly_area(ioap) &
                                /max(area_eps,o_oa_poly_area(iop))
      enddo
!$OMP END PARALLEL DO

      deallocate(o_poly_area)
      deallocate(a_poly_area)
      deallocate(oa_poly_area)
      deallocate(o_oa_poly_area)
      deallocate(a_oa_poly_area)
      
   end subroutine remap_init

   subroutine remap_o2a_cnsrv(o_area, a_area, o_field, a_field)
      ! Remap a field on an ocean grid to an atmosphere grid, conserving the
      ! area weighted sum. Suitable for flux fields.

      use types

      real (r8), dimension(nop), intent(in) :: o_area, o_field
      real (r8), dimension(nap), intent(in) :: a_area
      real (r8), dimension(nap), intent(out) :: a_field

      integer :: ioap, iop, iap

      a_field(:) = 0._r8
      do ioap = 1, noap
         iop = o_weight_index(ioap)
         iap = a_weight_index(ioap)
         a_field(iap) = a_field(iap) &
                      + o_field(iop)*o_area(iop)*o_weight_cnsrv(ioap)
      enddo
!$OMP PARALLEL DO
      do iap = 1, nap
         a_field(iap) = a_field(iap)/max(area_eps, a_area(iap))
      enddo
!$OMP END PARALLEL DO

   end subroutine remap_o2a_cnsrv

   subroutine remap_o2a_scalr(o_field, a_field)
      ! Remap a field on an ocean grid to an atmosphere grid. Suitable for
      ! scalar fields

      use types

      real (r8), dimension(nop), intent(in) :: o_field
      real (r8), dimension(nap), intent(out) :: a_field

      integer :: ioap, iop, iap

      a_field(:) = 0._r8
      do ioap = 1, noap
         iop = o_weight_index(ioap)
         iap = a_weight_index(ioap)
         a_field(iap) = a_field(iap) + o_field(iop)*o_weight_scalr(ioap)
      enddo

   end subroutine remap_o2a_scalr

   subroutine remap_a2o_cnsrv(a_area, o_area, a_field, o_field)
      ! Remap a field on an atmosphere grid to an ocean grid, conserving the
      ! area weighted sum. Suitable for flux fields.

      use types

      real (r8), dimension(nap), intent(in) :: a_area, a_field
      real (r8), dimension(nop), intent(in) :: o_area
      real (r8), dimension(nop), intent(out) :: o_field

      integer :: ioap, iop, iap

      o_field(:) = 0._r8
      do ioap = 1, noap
         iop = o_weight_index(ioap)
         iap = a_weight_index(ioap)
         o_field(iop) = o_field(iop) &
                      + a_field(iap)*a_area(iap)*a_weight_cnsrv(ioap)
      enddo
!$OMP PARALLEL DO
      do iop = 1, nop
         o_field(iop) = o_field(iop)/max(area_eps, o_area(iop))
      enddo
!$OMP END PARALLEL DO

   end subroutine remap_a2o_cnsrv

   subroutine remap_a2o_scalr(a_field, o_field)
      ! Remap a field on an atmosphere grid to an ocean grid. Suitable for
      ! scalar fields.

      use types

      real (r8), dimension(nap), intent(in) :: a_field
      real (r8), dimension(nop), intent(out) :: o_field

      integer :: ioap, iop, iap

      o_field(:) = 0._r8
      do ioap = 1, noap
         iop = o_weight_index(ioap)
         iap = a_weight_index(ioap)
         o_field(iop) = o_field(iop) + a_field(iap)*a_weight_scalr(ioap)
      enddo

   end subroutine remap_a2o_scalr

end module mo_gridremap
