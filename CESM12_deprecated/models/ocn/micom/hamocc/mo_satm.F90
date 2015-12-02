module mo_satm

   use types
   use ncutils
   use mod_xc
   use mo_gridremap
#ifdef __c_isotopes
   use mo_param1_bgc, only: natm, iatmco2, iatmo2, iatmn2, iatmc13, iatmc14
#else
   use mo_param1_bgc, only: natm, iatmco2, iatmo2, iatmn2
#endif

   implicit none

   integer, parameter :: &
      nax = 120, &
      nay =  60

   real (r8), parameter :: &
      difmer =  6666._r8, & ! Meridional diffusivity [m^2/s]
      difzon = 20000._r8    ! Zonal diffusivity [m^2/s]

   real (r8), dimension(itdm,jtdm) :: o_area
   real (r8), dimension(nax,nay) :: a_area
   real (r8), dimension(nax,nay,natm) :: acon_a
   real (r8), dimension(nay) :: a_area_i, udsc, vdsc, udrs, vdrs
   real (r8) :: dt
   integer :: nsubcyc

contains

   subroutine satm_init(o_area_patch, dt_in)
      ! Initialize slab atmosphere

      implicit none

      real (r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), intent(in) :: &
         o_area_patch
      real (r8), intent(in) :: dt_in

      real (r8), parameter :: &
         pi     = 3.14159265358979324e+00_r8, &  ! pi []
         rearth = 6.37122e6_r8                   ! NCAR CCCM constant [m]
      

      real (r8), dimension(nax*nay) :: a_area_1d
      real (r8), dimension(nay+1) :: scxv, scxp
      real (r8) :: scy, q, dx2, dy2, facsc, facrs
      integer ncid, i, j

      ! Gather the global ocean area
      call xcaget(o_area, o_area_patch, 1)

      ! Convert unit of ocean area from cm^2 to radians^2
      q = 1.e-4/(rearth*rearth)
      do j = 1, jtdm
         do i = 1, itdm
            o_area(i,j) = q*o_area(i,j)
         enddo
      enddo

      if (mnproc == 1) then

         ! Read atmosphere grid area
         call ncopen('gridremap.nc', ncid)
         call ncgetvar(ncid, 'b_poly_area', a_area_1d)
         call ncclose(ncid)
!$OMP PARALLEL DO
         do j = 1, nay
            do i = 1, nax
               a_area(i,j) = a_area_1d(i + (j - 1)*nax)
            enddo
         enddo
!$OMP END PARALLEL DO
         do j = 1, nay
            a_area_i(j) = 1._r8/a_area(1,j)
         enddo

         ! Initialize remapping between ocean and atmosphere grids
         call remap_init('gridremap.nc', itdm*jtdm, nax*nay)

         ! Time step
         dt = dt_in

         ! Compute scale factors and diffusion related coefficients
         scy = 2._r8*pi/nax
         do j = 2, nay
            scxv(j) = cos(.5_r8*(2*(j - 1) - nay)/real(nay, r8)*pi)*scy
            scxp(j) = cos(.5_r8*(2*(j - .5_r8) - nay)/real(nay, r8)*pi)*scy
         enddo
         q = 1._r8
         do j = 2, nay - 1
            dx2 = min(scxv(j), scxv(j+1))**2
            dy2 = scy*scy
            q = min(q, .5*dx2*dy2/(dx2+dy2))
         enddo
         facsc = q/(dt*max(difmer, difzon)/(rearth*rearth))
         nsubcyc = int(1._r8/facsc)
         facrs = 1._r8 - nsubcyc*facsc
         write (lp,*) 'satm_init: number of diffusion subcycles:', nsubcyc
         if (nsubcyc > 0) then
            write (lp,*) 'satm_init: diffusion factor for subcycles:', facsc
            write (lp,*) 'satm_init: diffusion factor for residual step:', facrs
         endif
         q = dt/(rearth*rearth)
         do j = 2, nay
            udsc(j) = facsc*q*difzon*scy/scxp(j)
            vdsc(j) = facsc*q*difmer*scxv(j)/scy
            udrs(j) = facrs*q*difzon*scy/scxp(j)
            vdrs(j) = facrs*q*difmer*scxv(j)/scy
         enddo

         ! Initialize concentration of atmospheric constituents
         do j = 1, nay
            do i = 1, nax
               acon_a(i,j,iatmco2) = 278.
               acon_a(i,j,iatmo2 ) = 196800.
               acon_a(i,j,iatmn2 ) = 802000.
#ifdef __c_isotopes
               acon_a(i,j,iatmc13) = 278.-(278.*0.0065)
               acon_a(i,j,iatmc14) = 278.-(278.*0.0065)**2
#endif
            enddo
         enddo

      endif

   end subroutine satm_init

   subroutine satm_step(aconflx_o, acon_o)
      ! Advance the slab atmosphere state one time step

      implicit none

      real (r8), dimension(idm,jdm,natm), intent(in) :: aconflx_o
      real (r8), dimension(idm,jdm,natm), intent(out) :: acon_o

      real (r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp_patch
      real (r8), dimension(itdm,jtdm) :: tmp_glob
      real (r8), dimension(nax,nay) :: flx_a, udf, vdf
      real (r8) :: sum_sp, sum_np, ave_sp, ave_np
      integer i, j, n, nsc, im1, ip1

!     write(lp,*) 'IN SATM_STEP',natm
      
      ! Loop trough all constituents
      do n = 1, natm

         ! Collect surface fluxes from grid patch to a global surface flux field
!$OMP PARALLEL DO
         do j = 1, jj
            do i = 1, ii
               tmp_patch(i,j) = aconflx_o(i,j,n)
            enddo
         enddo
!$OMP END PARALLEL DO
         call xcaget(tmp_glob, tmp_patch, 1)

         ! Update the atmospheric state on the master processor only
         if (mnproc == 1) then

            ! Remap the surface flux from the ocean to the atmospheric grid
            call remap_o2a_cnsrv(o_area, a_area, tmp_glob, flx_a)

            ! Update the concentration with the surface flux
!$OMP PARALLEL DO
            do j = 1, nay
               do i = 1, nax
                  acon_a(i,j,n) = acon_a(i,j,n) + flx_a(i,j)
!                  acon_a(i,j,n) = acon_a(i,j,n) + dt*flx_a(i,j)
               enddo
            enddo
!$OMP END PARALLEL DO

            ! Update the concentration with emissions

            ! Diffuse the concentration field
!            write (lp,*) n, sum(acon_a(:,:,n)*a_area(:,:))
            do nsc = 1, nsubcyc
!$OMP PARALLEL DO PRIVATE(im1)
               do j = 2, nay
                  do i = 1, nax
                     im1 = mod(i - 2 + nax, nax) + 1
                     udf(i,j) = udsc(j)*(acon_a(im1,j,n) - acon_a(i,j,n))
                     vdf(i,j) = vdsc(j)*(acon_a(i,j-1,n) - acon_a(i,j,n))
                  enddo
               enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(ip1)
               do j = 2, nay - 1
                  do i = 1, nax
                     ip1 = mod(i, nax) + 1 
                     acon_a(i,j,n) = acon_a(i,j,n) &
                                   - a_area_i(j)*( udf(ip1,j) - udf(i,j) &
                                                 + vdf(i,j+1) - vdf(i,j))
                  enddo
               enddo
!$OMP END PARALLEL DO
               sum_sp = 0._r8
               sum_np = 0._r8
               do i = 1, nax
                  sum_sp = sum_sp + acon_a(i,1  ,n) - a_area_i(1  )*vdf(i,2  )
                  sum_np = sum_np + acon_a(i,nay,n) + a_area_i(nay)*vdf(i,nay)
               enddo
               ave_sp = sum_sp/nax
               ave_np = sum_np/nax
               do i = 1, nax
                  acon_a(i,1  ,n) = ave_sp
                  acon_a(i,nay,n) = ave_np 
               enddo
            enddo
!$OMP PARALLEL DO PRIVATE(im1)
            do j = 2, nay
               do i = 1, nax
                  im1 = mod(i - 2 + nax, nax) + 1
                  udf(i,j) = udrs(j)*(acon_a(im1,j,n) - acon_a(i,j,n))
                  vdf(i,j) = vdrs(j)*(acon_a(i,j-1,n) - acon_a(i,j,n))
               enddo
            enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(ip1)
            do j = 2, nay - 1
               do i = 1, nax
                  ip1 = mod(i, nax) + 1 
                  acon_a(i,j,n) = acon_a(i,j,n) &
                                - a_area_i(j)*( udf(ip1,j) - udf(i,j) &
                                              + vdf(i,j+1) - vdf(i,j))
               enddo
            enddo
!$OMP END PARALLEL DO
            sum_sp = 0._r8
            sum_np = 0._r8
            do i = 1, nax
               sum_sp = sum_sp + acon_a(i,1  ,n) - a_area_i(1  )*vdf(i,2  )
               sum_np = sum_np + acon_a(i,nay,n) + a_area_i(nay)*vdf(i,nay)
            enddo
            ave_sp = sum_sp/nax
            ave_np = sum_np/nax
            do i = 1, nax
               acon_a(i,1  ,n) = ave_sp
               acon_a(i,nay,n) = ave_np 
            enddo
!            write (lp,*) n, sum(acon_a(:,:,n)*a_area(:,:))

            ! Remap the atmospheric concentration from the atmospheric to the
            ! ocean grid
            call remap_a2o_scalr(acon_a(1,1,n), tmp_glob)

         endif

         ! Distribute the global concentration field to the grid patches
         call xcaput(tmp_glob, tmp_patch, 1)
!$OMP PARALLEL DO
         do j = 1, jj
            do i = 1, ii
               acon_o(i,j,n) = tmp_patch(i,j)
            enddo
         enddo
!$OMP END PARALLEL DO

      enddo

      if (mnproc.eq.1) then
         open(10,file='satm.uf',form='unformatted')
         write(10) acon_a
         close(10)
      endif

   end subroutine satm_step

end module mo_satm
