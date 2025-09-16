! ------------------------------------------------------------------------------
! Copyright (C) 2017-2025 Mats Bentsen, Mariana Vertenstein, Mehmet Ilicak
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_swabs
! ----------------------------------------------------------------------------
! This module contains routines and specifies arrays related to shortwave
! radiation absorption.
!
! It is assumed that the vertical profile of shortwave radiation flux is
!
!   E(z) = E(0)*(swfc1*exp(-z/swal1) + swfc2*exp(-z/swal2))
!
! where E(0) is the shortwave radiation flux immediately below the ocean
! surface, z is the depth, swfc1 and swal1 are the fraction and attenuation
! length, respectively, of long visible visible wavelengths (also infrared
! wavelengths in case of Jerlov (1968) water types) and swfc2 and swal2 are the
! fraction and attenuation length of and short visible and ultraviolet
! wavelengths. For vcoord_type = 'isopyc_bulkml' it is assumed that the infrared
! and long visible wavelengths of the radiation is absorbed in the uppermost
! model layer.
!
! References:
!   Jerlov, N. G., 1968: Optical Oceanography. Elsevier, 194 pp.
!   Paulson, C. A., and J. J. Simpson, 1977: Irradiance Measurements in the
!     Upper Ocean. J. Phys. Oceanogr., 7, 952-956.
!   Morel, A., and D. Antoine, 1994: Heating Rate within the Upper Ocean in
!     Relation to Its Bio-Optical State. J. Phys. Oceanogr., 24, 1652-1665.
!   Sweeney, C., A. Gnanadesikan, S. M. Griffies, M. J. Harrison, A. J. Rosati,
!     and B. L. Samuels, 2005: Impacts of Shortwave Penetration Depth on
!     Large-Scale Ocean Circulation and Heat Transport. J. Phys. Oceanogr., 35,
!     1103-1118.
! ------------------------------------------------------------------------------

   use mod_types,    only: r8
   use dimensions,   only: idm, jdm, itdm, jtdm
   use mod_xc,       only: xcstop, xchalt, xcaput, xctilr, &
                           nbdy, ii, jj, ip, isp, ifp, ilp, lp, &
                           halo_ps, mnproc
   use mod_time,     only: xmi, l1mi, l2mi, l3mi, l4mi, l5mi
   use mod_checksum, only: csdiag, chksum
   use mod_intp1d,   only: intp1d
   use mod_utility,  only: fnmlen
   use netcdf

   implicit none
   private

   ! Variables to be set in namelist:
   !    swamth: Shortwave radiation absorption method. Valid methods:
   !            'top-layer'          : All radiation absorbed in top model
   !                                   layer.
   !            'jerlov'             : Absorption using modified Paulson and
   !                                   Simpson (1977) transmission
   !                                   parameterization of Jerlov (1968) water
   !                                   types.
   !            'chlorophyll'        : Absorption using modified Morel and
   !                                   Antoine (1994) chlorophyll concentration
   !                                   dependent transmission parameterization.
   !            'spatial_frac_attlen': Read spatially varying spectral band
   !                                   fractions and attenuation lengths from
   !                                   file.
   !   jwtype: Number indicating the Jerlov (1968) water type.
   !   chlopt: Chlorophyll concentration option. Valid options:
   !           'climatology': Monthly chlorophyll concentration climatology
   !                          obtained from SeaWiFS observations during
   !                          1997-2010.
   !   ccfile: Name of file containing chlorophyll concentration climatology.
   !   svfile: Name of file containing spatially varying spectral band fractions
   !           and attenuation lengths.
   character (len = 80)     :: swamth, chlopt
   character (len = fnmlen) :: ccfile, svfile
   integer                  :: jwtype

  ! Parameter arrays related to shortwave radiation absorption following Paulson
  ! and Simpson's (1977) fit to the data of Jerlov (1968):
  !    ps77_irfc: Fraction of infrared and long visible wavelengths.
  !    ps77_al1 : Attenuation length for infrared and long visible wavelengths.
  !    ps77_al2 : Attenuation length for short visible and ultraviolet
  !               wavelengths.
  ! Array index 1 through 5 correspond to Jerlov's classification of water types
  ! I, IA, IB, II and III, respectively.
  real(r8), dimension(5), parameter :: &
     ps77_irfc = [  .58_r8,   .62_r8,   .67_r8,   .77_r8,   .78_r8], &
     ps77_al1  = [  .35_r8,   .60_r8,  1.00_r8,  1.50_r8,  1.40_r8], &
     ps77_al2  = [23.00_r8, 20.00_r8, 17.00_r8, 14.00_r8,  7.90_r8]

  ! Parameters related to a modified Morel and Antoine (1994) transmission
  ! parameterization:
  !    chl10_min: Minimum value of log10 of chlorophyll concentration.
  !    chl10_max: Maximum value of log10 of chlorophyll concentration.
  !    ma94_irfc: Infrared fraction of shortwave radiation that is absorbed near
  !               the surface. The value 0.43 suggested by Sweeney et al. (2005)
  !               is used here.
  !    ma94_v2  : Coefficients of polynomial determining the fraction of short
  !               visible and ultraviolet wavelengths.
  !    ma94_z1  : Coefficients of polynomial determining the attenuation length
  !               of long visible wavelengths.
  !    ma94_z2  : Coefficients of polynomial determining the attenuation length
  !               of short visible and ultraviolet wavelengths.
  real(r8), parameter :: &
       chl10_min = -2._r8, &
       chl10_max = 1._r8, &
       ma94_irfc = .43_r8
  real(r8), dimension(6), parameter :: &
       ma94_v2 = [  .679_r8,  -.008_r8,  -.132_r8, &
                   -.038_r8,   .017_r8,   .007_r8], &
       ma94_z1 = [ 1.540_r8,  -.197_r8,   .166_r8, &
                   -.252_r8,  -.055_r8,   .042_r8], &
       ma94_z2 = [ 7.925_r8, -6.644_r8,  3.662_r8, &
                  -1.815_r8,  -.218_r8,   .502_r8]

  ! Other parameters:
  !    swamxd: Maximum depth of shortwave radiation penetration [m].
  real(r8), parameter ::  swamxd = 200._r8

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) :: chl10c
  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
     chl10, swfc1, swfc2, swal1, swal2

  public :: swamth, chlopt, ccfile, svfile, jwtype, &
            swamxd, swfc1, swfc2, swal1, swal2, &
            iniswa, updswa

contains

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine iniswa()
   ! ---------------------------------------------------------------------------
   ! Initialize shortwave radiation absorption functionality.
   ! ---------------------------------------------------------------------------

      ! Local variables
      real(r8), dimension(itdm,jtdm) :: tmp2d
      integer, dimension(3) :: istart, icount
      integer :: i, j, l, k, errstat, ncid, dimid, varid

      select case (trim(swamth))

         case ('top-layer')

            ! Set penetrative shortwave radiation fractions to zero to ensure
            ! all shortwave radiation is absorbed in the top layer.

            !$omp parallel do private(l, i)
            do j = 1, jj
               do l = 1, isp(j)
               do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                  swfc1(i,j) = 0._r8
                  swfc2(i,j) = 0._r8
                  swal1(i,j) = swamxd
                  swal2(i,j) = swamxd
               enddo
               enddo
            enddo
            !$omp end parallel do

         case ('jerlov')

            ! Set penetrative shortwave radiation fractions and attenuation
            ! lengths according to Jerlov water type.

            if (jwtype < 1 .or. jwtype > 5) then
               if (mnproc == 1) then
                  write (lp,'(a,i11,a)') ' jwtype = ', jwtype, &
                     ' is outside the valid interval of [1,5]!'
               endif
               call xcstop('(iniswa)')
                      stop '(iniswa)'
            endif
            !$omp parallel do private(l, i)
            do j = 1, jj
               do l = 1, isp(j)
               do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                  swfc1(i,j) = ps77_irfc(jwtype)
                  swfc2(i,j) = 1._r8 - ps77_irfc(jwtype)
                  swal1(i,j) = ps77_al1(jwtype)
                  swal2(i,j) = ps77_al2(jwtype)
               enddo
               enddo
            enddo
            !$omp end parallel do

         case ('chlorophyll')

            ! Initialize functionality for chlorophyll concentration dependent
            ! shortwave radiation absorption.

            if     (chlopt == 'climatology') then

               ! Read monthly chlorophyll concentration climatology.

               if (mnproc == 1) then
                  write (lp,*) &
                     'reading chlorophyll concentration climatology from '// &
                     trim(ccfile)
                  call flush(lp)

                  ! Open netCDF file.
                  errstat = nf90_open(ccfile, nf90_nowrite, ncid)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_open: '//trim(ccfile)//': '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif

                  ! Check dimensions.
                  errstat = nf90_inq_dimid(ncid, 'x', dimid)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inq_dimid: x: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
                  errstat=nf90_inquire_dimension(ncid, dimid, len = i)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inquire_dimension: x: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
                  errstat = nf90_inq_dimid(ncid, 'y', dimid)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inq_dimid: y: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
                  errstat=nf90_inquire_dimension(ncid, dimid, len = j)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inquire_dimension: y: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
                  errstat = nf90_inq_dimid(ncid, 'time', dimid)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inq_dimid: time: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
                  errstat=nf90_inquire_dimension(ncid, dimid, len = k)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inquire_dimension: time: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
                  if (i /= itdm .or. j /= jtdm .or. k /= 12) then
                     write (lp,*) 'wrong dimensions in '//trim(ccfile)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif

                  ! Get variable ID.
                  errstat = nf90_inq_varid(ncid, 'chlor_a', varid)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_inq_varid: chlor_a: '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif

                  ! Set start and count vectors for reading monthly slices of
                  ! data.
                  istart(1) = 1
                  istart(2) = 1
                  icount(1) = itdm
                  icount(2) = jtdm
                  icount(3) = 1

               endif

               ! Read data on master process and distribute to all processes.
               do k = 1, 12
                  if (mnproc == 1) then
                     istart(3) = k
                     errstat = nf90_get_var(ncid, varid, tmp2d, istart, icount)
                     if (errstat /= nf90_noerr) then
                        write(lp,*) 'nf90_get_var: chlor_a: '// &
                                    nf90_strerror(errstat)
                        call xchalt('(iniswa)')
                               stop '(iniswa)'
                     endif
                  endif
                  call xcaput(tmp2d, chl10c(1-nbdy,1-nbdy,k), 1)
               enddo

               ! Close file.
               if (mnproc == 1) then
                  errstat = nf90_close(ncid)
                  if (errstat /= nf90_noerr) then
                     write(lp,*) 'nf90_close: '//trim(ccfile)//': '// &
                                 nf90_strerror(errstat)
                     call xchalt('(iniswa)')
                            stop '(iniswa)'
                  endif
               endif

               ! Convert to log10 of chlorophyll concentration.
               !$omp parallel do private(k, l, i)
               do j = 1, jj
                  do k = 1, 12
                     do l = 1, isp(j)
                     do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                        chl10c(i,j,k) = log10(chl10c(i,j,k))
                     enddo
                     enddo
                  enddo
               enddo
               !$omp end parallel do

               ! Make sure halos are updated.
               call xctilr(chl10c, 1, 12, nbdy, nbdy, halo_ps)

            else
               if (mnproc == 1) then
                  write (lp,*) 'iniswa: chlopt = '//trim(chlopt)// &
                               ' is unsupported!'
               endif
               call xcstop('(iniswa)')
                      stop '(iniswa)'
            endif

         case ('spatial_frac_attlen')

            ! Initialize functionality for spatially varying spectral band
            ! fractions and attenuation lengths from file.

            if (mnproc == 1) then
               write (lp,*) &
                  'reading spatially varying spectral band fractions and '// &
                  'attenuation lengths from '//trim(svfile)
               call flush(lp)

               ! Open netCDF file.
               errstat = nf90_open(svfile, nf90_nowrite, ncid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_open: '//trim(svfile)//': '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif

               ! Check dimensions.
               errstat = nf90_inq_dimid(ncid, 'x', dimid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inq_dimid: x: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat=nf90_inquire_dimension(ncid, dimid, len = i)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inquire_dimension: x: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat = nf90_inq_dimid(ncid, 'y', dimid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inq_dimid: y: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat=nf90_inquire_dimension(ncid, dimid, len = j)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inquire_dimension: y: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               if (i /= itdm .or. j /= jtdm) then
                  write (lp,*) 'wrong dimensions in '//trim(svfile)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif

            endif

            ! Read fraction of infrared and long visible wavelengths.
            if (mnproc == 1) then
               errstat = nf90_inq_varid(ncid, 'swfc1', varid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inq_varid: swfc1: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat = nf90_get_var(ncid, varid, tmp2d)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_get_var: swfc1: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
            endif
            call xcaput(tmp2d, swfc1, 1)
            call xctilr(swfc1, 1, 1, nbdy, nbdy, halo_ps)

            ! Read fraction of short visible and ultraviolet wavelengths.
            if (mnproc == 1) then
               errstat = nf90_inq_varid(ncid, 'swfc2', varid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inq_varid: swfc2: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat = nf90_get_var(ncid, varid, tmp2d)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_get_var: swfc2: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
            endif
            call xcaput(tmp2d, swfc2, 1)
            call xctilr(swfc2, 1, 1, nbdy, nbdy, halo_ps)

            ! Read attenuation length of infrared and long visible wavelengths.
            if (mnproc == 1) then
               errstat = nf90_inq_varid(ncid, 'swal1', varid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inq_varid: swal1: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat = nf90_get_var(ncid, varid, tmp2d)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_get_var: swal1: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
            endif
            call xcaput(tmp2d, swal1, 1)
            call xctilr(swal1, 1, 1, nbdy, nbdy, halo_ps)

            ! Read attenuation length of short visible and ultraviolet
            ! wavelengths.
            if (mnproc == 1) then
               errstat = nf90_inq_varid(ncid, 'swal2', varid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_inq_varid: swal2: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
               errstat = nf90_get_var(ncid, varid, tmp2d)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_get_var: swal2: '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
            endif
            call xcaput(tmp2d, swal2, 1)
            call xctilr(swal2, 1, 1, nbdy, nbdy, halo_ps)

            ! Close file.
            if (mnproc == 1) then
               errstat = nf90_close(ncid)
               if (errstat /= nf90_noerr) then
                  write(lp,*) 'nf90_close: '//trim(svfile)//': '// &
                              nf90_strerror(errstat)
                  call xchalt('(iniswa)')
                         stop '(iniswa)'
               endif
            endif

         case default
            if (mnproc == 1) &
               write (lp,'(3a)') ' iniswa: swamth = ', trim(swamth), &
                                 ' is unsupported!'
            call xcstop('(iniswa)')
                   stop '(iniswa)'
      end select

   end subroutine iniswa

   subroutine updswa()
   ! ---------------------------------------------------------------------------
   ! Update arrays related to shortwave radiation absorption.
   ! ---------------------------------------------------------------------------

      ! Local variables
      integer :: i, j, l
      real(r8) :: q, v2

      select case (trim(swamth))

         case ('top-layer', 'jerlov', 'spatial_frac_attlen')

            ! Nothing to be done for methods 'top-layer', 'jerlov' or
            ! 'spatial_frac_attlen'.

            return

         case ('chlorophyll')

            if     (chlopt == 'climatology') then
               ! Time interpolation of chlorophyll concentration climatology.
               !$omp parallel do private(l, i)
               do j = 1, jj
                  do l = 1, isp(j)
                  do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                     chl10(i,j) = intp1d(chl10c(i,j,l1mi), chl10c(i,j,l2mi), &
                                         chl10c(i,j,l3mi), chl10c(i,j,l4mi), &
                                         chl10c(i,j,l5mi), xmi)
                  enddo
                  enddo
               enddo
               !$omp end parallel do
            else
               if (mnproc == 1) then
                  write (lp,'(3a)') ' updswa: chlopt = ', trim(chlopt), &
                                    ' is unsupported!'
               endif
               call xcstop('(updswa)')
                      stop '(updswa)'
            endif

            ! Compute penetrative shortwave radiation fractions and attenuation
            ! lengths according to modified Morel and Antoine (1994) scheme.
            !$omp parallel do private(l, i, q, v2)
            do j = 1, jj
              do l = 1, isp(j)
                do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                  q = max(chl10_min, min(chl10_max, chl10(i,j)))
                  v2 = ((((( ma94_v2(6) *q &
                           + ma94_v2(5))*q &
                           + ma94_v2(4))*q &
                           + ma94_v2(3))*q &
                           + ma94_v2(2))*q &
                           + ma94_v2(1))
                  swfc1(i,j) = (1._r8 - ma94_irfc)*(1._r8 - v2)
                  swfc2(i,j) = (1._r8 - ma94_irfc)*v2
                  swal1(i,j) = (((( ma94_z1(6) *q &
                                  + ma94_z1(5))*q &
                                  + ma94_z1(4))*q &
                                  + ma94_z1(3))*q &
                                  + ma94_z1(2))*q &
                                  + ma94_z1(1)
                  swal2(i,j) = (((( ma94_z2(6) *q &
                                  + ma94_z2(5))*q &
                                  + ma94_z2(4))*q &
                                  + ma94_z2(3))*q &
                                  + ma94_z2(2))*q &
                                  + ma94_z2(1)
                enddo
              enddo
            enddo
            !$omp end parallel do

         case default
            if (mnproc == 1) &
               write (lp,'(3a)') ' updswa: swamth = ', trim(swamth), &
                                 ' is unsupported!'
            call xcstop('(updswa)')
                   stop '(updswa)'
      end select

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'updswa:'
         endif
         call chksum(swfc1, 1, halo_ps, 'swfc1')
         call chksum(swfc2, 1, halo_ps, 'swfc2')
         call chksum(swal1, 1, halo_ps, 'swal1')
         call chksum(swal2, 1, halo_ps, 'swal2')
      endif

   end subroutine updswa

end module mod_swabs
