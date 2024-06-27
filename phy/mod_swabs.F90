! ------------------------------------------------------------------------------
! Copyright (C) 2017-2024 Mats Bentsen, Mariana Vertenstein
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

  ! ------------------------------------------------------------------
  ! This module contains routines and specifies arrays related to
  ! shortwave radiation absorption.
  !
  ! It is assumed that the vertical profile of shortwave radiation
  ! flux is
  !
  !   E(z) = E(0)*swbgfc*exp(-z/swbgal)
  !
  ! where E(0) is the shortwave radiation flux immediately below the
  ! ocean surface, z is the depth, swbgfc and swbgal is the fraction
  ! and attenuation length, respectively, of ultraviolet and short
  ! visible wavelengths (blue/green). Thus it is assumed that the
  ! infrared and long visible (red) wavelengths of the radiation is
  ! absorbed in the uppermost model layer.
  !
  ! References:
  !   Jerlov, N. G., 1968: Optical Oceanography. Elsevier, 194 pp.
  !   Paulson, C. A., and J. J. Simpson, 1977: Irradiance Measurements
  !     in the Upper Ocean. J. Phys. Oceanogr., 7, 952-956.
  !   Morel, A., and D. Antoine, 1994: Heating Rate within the Upper
  !     Ocean in Relation to Its Bio-Optical State. J. Phys.
  !     Oceanogr., 24, 1652-1665.
  !   Sweeney, C., A. Gnanadesikan, S. M. Griffies, M. J. Harrison, A.
  !     J. Rosati, and B. L. Samuels, 2005: Impacts of Shortwave
  !     Penetration Depth on Large-Scale Ocean Circulation and Heat
  !     Transport. J. Phys. Oceanogr., 35, 1103-1118.
  ! ------------------------------------------------------------------

  use dimensions,   only: idm, jdm, itdm, jtdm
  use mod_xc,       only: xcstop, xchalt, xcaput, xctilr, &
                          nbdy, ii, jj, ip, isp, ifp, ilp, lp, &
                          halo_ps, mnproc
  use mod_time,     only: xmi, l1mi, l2mi, l3mi, l4mi, l5mi
  use mod_checksum, only: csdiag, chksummsk
  use mod_intp1d,   only: intp1d
  use netcdf

  implicit none
  private

  ! Variables to be set in namelist:
  !   swamth: Shortwave radiation absorption method. Valid methods:
  !           'top-layer'  : All radiation absorbed in top model
  !                          layer.
  !           'jerlov'     : Absorption using modified Paulson and
  !                          Simpson (1977) transmission
  !                          parameterization of Jerlov (1968) water
  !                          types.
  !           'chlorophyll': Absorption using modified Morel and
  !                          Antoine (1994) chlorophyll concentration
  !                          dependent transmission parameterization.
  !  jwtype: Number indicating the Jerlov (1968) water type.
  !  chlopt: Chlorophyll concentration option. Valid options:
  !          'climatology': Monthly chlorophyll concentration
  !                         climatology obtained from SeaWiFS
  !                         observations during 1997-2010.
  !  ccfile: Name of file containing chlorophyll concentration
  !          climatology.
  character (len = 80),  public :: swamth,chlopt
  character (len = 256), public :: ccfile
  integer,               public :: jwtype

  ! Parameter arrays related to shortwave radiation absorption
  ! following Paulson and Simpson's (1977) fit to the data of Jerlov
  ! (1968):
  !   ps77rf: The infrared/red fraction absorbed near surface.
  !   ps77al: Attenuation length for blue/green spectral band.
  ! Array index 1 through 5 correspond to Jerlov's classification of
  ! water types I, IA, IB, II and III, respectively.
  real, dimension(5), parameter :: &
       ps77rf = (/.58,.62,.67,.77,.78/), &
       ps77al = (/23.,20.,17.,14.,7.9/)

  ! Parameters related to a modified Morel and Antoine (1994)
  ! transmission parameterization:
  !   cl10mn: Minimum value of log10 of chlorophyll concentration.
  !   cl10mx: Maximum value of log10 of chlorophyll concentration.
  !   ms94rf: Infrared (750-2500 nm) fraction of shortwave radiation
  !           that is absorbed in the upper few meter of the ocean.
  !           The value 0.43 suggested by Sweeney et al. (2005) is
  !           used here.
  !   ma94v2: Coefficients of polynomial determining the fraction of
  !           ultraviolet and visible light (300-750 nm) that belongs
  !           to ultraviolet and short visible wavelengths
  !           (blue/green).
  !   ma94z2: Coefficients of polynomial determining the attenuation
  !           length of ultraviolet and short visible wavelengths
  !           (blue/green).
  real, parameter :: &
       cl10mn = -2., &
       cl10mx = 1., &
       ma94rf = .43
  real, dimension(6), parameter :: &
       ma94v2 = (/  .679, -.008, -.132, -.038,  .017,  .007/), &
       ma94z2 = (/ 7.925,-6.644, 3.662,-1.815, -.218,  .502/)

  ! Other parameters:
  !----   swamxd: Maximum depth of shortwave radiation penetration [m].
  real, parameter, public ::  swamxd = 200.

  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12)      :: chl10c
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)         :: chl10
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), public :: swbgal
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), public :: swbgfc

  ! Public routines
  public :: iniswa, updswa

contains

  ! ------------------------------------------------------------------

  subroutine iniswa()

    ! ------------------------------------------------------------------
    ! Initialize shortwave radiation absorption functionality.
    ! ------------------------------------------------------------------

    ! Local variables
    real, dimension(itdm,jtdm) :: tmp2d
    integer, dimension(3) :: istart,icount
    integer :: i,j,l,k,istat,ncid,dimid,varid

    if     (swamth == 'top-layer') then

      ! Set penetrative blue/green shortwave radiation fraction to zero
      ! to ensure all shortwave radiation is absorbed in the top layer.
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            swbgfc(i,j) = 0.
            swbgal(i,j) = swamxd
          end do
        end do
      end do
      !$omp end parallel do

    else if (swamth == 'jerlov') then

      ! Set blue/green penetrative shortwave radiation fraction and
      ! attenuation length according to Jerlov water type.
      if (jwtype < 1.or.jwtype > 5) then
        if (mnproc == 1) then
          write (lp,'(a,i11,a)') ' jwtype = ',jwtype, &
               ' is outside the valid interval of [1,5]!'
        end if
        call xcstop('(iniswa)')
        stop '(iniswa)'
      end if
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            swbgfc(i,j) = 1.-ps77rf(jwtype)
            swbgal(i,j) = ps77al(jwtype)
          end do
        end do
      end do
      !$omp end parallel do

    else if (swamth == 'chlorophyll') then

      ! Initialize functionality for chlorophyll concentration dependent
      ! shortwave radiation absorption.
      if     (chlopt == 'climatology') then

        ! Read monthly chlorophyll concentration climatology.
        if (mnproc == 1) then
          write (lp,'(2a)') &
               ' reading chlorophyll concentration climatology from ', &
               trim(ccfile)
          call flush(lp)

          ! ----- Open netCDF file.
          istat = nf90_open(ccfile,nf90_nowrite,ncid)
          if (istat /= nf90_noerr) then
            write(lp,'(4a)') ' nf90_open: ',trim(ccfile),': ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if

          ! Check dimensions.
          istat = nf90_inq_dimid(ncid,'x',dimid)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inq_dimid: x: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
          istat=nf90_inquire_dimension(ncid,dimid,len = i)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inquire_dimension: x: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
          istat = nf90_inq_dimid(ncid,'y',dimid)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inq_dimid: y: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
          istat=nf90_inquire_dimension(ncid,dimid,len = j)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inquire_dimension: y: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
          istat = nf90_inq_dimid(ncid,'time',dimid)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inq_dimid: time: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
          istat=nf90_inquire_dimension(ncid,dimid,len = k)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inquire_dimension: time: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
          if (i /= itdm.or.j /= jtdm.or.k /= 12) then
            write (lp,'(2a)') &
                 ' wrong dimensions in ',trim(ccfile)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if

          ! Get variable ID.
          istat = nf90_inq_varid(ncid,'chlor_a',varid)
          if (istat /= nf90_noerr) then
            write(lp,'(2a)') ' nf90_inq_varid: chlor_a: ', &
                 nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if

          ! Set start and count vectors for reading monthly slices of
          ! data.
          istart(1) = 1
          istart(2) = 1
          icount(1) = itdm
          icount(2) = jtdm
          icount(3) = 1

        end if

        ! Read data on master process and distribute to all processes.
        do k = 1,12
          if (mnproc == 1) then
            istart(3) = k
            istat = nf90_get_var(ncid,varid,tmp2d,istart,icount)
            if (istat /= nf90_noerr) then
              write(lp,'(2a)') ' nf90_get_var: chlor_a: ', &
                   nf90_strerror(istat)
              call xchalt('(iniswa)')
              stop '(iniswa)'
            end if
          end if
          call xcaput(tmp2d,chl10c(1-nbdy,1-nbdy,k),1)
        end do

        ! Close file.
        if (mnproc == 1) then
          istat = nf90_close(ncid)
          if (istat /= nf90_noerr) then
            write(lp,'(4a)') &
                 ' nf90_close: ',trim(ccfile),': ',nf90_strerror(istat)
            call xchalt('(iniswa)')
            stop '(iniswa)'
          end if
        end if

        ! Convert to log10 of chlorophyll concentration.
        !$omp parallel do private(k,l,i)
        do j = 1,jj
          do k = 1,12
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                chl10c(i,j,k) = log10(chl10c(i,j,k))
              end do
            end do
          end do
        end do
        !$omp end parallel do

        ! --- Make sure halos are updated.
        call xctilr(chl10c, 1,12, nbdy,nbdy, halo_ps)

      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' chlopt = ',trim(chlopt),' is unsupported!'
        end if
        call xcstop('(iniswa)')
        stop '(iniswa)'
      end if

    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' swamth = ',trim(swamth),' is unsupported!'
      end if
      call xcstop('(iniswa)')
      stop '(iniswa)'
    end if

  end subroutine iniswa

  ! --- ------------------------------------------------------------------

  subroutine updswa()

    ! ------------------------------------------------------------------
    ! Update arrays related to shortwave radiation absorption.
    ! ------------------------------------------------------------------

    ! Local variables
    integer :: i,j,l
    real :: q

    if     (swamth == 'top-layer'.or.swamth == 'jerlov') then

      ! Nothing to be done for methods 'top-layer' or 'jerlov'.
      return

    else if (swamth == 'chlorophyll') then

      ! Time interpolation of chlorophyll concentration climatology.
      if     (chlopt == 'climatology') then
        !$OMP PARALLEL DO PRIVATE(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              chl10(i,j) = intp1d(chl10c(i,j,l1mi),chl10c(i,j,l2mi), &
                   chl10c(i,j,l3mi),chl10c(i,j,l4mi), &
                   chl10c(i,j,l5mi),xmi)
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      else
        if (mnproc == 1) then
          write (lp,'(3a)') ' chlopt = ',trim(chlopt),' is unsupported!'
        end if
        call xcstop('(updswa)')
        stop '(updswa)'
      end if

      !  Compute blue/green penetrative shortwave radiation fraction and
      !  attenuation length according to modified Morel and Antoine
      !  (1994) scheme.
      !$OMP PARALLEL DO PRIVATE(l,i,q)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = max(cl10mn,min(cl10mx,chl10(i,j)))
            swbgfc(i,j) = (1.-ma94rf) &
                 *(   ma94v2(1)+q*(ma94v2(2)+q*(ma94v2(3) &
                 +q*(ma94v2(4)+q*(ma94v2(5)+q* ma94v2(6))))))
            swbgal(i,j)=   ma94z2(1)+q*(ma94z2(2)+q*(ma94z2(3) &
                 +q*(ma94z2(4)+q*(ma94z2(5)+q* ma94z2(6)))))
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' swamth = ',trim(swamth),' is unsupported!'
      end if
      call xcstop('(updswa)')
      stop '(updswa)'
    end if

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'updswa:'
      end if
      call chksummsk(swbgfc,ip,1,'swbgfc')
      call chksummsk(swbgal,ip,1,'swbgal')
    end if

  end subroutine updswa

end module mod_swabs
