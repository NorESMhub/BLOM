! ------------------------------------------------------------------------------
! Copyright (C) 2017-2020 Mats Bentsen
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
c
c --- ------------------------------------------------------------------
c --- This module contains routines and specifies arrays related to
c --- shortwave radiation absorption. 
c ---
c --- It is assumed that the vertical profile of shortwave radiation
c --- flux is
c ---
c ---   E(z) = E(0)*swbgfc*exp(-z/swbgal)
c ---
c --- where E(0) is the shortwave radiation flux immediately below the
c --- ocean surface, z is the depth, swbgfc and swbgal is the fraction
c --- and attenuation length, respectively, of ultraviolet and short
c --- visible wavelengths (blue/green). Thus it is assumed that the
c --- infrared and long visible (red) wavelengths of the radiation is
c --- absorbed in the uppermost model layer.
c ---
c --- References:
c ---   Jerlov, N. G., 1968: Optical Oceanography. Elsevier, 194 pp.
c ---   Paulson, C. A., and J. J. Simpson, 1977: Irradiance Measurements
c ---     in the Upper Ocean. J. Phys. Oceanogr., 7, 952-956.
c ---   Morel, A., and D. Antoine, 1994: Heating Rate within the Upper
c ---     Ocean in Relation to Its Bio-Optical State. J. Phys.
c ---     Oceanogr., 24, 1652-1665.
c ---   Sweeney, C., A. Gnanadesikan, S. M. Griffies, M. J. Harrison, A.
c ---     J. Rosati, and B. L. Samuels, 2005: Impacts of Shortwave
c ---     Penetration Depth on Large-Scale Ocean Circulation and Heat
c ---     Transport. J. Phys. Oceanogr., 35, 1103-1118.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      private
c
c --- Variables to be set in namelist:
c ---   swamth: Shortwave radiation absorption method. Valid methods:
c ---           'top-layer'  : All radiation absorbed in top model
c ---                          layer.
c ---           'jerlov'     : Absorption using modified Paulson and
c ---                          Simpson (1977) transmission
c ---                          parameterization of Jerlov (1968) water
c ---                          types.
c ---           'chlorophyll': Absorption using modified Morel and
c ---                          Antoine (1994) chlorophyll concentration
c ---                          dependent transmission parameterization.
c ---  jwtype: Number indicating the Jerlov (1968) water type.
c ---  chlopt: Chlorophyll concentration option. Valid options:
c ---          'climatology': Monthly chlorophyll concentration
c ---                         climatology obtained from SeaWiFS
c ---                         observations during 1997-2010.
c ---  ccfile: Name of file containing chlorophyll concentration
c ---          climatology.
      character (len=80) :: swamth,chlopt
      character (len=256) :: ccfile
      integer :: jwtype
c
c --- Parameter arrays related to shortwave radiation absorption
c --- following Paulson and Simpson's (1977) fit to the data of Jerlov
c --- (1968):
c ---   ps77rf: The infrared/red fraction absorbed near surface.
c ---   ps77al: Attenuation length for blue/green spectral band.
c --- Array index 1 through 5 correspond to Jerlov's classification of
c --- water types I, IA, IB, II and III, respectively.
      real, dimension(5), parameter ::
     .  ps77rf=(/.58,.62,.67,.77,.78/),
     .  ps77al=(/23.,20.,17.,14.,7.9/)
c
c --- Parameters related to a modified Morel and Antoine (1994)
c --- transmission parameterization:
c ---   cl10mn: Minimum value of log10 of chlorophyll concentration.
c ---   cl10mx: Maximum value of log10 of chlorophyll concentration.
c ---   ms94rf: Infrared (750-2500 nm) fraction of shortwave radiation
c ---           that is absorbed in the upper few meter of the ocean.
c ---           The value 0.43 suggested by Sweeney et al. (2005) is
c ---           used here.
c ---   ma94v2: Coefficients of polynomial determining the fraction of
c ---           ultraviolet and visible light (300-750 nm) that belongs
c ---           to ultraviolet and short visible wavelengths
c ---           (blue/green).
c ---   ma94z2: Coefficients of polynomial determining the attenuation
c ---           length of ultraviolet and short visible wavelengths
c ---           (blue/green).
      real, parameter ::
     .  cl10mn=-2.,
     .  cl10mx=1.,
     .  ma94rf=.43
      real, dimension(6), parameter ::
     .  ma94v2=(/  .679, -.008, -.132, -.038,  .017,  .007/),
     .  ma94z2=(/ 7.925,-6.644, 3.662,-1.815, -.218,  .502/)
c
c --- Other parameters:
c----   swamxd: Maximum depth of shortwave radiation penetration.
      real, parameter ::
     .  swamxd = 200.
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) :: chl10c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  swbgal,swbgfc,chl10
c
      public :: swamth,jwtype,chlopt,ccfile,iniswa,updswa,swbgal,swbgfc,
     .          swamxd
c
      contains
c
c --- ------------------------------------------------------------------
c
      subroutine iniswa
c
c --- ------------------------------------------------------------------
c --- Initialize shortwave radiation absorption functionality.
c --- ------------------------------------------------------------------
c
      use netcdf
c
      implicit none
c
#include "common_blocks.h"
c
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(3) :: istart,icount
      integer :: i,j,l,k,istat,ncid,dimid,varid
c
      if     (swamth.eq.'top-layer') then
c
c --- - Set penetrative blue/green shortwave radiation fraction to zero
c --- - to ensure all shortwave radiation is absorbed in the top layer.
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            swbgfc(i,j)=0.
            swbgal(i,j)=swamxd
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (swamth.eq.'jerlov') then
c
c --- - Set blue/green penetrative shortwave radiation fraction and
c --- - attenuation length according to Jerlov water type.
        if (jwtype.lt.1.or.jwtype.gt.5) then
          if (mnproc.eq.1) then
            write (lp,'(a,i11,a)') ' jwtype=',jwtype,
     .        ' is outside the valid interval of [1,5]!'
          endif
          call xcstop('(iniswa)')
                 stop '(iniswa)'
        endif
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            swbgfc(i,j)=1.-ps77rf(jwtype)
            swbgal(i,j)=ps77al(jwtype)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (swamth.eq.'chlorophyll') then
c
c --- - Initialize functionality for chlorophyll concentration dependent
c --- - shortwave radiation absorption.
        if     (chlopt.eq.'climatology') then
c
c --- --- Read monthly chlorophyll concentration climatology.
          if (mnproc.eq.1) then
            write (lp,'(2a)')
     .        ' reading chlorophyll concentration climatology from ',
     .        trim(ccfile)
            call flush(lp)
c
c --- ----- Open netCDF file.
            istat=nf90_open(ccfile,nf90_nowrite,ncid)
            if (istat.ne.nf90_noerr) then
              write(lp,'(4a)') ' nf90_open: ',trim(ccfile),': ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
c
c --- ----- Check dimensions.
            istat=nf90_inq_dimid(ncid,'x',dimid)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inq_dimid: x: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
            istat=nf90_inquire_dimension(ncid,dimid,len=i)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inquire_dimension: x: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
            istat=nf90_inq_dimid(ncid,'y',dimid)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inq_dimid: y: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
            istat=nf90_inquire_dimension(ncid,dimid,len=j)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inquire_dimension: y: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
            istat=nf90_inq_dimid(ncid,'time',dimid)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inq_dimid: time: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
            istat=nf90_inquire_dimension(ncid,dimid,len=k)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inquire_dimension: time: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
            if (i.ne.itdm.or.j.ne.jtdm.or.k.ne.12) then
              write (lp,'(2a)')
     .          ' wrong dimensions in ',trim(ccfile)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
c
c --- ----- Get variable ID.
            istat=nf90_inq_varid(ncid,'chlor_a',varid)
            if (istat.ne.nf90_noerr) then
              write(lp,'(2a)') ' nf90_inq_varid: chlor_a: ',
     .                         nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
c
c --- ----- Set start and count vectors for reading monthly slices of
c --- ----- data.
            istart(1)=1
            istart(2)=1
            icount(1)=itdm
            icount(2)=jtdm
            icount(3)=1
c
          endif
c
c --- --- Read data on master process and distribute to all processes.
          do k=1,12
            if (mnproc.eq.1) then
              istart(3)=k
              istat=nf90_get_var(ncid,varid,tmp2d,istart,icount)
              if (istat.ne.nf90_noerr) then
                write(lp,'(2a)') ' nf90_get_var: chlor_a: ',
     .                           nf90_strerror(istat)
                call xchalt('(iniswa)')
                       stop '(iniswa)'
              endif
            endif
            call xcaput(tmp2d,chl10c(1-nbdy,1-nbdy,k),1)
          enddo
c
c --- --- Close file.
          if (mnproc.eq.1) then
            istat=nf90_close(ncid)
            if (istat.ne.nf90_noerr) then
              write(lp,'(4a)')
     .          ' nf90_close: ',trim(ccfile),': ',nf90_strerror(istat)
              call xchalt('(iniswa)')
                     stop '(iniswa)'
            endif
          endif
c
c --- --- Convert to log10 of chlorophyll concentration.
c$OMP PARALLEL DO PRIVATE(k,l,i)
          do j=1,jj
            do k=1,12
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                chl10c(i,j,k)=log10(chl10c(i,j,k))
              enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
c
c --- --- Make sure halos are updated.
          call xctilr(chl10c, 1,12, nbdy,nbdy, halo_ps)
c
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' chlopt=',trim(chlopt),' is unsupported!'
          endif
          call xcstop('(iniswa)')
                 stop '(iniswa)'
        endif
c
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' swamth=',trim(swamth),' is unsupported!'
        endif
        call xcstop('(iniswa)')
               stop '(iniswa)'
      endif
c
      end subroutine iniswa
c
c --- ------------------------------------------------------------------
c
      subroutine updswa
c
c --- ------------------------------------------------------------------
c --- Update arrays related to shortwave radiation absorption.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c
      integer :: i,j,l
      real :: q
c
      real :: intp1d
      external :: intp1d
c
      if     (swamth.eq.'top-layer'.or.swamth.eq.'jerlov') then
c
c --- - Nothing to be done for methods 'top-layer' or 'jerlov'.
        return
c
      elseif (swamth.eq.'chlorophyll') then
c
c --- - Time interpolation of chlorophyll concentration climatology.
        if     (chlopt.eq.'climatology') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              chl10(i,j)=intp1d(chl10c(i,j,l1mi),chl10c(i,j,l2mi),
     .                          chl10c(i,j,l3mi),chl10c(i,j,l4mi),
     .                          chl10c(i,j,l5mi),xmi)
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
c
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' chlopt=',trim(chlopt),' is unsupported!'
          endif
          call xcstop('(updswa)')
                 stop '(updswa)'
        endif
c
c --- - Compute blue/green penetrative shortwave radiation fraction and
c --- - attenuation length according to modified Morel and Antoine
c --- - (1994) scheme.
c$OMP PARALLEL DO PRIVATE(l,i,q)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=max(cl10mn,min(cl10mx,chl10(i,j)))
            swbgfc(i,j)=(1.-ma94rf)
     .                *(   ma94v2(1)+q*(ma94v2(2)+q*(ma94v2(3)
     .                 +q*(ma94v2(4)+q*(ma94v2(5)+q* ma94v2(6))))))
            swbgal(i,j)=   ma94z2(1)+q*(ma94z2(2)+q*(ma94z2(3)
     .                 +q*(ma94z2(4)+q*(ma94z2(5)+q* ma94z2(6)))))
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' swamth=',trim(swamth),' is unsupported!'
        endif
        call xcstop('(updswa)')
               stop '(updswa)'
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'updswa:'
        endif
        call chksummsk(swbgfc,ip,1,'swbgfc')
        call chksummsk(swbgal,ip,1,'swbgal')
      endif
c
      end subroutine updswa
c
      end module mod_swabs
