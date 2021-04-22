! ------------------------------------------------------------------------------
! Copyright (C) 2006-2020 Mats Bentsen
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

      module mod_temmin
c
c --- ------------------------------------------------------------------
c --- This module contains variables and procedures related to defining
c --- the lower bound of physical temperature values in isopycnic
c --- layers.
c --- ------------------------------------------------------------------
c
      use mod_types, only: r8
      use mod_config, only: expcnf
      use mod_grid, only: sigmar
      use mod_xc
      use mod_eos, only: ap11, ap12, ap13, ap14, ap15, ap16,
     .                   ap21, ap22, ap23, ap24, ap25, ap26, atf
      use mod_pointtest, only: itest, jtest, ptest
c
      implicit none
c
      private
c
      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     .  temmin ! Minimum temperature allowed in an isopycnic layer
               ! [deg C].
c
      public :: temmin, settemmin
c
      contains
c
c --- ------------------------------------------------------------------
c
      subroutine settemmin
c
c --- ------------------------------------------------------------------
c --- Set minimum physical temperature values in isopycnic layers
c --- ------------------------------------------------------------------
c
      integer i,j,k,l
      real salfrz,a,b,c
c
      if     (expcnf.eq.'cesm' .or. expcnf.eq.'single_column') then
c
c --- - Set temmin to a constant freezing temperature for all layers
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              temmin(i,j,k)=-3.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn'.or.
     .        expcnf.eq.'fuk95'.or.expcnf.eq.'channel') then
c
c --- - Let temmin be the freezing temperature of a given potential
c --- - density. This can be achieved by using potential density given
c --- - in the function sig and the salinity dependent freezing
c --- - temperature given in the function swtfrz.
c
c$OMP PARALLEL DO PRIVATE(k,l,i,a,b,c,salfrz)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              a=((ap14-ap24*sigmar(i,j,k))*atf
     .          + ap15-ap25*sigmar(i,j,k) )*atf
     .         +ap16-ap26*sigmar(i,j,k)
              b=(ap12-ap22*sigmar(i,j,k))*atf+ap13-ap23*sigmar(i,j,k)
              c=ap11-ap21*sigmar(i,j,k)
              salfrz=(-b+sqrt(b*b-4.*a*c))/(2.*a)
              temmin(i,j,k)=atf*salfrz
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (expcnf.eq.'isomip1'.or.expcnf.eq.'isomip2') then
c
c --- - Set temmin to a low value.
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              temmin(i,j,k)=-10.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(settemmin)')
               stop '(settemmin)'
      endif
c
      if (mnproc.eq.ptest) then
        write (lp,'(a/(6(i5,f8.3)))') 'minimum temperature values:',
     .    (k,temmin(itest,jtest,k),k=2,kk)
      endif
c
      end subroutine settemmin
c
      end module mod_temmin
