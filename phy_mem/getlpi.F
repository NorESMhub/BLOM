! ------------------------------------------------------------------------------
! Copyright (C) 2008 Mats Bentsen
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

      real function getlpi(temp,saln,phiu,phil,pu)
c
c --- get lower pressure interface of a layer knowing the temperature,
c --- salinity of the layer and the geopotential at upper and lower
c --- interface
c
      implicit none
c
      real temp,saln,phiu,phil,pu
c
      real pl,q,dphi,alpu,alpl
c
      real rho
      external rho
c
c --- first guess on pressure interface
      pl=pu-rho(pu,temp,saln)*(phil-phiu)
c
c --- improve the accuracy of the pressure interface by an
c --- iterative procedure
      q=1.
      do while (abs(q).gt.1.e-4)
        call delphi(pu,pl,temp,saln,dphi,alpu,alpl)
        q=(phil-phiu-dphi)/alpl
        pl=pl-q
      enddo
c
      getlpi=pl
c
      return
      end
