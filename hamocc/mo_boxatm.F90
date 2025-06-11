! Copyright (C) 2018-2019  A. Moree
! Copyright (C) 2023  J. Schwinger
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.


module mo_boxatm
  !*************************************************************************************************
  ! This module contains the routine update_boxatm for updating a 1-D/scalar/box atmosphere
  !
  ! The global sum of the air-sea C fluxes is calculated, then converted to ppm
  ! and added to the global atmospheric concentration. For C14, an atmospheric
  ! production term corresponding to the total decay in the ocean (plus sediment
  ! if activated) is assumed.
  !
  ! A. Moree,            *GFI, Bergen*      Oct 2019
  !
  ! Modified
  ! A. Moree,            *GFI, Bergen*      2019-10
  ! - 14C source added to atmosphere as the sum of all 14C loss (decay)
  ! J. Schwinger,        *NORCE, Bergen*    2023-08-02
  ! - ported into NorESM2 code, no functional changes
  !*************************************************************************************************

  implicit none
  private

  public :: update_boxatm

contains

  subroutine update_boxatm(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask)

    use mod_xc,         only: mnproc,nbdy,ips,xcsum
    use mo_kind,        only: rp
    use mo_control_bgc, only: io_stdo_bgc, use_cisonew, use_sedbypass
    use mo_carbch,      only: atmflx, atm, ocetra
    use mo_param_bgc,   only: rcar,c14dec
    use mo_param1_bgc,  only: iatmco2,iatmc13,iatmc14,isco214,idet14,icalc14,idoc14,iphy14,izoo14, &
                              ipowc14,issso14,isssc14
    use mo_sedmnt,      only: powtra,sedlay,seddw,porwat,porsol

    ! Arguments
    integer,intent(in) :: kpie,kpje,kpke
    real,   intent(in) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
    real,   intent(in) :: pddpo(kpie,kpje,kpke),omask(kpie,kpje)

    ! Local variables
    real, parameter :: pg2ppm = 1.0_rp/2.13_rp  ! conversion factor PgC -> ppm CO2
    integer :: i,j,k
    real    :: ztmp1(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)
    real    :: co2flux, co2flux_ppm
    real    :: ztmp2(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy) ! cisonew
    real    :: co213flux, co213flux_ppm ! cisonew
    real    :: co214flux, co214flux_ppm ! cisonew
    real    :: totc14dec, vol ! cisonew

    co2flux      = 0.0_rp

    ! Calculate global total air-sea flux [kmol]
    ztmp1(:,:)   = 0.0_rp
    do j=1,kpje
      do i=1,kpie
        ztmp1(i,j) = atmflx(i,j,iatmco2)*pdlxp(i,j)*pdlyp(i,j) ![kmol CO2/ m2] * [m] * [m]
      enddo
    enddo

    call xcsum(co2flux,ztmp1,ips)

    ! Convert global CO2 flux to ppm
    co2flux_ppm  = co2flux*12._rp*1.e-12_rp*pg2ppm ! [kmol C] -> [ppm]

    ! Update atmospheric pCO2
    do  j=1,kpje
      do  i=1,kpie
        atm(i,j,iatmco2)=atm(i,j,iatmco2) + co2flux_ppm
      enddo
    enddo

    if (use_cisonew) then
      co213flux    = 0.0_rp
      co214flux    = 0.0_rp

      ! Calculate global total air-sea flux for C isotopes [kmol]
      ztmp1(:,:)   = 0.0_rp
      ztmp2(:,:)   = 0.0_rp
      do j=1,kpje
        do i=1,kpie
          ztmp1(i,j) = atmflx(i,j,iatmc13)*pdlxp(i,j)*pdlyp(i,j) ![kmol 13CO2/ m2] * [m] * [m]
          ztmp2(i,j) = atmflx(i,j,iatmc14)*pdlxp(i,j)*pdlyp(i,j) ![kmol 14CO2/ m2] * [m] * [m]
        enddo
      enddo

      call xcsum(co213flux,ztmp1,ips)
      call xcsum(co214flux,ztmp2,ips)

      ! Convert global CO2 isotope fluxes to ppm isotope fluxes
      co213flux_ppm  = co213flux*13._rp*1.e-12_rp*pg2ppm*12._rp/13._rp ! [kmol 13CO2] -> [ppm]
      co214flux_ppm  = co214flux*14._rp*1.e-12_rp*pg2ppm*12._rp/14._rp ! [kmol 14CO2] -> [ppm]

      ! Calculate sum of 14C decay. Only decay in ocean, so only ocean tracers.
      totc14dec    = 0.0_rp
      ztmp1(:,:)   = 0.0_rp
      do k=1,kpke
        do j=1,kpje
          do i=1,kpie
            vol        = pdlxp(i,j)*pdlyp(i,j)*pddpo(i,j,k)*omask(i,j) ! ocean volume
            ztmp1(i,j) = ztmp1(i,j)+ocetra(i,j,k,isco214)*vol*(1.0_rp-c14dec)
            ztmp1(i,j) = ztmp1(i,j)+ocetra(i,j,k,idet14) *vol*(1.0_rp-c14dec)*rcar
            ztmp1(i,j) = ztmp1(i,j)+ocetra(i,j,k,icalc14)*vol*(1.0_rp-c14dec)
            ztmp1(i,j) = ztmp1(i,j)+ocetra(i,j,k,idoc14) *vol*(1.0_rp-c14dec)*rcar
            ztmp1(i,j) = ztmp1(i,j)+ocetra(i,j,k,iphy14) *vol*(1.0_rp-c14dec)*rcar
            ztmp1(i,j) = ztmp1(i,j)+ocetra(i,j,k,izoo14) *vol*(1.0_rp-c14dec)*rcar
            if (.not. use_sedbypass) then
              vol        = seddw(k)*pdlxp(i,j)*pdlyp(i,j)*porwat(i,j,k)*omask(i,j) ! porewater volume
              ztmp1(i,j) = ztmp1(i,j)+powtra(i,j,k,ipowc14) *vol*(1.0_rp-c14dec)
              vol        = seddw(k)*pdlxp(i,j)*pdlyp(i,j)*porsol(i,j,k)*omask(i,j) ! sediment volume
              ztmp1(i,j) = ztmp1(i,j)+sedlay(i,j,k,issso14) *vol*(1.0_rp-c14dec)*rcar
              ztmp1(i,j) = ztmp1(i,j)+sedlay(i,j,k,isssc14) *vol*(1.0_rp-c14dec)
            endif
          enddo
        enddo
      enddo

      call xcsum(totc14dec,ztmp1,ips)

      ! Update atmospheric p13CO2 and p14CO2
      do  j=1,kpje
        do  i=1,kpie
          atm(i,j,iatmc13)=atm(i,j,iatmc13) + co213flux_ppm
          atm(i,j,iatmc14)=atm(i,j,iatmc14) + co214flux_ppm
          atm(i,j,iatmc14)=atm(i,j,iatmc14) + totc14dec*14._rp*1.e-12_rp*pg2ppm*12._rp/14._rp ! add 14C decay (ppm)
        enddo
      enddo

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ' '
        write(io_stdo_bgc,*) 'Boxatm fluxes (ppm)'
        write(io_stdo_bgc,*) ' co213flux_ppm: ',co213flux_ppm
        write(io_stdo_bgc,*) ' co214flux_ppm: ',co214flux_ppm
        write(io_stdo_bgc,*) ' totc14dec (ppm): ',(totc14dec*14._rp*1.e-12_rp*pg2ppm*12._rp/14._rp)
        write(io_stdo_bgc,*) ' '
      endif

    endif ! end of use_cisonew

  end subroutine update_boxatm
  !*************************************************************************************************

end module mo_boxatm
