! ------------------------------------------------------------------------------
! Copyright (C) 2015-2024 Mats Bentsen, Alok Kumar Gupta, Mehmet Ilicak,
!                         Aleksi Nummelin, Mariana Vertenstein

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_inigeo

  use dimensions,           only: nreg
  use mod_config,           only: expcnf
  use mod_fuk95,            only: geoenv_fuk95
  use mod_channel,          only: geoenv_channel
  use mod_single_column,    only: geoenv_single_column
  use mod_geoenv,           only: geoenv_file, geoenv_test
  use mod_geoenv_cesmextra, only: geoenv_cesmextra
  use mod_xc,               only: xcstop, xctilr, xcsum, xchalt, &
                                  ii, jj, ips
  use mod_xc,               only: xcstop, xctilr, xcsum, xchalt, &
                                  ip, iq, iu, iv, ipwocn, halo_vs, cplmsk, &
                                  nproc, jpr, halo_ps, halo_us, halo_qs
  use mod_grid,             only: qclon, qclat, pclon, pclat, uclon, uclat, &
                                  vclon, vclat, scqx, scqy, scpx, scpy, scux, &
                                  scuy, scvx, scvy, scq2, scp2, scu2, scv2, &
                                  scq2i, scp2i, scq2i, scuxi, scuyi, scvxi, &
                                  scvyi, qlon, qlat, plon, plat, ulon, ulat, &
                                  vlon, vlat, depths, corioq, coriop, betafp, &
                                  area, nwp, inivar_grid
  use mod_utility,          only: util1, util2
  use mod_dia,              only: iotype
  use mod_checksum,         only: csdiag, chksummsk
  use mod_bigrid,           only: bigrid
  use mod_fill_global,      only: fill_global
  use mod_nctools

  implicit none
  private

  public :: inigeo

contains

  subroutine inigeo

    ! ------------------------------------------------------------------
    ! Initialize the geographic environment
    ! ------------------------------------------------------------------

    ! Local variables
    real :: avgbot
    real, parameter :: mval = -1.e12
    real, parameter :: fval = -1.e13
    real, dimension(itdm,jtdm) :: tmpg
    real :: rnwp,rmxnbp,rtnbp,rnbp
    integer :: i,j,k,l,kmax

    ! ------------------------------------------------------------------
    ! Initialize grid variables.
    ! ------------------------------------------------------------------

    call inivar_grid

    ! ------------------------------------------------------------------
    ! Define bathymetry, grid specification and Coriolis parameter
    ! ------------------------------------------------------------------

    select case (trim(expcnf))
      case ('cesm')
        call geoenv_file
        call geoenv_cesmextra
      case ('ben02clim', 'ben02syn', 'isomip1', 'isomip2')
        call geoenv_file
      case ('fuk95')
        call geoenv_fuk95
      case ('channel')
        call geoenv_channel
      case ('single_column')
        call geoenv_single_column
      case ('test')
        call geoenv_test
      case default
        if (mnproc == 1) then
          write (lp,'(3a)') ' inigeo: expcnf = ', trim(expcnf),' is unsupported!'
        end if
        call xcstop('(inigeo)')
        stop '(inigeo)'
    end select

    ! ------------------------------------------------------------------
    ! Compute auxilary grid parameters
    ! ------------------------------------------------------------------

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        scq2i(i,j) = 1./max(1.,scq2(i,j))
        scp2i(i,j) = 1./max(1.,scp2(i,j))
        scuxi(i,j) = 1./max(1.,scux(i,j))
        scvyi(i,j) = 1./max(1.,scvy(i,j))
        scuyi(i,j) = 1./max(1.,scuy(i,j))
        scvxi(i,j) = 1./max(1.,scvx(i,j))
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------------
    ! Determine do-loop limits for u,v,p,q points
    ! ------------------------------------------------------------------

    call bigrid(depths)

    ! ------------------------------------------------------------------
    ! Update halos for parameters related to the geographic environment
    ! ------------------------------------------------------------------

    call xctilr(qlat, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(qlon, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(plat, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(plon, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(ulat, 1,1, nbdy,nbdy, halo_us)
    call xctilr(ulon, 1,1, nbdy,nbdy, halo_us)
    call xctilr(vlat, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(vlon, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(scqx, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(scqy, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(scpx, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(scpy, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(scux, 1,1, nbdy,nbdy, halo_us)
    call xctilr(scuy, 1,1, nbdy,nbdy, halo_us)
    call xctilr(scvx, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(scvy, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(scq2, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(scp2, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(scu2, 1,1, nbdy,nbdy, halo_us)
    call xctilr(scv2, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(scq2i, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(scp2i, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(scuxi, 1,1, nbdy,nbdy, halo_us)
    call xctilr(scvyi, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(scuyi, 1,1, nbdy,nbdy, halo_us)
    call xctilr(scvxi, 1,1, nbdy,nbdy, halo_vs)
    call xctilr(corioq, 1,1, nbdy,nbdy, halo_qs)
    call xctilr(coriop, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(betafp, 1,1, nbdy,nbdy, halo_ps)
    call xctilr(qclat, 1,4, nbdy,nbdy, halo_qs)
    call xctilr(qclon, 1,4, nbdy,nbdy, halo_qs)
    call xctilr(pclat, 1,4, nbdy,nbdy, halo_ps)
    call xctilr(pclon, 1,4, nbdy,nbdy, halo_ps)
    call xctilr(uclat, 1,4, nbdy,nbdy, halo_us)
    call xctilr(uclon, 1,4, nbdy,nbdy, halo_us)
    call xctilr(vclat, 1,4, nbdy,nbdy, halo_vs)
    call xctilr(vclon, 1,4, nbdy,nbdy, halo_vs)
    if (expcnf == 'cesm') then
      !$omp parallel do private(i)
      do j = 1,jj
        do i = 1,ii
          util1(i,j) = cplmsk(i,j)
        end do
      end do
      !$omp end parallel do

      call xctilr(util1, 1,1, nbdy,nbdy, halo_ps)

      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          cplmsk(i,j) = nint(util1(i,j))
        end do
      end do
      !$omp end parallel do
    end if

    ! ------------------------------------------------------------------
    ! Set mask used for global sums
    ! ------------------------------------------------------------------

    if (nreg == 2.and.nproc == jpr) then
      !$omp parallel do private(i)
      do j = 1-nbdy,jj-1
        do i = 1-nbdy,ii+nbdy
          ips(i,j) = ip(i,j)
        end do
      end do
      !$omp end parallel do
      do j = jj,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          ips(i,j) = 0
        end do
      end do
    else
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          ips(i,j) = ip(i,j)
        end do
      end do
      !$omp end parallel do
    end if
    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        util1(i,j) = ip(i,j)
      end do
    end do
    !$omp end parallel do

    call xcsum(rnwp,util1,ip)

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        util1(i,j) = depths(i,j)*scp2(i,j)
        util2(i,j) = scp2(i,j)
      end do
    end do
    !$omp end parallel do

    call xcsum(avgbot,util1,ips)
    call xcsum(area,  util2,ips)
    avgbot = avgbot/area
    if (mnproc == 1) then
      if (nwp /= nint(rnwp)) then
        write (lp,'(a)') ' xcsum test failed!'
        write (lp,'(a,i7)') ' number of wet points:',nwp
        write (lp,'(a,i7)') ' xcsum on ocean mask: ',nint(rnwp)
        call xchalt('(inigeo)')
        stop '(inigeo)'
      end if
      write (lp,100) avgbot,area
100   format(' mean basin depth (m) and area (10^6 km^2):',f9.1,-16p,f9.1)
    end if

    ! ------------------------------------------------------------------
    ! Set mask for grid cells connected to the world ocean
    ! ------------------------------------------------------------------

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        util1(i,j) = ips(i,j)
      end do
    end do
    !$omp end parallel do
    call xcsum(rnwp,util1,ips)
    if (mnproc == 1) then
      write (lp,*) 'Number of wet points',nint(rnwp)
    end if

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        if (ips(i,j) == 1) then
          util1(i,j) = fval
        else
          util1(i,j) = mval
        end if
      end do
    end do
    !$omp end parallel do

    k = 0
    rmxnbp = 0.
    rtnbp = 0.
    do
      k = k+1
      call xcaget(tmpg,util1,1)
      if (mnproc == 1) then
        do l = 1,itdm*jtdm
          j = (l-1)/itdm+1
          i = l-(j-1)*itdm
          if (tmpg(i,j) == fval) then
            tmpg(i,j) = k
            exit
          end if
        end do
      end if
      call xcaput(tmpg,util1,1)
      call fill_global(mval,fval,halo_ps,util1)
      !$omp parallel do private(i)

      do j = 1,jj
        do i = 1,ii
          if (util1(i,j) == mval.or.util1(i,j) == fval) then
            util2(i,j) = 0.
          else
            if (nint(util1(i,j)) == k) then
              util2(i,j) = 1.
            else
              util2(i,j) = 0.
            end if
          end if
        end do
      end do
      !$omp end parallel do

      call xcsum(rnbp,util2,ips)
      if (mnproc == 1) then
        write (lp,*) 'Number of basin points',nint(rnbp)
      end if
      if (rnbp > rmxnbp) then
        rmxnbp = rnbp
        kmax = k
      end if
      rtnbp = rtnbp+rnbp
      if (nint(rtnbp-rnwp) == 0) exit
    end do

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        if (util1(i,j) == mval) then
          util1(i,j) = 0.
        else
          if (nint(util1(i,j)) == kmax) then
            util1(i,j) = 1.
          else
            util1(i,j) = 0.
          end if
        end if
      end do
    end do
    !$omp end parallel do

    call xctilr(util1, 1,1, nbdy,nbdy, halo_ps)

    !$omp parallel do private(i)
    do j = 1-nbdy,jj+nbdy
      do i = 1-nbdy,ii+nbdy
        ipwocn(i,j) = nint(util1(i,j))
      end do
    end do
    !$omp end parallel do

    call ncfopn('ipwocn.nc','w','c',1,iotype)
    call ncdims('x',itdm)
    call ncdims('y',jtdm)
    call ncdefvar('ipwocn','x y',nfint,2)
    call ncedef
    call ncwrti('ipwocn','x y',ipwocn,ip,1)
    call ncfcls

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'inigeo:'
      end if
      call chksummsk(depths,ip,1,'depths')
      call chksummsk(plat,ip,1,'plat')
      call chksummsk(plon,ip,1,'plon')
      call chksummsk(pclat,ip,4,'pclat')
      call chksummsk(pclon,ip,4,'pclon')
      call chksummsk(corioq,iq,1,'corioq')
      call chksummsk(coriop,ip,1,'coriop')
      call chksummsk(betafp,ip,1,'betafp')
      call chksummsk(scqx,iq,1,'scqx')
      call chksummsk(scqy,iq,1,'scqy')
      call chksummsk(scpx,ip,1,'scpx')
      call chksummsk(scpy,ip,1,'scpy')
      call chksummsk(scux,iu,1,'scux')
      call chksummsk(scuy,iu,1,'scuy')
      call chksummsk(scvx,iv,1,'scvx')
      call chksummsk(scvy,iv,1,'scvy')
      call chksummsk(scq2,iq,1,'scq2')
      call chksummsk(scp2,ip,1,'scp2')
      call chksummsk(scu2,iu,1,'scu2')
      call chksummsk(scv2,iv,1,'scv2')
      call chksummsk(scp2i,ip,1,'scp2i')
      call chksummsk(scq2i,iq,1,'scq2i')
      call chksummsk(scuxi,iu,1,'scuxi')
      call chksummsk(scvyi,iv,1,'scvyi')
      call chksummsk(scuyi,iu,1,'scuyi')
      call chksummsk(scvxi,iv,1,'scvxi')
    end if

  end subroutine inigeo

end module mod_inigeo
