! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2022 Mats Bentsen, Lars Inge Enstad

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

module mod_bigrid

  use dimensions, only: idm, jdm, itdm, jtdm
  use mod_xc,     only: xcmax, xcstop, xctilr, xchalt, xcsync, &  ! HYCOM communication interface
                        i0, j0, ii, jj, iq, ip, iu, iv, &
                        nbdy, ifq, ilq, isq, jfq, jlq, jsq, &
                        ifp, ilp, isp, jfp, jlp, jsp, ifu, ilu, &
                        isu, jfu, jlu, jsu, ifv, ilv, jfv, jsv, isv, jlv, &
                        nreg, mnproc, lp, ms, no_flush, &
                        halo_ps, halo_us, halo_vs, halo_qs

  implicit none
  private

  ! Public routines
  public  :: bigrid

  ! Private routines
  private :: indxi
  private :: indxj

contains

  subroutine bigrid(depth)

    ! --- set loop bounds for irregular basin in c-grid configuration
    ! --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
    ! --- 'depth' = basin depth array, zero values indicate land

    ! Arguments
    real, intent(inout), dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: depth

    ! Local variables
    real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: util1,util2,util3
    logical :: lperiodi,lperiodj,larctic
    integer :: i,j,nfill,nzero
    real :: depmax

    ! --- is the domain periodic in i-index?
    depmax = 0.0
    if (i0+ii == itdm) then
      do j= 1,jj
        depmax = max(depmax,depth(ii,j))
      end do
    end if
    call xcmax(depmax)
    lperiodi = depmax > 0.0

    ! --- is the domain periodic in j-index?
    depmax = 0.0
    if     (j0+jj == jtdm) then
      do i= 1,ii
        depmax = max(depmax,depth(i,jj))
      end do
    end if
    call xcmax(depmax)
    larctic = depmax > 0.0 .and. nreg == 2
    lperiodj = depmax > 0.0 .and. nreg /= 2

    ! --- is this consistent with nreg (from mod_xc)?
    if     (.not.lperiodi.and..not.lperiodj.and. &
         (nreg == 0.or.nreg == -1)) then
      nreg = 0 ! closed domain
    else if (lperiodi.and..not.lperiodj.and. &
         (nreg == 1.or.nreg == -1)) then
      nreg = 1 ! periodic domain in i-index
    else if (lperiodi.and.larctic.and. &
         (nreg == 2.or.nreg == -1)) then
      nreg = 2 ! global domain with arctic patch
    else if (lperiodi.and.lperiodj.and. &
         (nreg == 3.or.nreg == -1)) then
      nreg = 3 ! doubly periodic domain
    else if (.not.lperiodi.and.lperiodj.and. &
         (nreg == 4.or.nreg == -1)) then
      nreg = 4 ! periodic domain in j-index
    else
      if (mnproc == 1) then
        write(lp,'(/a,i2)') 'bigrid: nreg     =',nreg
        write(lp,'(a,l1)')  'bigrid: lperiodi =',lperiodi
        write(lp,'(a,l1)')  'bigrid: larctic  =',larctic
        write(lp,'(a,l1)')  'bigrid: lperiodj =',lperiodj
        write(lp,'(a/)')    'basin depth array inconsistent with nreg'
        call flush(lp)
      end if
      call xcstop('(bigrid)')
      stop '(bigrid)'
    end if

    if     (mnproc == 1) then
      write(lp,'(/a,i2)') 'bigrid: nreg =',nreg
      if     (nreg == 0) then
        write(lp,'(a/)') 'bigrid: closed domain'
      else if (nreg == 1) then
        write(lp,'(a/)') 'bigrid: periodic domain in i-index'
      else if (nreg == 2) then
        write(lp,'(a/)') 'bigrid: global domain with arctic patch'
      else if (nreg == 3) then
        write(lp,'(a/)') 'bigrid: doubly periodic domain'
      else if (nreg == 4) then
        write(lp,'(a/)') 'bigrid: periodic domain in j-index'
      end if
      call flush(lp)
    end if

    ! --- nreg is defined, so now safe to update halo
    call xctilr(depth, 1, 1, nbdy,nbdy, halo_ps)

    ! --- allow for non-periodic and non-arctic boundaries (part I).
    if     (.not.lperiodj .and. j0 == 0) then
      ! ---   south boundary is all land.
      do j = 1-nbdy,0
        do i = 1-nbdy,ii+nbdy
          depth(i,j) = 0.0
        end do
      end do
    end if

    if     (.not.lperiodj .and. .not.larctic .and. j0+jj == jtdm) then
      ! ---   north boundary is all land.
      do j = jj+1,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          depth(i,j) = 0.0
        end do
      end do
    end if

    if     (.not.lperiodi .and. i0 == 0) then
      ! ---   west boundary is all land.
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,0
          depth(i,j) = 0.0
        end do
      end do
    end if

    if     (.not.lperiodi .and. i0+ii == itdm) then
      ! ---   east boundary is all land.
      do j = 1-nbdy,jj+nbdy
        do i = ii+1,ii+nbdy
          depth(i,j) = 0.0
        end do
      end do
    end if

    ! --- detect (and abort on) single-width inlets and 1-point seas.
    nfill = 0
    do j = 1,jj
      do i = 1,ii
        nzero = 0
        if (depth(i,j) > 0.0) then
          if (depth(i-1,j) <= 0.0) nzero = nzero+1
          if (depth(i+1,j) <= 0.0) nzero = nzero+1
          if (depth(i,j-1) <= 0.0) nzero = nzero+1
          if (depth(i,j+1) <= 0.0) nzero = nzero+1
          if (nzero >= 3) then
            write (lp,'(a,i4,a,i4,a,i1,a)') &
                 'error - dh(',i0+i,',',j0+j,') has ', &
                 nzero,' land nieghbours'
            nfill = nfill+1
          end if
        end if
      end do
    end do
    call xcmax(nfill)
    if (nfill > 0) then
      if (mnproc == 1) then
        write(lp,'(/a/)') &
             'Must correct bathymetry before running BLOM'
        call flush(lp)
      end if
      call xcstop('(bigrid)')
      stop '(bigrid)'
    end if

    ! --- start out with masks as land everywhere
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1-nbdy,jdm+nbdy
      do i = 1-nbdy,idm+nbdy
        ip(i,j) = 0
        iq(i,j) = 0
        iu(i,j) = 0
        iv(i,j) = 0
        util1(i,j) = 0.
        util2(i,j) = 0.
        util3(i,j) = 0.
      end do
    end do
    !$OMP END PARALLEL DO

    ! --- mass points are defined where water depth is greater than zero
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1-nbdy,jj+nbdy
      do i = 1-nbdy,ii+nbdy
        if (depth(i,j) > 0.) then
          ip(i,j) = 1
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    ! --- u,v points are located halfway between any 2 adjoining mass points
    ! --- 'interior' q points require water on all 4 sides.
    ! --- 'promontory' q points require water on 3 (or at least 2
    ! --- diametrically opposed) sides
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1,jj
      do i = 1,ii
        if (ip(i-1,j) > 0.and.ip(i,j) > 0) then
          iu(i,j) = 1
        end if
        if (ip(i,j-1) > 0.and.ip(i,j) > 0) then
          iv(i,j) = 1
        end if
        if (min(ip(i,j),ip(i-1,j),ip(i,j-1),ip(i-1,j-1)) > 0) then
          iq(i,j) = 1
        else if ((ip(i  ,j) > 0.and.ip(i-1,j-1) > 0).or. &
             (ip(i-1,j) > 0.and.ip(i  ,j-1) > 0)    ) then
          iq(i,j) = 1
        end if
        util1(i,j) = iu(i,j)
        util2(i,j) = iv(i,j)
        util3(i,j) = iq(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    call xctilr(util1,1,1, nbdy,nbdy, halo_us)
    call xctilr(util2,1,1, nbdy,nbdy, halo_vs)
    call xctilr(util3,1,1, nbdy,nbdy, halo_qs)
    !$OMP PARALLEL DO PRIVATE(i)
    do j= 1-nbdy,jj+nbdy
      do i= 1-nbdy,ii+nbdy
        iu(i,j) = nint(util1(i,j))
        iv(i,j) = nint(util2(i,j))
        iq(i,j) = nint(util3(i,j))
      end do
    end do
    !$OMP END PARALLEL DO

    ! --- allow for non-periodic and non-arctic boundaries (part II).
    if     (.not.lperiodj .and. j0 == 0) then
      ! ---   south boundary is all land.
      do j = 1-nbdy,0
        do i = 1-nbdy,ii+nbdy
          iq(i,j) = 0
          iu(i,j) = 0
          iv(i,j) = 0
        end do
      end do
    end if

    if     (.not.lperiodj .and. .not.larctic .and. j0+jj == jtdm) then
      ! ---   north boundary is all land.
      do j = jj+1,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          iq(i,j) = 0
          iu(i,j) = 0
          iv(i,j) = 0
        end do
      end do
    end if

    if     (.not.lperiodi .and. i0 == 0) then
      ! ---   west boundary is all land.
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,0
          iq(i,j) = 0
          iu(i,j) = 0
          iv(i,j) = 0
        end do
      end do
    end if

    if     (.not.lperiodi .and. i0+ii == itdm) then
      ! ---   east boundary is all land.
      do j = 1-nbdy,jj+nbdy
        do i = ii+1,ii+nbdy
          iq(i,j) = 0
          iu(i,j) = 0
          iv(i,j) = 0
        end do
      end do
    end if

    ! --- determine loop bounds for vorticity points, including interior and
    ! --- promontory points
    call indxi(iq,ifq,ilq,isq)
    call indxj(iq,jfq,jlq,jsq)

    ! --- determine loop indices for mass and velocity points
    call indxi(ip,ifp,ilp,isp)
    call indxj(ip,jfp,jlp,jsp)
    call indxi(iu,ifu,ilu,isu)
    call indxj(iu,jfu,jlu,jsu)
    call indxi(iv,ifv,ilv,isv)
    call indxj(iv,jfv,jlv,jsv)

  end subroutine bigrid

  ! ------------------------------------------------------------------------------
  subroutine indxi(ipt,if,il,is)

    ! --- input array ipt contains 1 at grid point locations, 0 elsewhere
    ! --- output is arrays if, il, is  where
    ! --- if(j,k) gives row index of first point in column j for k-th section
    ! --- il(j,k) gives row index of last point
    ! --- is(j) gives number of sections in column j (maximum: ms)

    ! Arguments
    integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ipt
    integer, dimension (1-nbdy:jdm+nbdy,ms) :: if,il
    integer, dimension (1-nbdy:jdm+nbdy) :: is

    ! Local variables
    integer :: i,j,k,last

    do j = 1-nbdy,jj+nbdy
      is(j) = 0
      do k = 1,ms
        if(j,k) = 0
        il(j,k) = 0
      end do

      k = 1
      last = ipt(1-nbdy,j)
      if     (last  ==  1) then
        if(j,k) = 1-nbdy
      end if
      do i = 2-nbdy,ii+nbdy
        if      (last  ==  1 .and. ipt(i,j)  ==  0) then
          il(j,k) = i-1
          k = k+1
        else if (last  ==  0 .and. ipt(i,j)  ==  1) then
          if     (k  >  ms) then
            write(lp,'(a,i5)')  'indxi problem on proc ',mnproc
            write(lp,'(a,2i5)') &
                 ' error in indxi -- ms too small at i,j =',i0+i,j0+j
            call xchalt('(indxi)')
            stop '(indxi)'
          end if
          if(j,k) = i
        end if
        last = ipt(i,j)
      end do
      if     (last  ==  1) then
        il(j,k) = ii+nbdy
        is(j) = k
      else
        is(j) = k-1
      end if
    end do
    call xcsync(no_flush)

  end subroutine indxi

  ! ------------------------------------------------------------------------------
  subroutine indxj(jpt,jf,jl,js)

    ! --- input array jpt contains 1 at grid point locations, 0 elsewhere
    ! --- output is arrays jf, jl, js  where
    ! --- jf(i,k) gives column index of first point in row i for k-th section
    ! --- jl(i,k) gives column index of last point
    ! --- js(i) gives number of sections in row i (maximum: ms)

    ! Arguments
    integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: jpt
    integer, dimension (1-nbdy:idm+nbdy,ms) :: jf,jl
    integer, dimension (1-nbdy:idm+nbdy) :: js

    ! Local variables
    integer :: i,j,k,last

    do i = 1-nbdy,ii+nbdy
      js(i) = 0
      do k = 1,ms
        jf(i,k) = 0
        jl(i,k) = 0
      end do

      k = 1
      last = jpt(i,1-nbdy)
      if     (last  ==  1) then
        jf(i,k) = 1-nbdy
      end if
      do j = 2-nbdy,jj+nbdy
        if      (last  ==  1 .and. jpt(i,j)  ==  0) then
          jl(i,k) = j-1
          k = k+1
        else if (last  ==  0 .and. jpt(i,j)  ==  1) then
          if     (k  >  ms) then
            write(lp,'(a,i5)')  'indxj problem on proc ',mnproc
            write(lp,'(a,2i5)') &
                 ' error in indxj -- ms too small at i,j =',i0+i,j0+j
            call xchalt('(indxj)')
            stop '(indxj)'
          end if
          jf(i,k) = j
        end if
        last = jpt(i,j)
      end do
      if     (last  ==  1) then
        jl(i,k) = jj+nbdy
        js(i) = k
      else
        js(i) = k-1
      end if
    end do
    call xcsync(no_flush)

  end subroutine indxj

end module mod_bigrid
