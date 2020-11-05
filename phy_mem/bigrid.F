! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2020 Mats Bentsen, Lars Inge Enstad
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

      subroutine bigrid(depth)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &        depth
c
c --- set loop bounds for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &        util1,util2,util3
c
      integer    nchar
      parameter (nchar=120)
      logical   lperiodi,lperiodj,larctic
c
      integer   i,j,nfill,nzero
      real      depmax
c
      character fmt*13
      data fmt/'(i4,1x,120i1)'/
c
c --- is the domain periodic in i-index?
      depmax=0.0
      if     (i0+ii.eq.itdm) then
        do j= 1,jj
          depmax=max(depmax,depth(ii,j))
        enddo
      endif
      call xcmax(depmax)
      lperiodi=depmax.gt.0.0
c
c --- is the domain periodic in j-index?
      depmax=0.0
      if     (j0+jj.eq.jtdm) then
        do i= 1,ii
          depmax=max(depmax,depth(i,jj))
        enddo
      endif
      call xcmax(depmax)
      larctic=depmax.gt.0.0 .and. nreg.eq.2
      lperiodj=depmax.gt.0.0 .and. nreg.ne.2
c
c --- is this consistent with nreg (from mod_xc)?
      if     (.not.lperiodi.and..not.lperiodj.and.
     &        (nreg.eq.0.or.nreg.eq.-1)) then
        nreg=0 ! closed domain
      elseif (lperiodi.and..not.lperiodj.and.
     &        (nreg.eq.1.or.nreg.eq.-1)) then
        nreg=1 ! periodic domain in i-index
      elseif (lperiodi.and.larctic.and.
     &        (nreg.eq.2.or.nreg.eq.-1)) then
        nreg=2 ! global domain with arctic patch
      elseif (lperiodi.and.lperiodj.and.
     &        (nreg.eq.3.or.nreg.eq.-1)) then
        nreg=3 ! doubly periodic domain
      elseif (.not.lperiodi.and.lperiodj.and.
     &        (nreg.eq.4.or.nreg.eq.-1)) then
        nreg=4 ! periodic domain in j-index
      else  
        if (mnproc.eq.1) then
          write(lp,'(/a,i2)') 'bigrid: nreg     =',nreg
          write(lp,'(a,l1)')  'bigrid: lperiodi =',lperiodi
          write(lp,'(a,l1)')  'bigrid: larctic  =',larctic
          write(lp,'(a,l1)')  'bigrid: lperiodj =',lperiodj
          write(lp,'(a/)')    'basin depth array inconsistent with nreg'
          call flush(lp)
        endif
        call xcstop('(bigrid)')
               stop '(bigrid)'
      endif
c
      if     (mnproc.eq.1) then
        write(lp,'(/a,i2)') 'bigrid: nreg =',nreg
        if     (nreg.eq.0) then
          write(lp,'(a/)') 'bigrid: closed domain'
        elseif (nreg.eq.1) then
          write(lp,'(a/)') 'bigrid: periodic domain in i-index'
        elseif (nreg.eq.2) then
          write(lp,'(a/)') 'bigrid: global domain with arctic patch'
        elseif (nreg.eq.3) then
          write(lp,'(a/)') 'bigrid: doubly periodic domain'
        elseif (nreg.eq.4) then
          write(lp,'(a/)') 'bigrid: periodic domain in j-index'
        endif
        call flush(lp)
      endif
c
c --- nreg is defined, so now safe to update halo
      call xctilr(depth,1,1, nbdy,nbdy, halo_ps)
c
c --- allow for non-periodic and non-arctic boundaries (part I).
      if     (.not.lperiodj .and. j0.eq.0) then
c ---   south boundary is all land.
        do j=1-nbdy,0
          do i=1-nbdy,ii+nbdy
            depth(i,j) = 0.0
          enddo
        enddo
      endif
c
      if     (.not.lperiodj .and. .not.larctic .and. j0+jj.eq.jtdm) then
c ---   north boundary is all land.
        do j=jj+1,jj+nbdy
          do i=1-nbdy,ii+nbdy
            depth(i,j) = 0.0
          enddo
        enddo
      endif
c
      if     (.not.lperiodi .and. i0.eq.0) then
c ---   west boundary is all land.
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,0
            depth(i,j) = 0.0
          enddo
        enddo
      endif
c
      if     (.not.lperiodi .and. i0+ii.eq.itdm) then
c ---   east boundary is all land.
        do j=1-nbdy,jj+nbdy
          do i=ii+1,ii+nbdy
            depth(i,j) = 0.0
          enddo
        enddo
      endif
c
c --- detect (and abort on) single-width inlets and 1-point seas.
      nfill=0
      do j=1,jj
        do i=1,ii
          nzero=0
          if (depth(i,j).gt.0.0) then
            if (depth(i-1,j).le.0.0) nzero=nzero+1
            if (depth(i+1,j).le.0.0) nzero=nzero+1
            if (depth(i,j-1).le.0.0) nzero=nzero+1
            if (depth(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (lp,'(a,i4,a,i4,a,i1,a)')
     &          'error - dh(',i0+i,',',j0+j,') has ',
     &          nzero,' land nieghbours'
              nfill=nfill+1
            endif
          endif
        enddo
      enddo
      call xcmax(nfill)
      if (nfill.gt.0) then
        if (mnproc.eq.1) then
          write(lp,'(/a/)') 
     &     'Must correct bathymetry before running BLOM'
          call flush(lp)
        endif
        call xcstop('(bigrid)')
               stop '(bigrid)'
      endif
c
c --- start out with masks as land everywhere
c$OMP PARALLEL DO PRIVATE(i)
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          ip(i,j)=0
          iq(i,j)=0
          iu(i,j)=0
          iv(i,j)=0
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- mass points are defined where water depth is greater than zero
c$OMP PARALLEL DO PRIVATE(i)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          if (depth(i,j).gt.0.) then
            ip(i,j)=1
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- u,v points are located halfway between any 2 adjoining mass points
c --- 'interior' q points require water on all 4 sides.
c --- 'promontory' q points require water on 3 (or at least 2
c --- diametrically opposed) sides
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
            iu(i,j)=1
          endif
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
            iv(i,j)=1
          endif
          if (min(ip(i,j),ip(i-1,j),ip(i,j-1),ip(i-1,j-1)).gt.0) then
            iq(i,j)=1
          elseif ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or.
     &            (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) then
            iq(i,j)=1
          endif
          util1(i,j)=iu(i,j)
          util2(i,j)=iv(i,j)
          util3(i,j)=iq(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1,1,1, nbdy,nbdy, halo_us)
      call xctilr(util2,1,1, nbdy,nbdy, halo_vs)
      call xctilr(util3,1,1, nbdy,nbdy, halo_qs)
c$OMP PARALLEL DO PRIVATE(i)
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
          iu(i,j)=nint(util1(i,j))
          iv(i,j)=nint(util2(i,j))
          iq(i,j)=nint(util3(i,j))
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- allow for non-periodic and non-arctic boundaries (part II).
      if     (.not.lperiodj .and. j0.eq.0) then
c ---   south boundary is all land.
        do j=1-nbdy,0
          do i=1-nbdy,ii+nbdy
            iq(i,j) = 0
            iu(i,j) = 0
            iv(i,j) = 0
          enddo
        enddo
      endif
c
      if     (.not.lperiodj .and. .not.larctic .and. j0+jj.eq.jtdm) then
c ---   north boundary is all land.
        do j=jj+1,jj+nbdy
          do i=1-nbdy,ii+nbdy
            iq(i,j) = 0
            iu(i,j) = 0
            iv(i,j) = 0
          enddo
        enddo
      endif
c
      if     (.not.lperiodi .and. i0.eq.0) then
c ---   west boundary is all land.
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,0
            iq(i,j) = 0
            iu(i,j) = 0
            iv(i,j) = 0
          enddo
        enddo
      endif
c
      if     (.not.lperiodi .and. i0+ii.eq.itdm) then
c ---   east boundary is all land.
        do j=1-nbdy,jj+nbdy
          do i=ii+1,ii+nbdy
            iq(i,j) = 0
            iu(i,j) = 0
            iv(i,j) = 0
          enddo
        enddo
      endif
c
c --- determine loop bounds for vorticity points, including interior and
c --- promontory points
      call indxi(iq,ifq,ilq,isq)
      call indxj(iq,jfq,jlq,jsq)
c
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
      call indxj(ip,jfp,jlp,jsp)
      call indxi(iu,ifu,ilu,isu)
      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
      call indxj(iv,jfv,jlv,jsv)
c
      return
      end
c
c
      subroutine indxi(ipt,if,il,is)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &         ipt
      integer, dimension (1-nbdy:jdm+nbdy,ms) ::
     &         if,il
      integer, dimension (1-nbdy:jdm+nbdy) ::
     &         is
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      integer i,j,k,last
c
      do j=1-nbdy,jj+nbdy
        is(j) = 0
        do k=1,ms
          if(j,k) = 0
          il(j,k) = 0
        enddo
c
        k=1
        last = ipt(1-nbdy,j)
        if     (last .eq. 1) then
          if(j,k) = 1-nbdy
        endif
        do i=2-nbdy,ii+nbdy
          if      (last .eq. 1 .and. ipt(i,j) .eq. 0) then
            il(j,k) = i-1
            k = k+1
          elseif (last .eq. 0 .and. ipt(i,j) .eq. 1) then
            if     (k .gt. ms) then
              write(lp,'(a,i5)')  'indxi problem on proc ',mnproc
              write(lp,'(a,2i5)') 
     &          ' error in indxi -- ms too small at i,j =',i0+i,j0+j
              call xchalt('(indxi)')
                     stop '(indxi)'
            endif
            if(j,k) = i
          endif
          last = ipt(i,j)
        enddo
        if     (last .eq. 1) then
          il(j,k) = ii+nbdy
          is(j) = k
        else
          is(j) = k-1
        endif
      enddo
      call xcsync(no_flush)
      return
      end
c
      subroutine indxj(jpt,jf,jl,js)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &         jpt
      integer, dimension (1-nbdy:idm+nbdy,ms) ::
     &         jf,jl
      integer, dimension (1-nbdy:idm+nbdy) ::
     &         js
c
c --- input array jpt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays jf, jl, js  where
c --- jf(i,k) gives column index of first point in row i for k-th section
c --- jl(i,k) gives column index of last point
c --- js(i) gives number of sections in row i (maximum: ms)
c
      integer i,j,k,last
c
      do i=1-nbdy,ii+nbdy
        js(i) = 0
        do k=1,ms
          jf(i,k) = 0
          jl(i,k) = 0
        enddo
c
        k=1
        last = jpt(i,1-nbdy)
        if     (last .eq. 1) then
          jf(i,k) = 1-nbdy
        endif
        do j=2-nbdy,jj+nbdy
          if      (last .eq. 1 .and. jpt(i,j) .eq. 0) then
            jl(i,k) = j-1
            k = k+1
          elseif (last .eq. 0 .and. jpt(i,j) .eq. 1) then
            if     (k .gt. ms) then
              write(lp,'(a,i5)')  'indxj problem on proc ',mnproc
              write(lp,'(a,2i5)')
     &          ' error in indxj -- ms too small at i,j =',i0+i,j0+j
              call xchalt('(indxj)')
                     stop '(indxj)'
            endif
            jf(i,k) = j
          endif
          last = jpt(i,j)
        enddo
        if     (last .eq. 1) then
          jl(i,k) = jj+nbdy
          js(i) = k
        else
          js(i) = k-1
        endif
      enddo
      call xcsync(no_flush)
      return
      end
