! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen
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

      module mod_inicon
c
c --- ------------------------------------------------------------------
c --- This module contains variables and procedures related to advection
c --- of layer pressure thickness and tracers by calling incremental
c --- remapping routines.
c --- ------------------------------------------------------------------
c
      use mod_types, only: r8
      use mod_config, only: expcnf
      use mod_constants, only: g, epsil, onem
      use mod_time, only: nstep, delt1, dlt
      use mod_xc
      use mod_grid, only: sigmar, scuy, scvx, scuyi, scvxi, depths,
     .                    corioq
      use mod_state, only: u, v, dp, dpu, dpv, temp, saln, sigma, p, pu,
     .                     pv, phi, ubflxs, vbflxs, ub, vb, pb, pbu,
     .                     pbv, ubflxs_p, vbflxs_p, pb_p, pbu_p, pbv_p,
     .                     ubcors_p, vbcors_p, kfpla
      use mod_pgforc, only: pgfx, pgfy, pgfxm, pgfym,
     .                      xixp, xixm, xiyp, xiym, pgforc
      use mod_barotp, only: ubflx, vbflx, pb_mn, ubflx_mn, vbflx_mn,
     .                      pvtrop
      use mod_tmsmt, only: dpold
      use mod_temmin, only: settemmin
      use mod_cesm, only: inicon_cesm
      use mod_fuk95, only: ictsz_fuk95
      use mod_channel, only: ictsz_channel
      use mod_eos, only: rho, sig, sofsig, delphi
      use mod_swtfrz, only: swtfrz
      use mod_pointtest, only: itest, jtest, ptest
      use mod_checksum, only: csdiag, chksummsk
      use netcdf
c
      implicit none
c
      private
c
c --- Variables to be set in namelist:
      character(len = 256) ::
     .  icfile ! Name of file containing initial conditions, that is
               ! either a valid restart file or 'inicon.nc' if
               ! climatological based initial conditions are desired.
c
      public :: icfile, inicon
c
      contains
c
c --- ------------------------------------------------------------------
c --- Private procedures.
c --- ------------------------------------------------------------------
c
      function getpl(th,s,phiu,phil,pup) result(plo)
c
c --- ------------------------------------------------------------------
c --- Get lower pressure interface of a layer knowing the temperature,
c --- salinity of the layer and the geopotential at upper and lower
c --- interface.
c --- ------------------------------------------------------------------
c
      real(r8), intent(in) ::
     .  th,           ! Layer potential temperature [deg C].
     .  s,            ! Layer salinity [g kg-1].
     .  phiu,         ! Geopotential at upper interface [cm2 s-2].
     .  phil,         ! Geopotential at lower interface [cm2 s-2].
     .  pup           ! Pressure at upper interface [g cm-1 s-2].
c
      real(r8) :: plo ! Pressure at lower interface [g cm-1 s-2].
c
      real(r8) :: q,dphi,alpu,alpl
c
c --- first guess on pressure interface
      plo=pup-rho(pup,th,s)*(phil-phiu)
c
c --- improve the accuracy of the pressure interface by an
c --- iterative procedure
      q=1._r8
      do while (abs(q).gt.1.e-4_r8)
        call delphi(pup,plo,th,s,dphi,alpu,alpl)
        q=(phil-phiu-dphi)/alpl
        plo=plo-q
      enddo
c
      end function getpl
c
c --- ------------------------------------------------------------------
c
      subroutine ictsz_file
c
c --- ------------------------------------------------------------------
c --- Read initial conditions from file to define layer temperature and
c --- salinity and geopotential at layer interfaces.
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dz
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(3) :: start,count
      integer i,j,k,l,status,ncid,dimid,varid
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading initial condition from ',
     .                    trim(icfile)
        call flush(lp)
c
c --- - Open netcdf file
        status=nf90_open(icfile,nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_open: ',trim(icfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
c
c --- - Check dimensions
        status=nf90_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: x: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=i)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: x: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        status=nf90_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_dimid: y: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        status=nf90_inquire_dimension(ncid,dimid,len=j)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inquire_dimension: y: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in '//trim(icfile)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
c
      endif
c
      start(1)=1
      start(2)=1
      count(1)=itdm
      count(2)=jtdm
      count(3)=1
c
c --- Read reference potential density
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'sigma',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: sigma: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: sigma: ',
     .                       nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,sigmar(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- Read potential temperature
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'temp',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: temp: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: temp: ',
     .                       nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,temp(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- Read salinity
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'saln',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: saln: ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: saln: ',
     .                       nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,saln(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- Read layer thickness
      if (mnproc.eq.1) then
        status=nf90_inq_varid(ncid,'dz',varid)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') ' nf90_inq_varid: dz: ',nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
      do k=1,kk
        if (mnproc.eq.1) then
          start(3)=k
          status=nf90_get_var(ncid,varid,tmp2d,start,count)
          if (status.ne.nf90_noerr) then
            write(lp,'(2a)') ' nf90_get_var: dz: ',nf90_strerror(status)
            call xchalt('(ictsz_file)')
                   stop '(ictsz_file)'
          endif
        endif
        call xcaput(tmp2d,dz(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf90_close(ncid)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') ' nf90_close: ',trim(icfile),': ',
     .                     nf90_strerror(status)
          call xchalt('(ictsz_file)')
                 stop '(ictsz_file)'
        endif
      endif
c
c --- Construct interface depths [cm] from layer thicknesses [m] and
c --- convert unit of reference potential density from [kg/m^3] to
c --- [g/cm^3]
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c         z(i,j,1)=z(i,j,1)*1.e2
          z(i,j,1)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      do k=1,kk
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            z(i,j,k+1)=min(depths(i,j)*1.e2,z(i,j,k)+dz(i,j,k)*1.e2)
            sigmar(i,j,k)=sigmar(i,j,k)*1.e-3
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=2,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (z(i,j,kk+1)-z(i,j,k).lt.1.e-4)
     .        z(i,j,k)=depths(i,j)*1.e2
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          z(i,j,kk+1)=depths(i,j)*1.e2
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute layer interface geopotential
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk+1
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            phi(i,j,k)=-g*z(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end subroutine ictsz_file
c
c --- ------------------------------------------------------------------
c --- Public procedures.
c --- ------------------------------------------------------------------
c
      subroutine inicon
c
c --- ------------------------------------------------------------------
c --- Define initial conditions
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tfrz
c
      integer i,j,k,l
      real q,tsfac,dps
c
c --- ------------------------------------------------------------------
c --- Define layer interface heights and layer temperature and salinity
c --- ------------------------------------------------------------------
c
      select case (trim(expcnf))
        case ('cesm', 'ben02clim', 'ben02syn', 'single_column')
          call ictsz_file
        case ('fuk95')
          call ictsz_fuk95
        case ('channel')
          call ictsz_channel
        case ('isomip1')
c         call ictsz_isomip1
        case ('isomip2')
c         call ictsz_isomip2
        case default
          if (mnproc.eq.1) then
             write (lp,'(3a)') ' inicon: expcnf = ', trim(expcnf),
     .                         ' is unsupported!'
          endif
          call xcstop('(inicon)')
                 stop '(inicon)'
      end select
c
c --- ------------------------------------------------------------------
c --- Set minimum physical temperature for each isopycnic layer
c --- ------------------------------------------------------------------
c
      call settemmin
c
c --- ------------------------------------------------------------------
c --- Initialize configuration specific variables
c --- ------------------------------------------------------------------
c
      select case (trim(expcnf))
        case ('cesm')
          call inicon_cesm
        case ('ben02clim', 'ben02syn', 'single_column')
          call inicon_ben02
      end select
c
c --- ------------------------------------------------------------------
c --- Make sure layer temperature is greater than the lower physical
c --- bound and make temperature, salinity, and potential density
c --- variables consistent.
c --- ------------------------------------------------------------------
c
      do k=1,2
        tfrz(:,:)=swtfrz(p(:,:,1),saln(:,:,k))
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            temp(i,j,k)=max(tfrz(i,j),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
      do k=3,kk
        tfrz(:,:)=swtfrz(p(:,:,1),saln(:,:,k))
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            temp(i,j,k)=max(tfrz(i,j),temp(i,j,k))
            saln(i,j,k)=sofsig(sigmar(i,j,k),temp(i,j,k))
            temp(i,j,k+kk)=temp(i,j,k)
            saln(i,j,k+kk)=saln(i,j,k)
            sigma(i,j,k)=sig(temp(i,j,k),saln(i,j,k))
            sigma(i,j,k+kk)=sigma(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      if (mnproc.eq.ptest) then
        write (lp,'('' sigmar(k)    :'',7f9.5/(15x,7f9.5))')
     .     (sigmar(itest,jtest,k),k=1,kk)
      endif
c
c --- ------------------------------------------------------------------
c --- Find layer interface pressure
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          p(i,j,1)=getpl(temp(i,j,1),saln(i,j,1),0.,phi(i,j,1),0.)
        enddo
        enddo
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            p(i,j,k+1)=getpl(temp(i,j,k),saln(i,j,k),
     .                       phi(i,j,k),phi(i,j,k+1),p(i,j,k))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(p, 1,kk+1, 2,2, halo_ps)
      call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 1,1, halo_ps)
      call xctilr(temp, 1,kk, 1,1, halo_ps)
      call xctilr(saln, 1,kk, 1,1, halo_ps)
c
c --- ------------------------------------------------------------------
c --- Set layer thickness and bottom pressure
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=0,jj+1
        do k=1,kk
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
            dp(i,j,k+kk)=dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=0,jj+1
        do k=1,kk
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          pb(i,j,1)=p(i,j,kk+1)
          pb(i,j,2)=pb(i,j,1)
          pb_mn(i,j,1)=pb(i,j,1)
          pb_mn(i,j,2)=pb(i,j,1)
          pb_p(i,j)=pb(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pbu(i,j,1)=min(pb(i,j,1),pb(i-1,j,1))
          pbu(i,j,2)=pbu(i,j,1)
          pbu_p(i,j)=pbu(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          pbv(i,j,1)=min(pb(i,j,1),pb(i,j-1,1))
          pbv(i,j,2)=pbv(i,j,1)
          pbv_p(i,j)=pbv(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i,q)
      do j=-1,jj+2
        do k=1,kk
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,k)=
     .        .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k)))
     .           +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,k)=
     .        .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k)))
     .           +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            dpu(i,j,k+kk)=dpu(i,j,k)
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            dpv(i,j,k+kk)=dpv(i,j,k)
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Initialize potential vorticity of barotropic flow
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,q)
      do j=0,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          q=2./(pb_p(i,j)+pb_p(i-1,j))
          pvtrop(i,j  ,1)=corioq(i,j  )*q
          pvtrop(i,j+1,1)=corioq(i,j+1)*q
          pvtrop(i,j  ,2)=pvtrop(i,j  ,1)
          pvtrop(i,j+1,2)=pvtrop(i,j+1,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q)
      do j=1,jj
        do l=1,isv(j)
        do i=max(0,ifv(j,l)),min(ii,ilv(j,l))
          q=2./(pb_p(i,j)+pb_p(i,j-1))
          pvtrop(i  ,j,1)=corioq(i  ,j)*q
          pvtrop(i+1,j,1)=corioq(i+1,j)*q
          pvtrop(i  ,j,2)=pvtrop(i  ,j,1)
          pvtrop(i+1,j,2)=pvtrop(i+1,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isq(j)
        do i=max(1,ifq(j,l)),min(ii,ilq(j,l))
          pvtrop(i,j,1)=corioq(i,j)*4./(pb_p(i,j  )+pb_p(i-1,j  )
     .                                 +pb_p(i,j-1)+pb_p(i-1,j-1))
          pvtrop(i,j,2)=pvtrop(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Separate baroclinic and barotropic velocity components
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
c
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ub(i,j,1)=0.
        enddo
        do k=1,kk
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            ub(i,j,1)=ub(i,j,1)+u(i,j,k)*dpu(i,j,k)
          enddo
        enddo
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ub(i,j,1)=ub(i,j,1)/pbu_p(i,j)
          ub(i,j,2)=ub(i,j,1)
        enddo
        enddo
c
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vb(i,j,1)=0.
        enddo
        do k=1,kk
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vb(i,j,1)=vb(i,j,1)+v(i,j,k)*dpv(i,j,k)
          enddo
        enddo
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vb(i,j,1)=vb(i,j,1)/pbv_p(i,j)
          vb(i,j,2)=vb(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            u(i,j,k)=u(i,j,k)-ub(i,j,1)
            u(i,j,k+kk)=u(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            v(i,j,k)=v(i,j,k)-vb(i,j,1)
            v(i,j,k+kk)=v(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      tsfac=delt1/dlt
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ubflx_mn(i,j,1)=ub(i,j,1)*pbu_p(i,j)*scuy(i,j)
          ubflx_mn(i,j,2)=ubflx_mn(i,j,1)
          ubflx(i,j,1)=ubflx_mn(i,j,1)
          ubflx(i,j,2)=ubflx_mn(i,j,1)
          ubflxs(i,j,1)=ubflx_mn(i,j,1)*tsfac
          ubflxs(i,j,2)=ubflxs(i,j,1)
          ubflxs(i,j,3)=ubflxs(i,j,1)
          ubflxs_p(i,j,1)=ubflxs(i,j,1)
          ubflxs_p(i,j,2)=ubflxs(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vbflx_mn(i,j,1)=vb(i,j,1)*pbv_p(i,j)*scvx(i,j)
          vbflx_mn(i,j,2)=vbflx_mn(i,j,1)
          vbflx_mn(i,j,2)=vbflx_mn(i,j,1)
          vbflx(i,j,1)=vbflx_mn(i,j,1)
          vbflx(i,j,2)=vbflx_mn(i,j,1)
          vbflxs(i,j,1)=vbflx_mn(i,j,1)*tsfac
          vbflxs(i,j,2)=vbflxs(i,j,1)
          vbflxs(i,j,3)=vbflxs(i,j,1)
          vbflxs_p(i,j,1)=vbflxs(i,j,1)
          vbflxs_p(i,j,2)=vbflxs(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ubcors_p(i,j)= (vbflx_mn(i  ,j  ,1)*scvxi(i  ,j  )
     .                   +vbflx_mn(i  ,j+1,1)*scvxi(i  ,j+1)
     .                   +vbflx_mn(i-1,j  ,1)*scvxi(i-1,j  )
     .                   +vbflx_mn(i-1,j+1,1)*scvxi(i-1,j+1))
     .                   *(pvtrop(i,j,1)+pvtrop(i,j+1,1))*.125*tsfac
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          ubcors_p(i,j)=-(ubflx_mn(i  ,j  ,1)*scuyi(i  ,j  )
     .                   +ubflx_mn(i+1,j  ,1)*scuyi(i+1,j  )
     .                   +ubflx_mn(i  ,j-1,1)*scuyi(i  ,j-1)
     .                   +ubflx_mn(i+1,j-1,1)*scuyi(i+1,j-1))
     .                   *(pvtrop(i,j,1)+pvtrop(i+1,j,1))*.125*tsfac
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Initialize fields related to the pressure gradient force
c --- ------------------------------------------------------------------
c
      call pgforc(2,1,kk,0,kk+1,1)
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            pgfx(i,j,k+kk)=pgfx(i,j,k)
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            pgfy(i,j,k+kk)=pgfy(i,j,k)
          enddo
          enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pgfxm(i,j,2)=pgfxm(i,j,1)
          xixp(i,j,2)=xixp(i,j,1)
          xixm(i,j,2)=xixm(i,j,1)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          pgfym(i,j,2)=pgfym(i,j,1)
          xiyp(i,j,2)=xiyp(i,j,1)
          xiym(i,j,2)=xiym(i,j,1)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Define first physical interior layer
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(l,i,k,dps)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          k=3
          dps=0.
          do while (dp(i,j,k).lt.epsil)
            dps=dps+dp(i,j,k)
            dp(i,j,k)=0.
            dp(i,j,k+kk)=0.
            k=k+1
            if (k.gt.kk) exit
          enddo
          if (k.gt.kk) then
            dp(i,j,2)=dp(i,j,2)+dps
            dp(i,j,2+kk)=dp(i,j,2)
          else
            dp(i,j,k)=dp(i,j,k)+dps
            dp(i,j,k+kk)=dp(i,j,k)
          endif
          kfpla(i,j,1)=k
          kfpla(i,j,2)=k
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Set other time level layer thicknesses
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,k+kk)=dp(i,j,k)
            dpold(i,j,k)=dp(i,j,k)
            dpold(i,j,k+kk)=dp(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (mnproc.eq.ptest) then
        i=itest
        j=jtest
        write (lp,103) nstep,i0+i,j0+j,
     .    '  init.profile  temp    saln    dens   thkns    dpth',
     .    (k,temp(i,j,k),saln(i,j,k),1000.*sig(temp(i,j,k),saln(i,j,k)),
     .    dp(i,j,k)/onem,p(i,j,k+1)/onem,k=1,kk)
 103    format (i9,2i5,a/(28x,i3,3f8.2,2f8.1))
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'inicon:'
        endif
        call chksummsk(p,ip,kk+1,'p')
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(pb,ip,3,'pb')
        call chksummsk(pbu,iu,2,'pbu')
        call chksummsk(pbv,iv,2,'pbv')
        call chksummsk(pvtrop,iq,2,'pvtrop')
        call chksummsk(pu,iu,kk+1,'pu')
        call chksummsk(pv,iv,kk+1,'pv')
      endif
c
      end subroutine inicon
c
      end module mod_inicon
