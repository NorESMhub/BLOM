! ------------------------------------------------------------------------------
! Copyright (C) 2008-2024 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_blom_init

  implicit none
  private

  public :: blom_init

  private :: numerical_bounds

contains

  subroutine blom_init()
  ! ------------------------------------------------------------------
  ! initialize the model
  ! ------------------------------------------------------------------

    use dimensions,          only: itdm, nreg
    use mod_config,          only: expcnf
    use mod_time,            only: date, nday1, nday2, nstep1, nstep2, nstep, delt1, &
                                   time0, baclin
    use mod_timing,          only: init_timing, get_time
    use mod_xc,              only: xcspmd, xcbcst, xctilr, mnproc, nproc, &
                                   lp, ii, jj, kk, isp, ifp, isu, ifu, ilp, isv, ifv, &
                                   ilu, ilv, jpr, i0, nbdy, &
                                   halo_ps, halo_us, halo_vs, halo_uv, halo_vv, halo_qs
    use mod_pointtest,       only: init_ptest
    use mod_inicon,          only: icfile
    use mod_state,           only: dp, dpu, dpv, uflx, vflx, p, pu, pv, phi
    use mod_barotp,          only: pvtrop
    use mod_pgforc,          only: pgfxm, pgfym, xixp, xixm, xiyp, xiym
    use mod_niw,             only: uml, vml, umlres, vmlres
    use mod_eos,             only: inieos
    use mod_swabs,           only: iniswa
    use mod_tmsmt,           only: initms
    use mod_dia,             only: diaini, diasg1
    use mod_inicon,          only: inicon
    use mod_budget,          only: budget_init
    use mod_cmnfld_routines, only: cmnfld1
    use mod_tke,             only: initke
    use mod_rdlim,           only: rdlim
    use mod_inifrc,          only: inifrc
    use mod_inivar,          only: inivar
    use mod_vcoord,          only: vcoord_type_tag, isopyc_bulkml, sigmar
    use mod_inigeo,          only: inigeo
    use mod_iniphy,          only: iniphy
    use mod_restart,         only: restart_read
    use mod_ifdefs,          only: use_TRC, use_TKE
    use mod_tracers_update,  only: initrc
    use netcdf

    ! Local variables
    integer :: istat,ncid,varid,i,j,k,l,m,n,mm,nn,k1m,k1n,mt,mmt,kn,km
    real    :: q
    logical :: icrest,fexist
    integer :: icrest_int

    ! ---------------------------------------------------------------
    ! Initialize SPMD processing
    ! ------------------------------------------------------------------

    call xcspmd

    ! ------------------------------------------------------------------
    ! Initialize timing
    ! ------------------------------------------------------------------

    call init_timing

    ! print seconds elapsed since startup (should be almost zero)
    if (mnproc == 1) then
      write (lp,'(f12.4,a,i8)') get_time(),' Time 0 BLOM starting up'
      call flush(lp)
    end if

    ! ------------------------------------------------------------------
    ! Initialize check sum algorithm
    ! ------------------------------------------------------------------

    call crcinit

    ! ------------------------------------------------------------------
    ! Read limits file
    ! ------------------------------------------------------------------

    call rdlim

    ! ------------------------------------------------------------------
    ! Identify processor and horizontal indexes where detailed
    ! diagnostics are desired
    ! ------------------------------------------------------------------

    call init_ptest

    ! ------------------------------------------------------------------
    ! Initialize the geographic environment
    ! ------------------------------------------------------------------

    call inigeo

    ! ------------------------------------------------------------------
    ! Initialize various arrays
    ! ------------------------------------------------------------------

    call inivar

    ! ------------------------------------------------------------------
    ! Set various numerical bounds
    ! ------------------------------------------------------------------

    call numerical_bounds

    ! ------------------------------------------------------------------
    ! Initialize physical parameterizations
    ! ------------------------------------------------------------------

    call iniphy

    ! ------------------------------------------------------------------
    ! Initialize forcing
    ! ------------------------------------------------------------------

    call inifrc

    ! ------------------------------------------------------------------
    ! Define coefficients for equation of state functions
    ! ------------------------------------------------------------------

    call inieos

    ! ------------------------------------------------------------------
    ! Initialize shortwave radiation absorption
    ! ------------------------------------------------------------------

    call iniswa

    ! ------------------------------------------------------------------
    ! Initialize second order turbulence closure closure
    ! ------------------------------------------------------------------

    if (use_TRC .and. use_TKE) then
      call initke
    end if

    ! ------------------------------------------------------------------
    ! Initialize diagnostic accumulation fields
    ! ------------------------------------------------------------------

    call diaini

    ! ------------------------------------------------------------------
    ! Set up initial conditions or start from restart file
    ! ------------------------------------------------------------------

    ! check whether initial condition file given in namelist is a
    ! restart file
    icrest = .false.
    icrest_int = 0
    if (mnproc == 1) then
      inquire(file=icfile,exist = fexist)
      if (fexist) then
        istat = nf90_open(icfile,nf90_nowrite,ncid)
        if (istat == nf90_noerr) then
          istat = nf90_inq_varid(ncid,'dp',varid)
          if (istat == nf90_noerr) then
            icrest = .true.
          end if
        end if
      end if
      if (icrest) icrest_int = 1
    end if
    call xcbcst(icrest_int)
    icrest = (icrest_int == 1)

    if (nday1+nint(time0) == 0.and..not.icrest) then

      ! ----------------------------------------------------------------
      ! start from initial conditions derived from climatology
      ! ----------------------------------------------------------------

      if (date%month /= 1.or.date%day /= 1) then
        if (mnproc == 1) then
          write (lp,*) &
               'Warning! date is inconsistent with ini. cond. (Jan 1st)!'
          call flush(lp)
        end if
      end if

      delt1 = baclin

      call inicon
      if (use_TRC) then
        call initrc
      end if

    else ! nday1+nint(time0) > 0 .or. icrest

      ! ------------------------------------------------------------------
      ! start from restart file
      ! ------------------------------------------------------------------

      delt1 = baclin+baclin

      call restart_read()

    end if

    ! ------------------------------------------------------------------
    ! Initialize model time step and set time level indices consistent
    ! with starting state
    ! ------------------------------------------------------------------

    nstep = nstep1
    m = mod(nstep+1,2)+1
    n = mod(nstep  ,2)+1
    mm = (m-1)*kk
    nn = (n-1)*kk
    k1m = 1+mm
    k1n = 1+nn

    ! ------------------------------------------------------------------
    ! Initialize layer thicknesses
    ! ------------------------------------------------------------------

    call xctilr(dp, 1,2*kk, 3,3, halo_ps)

    if (vcoord_type_tag == isopyc_bulkml) then

      do mt = n,3-n,3-2*n
        mmt = (mt-1)*kk

        !$omp parallel do private(k,l,i)
        do j = -2,jj+2
          do k = 1,kk
            do l = 1,isp(j)
              do i = max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
                p(i,j,k+1) = p(i,j,k)+dp(i,j,k+mmt)
              end do
            end do
          end do
        end do
        !$omp end parallel do

        !$omp parallel do private(k,km,l,i,q)
        do j = -1,jj+2
          do k = 1,kk
            km = k+mmt
            do l = 1,isu(j)
              do i = max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
                q = min(p(i,j,kk+1),p(i-1,j,kk+1))
                dpu(i,j,km)= &
                     .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k))) &
                     +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
              end do
            end do
            do l = 1,isv(j)
              do i = max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
                q = min(p(i,j,kk+1),p(i,j-1,kk+1))
                dpv(i,j,km)= &
                     .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k))) &
                     +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
              end do
            end do
          end do
        end do
        !$omp end parallel do

      end do

    else

      call xctilr(dpu, 1,2*kk, 3,3, halo_us)
      call xctilr(dpv, 1,2*kk, 3,3, halo_vs)

      !$omp parallel do private(k,km,l,i)
      do j = -2,jj+2
        do k = 1,kk
          km = k+mm
          do l = 1,isp(j)
            do i = max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
              p(i,j,k+1) = p(i,j,k)+dp(i,j,km)
            end do
          end do
        end do
      end do
      !$omp end parallel do

      !$omp parallel do private(k,kn,l,i)
      do j = 1,jj
        do k = 1,kk
          kn = k+nn
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              pu(i,j,k+1) = pu(i,j,k)+dpu(i,j,kn)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              pv(i,j,k+1) = pv(i,j,k)+dpv(i,j,kn)
            end do
          end do
        end do
      end do
      !$omp end parallel do

    end if

    ! ------------------------------------------------------------------
    ! initialize budget calculations
    ! ------------------------------------------------------------------

    call budget_init

    ! ------------------------------------------------------------------
    ! update some halos
    ! ------------------------------------------------------------------

    call xctilr(sigmar, 1,kk, 2,2, halo_ps)
    call xctilr(uflx, 1,2*kk, 1,1, halo_uv)
    call xctilr(vflx, 1,2*kk, 1,1, halo_vv)
    call xctilr(phi(1-nbdy,1-nbdy,kk+1), 1,1, 2,2, halo_ps)
    call xctilr(pvtrop, 1,2, 1,3, halo_qs)
    call xctilr(pgfxm, 1,2, 1,2, halo_uv)
    call xctilr(xixp, 1,2, 1,2, halo_us)
    call xctilr(xixm, 1,2, 1,2, halo_us)
    call xctilr(pgfym, 1,2, 1,2, halo_vv)
    call xctilr(xiyp, 1,2, 1,2, halo_vs)
    call xctilr(xiym, 1,2, 1,2, halo_vs)
    if (vcoord_type_tag == isopyc_bulkml) then
       call xctilr(uml, 1,4, 1,0, halo_uv)
       call xctilr(vml, 1,4, 0,1, halo_vv)
       call xctilr(umlres, 1,2, 1,0, halo_uv)
       call xctilr(vmlres, 1,2, 0,1, halo_vv)
    end if

    ! with arctic patch, switch xixp and xixm and xiyp and xiym in the
    ! halo region adjacent to the arctic grid intersection
    if (nreg == 2.and.nproc == jpr) then
      do j = jj,jj+2
        do i = 0,ii+1
          q = xixp(i,j,1)
          xixp(i,j,1) = xixm(i,j,1)
          xixm(i,j,1) = q
          q = xixp(i,j,2)
          xixp(i,j,2) = xixm(i,j,2)
          xixm(i,j,2) = q
        end do
      end do
      do i = max(0,itdm/2-i0+1),ii+1
        q = xiyp(i,jj,1)
        xiyp(i,jj,1) = xiym(i,jj,1)
        xiym(i,jj,1) = q
        q = xiyp(i,jj,2)
        xiyp(i,jj,2) = xiym(i,jj,2)
        xiym(i,jj,2) = q
      end do
      do j = jj+1,jj+2
        do i = 0,ii+1
          q = xiyp(i,j,1)
          xiyp(i,j,1) = xiym(i,j,1)
          xiym(i,j,1) = q
          q = xiyp(i,j,2)
          xiyp(i,j,2) = xiym(i,j,2)
          xiym(i,j,2) = q
        end do
      end do
    end if

    ! ------------------------------------------------------------------
    ! Initialize time smoothing variables and some common fields.
    ! ------------------------------------------------------------------

    call initms(mm)
    call cmnfld1(m,n,mm,nn,k1m,k1n)

    ! ------------------------------------------------------------------
    ! Extract reference potential density vector representative of the
    ! dominating ocean domain
    ! ------------------------------------------------------------------

    call diasg1

    ! ------------------------------------------------------------------

    if (mnproc == 1.and.expcnf /= 'cesm') then
      write (lp,'(/2(a,i6),2(a,i9),a/)') &
           'model starts at day',nday1,', goes to day',nday2,'   (steps', &
           nstep1,' --',nstep2,')'
      call flush(lp)
    end if

    ! print seconds elapsed since last call to system_clock (Time 0)
    if (mnproc == 1) then
      write (lp,'(f12.4,a,i8)') get_time(),' Time 1 Just before main loop'
      call flush(lp)
    end if

  end subroutine blom_init

  subroutine numerical_bounds
  !------------------------------------------------------------------------
  ! Set various numerical bounds.
  !------------------------------------------------------------------------

    use mod_types,     only: r8
    use mod_constants, only: g, spval, L_mks2cgs
    use mod_time,      only: baclin
    use mod_xc
    use mod_grid,      only: scqx, scqy, scpx, scpy, scuy, scvx, scp2, depths
    use mod_diffusion, only: difmxp, difmxq
    use mod_utility,   only: umax, vmax
    use mod_checksum,  only: csdiag, chksummsk

    ! Local variables
    real(r8) :: dx2, dy2, btdtmx, umaxmin, vmaxmin, umaxmax, vmaxmax
    integer :: i, j, l

    ! Determine upper bound of lateral diffusivity based on numerical stability
    ! concerns.

    !$omp parallel do private(i, dx2, dy2)
    do j = 1 - nbdy, jj + nbdy
      do i = 1 - nbdy, ii + nbdy
        dx2 = scpx(i, j)*scpx(i, j)
        dy2 = scpy(i, j)*scpy(i, j)
        difmxp(i, j) = .9_r8*.5_r8*dx2*dy2 &
                       /max(1._r8,(dx2 + dy2)*(baclin + baclin))
        dx2 = scqx(i, j)*scqx(i, j)
        dy2 = scqy(i, j)*scqy(i, j)
        difmxq(i, j) = .9_r8*.5_r8*dx2*dy2 &
                       /max(1._r8,(dx2 + dy2)*(baclin + baclin))
      enddo
    enddo
    !$omp end parallel do

    ! Estimate maximum barotropic time step.
    btdtmx = 86400._r8

    !$omp parallel do private(l, i) reduction(min:btdtmx)
    do j = 1, jj
      do l = 1, isp(j)
        do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
          btdtmx = min(btdtmx, &
                      scpx(i, j)*scpy(i, j) &
                      /sqrt(g*depths(i, j)*L_mks2cgs*( scpx(i, j)*scpx(i, j) &
                      + scpy(i, j)*scpy(i, j))))
        enddo
      enddo
    enddo
    !$omp end parallel do
    call xcmin(btdtmx)
    if (mnproc == 1) then
      write (lp, *) 'estimated max. barotropic time step:', btdtmx/sqrt(2._r8)
      call flush(lp)
    endif

    ! Set maximum velocities allowable ensuring stability of the upwind scheme.

    umaxmin = spval
    vmaxmin = spval
    umaxmax = 0._r8
    vmaxmax = 0._r8
    !$omp parallel do private(l, i) &
    !$omp reduction(min:umaxmin, vmaxmin) reduction(max:umaxmax, vmaxmax)
    do j = 1, jj
      do l = 1, isu(j)
        do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
          umax(i, j) = .9_r8*.125_r8*min(scp2(i - 1, j), scp2(i, j)) &
                       /(scuy(i, j)*baclin)
          umaxmin = min(umaxmin, umax(i, j))
          umaxmax = max(umaxmax, umax(i, j))
        enddo
      enddo
      do l = 1, isv(j)
        do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
          vmax(i, j) = .9_r8*.125_r8*min(scp2(i, j - 1), scp2(i, j)) &
                       /(scvx(i, j)*baclin)
          vmaxmin = min(vmaxmin, vmax(i, j))
          vmaxmax = max(vmaxmax, vmax(i, j))
        enddo
      enddo
    enddo
    !$omp end parallel do

    call xctilr(umax, 1, 1, nbdy, nbdy, halo_us)
    call xctilr(vmax, 1, 1, nbdy, nbdy, halo_vs)

    call xcmin(umaxmin)
    call xcmax(umaxmax)
    call xcmin(vmaxmin)
    call xcmax(vmaxmax)

    if (mnproc == 1) then
      write (lp, *) 'min/max umax:', umaxmin, umaxmax
      write (lp, *) 'min/max vmax:', vmaxmin, vmaxmax
      call flush(lp)
    endif

    if (csdiag) then
      if (mnproc == 1) then
        write (lp, *) 'numerical_bounds:'
      endif
      call chksummsk(difmxp, ip, 1,'difmxp')
      call chksummsk(difmxq, iq, 1,'difmxq')
      call chksummsk(umax, iu, 1,'umax')
      call chksummsk(vmax, iv, 1,'vmax')
    endif

  end subroutine numerical_bounds

end module mod_blom_init
