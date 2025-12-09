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

module mod_blom_step

  use dimensions,          only: nreg
  use mod_config,          only: expcnf
  use mod_time,            only: date, nday_of_year, nstep1, &
                                 nstep, nstep_in_day, delt1, step_time, &
                                 baclin
  use mod_timing,          only: timer_start, timer_stop, timer_reset, &
                                 timer_statistics, timer_group_time
  use mod_xc,              only: lp, kk, mnproc, xctilr, xcsum
  use mod_vcoord,          only: vcoord_tag, vcoord_isopyc_bulkml, sigref_adapt
  use mod_ale_regrid_remap, only: ale_regrid_remap, ale_remap_diazlv
  use mod_ale_vdiff,       only: ale_vdifft, ale_vdiffm
  use mod_swabs,           only: updswa
  use mod_tmsmt,           only: tmsmt1, tmsmt2
  use mod_eddtra,          only: eddtra
  use mod_advect,          only: advect
  use mod_pbcor,           only: pbcor1, pbcor2
  use mod_pgforc,          only: pgforc
  use mod_momtum,          only: momtum
  use mod_mxlayr,          only: mxlayr
  use mod_barotp,          only: barotp
  use mod_cmnfld_routines, only: cmnfld_bfsqi_ale, cmnfld1, cmnfld2
  use mod_forcing,         only: fwbbal
  use mod_budget,          only: budget_sums, budget_output
  use mod_eddtra,          only: eddtra
  use mod_momtum,          only: momtum
  use mod_difest,          only: difest_isobml, difest_lateral_hybrid, &
                                 difest_vertical_hybrid
  use mod_chkvar,          only: chkvar
  use mod_dia,             only: diaacc, diaout_alarms, diaout, &
                                 rstann, rstmon, rstfrq, diagmon_phy, &
                                 alarm_phy, nphy
  use mod_diapfl,          only: diapfl
  use mod_state,           only: temp, saln, dp, init_fluxes
  use mod_thermf,          only: thermf
  use mod_getfrc,          only: getfrc
  use mod_diffus,          only: diffus
  use mod_sfcstr,          only: sfcstr
  use mod_convec,          only: convec
  use mod_wdiflx,          only: wdiflx
  use mod_restart,         only: restart_write
  use mod_tracers_update,  only: updtrc
  use mod_ifdefs,          only: use_TRC
  use mod_ale_forcing,     only: ale_forcing

  implicit none
  private

  ! Public routines
  public :: blom_step

contains

  subroutine blom_step()
  ! ------------------------------------------------------------------
  ! integrate a model time step
  ! ------------------------------------------------------------------

    ! Local variables
    real    :: total_step_time
    integer :: i,m,n,mm,nn,k1m,k1n,iogrp
    logical :: update_flux_halos,wrtrst

    call timer_start('blom_step')

    ! letter 'm' refers to mid-time level (example: dp(i,j,km) )
    ! letter 'n' refers to old and new time level

    m = mod(nstep  ,2)+1
    n = mod(nstep+1,2)+1
    mm = (m-1)*kk
    nn = (n-1)*kk
    k1m = 1+mm
    k1n = 1+nn

    call budget_sums(1,n,nn)
    call timer_stop('blom_step','budget_sg1')

    call step_time
#ifndef OFFLINE_SEDIMENT_SPINUP
    ! ------------------------------------------------------------------
    ! Reset fluxes to be accumulated over a model time step and update
    ! flux halos the first time step of a day to reproduce results after
    ! restart with tripolar grid.
    ! ------------------------------------------------------------------

    update_flux_halos = nreg == 2 .and. mod(nstep,nstep_in_day) == 1
    call init_fluxes(m,n,mm,nn,k1m,k1n,update_flux_halos)

    call timer_stop('blom_step','init_fluxes')

    ! ------------------------------------------------------------------
    ! Get forcing
    ! ------------------------------------------------------------------

    call getfrc

    ! ------------------------------------------------------------------
    ! Update arrays related to shortwave radiation absorption.
    ! ------------------------------------------------------------------

    call updswa

    call timer_stop('blom_step','getfrc')

    call tmsmt1(nn)
    call timer_stop('blom_step','tmsmt1')

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      call ale_regrid_remap(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','ale_regrid_remap')
      call budget_sums(2,n,nn)
      call timer_stop('blom_step','budget_sg2')
    end if

    call cmnfld2(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','cmnfld2')

    if (vcoord_tag == vcoord_isopyc_bulkml) then
      call difest_isobml(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','difest_isobml')
    else
      call difest_lateral_hybrid(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','difest_lateral_hybrid')
    end if
    call eddtra(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','eddtra')
    call advect(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','advect')
    call pbcor1(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','pbcor1')
    call diffus(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','diffus')

    if (vcoord_tag == vcoord_isopyc_bulkml) then
      call budget_sums(2,n,nn)
      call timer_stop('blom_step','budget_sg3')
    else
      call budget_sums(3,n,nn)
      call timer_stop('blom_step','budget_sg4')
    end if

    call sfcstr(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','sfcstr')

    call pgforc(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','pgforc')

    call momtum(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','momtum')

    if (vcoord_tag == vcoord_isopyc_bulkml) then

      call convec(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','convec')

      call budget_sums(3,n,nn)
      call timer_stop('blom_step','budget_sg5')

      call diapfl(n,nn,k1n)
      call timer_stop('blom_step','diapfl')

      call budget_sums(4,n,nn)
      call timer_stop('blom_step','budget_sg6')

    end if

    call thermf(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','thermf')

    if (vcoord_tag == vcoord_isopyc_bulkml) then
      call mxlayr(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','mxlayr')
    else
      call cmnfld_bfsqi_ale(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','cmnfld_bfsqi_ale')
      call ale_forcing(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','ale_forcing')
      call difest_vertical_hybrid(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','difest_vertical_hybrid')
      call ale_vdifft(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','ale_vdifft')
      call ale_vdiffm(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','ale_vdiffm')
      call budget_sums(4,n,nn)
      call timer_stop('blom_step','budget_sg7')
    end if
#endif

    if (use_TRC) then
      ! update tracer due to non-passive processes
      call updtrc(m,n,mm,nn,k1m,k1n)
      call timer_stop('blom_step','updtrc')
    end if

#ifndef OFFLINE_SEDIMENT_SPINUP
    call budget_sums(5,n,nn)
    call timer_stop('blom_step','budget_sg8')

    call barotp(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','barotp')

    call pbcor2(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','pbcor2')

    call budget_sums(6,m,mm)
    call timer_stop('blom_step','budget_sg9')

    call tmsmt2(m,mm,nn,k1m)
    call timer_stop('blom_step','tmsmt2')

    call budget_sums(7,m,mm)
    call timer_stop('blom_step','budget_sg10')

    call cmnfld1(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','cmnfld1')

    call sigref_adapt(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','sigref_adapt')

    call diaacc(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','diaacc')

    call fwbbal(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','fwbbal')

    call budget_output(m)
    call timer_stop('blom_step','budget_sg11')

    !-------------------------------------------------------------------
    ! output and diagnostic calculations
    !-------------------------------------------------------------------

    call chkvar(m,n,mm,nn,k1m,k1n)
    call timer_stop('blom_step','chkvar')

    if (mod(nstep,nstep_in_day) == 0.and.nday_of_year == 1) then

      ! ------------------------------------------------------------------
      ! output diagnosed heat and salt flux
      ! ------------------------------------------------------------------

      call wdiflx

    end if

    ! ------------------------------------------------------------------
    ! output of BLOM diagnostics
    ! ------------------------------------------------------------------

    call diaout_alarms

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      call ale_remap_diazlv(m,n,mm,nn,k1m,k1n)
    end if

    call diaout(m,n,mm,nn,k1m,k1n)

    call timer_stop('blom_step','diaout')

    wrtrst = .false.
    if (((rstann.and.nday_of_year == 1.or.rstmon.and.date%day == 1) &
         .and.mod(nstep,nstep_in_day) == 0).or. &
         .not.(rstann.or.rstmon).and.mod(nstep+.5,rstfrq) < 1.) then

      if (expcnf /= 'cesm') then

        ! ------------------------------------------------------------------
        ! output restart files
        ! ------------------------------------------------------------------

        wrtrst = .true.
        call restart_write()

      end if

    else
    end if
#endif

    if (expcnf /= 'cesm') call timer_stop('blom_step','restart_write')

    delt1 = baclin + baclin

    ! --------------------------------------------------------------------
    ! write timing diagnostics to standard out
    ! --------------------------------------------------------------------

#ifndef OFFLINE_SEDIMENT_SPINUP
    call timer_stop('total_step_time')
    call timer_group_time('total_step_time',total_step_time)
    call timer_reset('total_step_time')
    call timer_start('total_step_time')

    if (mnproc == 1) then
      write (lp,'(f12.4,a,i8)') total_step_time, ' sec for step ', nstep
    end if

    if (expcnf /= 'cesm') then
      do iogrp = 1, nphy
        if (diagmon_phy(iogrp) .and. (alarm_phy(iogrp) == 1 .or. wrtrst)) then
          call timer_statistics('blom_step')
          exit
        endif
      enddo
    endif
#endif

  end subroutine blom_step

end module mod_blom_step
