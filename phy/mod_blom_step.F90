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
  use mod_timing,          only: total_time, total_xio_time, &
                                 auxil_total_time, getfrc_total_time, &
                                 tmsmt1_total_time, advdif_total_time, &
                                 sfcstr_total_time, momtum_total_time, &
                                 pgforc_total_time, barotp_total_time, &
                                 pbcor2_total_time, convec_total_time, &
                                 diapfl_total_time, thermf_total_time, &
                                 mxlayr_total_time, tmsmt2_total_time, &
                                 diaacc_total_time, io_total_time, &
                                 get_time
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
                                 rstann, rstmon, rstfrq
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
    real    :: q
    integer :: i,m,n,mm,nn,k1m,k1n
    logical :: update_flux_halos
    real    :: total_step_time
    real    :: auxil_time
    real    :: getfrc_time
    real    :: tmsmt1_time
    real    :: advdif_time
    real    :: sfcstr_time
    real    :: momtum_time
    real    :: pgforc_time
    real    :: barotp_time
    real    :: pbcor2_time
    real    :: convec_time
    real    :: diapfl_time
    real    :: thermf_time
    real    :: mxlayr_time
    real    :: tmsmt2_time
    real    :: diaacc_time
    real    :: io_time

    ! letter 'm' refers to mid-time level (example: dp(i,j,km) )
    ! letter 'n' refers to old and new time level

    m = mod(nstep  ,2)+1
    n = mod(nstep+1,2)+1
    mm = (m-1)*kk
    nn = (n-1)*kk
    k1m = 1+mm
    k1n = 1+nn

    call budget_sums(1,n,nn)

    call step_time

    ! ------------------------------------------------------------------
    ! Reset fluxes to be accumulated over a model time step and update
    ! flux halos the first time step of a day to reproduce results after
    ! restart with tripolar grid.
    ! ------------------------------------------------------------------

    update_flux_halos = nreg == 2 .and. mod(nstep,nstep_in_day) == 1
    call init_fluxes(m,n,mm,nn,k1m,k1n,update_flux_halos)

    auxil_time = get_time()

    ! ------------------------------------------------------------------
    ! Get forcing
    ! ------------------------------------------------------------------

    call getfrc

    ! ------------------------------------------------------------------
    ! Update arrays related to shortwave radiation absorption.
    ! ------------------------------------------------------------------

    call updswa

    getfrc_time = get_time()

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      call ale_regrid_remap(m,n,mm,nn,k1m,k1n)
      convec_time = get_time()
      call budget_sums(2,n,nn)
    end if

    call cmnfld2(m,n,mm,nn,k1m,k1n)

    !diag write (lp,*) 'tmsmt1...'
    call tmsmt1(nn)
    tmsmt1_time = get_time()

    !diag write (lp,*) 'advdif...'
    if (vcoord_tag == vcoord_isopyc_bulkml) then
      call difest_isobml(m,n,mm,nn,k1m,k1n)
    else
      call difest_lateral_hybrid(m,n,mm,nn,k1m,k1n)
    end if
    call eddtra(m,n,mm,nn,k1m,k1n)
    call advect(m,n,mm,nn,k1m,k1n)
    call pbcor1(m,n,mm,nn,k1m,k1n)
    call diffus(m,n,mm,nn,k1m,k1n)
    advdif_time = get_time()

    if (vcoord_tag == vcoord_isopyc_bulkml) then
      call budget_sums(2,n,nn)
    else
      call budget_sums(3,n,nn)
    end if
    auxil_time = auxil_time+get_time()

    !diag write (lp,*) 'sfcstr...'
    call sfcstr(m,n,mm,nn,k1m,k1n)
    sfcstr_time = get_time()

    !diag write (lp,*) 'pgforc...'
    call pgforc(m,n,mm,nn,k1m,k1n)
    pgforc_time = get_time()

    !diag write (lp,*) 'momtum...'
    call momtum(m,n,mm,nn,k1m,k1n)
    momtum_time = get_time()

    if (vcoord_tag == vcoord_isopyc_bulkml) then

      !diag   write (lp,*) 'convec...'
      call convec(m,n,mm,nn,k1m,k1n)
      convec_time = get_time()

      call budget_sums(3,n,nn)
      auxil_time = auxil_time+get_time()

      !diag   write (lp,*) 'diapfl...'
      call diapfl(n,nn,k1n)
      diapfl_time = get_time()

      call budget_sums(4,n,nn)
      auxil_time = auxil_time+get_time()

    end if

    !diag write (lp,*) 'thermf...'
    call thermf(m,n,mm,nn,k1m,k1n)
    thermf_time = get_time()

    if (vcoord_tag == vcoord_isopyc_bulkml) then
      !diag   write (lp,*) 'mxlayr...'
      call mxlayr(m,n,mm,nn,k1m,k1n)
      mxlayr_time = get_time()
    else
      call cmnfld_bfsqi_ale(m,n,mm,nn,k1m,k1n)
      call ale_forcing(m,n,mm,nn,k1m,k1n)
      call difest_vertical_hybrid(m,n,mm,nn,k1m,k1n)
      mxlayr_time = get_time()
      call ale_vdifft(m,n,mm,nn,k1m,k1n)
      call ale_vdiffm(m,n,mm,nn,k1m,k1n)
      call budget_sums(4,n,nn)
      diapfl_time = get_time()
    end if

    if (use_TRC) then
      ! update tracer due to non-passive processes
      call updtrc(m,n,mm,nn,k1m,k1n)
    end if

    call budget_sums(5,n,nn)
    auxil_time = auxil_time+get_time()

    !diag write (lp,*) 'barotp...'
    call barotp(m,n,mm,nn,k1m,k1n)
    barotp_time = get_time()

    !diag write (lp,*) 'pbcor2...'
    call pbcor2(m,n,mm,nn,k1m,k1n)
    pbcor2_time = get_time()

    call budget_sums(6,m,mm)
    auxil_time = auxil_time+get_time()

    !diag write (lp,*) 'tmsmt2...'
    call tmsmt2(m,mm,nn,k1m)
    tmsmt2_time = get_time()

    call budget_sums(7,m,mm)

    call cmnfld1(m,n,mm,nn,k1m,k1n)

    call sigref_adapt(m,n,mm,nn,k1m,k1n)

    call diaacc(m,n,mm,nn,k1m,k1n)
    diaacc_time = get_time()

    call fwbbal(m,n,mm,nn,k1m,k1n)

    call budget_output(m)

    auxil_time = auxil_time+get_time()

    !-------------------------------------------------------------------
    ! output and diagnostic calculations
    !-------------------------------------------------------------------

    call chkvar(m,n,mm,nn,k1m,k1n)

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

    ! update total time spent by various tasks
    auxil_total_time = auxil_total_time+auxil_time
    getfrc_total_time = getfrc_total_time+getfrc_time
    tmsmt1_total_time = tmsmt1_total_time+tmsmt1_time
    advdif_total_time = advdif_total_time+advdif_time
    sfcstr_total_time = sfcstr_total_time+sfcstr_time
    momtum_total_time = momtum_total_time+momtum_time
    pgforc_total_time = pgforc_total_time+pgforc_time
    barotp_total_time = barotp_total_time+barotp_time
    pbcor2_total_time = pbcor2_total_time+pbcor2_time
    convec_total_time = convec_total_time+convec_time
    diapfl_total_time = diapfl_total_time+diapfl_time
    thermf_total_time = thermf_total_time+thermf_time
    mxlayr_total_time = mxlayr_total_time+mxlayr_time
    tmsmt2_total_time = tmsmt2_total_time+tmsmt2_time
    diaacc_total_time = diaacc_total_time+diaacc_time

    if (((rstann.and.nday_of_year == 1.or.rstmon.and.date%day == 1) &
         .and.mod(nstep,nstep_in_day) == 0).or. &
         .not.(rstann.or.rstmon).and.mod(nstep+.5,rstfrq) < 1.) then

      if (expcnf /= 'cesm') then

        ! ------------------------------------------------------------------
        ! output restart files
        ! ------------------------------------------------------------------

        call restart_write()

      end if

      io_time = get_time()

      ! ------------------------------------------------------------------
      ! write timing diagnostics to standard out
      ! ------------------------------------------------------------------

      io_total_time = io_total_time+io_time
      total_step_time = auxil_time +getfrc_time+tmsmt1_time+advdif_time &
                       +sfcstr_time+momtum_time+pgforc_time+barotp_time &
                       +pbcor2_time+convec_time+diapfl_time+thermf_time &
                       +mxlayr_time+tmsmt2_time+diaacc_time+io_time
      total_time = total_time+total_step_time
      total_xio_time = total_xio_time+total_step_time-io_time

      if (mnproc == 1) then
        write (lp,'(f12.4,a,i8)') total_step_time, '  sec for step ', nstep
        write (lp,'(f12.4,a,i8)') total_time/(nstep-nstep1),' Avg Time'
        write (lp,'(f12.4,a,i8)') total_xio_time/(nstep-nstep1),' Avg Time excluding IO'
        write (lp,'(f12.4,a,i8)') total_time,' Tot Time with contributions:'
        q = 100./total_time
        write (lp,'(f12.4,a,i8)') auxil_total_time*q ,'% auxil '
        write (lp,'(f12.4,a,i8)') getfrc_total_time*q,'% getfrc'
        write (lp,'(f12.4,a,i8)') tmsmt1_total_time*q,'% tmsmt1'
        write (lp,'(f12.4,a,i8)') advdif_total_time*q,'% advdif'
        write (lp,'(f12.4,a,i8)') sfcstr_total_time*q,'% sfcstr'
        write (lp,'(f12.4,a,i8)') momtum_total_time*q,'% momtum'
        write (lp,'(f12.4,a,i8)') pgforc_total_time*q,'% pgforc'
        write (lp,'(f12.4,a,i8)') barotp_total_time*q,'% barotp'
        write (lp,'(f12.4,a,i8)') pbcor2_total_time*q,'% pbcor2'
        write (lp,'(f12.4,a,i8)') convec_total_time*q,'% convec'
        write (lp,'(f12.4,a,i8)') diapfl_total_time*q,'% diapfl'
        write (lp,'(f12.4,a,i8)') thermf_total_time*q,'% thermf'
        write (lp,'(f12.4,a,i8)') mxlayr_total_time*q,'% mxlayr'
        write (lp,'(f12.4,a,i8)') tmsmt2_total_time*q,'% tmsmt2'
        write (lp,'(f12.4,a,i8)') diaacc_total_time*q,'% diaacc'
        write (lp,'(f12.4,a,i8)') io_total_time*q    ,'% IO'
      end if

    else

      ! ------------------------------------------------------------------
      ! write time spent for current time step
      ! ------------------------------------------------------------------

      io_time = get_time()
      io_total_time = io_total_time+io_time
      total_step_time = auxil_time  + getfrc_time + tmsmt1_time + advdif_time &
                        + sfcstr_time + momtum_time + pgforc_time + barotp_time &
                        + pbcor2_time + convec_time + diapfl_time + thermf_time &
                        + mxlayr_time + tmsmt2_time + diaacc_time + io_time
      total_time = total_time + total_step_time
      total_xio_time = total_xio_time + total_step_time-io_time

      if (mnproc == 1) then
        write (lp,'(f12.4,a,i8)') total_step_time, '  sec for step ', nstep
      end if

    end if

    delt1 = baclin + baclin

  end subroutine blom_step

end module mod_blom_step
