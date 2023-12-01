! ------------------------------------------------------------------------------
! Copyright (C) 2006-2023 Mats Bentsen, Mehmet Ilicak, Alok Kumar Gupta,
!                         Jerry Tjiputra, Ping-Gin Chiu, Aleksi Nummelin,
!                         Jörg Schwinger

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

module mod_restart_rd

  use dimensions,    only: idm, jdm
  use mod_config,    only: expcnf, runid, inst_suffix, resume_flag
  use mod_calendar,  only: date_type, daynum_diff, operator(/=), &
                           calendar_noerr, calendar_errstr
  use mod_time,      only: date0, date, nday1, nstep0, nstep1, &
                           time0, calendar, time
  use mod_xc,        only: nbdy, nfu, ii, jj, kk, ip, lp, iu, iv, &
                           halo_uv, halo_vv, mnproc, &
                           xcbcst, xctilr, xcstop
  use mod_vcoord,    only: vcoord_type_tag, isopyc_bulkml, &
                           cntiso_hybrid, sigmar
  use mod_inicon,    only: icfile
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                           uflx, vflx, utflx, vtflx, usflx, vsflx, &
                           phi, ubflxs, vbflxs, &
                           ub, vb, pb, pbu, pbv, ubflxs_p, vbflxs_p, &
                           pb_p, pbu_p, pbv_p, ubcors_p, vbcors_p, &
                           sealv, kfpla
  use mod_pgforc,    only: pgfx, pgfy, pgfxm, pgfym, &
                           xixp, xixm, xiyp, xiym
  use mod_barotp,    only: ubflx, vbflx, pb_mn, ubflx_mn, vbflx_mn, &
                           pvtrop
  use mod_dia,       only: iotype, phyh2d, phylyr, phylvl, &
                           acc_abswnd ,acc_alb ,acc_brnflx ,acc_brnpd ,acc_dfl, &
                           acc_eva ,acc_fice ,acc_fmltfz ,acc_hice ,acc_hmltfz, &
                           acc_hsnw ,acc_iage ,acc_idkedt ,acc_lamult ,acc_lasl, &
                           acc_lip ,acc_maxmld ,acc_mld ,acc_mlts ,acc_mltsmn, &
                           acc_mltsmx ,acc_mltssq ,acc_mtkeus ,acc_mtkeni ,acc_mtkebf, &
                           acc_mtkers ,acc_mtkepe ,acc_mtkeke ,acc_mty ,acc_nsf, &
                           acc_pbot ,acc_psrf ,acc_rfiflx ,acc_rnfflx ,acc_salflx, &
                           acc_salrlx ,acc_sbot ,acc_sealv ,acc_slvsq ,acc_sfl, &
                           acc_sop ,acc_sigmx ,acc_sss ,acc_ssssq ,acc_sst, &
                           acc_sstsq ,acc_surflx ,acc_surrlx ,acc_swa ,acc_t20d, &
                           acc_taux ,acc_tauy ,acc_tbot ,acc_tice ,acc_tsrf, &
                           acc_ub ,acc_ubflxs ,acc_uice ,acc_ustar ,acc_ustar3, &
                           acc_ustokes,acc_vb ,acc_vbflxs ,acc_vice ,acc_vstokes, &
                           acc_ztx ,acc_ivolu ,acc_ivolv ,acc_utilh2d, &
                           acc_bfsq ,acc_difdia ,acc_difvmo ,acc_difvho ,acc_difvso, &
                           acc_difint ,acc_difiso ,acc_dp ,acc_dpu ,acc_dpv, &
                           acc_dz ,acc_saln ,acc_temp ,acc_uflx ,acc_utflx, &
                           acc_usflx ,acc_umfltd ,acc_umflsm ,acc_utfltd ,acc_utflsm, &
                           acc_utflld ,acc_usfltd ,acc_usflsm ,acc_usflld ,acc_uvel, &
                           acc_vflx ,acc_vtflx ,acc_vsflx ,acc_vmfltd ,acc_vmflsm, &
                           acc_vtfltd ,acc_vtflsm ,acc_vtflld ,acc_vsfltd ,acc_vsflsm, &
                           acc_vsflld ,acc_vvel ,acc_wflx ,acc_wflx2 ,acc_avdsg, &
                           acc_dpvor ,acc_tke ,acc_gls_psi,acc_utillyr, &
                           acc_bfsqlvl ,acc_difdialvl ,acc_difvmolvl ,acc_difvholvl, &
                           acc_difvsolvl ,acc_difintlvl ,acc_difisolvl ,acc_dzlvl, &
                           acc_salnlvl ,acc_templvl ,acc_uflxlvl ,acc_utflxlvl, &
                           acc_usflxlvl ,acc_umfltdlvl ,acc_umflsmlvl ,acc_utfltdlvl, &
                           acc_utflsmlvl ,acc_utflldlvl ,acc_usfltdlvl ,acc_usflsmlvl, &
                           acc_usflldlvl ,acc_uvellvl ,acc_vflxlvl ,acc_vtflxlvl, &
                           acc_vsflxlvl ,acc_vmfltdlvl ,acc_vmflsmlvl ,acc_vtfltdlvl, &
                           acc_vtflsmlvl ,acc_vtflldlvl ,acc_vsfltdlvl ,acc_vsflsmlvl, &
                           acc_vsflldlvl ,acc_vvellvl ,acc_wflxlvl ,acc_wflx2lvl, &
                           acc_pvlvl ,acc_tkelvl ,acc_gls_psilvl,acc_uflxold, &
                           acc_vflxold ,acc_utillvl, &
                           acc_mmflxl,acc_mmflxd,acc_mmftdl,acc_mmfsml,acc_mmftdd, &
                           acc_mmfsmd,acc_mhflx ,acc_mhftd ,acc_mhfsm ,acc_mhfld, &
                           acc_msflx ,acc_msftd ,acc_msfsm ,acc_msfld ,acc_voltr, &
                           nphy, nacc_phy, ddm
  use mod_forcing,   only: ditflx, disflx, sprfac, &
                           tflxdi, sflxdi, nflxdi, &
                           prfac, eiacc, pracc, &
                           flxco2, flxdms, flxbrf, ustarb, buoyfl, &
                           ustar
  use mod_nctools,   only: ncinqa, ncinqv, ncread, ncgeti, ncgetr,  &
                           ncfopn, ncputi, ncfcls, &
                           ncdimc, ncdims, ncputi, ncputr
  use mod_niw,       only: uml, vml, umlres, vmlres
  use mod_difest,    only: OBLdepth
  use mod_diffusion, only: difiso, Kvisc_m, Kdiff_t, Kdiff_s, &
                           t_ns_nonloc, s_nb_nonloc, &
                           mu_nonloc, mv_nonloc, &
                           umfltd, vmfltd, umflsm, vmflsm, &
                           utfltd, vtfltd, utflsm, vtflsm, &
                           utflld, vtflld, usfltd, vsfltd, &
                           usflsm, vsflsm, usflld, vsflld, &
                           difdia
  use mod_tmsmt,     only: dpold
  use mod_cesm,      only: frzpot, mltpot, swa_da, nsf_da, hmlt_da, &
                           lip_da, sop_da, eva_da, rnf_da, rfi_da, &
                           fmltfz_da, sfl_da, ztx_da, mty_da, ustarw_da, &
                           slp_da, abswnd_da, atmco2_da, atmbrf_da, &
                           ficem_da, l1ci, l2ci
  use mod_ben02,     only: cd_d, ch_d, ce_d, wg2_d, cd_m, ch_m, ce_m, &
                           wg2_m, rhoa, tsi_tda, tml_tda, sml_tda, &
                           alb_tda, fice_tda, ntda, rnfres
  use mod_thdysi,    only: tsrfm, ticem
  use mod_seaice,    only: ficem, hicem, hsnwm, iagem
  use mod_temmin,    only: settemmin
  use mod_tracers,   only: itrtke, itrgls, itriag, trc
  use mod_tke,       only: L_scale
  use mod_ifdefs,    only: use_TRC, use_TKE, use_IDLAGE
  use mod_idlage,    only: idlage_init
  use mod_tracers_update, only: restart_trcrd

  implicit none
  private

  public :: restart_rd

contains

  subroutine restart_rd

    ! --- ------------------------------------------------------------------
    ! --- Read initial conditions from restart file
    ! --- ------------------------------------------------------------------

    ! Local variables
    type(date_type) :: date_rest
    integer :: errstat,dndiff,i,j,n
    character(len = 256) :: rstfnm,fnm
    character(len = 2) :: c2
    real   , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: rkfpla
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: iuu
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: ivv
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: iqq
    logical :: first = .true.
    logical :: fexist,vexist

    ! --- open restart file and adjust time information if needed

    if (nday1+nint(time0) == 0 .and. (.not.resume_flag)) then

      ! --- - open restart file for initial conditions and adjust integration
      ! --- - time corresponding to start date
      rstfnm = icfile
      if (mnproc == 1) inquire(file=rstfnm,exist = fexist)
      call xcbcst(fexist)
      if (fexist) then
        call ncfopn(rstfnm,'r',' ',1,iotype)
        call ncgetr('time',time)
        time0 = time
      else
        if (mnproc == 1) then
          write (lp,*) &
               'restart_rd: could not find restart file for initial conditions!'
        end if
        call xcstop('(restart_rd)')
        stop '(restart_rd)'
      end if

    else if (expcnf == 'cesm'.or.expcnf == 'channel') then

      ! --- - get restart file name from rpointer.ocn
      if (mnproc == 1) &
           inquire(file='rpointer.ocn'//trim(inst_suffix),exist = fexist)
      call xcbcst(fexist)
      if (fexist) then
        if (mnproc == 1) then
          open (unit=nfu,file = 'rpointer.ocn'//trim(inst_suffix))
          read (nfu,'(a)') rstfnm
          close (unit = nfu)
        end if
        call xcbcst(rstfnm)
      else
        if (mnproc == 1) then
          write (lp,*) &
               'restart_rd: could not find file rpointer.ocn'//trim(inst_suffix) &
               //'!'
        end if
        call xcstop('(restart_rd)')
        stop '(restart_rd)'
      end if

      ! --- - open restart file
      if (mnproc == 1) inquire(file=rstfnm,exist = fexist)
      call xcbcst(fexist)
      if (.not.fexist) then
        if (mnproc == 1) then
          write (lp,*) 'restart_rd: could not find restart file!'
        end if
        call xcstop('(restart_rd)')
        stop '(restart_rd)'
      end if
      call ncfopn(rstfnm,'r',' ',1,iotype)

    else

      ! --- - first try file name:
      ! --- -  <experiment name>_restphy_<year>.<month>.<day>_<integration day>.nc
      write (rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)') &
           trim(runid),'_restphy_', &
           date%year,'.',date%month,'.',date%day, &
           '_',nday1,'.nc'

      if (mnproc == 1) inquire(file=rstfnm,exist = fexist)
      call xcbcst(fexist)
      if (fexist) then
        call ncfopn(rstfnm,'r',' ',1,iotype)
        call ncgeti('nday0',date_rest%day)
        call ncgeti('nmonth0',date_rest%month)
        call ncgeti('nyear0',date_rest%year)
        call ncgetr('time0',time0)
        call ncgetr('time',time)
        if (date_rest /= date0) then
          if (mnproc == 1) then
            write (lp,'(2a)') &
                 ' restart_rd: expected identical initial experiment date in', &
                 ' namelist and restart but found:'
            write (lp,'(a,i4.4,2(i2.2))') &
                 ' restart_rd: initial date namelist: ',date0
            write (lp,'(a,i4.4,2(i2.2))') &
                 ' restart_rd: initial date restart:  ',date_rest
          end if
          call xcstop('(restart_rd)')
          stop '(restart_rd)'
        end if

      else

        ! --- --- then try automatic selection of file with consistent
        ! --- --- integration day and date among files named:
        ! --- ---   <experiment name>_restphy_1.nc
        ! --- ---   <experiment name>_restphy_2.nc
        ! --- ---   <experiment name>_restphy_3.nc
        do i = 1,4
          write (rstfnm,'(2a,i1,a)') &
               trim(runid),'_restphy_',i,'.nc'
          inquire(file=rstfnm,exist = fexist)
          if (fexist) then
            call ncfopn(rstfnm,'r',' ',1,iotype)
            call ncgeti('nday0',date_rest%day)
            call ncgeti('nmonth0',date_rest%month)
            call ncgeti('nyear0',date_rest%year)
            call ncgetr('time0',time0)
            call ncgetr('time',time)
            errstat = daynum_diff(calendar,date_rest,date,dndiff)
            if (errstat /= calendar_noerr) then
              if (mnproc == 1) then
                write (lp, '(2a)') ' restart_rd: daynum_diff error: ', &
                     trim(calendar_errstr(errstat))
              end if
              call xcstop('(restart_rd)')
              stop '(restart_rd)'
            end if
            if (nint(time) == nday1.and. &
                 dndiff == nint(time-time0)) exit
          end if
        end do
        if (i > 3) then
          if (mnproc == 1) then
            write (lp,*) &
                 'restart_rd: Could not find proper restart file!'
          end if
          call xcstop('(restart_rd)')
          stop '(restart_rd)'
        end if
      end if

    end if

    if (mnproc == 1) then
      write (lp,'(2a)') ' restart_rd: reading restart file ', &
           trim(rstfnm)
    end if

    ! --- Compute extended uv masks
    if (first) then
      first = .false.
      !$omp parallel do private(i)
      do j = 1,jj
        do i = 1,ii
          if ((ip(i,j)+ip(i-1,j)) >= 1) then
            iuu(i,j) = 1
          else
            iuu(i,j) = 0
          end if
          if ((ip(i,j)+ip(i,j-1)) >= 1) then
            ivv(i,j) = 1
          else
            ivv(i,j) = 0
          end if
          if ((iu(i,j)+iv(i,j)+iu(i,j-1)+iv(i-1,j)) >= 1) then
            iqq(i,j) = 1
          else
            iqq(i,j) = 0
          end if
        end do
      end do
      !$omp end parallel do
    end if

    call ncread('u',u,iuu,1,0.)
    call ncread('v',v,ivv,1,0.)
    call ncread('dp',dp,ip,1,0.)
    call ncread('dpold',dpold,ip,1,0.)
    call ncread('temp',temp,ip,1,0.)
    call ncread('saln',saln,ip,1,0.)
    call ncread('sigma',sigma,ip,1,0.)
    call ncread('sigmar',sigmar,ip,1,0.)
    call ncread('pgfx',pgfx,iuu,1,0.)
    call ncread('pgfy',pgfy,ivv,1,0.)
    call ncread('pb',pb,ip,1,0.)
    call ncread('pb_mn',pb_mn,ip,1,0.)
    call ncread('pb_p',pb_p,ip,1,0.)
    call ncread('pbu',pbu,iuu,1,0.)
    call ncread('pbv',pbv,ivv,1,0.)
    call ncread('pbu_p',pbu_p,iuu,1,0.)
    call ncread('pbv_p',pbv_p,ivv,1,0.)
    call ncread('ub',ub,iuu,1,0.)
    call ncread('vb',vb,ivv,1,0.)
    call ncread('uflx',uflx,iuu,1,0.)
    call ncread('utflx',utflx,iuu,1,0.)
    call ncread('usflx',usflx,iuu,1,0.)

    vexist = ncinqv('umfltd')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('umfltd',umfltd,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''umfltd'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('umflsm')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('umflsm',umflsm,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''umflsm'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('utfltd')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('utfltd',utfltd,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''utfltd'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('utflsm')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('utflsm',utflsm,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''utflsm'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('utflld')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('utflld',utflld,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''utflld'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('usfltd')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('usfltd',usfltd,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''usfltd'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('usflsm')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('usflsm',usflsm,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''usflsm'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('usflld')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('usflld',usflld,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''usflld'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    call ncread('vflx',vflx,ivv,1,0.)
    call ncread('vtflx',vtflx,ivv,1,0.)
    call ncread('vsflx',vsflx,ivv,1,0.)

    vexist = ncinqv('vmfltd')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vmfltd',vmfltd,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vmfltd'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vmflsm')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vmflsm',vmflsm,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vmflsm'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vtfltd')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vtfltd',vtfltd,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vtfltd'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vtflsm')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vtflsm',vtflsm,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vtflsm'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vtflld')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vtflld',vtflld,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vtflld'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vsfltd')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vsfltd',vsfltd,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vsfltd'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vsflsm')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vsflsm',vsflsm,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vsflsm'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vsflld')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vsflld',vsflld,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vsflld'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    call ncread('ubflx',ubflx,iuu,1,0.)
    call ncread('vbflx',vbflx,ivv,1,0.)
    call ncread('ubflx_mn',ubflx_mn,iuu,1,0.)
    call ncread('vbflx_mn',vbflx_mn,ivv,1,0.)
    call ncread('ubflxs',ubflxs,iuu,1,0.)
    call ncread('vbflxs',vbflxs,ivv,1,0.)
    call ncread('ubflxs_p',ubflxs_p,iuu,1,0.)
    call ncread('vbflxs_p',vbflxs_p,ivv,1,0.)
    call ncread('ubcors_p',ubcors_p,iuu,1,0.)
    call ncread('vbcors_p',vbcors_p,ivv,1,0.)
    call ncread('pvtrop',pvtrop,iqq,1,0.)
    call ncread('pgfxm',pgfxm,iuu,1,0.)
    call ncread('pgfym',pgfym,ivv,1,0.)
    call ncread('xixp',xixp,iuu,1,0.)
    call ncread('xixm',xixm,iuu,1,0.)
    call ncread('xiyp',xiyp,ivv,1,0.)
    call ncread('xiym',xiym,ivv,1,0.)
    call ncread('phi',phi(1-nbdy,1-nbdy,kk+1),ip,1,0.)
    call ncread('sealv',sealv,ip,1,0.)
    call ncread('ustar',ustar,ip,1,0.)
    call ncread('kfpla',rkfpla,ip,1,0.)
    call ncread('ficem',ficem,ip,1,0.)

    vexist = ncinqv('uml')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('uml',uml,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''uml'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vml')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vml',vml,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vml'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('umlres')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('umlres',umlres,iuu,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''umlres'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    vexist = ncinqv('vmlres')
    !      call xcbcst(vexist)
    if (vexist) then
      call ncread('vmlres',vmlres,ivv,1,0.)
    else if (mnproc == 1) then
      write (lp,*) &
           'Warning: field ''vmlres'' is not read from restart file and'
      write (lp,*) &
           'will be initialized to zero.'
    end if

    if (vcoord_type_tag == isopyc_bulkml) then
      call ncread('buoyfl',buoyfl,ip,1,0.)
    end if

    if (vcoord_type_tag == cntiso_hybrid) then
      call ncread('dpu',dpu,iu,1,0.)
      call ncread('dpv',dpv,iv,1,0.)
      call ncread('difiso',difiso,ip,1,0.)
      call ncread('OBLdepth',OBLdepth,ip,1,0.)
      call ncread('t_ns_nonloc',t_ns_nonloc,ip,1,0.)
      call ncread('s_nb_nonloc',s_nb_nonloc,ip,1,0.)
      call ncread('mu_nonloc',mu_nonloc,iu,1,0.)
      call ncread('mv_nonloc',mv_nonloc,iv,1,0.)
    end if

    if (sprfac) then
      vexist = ncinqa('prfac')
      !        call xcbcst(vexist)
      if (vexist) then
        call ncgetr('prfac',prfac)
        call xcbcst(prfac)
        call ncread('eiacc',eiacc,ip,1,0.)
        call ncread('pracc',pracc,ip,1,0.)
      else if (mnproc == 1) then
        write (lp,*) &
             'Warning: fields needed for balancing fresh water '// &
             'budget are not read'
        write (lp,*) 'from restart file and will be initialized.'
      end if
    end if

    if (expcnf == 'ben02clim'.or.expcnf == 'ben02syn') then
      call ncread('cd_d',cd_d,ip,1,0.)
      call ncread('ch_d',ch_d,ip,1,0.)
      call ncread('ce_d',ce_d,ip,1,0.)
      call ncread('wg2_d',wg2_d,ip,1,0.)
      call ncread('cd_m',cd_m,ip,1,0.)
      call ncread('ch_m',ch_m,ip,1,0.)
      call ncread('ce_m',ce_m,ip,1,0.)
      call ncread('wg2_m',wg2_m,ip,1,0.)
      call ncread('rhoa',rhoa,ip,1,0.)
      call ncread('tsi_tda',tsi_tda,ip,1,0.)
      call ncread('tml_tda',tml_tda,ip,1,0.)
      call ncread('sml_tda',sml_tda,ip,1,0.)
      call ncread('alb_tda',alb_tda,ip,1,0.)
      call ncread('fice_tda',fice_tda,ip,1,0.)

      call ncgeti('ntda',ntda)
      call xcbcst(ntda)

      call ncread('hicem',hicem,ip,1,0.)
      call ncread('tsrfm',tsrfm,ip,1,0.)
      call ncread('hsnwm',hsnwm,ip,1,0.)
      call ncread('ticem',ticem,ip,1,0.)
      call ncread('iagem',iagem,ip,1,0.)
      call ncread('rnfres',rnfres,ip,1,0.)
    end if

    if (expcnf == 'channel') then
      call ncgeti('ntda',ntda)
      call xcbcst(ntda)
      !        call ncread('rnfres',rnfres,ip,1,0.)
    end if

    if (expcnf == 'cesm') then
      vexist = ncinqv('ustarw_da')
      !        call xcbcst(vexist)
      if (vexist) then
        call ncread('ustarw_da',ustarw_da,ip,1,0.)
        call ncread('ztx_da',ztx_da,ip,1,0.)
        call ncread('mty_da',mty_da,ip,1,0.)
        call ncread('lip_da',lip_da,ip,1,0.)
        call ncread('sop_da',sop_da,ip,1,0.)
        call ncread('eva_da',eva_da,ip,1,0.)
        call ncread('rnf_da',rnf_da,ip,1,0.)
        call ncread('rfi_da',rfi_da,ip,1,0.)
        call ncread('fmltfz_da',fmltfz_da,ip,1,0.)
        call ncread('sfl_da',sfl_da,ip,1,0.)
        call ncread('swa_da',swa_da,ip,1,0.)
        call ncread('nsf_da',nsf_da,ip,1,0.)
        call ncread('hmlt_da',hmlt_da,ip,1,0.)
        call ncread('slp_da',slp_da,ip,1,0.)
        call ncread('ficem_da',ficem_da,ip,1,0.)
        call ncread('abswnd_da',abswnd_da,ip,1,0.)
        call ncread('atmco2_da',atmco2_da,ip,1,0.)
        call ncgeti('l2ci',l2ci)
        call xcbcst(l2ci)
        l1ci = 3-l2ci
      else
        if (mnproc == 1) then
          write (lp,*) &
               'Warning: time levels for interpolation of forcing fields is not'
          write (lp,*) &
               'read from restart file.'
        end if
        l1ci = 1
        l2ci = 1
      end if
      call ncread('frzpot',frzpot,ip,1,0.)
      call ncread('mltpot',mltpot,ip,1,0.)
      vexist = ncinqv('flxco2')
      !        call xcbcst(vexist)
      if (vexist) then
        call ncread('flxco2',flxco2,ip,1,0.)
      else if (mnproc == 1) then
        write (lp,*) &
             'Warning: Air-sea CO2 flux is not read from restart file and'
        write (lp,*) &
             'will be initialized to zero.'
      end if
      vexist = ncinqv('flxdms')
      !        call xcbcst(vexist)
      if (vexist) then
        call ncread('flxdms',flxdms,ip,1,0.)
      else if (mnproc == 1) then
        write (lp,*) &
             'Warning: DMS flux is not read from restart file and'
        write (lp,*) &
             'will be initialized to zero.'
      end if
      vexist = ncinqv('flxbrf')
      !        call xcbcst(vexist)
      if (vexist) then
        call ncread('flxbrf',flxbrf,ip,1,0.)
      else if (mnproc == 1) then
        write (lp,*) &
             'Warning: bromoform flux is not read from restart file and'
        write (lp,*) &
             'will be initialized to zero.'
      end if
    end if

    if (use_TRC) then
      if (use_TKE) then
        vexist = ncinqv('tke')
        !      call xcbcst(vexist)
        if (vexist) then
          call ncread('tke',trc(1-nbdy,1-nbdy,1,itrtke),ip,1,0.)
        else if (mnproc == 1) then
          write (lp,*) &
               'Warning: TKE is not read from restart file and'
          write (lp,*) &
               'will be initialized to tke_min.'
        end if

        vexist = ncinqv('gls_psi')
        !      call xcbcst(vexist)
        if (vexist) then
          call ncread('gls_psi',trc(1-nbdy,1-nbdy,1,itrgls),ip,1,0.)
        else if (mnproc == 1) then
          write (lp,*) &
               'Warning: gls_psi is not read from restart file and'
          write (lp,*) &
               'will be initialized to gls_psi_min.'
        end if

        vexist = ncinqv('L_scale')
        !      call xcbcst(vexist)
        if (vexist) then
          call ncread('L_scale',L_scale,ip,1,0.)
        else if (mnproc == 1) then
          write (lp,*) &
               'Warning: L_scale is not read from restart file and'
          write (lp,*) &
               'will be initialized to Ls_unlmt_min.'
        end if

        vexist = ncinqv('difdia')
        !      call xcbcst(vexist)
        if (vexist) then
          call ncread('difdia',difdia,ip,1,0.)
        else if (mnproc == 1) then
          write (lp,*) &
               'Warning: difdia is not read from restart file and'
          write (lp,*) &
               'will be initialized to zero.'
        end if

        vexist = ncinqv('ustarb')
        !      call xcbcst(vexist)
        if (vexist) then
          call ncread('ustarb',ustarb,ip,1,0.)
        else if (mnproc == 1) then
          write (lp,*) &
               'Warning: ustarb is not read from restart file and'
          write (lp,*) &
               'will be initialized to zero.'
        end if
      end if
      if (use_IDLAGE) then
        vexist = ncinqv('idlage')
        !      call xcbcst(vexist)
        if (vexist) then
          call ncread('idlage',trc(1-nbdy,1-nbdy,1,itriag),ip,1,0.)
        else
          if (mnproc == 1) then
            write (lp,*) &
                 'Warning: Ideal age tracer is not read from restart file and'
            write (lp,*) &
                 'will be initialized to zero.'
          end if
          call idlage_init
        end if
      end if
    end if

    ! --- read accumulated fields
    if (nstep1 > nstep0) then
      do n = 1,nphy
        write(c2,'(i2.2)') n
        if (ncinqa('nacc_phy'//c2)) then
          call ncgeti('nacc_phy'//c2,nacc_phy(n))
        else
          nacc_phy(n) = 0
        end if
        call xcbcst(nacc_phy(n))
        if (nacc_phy(n) > 0) then
          if (acc_ub(n) /= 0) call ncread('ub_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_UB(n)),iuu,1,0.)
          if (acc_vb(n) /= 0) call ncread('vb_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_VB(n)),ivv,1,0.)
          if (acc_ubflxs(n) /= 0) call ncread('ubflxs_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(n)),iuu,1,0.)
          if (acc_vbflxs(n) /= 0) call ncread('vbflxs_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(n)),ivv,1,0.)
          if (acc_ztx(n) /= 0) call ncread('ztx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_ZTX(n)),iuu,1,0.)
          if (acc_mty(n) /= 0) call ncread('mty_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTY(n)),ivv,1,0.)
          if (acc_taux(n) /= 0) call ncread('taux_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_TAUX(n)),iuu,1,0.)
          if (acc_tauy(n) /= 0) call ncread('tauy_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_TAUY(n)),ivv,1,0.)
          if (acc_uice(n) /= 0) call ncread('uice_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_UICE(n)),iuu,1,0.)
          if (acc_vice(n) /= 0) call ncread('vice_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_VICE(n)),ivv,1,0.)
          if (acc_ivolu(n) /= 0) call ncread('ivolu_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_IVOLU(n)),iuu,1,0.)
          if (acc_ivolv(n) /= 0) call ncread('ivolv_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_IVOLV(n)),ivv,1,0.)
          if (acc_psrf(n) /= 0) call ncread('psrf_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_PSRF(n)),ip,1,0.)
          if (acc_pbot(n) /= 0) call ncread('pbot_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_PBOT(n)),ip,1,0.)
          if (acc_sealv(n) /= 0) call ncread('sealv_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SEALV(n)),ip,1,0.)
          if (acc_slvsq(n) /= 0) call ncread('slvsq_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SLVSQ(n)),ip,1,0.)
          if (acc_sss(n) /= 0) call ncread('sss_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SSS(n)),ip,1,0.)
          if (acc_ssssq(n) /= 0) call ncread('ssssq_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SSSSQ(n)),ip,1,0.)
          if (acc_sbot(n) /= 0) call ncread('sbot_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SBOT(n)),ip,1,0.)
          if (acc_sst(n) /= 0) call ncread('sst_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SST(n)),ip,1,0.)
          if (acc_sstsq(n) /= 0) call ncread('sstsq_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SSTSQ(n)),ip,1,0.)
          if (acc_tbot(n) /= 0) call ncread('tbot_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_TBOT(n)),ip,1,0.)
          if (acc_sigmx(n) /= 0) call ncread('sigmx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SIGMX(n)),ip,1,0.)
          if (acc_mld(n) /= 0) call ncread('mld_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MLD(n)),ip,1,0.)
          if (acc_maxmld(n) /= 0) call ncread('maxmld_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MAXMLD(n)),ip,1,0.)
          if (acc_mlts(n) /= 0) call ncread('mlts_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MLTS(n)),ip,1,0.)
          if (acc_mltsmn(n) /= 0) call ncread('mltsmn_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MLTSMN(n)),ip,1,0.)
          if (acc_mltsmx(n) /= 0) call ncread('mltsmx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MLTSMX(n)),ip,1,0.)
          if (acc_mltssq(n) /= 0) call ncread('mltssq_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MLTSSQ(n)),ip,1,0.)
          if (acc_t20d(n) /= 0) call ncread('t20d_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_T20D(n)),ip,1,0.)
          if (acc_alb(n) /= 0) call ncread('alb_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_ALB(n)),ip,1,0.)
          if (acc_swa(n) /= 0) call ncread('swa_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SWA(n)),ip,1,0.)
          if (acc_nsf(n) /= 0) call ncread('nsf_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_NSF(n)),ip,1,0.)
          if (acc_dfl(n) /= 0) call ncread('dfl_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_DFL(n)),ip,1,0.)
          if (acc_lip(n) /= 0) call ncread('lip_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_LIP(n)),ip,1,0.)
          if (acc_sop(n) /= 0) call ncread('sop_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SOP(n)),ip,1,0.)
          if (acc_eva(n) /= 0) call ncread('eva_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_EVA(n)),ip,1,0.)
          if (acc_rnfflx(n) /= 0) call ncread('rnfflx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_RNFFLX(n)),ip,1,0.)
          if (acc_rfiflx(n) /= 0) call ncread('rfiflx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_RFIFLX(n)),ip,1,0.)
          if (acc_sfl(n) /= 0) call ncread('sfl_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SFL(n)),ip,1,0.)
          if (acc_brnflx(n) /= 0) call ncread('brnflx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_BRNFLX(n)),ip,1,0.)
          if (acc_brnpd(n) /= 0) call ncread('brnpd_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_BRNPD(n)),ip,1,0.)
          if (acc_surflx(n) /= 0) call ncread('surflx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SURFLX(n)),ip,1,0.)
          if (acc_surrlx(n) /= 0) call ncread('surrlx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SURRLX(n)),ip,1,0.)
          if (acc_salflx(n) /= 0) call ncread('salflx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SALFLX(n)),ip,1,0.)
          if (acc_salrlx(n) /= 0) call ncread('salrlx_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_SALRLX(n)),ip,1,0.)
          if (acc_abswnd(n) /= 0) call ncread('abswnd_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_ABSWND(n)),ip,1,0.)
          if (acc_ustar(n) /= 0) call ncread('ustar_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_USTAR(n)),ip,1,0.)
          if (acc_ustar3(n) /= 0) call ncread('ustar3_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_USTAR3(n)),ip,1,0.)
          if (acc_idkedt(n) /= 0) call ncread('idkedt_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_IDKEDT(n)),ip,1,0.)
          if (acc_mtkeus(n) /= 0) call ncread('mtkeus_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTKEUS(n)),ip,1,0.)
          if (acc_mtkeni(n) /= 0) call ncread('mtkeni_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTKENI(n)),ip,1,0.)
          if (acc_mtkebf(n) /= 0) call ncread('mtkebf_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTKEBF(n)),ip,1,0.)
          if (acc_mtkers(n) /= 0) call ncread('mtkers_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTKERS(n)),ip,1,0.)
          if (acc_mtkepe(n) /= 0) call ncread('mtkepe_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTKEPE(n)),ip,1,0.)
          if (acc_mtkeke(n) /= 0) call ncread('mtkeke_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_MTKEKE(n)),ip,1,0.)
          if (acc_lamult(n) /= 0) call ncread('lamult_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_LAMULT(n)),ip,1,0.)
          if (acc_lasl(n) /= 0) call ncread('lasl_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_LASL(n)),ip,1,0.)
          if (acc_ustokes(n) /= 0) call ncread('ustokes_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_USTOKES(n)),iuu,1,0.)
          if (acc_vstokes(n) /= 0) call ncread('vstokes_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_VSTOKES(n)),ivv,1,0.)
          if (acc_fmltfz(n) /= 0) call ncread('fmltfz_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_FMLTFZ(n)),ip,1,0.)
          if (acc_hmltfz(n) /= 0) call ncread('hmltfz_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_HMLTFZ(n)),ip,1,0.)
          if (acc_hice(n) /= 0) call ncread('hice_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_HICE(n)),ip,1,0.)
          if (acc_hsnw(n) /= 0) call ncread('hsnw_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_HSNW(n)),ip,1,0.)
          if (acc_fice(n) /= 0) call ncread('fice_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_FICE(n)),ip,1,0.)
          if (acc_tsrf(n) /= 0) call ncread('tsrf_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_TSRF(n)),ip,1,0.)
          if (acc_tice(n) /= 0) call ncread('tice_phy'//c2, &
               phyh2d(1-nbdy,1-nbdy,ACC_TICE(n)),ip,1,0.)
          if (acc_uvel(n) /= 0) call ncread('uvel_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)),iuu,1,0.)
          if (acc_vvel(n) /= 0) call ncread('vvel_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)),ivv,1,0.)
          if (acc_dpu(n) /= 0) call ncread('dpu_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)),iuu,1,0.)
          if (acc_dpv(n) /= 0) call ncread('dpv_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)),ivv,1,0.)
          if (acc_uflx(n) /= 0) call ncread('uflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)),iuu,1,0.)
          if (acc_vflx(n) /= 0) call ncread('vflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)),ivv,1,0.)
          if (acc_utflx(n) /= 0) call ncread('utflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)),iuu,1,0.)
          if (acc_vtflx(n) /= 0) call ncread('vtflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)),ivv,1,0.)
          if (acc_usflx(n) /= 0) call ncread('usflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)),iuu,1,0.)
          if (acc_vsflx(n) /= 0) call ncread('vsflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)),ivv,1,0.)
          if (acc_umfltd(n) /= 0) call ncread('umfltd_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)),iuu,1,0.)
          if (acc_vmfltd(n) /= 0) call ncread('vmfltd_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)),ivv,1,0.)
          if (acc_umflsm(n) /= 0) call ncread('umflsm_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UMFLSM(n)),iuu,1,0.)
          if (acc_vmflsm(n) /= 0) call ncread('vmflsm_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VMFLSM(n)),ivv,1,0.)
          if (acc_utfltd(n) /= 0) call ncread('utfltd_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)),iuu,1,0.)
          if (acc_vtfltd(n) /= 0) call ncread('vtfltd_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)),ivv,1,0.)
          if (acc_utflsm(n) /= 0) call ncread('utflsm_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLSM(n)),iuu,1,0.)
          if (acc_vtflsm(n) /= 0) call ncread('vtflsm_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLSM(n)),ivv,1,0.)
          if (acc_utflld(n) /= 0) call ncread('utflld_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)),iuu,1,0.)
          if (acc_vtflld(n) /= 0) call ncread('vtflld_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)),ivv,1,0.)
          if (acc_usfltd(n) /= 0) call ncread('usfltd_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)),iuu,1,0.)
          if (acc_vsfltd(n) /= 0) call ncread('vsfltd_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)),ivv,1,0.)
          if (acc_usflsm(n) /= 0) call ncread('usflsm_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_USFLSM(n)),iuu,1,0.)
          if (acc_vsflsm(n) /= 0) call ncread('vsflsm_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLSM(n)),ivv,1,0.)
          if (acc_usflld(n) /= 0) call ncread('usflld_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)),iuu,1,0.)
          if (acc_vsflld(n) /= 0) call ncread('vsflld_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)),ivv,1,0.)
          if (acc_saln(n) /= 0) call ncread('saln_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)),ip,1,0.)
          if (acc_temp(n) /= 0) call ncread('temp_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)),ip,1,0.)
          if (acc_dp(n) /= 0) call ncread('dp_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)),ip,1,0.)
          if (acc_dz(n) /= 0) call ncread('dz_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)),ip,1,0.)
          if (acc_bfsq(n) /= 0) call ncread('bfsq_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_BFSQ(n)),ip,1,0.)
          if (acc_difdia(n) /= 0) call ncread('difdia_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)),ip,1,0.)
          if (acc_difvmo(n) /= 0) call ncread('difvmo_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DIFVMO(n)),ip,1,0.)
          if (acc_difvho(n) /= 0) call ncread('difvho_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DIFVHO(n)),ip,1,0.)
          if (acc_difvso(n) /= 0) call ncread('difvso_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DIFVSO(n)),ip,1,0.)
          if (acc_difint(n) /= 0) call ncread('difint_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)),ip,1,0.)
          if (acc_difiso(n) /= 0) call ncread('difiso_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)),ip,1,0.)
          if (acc_wflx(n) /= 0) call ncread('wflx_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)),ip,1,0.)
          if (acc_wflx2(n) /= 0) call ncread('wflx2_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)),ip,1,0.)
          if (acc_avdsg(n) /= 0) call ncread('avdsg_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_AVDSG(n)),ip,1,0.)
          if (acc_dpvor(n) /= 0) call ncread('dpvor_phy'//c2, &
               phylyr(1-nbdy,1-nbdy,1,ACC_DPVOR(n)),ip,1,0.)
          if (use_TRC .and. use_TKE) then
            if (acc_tke(n) /= 0) call ncread('tke_phy'//c2, &
                 phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)),ip,1,0.)
            if (acc_gls_psi(n) /= 0) call ncread('gls_psi_phy'//c2, &
                 phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)),ip,1,0.)
          end if
          if (acc_uvellvl(n) /= 0) call ncread('uvellvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UVELLVL(n)),iuu,1,0.)
          if (acc_vvellvl(n) /= 0) call ncread('vvellvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VVELLVL(n)),ivv,1,0.)
          if (acc_uflxlvl(n) /= 0) then
            call ncread('uflxlvl_phy'//c2, &
                 phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)), &
                 iuu,1,0.)
            call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)), &
                 1,ddm, 1,1, halo_uv)
          end if
          if (acc_vflxlvl(n) /= 0) then
            call ncread('vflxlvl_phy'//c2, &
                 phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)), &
                 ivv,1,0.)
            call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)), &
                 1,ddm, 1,1, halo_vv)
          end if
          if (acc_utflxlvl(n) /= 0) call ncread('utflxlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLXLVL(n)),iuu,1,0.)
          if (acc_vtflxlvl(n) /= 0) call ncread('vtflxlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLXLVL(n)),ivv,1,0.)
          if (acc_usflxlvl(n) /= 0) call ncread('usflxlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_USFLXLVL(n)),iuu,1,0.)
          if (acc_vsflxlvl(n) /= 0) call ncread('vsflxlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLXLVL(n)),ivv,1,0.)
          if (acc_umfltdlvl(n) /= 0) call ncread('umfltdlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UMFLTDLVL(n)),iuu,1,0.)
          if (acc_vmfltdlvl(n) /= 0) call ncread('vmfltdlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VMFLTDLVL(n)),ivv,1,0.)
          if (acc_umflsmlvl(n) /= 0) call ncread('umflsmlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UMFLSMLVL(n)),iuu,1,0.)
          if (acc_vmflsmlvl(n) /= 0) call ncread('vmflsmlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VMFLSMLVL(n)),ivv,1,0.)
          if (acc_utfltdlvl(n) /= 0) call ncread('utfltdlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLTDLVL(n)),iuu,1,0.)
          if (acc_vtfltdlvl(n) /= 0) call ncread('vtfltdlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLTDLVL(n)),ivv,1,0.)
          if (acc_utflsmlvl(n) /= 0) call ncread('utflsmlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLSMLVL(n)),iuu,1,0.)
          if (acc_vtflsmlvl(n) /= 0) call ncread('vtflsmlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLSMLVL(n)),ivv,1,0.)
          if (acc_utflldlvl(n) /= 0) call ncread('utflldlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLLDLVL(n)),iuu,1,0.)
          if (acc_vtflldlvl(n) /= 0) call ncread('vtflldlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLLDLVL(n)),ivv,1,0.)
          if (acc_usfltdlvl(n) /= 0) call ncread('usfltdlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_USFLTDLVL(n)),iuu,1,0.)
          if (acc_vsfltdlvl(n) /= 0) call ncread('vsfltdlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLTDLVL(n)),ivv,1,0.)
          if (acc_usflsmlvl(n) /= 0) call ncread('usflsmlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_USFLSMLVL(n)),iuu,1,0.)
          if (acc_vsflsmlvl(n) /= 0) call ncread('vsflsmlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLSMLVL(n)),ivv,1,0.)
          if (acc_usflldlvl(n) /= 0) call ncread('usflldlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_USFLLDLVL(n)),iuu,1,0.)
          if (acc_vsflldlvl(n) /= 0) call ncread('vsflldlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLLDLVL(n)),ivv,1,0.)
          if (acc_salnlvl(n) /= 0) call ncread('salnlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_SALNLVL(n)),ip,1,0.)
          if (acc_templvl(n) /= 0) call ncread('templvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_TEMPLVL(n)),ip,1,0.)
          if (acc_dzlvl(n) /= 0) call ncread('dzlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DZLVL(n)),ip,1,0.)
          if (acc_bfsqlvl(n) /= 0) call ncread('bfsqlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_BFSQLVL(n)),ip,1,0.)
          if (acc_difdialvl(n) /= 0) call ncread('difdialvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DIFDIALVL(n)),ip,1,0.)
          if (acc_difvmolvl(n) /= 0) call ncread('difvmolvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DIFVMOLVL(n)),ip,1,0.)
          if (acc_difvholvl(n) /= 0) call ncread('difvholvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DIFVHOLVL(n)),ip,1,0.)
          if (acc_difvsolvl(n) /= 0) call ncread('difvsolvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DIFVSOLVL(n)),ip,1,0.)
          if (acc_difintlvl(n) /= 0) call ncread('difintlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DIFINTLVL(n)),ip,1,0.)
          if (acc_difisolvl(n) /= 0) call ncread('difisolvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_DIFISOLVL(n)),ip,1,0.)
          if (acc_wflxlvl(n) /= 0) call ncread('wflxlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)),ip,1,0.)
          if (acc_wflx2lvl(n) /= 0) call ncread('wflx2lvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)),ip,1,0.)
          if (acc_pvlvl(n) /= 0) call ncread('pvlvl_phy'//c2, &
               phylvl(1-nbdy,1-nbdy,1,ACC_PVLVL(n)),ip,1,0.)
          if (use_TRC .and. use_TKE) then
            if (acc_tkelvl(n) /= 0) call ncread('tkelvl_phy'//c2, &
                 phylvl(1-nbdy,1-nbdy,1,ACC_TKELVL(n)),ip,1,0.)
            if (acc_gls_psilvl(n) /= 0) &
                 call ncread('gls_psilvl_phy'//c2, &
                 phylvl(1-nbdy,1-nbdy,1,ACC_GLS_PSILVL(n)), &
                 ip,1,0.)
          end if
        end if
      end do
    end if

    call ncfcls

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        if (ip(i,j) == 1) then
          kfpla(i,j,1) = nint(rkfpla(i,j,1))
          kfpla(i,j,2) = nint(rkfpla(i,j,2))
        end if
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- set minimum physical temperature for each isopycnic layer
    ! --- ------------------------------------------------------------------

    call settemmin

    if (use_TRC) then
      if (.not.resume_flag) call restart_trcrd(rstfnm)
    end if

    if (ditflx) then

      ! --- - read diag. heat flux restart file if available
      if (expcnf == 'cesm') then
        fnm = trim(runid)//'.blom.rtflx.'//rstfnm(len_trim(runid)+10:)
      else
        fnm = trim(runid)//'_resttflx_'//rstfnm(len_trim(runid)+10:)
      end if
      if(mnproc  == 1) inquire(file=fnm,exist = fexist)
      call xcbcst(fexist)
      if (fexist) then
        call ncfopn(fnm,'r',' ',1,iotype)
        if (mnproc == 1) &
             write (lp,'(a,a)') ' reading diag. heat flux restart file ', &
             trim(fnm)
        call ncgetr('time',time)
        if (nint(time) /= nday1.and.mnproc == 1) then
          write (lp,'(a,i6.6,a)') &
               ' Integration day ',nint(time), &
               ' in diag. heat flux restart file differs from'
          write (lp,'(a,i6.6,a)') &
               ' start day ',nday1,' in limits file!'
          call ncfcls
          call xcstop('(restart_rd)')
          stop '(restart_rd)'
        end if
        call ncread('tflxdi',tflxdi,ip,1,0.)
        call ncgeti('nflxdi',nflxdi)
      else
        if (mnproc == 1) write (lp,*) &
             'Warning! No diag. heat flux restart file found'
      end if

    end if

    if (disflx) then

      ! --- - read diag. salt flux restart file if available
      if (expcnf == 'cesm') then
        fnm = trim(runid)//'.blom.rsflx.'//rstfnm(len_trim(runid)+10:)
      else
        fnm = trim(runid)//'_restsflx_'//rstfnm(len_trim(runid)+10:)
      end if
      if (mnproc == 1) inquire(file=fnm,exist = fexist)
      call xcbcst(fexist)
      if (fexist) then
        call ncfopn(fnm,'r',' ',1,iotype)
        if (mnproc == 1) &
             write (lp,'(a,a)') ' reading diag. salt flux restart file ', &
             trim(fnm)
        call ncgetr('time',time)
        if (nint(time) /= nday1.and.mnproc == 1) then
          write (lp,'(a,i6.6,a)') &
               ' Integration day ',nint(time), &
               ' in diag. salt flux restart file differs from'
          write (lp,'(a,i6.6,a)') &
               ' start day ',nday1,' in limits file!'
          call ncfcls
          call xcstop('(restart_rd)')
          stop '(restart_rd)'
        end if
        call ncread('sflxdi',sflxdi,ip,1,0.)
        call ncgeti('nflxdi',nflxdi)
        call ncfcls
      else
        if (mnproc == 1) write (lp,*) &
             'Warning! No diag. salt flux restart file found'
      end if

    end if

  end subroutine restart_rd

end module mod_restart_rd
