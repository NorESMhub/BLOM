! ------------------------------------------------------------------------------
! Copyright (C) 2006-2023 Mats Bentsen, Mehmet Ilicak, Alok Kumar Gupta,
!                         Ingo Bethke, Jerry Tjiputra, Ping-Gin Chiu,
!                         Aleksi Nummelin, Jörg Schwinger

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

module mod_restart_wt

  use dimensions,     only: idm, jdm, kdm, itdm, jtdm
  use mod_config,     only: expcnf, runid, inst_suffix
  use mod_time,       only: date0, date, nstep, nstep_in_day, nday_of_year, &
                            time, time0
  use mod_xc,         only: lp, ii, jj, kk, iu, iv, ip, &
                            ilu, isv, ifv, ilv, isq, ifq, ilq, &
                            nfu, nbdy, mnproc
  use mod_vcoord,     only: vcoord_type_tag, isopyc_bulkml, &
                            cntiso_hybrid, sigmar
  use mod_state,      only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                            uflx, vflx, utflx, vtflx, usflx, vsflx, &
                            phi, ubflxs, vbflxs, &
                            ub, vb, pb, pbu, pbv, ubflxs_p, vbflxs_p, &
                            pb_p, pbu_p, pbv_p, ubcors_p, vbcors_p, &
                            sealv, kfpla
  use mod_pgforc,     only: pgfx, pgfy, pgfxm, pgfym, &
                            xixp, xixm, xiyp, xiym
  use mod_barotp,     only: ubflx, vbflx, pb_mn, ubflx_mn, vbflx_mn, &
                            pvtrop
  use mod_dia,        only: rstfrq, iotype, rstfmt, rstcmp, rstmon, &
                            ddm, depthslev, phyh2d, phylyr, phylvl, &
                            nphy, nacc_phy, &
                            acc_abswnd ,acc_alb ,acc_brnflx ,acc_brnpd ,acc_dfl , &
                            acc_eva ,acc_fice ,acc_fmltfz ,acc_hice ,acc_hmltfz , &
                            acc_hsnw ,acc_iage ,acc_idkedt ,acc_lamult ,acc_lasl , &
                            acc_lip ,acc_maxmld ,acc_mld ,acc_mlts ,acc_mltsmn , &
                            acc_mltsmx ,acc_mltssq ,acc_mtkeus ,acc_mtkeni ,acc_mtkebf , &
                            acc_mtkers ,acc_mtkepe ,acc_mtkeke ,acc_mty ,acc_nsf , &
                            acc_pbot ,acc_psrf ,acc_rfiflx ,acc_rnfflx ,acc_salflx , &
                            acc_salrlx ,acc_sbot ,acc_sealv ,acc_slvsq ,acc_sfl , &
                            acc_sop ,acc_sigmx ,acc_sss ,acc_ssssq ,acc_sst , &
                            acc_sstsq ,acc_surflx ,acc_surrlx ,acc_swa ,acc_t20d , &
                            acc_taux ,acc_tauy ,acc_tbot ,acc_tice ,acc_tsrf , &
                            acc_ub ,acc_ubflxs ,acc_uice ,acc_ustar ,acc_ustar3 , &
                            acc_ustokes,acc_vb ,acc_vbflxs ,acc_vice ,acc_vstokes, &
                            acc_ztx ,acc_ivolu ,acc_ivolv ,acc_utilh2d, &
                            acc_bfsq ,acc_difdia ,acc_difvmo ,acc_difvho ,acc_difvso , &
                            acc_difint ,acc_difiso ,acc_dp ,acc_dpu ,acc_dpv , &
                            acc_dz ,acc_saln ,acc_temp ,acc_uflx ,acc_utflx , &
                            acc_usflx ,acc_umfltd ,acc_umflsm ,acc_utfltd ,acc_utflsm , &
                            acc_utflld ,acc_usfltd ,acc_usflsm ,acc_usflld ,acc_uvel , &
                            acc_vflx ,acc_vtflx ,acc_vsflx ,acc_vmfltd ,acc_vmflsm , &
                            acc_vtfltd ,acc_vtflsm ,acc_vtflld ,acc_vsfltd ,acc_vsflsm , &
                            acc_vsflld ,acc_vvel ,acc_wflx ,acc_wflx2 ,acc_avdsg , &
                            acc_dpvor ,acc_tke ,acc_gls_psi,acc_utillyr, &
                            acc_bfsqlvl ,acc_difdialvl ,acc_difvmolvl ,acc_difvholvl , &
                            acc_difvsolvl ,acc_difintlvl ,acc_difisolvl ,acc_dzlvl , &
                            acc_salnlvl ,acc_templvl ,acc_uflxlvl ,acc_utflxlvl , &
                            acc_usflxlvl ,acc_umfltdlvl ,acc_umflsmlvl ,acc_utfltdlvl , &
                            acc_utflsmlvl ,acc_utflldlvl ,acc_usfltdlvl ,acc_usflsmlvl , &
                            acc_usflldlvl ,acc_uvellvl ,acc_vflxlvl ,acc_vtflxlvl , &
                            acc_vsflxlvl ,acc_vmfltdlvl ,acc_vmflsmlvl ,acc_vtfltdlvl , &
                            acc_vtflsmlvl ,acc_vtflldlvl ,acc_vsfltdlvl ,acc_vsflsmlvl , &
                            acc_vsflldlvl ,acc_vvellvl ,acc_wflxlvl ,acc_wflx2lvl , &
                            acc_pvlvl ,acc_tkelvl ,acc_gls_psilvl,acc_uflxold , &
                            acc_vflxold ,acc_utillvl , &
                            acc_mmflxl,acc_mmflxd,acc_mmftdl,acc_mmfsml,acc_mmftdd, &
                            acc_mmfsmd,acc_mhflx ,acc_mhftd ,acc_mhfsm ,acc_mhfld , &
                            acc_msflx ,acc_msftd ,acc_msfsm ,acc_msfld ,acc_voltr
  use mod_forcing,    only: ditflx, disflx, sprfac, &
                            tflxdi, sflxdi, nflxdi, &
                            prfac, eiacc, pracc, &
                            flxco2, flxdms, ustarb, buoyfl,flxbrf, ustar
  use mod_niw,        only: uml, vml, umlres, vmlres
  use mod_difest,     only: OBLdepth
  use mod_diffusion,  only: difiso, Kvisc_m, Kdiff_t, Kdiff_s, &
                            t_ns_nonloc, s_nb_nonloc, &
                            mu_nonloc, mv_nonloc, &
                            umfltd, vmfltd, umflsm, vmflsm, &
                            utfltd, vtfltd, utflsm, vtflsm, &
                            utflld, vtflld, usfltd, vsfltd, &
                            usflsm, vsflsm, usflld, vsflld, &
                            difdia
  use mod_nctools,    only: ncfopn, ncdimc, ncdims, ncputr, ncputi, &
                            ncedef, ncfcls, ndouble, ncwrtr, ncdefvar, nccomp
  use mod_cesm,       only: frzpot, mltpot, swa_da, nsf_da, hmlt_da, &
                            lip_da, sop_da, eva_da, rnf_da, rfi_da, &
                            fmltfz_da, sfl_da, ztx_da, mty_da, ustarw_da, &
                            slp_da, abswnd_da, atmco2_da, atmbrf_da, &
                            ficem_da, l2ci
  use mod_ben02,      only: cd_d, ch_d, ce_d, wg2_d, cd_m, ch_m, ce_m, &
                            wg2_m, rhoa, tsi_tda, tml_tda, sml_tda, &
                            alb_tda, fice_tda, ntda, rnfres
  use mod_thdysi,     only: tsrfm, ticem
  use mod_seaice,     only: ficem, hicem, hsnwm, iagem
  use mod_tmsmt,      only: dpold
  use mod_tracers,    only: itrtke, itrgls, itriag, trc
  use mod_tke,        only: L_scale
#ifdef HAMOCC
  use mo_control_bgc, only: use_BROMO
#endif
  use mod_ifdefs,     only: use_TRC, use_TKE, use_IDLAGE
  use mod_tracers_update, only: restart_trcwt

  implicit none
  private

  public :: restart_wt
  public :: defvar_restart

  private :: wrtrst
  private :: defvarrst

contains

  subroutine restart_wt

    ! --- ------------------------------------------------------------------
    ! --- Write model state to restart files
    ! --- ------------------------------------------------------------------

    ! Local variables
    integer :: i,j,n
    character(len = 256), dimension(4) :: rstdate_str
    character(len = 256) :: rstfnm,fnm
    character(len = 2) :: c2
    character(len = 5) :: c5p,c5u,c5v,c5q
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: rkfpla
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: iuu
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: ivv
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: iqq
    logical :: first = .true.

    ! --- formulate restart name and open restart file

    if (expcnf == 'cesm') then
      write (rstfnm,'(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') &
           trim(runid),'.blom',trim(inst_suffix),'.r.', &
           date%year,'-',date%month,'-',date%day,'-', &
           86400*mod(nstep,nstep_in_day)/nstep_in_day,'.nc'
    else
      if (nday_of_year-max(1,nint(rstfrq/nstep_in_day)) <= 0) then
        write(rstfnm,'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)') &
             trim(runid),'_restphy_', &
             date%year,'.',date%month,'.',date%day, &
             '_',nint(time),'.nc'
      else
        if (rstmon) then
          write (rstfnm,'(2a,i1,a)') &
               trim(runid),'_restphy_', &
               mod(date%month+10,3)+1,'.nc'
        else
          write (rstfnm,'(2a,i1,a)') &
               trim(runid),'_restphy_', &
               mod(nint(min(nstep/rstfrq,time))-1,3)+1,'.nc'
        end if
        open (unit=nfu,file = 'rstdate.txt')
        i = 1
300     read (nfu,'(a)',end = 301) rstdate_str(i)
        i = i+1
        goto 300
301     close (unit = nfu)
        write(rstdate_str(i),'(2a,i4.4,a,i2.2,a,i2.2,a,i6.6)') &
             rstfnm(1:len_trim(runid)+13), &
             ': date ',date%year,'.',date%month,'.',date%day, &
             ', integration day ',nint(time)
        if (mnproc == 1) then
          if (i == 1) then
            open (unit=nfu,file = 'rstdate.txt')
            write (nfu,'(a)') rstdate_str(1)(1:len_trim(runid)+54)
            close (unit = nfu)
          else if (rstdate_str(max(1,i-2)) /= rstdate_str(i).and. &
               rstdate_str(i-1       ) /= rstdate_str(i)) then
            open (unit=nfu,file = 'rstdate.txt')
            do j = max(1,i-2),i
              write (nfu,'(a)') rstdate_str(j)(1:len_trim(runid)+54)
            end do
            close (unit = nfu)
          end if
        end if
      end if
      if (mnproc == 1) then
        write (lp,'(a,a)') ' saving restart file ',trim(rstfnm)
      end if
    end if

    if (rstfmt == 1) then
      call ncfopn(rstfnm,'w','6',1,iotype)
    else if (rstfmt == 2) then
      call ncfopn(rstfnm,'w','h',1,iotype)
    else
      call ncfopn(rstfnm,'w','c',1,iotype)
    end if
    call ncputi('nday0',date0%day)
    call ncputi('nmonth0',date0%month)
    call ncputi('nyear0',date0%year)
    call ncputr('time0',time0)
    call ncputr('time',time)

    ! --- define spatial and time dimensions
    if (first) then
      first = .false.
      !$OMP PARALLEL DO PRIVATE(i)
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
      !$OMP END PARALLEL DO
    end if
    if (rstcmp == 1) then
      call ncdimc('pcomp',ip,0)
      call ncdimc('qcomp',iqq,0)
      call ncdimc('ucomp',iuu,0)
      call ncdimc('vcomp',ivv,0)
      c5p = 'pcomp'
      c5u = 'ucomp'
      c5v = 'vcomp'
      c5q = 'qcomp'
    else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      c5p = 'x y'
      c5u = 'x y'
      c5v = 'x y'
      c5q = 'x y'
    end if
    call ncdims('k2',2)
    call ncdims('k3',3)
    call ncdims('k4',4)
    call ncdims('kk',kk)
    call ncdims('kkp1',kk+1)
    call ncdims('kk2',2*kk)
    call ncdims('plev',ddm)
    call ncputr('plev',depthslev)
    call ncdims('time',1)

    ! --- output model fields to restart file
    if (sprfac) then
      call ncputr('prfac',prfac)
    end if
    if (expcnf == 'ben02clim'.or.expcnf == 'ben02syn'.or. &
         expcnf == 'channel') then
      call ncputi('ntda',ntda)
    end if
    do n = 1,nphy
      write(c2,'(i2.2)') n
      call ncputi('nacc_phy'//c2,nacc_phy(n))
    end do
    if (expcnf == 'cesm') then
      call ncputi('l2ci',l2ci)
    end if

    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1,jj
      do i = 1,ii
        if (ip(i,j) == 1) then
          rkfpla(i,j,1) = real(kfpla(i,j,1))
          rkfpla(i,j,2) = real(kfpla(i,j,2))
        else
          rkfpla(i,j,1) = 0.
          rkfpla(i,j,2) = 0.
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    call defvar_restart(c5p,c5u,c5v,c5q)
    call wrtrst('u',trim(c5u)//' kk2 time',u,iuu)
    call wrtrst('v',trim(c5v)//' kk2 time',v,ivv)
    call wrtrst('dp',trim(c5p)//' kk2 time',dp,ip)
    call wrtrst('dpold',trim(c5p)//' kk2 time',dpold,ip)
    call wrtrst('temp',trim(c5p)//' kk2 time',temp,ip)
    call wrtrst('saln',trim(c5p)//' kk2 time',saln,ip)
    call wrtrst('sigma',trim(c5p)//' kk2 time',sigma,ip)
    call wrtrst('sigmar',trim(c5p)//' kk time',sigmar,ip)
    call wrtrst('pgfx',trim(c5u)//' kk2 time',pgfx,iuu)
    call wrtrst('pgfy',trim(c5v)//' kk2 time',pgfy,ivv)
    call wrtrst('pb',trim(c5p)//' k2 time',pb,ip)
    call wrtrst('pb_mn',trim(c5p)//' k2 time',pb_mn,ip)
    call wrtrst('pb_p',trim(c5p)//' time',pb_p,ip)
    call wrtrst('pbu',trim(c5u)//' k2 time',pbu,iuu)
    call wrtrst('pbv',trim(c5v)//' k2 time',pbv,ivv)
    call wrtrst('pbu_p',trim(c5u)//' time',pbu_p,iuu)
    call wrtrst('pbv_p',trim(c5v)//' time',pbv_p,ivv)
    call wrtrst('ub',trim(c5u)//' k2 time',ub,iuu)
    call wrtrst('vb',trim(c5v)//' k2 time',vb,ivv)
    call wrtrst('uflx',trim(c5u)//' kk2 time',uflx,iuu)
    call wrtrst('utflx',trim(c5u)//' kk2 time',utflx,iuu)
    call wrtrst('usflx',trim(c5u)//' kk2 time',usflx,iuu)
    call wrtrst('umfltd',trim(c5u)//' kk2 time',umfltd,iuu)
    call wrtrst('utfltd',trim(c5u)//' kk2 time',utfltd,iuu)
    call wrtrst('utflld',trim(c5u)//' kk2 time',utflld,iuu)
    call wrtrst('usfltd',trim(c5u)//' kk2 time',usfltd,iuu)
    call wrtrst('usflld',trim(c5u)//' kk2 time',usflld,iuu)
    call wrtrst('vflx',trim(c5v)//' kk2 time',vflx,ivv)
    call wrtrst('vtflx',trim(c5v)//' kk2 time',vtflx,ivv)
    call wrtrst('vsflx',trim(c5v)//' kk2 time',vsflx,ivv)
    call wrtrst('vmfltd',trim(c5v)//' kk2 time',vmfltd,ivv)
    call wrtrst('vtfltd',trim(c5v)//' kk2 time',vtfltd,ivv)
    call wrtrst('vtflld',trim(c5v)//' kk2 time',vtflld,ivv)
    call wrtrst('vsfltd',trim(c5v)//' kk2 time',vsfltd,ivv)
    call wrtrst('vsflld',trim(c5v)//' kk2 time',vsflld,ivv)
    call wrtrst('ubflx',trim(c5u)//' k2 time',ubflx,iuu)
    call wrtrst('vbflx',trim(c5v)//' k2 time',vbflx,ivv)
    call wrtrst('ubflx_mn',trim(c5u)//' k2 time',ubflx_mn,iuu)
    call wrtrst('vbflx_mn',trim(c5v)//' k2 time',vbflx_mn,ivv)
    call wrtrst('ubflxs',trim(c5u)//' k3 time',ubflxs,iuu)
    call wrtrst('vbflxs',trim(c5v)//' k3 time',vbflxs,ivv)
    call wrtrst('ubflxs_p',trim(c5u)//' k2 time',ubflxs_p,iuu)
    call wrtrst('vbflxs_p',trim(c5v)//' k2 time',vbflxs_p,ivv)
    call wrtrst('ubcors_p',trim(c5u)//' time',ubcors_p,iuu)
    call wrtrst('vbcors_p',trim(c5v)//' time',vbcors_p,ivv)
    call wrtrst('pvtrop',trim(c5q)//' k2 time',pvtrop,iqq)
    call wrtrst('pgfxm',trim(c5u)//' k2 time',pgfxm,iuu)
    call wrtrst('pgfym',trim(c5v)//' k2 time',pgfym,ivv)
    call wrtrst('xixp',trim(c5u)//' k2 time',xixp,iuu)
    call wrtrst('xixm',trim(c5u)//' k2 time',xixm,iuu)
    call wrtrst('xiyp',trim(c5v)//' k2 time',xiyp,ivv)
    call wrtrst('xiym',trim(c5v)//' k2 time',xiym,ivv)
    call wrtrst('phi',trim(c5p)//' time',phi(1-nbdy,1-nbdy,kk+1),ip)
    call wrtrst('sealv',trim(c5p)//' time',sealv,ip)
    call wrtrst('ustar',trim(c5p)//' time',ustar,ip)
    call wrtrst('kfpla',trim(c5p)//' k2 time',rkfpla,ip)
    call wrtrst('ficem',trim(c5p)//' time',ficem,ip)
    call wrtrst('uml',trim(c5u)//' k4 time',uml,iuu)
    call wrtrst('vml',trim(c5v)//' k4 time',vml,ivv)
    call wrtrst('umlres',trim(c5u)//' k2 time',umlres,iuu)
    call wrtrst('vmlres',trim(c5v)//' k2 time',vmlres,ivv)

    if (vcoord_type_tag == isopyc_bulkml) then
      call wrtrst('buoyfl',trim(c5p)//' time',buoyfl,ip)
    end if

    if (vcoord_type_tag == cntiso_hybrid) then
      call wrtrst('dpu',trim(c5p)//' kk2 time',dpu,iu)
      call wrtrst('dpv',trim(c5p)//' kk2 time',dpv,iv)
      call wrtrst('difiso',trim(c5p)//' kk time',difiso,ip)
      call wrtrst('OBLdepth',trim(c5p)//' time',OBLdepth,ip)
      call wrtrst('t_ns_nonloc',trim(c5p)//' kkp1 time', &
           t_ns_nonloc,ip)
      call wrtrst('s_nb_nonloc',trim(c5p)//' kkp1 time', &
           s_nb_nonloc,ip)
      call wrtrst('mu_nonloc',trim(c5u)//' kkp1 time',mu_nonloc,iu)
      call wrtrst('mv_nonloc',trim(c5v)//' kkp1 time',mv_nonloc,iv)
      call wrtrst('umflsm',trim(c5u)//' kk2 time',umflsm,iuu)
      call wrtrst('utflsm',trim(c5u)//' kk2 time',utflsm,iuu)
      call wrtrst('usflsm',trim(c5u)//' kk2 time',usflsm,iuu)
      call wrtrst('vmflsm',trim(c5v)//' kk2 time',vmflsm,ivv)
      call wrtrst('vtflsm',trim(c5v)//' kk2 time',vtflsm,ivv)
      call wrtrst('vsflsm',trim(c5v)//' kk2 time',vsflsm,ivv)
    end if

    if (sprfac) then
      call wrtrst('eiacc',trim(c5p)//' time',eiacc,ip)
      call wrtrst('pracc',trim(c5p)//' time',pracc,ip)
    end if

    if (expcnf == 'ben02clim'.or.expcnf == 'ben02syn') then
      call wrtrst('cd_d',trim(c5p)//' time',cd_d,ip)
      call wrtrst('ch_d',trim(c5p)//' time',ch_d,ip)
      call wrtrst('ce_d',trim(c5p)//' time',ce_d,ip)
      call wrtrst('wg2_d',trim(c5p)//' time',wg2_d,ip)
      call wrtrst('cd_m',trim(c5p)//' time',cd_m,ip)
      call wrtrst('ch_m',trim(c5p)//' time',ch_m,ip)
      call wrtrst('ce_m',trim(c5p)//' time',ce_m,ip)
      call wrtrst('wg2_m',trim(c5p)//' time',wg2_m,ip)
      call wrtrst('rhoa',trim(c5p)//' time',rhoa,ip)
      call wrtrst('tsi_tda',trim(c5p)//' time',tsi_tda,ip)
      call wrtrst('tml_tda',trim(c5p)//' time',tml_tda,ip)
      call wrtrst('sml_tda',trim(c5p)//' time',sml_tda,ip)
      call wrtrst('alb_tda',trim(c5p)//' time',alb_tda,ip)
      call wrtrst('fice_tda',trim(c5p)//' time',fice_tda,ip)
      call wrtrst('hicem',trim(c5p)//' time',hicem,ip)
      call wrtrst('tsrfm',trim(c5p)//' time',tsrfm,ip)
      call wrtrst('hsnwm',trim(c5p)//' time',hsnwm,ip)
      call wrtrst('ticem',trim(c5p)//' time',ticem,ip)
      call wrtrst('iagem',trim(c5p)//' time',iagem,ip)
      call wrtrst('rnfres',trim(c5p)//' time',rnfres,ip)
    end if

    !      if (expcnf.eq.'channel') then
    !        call wrtrst('rnfres',trim(c5p)//' time',rnfres,ip)
    !      endif

    if (expcnf == 'cesm') then
      call wrtrst('ustarw_da',trim(c5p)//' k2 time',ustarw_da,ip)
      call wrtrst('ztx_da',trim(c5p)//' k2 time',ztx_da,ip)
      call wrtrst('mty_da',trim(c5p)//' k2 time',mty_da,ip)
      call wrtrst('lip_da',trim(c5p)//' k2 time',lip_da,ip)
      call wrtrst('sop_da',trim(c5p)//' k2 time',sop_da,ip)
      call wrtrst('eva_da',trim(c5p)//' k2 time',eva_da,ip)
      call wrtrst('rnf_da',trim(c5p)//' k2 time',rnf_da,ip)
      call wrtrst('rfi_da',trim(c5p)//' k2 time',rfi_da,ip)
      call wrtrst('fmltfz_da',trim(c5p)//' k2 time',fmltfz_da,ip)
      call wrtrst('sfl_da',trim(c5p)//' k2 time',sfl_da,ip)
      call wrtrst('swa_da',trim(c5p)//' k2 time',swa_da,ip)
      call wrtrst('nsf_da',trim(c5p)//' k2 time',nsf_da,ip)
      call wrtrst('hmlt_da',trim(c5p)//' k2 time',hmlt_da,ip)
      call wrtrst('slp_da',trim(c5p)//' k2 time',slp_da,ip)
      call wrtrst('ficem_da',trim(c5p)//' k2 time',ficem_da,ip)
      call wrtrst('abswnd_da',trim(c5p)//' k2 time',abswnd_da,ip)
      call wrtrst('atmco2_da',trim(c5p)//' k2 time',atmco2_da,ip)
      call wrtrst('atmbrf_da',trim(c5p)//' k2 time',atmbrf_da,ip)  ! not read in restart_rd, necesarry?
      call wrtrst('frzpot',trim(c5p)//' time',frzpot,ip)
      call wrtrst('mltpot',trim(c5p)//' time',mltpot,ip)
      call wrtrst('flxco2',trim(c5p)//' time',flxco2,ip)
      call wrtrst('flxdms',trim(c5p)//' time',flxdms,ip)
#ifdef HAMOCC
      if (use_BROMO) then
        call wrtrst('flxbrf',trim(c5p)//' time',flxbrf,ip)
      end if
#endif
    end if

    if (use_TRC) then
      if (use_TKE) then
        call wrtrst('tke',trim(c5p)//' kk2 time', &
             trc(1-nbdy,1-nbdy,1,itrtke),ip)
        call wrtrst('gls_psi',trim(c5p)//' kk2 time', &
             trc(1-nbdy,1-nbdy,1,itrgls),ip)
        call wrtrst('L_scale',trim(c5p)//' kk time',L_scale,ip)
        call wrtrst('difdia',trim(c5p)//' kk time',difdia,ip)
        call wrtrst('ustarb',trim(c5p)//' time',ustarb,ip)
      end if
      if (use_IDLAGE) then
        call wrtrst('idlage',trim(c5p)//' kk2 time', &
             trc(1-nbdy,1-nbdy,1,itriag),ip)
      end if
    end if

    ! --- write accumulated fields
    do n = 1,nphy
      write(c2,'(i2.2)') n
      if (nacc_phy(n) /= 0) then
        if (acc_ub(n)     /= 0) call wrtrst('ub_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_UB(n)),iuu)
        if (acc_vb(n)     /= 0) call wrtrst('vb_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_VB(n)),ivv)
        if (acc_ubflxs(n)     /= 0) call wrtrst('ubflxs_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(n)),iuu)
        if (acc_vbflxs(n)     /= 0) call wrtrst('vbflxs_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(n)),ivv)
        if (acc_ztx(n)    /= 0) call wrtrst('ztx_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_ZTX(n)),iuu)
        if (acc_mty(n)    /= 0) call wrtrst('mty_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTY(n)),ivv)
        if (acc_taux(n)   /= 0) call wrtrst('taux_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_TAUX(n)),iuu)
        if (acc_tauy(n)   /= 0) call wrtrst('tauy_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_TAUY(n)),ivv)
        if (acc_uice(n)   /= 0) call wrtrst('uice_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_UICE(n)),iuu)
        if (acc_vice(n)   /= 0) call wrtrst('vice_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_VICE(n)),ivv)
        if (acc_ivolu(n)   /= 0) call wrtrst('ivolu_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_IVOLU(n)),iuu)
        if (acc_ivolv(n)   /= 0) call wrtrst('ivolv_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_IVOLV(n)),ivv)
        if (acc_psrf(n)  /= 0) call wrtrst('psrf_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_PSRF(n)),ip)
        if (acc_pbot(n)  /= 0) call wrtrst('pbot_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_PBOT(n)),ip)
        if (acc_sealv(n)  /= 0) call wrtrst('sealv_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SEALV(n)),ip)
        if (acc_slvsq(n)  /= 0) call wrtrst('slvsq_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SLVSQ(n)),ip)
        if (acc_sss(n)    /= 0) call wrtrst('sss_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SSS(n)),ip)
        if (acc_ssssq(n)    /= 0) call wrtrst('ssssq_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SSSSQ(n)),ip)
        if (acc_sbot(n)    /= 0) call wrtrst('sbot_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SBOT(n)),ip)
        if (acc_sst(n)    /= 0) call wrtrst('sst_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SST(n)),ip)
        if (acc_sstsq(n)    /= 0) call wrtrst('sstsq_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SSTSQ(n)),ip)
        if (acc_tbot(n)    /= 0) call wrtrst('tbot_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_TBOT(n)),ip)
        if (acc_sigmx(n)  /= 0) call wrtrst('sigmx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SIGMX(n)),ip)
        if (acc_mld(n)    /= 0) call wrtrst('mld_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MLD(n)),ip)
        if (acc_maxmld(n) /= 0) call wrtrst('maxmld_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MAXMLD(n)),ip)
        if (acc_mlts(n)    /= 0) call wrtrst('mlts_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MLTS(n)),ip)
        if (acc_mltsmn(n)    /= 0) call wrtrst('mltsmn_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MLTSMN(n)),ip)
        if (acc_mltsmx(n)    /= 0) call wrtrst('mltsmx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MLTSMX(n)),ip)
        if (acc_mltssq(n)    /= 0) call wrtrst('mltssq_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MLTSSQ(n)),ip)
        if (acc_t20d(n)    /= 0) call wrtrst('t20d_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_T20D(n)),ip)
        if (acc_alb(n)    /= 0) call wrtrst('alb_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_ALB(n)),ip)
        if (acc_swa(n)    /= 0) call wrtrst('swa_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SWA(n)),ip)
        if (acc_nsf(n)    /= 0) call wrtrst('nsf_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_NSF(n)),ip)
        if (acc_dfl(n)    /= 0) call wrtrst('dfl_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_DFL(n)),ip)
        if (acc_lip(n)    /= 0) call wrtrst('lip_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_LIP(n)),ip)
        if (acc_sop(n)    /= 0) call wrtrst('sop_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SOP(n)),ip)
        if (acc_eva(n)    /= 0) call wrtrst('eva_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_EVA(n)),ip)
        if (acc_rnfflx(n) /= 0) call wrtrst('rnfflx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_RNFFLX(n)),ip)
        if (acc_rfiflx(n) /= 0) call wrtrst('rfiflx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_RFIFLX(n)),ip)
        if (acc_sfl(n)    /= 0) call wrtrst('sfl_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SFL(n)),ip)
        if (acc_brnflx(n)    /= 0) call wrtrst('brnflx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_BRNFLX(n)),ip)
        if (acc_brnpd(n)    /= 0) call wrtrst('brnpd_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_BRNPD(n)),ip)
        if (acc_surflx(n) /= 0) call wrtrst('surflx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SURFLX(n)),ip)
        if (acc_surrlx(n) /= 0) call wrtrst('surrlx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SURRLX(n)),ip)
        if (acc_salflx(n) /= 0) call wrtrst('salflx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SALFLX(n)),ip)
        if (acc_salrlx(n) /= 0) call wrtrst('salrlx_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_SALRLX(n)),ip)
        if (acc_abswnd(n)    /= 0) call wrtrst('abswnd_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_ABSWND(n)),ip)
        if (acc_ustar(n)  /= 0) call wrtrst('ustar_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_USTAR(n)),ip)
        if (acc_ustar3(n)  /= 0) call wrtrst('ustar3_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_USTAR3(n)),ip)
        if (acc_idkedt(n)  /= 0) call wrtrst('idkedt_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_IDKEDT(n)),ip)
        if (acc_mtkeus(n)  /= 0) call wrtrst('mtkeus_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTKEUS(n)),ip)
        if (acc_mtkeni(n)  /= 0) call wrtrst('mtkeni_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTKENI(n)),ip)
        if (acc_mtkebf(n)  /= 0) call wrtrst('mtkebf_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTKEBF(n)),ip)
        if (acc_mtkers(n)  /= 0) call wrtrst('mtkers_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTKERS(n)),ip)
        if (acc_mtkepe(n)  /= 0) call wrtrst('mtkepe_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTKEPE(n)),ip)
        if (acc_mtkeke(n)  /= 0) call wrtrst('mtkeke_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_MTKEKE(n)),ip)
        if (acc_lamult(n)  /= 0) call wrtrst('lamult_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_LAMULT(n)),ip)
        if (acc_lasl(n)  /= 0) call wrtrst('lasl_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_LASL(n)),ip)
        if (acc_ustokes(n)  /= 0) call wrtrst('ustokes_phy'//c2, &
             trim(c5u)//' time',phyh2d(1-nbdy,1-nbdy,ACC_USTOKES(n)),iuu)
        if (acc_vstokes(n)  /= 0) call wrtrst('vstokes_phy'//c2, &
             trim(c5v)//' time',phyh2d(1-nbdy,1-nbdy,ACC_VSTOKES(n)),ivv)
        if (acc_fmltfz(n)    /= 0) call wrtrst('fmltfz_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_FMLTFZ(n)),ip)
        if (acc_hmltfz(n)    /= 0) call wrtrst('hmltfz_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_HMLTFZ(n)),ip)
        if (acc_hice(n)   /= 0) call wrtrst('hice_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_HICE(n)),ip)
        if (acc_hsnw(n)   /= 0) call wrtrst('hsnw_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_HSNW(n)),ip)
        if (acc_fice(n)   /= 0) call wrtrst('fice_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_FICE(n)),ip)
        if (acc_tsrf(n)   /= 0) call wrtrst('tsrf_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_TSRF(n)),ip)
        if (acc_tice(n)   /= 0) call wrtrst('tice_phy'//c2, &
             trim(c5p)//' time',phyh2d(1-nbdy,1-nbdy,ACC_TICE(n)),ip)
        if (acc_uvel(n)   /= 0) call wrtrst('uvel_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)), &
             iuu)
        if (acc_vvel(n)   /= 0) call wrtrst('vvel_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)), &
             ivv)
        if (acc_dpu(n)    /= 0) call wrtrst('dpu_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)), &
             iuu)
        if (acc_dpv(n)    /= 0) call wrtrst('dpv_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)), &
             ivv)
        if (acc_uflx(n)   /= 0) call wrtrst('uflx_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)), &
             iuu)
        if (acc_vflx(n)   /= 0) call wrtrst('vflx_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)), &
             ivv)
        if (acc_utflx(n)  /= 0) call wrtrst('utflx_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)), &
             iuu)
        if (acc_vtflx(n)  /= 0) call wrtrst('vtflx_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)), &
             ivv)
        if (acc_usflx(n)  /= 0) call wrtrst('usflx_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)), &
             iuu)
        if (acc_vsflx(n)  /= 0) call wrtrst('vsflx_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)), &
             ivv)
        if (acc_umfltd(n) /= 0) call wrtrst('umfltd_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)), &
             iuu)
        if (acc_vmfltd(n) /= 0) call wrtrst('vmfltd_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)), &
             ivv)
        if (acc_umflsm(n) /= 0) call wrtrst('umflsm_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UMFLSM(n)), &
             iuu)
        if (acc_vmflsm(n) /= 0) call wrtrst('vmflsm_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VMFLSM(n)), &
             ivv)
        if (acc_utfltd(n) /= 0) call wrtrst('utfltd_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)), &
             iuu)
        if (acc_vtfltd(n) /= 0) call wrtrst('vtfltd_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)), &
             ivv)
        if (acc_utflsm(n) /= 0) call wrtrst('utflsm_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLSM(n)), &
             iuu)
        if (acc_vtflsm(n) /= 0) call wrtrst('vtflsm_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLSM(n)), &
             ivv)
        if (acc_utflld(n) /= 0) call wrtrst('utflld_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)), &
             iuu)
        if (acc_vtflld(n) /= 0) call wrtrst('vtflld_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)), &
             ivv)
        if (acc_usfltd(n) /= 0) call wrtrst('usfltd_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)), &
             iuu)
        if (acc_vsfltd(n) /= 0) call wrtrst('vsfltd_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)), &
             ivv)
        if (acc_usflsm(n) /= 0) call wrtrst('usflsm_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLSM(n)), &
             iuu)
        if (acc_vsflsm(n) /= 0) call wrtrst('vsflsm_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLSM(n)), &
             ivv)
        if (acc_usflld(n) /= 0) call wrtrst('usflld_phy'//c2, &
             trim(c5u)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)), &
             iuu)
        if (acc_vsflld(n) /= 0) call wrtrst('vsflld_phy'//c2, &
             trim(c5v)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)), &
             ivv)
        if (acc_saln(n)   /= 0) call wrtrst('saln_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)), &
             ip)
        if (acc_temp(n)   /= 0) call wrtrst('temp_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)), &
             ip)
        if (acc_dp(n)     /= 0) call wrtrst('dp_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)),ip)
        if (acc_dz(n)     /= 0) call wrtrst('dz_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)),ip)
        if (acc_bfsq(n) /= 0) call wrtrst('bfsq_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_BFSQ(n)), &
             ip)
        if (acc_difdia(n) /= 0) call wrtrst('difdia_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)), &
             ip)
        if (acc_difvmo(n) /= 0) call wrtrst('difvmo_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFVMO(n)), &
             ip)
        if (acc_difvho(n) /= 0) call wrtrst('difvho_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFVHO(n)), &
             ip)
        if (acc_difvso(n) /= 0) call wrtrst('difvso_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFVSO(n)), &
             ip)
        if (acc_difint(n) /= 0) call wrtrst('difint_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)), &
             ip)
        if (acc_difiso(n) /= 0) call wrtrst('difiso_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)), &
             ip)
        if (acc_wflx(n)   /= 0) call wrtrst('wflx_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)), &
             ip)
        if (acc_wflx2(n)  /= 0) call wrtrst('wflx2_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)), &
             ip)
        if (acc_avdsg(n)  /= 0) call wrtrst('avdsg_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_AVDSG(n)), &
             ip)
        if (acc_dpvor(n)  /= 0) call wrtrst('dpvor_phy'//c2, &
             trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_DPVOR(n)), &
             ip)
        if (use_TRC .and. use_TKE) then
          if (acc_tke(n)   /= 0) call wrtrst('tke_phy'//c2, &
               trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)), &
               ip)
          if (acc_gls_psi(n) /= 0) call wrtrst('gls_psi_phy'//c2, &
               trim(c5p)//' kk time',phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)), &
               ip)
        end if
        if (acc_uvellvl(n)  /= 0) call wrtrst('uvellvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UVELLVL(n)),iuu)
        if (acc_vvellvl(n)  /= 0) call wrtrst('vvellvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VVELLVL(n)),ivv)
        if (acc_uflxlvl(n)  /= 0) call wrtrst('uflxlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UFLXLVL(n)),iuu)
        if (acc_vflxlvl(n)  /= 0) call wrtrst('vflxlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1 &
             ,ACC_VFLXLVL(n)),ivv)
        if (acc_utflxlvl(n) /= 0) call wrtrst('utflxlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UTFLXLVL(n)),iuu)
        if (acc_vtflxlvl(n) /= 0) call wrtrst('vtflxlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VTFLXLVL(n)),ivv)
        if (acc_usflxlvl(n) /= 0) call wrtrst('usflxlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_USFLXLVL(n)),iuu)
        if (acc_vsflxlvl(n) /= 0) call wrtrst('vsflxlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VSFLXLVL(n)),ivv)
        if (acc_umfltdlvl(n) /= 0) call wrtrst('umfltdlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UMFLTDLVL(n)),iuu)
        if (acc_vmfltdlvl(n) /= 0) call wrtrst('vmfltdlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VMFLTDLVL(n)),ivv)
        if (acc_umflsmlvl(n) /= 0) call wrtrst('umflsmlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UMFLSMLVL(n)),iuu)
        if (acc_vmflsmlvl(n) /= 0) call wrtrst('vmflsmlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VMFLSMLVL(n)),ivv)
        if (acc_utfltdlvl(n) /= 0) call wrtrst('utfltdlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UTFLTDLVL(n)),iuu)
        if (acc_vtfltdlvl(n) /= 0) call wrtrst('vtfltdlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VTFLTDLVL(n)),ivv)
        if (acc_utflsmlvl(n) /= 0) call wrtrst('utflsmlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UTFLSMLVL(n)),iuu)
        if (acc_vtflsmlvl(n) /= 0) call wrtrst('vtflsmlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VTFLSMLVL(n)),ivv)
        if (acc_utflldlvl(n) /= 0) call wrtrst('utflldlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_UTFLLDLVL(n)),iuu)
        if (acc_vtflldlvl(n) /= 0) call wrtrst('vtflldlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VTFLLDLVL(n)),ivv)
        if (acc_usfltdlvl(n) /= 0) call wrtrst('usfltdlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_USFLTDLVL(n)),iuu)
        if (acc_vsfltdlvl(n) /= 0) call wrtrst('vsfltdlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VSFLTDLVL(n)),ivv)
        if (acc_usflsmlvl(n) /= 0) call wrtrst('usflsmlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_USFLSMLVL(n)),iuu)
        if (acc_vsflsmlvl(n) /= 0) call wrtrst('vsflsmlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VSFLSMLVL(n)),ivv)
        if (acc_usflldlvl(n) /= 0) call wrtrst('usflldlvl_phy'//c2, &
             trim(c5u)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_USFLLDLVL(n)),iuu)
        if (acc_vsflldlvl(n) /= 0) call wrtrst('vsflldlvl_phy'//c2, &
             trim(c5v)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_VSFLLDLVL(n)),ivv)
        if (acc_salnlvl(n)  /= 0) call wrtrst('salnlvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_SALNLVL(n)),ip)
        if (acc_templvl(n)  /= 0) call wrtrst('templvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_TEMPLVL(n)),ip)
        if (acc_dzlvl(n)     /= 0) call wrtrst('dzlvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DZLVL(n)),ip)
        if (acc_bfsqlvl(n) /= 0) call wrtrst('bfsqlvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_BFSQLVL(n)),ip)
        if (acc_difdialvl(n) /= 0) call wrtrst('difdialvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DIFDIALVL(n)),ip)
        if (acc_difvmolvl(n) /= 0) call wrtrst('difvmolvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DIFVMOLVL(n)),ip)
        if (acc_difvholvl(n) /= 0) call wrtrst('difvholvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DIFVHOLVL(n)),ip)
        if (acc_difvsolvl(n) /= 0) call wrtrst('difvsolvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DIFVSOLVL(n)),ip)
        if (acc_difintlvl(n) /= 0) call wrtrst('difintlvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DIFINTLVL(n)),ip)
        if (acc_difisolvl(n) /= 0) call wrtrst('difisolvl_phy'//c2, &
             trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
             ACC_DIFISOLVL(n)),ip)
        if (acc_wflxlvl(n) /= 0) call wrtrst('wflxlvl_phy'//c2, &
             trim(c5p)//' plev time', &
             phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)),ip)
        if (acc_wflx2lvl(n) /= 0) call wrtrst('wflx2lvl_phy'//c2, &
             trim(c5p)//' plev time', &
             phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)),ip)
        if (acc_pvlvl(n) /= 0) call wrtrst('pvlvl_phy'//c2, &
             trim(c5p)//' plev time', &
             phylvl(1-nbdy,1-nbdy,1,ACC_PVLVL(n)),ip)
        if (use_TRC .and. use_TKE) then
          if (acc_tkelvl(n) /= 0) call wrtrst('tkelvl_phy'//c2, &
               trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
               ACC_TKELVL(n)),ip)
          if(acc_gls_psilvl(n) /= 0) call wrtrst('gls_psilvl_phy'//c2, &
               trim(c5p)//' plev time',phylvl(1-nbdy,1-nbdy,1, &
               ACC_GLS_PSILVL(n)),ip)
        end if
      end if
    end do

    call ncfcls

    if (use_TRC) then
      call restart_trcwt(rstfnm)
    end if

    if (ditflx) then

      ! --- - write diag. heat flux restart file
      if (expcnf == 'cesm') then
        fnm = trim(runid)//'.blom.rtflx.'//rstfnm(len_trim(runid)+10:)
      else
        fnm = trim(runid)//'_tflx_'//rstfnm(len_trim(runid)+10:)
      end if
      if (mnproc == 1) &
           write (lp,'(a,a)') ' saving diag. heat flux restart file ', &
           trim(fnm)
      if (rstfmt == 1) then
        call ncfopn(fnm,'w','6',1,iotype)
      else if (rstfmt == 2) then
        call ncfopn(fnm,'w','h',1,iotype)
      else
        call ncfopn(fnm,'w','c',1,iotype)
      end if

      if (rstcmp == 1) then
        call ncdimc('pcomp',ip,0)
      else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
      end if
      call ncdims('week',48)

      call ncputr('time',time)
      call ncputi('nflxdi',nflxdi)
      call defvarrst('tflxdi',trim(c5p)//' week')
      call ncedef

      call wrtrst('tflxdi',trim(c5p)//' week',tflxdi,ip)

      call ncfcls
    end if

    if (disflx) then

      ! --- - write diag. salt flux restart file
      if (expcnf == 'cesm') then
        fnm = trim(runid)//'.blom.rsflx.'//rstfnm(len_trim(runid)+10:)
      else
        fnm = trim(runid)//'_sflx_'//rstfnm(len_trim(runid)+10:)
      end if
      if (mnproc == 1) &
           write (lp,'(a,a)') ' saving diag. salt flux restart file ', &
           trim(fnm)
      if (rstfmt == 1) then
        call ncfopn(fnm,'w','6',1,iotype)
      else if (rstfmt == 2) then
        call ncfopn(fnm,'w','h',1,iotype)
      else
        call ncfopn(fnm,'w','c',1,iotype)
      end if

      if (rstcmp == 1) then
        call ncdimc('pcomp',ip,0)
      else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
      end if
      call ncdims('week',48)

      call ncputr('time',time)
      call ncputi('nflxdi',nflxdi)
      call defvarrst('sflxdi',trim(c5p)//' week')
      call ncedef

      call wrtrst('sflxdi',trim(c5p)//' week',sflxdi,ip)

      call ncfcls
    end if

    if (expcnf == 'cesm'.or.expcnf == 'channel') then

      ! --- - write restart filename to rpointer.ocn
      if (mnproc == 1) then
        open (unit=nfu,file = 'rpointer.ocn'//trim(inst_suffix))
        write (nfu,'(a)') rstfnm
        close (unit = nfu)
      end if
    end if

  end subroutine restart_wt

  ! ------------------------------------------------------------------
  subroutine wrtrst(vnm,dims,fld,msk)
    ! Arguments
    character(len=*), intent(in) :: vnm
    character(len=*), intent(in) :: dims
    real  , dimension(*), intent(in) :: fld
    integer,dimension(*), intent(in) :: msk

    ! --- Write data in compressed or uncompressed format
    if (dims(2:5) == 'comp') then
      call nccomp(vnm,dims,fld,msk,1.,0.,8)
    else
      call ncwrtr(vnm,dims,fld,msk,1,1.,0.,8)
    end if
  end subroutine wrtrst

  ! ------------------------------------------------------------------
  subroutine defvar_restart(c5p,c5u,c5v,c5q)
    ! Arguments
    character(len = 5), intent(in) :: c5p,c5u,c5v,c5q

    ! Local variables
    integer :: n
    character(len = 2) :: c2

    call defvarrst('u',trim(c5u)//' kk2 time')
    call defvarrst('v',trim(c5v)//' kk2 time')
    call defvarrst('dp',trim(c5p)//' kk2 time')
    call defvarrst('dpold',trim(c5p)//' kk2 time')
    call defvarrst('temp',trim(c5p)//' kk2 time')
    call defvarrst('saln',trim(c5p)//' kk2 time')
    call defvarrst('sigma',trim(c5p)//' kk2 time')
    call defvarrst('sigmar',trim(c5p)//' kk time')
    call defvarrst('pgfx',trim(c5u)//' kk2 time')
    call defvarrst('pgfy',trim(c5v)//' kk2 time')
    call defvarrst('pb',trim(c5p)//' k2 time')
    call defvarrst('pb_mn',trim(c5p)//' k2 time')
    call defvarrst('pb_p',trim(c5p)//' time')
    call defvarrst('pbu',trim(c5u)//' k2 time')
    call defvarrst('pbv',trim(c5v)//' k2 time')
    call defvarrst('pbu_p',trim(c5u)//' time')
    call defvarrst('pbv_p',trim(c5v)//' time')
    call defvarrst('ub',trim(c5u)//' k2 time')
    call defvarrst('vb',trim(c5v)//' k2 time')
    call defvarrst('uflx',trim(c5u)//' kk2 time')
    call defvarrst('utflx',trim(c5u)//' kk2 time')
    call defvarrst('usflx',trim(c5u)//' kk2 time')
    call defvarrst('umfltd',trim(c5u)//' kk2 time')
    call defvarrst('utfltd',trim(c5u)//' kk2 time')
    call defvarrst('utflld',trim(c5u)//' kk2 time')
    call defvarrst('usfltd',trim(c5u)//' kk2 time')
    call defvarrst('usflld',trim(c5u)//' kk2 time')
    call defvarrst('vflx',trim(c5v)//' kk2 time')
    call defvarrst('vtflx',trim(c5v)//' kk2 time')
    call defvarrst('vsflx',trim(c5v)//' kk2 time')
    call defvarrst('vmfltd',trim(c5v)//' kk2 time')
    call defvarrst('vtfltd',trim(c5v)//' kk2 time')
    call defvarrst('vtflld',trim(c5v)//' kk2 time')
    call defvarrst('vsfltd',trim(c5v)//' kk2 time')
    call defvarrst('vsflld',trim(c5v)//' kk2 time')
    call defvarrst('ubflx',trim(c5u)//' k2 time')
    call defvarrst('vbflx',trim(c5v)//' k2 time')
    call defvarrst('ubflx_mn',trim(c5u)//' k2 time')
    call defvarrst('vbflx_mn',trim(c5v)//' k2 time')
    call defvarrst('ubflxs',trim(c5u)//' k3 time')
    call defvarrst('vbflxs',trim(c5v)//' k3 time')
    call defvarrst('ubflxs_p',trim(c5u)//' k2 time')
    call defvarrst('vbflxs_p',trim(c5v)//' k2 time')
    call defvarrst('ubcors_p',trim(c5u)//' time')
    call defvarrst('vbcors_p',trim(c5v)//' time')
    call defvarrst('pvtrop',trim(c5q)//' k2 time')
    call defvarrst('pgfxm',trim(c5u)//' k2 time')
    call defvarrst('pgfym',trim(c5v)//' k2 time')
    call defvarrst('xixp',trim(c5u)//' k2 time')
    call defvarrst('xixm',trim(c5u)//' k2 time')
    call defvarrst('xiyp',trim(c5v)//' k2 time')
    call defvarrst('xiym',trim(c5v)//' k2 time')
    call defvarrst('phi',trim(c5p)//' time')
    call defvarrst('sealv',trim(c5p)//' time')
    call defvarrst('ustar',trim(c5p)//' time')
    call defvarrst('kfpla',trim(c5p)//' k2 time')
    call defvarrst('ficem',trim(c5p)//' time')
    call defvarrst('uml',trim(c5u)//' k4 time')
    call defvarrst('vml',trim(c5v)//' k4 time')
    call defvarrst('umlres',trim(c5u)//' k2 time')
    call defvarrst('vmlres',trim(c5v)//' k2 time')

    if (vcoord_type_tag == isopyc_bulkml) then
      call defvarrst('buoyfl',trim(c5p)//' time')
    end if

    if (vcoord_type_tag == cntiso_hybrid) then
      call defvarrst('dpu',trim(c5p)//' kk2 time')
      call defvarrst('dpv',trim(c5p)//' kk2 time')
      call defvarrst('difiso',trim(c5p)//' kk time')
      call defvarrst('OBLdepth',trim(c5p)//' time')
      call defvarrst('t_ns_nonloc',trim(c5p)//' kkp1 time')
      call defvarrst('s_nb_nonloc',trim(c5p)//' kkp1 time')
      call defvarrst('mu_nonloc',trim(c5u)//' kkp1 time')
      call defvarrst('mv_nonloc',trim(c5v)//' kkp1 time')
      call defvarrst('umflsm',trim(c5u)//' kk2 time')
      call defvarrst('utflsm',trim(c5u)//' kk2 time')
      call defvarrst('usflsm',trim(c5u)//' kk2 time')
      call defvarrst('vmflsm',trim(c5v)//' kk2 time')
      call defvarrst('vtflsm',trim(c5v)//' kk2 time')
      call defvarrst('vsflsm',trim(c5v)//' kk2 time')
    end if

    if (sprfac) then
      call defvarrst('eiacc',trim(c5p)//' time')
      call defvarrst('pracc',trim(c5p)//' time')
    end if

    if (expcnf == 'ben02clim'.or.expcnf == 'ben02syn') then
      call defvarrst('cd_d',trim(c5p)//' time')
      call defvarrst('ch_d',trim(c5p)//' time')
      call defvarrst('ce_d',trim(c5p)//' time')
      call defvarrst('wg2_d',trim(c5p)//' time')
      call defvarrst('cd_m',trim(c5p)//' time')
      call defvarrst('ch_m',trim(c5p)//' time')
      call defvarrst('ce_m',trim(c5p)//' time')
      call defvarrst('wg2_m',trim(c5p)//' time')
      call defvarrst('rhoa',trim(c5p)//' time')
      call defvarrst('tsi_tda',trim(c5p)//' time')
      call defvarrst('tml_tda',trim(c5p)//' time')
      call defvarrst('sml_tda',trim(c5p)//' time')
      call defvarrst('alb_tda',trim(c5p)//' time')
      call defvarrst('fice_tda',trim(c5p)//' time')
      call defvarrst('hicem',trim(c5p)//' time')
      call defvarrst('tsrfm',trim(c5p)//' time')
      call defvarrst('hsnwm',trim(c5p)//' time')
      call defvarrst('ticem',trim(c5p)//' time')
      call defvarrst('iagem',trim(c5p)//' time')
      call defvarrst('rnfres',trim(c5p)//' time')
    end if

    !      if (expcnf.eq.'channel') then
    !        call defvarrst('rnfres',trim(c5p)//' time')
    !      endif

    if (expcnf == 'cesm') then
      call defvarrst('ustarw_da',trim(c5p)//' k2 time')
      call defvarrst('ztx_da',trim(c5p)//' k2 time')
      call defvarrst('mty_da',trim(c5p)//' k2 time')
      call defvarrst('lip_da',trim(c5p)//' k2 time')
      call defvarrst('sop_da',trim(c5p)//' k2 time')
      call defvarrst('eva_da',trim(c5p)//' k2 time')
      call defvarrst('rnf_da',trim(c5p)//' k2 time')
      call defvarrst('rfi_da',trim(c5p)//' k2 time')
      call defvarrst('fmltfz_da',trim(c5p)//' k2 time')
      call defvarrst('sfl_da',trim(c5p)//' k2 time')
      call defvarrst('swa_da',trim(c5p)//' k2 time')
      call defvarrst('nsf_da',trim(c5p)//' k2 time')
      call defvarrst('hmlt_da',trim(c5p)//' k2 time')
      call defvarrst('slp_da',trim(c5p)//' k2 time')
      call defvarrst('ficem_da',trim(c5p)//' k2 time')
      call defvarrst('abswnd_da',trim(c5p)//' k2 time')
      call defvarrst('atmco2_da',trim(c5p)//' k2 time')
      call defvarrst('atmbrf_da',trim(c5p)//' k2 time')
      call defvarrst('frzpot',trim(c5p)//' time')
      call defvarrst('mltpot',trim(c5p)//' time')
      call defvarrst('flxco2',trim(c5p)//' time')
      call defvarrst('flxdms',trim(c5p)//' time')
#ifdef HAMOCC
      if (use_BROMO) then
        call defvarrst('flxbrf',trim(c5p)//' time')
      end if
#endif
    end if

    if (use_TRC) then
      if (use_TKE) then
        call defvarrst('tke',trim(c5p)//' kk2 time')
        call defvarrst('gls_psi',trim(c5p)//' kk2 time')
        call defvarrst('L_scale',trim(c5p)//' kk time')
        call defvarrst('difdia',trim(c5p)//' kk time')
        call defvarrst('ustarb',trim(c5p)//' time')
      end if
      if (use_IDLAGE) then
        call defvarrst('idlage',trim(c5p)//' kk2 time')
      end if
    end if

    ! --- write accumulated fields
    do n = 1,nphy
      write(c2,'(i2.2)') n
      if (nacc_phy(n) /= 0) then
        if (acc_ub(n)     /= 0) call defvarrst('ub_phy'//c2, &
             trim(c5u)//' time')
        if (acc_vb(n)     /= 0) call defvarrst('vb_phy'//c2, &
             trim(c5v)//' time')
        if (acc_ubflxs(n)     /= 0) call defvarrst('ubflxs_phy'//c2, &
             trim(c5u)//' time')
        if (acc_vbflxs(n)     /= 0) call defvarrst('vbflxs_phy'//c2, &
             trim(c5v)//' time')
        if (acc_ztx(n)    /= 0) call defvarrst('ztx_phy'//c2, &
             trim(c5u)//' time')
        if (acc_mty(n)    /= 0) call defvarrst('mty_phy'//c2, &
             trim(c5v)//' time')
        if (acc_taux(n)   /= 0) call defvarrst('taux_phy'//c2, &
             trim(c5u)//' time')
        if (acc_tauy(n)   /= 0) call defvarrst('tauy_phy'//c2, &
             trim(c5v)//' time')
        if (acc_uice(n)   /= 0) call defvarrst('uice_phy'//c2, &
             trim(c5u)//' time')
        if (acc_vice(n)   /= 0) call defvarrst('vice_phy'//c2, &
             trim(c5v)//' time')
        if (acc_ivolu(n)   /= 0) call defvarrst('ivolu_phy'//c2, &
             trim(c5u)//' time')
        if (acc_ivolv(n)   /= 0) call defvarrst('ivolv_phy'//c2, &
             trim(c5v)//' time')
        if (acc_psrf(n)  /= 0) call defvarrst('psrf_phy'//c2, &
             trim(c5p)//' time')
        if (acc_pbot(n)  /= 0) call defvarrst('pbot_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sealv(n)  /= 0) call defvarrst('sealv_phy'//c2, &
             trim(c5p)//' time')
        if (acc_slvsq(n)  /= 0) call defvarrst('slvsq_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sss(n)    /= 0) call defvarrst('sss_phy'//c2, &
             trim(c5p)//' time')
        if (acc_ssssq(n)    /= 0) call defvarrst('ssssq_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sbot(n)    /= 0) call defvarrst('sbot_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sst(n)    /= 0) call defvarrst('sst_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sstsq(n)    /= 0) call defvarrst('sstsq_phy'//c2, &
             trim(c5p)//' time')
        if (acc_tbot(n)    /= 0) call defvarrst('tbot_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sigmx(n)  /= 0) call defvarrst('sigmx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mld(n)    /= 0) call defvarrst('mld_phy'//c2, &
             trim(c5p)//' time')
        if (acc_maxmld(n) /= 0) call defvarrst('maxmld_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mlts(n)    /= 0) call defvarrst('mlts_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mltsmn(n)    /= 0) call defvarrst('mltsmn_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mltsmx(n)    /= 0) call defvarrst('mltsmx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mltssq(n)    /= 0) call defvarrst('mltssq_phy'//c2, &
             trim(c5p)//' time')
        if (acc_t20d(n)    /= 0) call defvarrst('t20d_phy'//c2, &
             trim(c5p)//' time')
        if (acc_alb(n)    /= 0) call defvarrst('alb_phy'//c2, &
             trim(c5p)//' time')
        if (acc_swa(n)    /= 0) call defvarrst('swa_phy'//c2, &
             trim(c5p)//' time')
        if (acc_nsf(n)    /= 0) call defvarrst('nsf_phy'//c2, &
             trim(c5p)//' time')
        if (acc_dfl(n)    /= 0) call defvarrst('dfl_phy'//c2, &
             trim(c5p)//' time')
        if (acc_lip(n)    /= 0) call defvarrst('lip_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sop(n)    /= 0) call defvarrst('sop_phy'//c2, &
             trim(c5p)//' time')
        if (acc_eva(n)    /= 0) call defvarrst('eva_phy'//c2, &
             trim(c5p)//' time')
        if (acc_rnfflx(n) /= 0) call defvarrst('rnfflx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_rfiflx(n) /= 0) call defvarrst('rfiflx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_sfl(n)    /= 0) call defvarrst('sfl_phy'//c2, &
             trim(c5p)//' time')
        if (acc_brnflx(n)    /= 0) call defvarrst('brnflx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_brnpd(n)    /= 0) call defvarrst('brnpd_phy'//c2, &
             trim(c5p)//' time')
        if (acc_surflx(n) /= 0) call defvarrst('surflx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_surrlx(n) /= 0) call defvarrst('surrlx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_salflx(n) /= 0) call defvarrst('salflx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_salrlx(n) /= 0) call defvarrst('salrlx_phy'//c2, &
             trim(c5p)//' time')
        if (acc_abswnd(n)    /= 0) call defvarrst('abswnd_phy'//c2, &
             trim(c5p)//' time')
        if (acc_ustar(n)  /= 0) call defvarrst('ustar_phy'//c2, &
             trim(c5p)//' time')
        if (acc_ustar3(n)  /= 0) call defvarrst('ustar3_phy'//c2, &
             trim(c5p)//' time')
        if (acc_idkedt(n)  /= 0) call defvarrst('idkedt_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mtkeus(n)  /= 0) call defvarrst('mtkeus_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mtkeni(n)  /= 0) call defvarrst('mtkeni_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mtkebf(n)  /= 0) call defvarrst('mtkebf_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mtkers(n)  /= 0) call defvarrst('mtkers_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mtkepe(n)  /= 0) call defvarrst('mtkepe_phy'//c2, &
             trim(c5p)//' time')
        if (acc_mtkeke(n)  /= 0) call defvarrst('mtkeke_phy'//c2, &
             trim(c5p)//' time')
        if (acc_lamult(n)  /= 0) call defvarrst('lamult_phy'//c2, &
             trim(c5p)//' time')
        if (acc_lasl(n)  /= 0) call defvarrst('lasl_phy'//c2, &
             trim(c5p)//' time')
        if (acc_ustokes(n)  /= 0) call defvarrst('ustokes_phy'//c2, &
             trim(c5u)//' time')
        if (acc_vstokes(n)  /= 0) call defvarrst('vstokes_phy'//c2, &
             trim(c5v)//' time')
        if (acc_fmltfz(n)    /= 0) call defvarrst('fmltfz_phy'//c2, &
             trim(c5p)//' time')
        if (acc_hmltfz(n)    /= 0) call defvarrst('hmltfz_phy'//c2, &
             trim(c5p)//' time')
        if (acc_hice(n)   /= 0) call defvarrst('hice_phy'//c2, &
             trim(c5p)//' time')
        if (acc_hsnw(n)   /= 0) call defvarrst('hsnw_phy'//c2, &
             trim(c5p)//' time')
        if (acc_fice(n)   /= 0) call defvarrst('fice_phy'//c2, &
             trim(c5p)//' time')
        if (acc_tsrf(n)   /= 0) call defvarrst('tsrf_phy'//c2, &
             trim(c5p)//' time')
        if (acc_tice(n)   /= 0) call defvarrst('tice_phy'//c2, &
             trim(c5p)//' time')
        if (acc_uvel(n)   /= 0) call defvarrst('uvel_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vvel(n)   /= 0) call defvarrst('vvel_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_dpu(n)    /= 0) call defvarrst('dpu_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_dpv(n)    /= 0) call defvarrst('dpv_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_uflx(n)   /= 0) call defvarrst('uflx_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vflx(n)   /= 0) call defvarrst('vflx_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_utflx(n)  /= 0) call defvarrst('utflx_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vtflx(n)  /= 0) call defvarrst('vtflx_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_usflx(n)  /= 0) call defvarrst('usflx_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vsflx(n)  /= 0) call defvarrst('vsflx_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_umfltd(n) /= 0) call defvarrst('umfltd_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vmfltd(n) /= 0) call defvarrst('vmfltd_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_umflsm(n) /= 0) call defvarrst('umflsm_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vmflsm(n) /= 0) call defvarrst('vmflsm_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_utfltd(n) /= 0) call defvarrst('utfltd_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vtfltd(n) /= 0) call defvarrst('vtfltd_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_utflsm(n) /= 0) call defvarrst('utflsm_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vtflsm(n) /= 0) call defvarrst('vtflsm_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_utflld(n) /= 0) call defvarrst('utflld_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vtflld(n) /= 0) call defvarrst('vtflld_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_usfltd(n) /= 0) call defvarrst('usfltd_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vsfltd(n) /= 0) call defvarrst('vsfltd_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_usflsm(n) /= 0) call defvarrst('usflsm_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vsflsm(n) /= 0) call defvarrst('vsflsm_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_usflld(n) /= 0) call defvarrst('usflld_phy'//c2, &
             trim(c5u)//' kk time')
        if (acc_vsflld(n) /= 0) call defvarrst('vsflld_phy'//c2, &
             trim(c5v)//' kk time')
        if (acc_saln(n)   /= 0) call defvarrst('saln_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_temp(n)   /= 0) call defvarrst('temp_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_dp(n)     /= 0) call defvarrst('dp_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_dz(n)     /= 0) call defvarrst('dz_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_bfsq(n)   /= 0) call defvarrst('bfsq_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_difdia(n) /= 0) call defvarrst('difdia_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_difvmo(n) /= 0) call defvarrst('difvmo_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_difvho(n) /= 0) call defvarrst('difvho_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_difvso(n) /= 0) call defvarrst('difvso_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_difint(n) /= 0) call defvarrst('difint_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_difiso(n) /= 0) call defvarrst('difiso_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_wflx(n)   /= 0) call defvarrst('wflx_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_wflx2(n)  /= 0) call defvarrst('wflx2_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_avdsg(n)  /= 0) call defvarrst('avdsg_phy'//c2, &
             trim(c5p)//' kk time')
        if (acc_dpvor(n)  /= 0) call defvarrst('dpvor_phy'//c2, &
             trim(c5p)//' kk time')
        if (use_TRC .and. use_TKE) then
          if (acc_tke(n)    /= 0) call defvarrst('tke_phy'//c2, &
               trim(c5p)//' kk time')
          if (acc_gls_psi(n) /= 0) call defvarrst('gls_psi_phy'//c2, &
               trim(c5p)//' kk time')
        end if
        if (acc_uvellvl(n)  /= 0) call defvarrst('uvellvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vvellvl(n)  /= 0) call defvarrst('vvellvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_uflxlvl(n)  /= 0) call defvarrst('uflxlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vflxlvl(n)  /= 0) call defvarrst('vflxlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_utflxlvl(n) /= 0) call defvarrst('utflxlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vtflxlvl(n) /= 0) call defvarrst('vtflxlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_usflxlvl(n) /= 0) call defvarrst('usflxlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vsflxlvl(n) /= 0) call defvarrst('vsflxlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_umfltdlvl(n) /= 0) call defvarrst('umfltdlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vmfltdlvl(n) /= 0) call defvarrst('vmfltdlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_umflsmlvl(n) /= 0) call defvarrst('umflsmlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vmflsmlvl(n) /= 0) call defvarrst('vmflsmlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_utfltdlvl(n) /= 0) call defvarrst('utfltdlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vtfltdlvl(n) /= 0) call defvarrst('vtfltdlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_utflsmlvl(n) /= 0) call defvarrst('utflsmlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vtflsmlvl(n) /= 0) call defvarrst('vtflsmlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_utflldlvl(n) /= 0) call defvarrst('utflldlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vtflldlvl(n) /= 0) call defvarrst('vtflldlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_usfltdlvl(n) /= 0) call defvarrst('usfltdlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vsfltdlvl(n) /= 0) call defvarrst('vsfltdlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_usflsmlvl(n) /= 0) call defvarrst('usflsmlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vsflsmlvl(n) /= 0) call defvarrst('vsflsmlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_usflldlvl(n) /= 0) call defvarrst('usflldlvl_phy'//c2, &
             trim(c5u)//' plev time')
        if (acc_vsflldlvl(n) /= 0) call defvarrst('vsflldlvl_phy'//c2, &
             trim(c5v)//' plev time')
        if (acc_salnlvl(n)  /= 0) call defvarrst('salnlvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_templvl(n)  /= 0) call defvarrst('templvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_dzlvl(n)     /= 0) call defvarrst('dzlvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_bfsqlvl(n) /= 0) call defvarrst('bfsqlvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_difdialvl(n) /= 0) call defvarrst('difdialvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_difvmolvl(n) /= 0) call defvarrst('difvmolvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_difvholvl(n) /= 0) call defvarrst('difvholvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_difvsolvl(n) /= 0) call defvarrst('difvsolvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_difintlvl(n) /= 0) call defvarrst('difintlvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_difisolvl(n) /= 0) call defvarrst('difisolvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_wflxlvl(n) /= 0) call defvarrst('wflxlvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_wflx2lvl(n) /= 0) call defvarrst('wflx2lvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (acc_pvlvl(n) /= 0) call defvarrst('pvlvl_phy'//c2, &
             trim(c5p)//' plev time')
        if (use_TRC .and. use_TKE) then
          if (acc_tkelvl(n) /= 0) call defvarrst('tkelvl_phy'//c2, &
               trim(c5p)//' plev time')
          if(acc_gls_psilvl(n) /= 0) call defvarrst('gls_psilvl_phy'//c2, &
               trim(c5p)//' plev time')
        end if
      end if
    end do

    call ncedef

  end subroutine defvar_restart

  ! ------------------------------------------------------------------
  subroutine defvarrst(vnm,dims)
    ! Arguments
    character(len=*), intent(in) :: vnm
    character(len=*), intent(in) :: dims

    call ncdefvar(vnm,dims,ndouble,8)
  end subroutine defvarrst

end module mod_restart_wt
