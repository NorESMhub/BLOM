! ------------------------------------------------------------------------------
! Copyright (C) 2006-2024 Mats Bentsen, Mehmet Ilicak, Alok Kumar Gupta,
!                         Ingo Bethke, Jerry Tjiputra, Ping-Gin Chiu,
!                         Aleksi Nummelin, Jörg Schwinger, Mariana Vertenstein, !                         Joeran Maerz
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

module mod_restart

   use mod_types,          only: r8
   use mod_config,         only: expcnf, runid, inst_suffix, resume_flag
   use mod_calendar,       only: date_type, daynum_diff, operator(/=), &
                                 calendar_noerr, calendar_errstr
   use mod_time,           only: date0, date, nday1, nstep0, nstep1, nstep, time, time0, &
                                 nstep_in_day, nday_of_year, calendar
   use mod_xc
   use mod_vcoord,         only: vcoord_tag, vcoord_isopyc_bulkml, sigmar
   use mod_inicon,         only: icfile
   use mod_state,          only: u, v, dp, dpu, dpv, temp, saln, sigma, uflx, vflx, &
                                 utflx, vtflx, usflx, vsflx, phi, ubflxs, vbflxs, &
                                 ub, vb, pb, pbu, pbv, ubflxs_p, vbflxs_p, &
                                 pb_p, pbu_p, pbv_p, ubcors_p, vbcors_p, sealv, kfpla
   use mod_pgforc,         only: pgfx, pgfy, pgfxm, pgfym, xixp, xixm, xiyp, xiym
   use mod_barotp,         only: ubflx, vbflx, pb_mn, ubflx_mn, vbflx_mn, pvtrop
   use mod_temmin,         only: settemmin
   use mod_nctools,        only: nccomp, ncwrtr, ncdefvar, ndouble, ncinqv, ncinqa, &
                                 ncread, ncfopn, ncgetr, ncfopn, ncgeti, ncgetr, ncfcls, &
                                 ncputi, ncputr, ncdims, ncdimc, ncedef
   use mod_dia,            only: rstfrq, iotype, rstfmt, rstcmp, rstmon, &
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
   use mod_forcing,        only: ditflx, disflx, sprfac, tflxdi, sflxdi, nflxdi, &
                                 prfac, eiacc, pracc, flxco2, flxdms, flxbrf, &
                                 flxn2o,flxnh3, &
                                 ustarb, wstar3, buoyfl, ustar
   use mod_niw,            only: uml, vml, umlres, vmlres
   use mod_difest,         only: OBLdepth
   use mod_diffusion,      only: difiso, Kvisc_m, Kdiff_t, Kdiff_s, &
                                 t_ns_nonloc, s_nb_nonloc, mu_nonloc, mv_nonloc, difdia,&
                                 umflsm, usfltd, utflld, usflsm, utfltd, &
                                 usflld, utflsm, usflld, utflld, umfltd, usflld, &
                                 vmflsm, vsfltd, vtflld, vsflsm, vtfltd, &
                                 vsflld, vtflsm, vsflld, vtflld, vmfltd, vsflld
   use mod_eddtra,         only: tau_growing_hbl, tau_decaying_hbl, &
                                 tau_growing_hml, tau_decaying_hml, &
                                 hbl_tf, wpup_tf, hml_tf1, hml_tf
   use mod_cesm,           only: frzpot, mltpot, swa_da, nsf_da, hmlt_da, lip_da, &
                                 sop_da, eva_da, rnf_da, rfi_da, fmltfz_da, sfl_da, &
                                 ztx_da, mty_da, ustarw_da, slp_da, abswnd_da, &
                                 atmco2_da, atmbrf_da, &
                                 atmco2_da, atmbrf_da, atmn2o_da, atmnh3_da, &
                                 ficem_da, l1ci, l2ci
   use mod_ben02,          only: cd_d, ch_d, ce_d, wg2_d, cd_m, ch_m, ce_m, wg2_m, &
                                 rhoa, tsi_tda, tml_tda, sml_tda, alb_tda, fice_tda, &
                                 ntda, rnfres
   use mod_thdysi,         only: tsrfm, ticem
   use mod_seaice,         only: ficem, hicem, hsnwm, iagem
   use mod_tmsmt,          only: dpold
   use mod_tracers,        only: itrtke, itrgls, itriag, trc
   use mod_tke,            only: L_scale
#ifdef HAMOCC
   use mo_control_bgc,     only: use_BROMO, use_extNcycle
#endif
   use mod_ifdefs,         only: use_TRC, use_TKE, use_IDLAGE, use_MKS
   use mod_idlage,         only: idlage_init
   use mod_tracers_update, only: restart_trcwt, restart_trcrd

   implicit none

   private

   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: rkfpla
   integer , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu, ivv, iqq
   character(len = 5) :: c5p, c5u, c5v, c5q
   logical :: first_call = .true.

   ! Unit conversion factors in case restart file and model use diffent unit
   ! systems. All set to 1 initially.
   real(r8), parameter :: &
      no_unitconv   = 1._r8    ! No conversion.
   real(r8) :: &
      l_unitconv    = 1._r8, & ! Length conversion.
      m_unitconv    = 1._r8, & ! Mass conversion.
      p_unitconv    = 1._r8, & ! Pressure conversion.
      r_unitconv    = 1._r8, & ! Density conversion.
      l2_unitconv   = 1._r8, & ! Length squared conversion.
      l3_unitconv   = 1._r8, & ! Length cubed conversion.
      lm_unitconv   = 1._r8, & ! Length times mass conversion.
      l2m2_unitconv = 1._r8, & ! Length squared times mass squared conversion.
      pi_unitconv   = 1._r8, & ! One over pressure conversion.
      ri_unitconv   = 1._r8, & ! One over density conversion.
      l2i_unitconv  = 1._r8, & ! One over length squared conversion.
      ml2i_unitconv = 1._r8    ! Mass over length squared conversion.

   public  :: restart_write, restart_read

   private :: extended_masks
   private :: defwrtfld
   private :: defwrtflds
   private :: readfld

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   subroutine extended_masks
   ! ---------------------------------------------------------------------------
   ! Compute extended masks at u-, v- and q-points.
   ! ---------------------------------------------------------------------------

      integer :: i, j

      if (first_call) then
         first_call = .false.
      !$omp parallel do private(i)
         do j = 1, jj
            do i = 1, ii
               if ((ip(i,j) + ip(i-1,j)) >= 1) then
                  iuu(i,j) = 1
               else
                  iuu(i,j) = 0
               endif
               if ((ip(i,j) + ip(i,j-1)) >= 1) then
                  ivv(i,j) = 1
               else
                  ivv(i,j) = 0
               endif
               if ((iu(i,j) + iv(i,j) + iu(i,j-1) + iv(i-1,j)) >= 1) then
                  iqq(i,j) = 1
               else
                  iqq(i,j) = 0
               endif
            enddo
         enddo
      !$omp end parallel do
      endif

   end subroutine extended_masks

   subroutine defwrtfld(vnm, dims, fld, msk, defmode)
   ! ---------------------------------------------------------------------------
   ! Depending on argument defmode, define or write field to netCDF file.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: vnm, dims
      real(r8), dimension(*), intent(in) :: fld
      integer, dimension(*), intent(in) :: msk
      logical, intent(in) :: defmode

      if (defmode) then
         ! Define field.
         call ncdefvar(vnm, dims, ndouble, 8)
      else
         if (dims(2:5) == 'comp') then
            ! Write field in compressed format.
            call nccomp(vnm, dims, fld, msk, 1._r8, 0._r8, 8)
         else
            ! Write field in uncompressed format.
            call ncwrtr(vnm, dims, fld, msk, 1, 1._r8, 0._r8, 8)
         endif
      endif

   end subroutine defwrtfld

   subroutine defwrtflds(defmode)
   ! ---------------------------------------------------------------------------
   ! Depending on argument defmode, define or write fields to netCDF file.
   ! ---------------------------------------------------------------------------

      logical, intent(in) :: defmode

      integer :: n
      character(len = 2) :: c2

      call defwrtfld('u', trim(c5u)//' kk2 time', &
                      u, iuu, defmode)
      call defwrtfld('v', trim(c5v)//' kk2 time', &
                      v, ivv, defmode)
      call defwrtfld('dp', trim(c5p)//' kk2 time', &
                      dp, ip, defmode)
      call defwrtfld('dpold', trim(c5p)//' kk2 time', &
                      dpold, ip, defmode)
      call defwrtfld('temp', trim(c5p)//' kk2 time', &
                      temp, ip, defmode)
      call defwrtfld('saln', trim(c5p)//' kk2 time', &
                      saln, ip, defmode)
      call defwrtfld('sigma', trim(c5p)//' kk2 time', &
                      sigma, ip, defmode)
      call defwrtfld('sigmar', trim(c5p)//' kk time', &
                      sigmar, ip, defmode)
      call defwrtfld('pgfx', trim(c5u)//' kk2 time', &
                      pgfx, iuu, defmode)
      call defwrtfld('pgfy', trim(c5v)//' kk2 time', &
                      pgfy, ivv, defmode)
      call defwrtfld('pb', trim(c5p)//' k2 time', &
                      pb, ip, defmode)
      call defwrtfld('pb_mn', trim(c5p)//' k2 time', &
                      pb_mn, ip, defmode)
      call defwrtfld('pb_p', trim(c5p)//' time', &
                      pb_p, ip, defmode)
      call defwrtfld('pbu', trim(c5u)//' k2 time', &
                      pbu, iuu, defmode)
      call defwrtfld('pbv', trim(c5v)//' k2 time', &
                      pbv, ivv, defmode)
      call defwrtfld('pbu_p', trim(c5u)//' time', &
                      pbu_p, iuu, defmode)
      call defwrtfld('pbv_p', trim(c5v)//' time', &
                      pbv_p, ivv, defmode)
      call defwrtfld('ub', trim(c5u)//' k2 time',&
                      ub, iuu, defmode)
      call defwrtfld('vb', trim(c5v)//' k2 time', &
                      vb, ivv, defmode)
      call defwrtfld('uflx', trim(c5u)//' kk2 time', &
                      uflx, iuu, defmode)
      call defwrtfld('utflx', trim(c5u)//' kk2 time', &
                      utflx, iuu, defmode)
      call defwrtfld('usflx', trim(c5u)//' kk2 time', &
                      usflx, iuu, defmode)
      call defwrtfld('umfltd', trim(c5u)//' kk2 time', &
                      umfltd, iuu, defmode)
      call defwrtfld('utfltd', trim(c5u)//' kk2 time', &
                      utfltd, iuu, defmode)
      call defwrtfld('utflld', trim(c5u)//' kk2 time', &
                      utflld, iuu, defmode)
      call defwrtfld('usfltd', trim(c5u)//' kk2 time', &
                      usfltd, iuu, defmode)
      call defwrtfld('usflld', trim(c5u)//' kk2 time', &
                      usflld, iuu, defmode)
      call defwrtfld('vflx', trim(c5v)//' kk2 time', &
                      vflx, ivv, defmode)
      call defwrtfld('vtflx', trim(c5v)//' kk2 time', &
                      vtflx, ivv, defmode)
      call defwrtfld('vsflx', trim(c5v)//' kk2 time', &
                      vsflx, ivv, defmode)
      call defwrtfld('vmfltd', trim(c5v)//' kk2 time', &
                      vmfltd, ivv, defmode)
      call defwrtfld('vtfltd', trim(c5v)//' kk2 time', &
                      vtfltd, ivv, defmode)
      call defwrtfld('vtflld', trim(c5v)//' kk2 time', &
                      vtflld, ivv, defmode)
      call defwrtfld('vsfltd', trim(c5v)//' kk2 time', &
                      vsfltd, ivv, defmode)
      call defwrtfld('vsflld', trim(c5v)//' kk2 time', &
                      vsflld, ivv, defmode)
      call defwrtfld('ubflx', trim(c5u)//' k2 time', &
                      ubflx, iuu, defmode)
      call defwrtfld('vbflx', trim(c5v)//' k2 time', &
                      vbflx, ivv, defmode)
      call defwrtfld('ubflx_mn', trim(c5u)//' k2 time', &
                      ubflx_mn, iuu, defmode)
      call defwrtfld('vbflx_mn', trim(c5v)//' k2 time', &
                      vbflx_mn, ivv, defmode)
      call defwrtfld('ubflxs', trim(c5u)//' k3 time', &
                      ubflxs, iuu, defmode)
      call defwrtfld('vbflxs', trim(c5v)//' k3 time', &
                      vbflxs, ivv, defmode)
      call defwrtfld('ubflxs_p', trim(c5u)//' k2 time', &
                      ubflxs_p, iuu, defmode)
      call defwrtfld('vbflxs_p', trim(c5v)//' k2 time', &
                      vbflxs_p, ivv, defmode)
      call defwrtfld('ubcors_p', trim(c5u)//' time', &
                      ubcors_p, iuu, defmode)
      call defwrtfld('vbcors_p', trim(c5v)//' time', &
                      vbcors_p, ivv, defmode)
      call defwrtfld('pvtrop', trim(c5q)//' k2 time', &
                      pvtrop, iqq, defmode)
      call defwrtfld('pgfxm', trim(c5u)//' k2 time', &
                      pgfxm, iuu, defmode)
      call defwrtfld('pgfym', trim(c5v)//' k2 time', &
                      pgfym, ivv, defmode)
      call defwrtfld('xixp', trim(c5u)//' k2 time', &
                      xixp, iuu, defmode)
      call defwrtfld('xixm', trim(c5u)//' k2 time', &
                      xixm, iuu, defmode)
      call defwrtfld('xiyp', trim(c5v)//' k2 time', &
                      xiyp, ivv, defmode)
      call defwrtfld('xiym', trim(c5v)//' k2 time', &
                      xiym, ivv, defmode)
      call defwrtfld('phi', trim(c5p)//' time', &
                      phi(1-nbdy,1-nbdy,kk+1), ip, defmode)
      call defwrtfld('sealv', trim(c5p)//' time', &
                      sealv, ip, defmode)
      call defwrtfld('ustar', trim(c5p)//' time', &
                      ustar, ip, defmode)
      call defwrtfld('kfpla', trim(c5p)//' k2 time', &
                      rkfpla, ip, defmode)
      call defwrtfld('ficem', trim(c5p)//' time', &
                      ficem, ip, defmode)

      if (vcoord_tag == vcoord_isopyc_bulkml) then
         call defwrtfld('buoyfl', trim(c5p)//' time', &
                         buoyfl, ip, defmode)
         call defwrtfld('uml', trim(c5u)//' k4 time', &
                         uml, iuu, defmode)
         call defwrtfld('vml', trim(c5v)//' k4 time', &
                         vml, ivv, defmode)
         call defwrtfld('umlres', trim(c5u)//' k2 time', &
                         umlres, iuu, defmode)
         call defwrtfld('vmlres', trim(c5v)//' k2 time', &
                         vmlres, ivv, defmode)
      endif

      if (vcoord_tag /= vcoord_isopyc_bulkml) then
         call defwrtfld('dpu', trim(c5u)//' kk2 time', &
                         dpu, iu, defmode)
         call defwrtfld('dpv', trim(c5v)//' kk2 time', &
                         dpv, iv, defmode)
         call defwrtfld('difiso', trim(c5p)//' kk time', &
                         difiso, ip, defmode)
         call defwrtfld('OBLdepth', trim(c5p)//' time', &
                         OBLdepth, ip, defmode)
         call defwrtfld('t_ns_nonloc', trim(c5p)//' kkp1 time', &
                         t_ns_nonloc, ip, defmode)
         call defwrtfld('s_nb_nonloc', trim(c5p)//' kkp1 time', &
                         s_nb_nonloc, ip, defmode)
         call defwrtfld('mu_nonloc', trim(c5u)//' kkp1 time', &
                         mu_nonloc, iu, defmode)
         call defwrtfld('mv_nonloc', trim(c5v)//' kkp1 time', &
                         mv_nonloc, iv, defmode)
         call defwrtfld('umflsm', trim(c5u)//' kk2 time', &
                         umflsm, iuu, defmode)
         call defwrtfld('utflsm', trim(c5u)//' kk2 time', &
                         utflsm, iuu, defmode)
         call defwrtfld('usflsm', trim(c5u)//' kk2 time', &
                         usflsm, iuu, defmode)
         call defwrtfld('vmflsm', trim(c5v)//' kk2 time', &
                         vmflsm, ivv, defmode)
         call defwrtfld('vtflsm', trim(c5v)//' kk2 time', &
                         vtflsm, ivv, defmode)
         call defwrtfld('vsflsm', trim(c5v)//' kk2 time', &
                         vsflsm, ivv, defmode)
         call defwrtfld('wstar3', trim(c5p)//' time', &
                         wstar3, ip, defmode)
         if (tau_growing_hbl > 0._r8 .or. tau_decaying_hbl > 0._r8) then
            call defwrtfld('hbl_tf', trim(c5p)//' time', &
                            hbl_tf, ip, defmode)
            call defwrtfld('wpup_tf', trim(c5p)//' time', &
                            wpup_tf, ip, defmode)
            call defwrtfld('hml_tf1', trim(c5p)//' time', &
                            hml_tf1, ip, defmode)
         endif
         if (tau_growing_hml > 0._r8 .or. tau_decaying_hml > 0._r8) then
            call defwrtfld('hml_tf', trim(c5p)//' time', &
                            hml_tf, ip, defmode)
         endif
      endif

      if (sprfac) then
         call defwrtfld('eiacc', trim(c5p)//' time', &
                         eiacc, ip, defmode)
         call defwrtfld('pracc', trim(c5p)//' time', &
                         pracc, ip, defmode)
      endif

      if (expcnf == 'ben02clim' .or. expcnf == 'ben02syn') then
         call defwrtfld('cd_d', trim(c5p)//' time', &
                         cd_d, ip, defmode)
         call defwrtfld('ch_d', trim(c5p)//' time', &
                         ch_d, ip, defmode)
         call defwrtfld('ce_d', trim(c5p)//' time', &
                         ce_d, ip, defmode)
         call defwrtfld('wg2_d', trim(c5p)//' time', &
                         wg2_d, ip, defmode)
         call defwrtfld('cd_m', trim(c5p)//' time', &
                         cd_m, ip, defmode)
         call defwrtfld('ch_m', trim(c5p)//' time', &
                        ch_m, ip, defmode)
         call defwrtfld('ce_m', trim(c5p)//' time', &
                         ce_m, ip, defmode)
         call defwrtfld('wg2_m', trim(c5p)//' time', &
                         wg2_m, ip, defmode)
         call defwrtfld('rhoa', trim(c5p)//' time', &
                         rhoa, ip, defmode)
         call defwrtfld('tsi_tda', trim(c5p)//' time', &
                         tsi_tda, ip, defmode)
         call defwrtfld('tml_tda', trim(c5p)//' time', &
                         tml_tda, ip, defmode)
         call defwrtfld('sml_tda', trim(c5p)//' time', &
                         sml_tda, ip, defmode)
         call defwrtfld('alb_tda', trim(c5p)//' time', &
                         alb_tda, ip, defmode)
         call defwrtfld('fice_tda', trim(c5p)//' time', &
                         fice_tda, ip, defmode)
         call defwrtfld('hicem', trim(c5p)//' time', &
                         hicem, ip, defmode)
         call defwrtfld('tsrfm', trim(c5p)//' time', &
                         tsrfm, ip, defmode)
         call defwrtfld('hsnwm', trim(c5p)//' time', &
                         hsnwm, ip, defmode)
         call defwrtfld('ticem', trim(c5p)//' time', &
                         ticem, ip, defmode)
         call defwrtfld('iagem', trim(c5p)//' time', &
                         iagem, ip, defmode)
         call defwrtfld('rnfres', trim(c5p)//' time', &
                         rnfres, ip, defmode)
      endif

!     if (expcnf == 'channel') then
!        call defwrtfld('rnfres', trim(c5p)//' time', &
!                        rnfres, ip, defmode)
!     endif

      if (expcnf == 'cesm') then
         call defwrtfld('ustarw_da', trim(c5p)//' k2 time', &
                         ustarw_da, ip, defmode)
         call defwrtfld('ztx_da', trim(c5p)//' k2 time', &
                         ztx_da, ip, defmode)
         call defwrtfld('mty_da', trim(c5p)//' k2 time', &
                         mty_da, ip, defmode)
         call defwrtfld('lip_da', trim(c5p)//' k2 time', &
                         lip_da, ip, defmode)
         call defwrtfld('sop_da', trim(c5p)//' k2 time', &
                         sop_da, ip, defmode)
         call defwrtfld('eva_da', trim(c5p)//' k2 time', &
                        eva_da, ip, defmode)
         call defwrtfld('rnf_da', trim(c5p)//' k2 time', &
                         rnf_da, ip, defmode)
         call defwrtfld('rfi_da', trim(c5p)//' k2 time', &
                         rfi_da, ip, defmode)
         call defwrtfld('fmltfz_da', trim(c5p)//' k2 time', &
                         fmltfz_da, ip, defmode)
         call defwrtfld('sfl_da', trim(c5p)//' k2 time', &
                         sfl_da, ip, defmode)
         call defwrtfld('swa_da', trim(c5p)//' k2 time', &
                         swa_da, ip, defmode)
         call defwrtfld('nsf_da', trim(c5p)//' k2 time', &
                         nsf_da, ip, defmode)
         call defwrtfld('hmlt_da', trim(c5p)//' k2 time', &
                         hmlt_da, ip, defmode)
         call defwrtfld('slp_da', trim(c5p)//' k2 time', &
                         slp_da, ip, defmode)
         call defwrtfld('ficem_da', trim(c5p)//' k2 time', &
                         ficem_da, ip, defmode)
         call defwrtfld('abswnd_da', trim(c5p)//' k2 time', &
                         abswnd_da, ip, defmode)
         call defwrtfld('atmco2_da', trim(c5p)//' k2 time', &
                         atmco2_da, ip, defmode)
         call defwrtfld('atmbrf_da', trim(c5p)//' k2 time', &
                         atmbrf_da, ip, defmode)  ! not read in restart_read, necessary?
         call defwrtfld('atmn2o_da', trim(c5p)//' k2 time', &
                         atmn2o_da, ip, defmode)  ! not read in restart_read, necessary?
         call defwrtfld('atmnh3_da', trim(c5p)//' k2 time', &
                         atmnh3_da, ip, defmode)  ! not read in restart_read, necessary?
         call defwrtfld('frzpot', trim(c5p)//' time', &
                         frzpot, ip, defmode)
         call defwrtfld('mltpot', trim(c5p)//' time', &
                         mltpot, ip, defmode)
         call defwrtfld('flxco2', trim(c5p)//' time', &
                         flxco2, ip, defmode)
         call defwrtfld('flxdms', trim(c5p)//' time', &
                         flxdms, ip, defmode)
#ifdef HAMOCC
         if (use_BROMO) then
            call defwrtfld('flxbrf', trim(c5p)//' time', flxbrf, ip, defmode)
         endif
         if (use_extNcycle) then
            call defwrtfld('flxn2o', trim(c5p)//' time', &
                            flxn2o, ip, defmode)
            call defwrtfld('flxnh3', trim(c5p)//' time', &
                            flxnh3, ip, defmode)
         endif
#endif
      endif

      if (use_TRC) then
         if (use_TKE) then
            call defwrtfld('tke', trim(c5p)//' kk2 time', trc(1-nbdy,1-nbdy,1,itrtke), ip, defmode)
            call defwrtfld('gls_psi', trim(c5p)//' kk2 time', trc(1-nbdy,1-nbdy,1,itrgls), ip, defmode)
            call defwrtfld('L_scale', trim(c5p)//' kk time', L_scale, ip, defmode)
            call defwrtfld('difdia', trim(c5p)//' kk time', difdia, ip, defmode)
            call defwrtfld('ustarb', trim(c5p)//' time', ustarb, ip, defmode)
         end if
         if (use_IDLAGE) then
            call defwrtfld('idlage', trim(c5p)//' kk2 time', trc(1-nbdy,1-nbdy,1,itriag), ip, defmode)
         end if
      end if

      ! Write accumulated fields.
      do n = 1, nphy
         write(c2, '(i2.2)') n
         if (nacc_phy(n) /= 0) then
            if (ACC_UB(n) /= 0) &
               call defwrtfld('ub_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_UB(n)), iuu, defmode)
            if (ACC_VB(n) /= 0) &
               call defwrtfld('vb_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_VB(n)), ivv, defmode)
            if (ACC_UBFLXS(n) /= 0) &
               call defwrtfld('ubflxs_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(n)), iuu, defmode)
            if (ACC_VBFLXS(n) /= 0) &
               call defwrtfld('vbflxs_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(n)), ivv, defmode)
            if (ACC_ZTX(n) /= 0) &
               call defwrtfld('ztx_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_ZTX(n)), iuu, defmode)
            if (ACC_MTY(n) /= 0) &
               call defwrtfld('mty_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTY(n)), ivv, defmode)
            if (ACC_TAUX(n) /= 0) &
               call defwrtfld('taux_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_TAUX(n)), iuu, defmode)
            if (ACC_TAUY(n) /= 0) &
               call defwrtfld('tauy_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_TAUY(n)), ivv, defmode)
            if (ACC_UICE(n) /= 0) &
               call defwrtfld('uice_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_UICE(n)), iuu, defmode)
            if (ACC_VICE(n) /= 0) &
               call defwrtfld('vice_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_VICE(n)), ivv, defmode)
            if (ACC_IVOLU(n) /= 0) &
               call defwrtfld('ivolu_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_IVOLU(n)), iuu, defmode)
            if (ACC_IVOLV(n) /= 0) &
               call defwrtfld('ivolv_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_IVOLV(n)), ivv, defmode)
            if (ACC_PSRF(n) /= 0) &
               call defwrtfld('psrf_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_PSRF(n)), ip, defmode)
            if (ACC_PBOT(n) /= 0) &
               call defwrtfld('pbot_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_PBOT(n)), ip, defmode)
            if (ACC_SEALV(n) /= 0) &
               call defwrtfld('sealv_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SEALV(n)), ip, defmode)
            if (ACC_SLVSQ(n) /= 0) &
               call defwrtfld('slvsq_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SLVSQ(n)), ip, defmode)
            if (ACC_SSS(n) /= 0) &
               call defwrtfld('sss_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SSS(n)), ip, defmode)
            if (ACC_SSSSQ(n) /= 0) &
               call defwrtfld('ssssq_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SSSSQ(n)), ip, defmode)
            if (ACC_SBOT(n) /= 0) &
               call defwrtfld('sbot_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SBOT(n)), ip, defmode)
            if (ACC_SST(n) /= 0) &
               call defwrtfld('sst_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SST(n)), ip, defmode)
            if (ACC_SSTSQ(n) /= 0) &
               call defwrtfld('sstsq_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SSTSQ(n)), ip, defmode)
            if (ACC_TBOT(n) /= 0) &
               call defwrtfld('tbot_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_TBOT(n)), ip, defmode)
            if (ACC_SIGMX(n) /= 0) &
               call defwrtfld('sigmx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SIGMX(n)), ip, defmode)
            if (ACC_MLD(n) /= 0) &
               call defwrtfld('mld_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MLD(n)), ip, defmode)
            if (ACC_MAXMLD(n) /= 0) &
               call defwrtfld('maxmld_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MAXMLD(n)), ip, defmode)
            if (ACC_MLTS(n) /= 0) &
               call defwrtfld('mlts_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MLTS(n)), ip, defmode)
            if (ACC_MLTSMN(n) /= 0) &
               call defwrtfld('mltsmn_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MLTSMN(n)), ip, defmode)
            if (ACC_MLTSMX(n) /= 0) &
               call defwrtfld('mltsmx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MLTSMX(n)), ip, defmode)
            if (ACC_MLTSSQ(n) /= 0) &
               call defwrtfld('mltssq_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MLTSSQ(n)), ip, defmode)
            if (ACC_T20D(n) /= 0) &
               call defwrtfld('t20d_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_T20D(n)), ip, defmode)
            if (ACC_ALB(n) /= 0) &
               call defwrtfld('alb_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_ALB(n)), ip, defmode)
            if (ACC_SWA(n) /= 0) &
               call defwrtfld('swa_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SWA(n)), ip, defmode)
            if (ACC_NSF(n) /= 0) &
               call defwrtfld('nsf_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_NSF(n)), ip, defmode)
            if (ACC_DFL(n) /= 0) &
               call defwrtfld('dfl_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_DFL(n)), ip, defmode)
            if (ACC_LIP(n) /= 0) &
               call defwrtfld('lip_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_LIP(n)), ip, defmode)
            if (ACC_SOP(n) /= 0) &
               call defwrtfld('sop_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SOP(n)), ip, defmode)
            if (ACC_EVA(n) /= 0) &
               call defwrtfld('eva_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_EVA(n)), ip, defmode)
            if (ACC_RNFFLX(n) /= 0) &
               call defwrtfld('rnfflx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_RNFFLX(n)), ip, defmode)
            if (ACC_RFIFLX(n) /= 0) &
               call defwrtfld('rfiflx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_RFIFLX(n)), ip, defmode)
            if (ACC_SFL(n) /= 0) &
               call defwrtfld('sfl_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SFL(n)), ip, defmode)
            if (ACC_BRNFLX(n) /= 0) &
               call defwrtfld('brnflx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_BRNFLX(n)), ip, defmode)
            if (ACC_BRNPD(n) /= 0) &
               call defwrtfld('brnpd_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_BRNPD(n)), ip, defmode)
            if (ACC_SURFLX(n) /= 0) &
               call defwrtfld('surflx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SURFLX(n)), ip, defmode)
            if (ACC_SURRLX(n) /= 0) &
               call defwrtfld('surrlx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SURRLX(n)), ip, defmode)
            if (ACC_SALFLX(n) /= 0) &
               call defwrtfld('salflx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SALFLX(n)), ip, defmode)
            if (ACC_SALRLX(n) /= 0) &
               call defwrtfld('salrlx_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_SALRLX(n)), ip, defmode)
            if (ACC_ABSWND(n) /= 0) &
               call defwrtfld('abswnd_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_ABSWND(n)), ip, defmode)
            if (ACC_USTAR(n) /= 0) &
               call defwrtfld('ustar_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_USTAR(n)), ip, defmode)
            if (ACC_USTAR3(n) /= 0) &
               call defwrtfld('ustar3_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_USTAR3(n)), ip, defmode)
            if (ACC_IDKEDT(n) /= 0) &
               call defwrtfld('idkedt_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_IDKEDT(n)), ip, defmode)
            if (ACC_MTKEUS(n) /= 0) &
               call defwrtfld('mtkeus_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTKEUS(n)), ip, defmode)
            if (ACC_MTKENI(n) /= 0) &
               call defwrtfld('mtkeni_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTKENI(n)), ip, defmode)
            if (ACC_MTKEBF(n) /= 0) &
               call defwrtfld('mtkebf_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTKEBF(n)), ip, defmode)
            if (ACC_MTKERS(n) /= 0) &
               call defwrtfld('mtkers_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTKERS(n)), ip, defmode)
            if (ACC_MTKEPE(n) /= 0) &
               call defwrtfld('mtkepe_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTKEPE(n)), ip, defmode)
            if (ACC_MTKEKE(n) /= 0) &
               call defwrtfld('mtkeke_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_MTKEKE(n)), ip, defmode)
            if (ACC_LAMULT(n) /= 0) &
               call defwrtfld('lamult_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_LAMULT(n)), ip, defmode)
            if (ACC_LASL(n) /= 0) &
               call defwrtfld('lasl_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_LASL(n)), ip, defmode)
            if (ACC_USTOKES(n) /= 0) &
               call defwrtfld('ustokes_phy'//c2, trim(c5u)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_USTOKES(n)), iuu, defmode)
            if (ACC_VSTOKES(n) /= 0) &
               call defwrtfld('vstokes_phy'//c2, trim(c5v)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_VSTOKES(n)), ivv, defmode)
            if (ACC_FMLTFZ(n) /= 0) &
               call defwrtfld('fmltfz_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_FMLTFZ(n)), ip, defmode)
            if (ACC_HMLTFZ(n) /= 0) &
               call defwrtfld('hmltfz_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_HMLTFZ(n)), ip, defmode)
            if (ACC_HICE(n) /= 0) &
               call defwrtfld('hice_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_HICE(n)), ip, defmode)
            if (ACC_HSNW(n) /= 0) &
               call defwrtfld('hsnw_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_HSNW(n)), ip, defmode)
            if (ACC_FICE(n) /= 0) &
               call defwrtfld('fice_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_FICE(n)), ip, defmode)
            if (ACC_TSRF(n) /= 0) &
               call defwrtfld('tsrf_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_TSRF(n)), ip, defmode)
            if (ACC_TICE(n) /= 0) &
               call defwrtfld('tice_phy'//c2, trim(c5p)//' time', &
                              phyh2d(1-nbdy,1-nbdy,ACC_TICE(n)), ip, defmode)
            if (ACC_UVEL(n) /= 0) &
               call defwrtfld('uvel_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)), iuu, defmode)
            if (ACC_VVEL(n) /= 0) &
                call defwrtfld('vvel_phy'//c2, trim(c5v)//' kk time', &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)), ivv, defmode)
            if (ACC_DPU(n) /= 0) &
               call defwrtfld('dpu_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)), iuu, defmode)
            if (ACC_DPV(n) /= 0) &
               call defwrtfld('dpv_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)), ivv, defmode)
            if (ACC_UFLX(n) /= 0) &
               call defwrtfld('uflx_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)), iuu, defmode)
            if (ACC_VFLX(n) /= 0) &
               call defwrtfld('vflx_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)), ivv, defmode)
            if (ACC_UTFLX(n) /= 0) &
               call defwrtfld('utflx_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)), iuu, defmode)
            if (ACC_VTFLX(n) /= 0) &
               call defwrtfld('vtflx_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)), ivv, defmode)
            if (ACC_USFLX(n) /= 0) &
               call defwrtfld('usflx_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)), iuu, defmode)
            if (ACC_VSFLX(n) /= 0) &
               call defwrtfld('vsflx_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)), ivv, defmode)
            if (ACC_UMFLTD(n) /= 0) &
               call defwrtfld('umfltd_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)), iuu, defmode)
            if (ACC_VMFLTD(n) /= 0) &
               call defwrtfld('vmfltd_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)), ivv, defmode)
            if (ACC_UMFLSM(n) /= 0) &
               call defwrtfld('umflsm_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UMFLSM(n)), iuu, defmode)
            if (ACC_VMFLSM(n) /= 0) &
               call defwrtfld('vmflsm_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VMFLSM(n)), ivv, defmode)
            if (ACC_UTFLTD(n) /= 0) &
               call defwrtfld('utfltd_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)), iuu, defmode)
            if (ACC_VTFLTD(n) /= 0) &
               call defwrtfld('vtfltd_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)), ivv, defmode)
            if (ACC_UTFLSM(n) /= 0) &
               call defwrtfld('utflsm_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UTFLSM(n)), iuu, defmode)
            if (ACC_VTFLSM(n) /= 0) &
               call defwrtfld('vtflsm_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VTFLSM(n)), ivv, defmode)
            if (ACC_UTFLLD(n) /= 0) &
               call defwrtfld('utflld_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)), iuu, defmode)
            if (ACC_VTFLLD(n) /= 0) &
               call defwrtfld('vtflld_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)), ivv, defmode)
            if (ACC_USFLTD(n) /= 0) &
               call defwrtfld('usfltd_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)), iuu, defmode)
            if (ACC_VSFLTD(n) /= 0) &
               call defwrtfld('vsfltd_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)), ivv, defmode)
            if (ACC_USFLSM(n) /= 0) &
               call defwrtfld('usflsm_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_USFLSM(n)), iuu, defmode)
            if (ACC_VSFLSM(n) /= 0) &
               call defwrtfld('vsflsm_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VSFLSM(n)), ivv, defmode)
            if (ACC_USFLLD(n) /= 0) &
               call defwrtfld('usflld_phy'//c2, trim(c5u)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)), iuu, defmode)
            if (ACC_VSFLLD(n) /= 0) &
               call defwrtfld('vsflld_phy'//c2, trim(c5v)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)), ivv, defmode)
            if (ACC_SALN(n) /= 0) &
               call defwrtfld('saln_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)), ip, defmode)
            if (ACC_TEMP(n) /= 0) &
               call defwrtfld('temp_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)), ip, defmode)
            if (ACC_DP(n) /= 0) &
               call defwrtfld('dp_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)), ip, defmode)
            if (ACC_DZ(n) /= 0) &
               call defwrtfld('dz_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)), ip, defmode)
            if (ACC_BFSQ(n) /= 0) &
               call defwrtfld('bfsq_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_BFSQ(n)), ip, defmode)
            if (ACC_DIFDIA(n) /= 0) &
               call defwrtfld('difdia_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)), ip, defmode)
            if (ACC_DIFVMO(n) /= 0) &
               call defwrtfld('difvmo_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DIFVMO(n)), ip, defmode)
            if (ACC_DIFVHO(n) /= 0) &
               call defwrtfld('difvho_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DIFVHO(n)), ip, defmode)
            if (ACC_DIFVSO(n) /= 0) &
               call defwrtfld('difvso_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DIFVSO(n)), ip, defmode)
            if (ACC_DIFINT(n) /= 0) &
               call defwrtfld('difint_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)), ip, defmode)
            if (ACC_DIFISO(n) /= 0) &
               call defwrtfld('difiso_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)), ip, defmode)
            if (ACC_WFLX(n) /= 0) &
               call defwrtfld('wflx_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)), ip, defmode)
            if (ACC_WFLX2(n) /= 0) &
               call defwrtfld('wflx2_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)), ip, defmode)
            if (ACC_AVDSG(n) /= 0) &
               call defwrtfld('avdsg_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_AVDSG(n)), ip, defmode)
            if (ACC_DPVOR(n) /= 0) &
               call defwrtfld('dpvor_phy'//c2, trim(c5p)//' kk time', &
                              phylyr(1-nbdy,1-nbdy,1,ACC_DPVOR(n)), ip, defmode)
            if (use_TRC .and. use_TKE) then
               if (ACC_TKE(n) /= 0) &
                    call defwrtfld('tke_phy'//c2, trim(c5p)//' kk time', &
                                   phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)), ip, defmode)
               if (ACC_GLS_PSI(n) /= 0) &
                    call defwrtfld('gls_psi_phy'//c2, trim(c5p)//' kk time', &
                                   phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)), ip, defmode)
            end if
            if (ACC_UVELLVL(n) /= 0) &
               call defwrtfld('uvellvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UVELLVL(n)), iuu, defmode)
            if (ACC_VVELLVL(n) /= 0) &
               call defwrtfld('vvellvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VVELLVL(n)), ivv, defmode)
            if (ACC_UFLXLVL(n) /= 0) &
               call defwrtfld('uflxlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)), iuu, defmode)
            if (ACC_VFLXLVL(n) /= 0) &
               call defwrtfld('vflxlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)), ivv, defmode)
            if (ACC_UTFLXLVL(n) /= 0) &
               call defwrtfld('utflxlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UTFLXLVL(n)), iuu, defmode)
            if (ACC_VTFLXLVL(n) /= 0) &
               call defwrtfld('vtflxlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VTFLXLVL(n)), ivv, defmode)
            if (ACC_USFLXLVL(n) /= 0) &
               call defwrtfld('usflxlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_USFLXLVL(n)), iuu, defmode)
            if (ACC_VSFLXLVL(n) /= 0) &
               call defwrtfld('vsflxlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VSFLXLVL(n)), ivv, defmode)
            if (ACC_UMFLTDLVL(n) /= 0) &
               call defwrtfld('umfltdlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UMFLTDLVL(n)), iuu, defmode)
            if (ACC_VMFLTDLVL(n) /= 0) &
               call defwrtfld('vmfltdlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VMFLTDLVL(n)), ivv, defmode)
            if (ACC_UMFLSMLVL(n) /= 0) &
               call defwrtfld('umflsmlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UMFLSMLVL(n)), iuu, defmode)
            if (ACC_VMFLSMLVL(n) /= 0) &
               call defwrtfld('vmflsmlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VMFLSMLVL(n)), ivv, defmode)
            if (ACC_UTFLTDLVL(n) /= 0) &
               call defwrtfld('utfltdlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UTFLTDLVL(n)), iuu, defmode)
            if (ACC_VTFLTDLVL(n) /= 0) &
               call defwrtfld('vtfltdlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VTFLTDLVL(n)), ivv, defmode)
            if (ACC_UTFLSMLVL(n) /= 0) &
               call defwrtfld('utflsmlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UTFLSMLVL(n)), iuu, defmode)
            if (ACC_VTFLSMLVL(n) /= 0) &
               call defwrtfld('vtflsmlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VTFLSMLVL(n)), ivv, defmode)
            if (ACC_UTFLLDLVL(n) /= 0) &
               call defwrtfld('utflldlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_UTFLLDLVL(n)), iuu, defmode)
            if (ACC_VTFLLDLVL(n) /= 0) &
               call defwrtfld('vtflldlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VTFLLDLVL(n)), ivv, defmode)
            if (ACC_USFLTDLVL(n) /= 0) &
               call defwrtfld('usfltdlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_USFLTDLVL(n)), iuu, defmode)
            if (ACC_VSFLTDLVL(n) /= 0) &
               call defwrtfld('vsfltdlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VSFLTDLVL(n)), ivv, defmode)
            if (ACC_USFLSMLVL(n) /= 0) &
               call defwrtfld('usflsmlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_USFLSMLVL(n)), iuu, defmode)
            if (ACC_VSFLSMLVL(n) /= 0) &
               call defwrtfld('vsflsmlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VSFLSMLVL(n)), ivv, defmode)
            if (ACC_USFLLDLVL(n) /= 0) &
               call defwrtfld('usflldlvl_phy'//c2, trim(c5u)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_USFLLDLVL(n)), iuu, defmode)
            if (ACC_VSFLLDLVL(n) /= 0) &
               call defwrtfld('vsflldlvl_phy'//c2, trim(c5v)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_VSFLLDLVL(n)), ivv, defmode)
            if (ACC_SALNLVL(n) /= 0) &
               call defwrtfld('salnlvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_SALNLVL(n)), ip, defmode)
            if (ACC_TEMPLVL(n) /= 0) &
               call defwrtfld('templvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_TEMPLVL(n)), ip, defmode)
            if (ACC_DZLVL(n) /= 0) &
               call defwrtfld('dzlvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DZLVL(n)), ip, defmode)
            if (ACC_BFSQLVL(n) /= 0) &
               call defwrtfld('bfsqlvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_BFSQLVL(n)), ip, defmode)
            if (ACC_DIFDIALVL(n) /= 0) &
               call defwrtfld('difdialvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DIFDIALVL(n)), ip, defmode)
            if (ACC_DIFVMOLVL(n) /= 0) &
               call defwrtfld('difvmolvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DIFVMOLVL(n)), ip, defmode)
            if (ACC_DIFVHOLVL(n) /= 0) &
               call defwrtfld('difvholvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DIFVHOLVL(n)), ip, defmode)
            if (ACC_DIFVSOLVL(n) /= 0) &
               call defwrtfld('difvsolvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DIFVSOLVL(n)), ip, defmode)
            if (ACC_DIFINTLVL(n) /= 0) &
               call defwrtfld('difintlvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DIFINTLVL(n)), ip, defmode)
            if (ACC_DIFISOLVL(n) /= 0) &
               call defwrtfld('difisolvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_DIFISOLVL(n)), ip, defmode)
            if (ACC_WFLXLVL(n) /= 0) &
               call defwrtfld('wflxlvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)), ip, defmode)
            if (ACC_WFLX2LVL(n) /= 0) &
               call defwrtfld('wflx2lvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)), ip, defmode)
            if (ACC_PVLVL(n) /= 0) &
               call defwrtfld('pvlvl_phy'//c2, trim(c5p)//' plev time', &
                              phylvl(1-nbdy,1-nbdy,1,ACC_PVLVL(n)), ip, defmode)
            if (use_TRC .and. use_TKE) then
               if (ACC_TKELVL(n) /= 0) &
                    call defwrtfld('tkelvl_phy'//c2, trim(c5p)//' plev time', &
                                   phylvl(1-nbdy,1-nbdy,1,ACC_TKELVL(n)), ip, defmode)
               if(ACC_GLS_PSILVL(n) /= 0) &
                    call defwrtfld('gls_psilvl_phy'//c2, trim(c5p)//' plev time', &
                                   phylvl(1-nbdy,1-nbdy,1,ACC_GLS_PSILVL(n)), ip, defmode)
            end if
         endif
      enddo

   end subroutine defwrtflds

   subroutine readfld(vnm, unitconv, fld, msk, required, fld_read)
   ! ---------------------------------------------------------------------------
   ! Read field from netCDF file.
   ! ---------------------------------------------------------------------------

      character(len = *), intent(in) :: vnm
      real(r8), intent(in) :: unitconv
      real(r8), dimension(*), intent(out) :: fld
      integer, dimension(*), intent(in) :: msk
      logical, intent(in), optional :: required
      logical, intent(out), optional :: fld_read

      logical :: required_loc

      ! Reading of a field is required by default.
      if (present(required)) then
         required_loc = required
      else
         required_loc = .true.
      endif

      if (ncinqv(vnm)) then
         call ncread(vnm, fld, msk, 1, 0._r8, unitconv)
         if (present(fld_read)) fld_read = .true.
      else
         if (required_loc) then
            if (mnproc == 1) &
               write(lp,*) 'readfld: could not read required field '// &
                           trim(vnm)//' from restart file!'
            call xcstop('(readfld)')
                   stop '(readfld)'
         else
            write(lp,*) 'readfld: field '//trim(vnm)// &
                        ' is not read from restart file.'
         endif
         if (present(fld_read)) fld_read = .false.
      endif

   end subroutine readfld

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine restart_write
   ! ---------------------------------------------------------------------------
   ! Write model state to restart files.
   ! ---------------------------------------------------------------------------

      integer :: i, j, n
      character(len = 256), dimension(4) :: rstdate_str
      character(len = 256) :: rstfnm, fnm
      character(len = 2) :: c2
      integer :: nfu

      ! Formulate restart filename.
      if (expcnf == 'cesm') then
         write(rstfnm, '(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') &
            trim(runid), '.blom', trim(inst_suffix), '.r.', &
            date%year, '-', date%month, '-', date%day, '-', &
            86400*mod(nstep, nstep_in_day)/nstep_in_day, '.nc'
      else
        if (nday_of_year-max(1, nint(rstfrq/nstep_in_day)) <= 0) then
           write(rstfnm, '(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)') &
              trim(runid), '_restphy_', &
              date%year, '.', date%month, '.', date%day, &
              '_', nint(time), '.nc'
        else
           if (rstmon) then
              write(rstfnm, '(2a,i1,a)') &
                 trim(runid), '_restphy_', &
                 mod(date%month + 10, 3) + 1, '.nc'
           else
              write(rstfnm, '(2a,i1,a)') &
                 trim(runid), '_restphy_', &
                 mod(nint(min(nstep/rstfrq, time)) - 1, 3) + 1, '.nc'
           endif
           open(newunit = nfu, file = 'rstdate.txt')
           i = 1
 300       read(nfu, '(a)', end = 301) rstdate_str(i)
           i = i + 1
           goto 300
 301       close(unit = nfu)
           write(rstdate_str(i), '(2a,i4.4,a,i2.2,a,i2.2,a,i6.6)') &
              rstfnm(1:len_trim(runid) + 13), &
              ': date ', date%year, '.', date%month, '.', date%day, &
              ', integration day ', nint(time)
           if (mnproc == 1) then
              if (i == 1) then
                 open(newunit = nfu, file = 'rstdate.txt')
                 write(nfu, '(a)') rstdate_str(1)(1:len_trim(runid) + 54)
                 close(unit = nfu)
              elseif (rstdate_str(max(1, i - 2)) /= rstdate_str(i) .and. &
                      rstdate_str(i - 1       ) /= rstdate_str(i)) then
                 open(newunit = nfu, file = 'rstdate.txt')
                 do j = max(1, i - 2), i
                    write(nfu, '(a)') rstdate_str(j)(1:len_trim(runid) + 54)
                 enddo
                 close(unit = nfu)
              endif
            endif
         endif
         if (mnproc == 1) &
            write(lp,*) 'restart_write: saving restart file '//trim(rstfnm)
      endif

      ! Compute extended masks at u-, v- and q-points.
      call extended_masks

      ! Open restart file.
      if (rstfmt == 1) then
         call ncfopn(rstfnm, 'w', '6', 1, iotype)
      elseif (rstfmt == 2) then
         call ncfopn(rstfnm, 'w', 'h', 1, iotype)
      else
         call ncfopn(rstfnm, 'w', 'c', 1, iotype)
      endif

      ! Write time variables.
      call ncputi('nday0', date0%day)
      call ncputi('nmonth0', date0%month)
      call ncputi('nyear0', date0%year)
      call ncputr('time0', time0)
      call ncputr('time', time)

      ! Define spatial and time dimensions.
      if (rstcmp == 1) then
         call ncdimc('pcomp', ip, 0)
         call ncdimc('qcomp', iqq, 0)
         call ncdimc('ucomp', iuu, 0)
         call ncdimc('vcomp', ivv, 0)
         c5p = 'pcomp'
         c5u = 'ucomp'
         c5v = 'vcomp'
         c5q = 'qcomp'
      else
         call ncdims('x', itdm)
         call ncdims('y', jtdm)
         c5p = 'x y'
         c5u = 'x y'
         c5v = 'x y'
         c5q = 'x y'
      endif
      call ncdims('k2', 2)
      call ncdims('k3', 3)
      call ncdims('k4', 4)
      call ncdims('kk', kk)
      call ncdims('kkp1', kk + 1)
      call ncdims('kk2', 2*kk)
      call ncdims('plev', ddm)
      call ncputr('plev', depthslev)
      call ncdims('time', 1)

      ! Output model fields to restart file.
      if (sprfac) then
         call ncputr('prfac', prfac)
      endif
      if (expcnf == 'ben02clim' .or. expcnf == 'ben02syn' .or. &
          expcnf == 'channel') then
         call ncputi('ntda', ntda)
      endif
      do n = 1, nphy
         write(c2, '(i2.2)') n
         call ncputi('nacc_phy'//c2, nacc_phy(n))
      enddo
      if (expcnf == 'cesm') then
         call ncputi('l2ci', l2ci)
      endif

   !$omp parallel do private(i)
      do j = 1, jj
         do i = 1, ii
            if (ip(i,j) == 1) then
               rkfpla(i,j,1) = real(kfpla(i,j,1))
               rkfpla(i,j,2) = real(kfpla(i,j,2))
            else
               rkfpla(i,j,1) = 0._r8
               rkfpla(i,j,2) = 0._r8
            endif
         enddo
      enddo
   !$omp end parallel do

      call defwrtflds(defmode = .true.)

      call ncedef

      call defwrtflds(defmode = .false.)

      call ncfcls

      if (use_TRC) then
         call restart_trcwt(rstfnm)
      end if

      if (ditflx) then

         ! Write diag. heat flux restart file.
         if (expcnf == 'cesm') then
            fnm = trim(runid)//'.blom.rtflx.'//rstfnm(len_trim(runid) + 10:)
         else
            fnm = trim(runid)//'_tflx_'//rstfnm(len_trim(runid) + 10:)
         endif
         if (mnproc == 1) &
            write(lp,*) &
               'restart_write: saving diag. heat flux restart file '//trim(fnm)
         if (rstfmt == 1) then
            call ncfopn(fnm, 'w', '6', 1, iotype)
         elseif (rstfmt == 2) then
            call ncfopn(fnm, 'w', 'h', 1, iotype)
         else
            call ncfopn(fnm, 'w', 'c', 1, iotype)
         endif

         if (rstcmp == 1) then
            call ncdimc('pcomp', ip, 0)
         else
            call ncdims('x', itdm)
            call ncdims('y', jtdm)
         endif
         call ncdims('week', 48)

         call ncputr('time', time)
         call ncputi('nflxdi', nflxdi)
         call defwrtfld('tflxdi', trim(c5p)//' week', tflxdi, ip, &
                        defmode = .true.)

         call ncedef

         call defwrtfld('tflxdi', trim(c5p)//' week', tflxdi, ip, &
                       defmode = .false.)

         call ncfcls
      endif

      if (disflx) then

         ! Write diag. salt flux restart file.
         if (expcnf == 'cesm') then
            fnm = trim(runid)//'.blom.rsflx.'//rstfnm(len_trim(runid) + 10:)
         else
            fnm = trim(runid)//'_sflx_'//rstfnm(len_trim(runid) + 10:)
         endif
         if (mnproc == 1) &
            write(lp,*) &
               'restart_write: saving diag. salt flux restart file '//trim(fnm)
         if (rstfmt == 1) then
            call ncfopn(fnm, 'w', '6', 1, iotype)
         elseif (rstfmt == 2) then
            call ncfopn(fnm, 'w', 'h', 1, iotype)
         else
            call ncfopn(fnm, 'w', 'c', 1, iotype)
         endif

         if (rstcmp == 1) then
            call ncdimc('pcomp', ip, 0)
         else
            call ncdims('x', itdm)
            call ncdims('y', jtdm)
         endif
         call ncdims('week', 48)

         call ncputr('time', time)
         call ncputi('nflxdi', nflxdi)
         call defwrtfld('sflxdi', trim(c5p)//' week', sflxdi, ip, &
                        defmode = .true.)

         call ncedef

         call defwrtfld('sflxdi', trim(c5p)//' week', sflxdi, ip, &
                        defmode = .false.)

         call ncfcls
      endif

      if (expcnf == 'cesm' .or. expcnf == 'channel') then
         ! Write restart filename to rpointer.ocn.
         if (mnproc == 1) then
            open(newunit = nfu, file = 'rpointer.ocn'//trim(inst_suffix))
            write(nfu, '(a)') rstfnm
            close(unit = nfu)
         endif
      endif

   end subroutine restart_write

   subroutine restart_read
   ! ---------------------------------------------------------------------------
   ! Read model state from restart files.
   ! ---------------------------------------------------------------------------

      type(date_type) :: date_rest
      integer errstat, dndiff, i, j, l, n
      character(len = 256) :: rstfnm, fnm
      character(len = 2) :: c2
      real(r8) :: pb_max, phi_min, rho_restart
      logical :: file_exist, fld_read
      integer :: nfu

      ! Open restart file and adjust time information if needed.
      if     (nday1 + nint(time0) == 0 .and. (.not.resume_flag)) then

         ! Open restart file for initial conditions and adjust integration time
         ! corresponding to start date.
         rstfnm = icfile
         if (mnproc == 1) inquire(file = rstfnm, exist = file_exist)
         call xcbcst(file_exist)
         if (file_exist) then
            call ncfopn(rstfnm, 'r', ' ', 1, iotype)
            call ncgetr('time', time)
            time0 = time
         else
            if (mnproc == 1) then
               write(lp,*) 'restart_read: could not find restart file for '// &
                           'initial conditions!'
            endif
            call xcstop('(restart_read)')
                   stop '(restart_read)'
         endif

      elseif (expcnf == 'cesm' .or. expcnf == 'channel') then

         ! Get restart file name from rpointer.ocn.
         if (mnproc == 1) &
            inquire(file = 'rpointer.ocn'//trim(inst_suffix), &
                    exist = file_exist)
         call xcbcst(file_exist)
         if (file_exist) then
            if (mnproc == 1) then
               open(newunit = nfu, file = 'rpointer.ocn'//trim(inst_suffix))
               read(nfu, '(a)') rstfnm
               close(unit = nfu)
            endif
            call xcbcst(rstfnm)
         else
            if (mnproc == 1) then
               write(lp,*) 'restart_read: could not find file rpointer.ocn'// &
                           trim(inst_suffix)//'!'
            endif
            call xcstop('(restart_read)')
                   stop '(restart_read)'
         endif

         ! Open restart file.
         if (mnproc == 1) inquire(file = rstfnm, exist = file_exist)
         call xcbcst(file_exist)
         if (.not.file_exist) then
            if (mnproc == 1) then
               write(lp,*) 'restart_read: could not find restart file!'
            endif
            call xcstop('(restart_read)')
                   stop '(restart_read)'
         endif
         call ncfopn(rstfnm, 'r', ' ', 1, iotype)

      else

         ! First try file name:
         !  <experiment name>_restphy_<year>.<month>.<day>_<integration day>.nc
         write(rstfnm, '(2a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)') &
           trim(runid), '_restphy_', &
           date%year, '.', date%month, '.', date%day, &
           '_', nday1, '.nc'

         if (mnproc == 1) inquire(file = rstfnm, exist = file_exist)
         call xcbcst(file_exist)
         if (file_exist) then
           call ncfopn(rstfnm, 'r', ' ', 1, iotype)
           call ncgeti('nday0', date_rest%day)
           call ncgeti('nmonth0', date_rest%month)
           call ncgeti('nyear0', date_rest%year)
           call ncgetr('time0', time0)
           call ncgetr('time', time)
           if (date_rest /= date0) then
             if (mnproc == 1) then
               write(lp,*) 'restart_read: expected identical initial '// &
                           'experiment date in namelist and restart but found:'
               write(lp, '(a,i4.4,2(i2.2))') &
                  ' restart_read: initial date namelist: ', date0
               write(lp, '(a,i4.4,2(i2.2))') &
                  ' restart_read: initial date restart:  ', date_rest
             endif
             call xcstop('(restart_read)')
                    stop '(restart_read)'
           endif

         else

           ! Then try automatic selection of file with consistent
           ! integration day and date among files named:
           !  <experiment name>_restphy_1.nc
           !  <experiment name>_restphy_2.nc
           !  <experiment name>_restphy_3.nc
           do i = 1, 4
              write(rstfnm, '(2a,i1,a)') trim(runid), '_restphy_', i, '.nc'
              inquire(file = rstfnm, exist = file_exist)
              if (file_exist) then
                 call ncfopn(rstfnm, 'r', ' ', 1, iotype)
                 call ncgeti('nday0', date_rest%day)
                 call ncgeti('nmonth0', date_rest%month)
                 call ncgeti('nyear0', date_rest%year)
                 call ncgetr('time0', time0)
                 call ncgetr('time', time)
                 errstat = daynum_diff(calendar, date_rest, date, dndiff)
                 if (errstat /= calendar_noerr) then
                    if (mnproc == 1) then
                       write(lp,*) 'restart_read: daynum_diff error: ', &
                                   trim(calendar_errstr(errstat))
                    endif
                    call xcstop('(restart_read)')
                           stop '(restart_read)'
                 endif
                 if (nint(time) == nday1 .and. &
                     dndiff == nint(time - time0)) exit
              endif
           enddo
           if (i > 3) then
              if (mnproc == 1) then
                 write(lp,*) 'restart_read: Could not find proper restart file!'
              endif
              call xcstop('(restart_read)')
                     stop '(restart_read)'
            endif
         endif

      endif

      if (mnproc == 1) then
         write(lp,*) 'restart_read: reading restart file '//trim(rstfnm)
      endif

      ! Compute extended masks at u-, v- and q-points.
      call extended_masks

      ! ------------------------------------------------------------------------
      ! Set unit conversion factors as needed. The ratio of bottom pressure to
      ! bottom geopotential should approximate sea water density that robustly
      ! can be tested if MKS or CGS units was used in the restart file. Convert
      ! units if needed.
      ! ------------------------------------------------------------------------

      call readfld('pb', no_unitconv, pb, ip)
      call readfld('phi', no_unitconv, phi(1-nbdy,1-nbdy,kk+1), ip)

      pb_max = 0._r8
      phi_min = 0._r8
   !$omp parallel do private(l, i) reduction(min:phi_min) reduction(max:pb_max)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            pb_max = max(pb_max, pb(i,j,1))
            phi_min = min(phi_min, phi(i,j,kk+1))
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xcmax(pb_max)
      call xcmin(phi_min)
      rho_restart = - pb_max/phi_min

      if (rho_restart > 1.e2_r8) then
         if (.not. use_MKS) then
            if (mnproc == 1) &
                 write(lp,*) 'restart_read: restart variables will be converted '// &
                             'from MKS to CGS units.'
            l_unitconv    = 1.e2_r8
            m_unitconv    = 1.e3_r8
            p_unitconv    = 1.e1_r8
            r_unitconv    = 1.e-3_r8
            l2_unitconv   = 1.e4_r8
            l3_unitconv   = 1.e6_r8
            lm_unitconv   = 1.e5_r8
            l2m2_unitconv = 1.e10_r8
            pi_unitconv   = 1.e-1_r8
            ri_unitconv   = 1.e3_r8
            l2i_unitconv  = 1.e-4_r8
            ml2i_unitconv = 1.e-1_r8
         end if
      else
         if (use_MKS) then
            if (mnproc == 1) &
                 write(lp,*) 'restart_read: restart variables will be converted '// &
                             'from CGS to MKS units.'
            l_unitconv    = 1.e-2_r8
            m_unitconv    = 1.e-3_r8
            p_unitconv    = 1.e-1_r8
            r_unitconv    = 1.e3_r8
            l2_unitconv   = 1.e-4_r8
            l3_unitconv   = 1.e-6_r8
            lm_unitconv   = 1.e-5_r8
            l2m2_unitconv = 1.e-10_r8
            pi_unitconv   = 1.e1_r8
            ri_unitconv   = 1.e-3_r8
            l2i_unitconv  = 1.e4_r8
            ml2i_unitconv = 1.e1_r8
         end if
      endif

      call readfld('u', l_unitconv, u, iuu)
      call readfld('v', l_unitconv, v, ivv)
      call readfld('dp', p_unitconv, dp, ip)
      call readfld('dpold', p_unitconv, dpold, ip)
      call readfld('temp', no_unitconv, temp, ip)
      call readfld('saln', no_unitconv, saln, ip)
      call readfld('sigma', r_unitconv, sigma, ip)
      call readfld('sigmar', r_unitconv, sigmar, ip)
      call readfld('pgfx', l2_unitconv, pgfx, iuu)
      call readfld('pgfy', l2_unitconv, pgfy, ivv)
      call readfld('pb', p_unitconv, pb, ip)
      call readfld('pb_mn', p_unitconv, pb_mn, ip)
      call readfld('pb_p', p_unitconv, pb_p, ip)
      call readfld('pbu', p_unitconv, pbu, iuu)
      call readfld('pbv', p_unitconv, pbv, ivv)
      call readfld('pbu_p', p_unitconv, pbu_p, iuu)
      call readfld('pbv_p', p_unitconv, pbv_p, ivv)
      call readfld('ub', l_unitconv, ub, iuu)
      call readfld('vb', l_unitconv, vb, ivv)
      call readfld('uflx', lm_unitconv, uflx, iuu)
      call readfld('utflx', lm_unitconv, utflx, iuu)
      call readfld('usflx', lm_unitconv, usflx, iuu)
      call readfld('umfltd', lm_unitconv, umfltd, iuu, required = .false.)
      call readfld('utfltd', lm_unitconv, utfltd, iuu, required = .false.)
      call readfld('utflld', lm_unitconv, utflld, iuu, required = .false.)
      call readfld('usfltd', lm_unitconv, usfltd, iuu, required = .false.)
      call readfld('usflld', lm_unitconv, usflld, iuu, required = .false.)
      call readfld('vflx', lm_unitconv, vflx, ivv)
      call readfld('vtflx', lm_unitconv, vtflx, ivv)
      call readfld('vsflx', lm_unitconv, vsflx, ivv)
      call readfld('vmfltd', lm_unitconv, vmfltd, ivv, required = .false.)
      call readfld('vtfltd', lm_unitconv, vtfltd, ivv, required = .false.)
      call readfld('vtflld', lm_unitconv, vtflld, ivv, required = .false.)
      call readfld('vsfltd', lm_unitconv, vsfltd, ivv, required = .false.)
      call readfld('vsflld', lm_unitconv, vsflld, ivv, required = .false.)

      call readfld('ubflx', lm_unitconv, ubflx, iuu)
      call readfld('vbflx', lm_unitconv, vbflx, ivv)
      call readfld('ubflx_mn', lm_unitconv, ubflx_mn, iuu)
      call readfld('vbflx_mn', lm_unitconv, vbflx_mn, ivv)
      call readfld('ubflxs', lm_unitconv, ubflxs, iuu)
      call readfld('vbflxs', lm_unitconv, vbflxs, ivv)
      call readfld('ubflxs_p', lm_unitconv, ubflxs_p, iuu)
      call readfld('vbflxs_p', lm_unitconv, vbflxs_p, ivv)
      call readfld('ubcors_p', l_unitconv, ubcors_p, iuu)
      call readfld('vbcors_p', l_unitconv, vbcors_p, ivv)
      call readfld('pvtrop', pi_unitconv, pvtrop, iqq)
      call readfld('pgfxm', l2_unitconv, pgfxm, iuu)
      call readfld('pgfym', l2_unitconv, pgfym, ivv)
      call readfld('xixp', ri_unitconv, xixp, iuu)
      call readfld('xixm', ri_unitconv, xixm, iuu)
      call readfld('xiyp', ri_unitconv, xiyp, ivv)
      call readfld('xiym', ri_unitconv, xiym, ivv)
      call readfld('phi', l2_unitconv, phi(1-nbdy,1-nbdy,kk+1), ip)
      call readfld('sealv', l_unitconv, sealv, ip)
      call readfld('ustar', l_unitconv, ustar, ip)
      call readfld('kfpla', no_unitconv, rkfpla, ip)
      call readfld('ficem', no_unitconv, ficem, ip)

      if (vcoord_tag == vcoord_isopyc_bulkml) then
         call readfld('buoyfl', l2_unitconv, buoyfl, ip)
         call readfld('uml', l_unitconv, uml, iuu, required = .false.)
         call readfld('vml', l_unitconv, vml, ivv, required = .false.)
         call readfld('umlres', l_unitconv, umlres, iuu, required = .false.)
         call readfld('vmlres', l_unitconv, vmlres, ivv, required = .false.)
      endif

      if (vcoord_tag /= vcoord_isopyc_bulkml) then
         call readfld('dpu', p_unitconv, dpu, iu)
         call readfld('dpv', p_unitconv, dpv, iv)
         call readfld('difiso', l2_unitconv, difiso, ip)
         call readfld('OBLdepth', no_unitconv, OBLdepth, ip)
         call readfld('t_ns_nonloc', no_unitconv, t_ns_nonloc, ip)
         call readfld('s_nb_nonloc', no_unitconv, s_nb_nonloc, ip)
         call readfld('mu_nonloc', no_unitconv, mu_nonloc, iu)
         call readfld('mv_nonloc', no_unitconv, mv_nonloc, iv)
         call readfld('umflsm', lm_unitconv, umflsm, iuu)
         call readfld('utflsm', lm_unitconv, utflsm, iuu)
         call readfld('usflsm', lm_unitconv, usflsm, iuu)
         call readfld('vmflsm', lm_unitconv, vmflsm, ivv)
         call readfld('vtflsm', lm_unitconv, vtflsm, ivv)
         call readfld('vsflsm', lm_unitconv, vsflsm, ivv)
         call readfld('wstar3', l3_unitconv, wstar3, ip, required = .false.)
         call readfld('hbl_tf', l_unitconv, hbl_tf, ip, required = .false.)
         call readfld('wpup_tf', l2_unitconv, wpup_tf, ip, required = .false.)
         call readfld('hml_tf1', l_unitconv, hml_tf1, ip, required = .false.)
         call readfld('hml_tf', l_unitconv, hml_tf, ip, required = .false.)
      endif

      if (sprfac) then
         call readfld('eiacc', no_unitconv, eiacc, ip, required = .false.)
         call readfld('pracc', no_unitconv, pracc, ip, required = .false., &
                      fld_read = fld_read)
         if (fld_read) then
            call ncgetr('prfac', prfac)
            call xcbcst(prfac)
         elseif (mnproc == 1) then
            write(lp,*) 'restart_read: warning: fields needed for '// &
                        'balancing fresh water budget are not read from '// &
                        'restart file and will be initialized.'
         endif
      endif

      if (expcnf == 'ben02clim' .or. expcnf == 'ben02syn') then

         call readfld('cd_d', no_unitconv, cd_d, ip)
         call readfld('ch_d', no_unitconv, ch_d, ip)
         call readfld('ce_d', no_unitconv, ce_d, ip)
         call readfld('wg2_d', no_unitconv, wg2_d, ip)
         call readfld('cd_m', no_unitconv, cd_m, ip)
         call readfld('ch_m', no_unitconv, ch_m, ip)
         call readfld('ce_m', no_unitconv, ce_m, ip)
         call readfld('wg2_m', no_unitconv, wg2_m, ip)
         call readfld('rhoa', no_unitconv, rhoa, ip)
         call readfld('tsi_tda', no_unitconv, tsi_tda, ip)
         call readfld('tml_tda', no_unitconv, tml_tda, ip)
         call readfld('sml_tda', no_unitconv, sml_tda, ip)
         call readfld('alb_tda', no_unitconv, alb_tda, ip)
         call readfld('fice_tda', no_unitconv, fice_tda, ip)

         call ncgeti('ntda', ntda)
         call xcbcst(ntda)

         call readfld('hicem', no_unitconv, hicem, ip)
         call readfld('tsrfm', no_unitconv, tsrfm, ip)
         call readfld('hsnwm', no_unitconv, hsnwm, ip)
         call readfld('ticem', no_unitconv, ticem, ip)
         call readfld('iagem', no_unitconv, iagem, ip)
         call readfld('rnfres', no_unitconv, rnfres, ip)

      endif

      if (expcnf == 'channel') then
         call ncgeti('ntda', ntda)
         call xcbcst(ntda)
!        call readfld('rnfres', no_unitconv, rnfres, ip)
      endif

      if (expcnf == 'cesm') then
         call readfld('ustarw_da', no_unitconv, ustarw_da, ip)
         call readfld('ztx_da', no_unitconv, ztx_da, ip, required = .false.)
         call readfld('mty_da', no_unitconv, mty_da, ip, required = .false.)
         call readfld('lip_da', no_unitconv, lip_da, ip, required = .false.)
         call readfld('sop_da', no_unitconv, sop_da, ip, required = .false.)
         call readfld('eva_da', no_unitconv, eva_da, ip, required = .false.)
         call readfld('rnf_da', no_unitconv, rnf_da, ip, required = .false.)
         call readfld('rfi_da', no_unitconv, rfi_da, ip, required = .false.)
         call readfld('fmltfz_da', no_unitconv, fmltfz_da, ip, &
                      required = .false.)
         call readfld('sfl_da', no_unitconv, sfl_da, ip, required = .false.)
         call readfld('swa_da', no_unitconv, swa_da, ip, required = .false.)
         call readfld('nsf_da', no_unitconv, nsf_da, ip, required = .false.)
         call readfld('hmlt_da', no_unitconv, hmlt_da, ip, required = .false.)
         call readfld('slp_da', no_unitconv, slp_da, ip, required = .false.)
         call readfld('ficem_da', no_unitconv, ficem_da, ip, required = .false.)
         call readfld('abswnd_da', no_unitconv, abswnd_da, ip, &
                      required = .false.)
         call readfld('atmco2_da', no_unitconv, atmco2_da, ip, &
                      required = .false., fld_read = fld_read)
         if (fld_read) then
            call ncgeti('l2ci', l2ci)
            call xcbcst(l2ci)
            l1ci = 3 - l2ci
         else
            if (mnproc == 1) then
               write(lp,*) 'restart_read: warning: time levels for '// &
                           'interpolation of forcing fields is not read '// &
                           'from restart file.'
            endif
            l1ci = 1
            l2ci = 1
         endif
         call readfld('frzpot', no_unitconv, frzpot, ip)
         call readfld('mltpot', no_unitconv, mltpot, ip)
         call readfld('flxco2', no_unitconv, flxco2, ip, &
                      required = .false., fld_read = fld_read)
         if (.not.fld_read .and. mnproc == 1) &
            write(lp,*) 'restart_read: warning: air-sea CO2 flux is not '// &
                        'read from restart file and will be initialized '// &
                        'to zero.'
         call readfld('flxdms', no_unitconv, flxdms, ip, &
                      required = .false., fld_read = fld_read)
         if (.not.fld_read .and. mnproc == 1) &
            write(lp,*) 'restart_read: warning: DMS flux is not read from '// &
                        'restart file and will be initialized to zero.'
         call readfld('flxbrf', no_unitconv, flxbrf, ip, &
                      required = .false., fld_read = fld_read)
         if (.not.fld_read .and. mnproc == 1) &
            write(lp,*) 'restart_read: warning: bromoform flux is not read '// &
                        'from restart file and will be initialized to zero.'
         call readfld('flxn2o', no_unitconv, flxn2o, ip, &
                      required = .false., fld_read = fld_read)
         if (.not.fld_read .and. mnproc == 1) &
            write(lp,*) 'restart_read: warning: N2O flux is not read '// &
                        'from restart file and will be initialized to zero.'
         call readfld('flxnh3', no_unitconv, flxnh3, ip, &
                      required = .false., fld_read = fld_read)
         if (.not.fld_read .and. mnproc == 1) &
            write(lp,*) 'restart_read: warning: NH3 flux is not read '// &
                        'from restart file and will be initialized to zero.'
      endif

      if (use_TRC) then
         if (use_TKE) then
            call readfld('tke', l2_unitconv, trc(1-nbdy,1-nbdy,1,itrtke), ip, required=.false.)
            call readfld('gls_psi', l2_unitconv, trc(1-nbdy,1-nbdy,1,itrgls), ip, required=.false.)
            call readfld('L_scale', l_unitconv, L_scale, ip, required = .false.)
            call readfld('difdia', l2_unitconv, difdia, ip, required = .false.)
            call readfld('ustarb', l_unitconv, ustarb, ip, required = .false.)
         end if
         if (use_IDLAGE) then
            call readfld('idlage', no_unitconv, trc(1-nbdy,1-nbdy,1,itriag), ip, &
                 required=.false., fld_read=fld_read)
            if (.not.fld_read) then
               if (mnproc == 1) then
                  write(lp,*) 'restart_read: warning: ideal age tracer is not '// &
                       'read from restart file and will be initialized to zero'
               end if
               call idlage_init
            endif
         end if
      end if

      ! Read accumulated fields.
      if (nstep1 > nstep0) then
         do n = 1, nphy
            write(c2, '(i2.2)') n
            if (ncinqa('nacc_phy'//c2)) then
               call ncgeti('nacc_phy'//c2, nacc_phy(n))
            else
               nacc_phy(n) = 0
            endif
            call xcbcst(nacc_phy(n))
            if (nacc_phy(n) > 0) then
               if (ACC_UB(n) /= 0) &
                  call readfld('ub_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_UB(n)), iuu)
               if (ACC_VB(n) /= 0) &
                  call readfld('vb_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_VB(n)), ivv)
               if (ACC_UBFLXS(n) /= 0) &
                  call readfld('ubflxs_phy'//c2, lm_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(n)), iuu)
               if (ACC_VBFLXS(n) /= 0) &
                  call readfld('vbflxs_phy'//c2, lm_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(n)), ivv)
               if (ACC_ZTX(n) /= 0) &
                  call readfld('ztx_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_ZTX(n)), iuu)
               if (ACC_MTY(n) /= 0) &
                  call readfld('mty_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTY(n)), ivv)
               if (ACC_TAUX(n) /= 0) &
                  call readfld('taux_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_TAUX(n)), iuu)
               if (ACC_TAUY(n) /= 0) &
                  call readfld('tauy_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_TAUY(n)), ivv)
               if (ACC_UICE(n) /= 0) &
                  call readfld('uice_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_UICE(n)), iuu)
               if (ACC_VICE(n) /= 0) &
                  call readfld('vice_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_VICE(n)), ivv)
               if (ACC_IVOLU(n) /= 0) &
                  call readfld('ivolu_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_IVOLU(n)), iuu)
               if (ACC_IVOLV(n) /= 0) &
                  call readfld('ivolv_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_IVOLV(n)), ivv)
               if (ACC_PSRF(n) /= 0) &
                  call readfld('psrf_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_PSRF(n)), ip)
               if (ACC_PBOT(n) /= 0) &
                  call readfld('pbot_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_PBOT(n)), ip)
               if (ACC_SEALV(n) /= 0) &
                  call readfld('sealv_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SEALV(n)), ip)
               if (ACC_SLVSQ(n) /= 0) &
                  call readfld('slvsq_phy'//c2, l2_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SLVSQ(n)), ip)
               if (ACC_SSS(n) /= 0) &
                  call readfld('sss_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SSS(n)), ip)
               if (ACC_SSSSQ(n) /= 0) &
                  call readfld('ssssq_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SSSSQ(n)), ip)
               if (ACC_SBOT(n) /= 0) &
                  call readfld('sbot_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SBOT(n)), ip)
               if (ACC_SST(n) /= 0) &
                  call readfld('sst_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SST(n)), ip)
               if (ACC_SSTSQ(n) /= 0) &
                  call readfld('sstsq_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SSTSQ(n)), ip)
               if (ACC_TBOT(n) /= 0) &
                  call readfld('tbot_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_TBOT(n)), ip)
               if (ACC_SIGMX(n) /= 0) &
                  call readfld('sigmx_phy'//c2, r_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SIGMX(n)), ip)
               if (ACC_MLD(n) /= 0) &
                  call readfld('mld_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MLD(n)), ip)
               if (ACC_MAXMLD(n) /= 0) &
                  call readfld('maxmld_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MAXMLD(n)), ip)
               if (ACC_MLTS(n) /= 0) &
                  call readfld('mlts_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MLTS(n)), ip)
               if (ACC_MLTSMN(n) /= 0) &
                  call readfld('mltsmn_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MLTSMN(n)), ip)
               if (ACC_MLTSMX(n) /= 0) &
                  call readfld('mltsmx_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MLTSMX(n)), ip)
               if (ACC_MLTSSQ(n) /= 0) &
                  call readfld('mltssq_phy'//c2, l2_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MLTSSQ(n)), ip)
               if (ACC_T20D(n) /= 0) &
                  call readfld('t20d_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_T20D(n)), ip)
               if (ACC_ALB(n) /= 0) &
                  call readfld('alb_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_ALB(n)), ip)
               if (ACC_SWA(n) /= 0) &
                  call readfld('swa_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SWA(n)), ip)
               if (ACC_NSF(n) /= 0) &
                  call readfld('nsf_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_NSF(n)), ip)
               if (ACC_DFL(n) /= 0) &
                  call readfld('dfl_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_DFL(n)), ip)
               if (ACC_LIP(n) /= 0) &
                  call readfld('lip_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_LIP(n)), ip)
               if (ACC_SOP(n) /= 0) &
                  call readfld('sop_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SOP(n)), ip)
               if (ACC_EVA(n) /= 0) &
                  call readfld('eva_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_EVA(n)), ip)
               if (ACC_RNFFLX(n) /= 0) &
                  call readfld('rnfflx_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_RNFFLX(n)), ip)
               if (ACC_RFIFLX(n) /= 0) &
                  call readfld('rfiflx_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_RFIFLX(n)), ip)
               if (ACC_SFL(n) /= 0) &
                  call readfld('sfl_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SFL(n)), ip)
               if (ACC_BRNFLX(n) /= 0) &
                  call readfld('brnflx_phy'//c2, ml2i_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_BRNFLX(n)), ip)
               if (ACC_BRNPD(n) /= 0) &
                  call readfld('brnpd_phy'//c2, p_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_BRNPD(n)), ip)
               if (ACC_SURFLX(n) /= 0) &
                  call readfld('surflx_phy'//c2, l2i_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SURFLX(n)), ip)
               if (ACC_SURRLX(n) /= 0) &
                  call readfld('surrlx_phy'//c2, l2i_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SURRLX(n)), ip)
               if (ACC_SALFLX(n) /= 0) &
                  call readfld('salflx_phy'//c2, ml2i_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SALFLX(n)), ip)
               if (ACC_SALRLX(n) /= 0) &
                  call readfld('salrlx_phy'//c2, ml2i_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_SALRLX(n)), ip)
               if (ACC_ABSWND(n) /= 0) &
                  call readfld('abswnd_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_ABSWND(n)), ip)
               if (ACC_USTAR(n) /= 0) &
                  call readfld('ustar_phy'//c2, l_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_USTAR(n)), ip)
               if (ACC_USTAR3(n) /= 0) &
                  call readfld('ustar3_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_USTAR3(n)), ip)
               if (ACC_IDKEDT(n) /= 0) &
                  call readfld('idkedt_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_IDKEDT(n)), ip)
               if (ACC_MTKEUS(n) /= 0) &
                  call readfld('mtkeus_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTKEUS(n)), ip)
               if (ACC_MTKENI(n) /= 0) &
                  call readfld('mtkeni_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTKENI(n)), ip)
               if (ACC_MTKEBF(n) /= 0) &
                  call readfld('mtkebf_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTKEBF(n)), ip)
               if (ACC_MTKERS(n) /= 0) &
                  call readfld('mtkers_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTKERS(n)), ip)
               if (ACC_MTKEPE(n) /= 0) &
                  call readfld('mtkepe_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTKEPE(n)), ip)
               if (ACC_MTKEKE(n) /= 0) &
                  call readfld('mtkeke_phy'//c2, l3_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_MTKEKE(n)), ip)
               if (ACC_LAMULT(n) /= 0) &
                  call readfld('lamult_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_LAMULT(n)), ip)
               if (ACC_LASL(n) /= 0) &
                  call readfld('lasl_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_LASL(n)), ip)
               if (ACC_USTOKES(n) /= 0) &
                  call readfld('ustokes_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_USTOKES(n)), iuu)
               if (ACC_VSTOKES(n) /= 0) &
                  call readfld('vstokes_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_VSTOKES(n)), ivv)
               if (ACC_FMLTFZ(n) /= 0) &
                  call readfld('fmltfz_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_FMLTFZ(n)), ip)
               if (ACC_HMLTFZ(n) /= 0) &
                  call readfld('hmltfz_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_HMLTFZ(n)), ip)
               if (ACC_HICE(n) /= 0) &
                  call readfld('hice_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_HICE(n)), ip)
               if (ACC_HSNW(n) /= 0) &
                  call readfld('hsnw_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_HSNW(n)), ip)
               if (ACC_FICE(n) /= 0) &
                  call readfld('fice_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_FICE(n)), ip)
               if (ACC_TSRF(n) /= 0) &
                  call readfld('tsrf_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_TSRF(n)), ip)
               if (ACC_TICE(n) /= 0) &
                  call readfld('tice_phy'//c2, no_unitconv, &
                               phyh2d(1-nbdy,1-nbdy,ACC_TICE(n)), ip)
               if (ACC_UVEL(n) /= 0) &
                  call readfld('uvel_phy'//c2, m_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UVEL(n)), iuu)
               if (ACC_VVEL(n) /= 0) &
                  call readfld('vvel_phy'//c2, m_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VVEL(n)), ivv)
               if (ACC_DPU(n) /= 0) &
                  call readfld('dpu_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DPU(n)), iuu)
               if (ACC_DPV(n) /= 0) &
                  call readfld('dpv_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DPV(n)), ivv)
               if (ACC_UFLX(n) /= 0) &
                  call readfld('uflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UFLX(n)), iuu)
               if (ACC_VFLX(n) /= 0) &
                  call readfld('vflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VFLX(n)), ivv)
               if (ACC_UTFLX(n) /= 0) &
                  call readfld('utflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLX(n)), iuu)
               if (ACC_VTFLX(n) /= 0) &
                  call readfld('vtflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLX(n)), ivv)
               if (ACC_USFLX(n) /= 0) &
                  call readfld('usflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_USFLX(n)), iuu)
               if (ACC_VSFLX(n) /= 0) &
                  call readfld('vsflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLX(n)), ivv)
               if (ACC_UMFLTD(n) /= 0) &
                  call readfld('umfltd_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UMFLTD(n)), iuu)
               if (ACC_VMFLTD(n) /= 0) &
                  call readfld('vmfltd_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VMFLTD(n)), ivv)
               if (ACC_UMFLSM(n) /= 0) &
                  call readfld('umflsm_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UMFLSM(n)), iuu)
               if (ACC_VMFLSM(n) /= 0) &
                  call readfld('vmflsm_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VMFLSM(n)), ivv)
               if (ACC_UTFLTD(n) /= 0) &
                  call readfld('utfltd_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLTD(n)), iuu)
               if (ACC_VTFLTD(n) /= 0) &
                  call readfld('vtfltd_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLTD(n)), ivv)
               if (ACC_UTFLSM(n) /= 0) &
                  call readfld('utflsm_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLSM(n)), iuu)
               if (ACC_VTFLSM(n) /= 0) &
                  call readfld('vtflsm_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLSM(n)), ivv)
               if (ACC_UTFLLD(n) /= 0) &
                  call readfld('utflld_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_UTFLLD(n)), iuu)
               if (ACC_VTFLLD(n) /= 0) &
                  call readfld('vtflld_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VTFLLD(n)), ivv)
               if (ACC_USFLTD(n) /= 0) &
                  call readfld('usfltd_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_USFLTD(n)), iuu)
               if (ACC_VSFLTD(n) /= 0) &
                  call readfld('vsfltd_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLTD(n)), ivv)
               if (ACC_USFLSM(n) /= 0) &
                  call readfld('usflsm_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_USFLSM(n)), iuu)
               if (ACC_VSFLSM(n) /= 0) &
                  call readfld('vsflsm_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLSM(n)), ivv)
               if (ACC_USFLLD(n) /= 0) &
                  call readfld('usflld_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_USFLLD(n)), iuu)
               if (ACC_VSFLLD(n) /= 0) &
                  call readfld('vsflld_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_VSFLLD(n)), ivv)
               if (ACC_SALN(n) /= 0) &
                  call readfld('saln_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_SALN(n)), ip)
               if (ACC_TEMP(n) /= 0) &
                  call readfld('temp_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_TEMP(n)), ip)
               if (ACC_DP(n) /= 0) &
                  call readfld('dp_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DP(n)), ip)
               if (ACC_DZ(n) /= 0) &
                  call readfld('dz_phy'//c2, l_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DZ(n)), ip)
               if (ACC_BFSQ(n) /= 0) &
                  call readfld('bfsq_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_BFSQ(n)), ip)
               if (ACC_DIFDIA(n) /= 0) &
                  call readfld('difdia_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DIFDIA(n)), ip)
               if (ACC_DIFVMO(n) /= 0) &
                  call readfld('difvmo_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DIFVMO(n)), ip)
               if (ACC_DIFVHO(n) /= 0) &
                  call readfld('difvho_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DIFVHO(n)), ip)
               if (ACC_DIFVSO(n) /= 0) &
                  call readfld('difvso_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DIFVSO(n)), ip)
               if (ACC_DIFINT(n) /= 0) &
                  call readfld('difint_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DIFINT(n)), ip)
               if (ACC_DIFISO(n) /= 0) &
                  call readfld('difiso_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DIFISO(n)), ip)
               if (ACC_WFLX(n) /= 0) &
                  call readfld('wflx_phy'//c2, lm_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_WFLX(n)), ip)
               if (ACC_WFLX2(n) /= 0) &
                  call readfld('wflx2_phy'//c2, l2m2_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_WFLX2(n)), ip)
               if (ACC_AVDSG(n) /= 0) &
                  call readfld('avdsg_phy'//c2, r_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_AVDSG(n)), ip)
               if (ACC_DPVOR(n) /= 0) &
                  call readfld('dpvor_phy'//c2, p_unitconv, &
                               phylyr(1-nbdy,1-nbdy,1,ACC_DPVOR(n)), ip)
               if (use_TRC .and. use_TKE) then
                  if (ACC_TKE(n) /= 0) &
                       call readfld('tke_phy'//c2, lm_unitconv, &
                                    phylyr(1-nbdy,1-nbdy,1,ACC_TKE(n)), ip)
                  if (ACC_GLS_PSI(n) /= 0) &
                       call readfld('gls_psi_phy'//c2, lm_unitconv, &
                                    phylyr(1-nbdy,1-nbdy,1,ACC_GLS_PSI(n)), ip)
               end if
               if (ACC_UVELLVL(n) /= 0) &
                  call readfld('uvellvl_phy'//c2, l_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UVELLVL(n)), iuu)
               if (ACC_VVELLVL(n) /= 0) &
                  call readfld('vvellvl_phy'//c2, l_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VVELLVL(n)), ivv)
               if (ACC_UFLXLVL(n) /= 0) then
                  call readfld('uflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)), iuu)
                  call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(n)), &
                              1, ddm, 1, 1, halo_uv)
               endif
               if (ACC_VFLXLVL(n) /= 0) then
                  call readfld('vflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)), ivv)
                  call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(n)), &
                              1, ddm, 1, 1, halo_vv)
               endif
               if (ACC_UTFLXLVL(n) /= 0) &
                  call readfld('utflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLXLVL(n)), iuu)
               if (ACC_VTFLXLVL(n) /= 0) &
                  call readfld('vtflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLXLVL(n)), ivv)
               if (ACC_USFLXLVL(n) /= 0) &
                  call readfld('usflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_USFLXLVL(n)), iuu)
               if (ACC_VSFLXLVL(n) /= 0) &
                  call readfld('vsflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLXLVL(n)), ivv)
               if (ACC_UMFLTDLVL(n) /= 0) &
                  call readfld('umfltdlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UMFLTDLVL(n)), iuu)
               if (ACC_VMFLTDLVL(n) /= 0) &
                  call readfld('vmfltdlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VMFLTDLVL(n)), ivv)
               if (ACC_UMFLSMLVL(n) /= 0) &
                  call readfld('umflsmlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UMFLSMLVL(n)), iuu)
               if (ACC_VMFLSMLVL(n) /= 0) &
                  call readfld('vmflsmlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VMFLSMLVL(n)), ivv)
               if (ACC_UTFLTDLVL(n) /= 0) &
                  call readfld('utfltdlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLTDLVL(n)), iuu)
               if (ACC_VTFLTDLVL(n) /= 0) &
                  call readfld('vtfltdlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLTDLVL(n)), ivv)
               if (ACC_UTFLSMLVL(n) /= 0) &
                  call readfld('utflsmlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLSMLVL(n)), iuu)
               if (ACC_VTFLSMLVL(n) /= 0) &
                  call readfld('vtflsmlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLSMLVL(n)), ivv)
               if (ACC_UTFLLDLVL(n) /= 0) &
                  call readfld('utflldlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_UTFLLDLVL(n)), iuu)
               if (ACC_VTFLLDLVL(n) /= 0) &
                  call readfld('vtflldlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VTFLLDLVL(n)), ivv)
               if (ACC_USFLTDLVL(n) /= 0) &
                  call readfld('usfltdlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_USFLTDLVL(n)), iuu)
               if (ACC_VSFLTDLVL(n) /= 0) &
                  call readfld('vsfltdlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLTDLVL(n)), ivv)
               if (ACC_USFLSMLVL(n) /= 0) &
                  call readfld('usflsmlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_USFLSMLVL(n)), iuu)
               if (ACC_VSFLSMLVL(n) /= 0) &
                  call readfld('vsflsmlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLSMLVL(n)), ivv)
               if (ACC_USFLLDLVL(n) /= 0) &
                  call readfld('usflldlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_USFLLDLVL(n)), iuu)
               if (ACC_VSFLLDLVL(n) /= 0) &
                  call readfld('vsflldlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_VSFLLDLVL(n)), ivv)
               if (ACC_SALNLVL(n) /= 0) &
                  call readfld('salnlvl_phy'//c2, no_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_SALNLVL(n)), ip)
               if (ACC_TEMPLVL(n) /= 0) &
                  call readfld('templvl_phy'//c2, no_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_TEMPLVL(n)), ip)
               if (ACC_DZLVL(n) /= 0) &
                  call readfld('dzlvl_phy'//c2, l_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DZLVL(n)), ip)
               if (ACC_BFSQLVL(n) /= 0) &
                  call readfld('bfsqlvl_phy'//c2, no_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_BFSQLVL(n)), ip)
               if (ACC_DIFDIALVL(n) /= 0) &
                  call readfld('difdialvl_phy'//c2, l2_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DIFDIALVL(n)), ip)
               if (ACC_DIFVMOLVL(n) /= 0) &
                  call readfld('difvmolvl_phy'//c2, l2_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DIFVMOLVL(n)), ip)
               if (ACC_DIFVHOLVL(n) /= 0) &
                  call readfld('difvholvl_phy'//c2, l2_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DIFVHOLVL(n)), ip)
               if (ACC_DIFVSOLVL(n) /= 0) &
                  call readfld('difvsolvl_phy'//c2, l2_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DIFVSOLVL(n)), ip)
               if (ACC_DIFINTLVL(n) /= 0) &
                  call readfld('difintlvl_phy'//c2, l2_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DIFINTLVL(n)), ip)
               if (ACC_DIFISOLVL(n) /= 0) &
                  call readfld('difisolvl_phy'//c2, l2_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_DIFISOLVL(n)), ip)
               if (ACC_WFLXLVL(n) /= 0) &
                  call readfld('wflxlvl_phy'//c2, lm_unitconv, &
                               phylvl(1-nbdy,1-nbdy,1,ACC_WFLXLVL(n)), ip)
               if (ACC_WFLX2LVL(n) /= 0) &
                  call readfld('wflx2lvl_phy'//c2, l2m2_unitconv, &
                                phylvl(1-nbdy,1-nbdy,1,ACC_WFLX2LVL(n)), ip)
               if (ACC_PVLVL(n) /= 0) &
                  call readfld('pvlvl_phy'//c2, l2i_unitconv, &
                                phylvl(1-nbdy,1-nbdy,1,ACC_PVLVL(n)), ip)
               if (use_TRC .and. use_TKE) then
                  if (ACC_TKELVL(n) /= 0) &
                       call readfld('tkelvl_phy'//c2, l2_unitconv, &
                                     phylvl(1-nbdy,1-nbdy,1,ACC_TKELVL(n)), ip)
                  if (ACC_GLS_PSILVL(n) /= 0) &
                       call readfld('gls_psilvl_phy'//c2, l2_unitconv, &
                                    phylvl(1-nbdy,1-nbdy,1,ACC_GLS_PSILVL(n)), ip)
               end if
            endif
         enddo
      endif

      call ncfcls

   !$omp parallel do private(i)
      do j = 1, jj
         do i = 1, ii
            if (ip(i,j) == 1) then
               kfpla(i,j,1) = nint(rkfpla(i,j,1))
               kfpla(i,j,2) = nint(rkfpla(i,j,2))
            endif
         enddo
      enddo
   !$omp end parallel do

      ! ------------------------------------------------------------------------
      ! Set minimum physical temperature for each isopycnic layer.
      ! ------------------------------------------------------------------------

      call settemmin

      if (use_TRC) then
         if (.not.resume_flag) call restart_trcrd(rstfnm)
      end if

      if (ditflx) then

         ! Read diag. heat flux restart file if available.
         if (expcnf == 'cesm') then
            fnm = trim(runid)//'.blom.rtflx.'//rstfnm(len_trim(runid) + 10:)
         else
            fnm = trim(runid)//'_resttflx_'//rstfnm(len_trim(runid) + 10:)
         endif
         if (mnproc  == 1) inquire(file = fnm, exist = file_exist)
         call xcbcst(file_exist)
         if (file_exist) then
            call ncfopn(fnm, 'r', ' ', 1, iotype)
            if (mnproc == 1) &
               write(lp,*) &
                  'restart_read: reading diag. heat flux restart file '// &
                  trim(fnm)
            call ncgetr('time', time)
            if (nint(time) /= nday1 .and. mnproc == 1) then
               write(lp, '(a,i6.6,a)') &
                  ' Integration day ', nint(time), &
                  ' in diag. heat flux restart file differs from'
               write(lp, '(a,i6.6,a)') ' start day ', nday1, ' in limits file!'
               call ncfcls
               call xcstop('(restart_read)')
                      stop '(restart_read)'
            endif
            call readfld('tflxdi', l2i_unitconv, tflxdi, ip)
            call ncgeti('nflxdi', nflxdi)
         else
            if (mnproc == 1) &
               write(lp,*) &
                  'restart_read: warning: No diag. heat flux restart file found'
         endif

      endif

      if (disflx) then

         ! Read diag. salt flux restart file if available.
         if (expcnf == 'cesm') then
            fnm = trim(runid)//'.blom.rsflx.'//rstfnm(len_trim(runid) + 10:)
         else
            fnm = trim(runid)//'_restsflx_'//rstfnm(len_trim(runid) + 10:)
         endif
         if (mnproc == 1) inquire(file = fnm, exist = file_exist)
         call xcbcst(file_exist)
         if (file_exist) then
            call ncfopn(fnm, 'r', ' ', 1, iotype)
            if (mnproc == 1) &
               write(lp,*) &
                  'restart_read: reading diag. salt flux restart file '// &
                  trim(fnm)
            call ncgetr('time', time)
            if (nint(time) /= nday1 .and. mnproc == 1) then
               write(lp, '(a,i6.6,a)') &
                  ' Integration day ', nint(time), &
                  ' in diag. salt flux restart file differs from'
               write(lp, '(a,i6.6,a)') ' start day ', nday1, ' in limits file!'
               call ncfcls
               call xcstop('(restart_read)')
                      stop '(restart_read)'
            endif
            call readfld('sflxdi', l2i_unitconv, sflxdi, ip)
            call ncgeti('nflxdi', nflxdi)
            call ncfcls
         else
            if (mnproc == 1) &
              write(lp,*) &
                 'restart_read: warning: No diag. salt flux restart file found'
         endif

      endif

   end subroutine restart_read

end module mod_restart
