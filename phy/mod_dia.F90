! ------------------------------------------------------------------------------
! Copyright (C) 2010-2025 Ingo Bethke, Mats Bentsen, Mehmet Ilicak,
!                         Alok Kumar Gupta, Jörg Schwinger, Ping-Gin Chi,
!                         Mariana Vertenstein
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

module mod_dia

  use dimensions,    only: kdm, nreg
  use mod_types,     only: i2
  use mod_config,    only: expcnf, runid, inst_suffix
  use mod_calendar,  only: date_type, date_offset, calendar_noerr
  use mod_time,      only: date0, date, calendar, nstep, nstep_in_day, &
                           nday_of_year, time, time0, baclin, dlt
  use mod_constants, only: grav, spcifh, t0deg, alpha0, epsilp, spval, &
                           onem, onecm, onemm, g2kg
  use mod_xc,        only: xcstop, xchalt, xcaget, xcbcst, xcaput, &
                           xceget, xctilr, xcsum, &
                           isp, ifp, ilp, isu, ifu, ilu, isv, ifv, ilv, &
                           isu, ifu, ilu, isv, ip, halo_ps, &
                           iu, iv, ips, halo_qs, halo_uv, halo_vv
  use mod_nctools
  use netcdf,        only: nf90_fill_double
  use mod_vcoord,    only: vcoord_type, vcoord_tag, vcoord_isopyc_bulkml, &
                           sigref_spec, sigmar, sigref, &
                           sigref_fun_spec, sigref_adaption
  use mod_grid,      only: scp2, depths, area
  use mod_eos,       only: rho, p_alpha
  use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, &
                           uflx, vflx, utflx, vtflx, usflx, vsflx, &
                           p, phi, ubflxs, vbflxs, ub, vb, sealv
  use mod_momtum,    only: absvor, dpvor
  use mod_tmsmt,     only: dpold
  use mod_mxlayr,    only: mtkeus, mtkeni, mtkebf, mtkers, mtkepe, &
                           mtkeke, pbrnda
  use mod_difest,    only: OBLdepth
  use mod_diffusion, only: difint, difiso, difdia, &
                           Kvisc_m, Kdiff_t, Kdiff_s, &
                           umfltd, vmfltd, umflsm, vmflsm, &
                           utfltd, vtfltd, utflsm, vtflsm, &
                           utflld, vtflld, usfltd, vsfltd, &
                           usflsm, vsflsm, usflld, vsflld
  use mod_cmnfld,    only: z, bfsql, dz, mlts
  use mod_seaice,    only: ficem, hicem, hsnwm, uicem, vicem, iagem
  use mod_forcing,   only: swa, nsf, hmat, hmltfz, lip, sop, eva, rnf, rfi, &
                           fmltfz, sfl, ztx, mty, abswnd, &
                           lamult, lasl, ustokes, vstokes, surflx, &
                           surrlx, salflx, brnflx, salrlx, taux, tauy, &
                           ustar, ustar3
  use mod_niw,       only: idkedt
  use mod_utility,   only: util1, util2, util3, util4, fnmlen
  use mod_ben02,     only: dfl, alb
  use mod_thdysi,    only: tsrfm,ticem
  use mod_tracers,   only: ntrocn, ntr, natr, itriag, itrtke, itrgls, trc
  use mod_ifdefs,    only: use_TRC, use_TKE, use_ATRC, use_IDLAGE
  implicit none
  private

  !---------------------------------------------------------------
  ! module variables related  to the accumulation and averaging of
  ! diagnostic fields
  !---------------------------------------------------------------

  ! Averaging and writing frequencies for diagnostics output
  integer                    , public :: nphy
  integer, parameter         , public :: nphymax = 10
  integer, dimension(nphymax), public :: nacc_phy
  real   , dimension(nphymax), public :: diagfq_phy
  logical, dimension(nphymax), public :: diagmon_phy
  logical, dimension(nphymax), public :: diagann_phy
  real   , dimension(nphymax), public :: filefq_phy
  logical, dimension(nphymax), public :: filemon_phy
  logical, dimension(nphymax), public :: fileann_phy
  integer, dimension(nphymax), public :: alarm_phy

  ! Restart parameters
  real,    public :: rstfrq
  integer, public :: iotype
  integer, public :: rstcmp,rstfmt
  logical, public :: rstmon,rstann

  ! Copies of BLOM variables that are used for HAMOCC diagnostics
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), public :: pbath,ubath,vbath
  integer, public :: nstepinday

  ! 2d and 3d diagnostic variables
  integer :: nphyh2d,nphylyr,nphylvl
  real, allocatable, dimension(:,:,:)  , public :: phyh2d
  real, allocatable, dimension(:,:,:,:), public :: phylyr
  real, allocatable, dimension(:,:,:,:), public :: phylvl

  ! Levitus levels
#ifndef LEVITUS2X
  integer, parameter, public :: ddm = 35
  integer, parameter :: k350 = 12
  real   , parameter :: w350 = 1.
  real,  parameter, dimension(ddm), public :: depthslev = (/ &
       0000.0,0010.0,0020.0,0030.0,0050.0,0075.0,0100.0,0125.0,0150.0, &
       0200.0,0250.0,0300.0,0400.0,0500.0,0600.0,0700.0,0800.0,0900.0, &
       1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1750.0,2000.0,2500.0, &
       3000.0,3500.0,4000.0,4500.0,5000.0,5500.0,6000.0,6500.0/)
  real,  parameter, dimension(2,ddm), public :: &
       depthslev_bnds = reshape((/ &
       0000.0,0005.0,0005.0,0015.0,0015.0,0025.0,0025.0,0040.0,0040.0, &
       0062.5,0062.5,0087.5,0087.5,0112.5,0112.5,0137.5,0137.5,0175.0, &
       0175.0,0225.0,0225.0,0275.0,0275.0,0350.0,0350.0,0450.0,0450.0, &
       0550.0,0550.0,0650.0,0650.0,0750.0,0750.0,0850.0,0850.0,0950.0, &
       0950.0,1050.0,1050.0,1150.0,1150.0,1250.0,1250.0,1350.0,1350.0, &
       1450.0,1450.0,1625.0,1625.0,1875.0,1875.0,2250.0,2250.0,2750.0, &
       2750.0,3250.0,3250.0,3750.0,3750.0,4250.0,4250.0,4750.0,4750.0, &
       5250.0,5250.0,5750.0,5750.0,6250.0,6250.0,8000.0/),(/2,ddm/))
#else
  integer, parameter, public :: ddm=70
  integer, parameter :: k350 = 25
  real   , parameter :: w350 = .5
  real,  parameter, dimension(ddm), public :: depthslev = (/ &
       0000.0,0005.0,0010.0,0015.0,0020.0,0025.0,0030.0,0040.0,0050.0, &
       0062.5,0075.0,0087.5,0100.0,0112.5,0125.0,0137.5,0150.0,0175.0, &
       0200.0,0225.0,0250.0,0275.0,0300.0,0350.0,0400.0,0450.0,0500.0, &
       0550.0,0600.0,0650.0,0700.0,0750.0,0800.0,0850.0,0900.0,0950.0, &
       1000.0,1050.0,1100.0,1150.0,1200.0,1250.0,1300.0,1350.0,1400.0, &
       1450.0,1500.0,1625.0,1750.0,1875.0,2000.0,2250.0,2500.0,2750.0, &
       3000.0,3250.0,3500.0,3750.0,4000.0,4250.0,4500.0,4750.0,5000.0, &
       5250.0,5500.0,5750.0,6000.0,6250.0,6500.0,6750.0/)
  real,  parameter, dimension(2,ddm), public :: &
       depthslev_bnds = reshape((/ &
       0000.0,0002.5,0002.5,0007.5,0007.5,0012.5,0012.5,0017.5,0017.5, &
       0022.5,0022.5,0027.5,0027.5,0035.0,0035.0,0045.0,0045.0,0056.2, &
       0056.2,0068.8,0068.8,0081.2,0081.2,0093.8,0093.8,0106.2,0106.2, &
       0118.8,0118.8,0131.2,0131.2,0143.8,0143.8,0162.5,0162.5,0187.5, &
       0187.5,0212.5,0212.5,0237.5,0237.5,0262.5,0262.5,0287.5,0287.5, &
       0325.0,0325.0,0375.0,0375.0,0425.0,0425.0,0475.0,0475.0,0525.0, &
       0525.0,0575.0,0575.0,0625.0,0625.0,0675.0,0675.0,0725.0,0725.0, &
       0775.0,0775.0,0825.0,0825.0,0875.0,0875.0,0925.0,0925.0,0975.0, &
       0975.0,1025.0,1025.0,1075.0,1075.0,1125.0,1125.0,1175.0,1175.0, &
       1225.0,1225.0,1275.0,1275.0,1325.0,1325.0,1375.0,1375.0,1425.0, &
       1425.0,1475.0,1475.0,1562.5,1562.5,1687.5,1687.5,1812.5,1812.5, &
       1937.5,1937.5,2125.0,2125.0,2375.0,2375.0,2625.0,2625.0,2875.0, &
       2875.0,3125.0,3125.0,3375.0,3375.0,3625.0,3625.0,3875.0,3875.0, &
       4125.0,4125.0,4375.0,4375.0,4625.0,4625.0,4875.0,4875.0,5125.0, &
       5125.0,5375.0,5375.0,5625.0,5625.0,5875.0,5875.0,6125.0,6125.0, &
       6375.0,6375.0,6625.0,6625.0,8000.0/),(/2,ddm/))
#endif

  ! Meridional overturning and flux diagnostics
  integer, parameter         :: ldm=itdm+jtdm
  integer, parameter         :: sdm=ldm
  integer, parameter, public :: odm=10
  integer, parameter         :: slenmax=50
  integer, parameter, public :: rflgdm=20
  character(len=slenmax), dimension(odm), public :: mer_regnam = ''
  character(len=fnmlen), public :: mer_orfile
  character(len=fnmlen), public :: mer_mifile
  integer, dimension(odm,rflgdm), public :: mer_regflg = -1
  integer, dimension(odm), public :: mer_nflg
  real, dimension(odm), public :: mer_minlat=-90.
  real, dimension(odm), public :: mer_maxlat = 90.
  integer, public :: mer_nreg
  integer :: lmax
  real, dimension(ldm) :: mtlat

  real, allocatable, dimension(:,:,:) :: &
       mmflxl,mmftdl,mmfsml,mmflxd,mmftdd,mmfsmd
  real, allocatable, dimension(:,:) :: &
       mhflx,mhftd,mhfsm,mhfld,msflx,msftd,msfsm,msfld

  ! Section transports
  character(len = fnmlen), public :: sec_sifile
  integer :: sec_num
  integer, parameter :: max_sec = 400
  character(len = slenmax) :: sec_name(max_sec)
  real, dimension(max_sec) :: voltr

  ! Global sums and averages
  real, dimension(1) :: massgs,volgs,salnga,tempga,sssga,sstga

  ! Pressure thickness [kg m-1 s-2] of region for bottom salinity and
  ! temperature diagnostics
  real, parameter :: dpbot = onem

  ! Namelist
  integer, dimension(nphymax), public :: &
       H2D_ABSWND ,H2D_ALB    ,H2D_BTMSTR ,H2D_BRNFLX ,H2D_BRNPD  , &
       H2D_DFL    ,H2D_EVA    ,H2D_FICE   ,H2D_FMLTFZ ,H2D_HICE   , &
       H2D_HMLTFZ ,H2D_HSNW   ,H2D_IAGE   ,H2D_IDKEDT ,H2D_LAMULT , H2D_HMAT, &
       H2D_LASL   ,H2D_LIP    ,H2D_MAXMLD ,H2D_MLD    ,H2D_MLTS   , &
       H2D_MLTSMN ,H2D_MLTSMX ,H2D_MLTSSQ ,H2D_MTKEUS ,H2D_MTKENI , &
       H2D_MTKEBF ,H2D_MTKERS ,H2D_MTKEPE ,H2D_MTKEKE ,H2D_MTY    , &
       H2D_NSF    ,H2D_PBOT   ,H2D_PSRF   ,H2D_RFIFLX ,H2D_RNFFLX , &
       H2D_SALFLX ,H2D_SALRLX ,H2D_SBOT   ,H2D_SEALV  ,H2D_SLVSQ  , &
       H2D_SFL    ,H2D_SOP    ,H2D_SIGMX  ,H2D_SSS    ,H2D_SSSSQ  , &
       H2D_SST    ,H2D_SSTSQ  ,H2D_SURFLX ,H2D_SURRLX ,H2D_SWA    , &
       H2D_T20D   ,H2D_TAUX   ,H2D_TAUY   ,H2D_TBOT   ,H2D_TICE   , &
       H2D_TSRF   ,H2D_UB     ,H2D_UICE   ,H2D_USTAR  ,H2D_USTAR3 , &
       H2D_USTOKES,H2D_VB     ,H2D_VICE   ,H2D_VSTOKES,H2D_ZTX    , &
       LYR_BFSQ   ,LYR_DIFDIA ,LYR_DIFVMO ,LYR_DIFVHO ,LYR_DIFVSO , &
       LYR_DIFINT ,LYR_DIFISO ,LYR_DP     ,LYR_DPU    ,LYR_DPV    , &
       LYR_DZ     ,LYR_SALN   ,LYR_TEMP   ,LYR_TRC    ,LYR_UFLX   , &
       LYR_UTFLX  ,LYR_USFLX  ,LYR_UMFLTD ,LYR_UMFLSM ,LYR_UTFLTD , &
       LYR_UTFLSM ,LYR_UTFLLD ,LYR_USFLTD ,LYR_USFLSM ,LYR_USFLLD , &
       LYR_UVEL   ,LYR_VFLX   ,LYR_VTFLX  ,LYR_VSFLX  ,LYR_VMFLTD , &
       LYR_VMFLSM ,LYR_VTFLTD ,LYR_VTFLSM ,LYR_VTFLLD ,LYR_VSFLTD , &
       LYR_VSFLSM ,LYR_VSFLLD ,LYR_VVEL   ,LYR_WFLX   ,LYR_WFLX2  , &
       LYR_PV     ,LYR_TKE    ,LYR_GLS_PSI,LYR_IDLAGE , &
       LVL_BFSQ   ,LVL_DIFDIA ,LVL_DIFVMO ,LVL_DIFVHO ,LVL_DIFVSO , &
       LVL_DIFINT ,LVL_DIFISO ,LVL_DZ     ,LVL_SALN   ,LVL_TEMP   , &
       LVL_TRC    ,LVL_UFLX   ,LVL_UTFLX  ,LVL_USFLX  ,LVL_UMFLTD , &
       LVL_UMFLSM ,LVL_UTFLTD ,LVL_UTFLSM ,LVL_UTFLLD ,LVL_USFLTD , &
       LVL_USFLSM ,LVL_USFLLD ,LVL_UVEL   ,LVL_VFLX   ,LVL_VTFLX  , &
       LVL_VSFLX  ,LVL_VMFLTD ,LVL_VMFLSM ,LVL_VTFLTD ,LVL_VTFLSM , &
       LVL_VTFLLD ,LVL_VSFLTD ,LVL_VSFLSM ,LVL_VSFLLD ,LVL_VVEL   , &
       LVL_WFLX   ,LVL_WFLX2  ,LVL_PV     ,LVL_TKE    ,LVL_GLS_PSI, &
       LVL_IDLAGE , &
       MSC_MMFLXL ,MSC_MMFLXD ,MSC_MMFTDL ,MSC_MMFSML ,MSC_MMFTDD , &
       MSC_MMFSMD ,MSC_MHFLX  ,MSC_MHFTD  ,MSC_MHFSM  ,MSC_MHFLD  , &
       MSC_MSFLX  ,MSC_MSFTD  ,MSC_MSFSM  ,MSC_MSFLD  ,MSC_VOLTR  , &
       MSC_MASSGS ,MSC_VOLGS  ,MSC_SALNGA ,MSC_TEMPGA ,MSC_SSSGA  , &
       MSC_SSTGA

  integer, dimension(nphymax), public :: &
       ACC_ABSWND ,ACC_ALB    ,ACC_BRNFLX ,ACC_BRNPD  ,ACC_DFL    , &
       ACC_EVA    ,ACC_FICE   ,ACC_FMLTFZ ,ACC_HICE   ,ACC_HMAT   , ACC_HMLTFZ , &
       ACC_HSNW   ,ACC_IAGE   ,ACC_IDKEDT ,ACC_LAMULT ,ACC_LASL   , &
       ACC_LIP    ,ACC_MAXMLD ,ACC_MLD    ,ACC_MLTS   ,ACC_MLTSMN , &
       ACC_MLTSMX ,ACC_MLTSSQ ,ACC_MTKEUS ,ACC_MTKENI ,ACC_MTKEBF , &
       ACC_MTKERS ,ACC_MTKEPE ,ACC_MTKEKE ,ACC_MTY    ,ACC_NSF    , &
       ACC_PBOT   ,ACC_PSRF   ,ACC_RFIFLX ,ACC_RNFFLX ,ACC_SALFLX , &
       ACC_SALRLX ,ACC_SBOT   ,ACC_SEALV  ,ACC_SLVSQ  ,ACC_SFL    , &
       ACC_SOP    ,ACC_SIGMX  ,ACC_SSS    ,ACC_SSSSQ  ,ACC_SST    , &
       ACC_SSTSQ  ,ACC_SURFLX ,ACC_SURRLX ,ACC_SWA    ,ACC_T20D   , &
       ACC_TAUX   ,ACC_TAUY   ,ACC_TBOT   ,ACC_TICE   ,ACC_TSRF   , &
       ACC_UB     ,ACC_UBFLXS ,ACC_UICE   ,ACC_USTAR  ,ACC_USTAR3 , &
       ACC_USTOKES,ACC_VB     ,ACC_VBFLXS ,ACC_VICE   ,ACC_VSTOKES, &
       ACC_ZTX    ,ACC_IVOLU  ,ACC_IVOLV  ,ACC_UTILH2D, &
       ACC_BFSQ   ,ACC_DIFDIA ,ACC_DIFVMO ,ACC_DIFVHO ,ACC_DIFVSO , &
       ACC_DIFINT ,ACC_DIFISO ,ACC_DP     ,ACC_DPU    ,ACC_DPV    , &
       ACC_DZ     ,ACC_SALN   ,ACC_TEMP   ,ACC_UFLX   ,ACC_UTFLX  , &
       ACC_USFLX  ,ACC_UMFLTD ,ACC_UMFLSM ,ACC_UTFLTD ,ACC_UTFLSM , &
       ACC_UTFLLD ,ACC_USFLTD ,ACC_USFLSM ,ACC_USFLLD ,ACC_UVEL   , &
       ACC_VFLX   ,ACC_VTFLX  ,ACC_VSFLX  ,ACC_VMFLTD ,ACC_VMFLSM , &
       ACC_VTFLTD ,ACC_VTFLSM ,ACC_VTFLLD ,ACC_VSFLTD ,ACC_VSFLSM , &
       ACC_VSFLLD ,ACC_VVEL   ,ACC_WFLX   ,ACC_WFLX2  ,ACC_AVDSG  , &
       ACC_DPVOR  ,ACC_TKE    ,ACC_GLS_PSI,ACC_UTILLYR, &
       ACC_BFSQLVL   ,ACC_DIFDIALVL ,ACC_DIFVMOLVL ,ACC_DIFVHOLVL , &
       ACC_DIFVSOLVL ,ACC_DIFINTLVL ,ACC_DIFISOLVL ,ACC_DZLVL     , &
       ACC_SALNLVL   ,ACC_TEMPLVL   ,ACC_UFLXLVL   ,ACC_UTFLXLVL  , &
       ACC_USFLXLVL  ,ACC_UMFLTDLVL ,ACC_UMFLSMLVL ,ACC_UTFLTDLVL , &
       ACC_UTFLSMLVL ,ACC_UTFLLDLVL ,ACC_USFLTDLVL ,ACC_USFLSMLVL , &
       ACC_USFLLDLVL ,ACC_UVELLVL   ,ACC_VFLXLVL   ,ACC_VTFLXLVL  , &
       ACC_VSFLXLVL  ,ACC_VMFLTDLVL ,ACC_VMFLSMLVL ,ACC_VTFLTDLVL , &
       ACC_VTFLSMLVL ,ACC_VTFLLDLVL ,ACC_VSFLTDLVL ,ACC_VSFLSMLVL , &
       ACC_VSFLLDLVL ,ACC_VVELLVL   ,ACC_WFLXLVL   ,ACC_WFLX2LVL  , &
       ACC_PVLVL     ,ACC_TKELVL    ,ACC_GLS_PSILVL,ACC_IDLAGELVL , &
       ACC_UFLXOLD   ,ACC_VFLXOLD   ,ACC_UTILLVL   , &
       ACC_MMFLXL,ACC_MMFLXD,ACC_MMFTDL,ACC_MMFSML,ACC_MMFTDD, &
       ACC_MMFSMD,ACC_MHFLX ,ACC_MHFTD ,ACC_MHFSM ,ACC_MHFLD , &
       ACC_MSFLX ,ACC_MSFTD ,ACC_MSFSM ,ACC_MSFLD ,ACC_VOLTR

  character(len=10) , dimension(nphymax), public :: GLB_FNAMETAG
  integer           , dimension(nphymax), public :: GLB_AVEPERIO
  integer           , dimension(nphymax), public :: GLB_FILEFREQ
  integer           , dimension(nphymax), public :: GLB_COMPFLAG
  integer           , dimension(nphymax), public :: GLB_NCFORMAT

  namelist /MERDIA/ &
       MER_ORFILE,MER_MIFILE,MER_REGNAM,MER_REGFLG,MER_MINLAT,MER_MAXLAT

  namelist /SECDIA/ &
       SEC_SIFILE

  namelist /DIAPHY/ &
       H2D_ABSWND ,H2D_ALB    ,H2D_BTMSTR ,H2D_BRNFLX ,H2D_BRNPD  , &
       H2D_DFL    ,H2D_EVA    ,H2D_FICE   ,H2D_FMLTFZ ,H2D_HICE   , H2D_HMAT, &
       H2D_HMLTFZ ,H2D_HSNW   ,H2D_IAGE   ,H2D_IDKEDT ,H2D_LAMULT , &
       H2D_LASL   ,H2D_LIP    ,H2D_MAXMLD ,H2D_MLD    ,H2D_MLTS   , &
       H2D_MLTSMN ,H2D_MLTSMX ,H2D_MLTSSQ ,H2D_MTKEUS ,H2D_MTKENI , &
       H2D_MTKEBF ,H2D_MTKERS ,H2D_MTKEPE ,H2D_MTKEKE ,H2D_MTY    , &
       H2D_NSF    ,H2D_PBOT   ,H2D_PSRF   ,H2D_RFIFLX ,H2D_RNFFLX , &
       H2D_SALFLX ,H2D_SALRLX ,H2D_SBOT   ,H2D_SEALV  ,H2D_SLVSQ  , &
       H2D_SFL    ,H2D_SOP    ,H2D_SIGMX  ,H2D_SSS    ,H2D_SSSSQ  , &
       H2D_SST    ,H2D_SSTSQ  ,H2D_SURFLX ,H2D_SURRLX ,H2D_SWA    , &
       H2D_T20D   ,H2D_TAUX   ,H2D_TAUY   ,H2D_TBOT   ,H2D_TICE   , &
       H2D_TSRF   ,H2D_UB     ,H2D_UICE   ,H2D_USTAR  ,H2D_USTAR3 , &
       H2D_USTOKES,H2D_VB     ,H2D_VICE   ,H2D_VSTOKES,H2D_ZTX    , &
       LYR_BFSQ   ,LYR_DIFDIA ,LYR_DIFVMO ,LYR_DIFVHO ,LYR_DIFVSO , &
       LYR_DIFINT ,LYR_DIFISO ,LYR_DP     ,LYR_DPU    ,LYR_DPV    , &
       LYR_DZ     ,LYR_SALN   ,LYR_TEMP   ,LYR_TRC    ,LYR_UFLX   , &
       LYR_UTFLX  ,LYR_USFLX  ,LYR_UMFLTD ,LYR_UMFLSM ,LYR_UTFLTD , &
       LYR_UTFLSM ,LYR_UTFLLD ,LYR_USFLTD ,LYR_USFLSM ,LYR_USFLLD , &
       LYR_UVEL   ,LYR_VFLX   ,LYR_VTFLX  ,LYR_VSFLX  ,LYR_VMFLTD , &
       LYR_VMFLSM ,LYR_VTFLTD ,LYR_VTFLSM ,LYR_VTFLLD ,LYR_VSFLTD , &
       LYR_VSFLSM ,LYR_VSFLLD ,LYR_VVEL   ,LYR_WFLX   ,LYR_WFLX2  , &
       LYR_PV     ,LYR_TKE    ,LYR_GLS_PSI,LYR_IDLAGE , &
       LVL_BFSQ   ,LVL_DIFDIA ,LVL_DIFVMO ,LVL_DIFVHO ,LVL_DIFVSO , &
       LVL_DIFINT ,LVL_DIFISO ,LVL_DZ     ,LVL_SALN   ,LVL_TEMP   , &
       LVL_TRC    ,LVL_UFLX   ,LVL_UTFLX  ,LVL_USFLX  ,LVL_UMFLTD , &
       LVL_UMFLSM ,LVL_UTFLTD ,LVL_UTFLSM ,LVL_UTFLLD ,LVL_USFLTD , &
       LVL_USFLSM ,LVL_USFLLD ,LVL_UVEL   ,LVL_VFLX   ,LVL_VTFLX  , &
       LVL_VSFLX  ,LVL_VMFLTD ,LVL_VMFLSM ,LVL_VTFLTD ,LVL_VTFLSM , &
       LVL_VTFLLD ,LVL_VSFLTD ,LVL_VSFLSM ,LVL_VSFLLD ,LVL_VVEL   , &
       LVL_WFLX   ,LVL_WFLX2  ,LVL_PV     ,LVL_TKE    ,LVL_GLS_PSI, &
       LVL_IDLAGE , &
       MSC_MMFLXL ,MSC_MMFLXD ,MSC_MMFTDL ,MSC_MMFSML ,MSC_MMFTDD , &
       MSC_MMFSMD ,MSC_MHFLX  ,MSC_MHFTD  ,MSC_MHFSM  ,MSC_MHFLD  , &
       MSC_MSFLX  ,MSC_MSFTD  ,MSC_MSFSM  ,MSC_MSFLD  ,MSC_VOLTR  , &
       MSC_MASSGS ,MSC_VOLGS  ,MSC_SALNGA ,MSC_TEMPGA ,MSC_SSSGA  , &
       MSC_SSTGA  , &
       GLB_AVEPERIO,GLB_FILEFREQ,GLB_COMPFLAG,GLB_NCFORMAT, &
       GLB_FNAMETAG

  ! public namelists
  public :: DIAPHY, MERDIA, SECDIA

  ! Public routines
  public :: diafnm
  public :: diaini
  public :: diaacc
  public :: diaout_alarms
  public :: diaout
  public :: diasec
  public :: diamer
  public :: diavfl
  public :: diazlv

  ! Private routines
  private :: inih2d
  private :: inilyr
  private :: inilvl
  private :: acch2d
  private :: maxh2d
  private :: minh2d
  private :: sqh2d
  private :: acclyr
  private :: accily
  private :: acclvl
  private :: accilv
  private :: inifld
  private :: finh2d
  private :: finlyr
  private :: wrth2d
  private :: wrtlyr
  private :: wrtlvl
  private :: logh2d
  private :: loglyr
  private :: loglvl
  private :: msklvl
  private :: definevar

contains

  subroutine diafnm(ctag,diagfq,diagmon,diagann,fname)
  !---------------------------------------------------------------
  ! Description: creates file name for the diagnostic output
  !---------------------------------------

    ! Arguments
    character(len = *), intent(in) :: ctag    ! string used in middle of file name
    real,               intent(in) :: diagfq  ! diagnostic frequency
    logical,            intent(in) :: diagmon ! switch to show whether diagfq=month
    logical,            intent(in) :: diagann ! switch to show whether diagfq=year
    character(len = *), intent(out):: fname   ! file name

    ! Local variables
    character(len = fnmlen) :: prefix
    character(len = 1) :: sep1,sep2
    type(date_type) :: date_tmp
    integer :: errstat,ns

    if (expcnf == 'cesm') then
      prefix = trim(runid)//'.blom'//trim(inst_suffix)
      sep1 = '.'
      sep2 = '-'
    else
      prefix = trim(runid)
      sep1 = '_'
      sep2 = '.'
    end if

    date_tmp = date

    if (diagfq+epsilp > 1.) then
      errstat = date_offset(calendar,date_tmp,-1)
      if (errstat /= calendar_noerr) then
        if (mnproc == 1) then
          write (lp, '(2a)') ' diafnm: date_offset error: ', &
               trim(calendar_errstr(errstat))
        end if
        call xcstop('(diafnm)')
        stop '(diafnm)'
      end if
      if     (diagann) then
        write(fname,'(4a,i4.4,a)') &
             trim(prefix),sep1,trim(ctag),sep1, &
             date_tmp%year,'.nc'
      else if (diagmon) then
        write(fname,'(4a,i4.4,a,i2.2,a)') &
             trim(prefix),sep1,trim(ctag),sep1, &
             date_tmp%year,sep2,date_tmp%month,'.nc'
      else
        write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)') &
             trim(prefix),sep1,trim(ctag),sep1, &
             date_tmp%year,sep2,date_tmp%month,sep2,date_tmp%day,'.nc'
      end if
    else
      if (mod(nstep,nstep_in_day) == 0) then
        errstat = date_offset(calendar,date_tmp,-1)
        if (errstat /= calendar_noerr) then
          if (mnproc == 1) then
            write (lp, '(2a)') ' diafnm: date_offset error: ', &
                 trim(calendar_errstr(errstat))
          end if
          call xcstop('(diafnm)')
          stop '(diafnm)'
        end if
      end if
      ns = 86400*(mod(nstep-1,nstep_in_day)+1)/nstep_in_day
      write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') &
           trim(prefix),sep1,trim(ctag),sep1, &
           date_tmp%year,sep2,date_tmp%month,sep2,date_tmp%day,sep2,ns, &
           '.nc'
    end if

  end subroutine diafnm



  subroutine diaini
  !---------------------------------------------------------------
  ! initialize diagnostic variables
  !---------------------------------------------------------------

    integer :: i,j,l,n,istat,istatsum
    integer, parameter :: imn=1-nbdy,imx=idm+nbdy,jmn=imn,jmx = jdm+nbdy
    logical :: fexist

    ! Check existence of data files for meridional and section transport
    ! diagnostics
    if (mnproc == 1) then
      if (sum(MSC_MMFLXL(1:nphy)+MSC_MMFLXD(1:nphy)+MSC_MMFTDL(1:nphy) &
           +MSC_MMFSML(1:nphy)+MSC_MMFTDD(1:nphy)+MSC_MMFSMD(1:nphy) &
           +MSC_MHFLX (1:nphy)+MSC_MHFTD (1:nphy)+MSC_MHFSM (1:nphy) &
           +MSC_MHFLD (1:nphy)+MSC_MSFLX (1:nphy)+MSC_MSFTD (1:nphy) &
           +msc_msfsm (1:nphy)+msc_msfld (1:nphy)) /= 0) then
        inquire(file=mer_orfile,exist = fexist)
        if (.not.fexist) then
          write (lp,'(3a)') ' Could not find file ', &
               trim(mer_orfile),'!'
          call xchalt('(diaini)')
          stop '(diaini)'
        end if
        inquire(file=mer_mifile,exist = fexist)
        if (.not.fexist) then
          write (lp,'(3a)') ' Could not find file ', &
               trim(mer_mifile),'!'
          call xchalt('(diaini)')
          stop '(diaini)'
        end if
      end if
      if (sum(msc_voltr(1:nphy)) /= 0) then
        inquire(file=sec_sifile,exist = fexist)
        if (.not.fexist) then
          write (lp,'(3a)') ' Could not find file ', &
               trim(sec_sifile),'!'
          call xchalt('(diaini)')
          stop '(diaini)'
        end if
      end if
    end if

    nphyh2d = 0
    nphylyr = 0
    nphylvl = 0

    ! Loop over io groups
    do n = 1,nphy
      nacc_phy(n) = 0

      ! Activate alarm for io group output to avoid accumulation of diagnostic
      ! variables in call to ALE regrid/remapping during initialization
      alarm_phy(n) = 1

      ! - Solve dependencies for diagnostic variables (0=skipped)
      ACC_ABSWND(n)   = H2D_ABSWND(n)
      ACC_ALB(n)      = H2D_ALB(n)
      ACC_BRNFLX(n)   = H2D_BRNFLX(n)
      ACC_BRNPD(n)    = H2D_BRNPD(n)
      ACC_DFL(n)      = H2D_DFL(n)
      ACC_EVA(n)      = H2D_EVA(n)
      ACC_FMLTFZ(n)   = H2D_FMLTFZ(n)
      ACC_FICE(n)     = H2D_FICE(n)   + H2D_HICE(n)   + H2D_UICE(n) &
           + H2D_VICE(n)   + H2D_HSNW(n)
      ACC_HICE(n)     = H2D_HICE(n)   + H2D_UICE(n)   + H2D_VICE(n)
      ACC_HMAT(n)     = H2D_HMAT(n)
      ACC_HMLTFZ(n)   = H2D_HMLTFZ(n)
      ACC_HSNW(n)     = H2D_HSNW(n)
      ACC_IAGE(n)     = H2D_IAGE(n)
      ACC_IDKEDT(n)   = H2D_IDKEDT(n)
      ACC_IVOLU(n)    = H2D_UICE(n)
      ACC_IVOLV(n)    = H2D_VICE(n)
      ACC_LAMULT(n)   = H2D_LAMULT(n)
      ACC_LASL(n)     = H2D_LASL(n)
      ACC_LIP(n)      = H2D_LIP(n)
      ACC_MAXMLD(n)   = H2D_MAXMLD(n)
      ACC_MLD(n)      = H2D_MLD(n)
      ACC_MLTS(n)     = H2D_MLTS(n)
      ACC_MLTSMN(n)   = H2D_MLTSMN(n)
      ACC_MLTSMX(n)   = H2D_MLTSMX(n)
      ACC_MLTSSQ(n)   = H2D_MLTSSQ(n)
      ACC_MTKEUS(n)   = H2D_MTKEUS(n)
      ACC_MTKENI(n)   = H2D_MTKENI(n)
      ACC_MTKEBF(n)   = H2D_MTKEBF(n)
      ACC_MTKERS(n)   = H2D_MTKERS(n)
      ACC_MTKEPE(n)   = H2D_MTKEPE(n)
      ACC_MTKEKE(n)   = H2D_MTKEKE(n)
      ACC_MTY(n)      = H2D_MTY(n)
      ACC_NSF(n)      = H2D_NSF(n)
      ACC_PBOT(n)     = H2D_PBOT(n)
      ACC_PSRF(n)     = H2D_PSRF(n)
      ACC_RFIFLX(n)   = H2D_RFIFLX(n)
      ACC_RNFFLX(n)   = H2D_RNFFLX(n)
      ACC_SURFLX(n)   = H2D_SURFLX(n)
      ACC_SURRLX(n)   = H2D_SURRLX(n)
      ACC_SALFLX(n)   = H2D_SALFLX(n)
      ACC_SALRLX(n)   = H2D_SALRLX(n)
      ACC_SBOT(n)     = H2D_SBOT(n)
      ACC_SEALV(n)    = H2D_SEALV(n)
      ACC_SLVSQ(n)    = H2D_SLVSQ(n)
      ACC_SFL(n)      = H2D_SFL(n)
      ACC_SIGMX(n)    = H2D_SIGMX(n)  + MSC_MMFLXL(n) + MSC_MMFTDL(n) &
           + MSC_MMFSML(n)
      ACC_SOP(n)      = H2D_SOP(n)
      ACC_SSS(n)      = H2D_SSS(n)    + MSC_SSSGA(n)
      ACC_SSSSQ(n)    = H2D_SSSSQ(n)
      ACC_SST(n)      = H2D_SST(n)    + MSC_SSTGA(n)
      ACC_SSTSQ(n)    = H2D_SSTSQ(n)
      ACC_SWA(n)      = H2D_SWA(n)
      ACC_T20D(n)     = H2D_T20D(n)
      ACC_TAUX(n)     = H2D_TAUX(n)
      ACC_TAUY(n)     = H2D_TAUY(n)
      ACC_TBOT(n)     = H2D_TBOT(n)
      ACC_TICE(n)     = H2D_TICE(n)
      ACC_TSRF(n)     = H2D_TSRF(n)
      ACC_UB(n)       = H2D_UB(n)
      ACC_UBFLXS(n)   = H2D_BTMSTR(n)
      ACC_UICE(n)     = H2D_UICE(n)
      ACC_USTAR(n)    = H2D_USTAR(n)
      ACC_USTAR3(n)   = H2D_USTAR3(n)
      ACC_USTOKES(n)  = H2D_USTOKES(n)
      ACC_VB(n)       = H2D_VB(n)
      ACC_VBFLXS(n)   = H2D_BTMSTR(n)
      ACC_VICE(n)     = H2D_VICE(n)
      ACC_VSTOKES(n)  = H2D_VSTOKES(n)
      ACC_ZTX(n)      = H2D_ZTX(n)
      ACC_BFSQ(n)     = LYR_BFSQ(n)
      ACC_BFSQLVL(n)  = LVL_BFSQ(n)
      ACC_DIFDIA(n)   = LYR_DIFDIA(n)
      ACC_DIFDIALVL(n)= LVL_DIFDIA(n)
      ACC_DIFVMO(n)   = LYR_DIFVMO(n)
      ACC_DIFVMOLVL(n)= LVL_DIFVMO(n)
      ACC_DIFVHO(n)   = LYR_DIFVHO(n)
      ACC_DIFVHOLVL(n)= LVL_DIFVHO(n)
      ACC_DIFVSO(n)   = LYR_DIFVSO(n)
      ACC_DIFVSOLVL(n)= LVL_DIFVSO(n)
      ACC_DIFINT(n)   = LYR_DIFINT(n)
      ACC_DIFINTLVL(n)= LVL_DIFINT(n)
      ACC_DIFISO(n)   = LYR_DIFISO(n)
      ACC_DIFISOLVL(n)= LVL_DIFISO(n)
      ACC_DP(n)       = LYR_DP(n)     + LYR_BFSQ(n)   + LYR_SALN(n) &
           + LYR_TEMP(n)   + LYR_DIFDIA(n) + LYR_DIFVMO(n) &
           + LYR_DIFVHO(n) + LYR_DIFVSO(n) + LYR_DIFINT(n) &
           + LYR_DIFISO(n) + LYR_TKE(n)    + LYR_GLS_PSI(n) &
           + LVL_BFSQ(n)   + LVL_SALN(n)   + LVL_TEMP(n) &
           + LVL_DIFDIA(n) + LVL_DIFVMO(n) + LVL_DIFVHO(n) &
           + LVL_DIFVSO(n) + LVL_DIFINT(n) + LVL_DIFISO(n) &
           + LVL_TKE(n)    + LVL_GLS_PSI(n) &
           + MSC_MASSGS(n) + MSC_SALNGA(n) + MSC_TEMPGA(n)
      ACC_DPU(n)      = LYR_DPU(n)    + LYR_UVEL(n)
      ACC_DPV(n)      = LYR_DPV(n)    + LYR_VVEL(n)
      ACC_DZ(n)       = LYR_DZ(n)     + MSC_VOLGS(n)
      ACC_DZLVL(n)    = LVL_DZ(n)
      ACC_SALN(n)     = LYR_SALN(n)   + MSC_SALNGA(n)
      ACC_SALNLVL(n)  = LVL_SALN(n)
      ACC_TEMP(n)     = LYR_TEMP(n)   + MSC_TEMPGA(n)
      ACC_TEMPLVL(n)  = LVL_TEMP(n)
      ACC_UFLX(n)     = LYR_UFLX(n)   + MSC_MMFLXL(n) + LYR_WFLX(n) &
           + LYR_WFLX2(n)
      ACC_UFLXLVL(n)  = LVL_UFLX(n)   + MSC_MMFLXD(n) + MSC_VOLTR(n) &
           + LVL_WFLX(n)   + LVL_WFLX2(n)
      ACC_UFLXOLD(n)  = LVL_WFLX(n)   + LVL_WFLX2(n)
      ACC_UTFLX(n)    = LYR_UTFLX(n)  + MSC_MHFLX(n)
      ACC_UTFLXLVL(n) = LVL_UTFLX(n)
      ACC_USFLX(n)    = LYR_USFLX(n)  + MSC_MSFLX(n)
      ACC_USFLXLVL(n) = LVL_USFLX(n)
      ACC_UMFLTD(n)   = LYR_UMFLTD(n) + MSC_MMFTDL(n)
      ACC_UMFLSM(n)   = LYR_UMFLSM(n) + MSC_MMFSML(n)
      ACC_UMFLTDLVL(n)= LVL_UMFLTD(n) + MSC_MMFTDD(n)
      ACC_UMFLSMLVL(n)= LVL_UMFLSM(n) + MSC_MMFSMD(n)
      ACC_UTFLTD(n)   = LYR_UTFLTD(n) + MSC_MHFTD(n)
      ACC_UTFLSM(n)   = LYR_UTFLSM(n) + MSC_MHFSM(n)
      ACC_UTFLTDLVL(n)= LVL_UTFLTD(n)
      ACC_UTFLSMLVL(n)= LVL_UTFLSM(n)
      ACC_UTFLLD(n)   = LYR_UTFLLD(n) + MSC_MHFLD(n)
      ACC_UTFLLDLVL(n)= LVL_UTFLLD(n)
      ACC_USFLTD(n)   = LYR_USFLTD(n) + MSC_MSFTD(n)
      ACC_USFLSM(n)   = LYR_USFLSM(n) + MSC_MSFSM(n)
      ACC_USFLTDLVL(n)= LVL_USFLTD(n)
      ACC_USFLSMLVL(n)= LVL_USFLSM(n)
      ACC_USFLLD(n)   = LYR_USFLLD(n) + MSC_MSFLD(n)
      ACC_USFLLDLVL(n)= LVL_USFLLD(n)
      ACC_UVEL(n)     = LYR_UVEL(n)
      ACC_UVELLVL(n)  = LVL_UVEL(n)
      ACC_VFLX(n)     = LYR_VFLX(n)   + MSC_MMFLXL(n) + LYR_WFLX(n) &
           + LYR_WFLX2(n)
      ACC_VFLXLVL(n)  = LVL_VFLX(n)   + MSC_MMFLXD(n) + MSC_VOLTR(n) &
           + LVL_WFLX(n)   + LVL_WFLX2(n)
      ACC_VFLXOLD(n)  = LVL_WFLX(n)   + LVL_WFLX2(n)
      ACC_VTFLX(n)    = LYR_VTFLX(n)  + MSC_MHFLX(n)
      ACC_VTFLXLVL(n) = LVL_VTFLX(n)
      ACC_VSFLX(n)    = LYR_VSFLX(n)  + MSC_MSFLX(n)
      ACC_VSFLXLVL(n) = LVL_VSFLX(n)
      ACC_VMFLTD(n)   = LYR_VMFLTD(n) + MSC_MMFTDL(n)
      ACC_VMFLSM(n)   = LYR_VMFLSM(n) + MSC_MMFSML(n)
      ACC_VMFLTDLVL(n)= LVL_VMFLTD(n) + MSC_MMFTDD(n)
      ACC_VMFLSMLVL(n)= LVL_VMFLSM(n) + MSC_MMFSMD(n)
      ACC_VTFLTD(n)   = LYR_VTFLTD(n) + MSC_MHFTD(n)
      ACC_VTFLSM(n)   = LYR_VTFLSM(n) + MSC_MHFSM(n)
      ACC_VTFLTDLVL(n)= LVL_VTFLTD(n)
      ACC_VTFLSMLVL(n)= LVL_VTFLSM(n)
      ACC_VTFLLD(n)   = LYR_VTFLLD(n) + MSC_MHFLD(n)
      ACC_VTFLLDLVL(n)= LVL_VTFLLD(n)
      ACC_VSFLTD(n)   = LYR_VSFLTD(n) + MSC_MSFTD(n)
      ACC_VSFLSM(n)   = LYR_VSFLSM(n) + MSC_MSFSM(n)
      ACC_VSFLTDLVL(n)= LVL_VSFLTD(n)
      ACC_VSFLSMLVL(n)= LVL_VSFLSM(n)
      ACC_VSFLLD(n)   = LYR_VSFLLD(n) + MSC_MSFLD(n)
      ACC_VSFLLDLVL(n)= LVL_VSFLLD(n)
      ACC_VVEL(n)     = LYR_VVEL(n)
      ACC_VVELLVL(n)  = LVL_VVEL(n)
      ACC_WFLX(n)     = LYR_WFLX(n)   + LYR_WFLX2(n)  + LVL_WFLX(n) &
           + LVL_WFLX2(n)
      ACC_WFLXLVL(n)  = LVL_WFLX(n)   + LVL_WFLX2(n)  + LYR_WFLX(n) &
           + LYR_WFLX2(n)
      ACC_WFLX2(n)    = LYR_WFLX2(n)  + LYR_WFLX(n)   + LVL_WFLX(n) &
           + LVL_WFLX2(n)
      ACC_WFLX2LVL(n) = LVL_WFLX2(n)  + LVL_WFLX(n)   + LYR_WFLX(n) &
           + LYR_WFLX2(n)
      ACC_AVDSG(n)    = LYR_PV(n)
      ACC_DPVOR(n)    = LYR_PV(n)
      ACC_PVLVL(n)    = LVL_PV(n)
      ACC_TKE(n)      = LYR_TKE(n)
      ACC_TKELVL(n)   = LVL_TKE(n)
      ACC_GLS_PSI(n)  = LYR_GLS_PSI(n)
      ACC_GLS_PSILVL(n) = LVL_GLS_PSI(n)
      ACC_IDLAGELVL(n)= LVL_IDLAGE(n)
      ACC_MMFLXL(n)   = MSC_MMFLXL(n)
      ACC_MMFLXD(n)   = MSC_MMFLXD(n)
      ACC_MMFTDL(n)   = MSC_MMFTDL(n)
      ACC_MMFSML(n)   = MSC_MMFSML(n)
      ACC_MMFTDD(n)   = MSC_MMFTDD(n)
      ACC_MMFSMD(n)   = MSC_MMFSMD(n)
      ACC_MHFLX(n)    = MSC_MHFLX(n)
      ACC_MHFTD(n)    = MSC_MHFTD(n)
      ACC_MHFSM(n)    = MSC_MHFSM(n)
      ACC_MHFLD(n)    = MSC_MHFLD(n)
      ACC_MSFLX(n)    = MSC_MSFLX(n)
      ACC_MSFTD(n)    = MSC_MSFTD(n)
      ACC_MSFSM(n)    = MSC_MSFSM(n)
      ACC_MSFLD(n)    = MSC_MSFLD(n)
      ACC_VOLTR(n)    = MSC_VOLTR(n)

      ! - Determine position in buffer
      if (acc_abswnd(n) /= 0) nphyh2d = nphyh2d+1
      ACC_ABSWND(n) = nphyh2d*min(1,ACC_ABSWND(n))
      if (acc_alb(n) /= 0) nphyh2d = nphyh2d+1
      ACC_ALB(n) = nphyh2d*min(1,ACC_ALB(n))
      if (acc_brnflx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_BRNFLX(n) = nphyh2d*min(1,ACC_BRNFLX(n))
      if (acc_brnpd(n) /= 0) nphyh2d = nphyh2d+1
      ACC_BRNPD(n) = nphyh2d*min(1,ACC_BRNPD(n))
      if (acc_dfl(n) /= 0) nphyh2d = nphyh2d+1
      ACC_DFL(n) = nphyh2d*min(1,ACC_DFL(n))
      if (acc_eva(n) /= 0) nphyh2d = nphyh2d+1
      ACC_EVA(n) = nphyh2d*min(1,ACC_EVA(n))
      if (acc_fmltfz(n) /= 0) nphyh2d = nphyh2d+1
      ACC_FMLTFZ(n) = nphyh2d*min(1,ACC_FMLTFZ(n))
      if (acc_fice(n) /= 0) nphyh2d = nphyh2d+1
      ACC_FICE(n) = nphyh2d*min(1,ACC_FICE(n))
      if (acc_hice(n) /= 0) nphyh2d = nphyh2d+1
      ACC_HICE(n) = nphyh2d*min(1,ACC_HICE(n))
      if (acc_hmat(n) /= 0) nphyh2d = nphyh2d+1
      ACC_HMAT(n) = nphyh2d*min(1,ACC_HMAT(n))
      if (acc_hmltfz(n) /= 0) nphyh2d = nphyh2d+1
      ACC_HMLTFZ(n) = nphyh2d*min(1,ACC_HMLTFZ(n))
      if (acc_hsnw(n) /= 0) nphyh2d = nphyh2d+1
      ACC_HSNW(n) = nphyh2d*min(1,ACC_HSNW(n))
      if (acc_iage(n) /= 0) nphyh2d = nphyh2d+1
      ACC_IAGE(n) = nphyh2d*min(1,ACC_IAGE(n))
      if (acc_idkedt(n) /= 0) nphyh2d = nphyh2d+1
      ACC_IDKEDT(n) = nphyh2d*min(1,ACC_IDKEDT(n))
      if (acc_ivolu(n) /= 0) nphyh2d = nphyh2d+1
      ACC_IVOLU(n) = nphyh2d*min(1,ACC_IVOLU(n))
      if (acc_ivolv(n) /= 0) nphyh2d = nphyh2d+1
      ACC_IVOLV(n) = nphyh2d*min(1,ACC_IVOLV(n))
      if (acc_lamult(n) /= 0) nphyh2d = nphyh2d+1
      ACC_LAMULT(n) = nphyh2d*min(1,ACC_LAMULT(n))
      if (acc_lasl(n) /= 0) nphyh2d = nphyh2d+1
      ACC_LASL(n) = nphyh2d*min(1,ACC_LASL(n))
      if (acc_lip(n) /= 0) nphyh2d = nphyh2d+1
      ACC_LIP(n) = nphyh2d*min(1,ACC_LIP(n))
      if (acc_maxmld(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MAXMLD(n) = nphyh2d*min(1,ACC_MAXMLD(n))
      if (acc_mld(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MLD(n) = nphyh2d*min(1,ACC_MLD(n))
      if (acc_mlts(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MLTS(n) = nphyh2d*min(1,ACC_MLTS(n))
      if (acc_mltsmn(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MLTSMN(n) = nphyh2d*min(1,ACC_MLTSMN(n))
      if (acc_mltsmx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MLTSMX(n) = nphyh2d*min(1,ACC_MLTSMX(n))
      if (acc_mltssq(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MLTSSQ(n) = nphyh2d*min(1,ACC_MLTSSQ(n))
      if (acc_mtkeus(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTKEUS(n) = nphyh2d*min(1,ACC_MTKEUS(n))
      if (acc_mtkeni(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTKENI(n) = nphyh2d*min(1,ACC_MTKENI(n))
      if (acc_mtkebf(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTKEBF(n) = nphyh2d*min(1,ACC_MTKEBF(n))
      if (acc_mtkers(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTKERS(n) = nphyh2d*min(1,ACC_MTKERS(n))
      if (acc_mtkepe(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTKEPE(n) = nphyh2d*min(1,ACC_MTKEPE(n))
      if (acc_mtkeke(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTKEKE(n) = nphyh2d*min(1,ACC_MTKEKE(n))
      if (acc_mty(n) /= 0) nphyh2d = nphyh2d+1
      ACC_MTY(n) = nphyh2d*min(1,ACC_MTY(n))
      if (acc_nsf(n) /= 0) nphyh2d = nphyh2d+1
      ACC_NSF(n) = nphyh2d*min(1,ACC_NSF(n))
      if (acc_pbot(n) /= 0) nphyh2d = nphyh2d+1
      ACC_PBOT(n) = nphyh2d*min(1,ACC_PBOT(n))
      if (acc_psrf(n) /= 0) nphyh2d = nphyh2d+1
      ACC_PSRF(n) = nphyh2d*min(1,ACC_PSRF(n))
      if (acc_rfiflx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_RFIFLX(n) = nphyh2d*min(1,ACC_RFIFLX(n))
      if (acc_rnfflx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_RNFFLX(n) = nphyh2d*min(1,ACC_RNFFLX(n))
      if (acc_surflx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SURFLX(n) = nphyh2d*min(1,ACC_SURFLX(n))
      if (acc_surrlx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SURRLX(n) = nphyh2d*min(1,ACC_SURRLX(n))
      if (acc_salflx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SALFLX(n) = nphyh2d*min(1,ACC_SALFLX(n))
      if (acc_salrlx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SALRLX(n) = nphyh2d*min(1,ACC_SALRLX(n))
      if (acc_sbot(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SBOT(n) = nphyh2d*min(1,ACC_SBOT(n))
      if (acc_sealv(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SEALV(n) = nphyh2d*min(1,ACC_SEALV(n))
      if (acc_slvsq(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SLVSQ(n) = nphyh2d*min(1,ACC_SLVSQ(n))
      if (acc_sfl(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SFL(n) = nphyh2d*min(1,ACC_SFL(n))
      if (acc_sigmx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SIGMX(n) = nphyh2d*min(1,ACC_SIGMX(n))
      if (acc_sop(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SOP(n) = nphyh2d*min(1,ACC_SOP(n))
      if (acc_sss(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SSS(n) = nphyh2d*min(1,ACC_SSS(n))
      if (acc_ssssq(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SSSSQ(n) = nphyh2d*min(1,ACC_SSSSQ(n))
      if (acc_sst(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SST(n) = nphyh2d*min(1,ACC_SST(n))
      if (acc_sstsq(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SSTSQ(n) = nphyh2d*min(1,ACC_SSTSQ(n))
      if (acc_swa(n) /= 0) nphyh2d = nphyh2d+1
      ACC_SWA(n) = nphyh2d*min(1,ACC_SWA(n))
      if (acc_t20d(n) /= 0) nphyh2d = nphyh2d+1
      ACC_T20D(n) = nphyh2d*min(1,ACC_T20D(n))
      if (acc_taux(n) /= 0) nphyh2d = nphyh2d+1
      ACC_TAUX(n) = nphyh2d*min(1,ACC_TAUX(n))
      if (acc_tauy(n) /= 0) nphyh2d = nphyh2d+1
      ACC_TAUY(n) = nphyh2d*min(1,ACC_TAUY(n))
      if (acc_tbot(n) /= 0) nphyh2d = nphyh2d+1
      ACC_TBOT(n) = nphyh2d*min(1,ACC_TBOT(n))
      if (acc_tice(n) /= 0) nphyh2d = nphyh2d+1
      ACC_TICE(n) = nphyh2d*min(1,ACC_TICE(n))
      if (acc_tsrf(n) /= 0) nphyh2d = nphyh2d+1
      ACC_TSRF(n) = nphyh2d*min(1,ACC_TSRF(n))
      if (acc_ub(n) /= 0) nphyh2d = nphyh2d+1
      ACC_UB(n) = nphyh2d*min(1,ACC_UB(n))
      if (acc_ubflxs(n) /= 0) nphyh2d = nphyh2d+1
      ACC_UBFLXS(n) = nphyh2d*min(1,ACC_UBFLXS(n))
      if (acc_uice(n) /= 0) nphyh2d = nphyh2d+1
      ACC_UICE(n) = nphyh2d*min(1,ACC_UICE(n))
      if (acc_ustar(n) /= 0) nphyh2d = nphyh2d+1
      ACC_USTAR(n) = nphyh2d*min(1,ACC_USTAR(n))
      if (acc_ustar3(n) /= 0) nphyh2d = nphyh2d+1
      ACC_USTAR3(n) = nphyh2d*min(1,ACC_USTAR3(n))
      if (acc_ustokes(n) /= 0) nphyh2d = nphyh2d+1
      ACC_USTOKES(n) = nphyh2d*min(1,ACC_USTOKES(n))
      if (acc_vb(n) /= 0) nphyh2d = nphyh2d+1
      ACC_VB(n) = nphyh2d*min(1,ACC_VB(n))
      if (acc_vbflxs(n) /= 0) nphyh2d = nphyh2d+1
      ACC_VBFLXS(n) = nphyh2d*min(1,ACC_VBFLXS(n))
      if (acc_vice(n) /= 0) nphyh2d = nphyh2d+1
      ACC_VICE(n) = nphyh2d*min(1,ACC_VICE(n))
      if (acc_vstokes(n) /= 0) nphyh2d = nphyh2d+1
      ACC_VSTOKES(n) = nphyh2d*min(1,ACC_VSTOKES(n))
      if (acc_ztx(n) /= 0) nphyh2d = nphyh2d+1
      ACC_ZTX(n) = nphyh2d*min(1,ACC_ZTX(n))

      if (acc_bfsq(n) /= 0) nphylyr = nphylyr+1
      ACC_BFSQ(n) = nphylyr*min(1,ACC_BFSQ(n))
      if (acc_difdia(n) /= 0) nphylyr = nphylyr+1
      ACC_DIFDIA(n) = nphylyr*min(1,ACC_DIFDIA(n))
      if (acc_difvmo(n) /= 0) nphylyr = nphylyr+1
      ACC_DIFVMO(n) = nphylyr*min(1,ACC_DIFVMO(n))
      if (acc_difvho(n) /= 0) nphylyr = nphylyr+1
      ACC_DIFVHO(n) = nphylyr*min(1,ACC_DIFVHO(n))
      if (acc_difvso(n) /= 0) nphylyr = nphylyr+1
      ACC_DIFVSO(n) = nphylyr*min(1,ACC_DIFVSO(n))
      if (acc_difint(n) /= 0) nphylyr = nphylyr+1
      ACC_DIFINT(n) = nphylyr*min(1,ACC_DIFINT(n))
      if (acc_difiso(n) /= 0) nphylyr = nphylyr+1
      ACC_DIFISO(n) = nphylyr*min(1,ACC_DIFISO(n))
      if (acc_dp(n) /= 0) nphylyr = nphylyr+1
      ACC_DP(n) = nphylyr*min(1,ACC_DP(n))
      if (acc_dpu(n) /= 0) nphylyr = nphylyr+1
      ACC_DPU(n) = nphylyr*min(1,ACC_DPU(n))
      if (acc_dpv(n) /= 0) nphylyr = nphylyr+1
      ACC_DPV(n) = nphylyr*min(1,ACC_DPV(n))
      if (acc_dz(n) /= 0) nphylyr = nphylyr+1
      ACC_DZ(n) = nphylyr*min(1,ACC_DZ(n))
      if (acc_saln(n) /= 0) nphylyr = nphylyr+1
      ACC_SALN(n) = nphylyr*min(1,ACC_SALN(n))
      if (acc_temp(n) /= 0) nphylyr = nphylyr+1
      ACC_TEMP(n) = nphylyr*min(1,ACC_TEMP(n))
      if (acc_uflx(n) /= 0) nphylyr = nphylyr+1
      ACC_UFLX(n) = nphylyr*min(1,ACC_UFLX(n))
      if (acc_utflx(n) /= 0) nphylyr = nphylyr+1
      ACC_UTFLX(n) = nphylyr*min(1,ACC_UTFLX(n))
      if (acc_usflx(n) /= 0) nphylyr = nphylyr+1
      ACC_USFLX(n) = nphylyr*min(1,ACC_USFLX(n))
      if (acc_umfltd(n) /= 0) nphylyr = nphylyr+1
      ACC_UMFLTD(n) = nphylyr*min(1,ACC_UMFLTD(n))
      if (acc_umflsm(n) /= 0) nphylyr = nphylyr+1
      ACC_UMFLSM(n) = nphylyr*min(1,ACC_UMFLSM(n))
      if (acc_utfltd(n) /= 0) nphylyr = nphylyr+1
      ACC_UTFLTD(n) = nphylyr*min(1,ACC_UTFLTD(n))
      if (acc_utflsm(n) /= 0) nphylyr = nphylyr+1
      ACC_UTFLSM(n) = nphylyr*min(1,ACC_UTFLSM(n))
      if (acc_utflld(n) /= 0) nphylyr = nphylyr+1
      ACC_UTFLLD(n) = nphylyr*min(1,ACC_UTFLLD(n))
      if (acc_usfltd(n) /= 0) nphylyr = nphylyr+1
      ACC_USFLTD(n) = nphylyr*min(1,ACC_USFLTD(n))
      if (acc_usflsm(n) /= 0) nphylyr = nphylyr+1
      ACC_USFLSM(n) = nphylyr*min(1,ACC_USFLSM(n))
      if (acc_usflld(n) /= 0) nphylyr = nphylyr+1
      ACC_USFLLD(n) = nphylyr*min(1,ACC_USFLLD(n))
      if (acc_uvel(n) /= 0) nphylyr = nphylyr+1
      ACC_UVEL(n) = nphylyr*min(1,ACC_UVEL(n))
      if (acc_vflx(n) /= 0) nphylyr = nphylyr+1
      ACC_VFLX(n) = nphylyr*min(1,ACC_VFLX(n))
      if (acc_vtflx(n) /= 0) nphylyr = nphylyr+1
      ACC_VTFLX(n) = nphylyr*min(1,ACC_VTFLX(n))
      if (acc_vsflx(n) /= 0) nphylyr = nphylyr+1
      ACC_VSFLX(n) = nphylyr*min(1,ACC_VSFLX(n))
      if (acc_vmfltd(n) /= 0) nphylyr = nphylyr+1
      ACC_VMFLTD(n) = nphylyr*min(1,ACC_VMFLTD(n))
      if (acc_vmflsm(n) /= 0) nphylyr = nphylyr+1
      ACC_VMFLSM(n) = nphylyr*min(1,ACC_VMFLSM(n))
      if (acc_vtfltd(n) /= 0) nphylyr = nphylyr+1
      ACC_VTFLTD(n) = nphylyr*min(1,ACC_VTFLTD(n))
      if (acc_vtflsm(n) /= 0) nphylyr = nphylyr+1
      ACC_VTFLSM(n) = nphylyr*min(1,ACC_VTFLSM(n))
      if (acc_vtflld(n) /= 0) nphylyr = nphylyr+1
      ACC_VTFLLD(n) = nphylyr*min(1,ACC_VTFLLD(n))
      if (acc_vsfltd(n) /= 0) nphylyr = nphylyr+1
      ACC_VSFLTD(n) = nphylyr*min(1,ACC_VSFLTD(n))
      if (acc_vsflsm(n) /= 0) nphylyr = nphylyr+1
      ACC_VSFLSM(n) = nphylyr*min(1,ACC_VSFLSM(n))
      if (acc_vsflld(n) /= 0) nphylyr = nphylyr+1
      ACC_VSFLLD(n) = nphylyr*min(1,ACC_VSFLLD(n))
      if (acc_vvel(n) /= 0) nphylyr = nphylyr+1
      ACC_VVEL(n) = nphylyr*min(1,ACC_VVEL(n))
      if (acc_wflx(n) /= 0) nphylyr = nphylyr+1
      ACC_WFLX(n) = nphylyr*min(1,ACC_WFLX(n))
      if (acc_wflx2(n) /= 0) nphylyr = nphylyr+1
      ACC_WFLX2(n) = nphylyr*min(1,ACC_WFLX2(n))
      if (acc_avdsg(n) /= 0) nphylyr = nphylyr+1
      ACC_AVDSG(n) = nphylyr*min(1,ACC_AVDSG(n))
      if (acc_dpvor(n) /= 0) nphylyr = nphylyr+1
      ACC_DPVOR(n) = nphylyr*min(1,ACC_DPVOR(n))
      if (acc_tke(n) /= 0) nphylyr = nphylyr+1
      ACC_TKE(n) = nphylyr*min(1,ACC_TKE(n))
      if (acc_gls_psi(n) /= 0) nphylyr = nphylyr+1
      ACC_GLS_PSI(n) = nphylyr*min(1,ACC_GLS_PSI(n))

      if (acc_bfsqlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_BFSQLVL(n) = nphylvl*min(1,ACC_BFSQLVL(n))
      if (acc_difdialvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DIFDIALVL(n) = nphylvl*min(1,ACC_DIFDIALVL(n))
      if (acc_difvmolvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DIFVMOLVL(n) = nphylvl*min(1,ACC_DIFVMOLVL(n))
      if (acc_difvholvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DIFVHOLVL(n) = nphylvl*min(1,ACC_DIFVHOLVL(n))
      if (acc_difvsolvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DIFVSOLVL(n) = nphylvl*min(1,ACC_DIFVSOLVL(n))
      if (acc_difintlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DIFINTLVL(n) = nphylvl*min(1,ACC_DIFINTLVL(n))
      if (acc_difisolvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DIFISOLVL(n) = nphylvl*min(1,ACC_DIFISOLVL(n))
      if (acc_dzlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_DZLVL(n) = nphylvl*min(1,ACC_DZLVL(n))
      if (acc_salnlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_SALNLVL(n) = nphylvl*min(1,ACC_SALNLVL(n))
      if (acc_templvl(n) /= 0) nphylvl = nphylvl+1
      ACC_TEMPLVL(n) = nphylvl*min(1,ACC_TEMPLVL(n))
      if (acc_uflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UFLXLVL(n) = nphylvl*min(1,ACC_UFLXLVL(n))
      if (acc_uflxold(n) /= 0) nphylvl = nphylvl+1
      ACC_UFLXOLD(n) = nphylvl*min(1,ACC_UFLXOLD(n))
      if (acc_utflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UTFLXLVL(n) = nphylvl*min(1,ACC_UTFLXLVL(n))
      if (acc_usflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_USFLXLVL(n) = nphylvl*min(1,ACC_USFLXLVL(n))
      if (acc_umfltdlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UMFLTDLVL(n) = nphylvl*min(1,ACC_UMFLTDLVL(n))
      if (acc_umflsmlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UMFLSMLVL(n) = nphylvl*min(1,ACC_UMFLSMLVL(n))
      if (acc_utfltdlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UTFLTDLVL(n) = nphylvl*min(1,ACC_UTFLTDLVL(n))
      if (acc_utflsmlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UTFLSMLVL(n) = nphylvl*min(1,ACC_UTFLSMLVL(n))
      if (acc_utflldlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UTFLLDLVL(n) = nphylvl*min(1,ACC_UTFLLDLVL(n))
      if (acc_usfltdlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_USFLTDLVL(n) = nphylvl*min(1,ACC_USFLTDLVL(n))
      if (acc_usflsmlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_USFLSMLVL(n) = nphylvl*min(1,ACC_USFLSMLVL(n))
      if (acc_usflldlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_USFLLDLVL(n) = nphylvl*min(1,ACC_USFLLDLVL(n))
      if (acc_uvellvl(n) /= 0) nphylvl = nphylvl+1
      ACC_UVELLVL(n) = nphylvl*min(1,ACC_UVELLVL(n))
      if (acc_vflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VFLXLVL(n) = nphylvl*min(1,ACC_VFLXLVL(n))
      if (acc_vflxold(n) /= 0) nphylvl = nphylvl+1
      ACC_VFLXOLD(n) = nphylvl*min(1,ACC_VFLXOLD(n))
      if (acc_vtflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VTFLXLVL(n) = nphylvl*min(1,ACC_VTFLXLVL(n))
      if (acc_vsflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VSFLXLVL(n) = nphylvl*min(1,ACC_VSFLXLVL(n))
      if (acc_vmfltdlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VMFLTDLVL(n) = nphylvl*min(1,ACC_VMFLTDLVL(n))
      if (acc_vmflsmlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VMFLSMLVL(n) = nphylvl*min(1,ACC_VMFLSMLVL(n))
      if (acc_vtfltdlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VTFLTDLVL(n) = nphylvl*min(1,ACC_VTFLTDLVL(n))
      if (acc_vtflsmlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VTFLSMLVL(n) = nphylvl*min(1,ACC_VTFLSMLVL(n))
      if (acc_vtflldlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VTFLLDLVL(n) = nphylvl*min(1,ACC_VTFLLDLVL(n))
      if (acc_vsfltdlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VSFLTDLVL(n) = nphylvl*min(1,ACC_VSFLTDLVL(n))
      if (acc_vsflsmlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VSFLSMLVL(n) = nphylvl*min(1,ACC_VSFLSMLVL(n))
      if (acc_vsflldlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VSFLLDLVL(n) = nphylvl*min(1,ACC_VSFLLDLVL(n))
      if (acc_vvellvl(n) /= 0) nphylvl = nphylvl+1
      ACC_VVELLVL(n) = nphylvl*min(1,ACC_VVELLVL(n))
      if (acc_wflxlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_WFLXLVL(n) = nphylvl*min(1,ACC_WFLXLVL(n))
      if (acc_wflx2lvl(n) /= 0) nphylvl = nphylvl+1
      ACC_WFLX2LVL(n) = nphylvl*min(1,ACC_WFLX2LVL(n))
      if (acc_pvlvl(n) /= 0) nphylvl = nphylvl+1
      ACC_PVLVL(n) = nphylvl*min(1,ACC_PVLVL(n))
      if (acc_tkelvl(n) /= 0) nphylvl = nphylvl+1
      ACC_TKELVL(n) = nphylvl*min(1,ACC_TKELVL(n))
      if (acc_gls_psilvl(n) /= 0) nphylvl = nphylvl+1
      ACC_GLS_PSILVL(n) = nphylvl*min(1,ACC_GLS_PSILVL(n))
      if (acc_idlagelvl(n) /= 0) nphylvl = nphylvl+1
      ACC_IDLAGELVL(n) = nphylvl*min(1,ACC_IDLAGELVL(n))

      ! End loop over io groups
    end do

    ! Assign buffer positions for utility fields
    ACC_UTILH2D = 0
    nphyh2d = nphyh2d+1
    ACC_UTILH2D(1) = nphyh2d

    ACC_UTILLYR = 0
    nphylyr = nphylyr+1
    ACC_UTILLYR(1) = nphylyr

    ACC_UTILLVL = 0
    nphylvl = nphylvl+1
    ACC_UTILLVL(1) = nphylvl

    ! Allocate buffers
    istatsum = 0
    istat = 0
    if (nphyh2d /= 0) &
         allocate(phyh2d(imn:imx,jmn:jmx,nphyh2d),stat = istat)
    istatsum = istatsum+istat
    if (nphylyr /= 0) &
         allocate(phylyr(imn:imx,jmn:jmx,kdm,nphylyr),stat = istat)
    istatsum = istatsum+istat
    if (nphylvl /= 0) &
         allocate(phylvl(imn:imx,jmn:jmx,ddm,nphylvl),stat = istat)
    istatsum = istatsum+istat
    if (istatsum /= 0) then
      write (lp,*) 'Cannot allocate enough memory!'
      call xchalt('(diaini)')
      stop '(diaini)'
    end if

    ! initialisation of h2h, lyr and lvl fields
    do n = 1,nphy
      call inifld(n)
    end do

    ! Load bathymetry into module mod_dia (used for vertical
    ! interpolation in BLOM and HAMOCC)
    nstepinday = nstep_in_day
    !$omp parallel do private(l,i)
    do j = 1,jj+1
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbath(i,j) = depths(i,j)
        end do
      end do
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          ubath(i,j) = min(depths(i,j),depths(i-1,j))
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          vbath(i,j) = min(depths(i,j),depths(i,j-1))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine diaini


  subroutine diaacc(m,n,mm,nn,k1m,k1n)
  !---------------------------------------------------------------
  ! accumulate diagnostic variables
  !---------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    integer :: i,j,k,l,km,kup,iogrp
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2, &
         ipsw,ipse,ipnw,ipne
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts, &
         wghtsflx
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: uvel,vvel, &
         avdsg_p,dpvor_p,pv_p,dummy
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         dpml,sbot,tbot,dps,t20d
    real :: dsig,q,zup,zlo,tup,tlo

    ! Increase counter
    do iogrp = 1,nphy
      nacc_phy(iogrp) = nacc_phy(iogrp)+1
    end do

    ! Define auxillary variables

    if (sum(acc_uice(1:nphy)+acc_vice(1:nphy)) /= 0) then
      call xctilr(hicem, 1,1, 1,1, halo_ps)
      call xctilr(ficem, 1,1, 1,1, halo_ps)
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            util1(i,j) = hicem(i-1,j)*ficem(i-1,j)+hicem(i,j)* &
                 ficem(i,j)
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            util3(i,j) = hicem(i,j-1)*ficem(i,j-1)+hicem(i,j)* &
                 ficem(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (sum(acc_mld(1:nphy)+acc_maxmld(1:nphy)) /= 0) then
      select case (vcoord_tag)
        case (vcoord_isopyc_bulkml)
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                dpml(i,j) = dp(i,j,1+mm)+dp(i,j,2+mm)
              end do
            end do
          end do
          !$omp end parallel do
        case default
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                dpml(i,j) = OBLdepth(i,j)*onem
              end do
            end do
          end do
          !$omp end parallel do
      end select
    end if

    if (sum(acc_uvel(1:nphy)+acc_uvellvl(1:nphy)) /= 0) then
      !$omp parallel do private(k,km,l,i)
      do j = 1,jj
        do k = 1,kk
          km = k+mm
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              uvel(i,j,k) = u(i,j,km)+ub(i,j,m)
            end do
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (sum(acc_vvel(1:nphy)+acc_vvellvl(1:nphy)) /= 0) then
      !$omp parallel do private(k,km,l,i)
      do j = 1,jj+1
        do k = 1,kk
          km = k+mm
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vvel(i,j,k) = v(i,j,km)+vb(i,j,m)
            end do
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (sum(acc_avdsg(1:nphy)+acc_pvlvl(1:nphy)) /= 0) then
      !$omp parallel do private(l,i,k,km,dsig)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            ipsw(i,j) = min(1,iu(i  ,j)+iv(i,j  ))
            ipse(i,j) = min(1,iu(i+1,j)+iv(i,j  ))
            ipnw(i,j) = min(1,iu(i  ,j)+iv(i,j+1))
            ipne(i,j) = min(1,iu(i+1,j)+iv(i,j+1))
          end do
        end do
        do k = 1,kk
          km = k+mm
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k <= 2) then
                dsig = max(0.,sigma(i,j,2+mm)-sigma(i,j,1+mm))
              else if (k == kk) then
                dsig = max(0.,sigma(i,j,km) &
                      -max(sigma(i,j,km-1),sigma(i,j,2+mm)))
              else
                dsig = .5*max(0.,sigma(i,j,km+1) &
                         -max(sigma(i,j,km-1),sigma(i,j,2+mm)))
              end if
              avdsg_p(i,j,k) = (absvor(i  ,j  ,k)*ipsw(i,j) &
                   +absvor(i+1,j  ,k)*ipse(i,j) &
                   +absvor(i  ,j+1,k)*ipnw(i,j) &
                   +absvor(i+1,j+1,k)*ipne(i,j))*dsig
              dpvor_p(i,j,k) = dpvor(i  ,j  ,k)*ipsw(i,j) &
                   +dpvor(i+1,j  ,k)*ipse(i,j) &
                   +dpvor(i  ,j+1,k)*ipnw(i,j) &
                   +dpvor(i+1,j+1,k)*ipne(i,j)
              pv_p(i,j,k) = avdsg_p(i,j,k)/dpvor_p(i,j,k)
            end do
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (sum(acc_sbot(1:nphy)+acc_tbot(1:nphy)) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            sbot(i,j) = 0.
            tbot(i,j) = 0.
            dps(i,j) = 0.
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(k,km,l,i,q)
      do j = 1,jj
        do k = 1,kk
          km = k+mm
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              q = max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpbot,p(i,j,k)))
              sbot(i,j) = sbot(i,j)+saln(i,j,km)*q
              tbot(i,j) = tbot(i,j)+temp(i,j,km)*q
              dps(i,j) = dps(i,j)+q
            end do
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(l,i,q)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (dps(i,j) > onemm) then
              q = 1./dps(i,j)
              sbot(i,j) = sbot(i,j)*q
              tbot(i,j) = tbot(i,j)*q
            else
              sbot(i,j) = saln(i,j,1+mm)
              tbot(i,j) = temp(i,j,1+mm)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    end if

    if (sum(acc_t20d(1:nphy)) /= 0) then
      !$omp parallel do private(l,i,k,km,kup,zup,zlo,tup,tlo)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            k = 1
            km = k+mm
            do
              if (dp(i,j,km) > onecm) then
                if (temp(i,j,km) > 20.) then
                  kup = k
                else
                  exit
                end if
              end if
              k = k+1
              km = k+mm
              if (k > kk) exit
            end do
            if     (k == 1) then
              t20d(i,j) = 0.
            else if (k > kk) then
              t20d(i,j) = z(i,j,kk+1)
            else
              zup = z(i,j,kup)+.5*dz(i,j,kup)
              zlo = z(i,j,k  )+.5*dz(i,j,k  )
              tup = temp(i,j,kup+mm)
              tlo = min(temp(i,j,km),tup-epsilp)
              t20d(i,j) = (zup*(tlo-20.)+zlo*(20.-tup))/(tlo-tup)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    end if

    !---------------------------------------------------------------
    ! accumulate 2d diagnostic variables
    !---------------------------------------------------------------

    ! u-component of barotropic velocity [m/s]
    call acch2d(ACC_UB,ub(1-nbdy,1-nbdy,m),dummy,0,'u')

    ! u-component of barotropic mass flux [kg*m/s^3]
    call acch2d(ACC_UBFLXS,ubflxs(1-nbdy,1-nbdy,n),dummy,0,'u')

    ! u-component of wind stress [N/m^2]
    call acch2d(ACC_ZTX,ztx,dummy,0,'u')

    ! u-component of momentum flux received by the ocean [kg/m/s^2]
    call acch2d(ACC_TAUX,taux,dummy,0,'u')

    ! weighted u-component of ice velocity [m^2/s]
    call acch2d(ACC_UICE,uicem,util1,1,'u')

    ! u-component of surface Stokes drift [m/s]
    call acch2d(ACC_USTOKES,ustokes,dummy,0,'u')

    ! v-component of barotropic velocity [m/s]
    call acch2d(ACC_VB,vb(1-nbdy,1-nbdy,m),dummy,0,'v')

    ! v-component of barotropic mass flux [kg*m/s^3]
    call acch2d(ACC_VBFLXS,vbflxs(1-nbdy,1-nbdy,n),dummy,0,'v')

    ! v-component of wind stress [N/m^2]
    call acch2d(ACC_MTY,mty,dummy,0,'v')

    ! v-component of momentum flux received by the ocean [kg/m/s^2]
    call acch2d(ACC_TAUY,tauy,dummy,0,'v')

    ! weighted v-component of ice velocity [m^2/s]
    call acch2d(ACC_VICE,vicem,util3,1,'v')

    ! v-component of surface Stokes drift [m/s]
    call acch2d(ACC_VSTOKES,vstokes,dummy,0,'v')

    ! surface pressure [kg/m/s^2]
    call acch2d(ACC_PSRF,p(1-nbdy,1-nbdy,1),dummy,0,'p')

    ! bottom pressure [kg/m/s^2]
    call acch2d(ACC_PBOT,p(1-nbdy,1-nbdy,kk+1),dummy,0,'p')

    ! sea level height [m]
    call acch2d(ACC_SEALV,z(1-nbdy,1-nbdy,1),dummy,0,'p')

    ! mixed layer density (sigma units)
    call acch2d(ACC_SIGMX,sigma(1-nbdy,1-nbdy,k1m),dummy,0,'p')

    ! weighted ice thickness [m]
    call acch2d(ACC_HICE,hicem,ficem,1,'p')

    ! weighted snow thickness [m]
    call acch2d(ACC_HSNW,hsnwm,ficem,1,'p')

    ! fractional ice cover
    call acch2d(ACC_FICE,ficem,dummy,0,'p')

    ! ice volume in u-points [m]
    call acch2d(ACC_IVOLU,util1,dummy,0,'u')

    ! ice volume in v-points [m]
    call acch2d(ACC_IVOLV,util3,dummy,0,'v')

    ! surface temperature [K]
    call acch2d(ACC_TSRF,tsrfm,dummy,0,'p')

    ! ice temperature [K]
    call acch2d(ACC_TICE,ticem,dummy,0,'p')

    ! short wave heat flux [W/m^2]
    call acch2d(ACC_SWA,swa,dummy,0,'p')

    ! non-solar heat flux [W/m^2]
    call acch2d(ACC_NSF,nsf,dummy,0,'p')

    ! heat flux due to material enthalpy flux [W/m^2]
    call acch2d(ACC_HMAT,hmat,dummy,0,'p')

    ! heat flux due to melting/freezing [W/m^2]
    call acch2d(ACC_HMLTFZ,hmltfz,dummy,0,'p')

    ! derivative of non-solar heat flux by surface temperature [W/m^2/K]
    if (allocated(dfl)) call acch2d(ACC_DFL,dfl,dummy,0,'p')

    ! liquid precipitation [mm/day]
    call acch2d(ACC_LIP,lip,dummy,0,'p')

    ! solid precipitation [mm/day]
    call acch2d(ACC_SOP,sop,dummy,0,'p')

    ! evaporation [mm/day]
    call acch2d(ACC_EVA,eva,dummy,0,'p')

    ! fresh water flux due to melting/freezing [kg/m^2/s]
    call acch2d(ACC_FMLTFZ,fmltfz,dummy,0,'p')

    ! salt flux [kg/m^2/s]
    call acch2d(ACC_SFL,sfl,dummy,0,'p')

    ! albedo
    if (allocated(alb)) call acch2d(ACC_ALB,alb,dummy,0,'p')

    ! liquid runoff [kg m-2 s-1]
    call acch2d(ACC_RNFFLX,rnf,dummy,0,'p')

    ! frozen runoff [kg m-2 s-1]
    call acch2d(ACC_RFIFLX,rfi,dummy,0,'p')

    ! friction velocity [m s-1]
    call acch2d(ACC_USTAR,ustar,dummy,0,'p')

    ! friction velocity cubed [m3 s-3]
    call acch2d(ACC_USTAR3,ustar3,dummy,0,'p')

    ! mixed layer integrated inertial kinetic energy tendency [m3 s-3]
    call acch2d(ACC_IDKEDT,idkedt,dummy,0,'p')

    ! absolute wind speed [m s-1]
    call acch2d(ACC_ABSWND,abswnd,dummy,0,'p')

    ! mixed layer TKE tendency related to friction velocity [m3 s-3]
    call acch2d(ACC_MTKEUS,mtkeus,dummy,0,'p')

    ! mixed layer TKE tendency related to near inertial motions [m3 s-3]
    call acch2d(ACC_MTKENI,mtkeni,dummy,0,'p')

    ! mixed layer TKE tendency related to buoyancy forcing [m3 s-3]
    call acch2d(ACC_MTKEBF,mtkebf,dummy,0,'p')

    ! mixed layer TKE tendency related to eddy restratification [m3 s-3]
    call acch2d(ACC_MTKERS,mtkers,dummy,0,'p')

    ! mixed layer TKE tendency related to potential energy change [m3 s-3]
    call acch2d(ACC_MTKEPE,mtkepe,dummy,0,'p')

    ! mixed layer TKE tendency related to kinetic energy change [m3 s-3]
    call acch2d(ACC_MTKEKE,mtkeke,dummy,0,'p')

    ! Langmuir enhancement factor []
    call acch2d(ACC_LAMULT,lamult,dummy,0,'p')

    ! Surface layer averaged Langmuir number []
    call acch2d(ACC_LASL,lasl,dummy,0,'p')

    ! sea surface salinity [g kg-1]
    call acch2d(ACC_SSS,saln(1-nbdy,1-nbdy,k1m),dummy,0,'p')

    ! sea surface temperature [degC]
    call acch2d(ACC_SST,temp(1-nbdy,1-nbdy,k1m),dummy,0,'p')

    ! bottom salinity [g kg-1]
    call acch2d(ACC_SBOT,sbot,dummy,0,'p')

    ! bottom temperature [degC]
    call acch2d(ACC_TBOT,tbot,dummy,0,'p')

    ! mixed layer pressure thickness [kg/m/s^2]
    call acch2d(ACC_MLD,dpml,dummy,0,'p')

    ! mixed layer thickness using "sigma-t" criterion [m]
    call acch2d(ACC_MLTS,mlts,dummy,0,'p')

    ! 20C isoterm depth [m]
    call acch2d(ACC_T20D,t20d,dummy,0,'p')

    ! heat flux given by the ocean [W/m^2]
    call acch2d(ACC_SURFLX,surflx,dummy,0,'p')

    ! salt flux given by the ocean [g/m^2/s]
    call acch2d(ACC_SALFLX,salflx,dummy,0,'p')

    ! restoring heat flux received by the ocean [W/m^2]
    call acch2d(ACC_SURRLX,surrlx,dummy,0,'p')

    ! restoring salt flux received by the ocean [g/m^2/s]
    call acch2d(ACC_SALRLX,salrlx,dummy,0,'p')

    ! brine flux received by the ocean [g/m^2/s]
    call acch2d(ACC_BRNFLX,brnflx,dummy,0,'p')

    ! brine plume pressure depth [kg/m/s^2]
    call acch2d(ACC_BRNPD,pbrnda,dummy,0,'p')

    !---------------------------------------------------------------
    ! store minimum or maximum of 2d diagnostic variables
    !---------------------------------------------------------------

    ! maximum mixed layer pressure thickness [kg/m/s^2]
    call maxh2d(ACC_MAXMLD,dpml,'p')

    ! minimum mixed layer thickness using "sigma-t" criterion [m]
    call minh2d(ACC_MLTSMN,mlts,'p')

    ! maximum mixed layer thickness using "sigma-t" criterion [m]
    call maxh2d(ACC_MLTSMX,mlts,'p')

    !---------------------------------------------------------------
    ! store squared of 2d diagnostic variables
    !---------------------------------------------------------------

    ! mixed layer thickness squared using "sigma-t" criterion [m^2]
    call sqh2d(ACC_MLTSSQ,mlts,'p')

    ! sea level height squared [m^2]
    call sqh2d(ACC_SLVSQ,z(1-nbdy,1-nbdy,1),'p')

    ! sea surface salinity squared [g^2/kg^2]
    call sqh2d(ACC_SSSSQ,saln(1-nbdy,1-nbdy,k1m),'p')

    ! sea surface temperature squared [degC2]
    call sqh2d(ACC_SSTSQ,temp(1-nbdy,1-nbdy,k1m),'p')

    !---------------------------------------------------------------
    ! accumulate 3d diagnostic variables
    !---------------------------------------------------------------

    ! weighted u-component of total velocity [kg/s^3]
    call acclyr(ACC_UVEL,uvel,dpu(1-nbdy,1-nbdy,k1m),1,'u')

    ! layer pressure thickness at u-point [kg/m/s^2]
    call acclyr(ACC_DPU,dpu(1-nbdy,1-nbdy,k1m),dummy,0,'u')

    ! u-component of mass flux [kg*m/s^2]
    call acclyr(ACC_UFLX,uflx(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of heat flux [K*kg*m/s^2]
    call acclyr(ACC_UTFLX,utflx(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of salt flux [g*m/s^2]
    call acclyr(ACC_USFLX,usflx(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of mass flux due to thickness diffusion [kg*m/s^2]
    call acclyr(ACC_UMFLTD,umfltd(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of mass flux due to submesoscale transport [kg*m/s^2]
    call acclyr(ACC_UMFLSM,umflsm(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of heat flux due to thickness diffusion [K*kg*m/s^2]
    call acclyr(ACC_UTFLTD,utfltd(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of heat flux due to submesoscale transport [K*kg*m/s^2]
    call acclyr(ACC_UTFLSM,utflsm(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of heat flux due to lateral diffusion [K*kg*m/s^2]
    call acclyr(ACC_UTFLLD,utflld(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of salt flux due to thickness diffusion [g*m/s^2]
    call acclyr(ACC_USFLTD,usfltd(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of salt flux due to submesoscale transport [g*m/s^2]
    call acclyr(ACC_USFLSM,usflsm(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! u-component of salt flux due to lateral diffusion [g*m/s^2]
    call acclyr(ACC_USFLLD,usflld(1-nbdy,1-nbdy,k1n),dummy,0,'u')

    ! weighted v-component of total velocity [kg/s^3]
    call acclyr(ACC_VVEL,vvel,dpv(1-nbdy,1-nbdy,k1m),1,'v')

    ! layer pressure thickness at v-point [kg/m/s^2]
    call acclyr(ACC_DPV,dpv(1-nbdy,1-nbdy,k1m),dummy,0,'v')

    ! v-component of mass flux [kg*m/s^2]
    call acclyr(ACC_VFLX,vflx(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of heat flux [K*kg*m/s^2]
    call acclyr(ACC_VTFLX,vtflx(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of salt flux [g*m/s^2]
    call acclyr(ACC_VSFLX,vsflx(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of mass flux due to thickness diffusion [kg*m/s^2]
    call acclyr(ACC_VMFLTD,vmfltd(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of mass flux due to submesoscale transport [kg*m/s^2]
    call acclyr(ACC_VMFLSM,vmflsm(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of heat flux due to thickness diffusion [K*kg*m/s^2]
    call acclyr(ACC_VTFLTD,vtfltd(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of heat flux due to submesoscale transport [K*kg*m/s^2]
    call acclyr(ACC_VTFLSM,vtflsm(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of heat flux due to lateral diffusion [K*kg*m/s^2]
    call acclyr(ACC_VTFLLD,vtflld(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of salt flux due to thickness diffusion [g*m/s^2]
    call acclyr(ACC_VSFLTD,vsfltd(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of salt flux due to submesoscale transport [g*m/s^2]
    call acclyr(ACC_VSFLSM,vsflsm(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! v-component of salt flux due to lateral diffusion [g*m/s^2]
    call acclyr(ACC_VSFLLD,vsflld(1-nbdy,1-nbdy,k1n),dummy,0,'v')

    ! weighted salinity [g/m/s^2]
    call acclyr(ACC_SALN,saln(1-nbdy,1-nbdy,k1m), dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted temperature [degC*kg/m/s^2]
    call acclyr(ACC_TEMP,temp(1-nbdy,1-nbdy,k1m), dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! layer pressure thickness [kg/m/s^2]
    call acclyr(ACC_DP,dp(1-nbdy,1-nbdy,k1m),dummy,0,'p')

    ! layer thickness [m]
    call acclyr(ACC_DZ,dz,dummy,0,'p')

    ! weighted buoyancy frequency squared [kg/m/s^4]
    call acclyr(ACC_BFSQ,bfsql,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted layer interface diffusivity [m*kg/s^3]
    call acclyr(ACC_DIFINT,difint,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted isopycnal diffusivity [m*kg/s^3]
    call acclyr(ACC_DIFISO,difiso,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted vertical diffusivity (vcoord == 'isopyc_bulkml')
    ! [m*kg/s^3]
    call acclyr(ACC_DIFDIA,difdia,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted vertical momentum diffusivity (vcoord /= 'isopyc_bulkml')
    ! [m*kg/s^3]
    call accily(ACC_DIFVMO,Kvisc_m,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted vertical heat diffusivity (vcoord /= 'isopyc_bulkml')
    ! [m*kg/s^3]
    call accily(ACC_DIFVHO,Kdiff_t,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! weighted vertical salt diffusivity (vcoord /= 'isopyc_bulkml')
    ! [m*kg/s^3]
    call accily(ACC_DIFVSO,Kdiff_s,dp(1-nbdy,1-nbdy,k1m),1,'p')

    ! absolute vorticity multiplied with potential density difference
    ! over layer [kg/m^3/s]
    call acclyr(ACC_AVDSG,avdsg_p,dummy,0,'p')

    ! layer pressure thickness used in vorticity computation [kg/m/s^2]
    call acclyr(ACC_DPVOR,dpvor_p,dummy,0,'p')

    if (use_TRC .and. use_TKE) then
      ! weighted tke [m*kg/s^4]
      call acclyr(ACC_TKE,trc(1-nbdy,1-nbdy,k1m,itrtke), &
           dp(1-nbdy,1-nbdy,k1m),1,'p')

      ! weighted gls_psi [m*kg/s^5]
      call acclyr(ACC_GLS_PSI,trc(1-nbdy,1-nbdy,k1m,itrgls), &
           dp(1-nbdy,1-nbdy,k1m),1,'p')

    end if
    !---------------------------------------------------------------
    ! accumulate 3d diagnostic variables on Levitus levels
    !---------------------------------------------------------------

    do iogrp = 1,nphy
      if (acc_wflxlvl(iogrp)+acc_wflx2lvl(iogrp) /= 0) then
        !$omp parallel do private(k,l,i)
        do j = 1,jj+1
          do k = 1,ddm
            do l = 1,isu(j)
              do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                phylvl(i,j,k,ACC_UFLXOLD(iogrp))= &
                     phylvl(i,j,k,ACC_UFLXLVL(iogrp))
              end do
            end do
            do l = 1,isv(j)
              do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                phylvl(i,j,k,ACC_VFLXOLD(iogrp))= &
                     phylvl(i,j,k,ACC_VFLXLVL(iogrp))
              end do
            end do
          end do
        end do
        !$omp end parallel do
      end if
    end do

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      if (sum(ACC_UFLXLVL  (1:nphy) &
             +ACC_UTFLXLVL (1:nphy)+ACC_USFLXLVL (1:nphy) &
             +ACC_UMFLTDLVL(1:nphy)+ACC_UMFLSMLVL(1:nphy) &
             +ACC_UTFLTDLVL(1:nphy)+ACC_UTFLSMLVL(1:nphy) &
             +ACC_UTFLLDLVL(1:nphy)+ACC_USFLTDLVL(1:nphy) &
             +ACC_USFLSMLVL(1:nphy)+ACC_USFLLDLVL(1:nphy)) /= 0) then
        do k = 1,kk
          call diazlv('u',k,mm,nn,ind1,ind2,wghts,wghtsflx)

          ! u-component of mass flux [kg*m/s^2]
          call acclvl(ACC_UFLXLVL,uflx(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux [K*kg*m/s^2]
          call acclvl(ACC_UTFLXLVL,utflx(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux [g*m/s^2]
          call acclvl(ACC_USFLXLVL,usflx(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of mass flux due to thickness diffusion [kg*m/s^2]
          call acclvl(ACC_UMFLTDLVL,umfltd(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of mass flux due to submesoscale transport [kg*m/s^2]
          call acclvl(ACC_UMFLSMLVL,umflsm(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux due to thickness diffusion [K*kg*m/s^2]
          call acclvl(ACC_UTFLTDLVL,utfltd(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux due to submesoscale transport [K*kg*m/s^2]
          call acclvl(ACC_UTFLSMLVL,utflsm(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux due to lateral diffusion [K*kg*m/s^2]
          call acclvl(ACC_UTFLLDLVL,utflld(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux due to thickness diffusion [g*m/s^2]
          call acclvl(ACC_USFLTDLVL,usfltd(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux due to submesoscale transport [g*m/s^2]
          call acclvl(ACC_USFLSMLVL,usflsm(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux due to lateral diffusion [g*m/s^2]
          call acclvl(ACC_USFLLDLVL,usflld(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)
        end do
      end if
    else
      if (sum(ACC_UVELLVL  (1:nphy)+ACC_UFLXLVL  (1:nphy) &
             +ACC_UTFLXLVL (1:nphy)+ACC_USFLXLVL (1:nphy) &
             +ACC_UMFLTDLVL(1:nphy)+ACC_UMFLSMLVL(1:nphy) &
             +ACC_UTFLTDLVL(1:nphy)+ACC_UTFLSMLVL(1:nphy) &
             +ACC_UTFLLDLVL(1:nphy)+ACC_USFLTDLVL(1:nphy) &
             +ACC_USFLSMLVL(1:nphy)+ACC_USFLLDLVL(1:nphy)) /= 0) then
        do k = 1,kk
          call diazlv('u',k,mm,nn,ind1,ind2,wghts,wghtsflx)

          ! u-component of total velocity [m/s]
          call acclvl(ACC_UVELLVL,uvel,'u',k,ind1,ind2,wghts)

          ! u-component of mass flux [kg*m/s^2]
          call acclvl(ACC_UFLXLVL,uflx(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux [K*kg*m/s^2]
          call acclvl(ACC_UTFLXLVL,utflx(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux [g*m/s^2]
          call acclvl(ACC_USFLXLVL,usflx(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of mass flux due to thickness diffusion [kg*m/s^2]
          call acclvl(ACC_UMFLTDLVL,umfltd(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of mass flux due to submesoscale transport [kg*m/s^2]
          call acclvl(ACC_UMFLSMLVL,umflsm(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux due to thickness diffusion [K*kg*m/s^2]
          call acclvl(ACC_UTFLTDLVL,utfltd(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux due to submesoscale transport [K*kg*m/s^2]
          call acclvl(ACC_UTFLSMLVL,utflsm(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of heat flux due to lateral diffusion [K*kg*m/s^2]
          call acclvl(ACC_UTFLLDLVL,utflld(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux due to thickness diffusion [g*m/s^2]
          call acclvl(ACC_USFLTDLVL,usfltd(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux due to submesoscale transport [g*m/s^2]
          call acclvl(ACC_USFLSMLVL,usflsm(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)

          ! u-component of salt flux due to lateral diffusion [g*m/s^2]
          call acclvl(ACC_USFLLDLVL,usflld(1-nbdy,1-nbdy,k1n), &
               'u',k,ind1,ind2,wghtsflx)
        end do
      end if
    end if

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      if (sum(ACC_VFLXLVL  (1:nphy) &
             +ACC_VTFLXLVL (1:nphy)+ACC_VSFLXLVL (1:nphy) &
             +ACC_VMFLTDLVL(1:nphy)+ACC_VMFLSMLVL(1:nphy) &
             +ACC_VTFLTDLVL(1:nphy)+ACC_VTFLSMLVL(1:nphy) &
             +ACC_VTFLLDLVL(1:nphy)+ACC_VSFLTDLVL(1:nphy) &
             +ACC_VSFLSMLVL(1:nphy)+ACC_VSFLLDLVL(1:nphy)) /= 0) then
        do k = 1,kk
          call diazlv('v',k,mm,nn,ind1,ind2,wghts,wghtsflx)

          ! v-component of mass flux [kg*m/s^2]
          call acclvl(ACC_VFLXLVL,vflx(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux [K*kg*m/s^2]
          call acclvl(ACC_VTFLXLVL,vtflx(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux [g*m/s^2]
          call acclvl(ACC_VSFLXLVL,vsflx(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of mass flux due to thickness diffusion [kg*m/s^2]
          call acclvl(ACC_VMFLTDLVL,vmfltd(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of mass flux due to submesoscale transport [kg*m/s^2]
          call acclvl(ACC_VMFLSMLVL,vmflsm(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux due to thickness diffusion [K*kg*m/s^2]
          call acclvl(ACC_VTFLTDLVL,vtfltd(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux due to submesoscale transport [K*kg*m/s^2]
          call acclvl(ACC_VTFLSMLVL,vtflsm(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux due to lateral diffusion [K*kg*m/s^2]
          call acclvl(ACC_VTFLLDLVL,vtflld(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux due to thickness diffusion [g*m/s^2]
          call acclvl(ACC_VSFLTDLVL,vsfltd(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux due to submesoscale transport [g*m/s^2]
          call acclvl(ACC_VSFLSMLVL,vsflsm(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux due to lateral diffusion [g*m/s^2]
          call acclvl(ACC_VSFLLDLVL,vsflld(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)
        end do
      end if
    else
      if (sum(ACC_VVELLVL  (1:nphy)+ACC_VFLXLVL  (1:nphy) &
             +ACC_VTFLXLVL (1:nphy)+ACC_VSFLXLVL (1:nphy) &
             +ACC_VMFLTDLVL(1:nphy)+ACC_VMFLSMLVL(1:nphy) &
             +ACC_VTFLTDLVL(1:nphy)+ACC_VTFLSMLVL(1:nphy) &
             +ACC_VTFLLDLVL(1:nphy)+ACC_VSFLTDLVL(1:nphy) &
             +ACC_VSFLSMLVL(1:nphy)+ACC_VSFLLDLVL(1:nphy)) /= 0) then
        do k = 1,kk
          call diazlv('v',k,mm,nn,ind1,ind2,wghts,wghtsflx)

          ! v-component of total velocity [m/s]
          call acclvl(ACC_VVELLVL,vvel,'v',k,ind1,ind2,wghts)

          ! v-component of mass flux [kg*m/s^2]
          call acclvl(ACC_VFLXLVL,vflx(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux [K*kg*m/s^2]
          call acclvl(ACC_VTFLXLVL,vtflx(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux [g*m/s^2]
          call acclvl(ACC_VSFLXLVL,vsflx(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of mass flux due to thickness diffusion [kg*m/s^2]
          call acclvl(ACC_VMFLTDLVL,vmfltd(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of mass flux due to submesoscale transport [kg*m/s^2]
          call acclvl(ACC_VMFLSMLVL,vmflsm(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux due to thickness diffusion [K*kg*m/s^2]
          call acclvl(ACC_VTFLTDLVL,vtfltd(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux due to submesoscale transport [K*kg*m/s^2]
          call acclvl(ACC_VTFLSMLVL,vtflsm(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of heat flux due to lateral diffusion [K*kg*m/s^2]
          call acclvl(ACC_VTFLLDLVL,vtflld(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux due to thickness diffusion [g*m/s^2]
          call acclvl(ACC_VSFLTDLVL,vsfltd(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux due to submesoscale transport [g*m/s^2]
          call acclvl(ACC_VSFLSMLVL,vsflsm(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)

          ! v-component of salt flux due to lateral diffusion [g*m/s^2]
          call acclvl(ACC_VSFLLDLVL,vsflld(1-nbdy,1-nbdy,k1n), &
               'v',k,ind1,ind2,wghtsflx)
        end do
      end if
    end if

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      if (sum(ACC_BFSQLVL(1:nphy)   +ACC_DIFDIALVL(1:nphy) &
             +ACC_DIFVMOLVL(1:nphy) +ACC_DIFVHOLVL(1:nphy) &
             +ACC_DIFVSOLVL(1:nphy) +ACC_DIFINTLVL(1:nphy) &
             +ACC_DIFISOLVL(1:nphy) +ACC_TKELVL(1:nphy) &
             +ACC_GLS_PSILVL(1:nphy)+ACC_PVLVL(1:nphy) &
             +ACC_DZLVL(1:nphy)     ) /= 0) then
        do k = 1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)

          ! buoyancy frequency squared [1/s^2]
          call acclvl(ACC_BFSQLVL,bfsql,'p',k,ind1,ind2,wghts)

          ! layer interface diffusivity [m^2/s]
          call acclvl(ACC_DIFINTLVL,difint,'p',k,ind1,ind2,wghts)

          ! isopycnal diffusivity [m^2/s]
          call acclvl(ACC_DIFISOLVL,difiso,'p',k,ind1,ind2,wghts)

          ! vertical diffusivity (vcoord == 'isopyc_bulkml') [m^2/s]
          call acclvl(ACC_DIFDIALVL,difdia,'p',k,ind1,ind2,wghts)

          ! vertical momentum diffusivity (vcoord /= 'isopyc_bulkml') [m^2/s]
          call accilv(ACC_DIFVMOLVL,Kvisc_m,'p',k,ind1,ind2,wghts)

          ! vertical heat diffusivity (vcoord /= 'isopyc_bulkml') [m^2/s]
          call accilv(ACC_DIFVHOLVL,Kdiff_t,'p',k,ind1,ind2,wghts)

          ! vertical salt diffusivity (vcoord /= 'isopyc_bulkml') [m^2/s]
          call accilv(ACC_DIFVSOLVL,Kdiff_s,'p',k,ind1,ind2,wghts)

          ! potential vorticity [s m-2]
          call acclvl(ACC_PVLVL,pv_p,'p',k,ind1,ind2,wghts)

          if (use_TRC .and. use_TKE) then
            ! tke [m^2/s^2]
            call acclvl(ACC_TKELVL,trc(1-nbdy,1-nbdy,k1m,itrtke),'p',k,ind1,ind2,wghts)

            ! gls_psi [m^2/s^3]
            call acclvl(ACC_GLS_PSILVL,trc(1-nbdy,1-nbdy,k1m,itrgls),'p',k,ind1,ind2,wghts)

          end if
          ! layer thickness [m]
          call acclvl(ACC_DZLVL,dz,'p',k,ind1,ind2,wghts)

        end do
      end if
    else
      if (sum(ACC_SALNLVL(1:nphy)   +ACC_TEMPLVL(1:nphy) &
             +ACC_BFSQLVL(1:nphy)   +ACC_DIFDIALVL(1:nphy) &
             +ACC_DIFVMOLVL(1:nphy) +ACC_DIFVHOLVL(1:nphy) &
             +ACC_DIFVSOLVL(1:nphy) +ACC_DIFINTLVL(1:nphy) &
             +ACC_DIFISOLVL(1:nphy) +ACC_TKELVL(1:nphy) &
             +ACC_GLS_PSILVL(1:nphy)+ACC_PVLVL(1:nphy) &
             +ACC_DZLVL(1:nphy)     ) /= 0) then
        do k = 1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)

          ! salinity [g/kg]
          call acclvl(ACC_SALNLVL,saln(1-nbdy,1-nbdy,k1m),'p',k,ind1,ind2,wghts)

          ! temperature [degC]
          call acclvl(ACC_TEMPLVL,temp(1-nbdy,1-nbdy,k1m),'p',k,ind1,ind2,wghts)

          ! buoyancy frequency squared [1/s^2]
          call acclvl(ACC_BFSQLVL,bfsql,'p',k,ind1,ind2,wghts)

          ! layer interface diffusivity [m^2/s]
          call acclvl(ACC_DIFINTLVL,difint,'p',k,ind1,ind2,wghts)

          ! isopycnal diffusivity [m^2/s]
          call acclvl(ACC_DIFISOLVL,difiso,'p',k,ind1,ind2,wghts)

          ! vertical diffusivity (vcoord == 'isopyc_bulkml') [m^2/s]
          call acclvl(ACC_DIFDIALVL,difdia,'p',k,ind1,ind2,wghts)

          ! vertical momentum diffusivity (vcoord /= 'isopyc_bulkml') [m^2/s]
          call accilv(ACC_DIFVMOLVL,Kvisc_m,'p',k,ind1,ind2,wghts)

          ! vertical heat diffusivity (vcoord /= 'isopyc_bulkml') [m^2/s]
          call accilv(ACC_DIFVHOLVL,Kdiff_t,'p',k,ind1,ind2,wghts)

          ! vertical salt diffusivity (vcoord /= 'isopyc_bulkml') [m^2/s]
          call accilv(ACC_DIFVSOLVL,Kdiff_s,'p',k,ind1,ind2,wghts)

          ! potential vorticity [s m-2]
          call acclvl(ACC_PVLVL,pv_p,'p',k,ind1,ind2,wghts)

          if (use_TRC .and. use_TKE) then
            ! tke [m^2/s^2]
            call acclvl(ACC_TKELVL,trc(1-nbdy,1-nbdy,k1m,itrtke),'p',k,ind1,ind2,wghts)

            ! gls_psi [m^2/s^3]
            call acclvl(ACC_GLS_PSILVL,trc(1-nbdy,1-nbdy,k1m,itrgls),'p',k,ind1,ind2,wghts)

          end if
          ! layer thickness [m]
          call acclvl(ACC_DZLVL,dz,'p',k,ind1,ind2,wghts)

        end do
      end if
    end if

    ! Accumulate vertical velocity
    do iogrp = 1,nphy
      if (ACC_WFLX(iogrp)+ACC_WFLX2(iogrp)+ACC_WFLXLVL(iogrp) + acc_wflx2lvl(iogrp) /= 0) then
         call diavfl(iogrp,m,n,mm,nn,k1m,k1n)
      end if
    end do

  end subroutine diaacc



  subroutine diaout_alarms
  !----------------------------------
  ! Set alarms for when output of io group is due.
  !----------------------------------

    ! Local variables
    integer :: iogrp

    do iogrp = 1,nphy
      if (((diagann_phy(iogrp).and.nday_of_year == 1.or.diagmon_phy(iogrp) &
           .and.date%day == 1).and.mod(nstep,nstep_in_day) == 0).or. &
           .not.(diagann_phy(iogrp).or.diagmon_phy(iogrp)).and. &
           mod(nstep+.5,diagfq_phy(iogrp)) < 1.) then
        alarm_phy(iogrp) = 1
      else
        alarm_phy(iogrp) = 0
      end if
    end do

  end subroutine diaout_alarms



  subroutine diaout(m,n,mm,nn,k1m,k1n)
  !----------------------------------
  ! Write diagnostic fields if io group alarms are activated.
  !----------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    integer :: iogrp

    do iogrp = 1, nphy
      if (alarm_phy(iogrp) == 1) call diaout_iogrp(iogrp,m,n,mm,nn,k1m,k1n)
    end do

  end subroutine diaout



  subroutine diaout_iogrp(iogrp,m,n,mm,nn,k1m,k1n)
  !----------------------------------
  ! Write diagnostic fields for io group
  !----------------------------------

    ! Arguments
    integer, intent(in) :: iogrp,m,n,mm,nn,k1m,k1n

    ! Local variables
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: iuu
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: ivv
    integer,               save :: irec(nphymax)
    character(len=fnmlen), save :: fname(nphymax)
    logical                     :: iniflg = .true.
    logical                     :: append2file(nphymax) = .false.
    integer                     :: i,j,k,l,cmpflg
    character(len=30)           :: timeunits
    character(len=20)           :: startdate
    real                        :: datenum,rnacc
    real, dimension(itdm,jtdm)  :: bflxg,strg
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghtsflx
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: tmp3d
    integer :: nt,nat,km
    character :: trcnm*80,trcnml*80,dimslyr*100
    real :: treps
    parameter (treps = 1.e-14)

    if (mnproc == 1) &
         write (lp,'(a,f6.2,a)') 'diaout_iogrp: fields averaged over ', &
         real(nacc_phy(iogrp))/nstep_in_day,' days'

    rnacc = 1./real(nacc_phy(iogrp))
    cmpflg = GLB_COMPFLAG(iogrp)

    ! compute meridional transports and transports through sections
    if (ACC_MMFLXL(iogrp)+ACC_MMFLXD(iogrp)+ACC_MMFTDL(iogrp) &
       +ACC_MMFSML(iogrp)+ACC_MMFTDD(iogrp)+ACC_MMFSMD(iogrp) &
       +ACC_MHFLX (iogrp)+ACC_MHFTD (iogrp)+ACC_MHFSM (iogrp) &
       +ACC_MHFLD (iogrp)+ACC_MSFLX (iogrp)+ACC_MSFTD (iogrp) &
       +acc_msfsm (iogrp) +acc_msfld(iogrp) /= 0) call diamer(iogrp)
    if (acc_voltr(iogrp) /= 0) call diasec(iogrp)

    ! compute barotropic mass streamfunction
    if (h2d_btmstr(iogrp) /= 0) then
      if     (nreg <= 2) then
        call xcaget(bflxg,phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(iogrp)),1)
        if (mnproc == 1) then
          do i = 1,itdm
            strg(i,1) = 0.
          end do
          do j = 1,jtdm-1
            do i = 1,itdm
              strg(i,j+1) = strg(i,j)-bflxg(i,j)
            end do
          end do
        end if
      else if (nreg == 4) then
        call xcaget(bflxg,phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(iogrp)),1)
        if (mnproc == 1) then
          do j = 1,jtdm
            strg(itdm,j) = 0.
          end do
          do j = 1,jtdm
            do i = itdm-1,1,-1
              strg(i,j) = strg(i+1,j)-bflxg(i,j)
            end do
          end do
        end if
      else
        if (mnproc == 1) then
          write (lp,'(a,i2,a)') &
               'mod_dia: cannot compute streamfunction for nreg =',nreg, &
               '!'
        end if
        call xcstop('(mod_dia)')
        stop '(mod_dia)'
      end if
      call xcaput(strg,util1,1)
      call xctilr(util1,1,1, 1,1, halo_qs)
      call inih2d(ACC_UTILH2D(1),'p',0.)
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            phyh2d(i,j,ACC_UTILH2D(1))= &
                 .25*(util1(i  ,j  )+util1(i+1,j  ) &
                     +util1(i  ,j+1)+util1(i+1,j+1))
          end do
        end do
      end do
      !$omp end parallel do
    end if

    ! compute global sums and averages
    if (MSC_MASSGS(iogrp)+MSC_SALNGA(iogrp) +msc_tempga(iogrp) /= 0) then

      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = phylyr(i,j,1,ACC_DP(iogrp))
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j) = util1(i,j)+phylyr(i,j,k,ACC_DP(iogrp))
            end do
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = util1(i,j)*scp2(i,j)
          end do
        end do
      end do
      !$omp end parallel do
      call xcsum(massgs(1),util1,ips)
    end if

    if (msc_volgs(iogrp) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = phylyr(i,j,1,ACC_DZ(iogrp))
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j) = util1(i,j)+phylyr(i,j,k,ACC_DZ(iogrp))
            end do
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = util1(i,j)*scp2(i,j)
          end do
        end do
      end do
      !$omp end parallel do
      call xcsum(volgs(1),util1,ips)
      volgs(1) = rnacc*volgs(1)/grav
    end if
    if (msc_salnga(iogrp) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = phylyr(i,j,1,ACC_SALN(iogrp))
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j) = util1(i,j)+phylyr(i,j,k,ACC_SALN(iogrp))
            end do
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = util1(i,j)*scp2(i,j)
          end do
        end do
      end do
      !$omp end parallel do
      call xcsum(salnga(1),util1,ips)
      salnga(1) = salnga(1)/massgs(1)
    end if
    if (msc_tempga(iogrp) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = phylyr(i,j,1,ACC_TEMP(iogrp))
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j) = util1(i,j)+phylyr(i,j,k,ACC_TEMP(iogrp))
            end do
          end do
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = util1(i,j)*scp2(i,j)
          end do
        end do
      end do
      !$omp end parallel do
      call xcsum(tempga(1),util1,ips)
      tempga(1) = tempga(1)/massgs(1)
    end if
    if (msc_massgs(iogrp) /= 0) then
      massgs(1) = rnacc*massgs(1)/grav
    end if
    if (msc_sssga(iogrp) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = phyh2d(i,j,ACC_SSS(iogrp))*scp2(i,j)
          end do
        end do
      end do
      !$omp end parallel do
      call xcsum(sssga(1),util1,ips)
      sssga(1) = rnacc*sssga(1)/area
    end if
    if (msc_sstga(iogrp) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = phyh2d(i,j,ACC_SST(iogrp))*scp2(i,j)
          end do
        end do
      end do
      !$omp end parallel do
      call xcsum(sstga(1),util1,ips)
      sstga(1) = rnacc*sstga(1)/area
    end if

    ! finalize accumulation of 2d fields
    call finh2d(ACC_HICE(iogrp),ACC_FICE(iogrp),'p')
    call finh2d(ACC_HSNW(iogrp),ACC_FICE(iogrp),'p')
    call finh2d(ACC_UICE(iogrp),ACC_IVOLU(iogrp),'u')
    call finh2d(ACC_VICE(iogrp),ACC_IVOLV(iogrp),'v')

    ! finalize accumulation of layer fields
    call finlyr(ACC_UVEL(iogrp),ACC_DPU(iogrp),'u')
    call finlyr(ACC_VVEL(iogrp),ACC_DPV(iogrp),'v')
    call finlyr(ACC_SALN(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_TEMP(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_BFSQ(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_DIFDIA(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_DIFVMO(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_DIFVHO(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_DIFVSO(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_DIFINT(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_DIFISO(iogrp),ACC_DP(iogrp),'p')
    call finlyr(ACC_AVDSG(iogrp),ACC_DPVOR(iogrp),'p')
    if (use_TRC .and. use_TKE) then
      call finlyr(ACC_TKE(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_GLS_PSI(iogrp),ACC_DP(iogrp),'p')
    end if

    ! compute log10 of diffusivities
    if (lyr_difdia(iogrp) == 2) &
         call loglyr(ACC_DIFDIA(iogrp),'p',1.,0.)
    if (lyr_difvmo(iogrp) == 2) &
         call loglyr(ACC_DIFVMO(iogrp),'p',1.,0.)
    if (lyr_difvho(iogrp) == 2) &
         call loglyr(ACC_DIFVHO(iogrp),'p',1.,0.)
    if (lyr_difvso(iogrp) == 2) &
         call loglyr(ACC_DIFVSO(iogrp),'p',1.,0.)
    if (lyr_difint(iogrp) == 2) &
         call loglyr(ACC_DIFINT(iogrp),'p',1.,0.)
    if (lyr_difiso(iogrp) == 2) &
         call loglyr(ACC_DIFISO(iogrp),'p',1.,0.)

    if (lvl_difdia(iogrp) == 2) &
         call loglvl(ACC_DIFDIALVL(iogrp),'p',rnacc,0.)
    if (lvl_difvmo(iogrp) == 2) &
         call loglvl(ACC_DIFVMOLVL(iogrp),'p',rnacc,0.)
    if (lvl_difvho(iogrp) == 2) &
         call loglvl(ACC_DIFVHOLVL(iogrp),'p',rnacc,0.)
    if (lvl_difvso(iogrp) == 2) &
         call loglvl(ACC_DIFVSOLVL(iogrp),'p',rnacc,0.)
    if (lvl_difint(iogrp) == 2) &
         call loglvl(ACC_DIFINTLVL(iogrp),'p',rnacc,0.)
    if (lvl_difiso(iogrp) == 2) &
         call loglvl(ACC_DIFISOLVL(iogrp),'p',rnacc,0.)

    ! mask sea floor of level fields
    call msklvl(ACC_BFSQLVL(iogrp),'p')
    call msklvl(ACC_DIFDIALVL(iogrp),'p')
    call msklvl(ACC_DIFVMOLVL(iogrp),'p')
    call msklvl(ACC_DIFVHOLVL(iogrp),'p')
    call msklvl(ACC_DIFVSOLVL(iogrp),'p')
    call msklvl(ACC_DIFINTLVL(iogrp),'p')
    call msklvl(ACC_DIFISOLVL(iogrp),'p')
    call msklvl(ACC_DZLVL(iogrp),'p')
    call msklvl(ACC_UVELLVL(iogrp),'u')
    call msklvl(ACC_VVELLVL(iogrp),'v')
    call msklvl(ACC_UFLXLVL(iogrp),'u')
    call msklvl(ACC_VFLXLVL(iogrp),'v')
    call msklvl(ACC_UTFLXLVL(iogrp),'u')
    call msklvl(ACC_VTFLXLVL(iogrp),'v')
    call msklvl(ACC_USFLXLVL(iogrp),'u')
    call msklvl(ACC_VSFLXLVL(iogrp),'v')
    call msklvl(ACC_UMFLTDLVL(iogrp),'u')
    call msklvl(ACC_VMFLTDLVL(iogrp),'v')
    call msklvl(ACC_UMFLSMLVL(iogrp),'u')
    call msklvl(ACC_VMFLSMLVL(iogrp),'v')
    call msklvl(ACC_UTFLTDLVL(iogrp),'u')
    call msklvl(ACC_VTFLTDLVL(iogrp),'v')
    call msklvl(ACC_UTFLSMLVL(iogrp),'u')
    call msklvl(ACC_VTFLSMLVL(iogrp),'v')
    call msklvl(ACC_UTFLLDLVL(iogrp),'u')
    call msklvl(ACC_VTFLLDLVL(iogrp),'v')
    call msklvl(ACC_USFLTDLVL(iogrp),'u')
    call msklvl(ACC_VSFLTDLVL(iogrp),'v')
    call msklvl(ACC_USFLSMLVL(iogrp),'u')
    call msklvl(ACC_VSFLSMLVL(iogrp),'v')
    call msklvl(ACC_USFLLDLVL(iogrp),'u')
    call msklvl(ACC_VSFLLDLVL(iogrp),'v')
    call msklvl(ACC_SALNLVL(iogrp),'p')
    call msklvl(ACC_TEMPLVL(iogrp),'p')
    call msklvl(ACC_WFLXLVL(iogrp),'p')
    call msklvl(ACC_WFLX2LVL(iogrp),'p')
    call msklvl(ACC_PVLVL(iogrp),'p')
    if (use_TRC .and. use_TKE) then
      call msklvl(ACC_TKELVL(iogrp),'p')
      call msklvl(ACC_GLS_PSILVL(iogrp),'p')
    end if

    ! get instantaneous values for ice age
    if (acc_iage(iogrp) /= 0) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            phyh2d(i,j,ACC_IAGE(iogrp)) = iagem(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    ! set time information
    timeunits = ' '
    startdate = ' '
    write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)') &
         'days since ',min(1800,date0%year),'-',1,'-',1,' 00:00'
    write(startdate,'(i4.4,a1,i2.2,a1,i2.2,a6)') &
         date0%year,'-',date0%month,'-',date0%day,' 00:00'
    datenum = time-time0-.5*diagfq_phy(iogrp)/nstep_in_day

    ! create file name
    if (.not.append2file(iogrp)) then
      call diafnm(GLB_FNAMETAG(IOGRP), &
           filefq_phy(iogrp)/real(nstep_in_day), &
           filemon_phy(iogrp),fileann_phy(iogrp),fname(iogrp))
      append2file(iogrp) = .true.
      irec(iogrp) = 1
    else
      irec(iogrp) = irec(iogrp)+1
    end if
    if (((fileann_phy(iogrp).and.nday_of_year == 1.or. &
         filemon_phy(iogrp).and.date%day == 1).and. &
         mod(nstep,nstep_in_day) == 0).or. &
         .not.(fileann_phy(iogrp).or.filemon_phy(iogrp)).and. &
         mod(nstep+.5,filefq_phy(iogrp)) < 1.) then
      append2file(iogrp) = .false.
    end if

    ! open output file
    if (mnproc == 1) &
         write (lp,'(2a)') 'Writing physical diagnostics to file: ', &
         trim(fname(iogrp))
    if (glb_ncformat(iogrp) == 1) then
      call ncfopn(fname(iogrp),'w','6',irec(iogrp),iotype)
    else if (glb_ncformat(iogrp) == 2) then
      call ncfopn(fname(iogrp),'w','h',irec(iogrp),iotype)
    else
      call ncfopn(fname(iogrp),'w','c',irec(iogrp),iotype)
    end if

    ! compute extended ocean masks
    if (iniflg) then
      iniflg = .false.
      !$omp parallel do private(i)
      do j = 1,jj
        do i = 1,ii
          if((ip(i,j)+ip(i-1,j)) >= 1) then
            iuu(i,j) = 1
          else
            iuu(i,j) = 0
          end if
          if((ip(i,j)+ip(i,j-1)) >= 1) then
            ivv(i,j) = 1
          else
            ivv(i,j) = 0
          end if
        end do
      end do
      !$omp end parallel do
    end if

    ! define output dimensions
    if (cmpflg /= 0) then
      call ncdimc('pcomp',ip,0)
      call ncdimc('ucomp',iuu,0)
      call ncdimc('vcomp',ivv,0)
    else
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
    end if
    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      call ncdims('layer',kdm)
    else
      call ncdims('sigma',kdm)
    endif
    call ncdims('depth',ddm)
    call ncdims('bounds',2)
    call ncdims('time',0)

    if (ACC_MMFLXL(iogrp)+ACC_MMFLXD(iogrp)+ACC_MMFTDL(iogrp) &
       +ACC_MMFSML(iogrp)+ACC_MMFTDD(iogrp)+ACC_MMFSMD(iogrp) &
       +ACC_MHFLX (iogrp)+ACC_MHFTD (iogrp)+ACC_MHFSM (iogrp) &
       +ACC_MHFLD (iogrp)+ACC_MSFLX (iogrp)+ACC_MSFTD (iogrp) &
       +acc_msfsm (iogrp)+acc_msfld (iogrp)+msc_voltr (iogrp) /= 0) then
      call ncdims('slenmax',slenmax)
    end if

    if (ACC_MMFLXL(iogrp)+ACC_MMFLXD(iogrp)+ACC_MMFTDL(iogrp) &
       +ACC_MMFSML(iogrp)+ACC_MMFTDD(iogrp)+ACC_MMFSMD(iogrp) &
       +ACC_MHFLX (iogrp)+ACC_MHFTD (iogrp)+ACC_MHFSM (iogrp) &
       +ACC_MHFLD (iogrp)+ACC_MSFLX (iogrp)+ACC_MSFTD (iogrp) &
       +acc_msfsm (iogrp)+acc_msfld (iogrp) /= 0) then
      if ((lmax > 0.and.lmax <= ldm)) then
        call ncdims('lat',lmax)
        call ncdims('region',mer_nreg)
      else
        write (lp,*) &
             'Illegal dimension of meridional diagnostics: lmax = ',lmax, &
             ' ldm = ',ldm
        call xchalt('(diaout_iogrp)')
        stop '(diaout_iogrp)'
      end if
    end if

    if (acc_voltr(iogrp) /= 0) then
      if ((sec_num > 0.and.sec_num <= max_sec)) then
        call ncdims('section',sec_num)
      else
        write (lp,*) 'Illegal number of sections: sec_num = ',sec_num, &
             ' max_sec = ',max_sec
        call xchalt('(diaout_iogrp)')
        stop '(diaout_iogrp)'
      end if
    end if

    ! Record vertical coordinate settings as global attributes.
    call ncattr('vcoord_type',trim(vcoord_type))
    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      call ncattr('sigref_spec',trim(sigref_spec))
      if (trim(sigref_spec) == 'function') then
        if (sigref_adaption) then
          call ncattr('sigref_adaption','.true.')
        else
          call ncattr('sigref_adaption','.false.')
        endif
        call ncputr('sigref_fun_spec_dsdz_bot', sigref_fun_spec%dsdz_bot)
        call ncputr('sigref_fun_spec_sp1'     , sigref_fun_spec%sp1)
        call ncputr('sigref_fun_spec_zp2'     , sigref_fun_spec%zp2)
        call ncputr('sigref_fun_spec_zp3'     , sigref_fun_spec%zp3)
        call ncputr('sigref_fun_spec_sp4'     , sigref_fun_spec%sp4)
        call ncputr('sigref_fun_spec_z_top'   , sigref_fun_spec%z_top)
        call ncputr('sigref_fun_spec_s_top'   , sigref_fun_spec%s_top)
        call ncputr('sigref_fun_spec_z_bot'   , sigref_fun_spec%z_bot)
        call ncputr('sigref_fun_spec_s_bot'   , sigref_fun_spec%s_bot)
      endif
    endif

    call definevar(irec(iogrp),iogrp,cmpflg,timeunits,calendar)

    call nctime(datenum,calendar,timeunits,startdate)

    ! write auxillary dimension information
    if (irec(iogrp) == 1) then
      if (vcoord_tag /= vcoord_isopyc_bulkml) then
        ! sigma levels
        call ncwrt1('sigma','layer',sigref(1:kk))
        call ncattr('long_name','Potential density')
        call ncattr('standard_name','sea_water_sigma_theta')
        call ncattr('units','kg m-3')
        call ncattr('positive','down')
      else
        ! sigma levels
        call ncwrt1('sigma','sigma',sigref(1:kk))
        call ncattr('long_name','Potential density')
        call ncattr('standard_name','sea_water_sigma_theta')
        call ncattr('units','kg m-3')
        call ncattr('positive','down')
      endif
      ! zlevel
      call ncwrt1('depth','depth',depthslev)
      call ncattr('long_name','z level')
      call ncattr('units','m')
      call ncattr('positive','down')
      call ncattr('bounds','depth_bnds')
      call ncwrt1('depth_bnds','bounds depth',depthslev_bnds)
      if (MSC_MMFLXL(iogrp)+MSC_MMFLXD(iogrp)+MSC_MMFTDL(iogrp) &
         +MSC_MMFSML(iogrp)+MSC_MMFTDD(iogrp)+MSC_MMFSMD(iogrp) &
         +MSC_MHFLX (iogrp)+MSC_MHFTD (iogrp)+MSC_MHFSM (iogrp) &
         +MSC_MHFLD (iogrp)+MSC_MSFLX (iogrp)+MSC_MSFTD (iogrp) &
         +msc_msfsm (iogrp)+msc_msfld (iogrp) /= 0) then
        call ncwrt1('lat','lat',mtlat)
        call ncattr('long_name','Latitude')
        call ncattr('standard_name','latitude')
        call ncattr('units','degree_north')
        call ncwrtc('region','slenmax region',mer_regnam)
        call ncattr('long_name','Region name')
      end if
      if (msc_voltr(iogrp) /= 0) then
        call ncwrtc('section','slenmax section',sec_name)
        call ncattr('long_name','Section name')
      end if
    end if

    ! write 2d fields
    call wrth2d(ACC_SIGMX(iogrp),H2D_SIGMX(iogrp),rnacc,0., &
         cmpflg,ip,'p','sigmx','Mixed layer density',' ','kg m-3')

    call wrth2d(ACC_UB(iogrp),H2D_UB(iogrp),rnacc,0., &
         cmpflg,iuu,'u','ubaro','Barotropic velocity x-component',' ','m s-1')

    call wrth2d(ACC_VB(iogrp),H2D_VB(iogrp),rnacc,0., &
         cmpflg,ivv,'v','vbaro','Barotropic velocity y-component',' ','m s-1')

    call wrth2d(ACC_PSRF(iogrp),H2D_PSRF(iogrp),rnacc,0., &
         cmpflg,ip,'p','psrf','Surface pressure',' ','Pa')

    call wrth2d(ACC_PBOT(iogrp),H2D_PBOT(iogrp),rnacc,0., &
         cmpflg,ip,'p','pbot','Bottom pressure',' ','Pa')

    call wrth2d(ACC_SEALV(iogrp),H2D_SEALV(iogrp),-rnacc,0., &
         cmpflg,ip,'p','sealv','Sea level',' ','m')

    call wrth2d(ACC_SLVSQ(iogrp),H2D_SLVSQ(iogrp),rnacc,0., &
         cmpflg,ip,'p','slvsq','Sea level squared',' ','m2')

    call wrth2d(ACC_UTILH2D(1),H2D_BTMSTR(iogrp), &
         rnacc*.5*dlt/(grav*baclin),0.,cmpflg,ip,'p','btmstr', &
         'Barotropic mass streamfunction',' ','kg s-1')

    call wrth2d(ACC_HICE(iogrp),H2D_HICE(iogrp),1.,0., &
         cmpflg,ip,'p','hice','Ice thickness',' ','m')

    call wrth2d(ACC_TICE(iogrp),H2D_TICE(iogrp),rnacc,-t0deg, &
         cmpflg,ip,'p','tice','Ice temperature',' ','degC')

    call wrth2d(ACC_HSNW(iogrp),H2D_HSNW(iogrp),1.,0., &
         cmpflg,ip,'p','hsnw','Snow depth',' ','m')

    call wrth2d(ACC_FICE(iogrp),H2D_FICE(iogrp),rnacc*1e2,0., &
         cmpflg,ip,'p','fice','Ice concentration',' ','%')

    call wrth2d(ACC_TSRF(iogrp),H2D_TSRF(iogrp),rnacc,-t0deg, &
         cmpflg,ip,'p','tsrf','Surface temperature',' ','degC')

    call wrth2d(ACC_IAGE(iogrp),H2D_IAGE(iogrp),1.,0., &
         cmpflg,ip,'p','iage','Ice age',' ','day')

    call wrth2d(ACC_UICE(iogrp),H2D_UICE(iogrp),1.,0., &
         cmpflg,iuu,'u','uice','Ice velocity x-component',' ','m s-1')

    call wrth2d(ACC_VICE(iogrp),H2D_VICE(iogrp),1.,0., &
         cmpflg,ivv,'v','vice','Ice velocity y-component',' ','m s-1')

    call wrth2d(ACC_SWA(iogrp),H2D_SWA(iogrp),rnacc,0., &
         cmpflg,ip,'p','swa','Short-wave heat flux',' ','W m-2')

    call wrth2d(ACC_NSF(iogrp),H2D_NSF(iogrp),rnacc,0., &
         cmpflg,ip,'p','nsf','Non-solar heat flux',' ','W m-2')

    call wrth2d(ACC_HMAT(iogrp),H2D_HMAT(iogrp),rnacc,0., &
         cmpflg,ip,'p','hmat','Heat flux due to material enthalpy flux', &
         ' ','W m-2')

    call wrth2d(ACC_HMLTFZ(iogrp),H2D_HMLTFZ(iogrp),rnacc,0., &
         cmpflg,ip,'p','hmltfz','Heat flux due to melting/freezing',' ','W m-2')

    call wrth2d(ACC_DFL(iogrp),H2D_DFL(iogrp),rnacc,0., &
         cmpflg,ip,'p','dfl','Non-solar heat flux derivative',' ','W m-2 K-1')

    call wrth2d(ACC_SURFLX(iogrp),H2D_SURFLX(iogrp),-rnacc,0., &
         cmpflg,ip,'p','hflx','Heat flux received by ocean',' ','W m-2')

    call wrth2d(ACC_SURRLX(iogrp),H2D_SURRLX(iogrp),-rnacc,0., &
         cmpflg,ip,'p','hrflx', &
         'Restoring heat flux received by ocean',' ','W m-2')

    call wrth2d(ACC_LIP(iogrp),H2D_LIP(iogrp),rnacc,0., &
         cmpflg,ip,'p','lip','Liquid precipitation',' ','kg m-2 s-1')

    call wrth2d(ACC_SOP(iogrp),H2D_SOP(iogrp),rnacc,0., &
         cmpflg,ip,'p','sop','Solid precipitation',' ','kg m-2 s-1')

    call wrth2d(ACC_EVA(iogrp),H2D_EVA(iogrp),rnacc,0., &
         cmpflg,ip,'p','eva','Evaporation',' ','kg m-2 s-1')

    call wrth2d(ACC_FMLTFZ(iogrp),H2D_FMLTFZ(iogrp),rnacc,0., &
         cmpflg,ip,'p','fmltfz', &
         'Fresh water flux due to melting/freezing',' ','kg m-2 s-1')

    call wrth2d(ACC_RNFFLX(iogrp),H2D_RNFFLX(iogrp),rnacc,0., &
         cmpflg,ip,'p','rnf','Liquid runoff',' ','kg m-2 s-1')

    call wrth2d(ACC_RFIFLX(iogrp),H2D_RFIFLX(iogrp),rnacc,0., &
         cmpflg,ip,'p','rfi','Frozen runoff',' ','kg m-2 s-1')

    call wrth2d(ACC_SALFLX(iogrp),H2D_SALFLX(iogrp),-rnacc*g2kg,0., &
         cmpflg,ip,'p','sflx','Salt flux received by ocean',' ','kg m-2 s-1')

    call wrth2d(ACC_SALRLX(iogrp),H2D_SALRLX(iogrp),-rnacc*g2kg,0., &
         cmpflg,ip,'p','srflx', &
         'Restoring salt flux received by ocean',' ','kg m-2 s-1')

    call wrth2d(ACC_BRNFLX(iogrp),H2D_BRNFLX(iogrp),-rnacc*g2kg,0., &
         cmpflg,ip,'p','bflx','Brine flux',' ','kg m-2 s-1')

    call wrth2d(ACC_ZTX(iogrp),H2D_ZTX(iogrp),rnacc,0., &
         cmpflg,iuu,'u','ztx','Wind stress x-component',' ','N m-2')

    call wrth2d(ACC_MTY(iogrp),H2D_MTY(iogrp),rnacc,0., &
         cmpflg,ivv,'v','mty','Wind stress y-component',' ','N m-2')

    call wrth2d(ACC_TAUX(iogrp),H2D_TAUX(iogrp),rnacc,0., &
         cmpflg,iuu,'u','taux', &
         'Momentum flux received by ocean x-component',' ','N m-2')

    call wrth2d(ACC_TAUY(iogrp),H2D_TAUY(iogrp),rnacc,0., &
         cmpflg,ivv,'v','tauy', &
         'Momentum flux received by ocean y-component',' ','N m-2')

    call wrth2d(ACC_IDKEDT(iogrp),H2D_IDKEDT(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','idkedt', &
         'Mixed layer inertial kinetic energy tendency per unit area', &
         ' ','kg s-3')

    call wrth2d(ACC_USTAR(iogrp),H2D_USTAR(iogrp),rnacc,0., &
         cmpflg,ip,'p','ustar','Friction velocity',' ','m s-1')

    call wrth2d(ACC_USTAR3(iogrp),H2D_USTAR3(iogrp),rnacc,0., &
         cmpflg,ip,'p','ustar3','Friction velocity cubed',' ','m3 s-3')

    call wrth2d(ACC_ABSWND(iogrp),H2D_ABSWND(iogrp),rnacc,0., &
         cmpflg,ip,'p','abswnd','Absolute wind speed',' ','m s-1')

    call wrth2d(ACC_MTKEUS(iogrp),H2D_MTKEUS(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','mtkeus', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to friction velocity', &
         ' ','kg s-3')

    call wrth2d(ACC_MTKENI(iogrp),H2D_MTKENI(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','mtkeni', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to near inertial motions', &
         ' ','kg s-3')

    call wrth2d(ACC_MTKEBF(iogrp),H2D_MTKEBF(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','mtkebf', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to buoyancy forcing', &
         ' ','kg s-3')

    call wrth2d(ACC_MTKERS(iogrp),H2D_MTKERS(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','mtkers', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to eddy restratification', &
         ' ','kg s-3')

    call wrth2d(ACC_MTKEPE(iogrp),H2D_MTKEPE(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','mtkepe', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to potential energy change', &
         ' ','kg s-3')

    call wrth2d(ACC_MTKEKE(iogrp),H2D_MTKEKE(iogrp),rnacc/alpha0,0., &
         cmpflg,ip,'p','mtkeke', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to kinetic energy change', &
         ' ','kg s-3')

    call wrth2d(ACC_LAMULT(iogrp),H2D_LAMULT(iogrp),rnacc,0., &
         cmpflg,ip,'p','lamult','Langmuir enhancement factor',' ','1')

    call wrth2d(ACC_LASL(iogrp),H2D_LASL(iogrp),rnacc,0., &
         cmpflg,ip,'p','lasl','Surface layer averaged Langmuir number',' ','1')

    call wrth2d(ACC_USTOKES(iogrp),H2D_USTOKES(iogrp),rnacc,0., &
         cmpflg,iuu,'u','ustokes', &
         'Surface Stokes drift x-component',' ','m s-1')

    call wrth2d(ACC_VSTOKES(iogrp),H2D_VSTOKES(iogrp),rnacc,0., &
         cmpflg,ivv,'v','vstokes', &
         'Surface Stokes drift y-component',' ','m s-1')

    call wrth2d(ACC_SFL(iogrp),H2D_SFL(iogrp),rnacc,0., &
         cmpflg,ip,'p','sfl','Salt flux',' ','kg m-2 s-1')

    call wrth2d(ACC_ALB(iogrp),H2D_ALB(iogrp),rnacc,0., &
         cmpflg,ip,'p','alb','Surface albedo',' ','1')

    call wrth2d(ACC_MLD(iogrp),H2D_MLD(iogrp),rnacc/onem,0., &
         cmpflg,ip,'p','mld','Mixed layer depth',' ','m')

    call wrth2d(ACC_MAXMLD(iogrp),H2D_MAXMLD(iogrp),1./onem,0., &
         cmpflg,ip,'p','maxmld','Maximum mixed layer depth',' ','m')

    call wrth2d(ACC_MLTS(iogrp),H2D_MLTS(iogrp),rnacc,0., &
         cmpflg,ip,'p','mlts', &
         'Mixed layer thickness defined by sigma t',' ','m')

    call wrth2d(ACC_MLTSMN(iogrp),H2D_MLTSMN(iogrp),1.,0., &
         cmpflg,ip,'p','mltsmn', &
         'Minimum mixed layer thickness defined by sigma t',' ','m')

    call wrth2d(ACC_MLTSMX(iogrp),H2D_MLTSMX(iogrp),1.,0., &
         cmpflg,ip,'p','mltsmx', &
         'Maximum mixed layer thickness defined by sigma t',' ','m')

    call wrth2d(ACC_MLTSSQ(iogrp),H2D_MLTSSQ(iogrp),rnacc,0., &
         cmpflg,ip,'p','mltssq', &
         'Mixed layer thickness squared defined by sigma t',' ','m2')

    call wrth2d(ACC_T20D(iogrp),H2D_T20D(iogrp),rnacc,0., &
         cmpflg,ip,'p','t20d','20C isoterm depth',' ','m')

    call wrth2d(ACC_BRNPD(iogrp),H2D_BRNPD(iogrp),rnacc/onem,0., &
         cmpflg,ip,'p','brnpd','Brine plume depth',' ','m')

    call wrth2d(ACC_SSS(iogrp),H2D_SSS(iogrp),rnacc,0., &
         cmpflg,ip,'p','sss','Ocean surface salinity',' ','g kg-1')

    call wrth2d(ACC_SSSSQ(iogrp),H2D_SSSSQ(iogrp),rnacc,0., &
         cmpflg,ip,'p','ssssq','Ocean surface salinity squared',' ','g2 kg-2')

    call wrth2d(ACC_SBOT(iogrp),H2D_SBOT(iogrp),rnacc,0., &
         cmpflg,ip,'p','sbot','Bottom salinity',' ','g kg-1')

    call wrth2d(ACC_SST(iogrp),H2D_SST(iogrp),rnacc,0., &
         cmpflg,ip,'p','sst','Ocean surface temperature',' ','degC')

    call wrth2d(ACC_SSTSQ(iogrp),H2D_SSTSQ(iogrp),rnacc,0., &
         cmpflg,ip,'p','sstsq','Ocean surface temperature squared',' ','degC2')

    call wrth2d(ACC_TBOT(iogrp),H2D_TBOT(iogrp),rnacc,0., &
         cmpflg,ip,'p','tbot','Bottom temperature',' ','degC')

    ! write 3d layer fields
    call wrtlyr(ACC_DP(iogrp),LYR_DP(iogrp),rnacc,0., &
         cmpflg,ip,'p','dp','Layer pressure thickness',' ','Pa')

    call wrtlyr(ACC_DZ(iogrp),LYR_DZ(iogrp),rnacc,0., &
         cmpflg,ip,'p','dz','Layer thickness',' ','m')

    call wrtlyr(ACC_TEMP(iogrp),LYR_TEMP(iogrp),1.,0., &
         cmpflg,ip,'p','temp','Temperature','Ocean temperature','degC')

    call wrtlyr(ACC_SALN(iogrp),LYR_SALN(iogrp),1.,0., &
         cmpflg,ip,'p','saln','Salinity','Ocean salinity','g kg-1')

    call wrtlyr(ACC_UVEL(iogrp),LYR_UVEL(iogrp),1.,0., &
         cmpflg,iuu,'u','uvel','Velocity x-component',' ','m s-1')

    call wrtlyr(ACC_VVEL(iogrp),LYR_VVEL(iogrp),1.,0., &
         cmpflg,ivv,'v','vvel','Velocity y-component',' ','m s-1')

    call wrtlyr(ACC_UFLX(iogrp),LYR_UFLX(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,iuu,'u','uflx', &
         'Mass flux in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VFLX(iogrp),LYR_VFLX(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ivv,'v','vflx', &
         'Mass flux in y-direction',' ','kg s-1')

    call wrtlyr(ACC_UTFLX(iogrp),LYR_UTFLX(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhflx', &
         'Heat flux in x-direction',' ','W')

    call wrtlyr(ACC_VTFLX(iogrp),LYR_VTFLX(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhflx', &
         'Heat flux in y-direction',' ','W')

    call wrtlyr(ACC_USFLX(iogrp),LYR_USFLX(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usflx', &
         'Salt flux in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VSFLX(iogrp),LYR_VSFLX(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsflx', &
         'Salt flux in y-direction',' ','kg s-1')

    call wrtlyr(ACC_UMFLTD(iogrp),LYR_UMFLTD(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,iuu,'u','umfltd', &
         'Mass flux due to thickness diffusion in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VMFLTD(iogrp),LYR_VMFLTD(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ivv,'v','vmfltd', &
         'Mass flux due to thickness diffusion in y-direction',' ','kg s-1')

    call wrtlyr(ACC_UMFLSM(iogrp),LYR_UMFLSM(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,iuu,'u','umflsm', &
         'Mass flux due to submesoscale transport in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VMFLSM(iogrp),LYR_VMFLSM(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ivv,'v','vmflsm', &
         'Mass flux due to submesoscale transport in y-direction',' ','kg s-1')

    call wrtlyr(ACC_UTFLTD(iogrp),LYR_UTFLTD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhfltd', &
         'Heat flux due to thickness diffusion in x-direction',' ','W')

    call wrtlyr(ACC_VTFLTD(iogrp),LYR_VTFLTD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhfltd', &
         'Heat flux due to thickness diffusion in y-direction',' ','W')

    call wrtlyr(ACC_UTFLSM(iogrp),LYR_UTFLSM(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhflsm', &
         'Heat flux due to submesoscale transport in x-direction',' ','W')

    call wrtlyr(ACC_VTFLSM(iogrp),LYR_VTFLSM(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhflsm', &
         'Heat flux due to submesoscale transport in y-direction',' ','W')

    call wrtlyr(ACC_UTFLLD(iogrp),LYR_UTFLLD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhflld', &
         'Heat flux due to lateral diffusion in x-direction',' ','W')

    call wrtlyr(ACC_VTFLLD(iogrp),LYR_VTFLLD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhflld', &
         'Heat flux due to lateral diffusion in y-direction',' ','W')

    call wrtlyr(ACC_USFLTD(iogrp),LYR_USFLTD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usfltd', &
         'Salt flux due to thickness diffusion in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VSFLTD(iogrp),LYR_VSFLTD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsfltd', &
         'Salt flux due to thickness diffusion in y-direction',' ','kg s-1')

    call wrtlyr(ACC_USFLSM(iogrp),LYR_USFLSM(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usflsm', &
         'Salt flux due to submesoscale transport in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VSFLSM(iogrp),LYR_VSFLSM(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsflsm', &
         'Salt flux due to submesoscale transport in y-direction',' ','kg s-1')

    call wrtlyr(ACC_USFLLD(iogrp),LYR_USFLLD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usflld', &
         'Salt flux due to lateral diffusion in x-direction',' ','kg s-1')

    call wrtlyr(ACC_VSFLLD(iogrp),LYR_VSFLLD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsflld', &
         'Salt flux due to lateral diffusion in y-direction',' ','kg s-1')

    call wrtlyr(ACC_WFLX(iogrp),LYR_WFLX(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ip,'p','wflx', &
         'Vertical mass flux',' ','kg s-1')

    call wrtlyr(ACC_WFLX2(iogrp),LYR_WFLX2(iogrp), &
         rnacc*(.5/(grav*baclin))**2,0.,cmpflg,ip,'p','wflx2', &
         'Vertical mass flux squared',' ','kg2 s-2')

    call wrtlyr(ACC_BFSQ(iogrp),LYR_BFSQ(iogrp),1.,0., &
         cmpflg,ip,'p','bfsq','Squared buoyancy frequency',' ','s-1')

    call wrtlyr(ACC_AVDSG(iogrp),LYR_PV(iogrp),grav,0., &
         cmpflg,ip,'p','pv','Potential vorticity',' ','m-1 s-1')

    if (lyr_difint(iogrp) == 2) then
      call wrtlyr(ACC_DIFINT(iogrp),LYR_DIFINT(iogrp),1.,0., &
           cmpflg,ip,'p','difint', &
           'Layer interface diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlyr(ACC_DIFINT(iogrp),LYR_DIFINT(iogrp),1.,0., &
           cmpflg,ip,'p','difint','Layer interface diffusivity',' ','m2 s-1')
    end if

    if (lyr_difiso(iogrp) == 2) then
      call wrtlyr(ACC_DIFISO(iogrp),LYR_DIFISO(iogrp),1.,0., &
           cmpflg,ip,'p','difiso', &
           'Isopycnal diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlyr(ACC_DIFISO(iogrp),LYR_DIFISO(iogrp),1.,0., &
           cmpflg,ip,'p','difiso','Isopycnal diffusivity',' ','m2 s-1')
    end if

    if (lyr_difdia(iogrp) == 2) then
      call wrtlyr(ACC_DIFDIA(iogrp),LYR_DIFDIA(iogrp),1.,0., &
           cmpflg,ip,'p','difdia', &
           'Vertical diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlyr(ACC_DIFDIA(iogrp),LYR_DIFDIA(iogrp),1.,0., &
           cmpflg,ip,'p','difdia','Vertical diffusivity',' ','m2 s-1')
    end if

    if (lyr_difvmo(iogrp) == 2) then
      call wrtlyr(ACC_DIFVMO(iogrp),LYR_DIFVMO(iogrp),1.,0., &
           cmpflg,ip,'p','difvmo', &
           'Vertical momentum diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlyr(ACC_DIFVMO(iogrp),LYR_DIFVMO(iogrp),1.,0., &
           cmpflg,ip,'p','difvmo','Vertical momentum diffusivity',' ','m2 s-1')
    end if

    if (lyr_difvho(iogrp) == 2) then
      call wrtlyr(ACC_DIFVHO(iogrp),LYR_DIFVHO(iogrp),1.,0., &
           cmpflg,ip,'p','difvho', &
           'Vertical heat diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlyr(ACC_DIFVHO(iogrp),LYR_DIFVHO(iogrp),1.,0., &
           cmpflg,ip,'p','difvho','Vertical heat diffusivity',' ','m2 s-1')
    end if

    if (lyr_difvso(iogrp) == 2) then
      call wrtlyr(ACC_DIFVSO(iogrp),LYR_DIFVSO(iogrp),1.,0., &
           cmpflg,ip,'p','difvso', &
           'Vertical salt diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlyr(ACC_DIFVSO(iogrp),LYR_DIFVSO(iogrp),1.,0., &
           cmpflg,ip,'p','difvso','Vertical salt diffusivity',' ','m2 s-1')
    end if

    if (use_TRC .and. use_TKE) then
      call wrtlyr(ACC_TKE(iogrp),LYR_TKE(iogrp),1.,0., &
           cmpflg,ip,'p','tke','TKE','Turbulent kinetic energy','m2 s-2')

      call wrtlyr(ACC_GLS_PSI(iogrp),LYR_GLS_PSI(iogrp),1.,0., &
           cmpflg,ip,'p','gls_psi','GLS_PSI','Generic length scale','m2 s-3')

    end if
    ! Write 3d depth fields
    call wrtlvl(ACC_DZLVL(iogrp),LVL_DZ(iogrp),rnacc,0., &
         cmpflg,ip,'p','dzlvl','Layer thickness',' ','m')

    call wrtlvl(ACC_TEMPLVL(iogrp),LVL_TEMP(iogrp),rnacc,0., &
         cmpflg,ip,'p','templvl','Temperature','Ocean temperature','degC')

    call wrtlvl(ACC_SALNLVL(iogrp),LVL_SALN(iogrp),rnacc,0., &
         cmpflg,ip,'p','salnlvl','Salinity','Ocean salinity','g kg-1')

    call wrtlvl(ACC_UVELLVL(iogrp),LVL_UVEL(iogrp),rnacc,0., &
         cmpflg,iuu,'u','uvellvl','Velocity x-component',' ','m s-1')

    call wrtlvl(ACC_VVELLVL(iogrp),LVL_VVEL(iogrp),rnacc,0., &
         cmpflg,ivv,'v','vvellvl','Velocity y-component',' ','m s-1')

    call wrtlvl(ACC_UFLXLVL(iogrp),LVL_UFLX(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,iuu,'u','uflxlvl', &
         'Mass flux in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VFLXLVL(iogrp),LVL_VFLX(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ivv,'v','vflxlvl', &
         'Mass flux in y-direction',' ','kg s-1')

    call wrtlvl(ACC_UTFLXLVL(iogrp),LVL_UTFLX(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhflxlvl', &
         'Heat flux in x-direction',' ','W')

    call wrtlvl(ACC_VTFLXLVL(iogrp),LVL_VTFLX(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhflxlvl', &
         'Heat flux in y-direction',' ','W')

    call wrtlvl(ACC_USFLXLVL(iogrp),LVL_USFLX(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u', &
         'usflxlvl','Salt flux in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VSFLXLVL(iogrp),LVL_VSFLX(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v', &
         'vsflxlvl','Salt flux in y-direction',' ','kg s-1')

    call wrtlvl(ACC_UMFLTDLVL(iogrp),LVL_UMFLTD(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,iuu,'u','umfltdlvl', &
         'Mass flux due to thickness diffusion in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VMFLTDLVL(iogrp),LVL_VMFLTD(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ivv,'v','vmfltdlvl', &
         'Mass flux due to thickness diffusion in y-direction',' ','kg s-1')

    call wrtlvl(ACC_UMFLSMLVL(iogrp),LVL_UMFLSM(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,iuu,'u','umflsmlvl', &
         'Mass flux due to submesoscale transport in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VMFLSMLVL(iogrp),LVL_VMFLSM(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ivv,'v','vmflsmlvl', &
         'Mass flux due to submesoscale transport in y-direction',' ','kg s-1')

    call wrtlvl(ACC_UTFLTDLVL(iogrp),LVL_UTFLTD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhfltdlvl', &
         'Heat flux due to thickness diffusion in x-direction',' ','W')

    call wrtlvl(ACC_VTFLTDLVL(iogrp),LVL_VTFLTD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhfltdlvl', &
         'Heat flux due to thickness diffusion in y-direction',' ','W')

    call wrtlvl(ACC_UTFLSMLVL(iogrp),LVL_UTFLSM(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhflsmlvl', &
         'Heat flux due to submesoscale transport in x-direction',' ','W')

    call wrtlvl(ACC_VTFLSMLVL(iogrp),LVL_VTFLSM(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhflsmlvl', &
         'Heat flux due to submesoscale transport in y-direction',' ','W')

    call wrtlvl(ACC_UTFLLDLVL(iogrp),LVL_UTFLLD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,iuu,'u','uhflldlvl', &
         'Heat flux due to lateral diffusion in x-direction',' ','W')

    call wrtlvl(ACC_VTFLLDLVL(iogrp),LVL_VTFLLD(iogrp), &
         rnacc*spcifh*.5/(grav*baclin),0.,cmpflg,ivv,'v','vhflldlvl', &
         'Heat flux due to lateral diffusion in y-direction',' ','W')

    call wrtlvl(ACC_USFLTDLVL(iogrp),LVL_USFLTD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usfltdlvl', &
         'Salt flux due to thickness diffusion in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VSFLTDLVL(iogrp),LVL_VSFLTD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsfltdlvl', &
         'Salt flux due to thickness diffusion in y-direction',' ','kg s-1')

    call wrtlvl(ACC_USFLSMLVL(iogrp),LVL_USFLSM(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usflsmlvl', &
         'Salt flux due to submesoscale transport in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VSFLSMLVL(iogrp),LVL_VSFLSM(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsflsmlvl', &
         'Salt flux due to submesoscale transport in y-direction',' ','kg s-1')

    call wrtlvl(ACC_USFLLDLVL(iogrp),LVL_USFLLD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,iuu,'u','usflldlvl', &
         'Salt flux due to lateral diffusion in x-direction',' ','kg s-1')

    call wrtlvl(ACC_VSFLLDLVL(iogrp),LVL_VSFLLD(iogrp), &
         rnacc*.5*g2kg/(grav*baclin),0.,cmpflg,ivv,'v','vsflldlvl', &
         'Salt flux due to lateral diffusion in y-direction',' ','kg s-1')

    call wrtlvl(ACC_WFLXLVL(iogrp),LVL_WFLX(iogrp), &
         rnacc*.5/(grav*baclin),0.,cmpflg,ip,'p','wflxlvl', &
         'Vertical mass flux',' ','kg s-1')

    call wrtlvl(ACC_WFLX2LVL(iogrp),LVL_WFLX2(iogrp), &
         rnacc*(.5/(grav*baclin))**2,0.,cmpflg,ip,'p','wflx2lvl', &
         'Vertical mass flux squared',' ','kg2 s-2')

    call wrtlvl(ACC_BFSQLVL(iogrp),LVL_BFSQ(iogrp),rnacc,0., &
         cmpflg,ip,'p','bfsqlvl','Squared buoyancy frequency',' ','s-1')

    call wrtlvl(ACC_PVLVL(iogrp),LVL_PV(iogrp),rnacc*grav,0., &
         cmpflg,ip,'p','pvlvl','Potential vorticity',' ','m-1 s-1')

    if (lvl_difint(iogrp) == 2) then
      call wrtlvl(ACC_DIFINTLVL(iogrp),LVL_DIFINT(iogrp),1.,0., &
           cmpflg,ip,'p','difintlvl', &
           'Layer interface diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlvl(ACC_DIFINTLVL(iogrp),LVL_DIFINT(iogrp),rnacc,0., &
           cmpflg,ip,'p','difintlvl', &
           'Layer interface diffusivity',' ','m2 s-1')
    end if

    if (lvl_difiso(iogrp) == 2) then
      call wrtlvl(ACC_DIFISOLVL(iogrp),LVL_DIFISO(iogrp),1.,0., &
           cmpflg,ip,'p','difisolvl', &
           'Isopycnal diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlvl(ACC_DIFISOLVL(iogrp),LVL_DIFISO(iogrp),rnacc,0., &
           cmpflg,ip,'p','difisolvl', &
           'Isopycnal diffusivity',' ','m2 s-1')
    end if

    if (lvl_difdia(iogrp) == 2) then
      call wrtlvl(ACC_DIFDIALVL(iogrp),LVL_DIFDIA(iogrp),1.,0., &
           cmpflg,ip,'p','difdialvl', &
           'Vertical diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlvl(ACC_DIFDIALVL(iogrp),LVL_DIFDIA(iogrp),rnacc,0., &
           cmpflg,ip,'p','difdialvl', &
           'Vertical diffusivity',' ','m2 s-1')
    end if

    if (lvl_difvmo(iogrp) == 2) then
      call wrtlvl(ACC_DIFVMOLVL(iogrp),LVL_DIFVMO(iogrp),1.,0., &
           cmpflg,ip,'p','difvmolvl', &
           'Vertical momentum diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlvl(ACC_DIFVMOLVL(iogrp),LVL_DIFVMO(iogrp),rnacc,0., &
           cmpflg,ip,'p','difvmolvl', &
           'Vertical momentum diffusivity',' ','m2 s-1')
    end if

    if (lvl_difvho(iogrp) == 2) then
      call wrtlvl(ACC_DIFVHOLVL(iogrp),LVL_DIFVHO(iogrp),1.,0., &
           cmpflg,ip,'p','difvholvl', &
           'Vertical heat diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlvl(ACC_DIFVHOLVL(iogrp),LVL_DIFVHO(iogrp),rnacc,0., &
           cmpflg,ip,'p','difvholvl', &
           'Vertical heat diffusivity',' ','m2 s-1')
    end if

    if (lvl_difvso(iogrp) == 2) then
      call wrtlvl(ACC_DIFVSOLVL(iogrp),LVL_DIFVSO(iogrp),1., 0., &
           cmpflg,ip,'p','difvsolvl', &
           'Vertical salt diffusivity',' ','log10(m2 s-1)')
    else
      call wrtlvl(ACC_DIFVSOLVL(iogrp),LVL_DIFVSO(iogrp),rnacc,0., &
           cmpflg,ip,'p','difvsolvl', &
           'Vertical salt diffusivity',' ','m2 s-1')
    end if

    if (use_TRC .and. use_TKE) then
      call wrtlvl(ACC_TKELVL(iogrp),LVL_TKE(iogrp),rnacc,0., &
           cmpflg,ip,'p','tkelvl','Turbulent kinetic energy',' ','m2 s-2')

      call wrtlvl(ACC_GLS_PSILVL(iogrp),LVL_GLS_PSI(iogrp),rnacc,0., &
           cmpflg,ip,'p','gls_psilvl','Generic length scale',' ','m2 s-3')
    end if

    ! store meridional transports
    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      dimslyr='lat layer region time'
    else
      dimslyr='lat sigma region time'
    end if
    if (msc_mmflxl(iogrp) /= 0) then
      call ncwrt1('mmflxl',dimslyr,mmflxl)
      call ncattr('long_name', &
           'Overturning stream-function on layers')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmflxd(iogrp) /= 0) then
      call ncwrt1('mmflxd','lat depth region time',mmflxd)
      call ncattr('long_name', &
           'Overturning stream-function on z-levels')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmftdl(iogrp) /= 0) then
      call ncwrt1('mmftdl',dimslyr,mmftdl)
      call ncattr('long_name', &
           'Overturning stream-function due to thickness diffusion '// &
           'on layers')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmfsml(iogrp) /= 0) then
      call ncwrt1('mmfsml',dimslyr,mmfsml)
      call ncattr('long_name', &
           'Overturning stream-function due to submesoscale transport '// &
           'on layers')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmftdd(iogrp) /= 0) then
      call ncwrt1('mmftdd','lat depth region time',mmftdd)
      call ncattr('long_name', &
           'Overturning stream-function due to thickness diffusion '// &
           'on z-levels')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmfsmd(iogrp) /= 0) then
      call ncwrt1('mmfsmd','lat depth region time',mmfsmd)
      call ncattr('long_name', &
           'Overturning stream-function due to submesoscale transport '// &
           'on z-levels')
      call ncattr('units','kg s-1')
    end if
    if (msc_mhflx(iogrp) /= 0) then
      call ncwrt1('mhflx','lat region time',mhflx)
      call ncattr('long_name','Meridional heat flux')
      call ncattr('units','W')
    end if
    if (msc_mhftd(iogrp) /= 0) then
      call ncwrt1('mhftd','lat region time',mhftd)
      call ncattr('long_name', &
           'Meridional heat flux due to thickness diffusion')
      call ncattr('units','W')
    end if
    if (msc_mhfsm(iogrp) /= 0) then
      call ncwrt1('mhfsm','lat region time',mhfsm)
      call ncattr('long_name', &
           'Meridional heat flux due to submesoscale transport')
      call ncattr('units','W')
    end if
    if (msc_mhfld(iogrp) /= 0) then
      call ncwrt1('mhfld','lat region time',mhfld)
      call ncattr('long_name', &
           'Meridional heat flux due to lateral diffusion')
      call ncattr('units','W')
    end if
    if (msc_msflx(iogrp) /= 0) then
      call ncwrt1('msflx','lat region time',msflx)
      call ncattr('long_name','Meridional salt flux')
      call ncattr('units','kg s-1')
    end if
    if (msc_msftd(iogrp) /= 0) then
      call ncwrt1('msftd','lat region time',msftd)
      call ncattr('long_name', &
           'Meridional salt flux due to thickness diffusion')
      call ncattr('units','kg s-1')
    end if
    if (msc_msfsm(iogrp) /= 0) then
      call ncwrt1('msfsm','lat region time',msfsm)
      call ncattr('long_name', &
           'Meridional salt flux due to submesoscale transport')
      call ncattr('units','kg s-1')
    end if
    if (msc_msfld(iogrp) /= 0) then
      call ncwrt1('msfld','lat region time',msfld)
      call ncattr('long_name', &
           'Meridional salt flux due to lateral diffusion')
      call ncattr('units','kg s-1')
    end if

    ! store section transports
    if (msc_voltr(iogrp) /= 0) then
      call ncwrt1('voltr','section time',voltr)
      call ncattr('long_name','Section transports')
      call ncattr('units','kg s-1')
    end if

    ! store global sums and averages
    if (msc_massgs(iogrp) /= 0) then
      call ncwrt1('massgs','time',massgs)
      call ncattr('long_name','Sea water mass')
      call ncattr('units','kg')
    end if
    if (msc_volgs(iogrp) /= 0) then
      call ncwrt1('volgs','time',volgs)
      call ncattr('long_name','Sea water volume')
      call ncattr('units','m3')
    end if
    if (msc_salnga(iogrp) /= 0) then
      call ncwrt1('salnga','time',salnga)
      call ncattr('long_name','Global average salinity')
      call ncattr('units','g kg-1')
    end if
    if (msc_tempga(iogrp) /= 0) then
      call ncwrt1('tempga','time',tempga)
      call ncattr('long_name','Global average temperature')
      call ncattr('units','degC')
    end if
    if (msc_sssga(iogrp) /= 0) then
      call ncwrt1('sssga','time',sssga)
      call ncattr('long_name','Global average sea surface salinity')
      call ncattr('units','g kg-1')
    end if
    if (msc_sstga(iogrp) /= 0) then
      call ncwrt1('sstga','time',sstga)
      call ncattr('long_name', &
           'Global average sea surface temperature')
      call ncattr('units','degC')
    end if

    if (use_TRC) then
      if (lyr_idlage(iogrp) /= 0.or.lyr_trc(iogrp) /= 0) then
        call inilyr(ACC_UTILLYR(1),'p',0.)
        call acclyr(ACC_UTILLYR,dp(1-nbdy,1-nbdy,k1m),tmp3d,0,'p')
        call wrtlyr(ACC_UTILLYR(1),max(LYR_IDLAGE(iogrp),LYR_TRC(iogrp)), &
                    1.,0.,cmpflg,ip,'p','dp_trc', &
                    'Layer pressure thickness',' ','Pa')
      end if
      if (use_IDLAGE) then

        ! ideal age tracer
        if (lyr_idlage(iogrp) /= 0) then
          call inilyr(ACC_UTILLYR(1),'p',0.)
          call acclyr(ACC_UTILLYR,trc(1-nbdy,1-nbdy,k1m,itriag),tmp3d,0, &
               'p')
          call wrtlyr(ACC_UTILLYR(1),LYR_IDLAGE(iogrp),1.,0., &
               cmpflg,ip,'p','idlage','Ideal age', &
               'sea_water_age_since_surface_contact','year')
        end if
        if (lvl_idlage(iogrp) /= 0) then
          if (vcoord_tag == vcoord_isopyc_bulkml) then
            call inilvl(ACC_IDLAGELVL(iogrp),'p',0.)
            do k = 1,kk
              call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
              call acclvl(ACC_IDLAGELVL,trc(1-nbdy,1-nbdy,k1m,itriag),'p', &
                   k,ind1,ind2,wghts)
            end do
          end if
          call msklvl(ACC_IDLAGELVL(iogrp),'p')
          call wrtlvl(ACC_IDLAGELVL(iogrp),LVL_IDLAGE(iogrp), &
               1.,0.,cmpflg,ip,'p','idlagelvl','Ideal age', &
               'sea_water_age_since_surface_contact','year')
        end if
      end if ! use_IDLAGE

      ! ocean tracers
      if (lyr_trc(iogrp) > 0.and.ntrocn > 0) then
        if (use_ATRC) then
          do nt = 1,ntrocn-natr
            call inilyr(ACC_UTILLYR(1),'p',0.)
            call acclyr(ACC_UTILLYR,trc(1-nbdy,1-nbdy,k1m,nt),tmp3d,0,'p')
            write (trcnm,'(a,i3.3)') 'trc',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call wrtlyr(ACC_UTILLYR(1),LYR_TRC(iogrp),1.,0., &
                 cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
          end do
          do nt = 1,natr
            nat = ntr-natr+nt
            call inilyr(ACC_UTILLYR(1),'p',0.)
            !$omp parallel do private(k,km,l,i)
            do j = 1,jj
              do k = 1,kk
                km = k+mm
                do l = 1,isp(j)
                  do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                    tmp3d(i,j,k) = trc(i,j,km,nat)/max(trc(i,j,km,nt),treps)
                  end do
                end do
              end do
            end do
            !$omp end parallel do
            call acclyr(ACC_UTILLYR,tmp3d,tmp3d,0,'p')
            write (trcnm,'(a,i3.3)') 'atrc',nt
            write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
            call wrtlyr(ACC_UTILLYR(1),LYR_TRC(iogrp),1.,0., &
                 cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
          end do
        else
          do nt = 1,ntrocn
            call inilyr(ACC_UTILLYR(1),'p',0.)
            call acclyr(ACC_UTILLYR,trc(1-nbdy,1-nbdy,k1m,nt),tmp3d,0,'p')
            write (trcnm,'(a,i3.3)') 'trc',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call wrtlyr(ACC_UTILLYR(1),LYR_TRC(iogrp),1.,0., &
                 cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
          end do
        end if
      end if
      if (lvl_trc(iogrp) > 0.and.ntrocn > 0) then
        if (use_ATRC) then
          do nt = 1,ntrocn-natr
            call inilvl(ACC_UTILLVL(1),'p',0.)
            do k = 1,kk
              call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
              call acclvl(ACC_UTILLVL,trc(1-nbdy,1-nbdy,k1m,nt),'p', &
                   k,ind1,ind2,wghts)
            end do
            call msklvl(ACC_UTILLVL(1),'p')
            write (trcnm,'(a,i3.3)') 'trclvl',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call wrtlvl(ACC_UTILLVL(1),LVL_TRC(iogrp), &
                 1.,0.,cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
          end do
          do nt = 1,natr
            nat = ntr-natr+nt
            call inilvl(ACC_UTILLVL(1),'p',0.)
            !$omp parallel do private(k,km,l,i)
            do j = 1,jj
              do k = 1,kk
                km = k+mm
                do l = 1,isp(j)
                  do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                    tmp3d(i,j,k) = trc(i,j,km,nat)/max(trc(i,j,km,nt),treps)
                  end do
                end do
              end do
            end do
            !$omp end parallel do
            do k = 1,kk
              call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
              call acclvl(ACC_UTILLVL,tmp3d,'p',k,ind1,ind2,wghts)
            end do
            call msklvl(ACC_UTILLVL(1),'p')
            write (trcnm,'(a,i3.3)') 'atrclvl',nt
            write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
            call wrtlvl(ACC_UTILLVL(1),LVL_TRC(iogrp), &
                 1.,0.,cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
          end do
        else
          do nt = 1,ntrocn
            call inilvl(ACC_UTILLVL(1),'p',0.)
            do k = 1,kk
              call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
              call acclvl(ACC_UTILLVL,trc(1-nbdy,1-nbdy,k1m,nt),'p', &
                   k,ind1,ind2,wghts)
            end do
            call msklvl(ACC_UTILLVL(1),'p')
            write (trcnm,'(a,i3.3)') 'trclvl',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call wrtlvl(ACC_UTILLVL(1),LVL_TRC(iogrp), &
                 1.,0.,cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
          end do
        end if
      end if
    end if

    ! close netcdf file
    call ncfcls

    ! initialisation of fields
    call inifld(iogrp)

    ! reset accumulation counter
    nacc_phy(iogrp) = 0

  end subroutine diaout_iogrp



  subroutine diasec(iogrp)

    ! Arguments
    integer :: iogrp

    ! Local variables
    integer       :: nfu,iostatus,n,i,j,k,s,l
    integer, save :: nsi(max_sec)
    integer, save :: isi(max_sec,sdm)
    integer, save :: jsi(max_sec,sdm)
    integer, save :: usi(max_sec,sdm)
    integer, save :: vsi(max_sec,sdm)
    integer, save :: equat_sec
    character(len = 120) :: char120
    logical :: iniflg = .true.
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         uflx_cum,vflx_cum,uflx_cum350,vflx_cum350
    real, dimension(itdm,jtdm) :: &
         uflx_cumt,vflx_cumt,uflx_cum350t,vflx_cum350t
    real(8) :: volu,volv

    !---------------------------------------------------------------
    ! read section information
    !---------------------------------------------------------------
    if (iniflg) then
      if (mnproc == 1) then
        equat_sec = -1
        open(newunit=nfu,file=sec_sifile,status = 'old')
        sec_num = 0
        do
          read(nfu,'(a120)',iostat = iostatus) char120
          if (iostatus /= 0) exit
          if (char120(1:4) == 'Name') then
            sec_num = sec_num+1
            sec_name(sec_num) = char120(7:120)
            if (index(char120,'equatorial_undercurrent') > 0) then
              equat_sec = sec_num
            end if
            nsi(sec_num) = 0
          else
            nsi(sec_num) = nsi(sec_num)+1
            read(char120,*) isi(sec_num,nsi(sec_num)), &
                 jsi(sec_num,nsi(sec_num)), &
                 usi(sec_num,nsi(sec_num)), &
                 vsi(sec_num,nsi(sec_num))
          end if
        end do
        close(nfu)
        write(lp,*) 'number of sections = ',sec_num
      end if
      call xcbcst(sec_num)
      iniflg = .false.
    end if

    ! Prepare 2d field
    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        uflx_cum(i,j) = 0.
        vflx_cum(i,j) = 0.
        uflx_cum350(i,j) = 0.
        vflx_cum350(i,j) = 0.
      end do
    end do
    !$omp end parallel do

    ! Compute accumulated transports
    !$omp parallel do private(k,l,i)
    do j = 1,jj
      do k = 1,ddm
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            uflx_cum(i,j) = uflx_cum(i,j)+ &
                 phylvl(i,j,k,ACC_UFLXLVL(iogrp)) &
                 *.5/(grav*baclin*nacc_phy(iogrp))
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vflx_cum(i,j) = vflx_cum(i,j)+ &
                 phylvl(i,j,k,ACC_VFLXLVL(iogrp)) &
                 *.5/(grav*baclin*nacc_phy(iogrp))
          end do
        end do

        ! the upper 350 m  for equatorial_undercurrent
        if (k == k350-1) then
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              uflx_cum350(i,j) = uflx_cum(i,j)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vflx_cum350(i,j) = vflx_cum(i,j)
            end do
          end do
        else if (k == k350) then
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              uflx_cum350(i,j) = uflx_cum350(i,j)+ &
                   w350*(uflx_cum(i,j)-uflx_cum350(i,j))
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vflx_cum350(i,j) = vflx_cum350(i,j)+ &
                   w350*(vflx_cum(i,j)-vflx_cum350(i,j))
            end do
          end do
        end if
      end do
    end do
    !$omp end parallel do

    ! Collect data on master node
    call xcaget(uflx_cumt,uflx_cum,1)
    call xcaget(vflx_cumt,vflx_cum,1)
    call xcaget(uflx_cum350t,uflx_cum350,1)
    call xcaget(vflx_cum350t,vflx_cum350,1)

    ! Compute section transports
    if (mnproc == 1) then
      do s = 1,sec_num
        voltr(s) = 0.
        if (s == equat_sec) then
          do n = 1,nsi(s)
            i = isi(s,n)
            j = jsi(s,n)
            volu = uflx_cum350t(i,j)*real(usi(s,n))
            volv = vflx_cum350t(i,j)*real(vsi(s,n))
            voltr(s) = voltr(s)+volu+volv
          end do
        else
          do n = 1,nsi(s)
            i = isi(s,n)
            j = jsi(s,n)
            volu = uflx_cumt(i,j)*real(usi(s,n))
            volv = vflx_cumt(i,j)*real(vsi(s,n))
            voltr(s) = voltr(s)+volu+volv
          end do
        end if
      end do
    end if

  end subroutine diasec


  subroutine diamer(iogrp)

    ! Arguments
    integer, intent(in) :: iogrp

    ! Local variables
    integer :: ncid,dimid,varid,i,j,k,l,m,n,o,s,ocn_nreg,nfu,iostatus
    integer :: istat,iind1,jind1,uflg1,vflg1,nind1
    integer :: nfld,ACC_UIND,ACC_VIND,nind(ldm),iind(sdm,ldm),jind(sdm,ldm)
    integer :: kmxl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer(i2) :: uflg(sdm,ldm),vflg(sdm,ldm),oflg(itdm,jtdm)
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ucum,vcum
    real, dimension(itdm,jtdm) :: depthst,ucumg,vcumg
    real, save, allocatable, dimension(:,:) ::  mflx_or,mflx_mr,mflx_last_mr
    integer, save, allocatable, dimension(:,:) :: mcnt_or,mcnt_mr,kmax
    real :: r
    character :: c20*20
    logical :: iniflg = .true.

    save nind,iind,jind,oflg,uflg,vflg,depthst,iniflg,ocn_nreg

    if (iniflg) then

      if (mnproc == 1) then

        ! Read ocean region flags from ocean_regions.nc
        call ncerro(nf90_open(mer_orfile,nf90_nowrite,ncid))
        call ncerro(nf90_inq_dimid(ncid,'x',dimid))
        call ncerro(nf90_inquire_dimension(ncid,dimid,len = i))
        call ncerro(nf90_inq_dimid(ncid,'y',dimid))
        call ncerro(nf90_inquire_dimension(ncid,dimid,len = j))
        if (i /= itdm.or.j /= jtdm) then
          write (lp,'(2a)') ' wrong dimensions in',trim(mer_orfile)
          call xchalt('(diamer)')
          stop '(diamer)'
        end if
        call ncerro(nf90_inq_dimid(ncid,'regions',dimid))
        call ncerro(nf90_inquire_dimension(ncid,dimid,len = ocn_nreg))
        call ncerro(nf90_inq_varid(ncid,'region',varid))
        call ncerro(nf90_get_var(ncid,varid,oflg))
        call ncerro(nf90_close(ncid))

        ! Read section file metra_index.dat
        open(newunit=nfu,file=mer_mifile,status = 'old')
        lmax = 0
        do l = 1,ldm
          c20 = ' '
          read(nfu,'(a)',iostat = iostatus) c20
          if (iostatus /= 0) exit
          if (c20(1:7) == 'Section') then
            lmax = l
            read(c20(9:15),'(f7.3)') mtlat(l)
            read(c20(16:20),'(i7)') nind1
          else
            write(lp,'(2a)') ' problem in ',trim(mer_mifile)
            call xchalt('(diamer)')
            stop '(diamer)'
          end if
          nind(l) = 0
          do s = 1,nind1
            read(nfu,'(a)',iostat = iostatus) c20
            if (iostatus /= 0) then
              write(lp,*) 'section too short?'
              call xchalt('(diamer)')
              stop '(diamer)'
            end if
            read(c20,*) iind1,jind1,uflg1,vflg1
            if (oflg(iind1,jind1) > 0) then
              nind(l) = nind(l)+1
              iind(nind(l),l) = iind1
              jind(nind(l),l) = jind1
              uflg(nind(l),l) = int(uflg1,i2)
              vflg(nind(l),l) = int(vflg1,i2)
              if (iind(nind(l),l) < 1.or. &
                   iind(nind(l),l) > itdm.or. &
                   jind(nind(l),l) < 1.or. &
                   jind(nind(l),l) > jtdm) then
                write(lp,*) 'iind=',iind(nind(l),l),' itdm = ',itdm
                write(lp,*) 'jind=',jind(nind(l),l),' jtdm = ',jtdm
                call flush(lp)
                write(lp,*) 'indices out of range!'
                call xchalt('(diamer)')
                stop '(diamer)'
              end if
            end if
          end do
          if (lmax > ldm) then
            write(lp,*) 'too many or too long sections!'
            call xchalt('(diamer)')
            stop '(diamer)'
          end if
        end do
        close(nfu)

        ! Allocate utility arrays for meridional fluxes
        allocate(mflx_or(lmax,ocn_nreg),mflx_mr(lmax,mer_nreg), &
             mflx_last_mr(lmax,mer_nreg), &
             mcnt_or(lmax,ocn_nreg),mcnt_mr(lmax,mer_nreg), &
             kmax(lmax,mer_nreg), &
             stat = istat)
        if (istat /= 0) then
          write (lp,*) 'Cannot allocate enough memory!'
          call xchalt('(diamer)')
          stop '(diamer)'
        end if

      end if

      call xcbcst(lmax)

      ! - Allocate arrays for meridional fluxes
      allocate(mmflxl(lmax,kdm,mer_nreg),mmftdl(lmax,kdm,mer_nreg), &
           mmfsml(lmax,kdm,mer_nreg),mmflxd(lmax,ddm,mer_nreg), &
           mmftdd(lmax,ddm,mer_nreg),mmfsmd(lmax,ddm,mer_nreg), &
           mhflx(lmax,mer_nreg),mhftd(lmax,mer_nreg), &
           mhfsm(lmax,mer_nreg),mhfld(lmax,mer_nreg), &
           msflx(lmax,mer_nreg),msftd(lmax,mer_nreg), &
           msfsm(lmax,mer_nreg),msfld(lmax,mer_nreg), &
           stat = istat)
      if (istat /= 0) then
        write (lp,*) 'Cannot allocate enough memory!'
        call xchalt('(diamer)')
        stop '(diamer)'
      end if

    end if

    ! Compute vertical integrated heat and salt transports

    !$omp parallel do private(i)
    do j = 1,jj
      do i = 1,ii
        ucum(i,j) = 0.
        vcum(i,j) = 0.
      end do
    end do
    !$omp end parallel do

    do nfld = 1,8

      if     (nfld == 1) then
        if (acc_mhflx(iogrp) == 0) cycle
        ACC_UIND = ACC_UTFLX(iogrp)
        ACC_VIND = ACC_VTFLX(iogrp)
        r = spcifh*.5/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 2) then
        if (acc_mhftd(iogrp) == 0) cycle
        ACC_UIND = ACC_UTFLTD(iogrp)
        ACC_VIND = ACC_VTFLTD(iogrp)
        r = spcifh*.5/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 3) then
        if (acc_mhfsm(iogrp) == 0) cycle
        ACC_UIND = ACC_UTFLSM(iogrp)
        ACC_VIND = ACC_VTFLSM(iogrp)
        r = spcifh*.5/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 4) then
        if (acc_mhfld(iogrp) == 0) cycle
        ACC_UIND = ACC_UTFLLD(iogrp)
        ACC_VIND = ACC_VTFLLD(iogrp)
        r = spcifh*.5/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 5) then
        if (acc_msflx(iogrp) == 0) cycle
        ACC_UIND = ACC_USFLX(iogrp)
        ACC_VIND = ACC_VSFLX(iogrp)
        r = .5*g2kg/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 6) then
        if (acc_msftd(iogrp) == 0) cycle
        ACC_UIND = ACC_USFLTD(iogrp)
        ACC_VIND = ACC_VSFLTD(iogrp)
        r = .5*g2kg/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 7) then
        if (acc_msfsm(iogrp) == 0) cycle
        ACC_UIND = ACC_USFLSM(iogrp)
        ACC_VIND = ACC_VSFLSM(iogrp)
        r = .5*g2kg/(grav*baclin*nacc_phy(iogrp))
      else if (nfld == 8) then
        if (acc_msfld(iogrp) == 0) cycle
        ACC_UIND = ACC_USFLLD(iogrp)
        ACC_VIND = ACC_VSFLLD(iogrp)
        r = .5*g2kg/(grav*baclin*nacc_phy(iogrp))
      else
        write(lp,*) 'field index out of range'
        call xchalt('(diamer)')
        stop '(diamer)'
      end if

      !$omp parallel do private(l,i,k)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            ucum(i,j) = phylyr(i,j,1,ACC_UIND)*r
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vcum(i,j) = phylyr(i,j,1,ACC_VIND)*r
          end do
        end do
        do k = 2,kk
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              ucum(i,j) = ucum(i,j)+phylyr(i,j,k,ACC_UIND)*r
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vcum(i,j) = vcum(i,j)+phylyr(i,j,k,ACC_VIND)*r
            end do
          end do
        end do
      end do
      !$omp end parallel do
      call xcaget(ucumg,ucum,1)
      call xcaget(vcumg,vcum,1)
      if (mnproc == 1) then
        do l = 1,lmax
          ! Accumulate meridional fluxes in seperate ocean regions
          do o = 1,ocn_nreg
            mflx_or(l,o) = 0.
            mcnt_or(l,o) = 0
          end do
          do s = 1,nind(l)
            i = iind(s,l)
            j = jind(s,l)
            o = oflg(i,j)
            mflx_or(l,o) = mflx_or(l,o) &
                 +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
            mcnt_or(l,o) = mcnt_or(l,o)+1
          end do
          ! Add together the ocean regions that belong to the regions of
          ! meridional flux diagnostics and apply a fill value outside
          ! the regions latitude bounds or when no values have been
          ! accumulated
          do m = 1,mer_nreg
            if (mtlat(l) < mer_minlat(m).or. &
                mtlat(l) > mer_maxlat(m)) then
              mflx_mr(l,m) = nf90_fill_double
            else
              mflx_mr(l,m) = 0.
              mcnt_mr(l,m) = 0
              if (mer_regflg(m,1) == 0) then
                do o = 1,ocn_nreg
                  mflx_mr(l,m) = mflx_mr(l,m)+mflx_or(l,o)
                  mcnt_mr(l,m) = mcnt_mr(l,m)+mcnt_or(l,o)
                end do
              else
                do n = 1,mer_nflg(m)
                  o = mer_regflg(m,n)
                  mflx_mr(l,m) = mflx_mr(l,m)+mflx_or(l,o)
                  mcnt_mr(l,m) = mcnt_mr(l,m)+mcnt_or(l,o)
                end do
              end if
              if (mcnt_mr(l,m) == 0) mflx_mr(l,m) = nf90_fill_double
            end if
          end do
        end do

        if     (nfld == 1) then
          do l = 1,lmax
            do m = 1,mer_nreg
              mhflx(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 2) then
          do l = 1,lmax
            do m = 1,mer_nreg
              mhftd(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 3) then
          do l = 1,lmax
            do m = 1,mer_nreg
              mhfsm(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 4) then
          do l = 1,lmax
            do m = 1,mer_nreg
              mhfld(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 5) then
          do l = 1,lmax
            do m = 1,mer_nreg
              msflx(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 6) then
          do l = 1,lmax
            do m = 1,mer_nreg
              msftd(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 7) then
          do l = 1,lmax
            do m = 1,mer_nreg
              msfsm(l,m) = mflx_mr(l,m)
            end do
          end do
        else if (nfld == 8) then
          do l = 1,lmax
            do m = 1,mer_nreg
              msfld(l,m) = mflx_mr(l,m)
            end do
          end do
        else
          write(lp,*) 'field index out of range'
          call xchalt('(diamer)')
          stop '(diamer)'
        end if

      end if
    end do

    ! Compute overturning stream function at layer interfaces

    r = 1./real(nacc_phy(iogrp))
    kmxl = 0
    !$omp parallel do private(l,i,k)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          k = 1
          do while (phyh2d(i,j,acc_sigmx(iogrp))*r >= sigmar(i,j,k))
            k = k+1
            if (k == kk) exit
          end do
          kmxl(i,j) = k
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j) = kmxl(i,j)
        end do
      end do
    end do
    !$omp end parallel do
    call xctilr(util1, 1,1, 1,1, halo_ps)
    !$omp parallel do private(l,i)
    do j = 0,jj+1
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmxl(i,j) = nint(util1(i,j))
        end do
      end do
    end do
    !$omp end parallel do

    r = .5/(grav*baclin*nacc_phy(iogrp))

    do nfld = 1,3

      if     (nfld == 1) then
        if (acc_mmflxl(iogrp) == 0) cycle
        ACC_UIND = ACC_UFLX(iogrp)
        ACC_VIND = ACC_VFLX(iogrp)
      else if (nfld == 2) then
        if (acc_mmftdl(iogrp) == 0) cycle
        ACC_UIND = ACC_UMFLTD(iogrp)
        ACC_VIND = ACC_VMFLTD(iogrp)
      else if (nfld == 3) then
        if (acc_mmfsml(iogrp) == 0) cycle
        ACC_UIND = ACC_UMFLSM(iogrp)
        ACC_VIND = ACC_VMFLSM(iogrp)
      else
        write(lp,*) 'field index out of range'
        call xchalt('(diamer)')
        stop '(diamer)'
      end if

      if (mnproc == 1) then
        do l = 1,lmax
          do m = 1,mer_nreg
            mflx_last_mr(l,m) = 0.
          end do
        end do
      end if
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            ucum(i,j) = 0.
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vcum(i,j) = 0.
          end do
        end do
      end do
      !$omp end parallel do
      do k = 1,kk
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (k > 2) then
                ucum(i,j) = ucum(i,j)+r*phylyr(i,j,k,ACC_UIND)
              end if
              if (k == min(kmxl(i-1,j),kmxl(i,j))) then
                ucum(i,j) = ucum(i,j)+r*(phylyr(i,j,1,ACC_UIND) &
                     +phylyr(i,j,2,ACC_UIND))
              end if
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (k > 2) then
                vcum(i,j) = vcum(i,j)+r*phylyr(i,j,k,ACC_VIND)
              end if
              if (k == min(kmxl(i,j-1),kmxl(i,j))) then
                vcum(i,j) = vcum(i,j)+r*(phylyr(i,j,1,ACC_VIND) &
                     +phylyr(i,j,2,ACC_VIND))
              end if
            end do
          end do
        end do
        !$omp end parallel do
        call xcaget(ucumg,ucum,1)
        call xcaget(vcumg,vcum,1)
        if (mnproc == 1) then
          do l = 1,lmax
            ! Accumulate meridional fluxes in seperate ocean regions
            do o = 1,ocn_nreg
              mflx_or(l,o) = 0.
            end do
            do s = 1,nind(l)
              i = iind(s,l)
              j = jind(s,l)
              o = oflg(i,j)
              mflx_or(l,o) = mflx_or(l,o)&
                            +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
            end do
            ! Add together the ocean regions that belong to the regions
            ! of meridional flux diagnostics and apply a fill value
            ! outside the regions latitude bounds or when no values have
            ! been accumulated
            do m = 1,mer_nreg
              if (mtlat(l) < mer_minlat(m).or. &
                  mtlat(l) > mer_maxlat(m)) then
                mflx_mr(l,m) = nf90_fill_double
              else
                mflx_mr(l,m) = 0.
                if (mer_regflg(m,1) == 0) then
                  do o = 1,ocn_nreg
                    mflx_mr(l,m) = mflx_mr(l,m)+mflx_or(l,o)
                  end do
                else
                  do n = 1,mer_nflg(m)
                    o = mer_regflg(m,n)
                    mflx_mr(l,m) = mflx_mr(l,m)+mflx_or(l,o)
                  end do
                end if
                if (abs(mflx_mr(l,m)-mflx_last_mr(l,m)) < &
                     1.e5*epsilp) then
                  mflx_last_mr(l,m) = mflx_mr(l,m)
                  mflx_mr(l,m) = nf90_fill_double
                else
                  mflx_last_mr(l,m) = mflx_mr(l,m)
                end if
              end if
            end do
          end do

          if     (nfld == 1) then
            do l = 1,lmax
              do m = 1,mer_nreg
                mmflxl(l,k,m) = mflx_mr(l,m)
              end do
            end do
          else if (nfld == 2) then
            do l = 1,lmax
              do m = 1,mer_nreg
                mmftdl(l,k,m) = mflx_mr(l,m)
              end do
            end do
          else if (nfld == 3) then
            do l = 1,lmax
              do m = 1,mer_nreg
                mmfsml(l,k,m) = mflx_mr(l,m)
              end do
            end do
          else
            write(lp,*) 'field index out of range'
            call xchalt('(diamer)')
            stop '(diamer)'
          end if

        end if
      end do
    end do

    ! Compute overturning stream function at levitus level interfaces
    ! Prepare depth mask

    if (iniflg) call xcaget(depthst,depths,1)
    if (iniflg.and.mnproc == 1) then
      do l = 1,lmax
        do m = 1,mer_nreg
          kmax(l,m) = 0
        end do
      end do
      do k = 1,ddm
        do l = 1,lmax
          do s = 1,nind(l)
            i = iind(s,l)
            j = jind(s,l)
            if (depthslev_bnds(1,k) < depthst(i,j)) then
              do m = 1,mer_nreg
                if (mer_regflg(m,1) == 0) then
                  kmax(l,m) = k
                else
                  do n = 1,mer_nflg(m)
                    if (mer_regflg(m,n) == oflg(i,j)) kmax(l,m) = k
                  end do
                end if
              end do
            end if
          end do
        end do
      end do
    end if

    r = .5/(grav*baclin*nacc_phy(iogrp))

    do nfld = 1,3

      if     (nfld == 1) then
        if (acc_mmflxd(iogrp) == 0) cycle
        ACC_UIND = ACC_UFLXLVL(iogrp)
        ACC_VIND = ACC_VFLXLVL(iogrp)
      else if (nfld == 2) then
        if (acc_mmftdd(iogrp) == 0) cycle
        ACC_UIND = ACC_UMFLTDLVL(iogrp)
        ACC_VIND = ACC_VMFLTDLVL(iogrp)
      else if (nfld == 3) then
        if (acc_mmfsmd(iogrp) == 0) cycle
        ACC_UIND = ACC_UMFLSMLVL(iogrp)
        ACC_VIND = ACC_VMFLSMLVL(iogrp)
      else
        write(lp,*) 'field index out of range'
        call xchalt('(diamer)')
        stop '(diamer)'
      end if

      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            ucum(i,j) = 0.
          end do
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            vcum(i,j) = 0.
          end do
        end do
      end do
      !$omp end parallel do
      do k = 1,ddm
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              ucum(i,j) = ucum(i,j)+r*phylvl(i,j,k,ACC_UIND)
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vcum(i,j) = vcum(i,j)+r*phylvl(i,j,k,ACC_VIND)
            end do
          end do
        end do
        !$omp end parallel do
        call xcaget(ucumg,ucum,1)
        call xcaget(vcumg,vcum,1)
        if (mnproc == 1) then
          do l = 1,lmax
            ! Accumulate meridional fluxes in seperate ocean regions
            do o = 1,ocn_nreg
              mflx_or(l,o) = 0.
            end do
            do s = 1,nind(l)
              i = iind(s,l)
              j = jind(s,l)
              o = oflg(i,j)
              mflx_or(l,o) = mflx_or(l,o) &
                   +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
            end do
            ! Add together the ocean regions that belong to the regions
            ! of meridional flux diagnostics and apply a fill value
            ! outside the regions latitude and depth bounds
            do m = 1,mer_nreg
              if (mtlat(l) < mer_minlat(m).or. &
                   mtlat(l) > mer_maxlat(m).or. &
                   kmax(l,m) < k) then
                mflx_mr(l,m) = nf90_fill_double
              else
                mflx_mr(l,m) = 0.
                if (mer_regflg(m,1) == 0) then
                  do o = 1,ocn_nreg
                    mflx_mr(l,m) = mflx_mr(l,m)+mflx_or(l,o)
                  end do
                else
                  do n = 1,mer_nflg(m)
                    o = mer_regflg(m,n)
                    mflx_mr(l,m) = mflx_mr(l,m)+mflx_or(l,o)
                  end do
                end if
              end if
            end do
          end do

          if     (nfld == 1) then
            do l = 1,lmax
              do m = 1,mer_nreg
                mmflxd(l,k,m) = mflx_mr(l,m)
              end do
            end do
          else if (nfld == 2) then
            do l = 1,lmax
              do m = 1,mer_nreg
                mmftdd(l,k,m) = mflx_mr(l,m)
              end do
            end do
          else if (nfld == 3) then
            do l = 1,lmax
              do m = 1,mer_nreg
                mmfsmd(l,k,m) = mflx_mr(l,m)
              end do
            end do
          else
            write(lp,*) 'field index out of range'
            call xchalt('(diamer)')
            stop '(diamer)'
          end if

        end if
      end do
    end do

    if (iniflg) iniflg = .false.

  end subroutine diamer



  subroutine diavfl(iogrp,m,n,mm,nn,k1m,k1n)
  !---------------------------------------------------------------
  ! computation of vertical mass flux at layer interfaces
  !---------------------------------------------------------------

    integer :: iogrp,m,n,mm,nn,k1m,k1n

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: wflx,ucum,vcum
    integer :: i,j,k,km,kn,l
    real :: q

    !
    ! Compute vertical mass flux at layer interfaces
    !
    if (acc_wflx(iogrp)+acc_wflx2(iogrp) /= 0) then

      !$omp parallel do private(i)
      do j = 1,jj
        do i = 1,ii
          wflx(i,j) = 0.
        end do
      end do
      !$omp end parallel do
      do k = kk,1,-1
        km = k+mm
        kn = k+nn
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              wflx(i,j) = wflx(i,j) &
                   -(uflx(i+1,j,kn)-uflx(i,j,kn) &
                   +vflx(i,j+1,kn)-vflx(i,j,kn)) &
                   -(dp(i,j,km)-dpold(i,j,km))*scp2(i,j)
              phylyr(i,j,k,ACC_WFLX(iogrp))= &
                   phylyr(i,j,k,ACC_WFLX(iogrp))+wflx(i,j)
              phylyr(i,j,k,ACC_WFLX2(iogrp))= &
                   phylyr(i,j,k,ACC_WFLX2(iogrp))+wflx(i,j)**2
            end do
          end do
        end do
        !$omp end parallel do
      end do
    end if

    ! Computation of vertical mass flux at levitus layer interfaces
    if (acc_wflxlvl(iogrp)+acc_wflx2lvl(iogrp) /= 0) then

      call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(iogrp)), &
           1,ddm, 1,1, halo_uv)
      call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(iogrp)), &
           1,ddm, 1,1, halo_vv)
      !$omp parallel do private(i)
      do j = 1,jj+1
        do i = 1,ii+1
          ucum(i,j) = 0.
          vcum(i,j) = 0.
        end do
      end do
      !$omp end parallel do
      do k = ddm,1,-1
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              ucum(i,j) = ucum(i,j) &
                   +phylvl(i,j,k,ACC_UFLXLVL(iogrp)) &
                   -phylvl(i,j,k,ACC_UFLXOLD(iogrp))
            end do
          end do
        end do
        !$omp end parallel do
        !$omp parallel do private(l,i)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              vcum(i,j) = vcum(i,j) &
                   +phylvl(i,j,k,ACC_VFLXLVL(iogrp)) &
                   -phylvl(i,j,k,ACC_VFLXOLD(iogrp))
            end do
          end do
        end do
        !$omp end parallel do
        !$omp parallel do private(l,i,q)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              q = -(ucum(i+1,j)-ucum(i,j) &
                   +vcum(i,j+1)-vcum(i,j))
              phylvl(i,j,k,ACC_WFLXLVL(iogrp))= &
                   phylvl(i,j,k,ACC_WFLXLVL(iogrp))+q
              phylvl(i,j,k,ACC_WFLX2LVL(iogrp))= &
                   phylvl(i,j,k,ACC_WFLX2LVL(iogrp))+q**2
            end do
          end do
        end do
        !$omp end parallel do
      end do
    end if

  end subroutine diavfl



  subroutine diazlv(gridid,k,mm,nn,ind1,ind2,weights,weightsflx)

    ! Arguments
    character(len=*), intent(in) :: gridid
    integer,          intent(in) :: k
    integer,          intent(in) :: mm
    integer,          intent(in) :: nn
    integer,          intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     :: ind1
    integer,          intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     :: ind2
    real,             intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: weights
    real,             intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: weightsflx

    ! Local variables
    integer :: d,i,j,l,kl,km,kn,kml,k1m
    real    :: r,dzeps,dpeps,flxeps
    logical :: iniflg = .true.
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk), save  :: ztop
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk), save  :: zbot
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm),save :: dlevp
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm),save :: dlevu
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm),save :: dlevv

    ! Define thresholds
    dzeps = 1e1*epsilp
    dpeps = 1e5*epsilp
    flxeps = 1e5*epsilp

    ! Sort out stuff related to time stepping
    km = k+mm
    kn = k+nn
    k1m = 1+mm

    ! Adjust bounds of levitus levels according to model bathymetry
    if (iniflg) then
      !$omp parallel do private(d,l,i)
      do j = 1,jj+1
        do d = 1,ddm
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              dlevp(i,j,d) = max(dzeps,min(pbath(i,j), &
                   depthslev_bnds(2,d))-depthslev_bnds(1,d))
            end do
          end do
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              dlevu(i,j,d) = max(dzeps,min(ubath(i,j), &
                   depthslev_bnds(2,d))-depthslev_bnds(1,d))
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              dlevv(i,j,d) = max(dzeps,min(vbath(i,j), &
                   depthslev_bnds(2,d))-depthslev_bnds(1,d))
            end do
          end do
        end do
      end do
      !$omp end parallel do
      iniflg = .false.
    end if


    ! Compute top and bottom depths of density layers
    if (k == 1) then
      if (gridid == 'p') then
        !$omp parallel do private(l,i,kl,kml)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,1) = dp(i,j,k1m)
            end do
          end do
          do kl = 2,kk
            kml = kl+mm
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,kl) = zbot(i,j,kl-1)+dp(i,j,kml)
              end do
            end do
          end do
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,1) = zbot(i,j,1)*pbath(i,j)/zbot(i,j,kk)
              ztop(i,j,1) = 0.
              ind1(i,j) = 1
            end do
          end do
          do kl = 2,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,kl) = zbot(i,j,kl)*pbath(i,j)/zbot(i,j,kk)
                ztop(i,j,kl) = zbot(i,j,kl-1)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'u') then
        !$omp parallel do private(l,i,kl,kml)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              zbot(i,j,1) = max(dpu(i,j,k1m),dpeps)
            end do
          end do
          do kl = 2,kk
            kml = kl+mm
            do l = 1,isu(j)
              do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                zbot(i,j,kl) = zbot(i,j,kl-1)+max(dpu(i,j,kml),dpeps)
              end do
            end do
          end do
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              zbot(i,j,1) = zbot(i,j,1)*ubath(i,j)/zbot(i,j,kk)
              ztop(i,j,1) = 0.
              ind1(i,j) = 1
            end do
          end do
          do kl = 2,kk
            do l = 1,isu(j)
              do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                zbot(i,j,kl) = zbot(i,j,kl)*ubath(i,j)/zbot(i,j,kk)
                ztop(i,j,kl) = zbot(i,j,kl-1)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'v') then
        !$omp parallel do private(l,i,kl,kml)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              zbot(i,j,1) = max(dpv(i,j,k1m),dpeps)
            end do
          end do
          do kl = 2,kk
            kml = kl+mm
            do l = 1,isv(j)
              do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                zbot(i,j,kl) = zbot(i,j,kl-1)+max(dpv(i,j,kml),dpeps)
              end do
            end do
          end do
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              zbot(i,j,1) = zbot(i,j,1)*vbath(i,j)/zbot(i,j,kk)
              ztop(i,j,1) = 0.
              ind1(i,j) = 1
            end do
          end do
          do kl = 2,kk
            do l = 1,isv(j)
              do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                zbot(i,j,kl) = zbot(i,j,kl)*vbath(i,j)/zbot(i,j,kk)
                ztop(i,j,kl) = zbot(i,j,kl-1)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else
        write (lp,'(a)') 'cannot identify grid'
        flush(lp)
        call xchalt('(diazlv)')
        stop '(diazlv)'
      end if
    end if

    ! Compute interpolation weights
    if (gridid == 'p') then
      !$omp parallel do private(l,i,d)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            ind2(i,j) = 0
            if (dp(i,j,km) > dpeps) then
              do d = ind1(i,j),ddm
                if (depthslev_bnds(2,d) <= ztop(i,j,k)) then
                  ind1(i,j) = d+1
                  cycle
                else if (depthslev_bnds(1,d) >= zbot(i,j,k)) then
                  exit
                end if
                ind2(i,j) = d
                weights(i,j,d) = (min(zbot(i,j,k), &
                     depthslev_bnds(2,d))-max(ztop(i,j,k), &
                     depthslev_bnds(1,d)))/dlevp(i,j,d)
              end do
            end if
          end do
        end do
      end do
      !$omp end parallel do

    else if (gridid == 'u') then
      !$omp parallel do private(l,i,d,r)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            ind2(i,j) = 0
            if (abs(uflx(i,j,kn)) > flxeps) then
              do d = ind1(i,j),ddm
                if (depthslev_bnds(2,d) < ztop(i,j,k)) then
                  ind1(i,j) = d+1
                  cycle
                else if (depthslev_bnds(1,d) > zbot(i,j,k)) then
                  exit
                end if
                ind2(i,j) = d
                r = (min(zbot(i,j,k),depthslev_bnds(2,d))- &
                     max(ztop(i,j,k),depthslev_bnds(1,d)))
                weights(i,j,d) = r/dlevu(i,j,d)
                weightsflx(i,j,d) = r/(zbot(i,j,k)-ztop(i,j,k))
              end do
            end if
          end do
        end do
      end do
      !$omp end parallel do

    else if (gridid == 'v') then
      !$omp parallel do private(l,i,d,r)
      do j = 1,jj+1
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            ind2(i,j) = 0
            if (abs(vflx(i,j,kn)) > flxeps) then
              do d = ind1(i,j),ddm
                if (depthslev_bnds(2,d) < ztop(i,j,k)) then
                  ind1(i,j) = d+1
                  cycle
                else if (depthslev_bnds(1,d) > zbot(i,j,k)) then
                  exit
                end if
                ind2(i,j) = d
                r = (min(zbot(i,j,k),depthslev_bnds(2,d))- &
                     max(ztop(i,j,k),depthslev_bnds(1,d)))
                weights(i,j,d) = r/dlevv(i,j,d)
                weightsflx(i,j,d) = r/(zbot(i,j,k)-ztop(i,j,k))
              end do
            end if
          end do
        end do
      end do
      !$omp end parallel do

    else
      write (lp,'(a)') 'cannot identify grid'
      flush(lp)
      call xchalt('(diazlv)')
      stop '(diazlv)'
    end if

  end subroutine diazlv


  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !---------------------------------------------------------------


  subroutine inih2d(pos,gridid,inival)
  !---------------------------------------------------------------
  ! Description: initialise 2d diagnostic field
  !---------------------------------------------------------------

    ! Arguments
    integer,   intent(in) :: pos    ! position in common buffer
    real,      intent(in) :: inival ! value used for initalisation
    character, intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j

    ! Check whether field should be initialised
    if (pos == 0) return

    if (gridid(1:1) == 'u') then
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          phyh2d(i,j,pos) = inival*iu(i,j)
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'v') then
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          phyh2d(i,j,pos) = inival*iv(i,j)
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'p') then
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          phyh2d(i,j,pos) = inival*ip(i,j)
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          phyh2d(i,j,pos) = inival
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine inih2d



  subroutine inilyr(pos,gridid,inival)
  !---------------------------------------------------------------
  ! Description: initialise layer diagnostic field
  !---------------------------------------------------------------

    ! Arguments
    integer,   intent(in) :: pos    ! position in common buffer
    character, intent(in) :: gridid ! grid identifier ('p','u' or 'v')
    real,      intent(in) :: inival ! value used for initalisation

    ! LOCAL VARIABLES
    integer :: i,j,k

    ! Check whether field should be initialised
    if (pos == 0) return

    if (gridid(1:1) == 'u') then
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,kk
          do i = 1-nbdy,ii+nbdy
            phylyr(i,j,k,pos) = inival*iu(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'v') then
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,kk
          do i = 1-nbdy,ii+nbdy
            phylyr(i,j,k,pos) = inival*iv(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'p') then
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,kk
          do i = 1-nbdy,ii+nbdy
            phylyr(i,j,k,pos) = inival*ip(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,kk
          do i = 1-nbdy,ii+nbdy
            phylyr(i,j,k,pos) = inival
          end do
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine inilyr



  subroutine inilvl(pos,gridid,inival)
  !---------------------------------------------------------------
  ! Description: initialise level diagnostic field
  !
  ! Arguments:
  !   int  pos      (in)     : position in common buffer
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !   real inival   (in)     : value used for initalisation
  !---------------------------------------------------------------

    integer :: pos
    real :: inival
    character :: gridid

    integer :: i,j,k

    ! Check whether field should be initialised
    if (pos == 0) return

    if (gridid(1:1) == 'u') then
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,ddm
          do i = 1-nbdy,ii+nbdy
            phylvl(i,j,k,pos) = inival*iu(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'v') then
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,ddm
          do i = 1-nbdy,ii+nbdy
            phylvl(i,j,k,pos) = inival*iv(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'p') then
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,ddm
          do i = 1-nbdy,ii+nbdy
            phylvl(i,j,k,pos) = inival*ip(i,j)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(k,i)
      do j = 1-nbdy,jj+nbdy
        do k = 1,ddm
          do i = 1-nbdy,ii+nbdy
            phylvl(i,j,k,pos) = inival
          end do
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine inilvl



  subroutine acch2d(pos,fld,wghts,wghtsflg,gridid)
  !---------------------------------------------------------------
  ! Description: accumulate 2d fields
  !
  ! Arguments:
  !   int  pos      (in)     : position in 2d buffer
  !   real fld      (in)     : input data used for accumulation
  !   real wghts    (in)     : weights used for accumulation
  !   int  wghtsflg (in)     : weights flag (0=no weighting)
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !---------------------------------------------------------------

    integer :: pos(nphymax),wghtsflg
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld,wghts
    character :: gridid

    integer :: i,j,l,o

    ! Check whether field should be accumulated
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        if (wghtsflg == 0) then
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isu(j)
              do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
                phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isu(j)
              do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
                phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)* &
                     wghts(i,j)
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else if (gridid == 'v') then
        if (wghtsflg == 0) then
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isv(j)
              do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isv(j)
              do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)* &
                     wghts(i,j)
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else if (gridid == 'p') then
        if (wghtsflg == 0) then
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(l,i)
          do j = 1,jj
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)* &
                     wghts(i,j)
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(acch2d)')
        stop '(acch2d)'
      end if

    end do

  end subroutine acch2d



  subroutine maxh2d(pos,fld,gridid)
  !---------------------------------------------------------------
  ! Description: store maximum of 2d fields
  !
  ! Arguments:
  !   int  pos      (in)     : position in 2d buffer
  !   real fld      (in)     : input data used for finding maximum
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !---------------------------------------------------------------

    integer :: pos(nphymax)
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
    character :: gridid

    integer :: i,j,l,o

    ! Check whether maximum of field should be stored
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              phyh2d(i,j,pos(o)) = max(phyh2d(i,j,pos(o)),fld(i,j))
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'v') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              phyh2d(i,j,pos(o)) = max(phyh2d(i,j,pos(o)),fld(i,j))
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'p') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              phyh2d(i,j,pos(o)) = max(phyh2d(i,j,pos(o)),fld(i,j))
            end do
          end do
        end do
        !$omp end parallel do
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(maxh2d)')
        stop '(maxh2d)'
      end if

    end do

  end subroutine maxh2d



  subroutine minh2d(pos,fld,gridid)
  !---------------------------------------------------------------
  ! Description: store minimum of 2d fields
  !
  ! Arguments:
  !   int  pos      (in)     : position in 2d buffer
  !   real fld      (in)     : input data used for finding minimum
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !---------------------------------------------------------------

    integer :: pos(nphymax)
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
    character :: gridid

    integer :: i,j,l,o

    ! Check whether minimum of field should be stored
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              phyh2d(i,j,pos(o)) = min(phyh2d(i,j,pos(o)),fld(i,j))
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'v') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              phyh2d(i,j,pos(o)) = min(phyh2d(i,j,pos(o)),fld(i,j))
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'p') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              phyh2d(i,j,pos(o)) = min(phyh2d(i,j,pos(o)),fld(i,j))
            end do
          end do
        end do
        !$omp end parallel do
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(minh2d)')
        stop '(minh2d)'
      end if

    end do

  end subroutine minh2d



  subroutine sqh2d(pos,fld,gridid)
  !---------------------------------------------------------------
  ! Description: accumulate square of 2d fields
  !
  ! Arguments:
  !   int  pos      (in)     : position in 2d buffer
  !   real fld      (in)     : input data used for accumulation
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !---------------------------------------------------------------

    integer :: pos(nphymax)
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
    character :: gridid

    integer :: i,j,l,o

    ! Check whether field should be accumulated
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)*fld(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'v') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)*fld(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'p') then
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              phyh2d(i,j,pos(o)) = phyh2d(i,j,pos(o))+fld(i,j)*fld(i,j)
            end do
          end do
        end do
        !$omp end parallel do
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(minh2d)')
        stop '(minh2d)'
      end if

    end do

  end subroutine sqh2d



  subroutine acclyr(pos,fld,wghts,wghtsflg,gridid)
  !---------------------------------------------------------------
  ! Description: accumulate layer fields
  !
  ! Arguments:
  !   int  pos      (in)     : position in 3d layer buffer
  !   real fld      (in)     : input data used for accumulation
  !   real wghts    (in)     : weights used for accumulation
  !   int  wghtsflg (in)     : weights flag (0=no weighting)
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !---------------------------------------------------------------

    integer :: pos(nphymax),wghtsflg
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: fld,wghts
    character :: gridid

    integer :: i,j,k,l,o

    ! Check whether field should be accumulated
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        if (wghtsflg == 0) then
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isu(j)
                do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o))+fld(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isu(j)
                do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +fld(i,j,k)*wghts(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else if (gridid == 'v') then
        if (wghtsflg == 0) then
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isv(j)
                do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o))+fld(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isv(j)
                do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +fld(i,j,k)*wghts(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else if (gridid == 'p') then
        if (wghtsflg == 0) then
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isp(j)
                do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o))+fld(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isp(j)
                do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +fld(i,j,k)*wghts(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(acclyr)')
        stop '(acclyr)'
      end if

    end do

  end subroutine acclyr



  subroutine accily(pos,fld,wghts,wghtsflg,gridid)
  !---------------------------------------------------------------
  ! Description: accumulate interface fields after interpolation to
  !              layers
  !
  ! Arguments:
  !   int  pos      (in)     : position in 3d layer buffer
  !   real fld      (in)     : input data used for accumulation
  !   real wghts    (in)     : weights used for accumulation
  !   int  wghtsflg (in)     : weights flag (0=no weighting)
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !---------------------------------------------------------------

    integer :: pos(nphymax),wghtsflg
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: fld
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: wghts
    character :: gridid

    integer :: i,j,k,l,o

    ! Check whether field should be accumulated
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        if (wghtsflg == 0) then
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isu(j)
                do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +.5*(fld(i,j,k)+fld(i,j,k+1))
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isu(j)
                do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +.5*(fld(i,j,k)+fld(i,j,k+1)) &
                       *wghts(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else if (gridid == 'v') then
        if (wghtsflg == 0) then
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isv(j)
                do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +.5*(fld(i,j,k)+fld(i,j,k+1))
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isv(j)
                do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +.5*(fld(i,j,k)+fld(i,j,k+1)) &
                       *wghts(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else if (gridid == 'p') then
        if (wghtsflg == 0) then
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isp(j)
                do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +.5*(fld(i,j,k)+fld(i,j,k+1))
                end do
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(k,l,i)
          do j = 1,jj
            do k = 1,kk
              do l = 1,isp(j)
                do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                  phylyr(i,j,k,pos(o)) = phylyr(i,j,k,pos(o)) &
                       +.5*(fld(i,j,k)+fld(i,j,k+1)) &
                       *wghts(i,j,k)
                end do
              end do
            end do
          end do
          !$omp end parallel do
        end if
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(accily)')
        stop '(accily)'
      end if

    end do

  end subroutine accily



  subroutine acclvl(pos,fld,gridid,k,ind1,ind2,wghts)
  !---------------------------------------------------------------
  ! Description: accumulate layer fields mapped to levels
  !
  ! Arguments:
  !   int  pos      (in)     : position in buffer
  !   real fld      (in)     : input data used for accumulation
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !   int  k        (in)     : layer index of fld
  !   int  ind1     (in)     : index field for first accumulated level
  !   int  ind2     (in)     : index field for last accumulated level
  !   real wghts    (in)     : weights used for accumulation
  !---------------------------------------------------------------

    integer :: pos(nphymax),k
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: fld
    character :: gridid

    integer :: d,i,j,l,o

    ! Check whether field should be accumulated
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        !$omp parallel do private(l,i,d)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              do d = ind1(i,j),ind2(i,j)
                phylvl(i,j,d,pos(o)) = phylvl(i,j,d,pos(o)) &
                                      +fld(i,j,k)*wghts(i,j,d)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'v') then
        !$omp parallel do private(l,i,d)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              do d = ind1(i,j),ind2(i,j)
                phylvl(i,j,d,pos(o)) = phylvl(i,j,d,pos(o)) &
                                      +fld(i,j,k)*wghts(i,j,d)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid(1:1) == 'p') then
        !$omp parallel do private(l,i,d)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              do d = ind1(i,j),ind2(i,j)
                phylvl(i,j,d,pos(o)) = phylvl(i,j,d,pos(o)) &
                                      +fld(i,j,k)*wghts(i,j,d)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(acclvl)')
        stop '(acclvl)'
      end if
    end do

  end subroutine acclvl



  subroutine accilv(pos,fld,gridid,k,ind1,ind2,wghts)
  !---------------------------------------------------------------
  ! Description: accumulate interface fields mapped to levels
  !
  ! Arguments:
  !   int  pos      (in)     : position in buffer
  !   real fld      (in)     : input data used for accumulation
  !   char gridid   (in)     : grid identifier ('p','u' or 'v')
  !   int  k        (in)     : layer index of fld
  !   int  ind1     (in)     : index field for first accumulated level
  !   int  ind2     (in)     : index field for last accumulated level
  !   real wghts    (in)     : weights used for accumulation
  !---------------------------------------------------------------

    integer :: pos(nphymax),k
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: fld
    character :: gridid

    integer :: d,i,j,l,o

    ! Check whether field should be accumulated
    do o = 1,nphy
      if (pos(o) == 0) cycle

      if (gridid == 'u') then
        !$omp parallel do private(l,i,d)
        do j = 1,jj
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              do d = ind1(i,j),ind2(i,j)
                phylvl(i,j,d,pos(o)) = phylvl(i,j,d,pos(o)) &
                     +.5*(fld(i,j,k)+fld(i,j,k+1)) &
                     *wghts(i,j,d)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid == 'v') then
        !$omp parallel do private(l,i,d)
        do j = 1,jj+1
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              do d = ind1(i,j),ind2(i,j)
                phylvl(i,j,d,pos(o)) = phylvl(i,j,d,pos(o)) &
                     +.5*(fld(i,j,k)+fld(i,j,k+1)) &
                     *wghts(i,j,d)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else if (gridid(1:1) == 'p') then
        !$omp parallel do private(l,i,d)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              do d = ind1(i,j),ind2(i,j)
                phylvl(i,j,d,pos(o)) = phylvl(i,j,d,pos(o)) &
                     +.5*(fld(i,j,k)+fld(i,j,k+1)) &
                     *wghts(i,j,d)
              end do
            end do
          end do
        end do
        !$omp end parallel do
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(accilv)')
        stop '(accilv)'
      end if
    end do

  end subroutine accilv



  subroutine inifld(iogrp)

    ! Arguments
    integer, intent(in) :: iogrp

    ! initialisation of 2d fields
    call inih2d(ACC_UB(iogrp),'u',0.)
    call inih2d(ACC_UBFLXS(iogrp),'u',0.)
    call inih2d(ACC_ZTX(iogrp),'u',0.)
    call inih2d(ACC_TAUX(iogrp),'u',0.)
    call inih2d(ACC_UICE(iogrp),'u',0.)
    call inih2d(ACC_IVOLU(iogrp),'u',0.)
    call inih2d(ACC_USTOKES(iogrp),'u',0.)

    call inih2d(ACC_VB(iogrp),'v',0.)
    call inih2d(ACC_VBFLXS(iogrp),'v',0.)
    call inih2d(ACC_MTY(iogrp),'v',0.)
    call inih2d(ACC_TAUY(iogrp),'v',0.)
    call inih2d(ACC_VICE(iogrp),'v',0.)
    call inih2d(ACC_IVOLV(iogrp),'v',0.)
    call inih2d(ACC_VSTOKES(iogrp),'v',0.)

    call inih2d(ACC_PSRF(iogrp),'p',0.)
    call inih2d(ACC_PBOT(iogrp),'p',0.)
    call inih2d(ACC_SEALV(iogrp),'p',0.)
    call inih2d(ACC_SLVSQ(iogrp),'p',0.)
    call inih2d(ACC_SSS(iogrp),'p',0.)
    call inih2d(ACC_SSSSQ(iogrp),'p',0.)
    call inih2d(ACC_SBOT(iogrp),'p',0.)
    call inih2d(ACC_SST(iogrp),'p',0.)
    call inih2d(ACC_SSTSQ(iogrp),'p',0.)
    call inih2d(ACC_TBOT(iogrp),'p',0.)
    call inih2d(ACC_SIGMX(iogrp),'p',0.)
    call inih2d(ACC_MLD(iogrp),'p',0.)
    call inih2d(ACC_MAXMLD(iogrp),'p',-spval)
    call inih2d(ACC_MLTS(iogrp),'p',0.)
    call inih2d(ACC_MLTSMN(iogrp),'p', spval)
    call inih2d(ACC_MLTSMX(iogrp),'p',-spval)
    call inih2d(ACC_MLTSSQ(iogrp),'p',0.)
    call inih2d(ACC_T20D(iogrp),'p',0.)
    call inih2d(ACC_ALB(iogrp),'p',0.)
    call inih2d(ACC_SWA(iogrp),'p',0.)
    call inih2d(ACC_NSF(iogrp),'p',0.)
    call inih2d(ACC_HMAT(iogrp),'p',0.)
    call inih2d(ACC_DFL(iogrp),'p',0.)
    call inih2d(ACC_LIP(iogrp),'p',0.)
    call inih2d(ACC_SOP(iogrp),'p',0.)
    call inih2d(ACC_EVA(iogrp),'p',0.)
    call inih2d(ACC_RNFFLX(iogrp),'p',0.)
    call inih2d(ACC_RFIFLX(iogrp),'p',0.)
    call inih2d(ACC_SFL(iogrp),'p',0.)
    call inih2d(ACC_BRNFLX(iogrp),'p',0.)
    call inih2d(ACC_BRNPD(iogrp),'p',0.)
    call inih2d(ACC_SURFLX(iogrp),'p',0.)
    call inih2d(ACC_SURRLX(iogrp),'p',0.)
    call inih2d(ACC_SALFLX(iogrp),'p',0.)
    call inih2d(ACC_SALRLX(iogrp),'p',0.)
    call inih2d(ACC_ABSWND(iogrp),'p',0.)
    call inih2d(ACC_USTAR(iogrp),'p',0.)
    call inih2d(ACC_USTAR3(iogrp),'p',0.)
    call inih2d(ACC_IDKEDT(iogrp),'p',0.)
    call inih2d(ACC_MTKEUS(iogrp),'p',0.)
    call inih2d(ACC_MTKENI(iogrp),'p',0.)
    call inih2d(ACC_MTKEBF(iogrp),'p',0.)
    call inih2d(ACC_MTKERS(iogrp),'p',0.)
    call inih2d(ACC_MTKEPE(iogrp),'p',0.)
    call inih2d(ACC_MTKEKE(iogrp),'p',0.)
    call inih2d(ACC_LAMULT(iogrp),'p',0.)
    call inih2d(ACC_LASL(iogrp),'p',0.)
    call inih2d(ACC_FMLTFZ(iogrp),'p',0.)
    call inih2d(ACC_HMLTFZ(iogrp),'p',0.)
    call inih2d(ACC_HICE(iogrp),'p',0.)
    call inih2d(ACC_HSNW(iogrp),'p',0.)
    call inih2d(ACC_FICE(iogrp),'p',0.)
    call inih2d(ACC_TSRF(iogrp),'p',0.)
    call inih2d(ACC_TICE(iogrp),'p',0.)

    ! initialisation of 3d layer fields
    call inilyr(ACC_UVEL(iogrp),'u',0.)
    call inilyr(ACC_DPU(iogrp),'u',0.)
    call inilyr(ACC_UFLX(iogrp),'u',0.)
    call inilyr(ACC_UTFLX(iogrp),'u',0.)
    call inilyr(ACC_USFLX(iogrp),'u',0.)
    call inilyr(ACC_UMFLTD(iogrp),'u',0.)
    call inilyr(ACC_UMFLSM(iogrp),'u',0.)
    call inilyr(ACC_UTFLTD(iogrp),'u',0.)
    call inilyr(ACC_UTFLSM(iogrp),'u',0.)
    call inilyr(ACC_UTFLLD(iogrp),'u',0.)
    call inilyr(ACC_USFLTD(iogrp),'u',0.)
    call inilyr(ACC_USFLSM(iogrp),'u',0.)
    call inilyr(ACC_USFLLD(iogrp),'u',0.)

    call inilyr(ACC_VVEL(iogrp),'v',0.)
    call inilyr(ACC_DPV(iogrp),'v',0.)
    call inilyr(ACC_VFLX(iogrp),'v',0.)
    call inilyr(ACC_VTFLX(iogrp),'v',0.)
    call inilyr(ACC_VSFLX(iogrp),'v',0.)
    call inilyr(ACC_VMFLTD(iogrp),'v',0.)
    call inilyr(ACC_VMFLSM(iogrp),'v',0.)
    call inilyr(ACC_VTFLTD(iogrp),'v',0.)
    call inilyr(ACC_VTFLSM(iogrp),'v',0.)
    call inilyr(ACC_VTFLLD(iogrp),'v',0.)
    call inilyr(ACC_VSFLTD(iogrp),'v',0.)
    call inilyr(ACC_VSFLSM(iogrp),'v',0.)
    call inilyr(ACC_VSFLLD(iogrp),'v',0.)

    call inilyr(ACC_SALN(iogrp),'p',0.)
    call inilyr(ACC_TEMP(iogrp),'p',0.)
    call inilyr(ACC_DP(iogrp),'p',0.)
    call inilyr(ACC_DZ(iogrp),'p',0.)
    call inilyr(ACC_BFSQ(iogrp),'p',0.)
    call inilyr(ACC_DIFDIA(iogrp),'p',0.)
    call inilyr(ACC_DIFVMO(iogrp),'p',0.)
    call inilyr(ACC_DIFVHO(iogrp),'p',0.)
    call inilyr(ACC_DIFVSO(iogrp),'p',0.)
    call inilyr(ACC_DIFINT(iogrp),'p',0.)
    call inilyr(ACC_DIFISO(iogrp),'p',0.)
    call inilyr(ACC_WFLX(iogrp),'p',0.)
    call inilyr(ACC_WFLX2(iogrp),'p',0.)
    call inilyr(ACC_AVDSG(iogrp),'p',0.)
    call inilyr(ACC_DPVOR(iogrp),'p',0.)
    if (use_TRC .and. use_TKE) then
      call inilyr(ACC_TKE(iogrp),'p',0.)
      call inilyr(ACC_GLS_PSI(iogrp),'p',0.)
    endif

    ! initialsation of 3d level fields
    call inilvl(ACC_UVELLVL(iogrp),'u',0.)
    call inilvl(ACC_UFLXLVL(iogrp),'u',0.)
    call inilvl(ACC_UTFLXLVL(iogrp),'u',0.)
    call inilvl(ACC_USFLXLVL(iogrp),'u',0.)
    call inilvl(ACC_UMFLTDLVL(iogrp),'u',0.)
    call inilvl(ACC_UMFLSMLVL(iogrp),'u',0.)
    call inilvl(ACC_UTFLTDLVL(iogrp),'u',0.)
    call inilvl(ACC_UTFLSMLVL(iogrp),'u',0.)
    call inilvl(ACC_UTFLLDLVL(iogrp),'u',0.)
    call inilvl(ACC_USFLTDLVL(iogrp),'u',0.)
    call inilvl(ACC_USFLSMLVL(iogrp),'u',0.)
    call inilvl(ACC_USFLLDLVL(iogrp),'u',0.)

    call inilvl(ACC_VVELLVL(iogrp),'v',0.)
    call inilvl(ACC_VFLXLVL(iogrp),'v',0.)
    call inilvl(ACC_VTFLXLVL(iogrp),'v',0.)
    call inilvl(ACC_VSFLXLVL(iogrp),'v',0.)
    call inilvl(ACC_VMFLTDLVL(iogrp),'v',0.)
    call inilvl(ACC_VMFLSMLVL(iogrp),'v',0.)
    call inilvl(ACC_VTFLTDLVL(iogrp),'v',0.)
    call inilvl(ACC_VTFLSMLVL(iogrp),'v',0.)
    call inilvl(ACC_VTFLLDLVL(iogrp),'v',0.)
    call inilvl(ACC_VSFLTDLVL(iogrp),'v',0.)
    call inilvl(ACC_VSFLSMLVL(iogrp),'v',0.)
    call inilvl(ACC_VSFLLDLVL(iogrp),'v',0.)

    call inilvl(ACC_BFSQLVL(iogrp),'p',0.)
    call inilvl(ACC_DIFDIALVL(iogrp),'p',0.)
    call inilvl(ACC_DIFVMOLVL(iogrp),'p',0.)
    call inilvl(ACC_DIFVHOLVL(iogrp),'p',0.)
    call inilvl(ACC_DIFVSOLVL(iogrp),'p',0.)
    call inilvl(ACC_DIFINTLVL(iogrp),'p',0.)
    call inilvl(ACC_DIFISOLVL(iogrp),'p',0.)
    call inilvl(ACC_DZLVL(iogrp),'p',0.)
    call inilvl(ACC_SALNLVL(iogrp),'p',0.)
    call inilvl(ACC_TEMPLVL(iogrp),'p',0.)
    call inilvl(ACC_WFLXLVL(iogrp),'p',0.)
    call inilvl(ACC_WFLX2LVL(iogrp),'p',0.)
    call inilvl(ACC_PVLVL(iogrp),'p',0.)
    if (use_TRC .and. use_TKE) then
      call inilvl(ACC_TKELVL(iogrp),'p',0.)
      call inilvl(ACC_GLS_PSILVL(iogrp),'p',0.)
    endif

  end subroutine inifld



  subroutine finh2d(posacc,poswgt,gridid)
  !---------------------------------------------------------------
  ! Description: finalise accumulation of weighted 2d fields
  !---------------------------------------------------------------

    ! Arguments
    integer,   intent(in) :: posacc ! position of accumulated field in buffer
    integer,   intent(in) :: poswgt ! position of accumulated weights
    character, intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j,l
    real, parameter :: epsil = 1e-11

    ! Check whether field should be initialised
    if (posacc == 0) return

    if (gridid == 'u') then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            phyh2d(i,j,posacc) = phyh2d(i,j,posacc)/ &
                 max(epsil,phyh2d(i,j,poswgt))
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid(1:1) == 'v') then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            phyh2d(i,j,posacc) = phyh2d(i,j,posacc)/ &
                 max(epsil,phyh2d(i,j,poswgt))
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'p') then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            phyh2d(i,j,posacc) = phyh2d(i,j,posacc)/ &
                 max(epsil,phyh2d(i,j,poswgt))
          end do
        end do
      end do
      !$omp end parallel do
    else
      write (lp,*) 'cannot identify grid '//gridid//'!'
      call xchalt('(finh2d)')
      stop '(finh2d)'
    end if

  end subroutine finh2d



  subroutine finlyr(posacc,poswgt,gridid)
  !---------------------------------------------------------------
  ! Description: finalise accumulation of weighted 3d layer fields
  !---------------------------------------------------------------

    ! Arguments
    integer,   intent(in) :: posacc ! position of accumulated field in buffer
    integer,   intent(in) :: poswgt ! position of accumulated weights
    character, intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j,k,l
    real, parameter :: epsil = 1e-11

    ! Check whether field should be initialised
    if (posacc == 0) return

    if (gridid == 'u') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,kk
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (phylyr(i,j,k,poswgt) > epsil) then
                phylyr(i,j,k,posacc) = phylyr(i,j,k,posacc)/phylyr(i,j,k,poswgt)
              else
                phylyr(i,j,k,posacc) = nf90_fill_double
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'v') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,kk
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (phylyr(i,j,k,poswgt) > epsil) then
                phylyr(i,j,k,posacc) = phylyr(i,j,k,posacc)/phylyr(i,j,k,poswgt)
              else
                phylyr(i,j,k,posacc) = nf90_fill_double
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'p') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (phylyr(i,j,k,poswgt) > epsil) then
                phylyr(i,j,k,posacc) = phylyr(i,j,k,posacc)/phylyr(i,j,k,poswgt)
              else
                phylyr(i,j,k,posacc) = nf90_fill_double
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else
      write (lp,*) 'cannot identify grid '//gridid//'!'
      call xchalt('(finlyr)')
      stop '(finlyr)'
    end if

  end subroutine finlyr



  subroutine wrth2d(pos,frmt,sfac,offs,cmpflg,msk,gridid, &
       vnm,vlngnm,vstdnm,vunits)
  !---------------------------------------------------------------
  ! Description: writes diagnostic 2d field to file
  !---------------------------------------------------------------

    ! Arguments
    real,                  intent(in) :: sfac   ! user defined scale factor to be applied
    real,                  intent(in) :: offs   ! user defined offset to be added
    integer,               intent(in) :: frmt   ! format/precision of output
    ! 0=field is not written
    ! 2=field is written as int2 with scale
    ! factor and offset
    ! 4=field is written as real4
    ! 8=field is written as real8
    integer,               intent(in) :: cmpflg ! compression flag; only wet points are
    ! written if flag is set to 1
    integer,               intent(in) :: pos    ! variable position in common buffer
    character(len=*),      intent(in) :: gridid ! grid identifier ('p','u' or 'v')
    character(len=*),      intent(in) :: vnm    ! variable name used in nc-file
    character(len=*),      intent(in) :: vlngnm ! variable long name (skipped if ' ')
    character(len=*),      intent(in) :: vstdnm ! variable standard name (skipped if ' ')
    character(len=*),      intent(in) :: vunits ! variable units (skipped if ' ')
    integer, dimension(*), intent(in) :: msk    ! ocean mask

    ! Local variables
    character(len = 100) :: dims

    ! Check whether field should be written
    if (frmt == 0) return

    ! Create dimension string
    if (cmpflg == 1) then
      dims = gridid(1:1)//'comp time'
    else
      dims = 'x y time'
    end if

    ! Check output format
    if (frmt == 2) then
      if (cmpflg == 1) then
        call nccopa(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,sfac, &
             offs)
      else
        call ncpack(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk, &
             1,sfac,offs)
      end if
    else if (frmt == 4) then
      if (cmpflg == 1) then
        call nccomp(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,sfac, &
             offs,4)
      else
        call ncwrtr(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk, &
             1,sfac,offs,4)
      end if
    else if (frmt == 8) then
      if (cmpflg == 1) then
        call nccomp(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,sfac, &
             offs,8)
      else
        call ncwrtr(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk, &
             1,sfac,offs,8)
      end if
    else
      write (lp,*) 'unknown output format!'
      call xchalt('(wrth2d)')
      stop '(wrth2d)'
    end if

    ! Define attributes
    !      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
    !      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
    !      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
    !      call ncattr('coordinates',
    !     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
    !      call ncattr('cell_measures','area: '//gridid(1:1)//'area')

  end subroutine wrth2d



  subroutine wrtlyr(pos,frmt,sfac,offs,cmpflg,msk,gridid, &
       vnm,vlngnm,vstdnm,vunits)
  !---------------------------------------------------------------
  ! Description: writes diagnostic layer field to file
  !---------------------------------------------------------------

    ! Arguments
    real,                  intent(in) :: sfac   ! user defined scale factor to be applied
    real,                  intent(in) :: offs   ! user defined offset to be added
    integer,               intent(in) :: frmt   ! format/precision of output
    ! 0=field is not written
    ! 2=field is written as int2 with scale
    ! factor and offset
    ! 4=field is written as real4
    ! 8=field is written as real8
    integer,               intent(in) :: cmpflg ! compression flag; only wet points are
    ! written if flag is set to 1
    integer,               intent(in) :: pos    ! variable position in common buffer
    character(len=*),      intent(in) :: gridid ! grid identifier ('p','u' or 'v')
    character(len=*),      intent(in) :: vnm    ! variable name used in nc-file
    character(len=*),      intent(in) :: vlngnm ! variable long name (skipped if ' ')
    character(len=*),      intent(in) :: vstdnm ! variable standard name (skipped if ' ')
    character(len=*),      intent(in) :: vunits ! variable units (skipped if ' ')
    integer, dimension(*), intent(in) :: msk    ! ocean mask

    ! Local variables
    character(len = 100) :: dims

    ! Check whether field should be written
    if (frmt == 0) return

    ! Create dimension string
    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      if (cmpflg == 1) then
        dims = gridid(1:1)//'comp layer time'
      else
        dims = 'x y layer time'
      end if
    else
      if (cmpflg == 1) then
        dims = gridid(1:1)//'comp sigma time'
      else
        dims = 'x y sigma time'
      end if
    end if

    ! Check output format
    if (frmt == 2) then
      if (cmpflg == 1) then
        call nccopa(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,sfac, &
             offs)
      else
        call ncpack(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk, &
             2,sfac,offs)
      end if
    else if (frmt == 4) then
      if (cmpflg == 1) then
        call nccomp(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,sfac, &
             offs,4)
      else
        call ncwrtr(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk, &
             1,sfac,offs,4)
      end if
    else if (frmt == 8) then
      if (cmpflg == 1) then
        call nccomp(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,sfac, &
             offs,8)
      else
        call ncwrtr(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk, &
             1,sfac,offs,8)
      end if
    else
      write (lp,*) 'unknown output format!'
      call xchalt('(wrtlyr)')
      stop '(wrtlyr)'
    end if

    ! Define attributes
    !      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
    !      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
    !      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
    !      call ncattr('coordinates',
    !     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
    !      call ncattr('cell_measures','area: '//gridid(1:1)//'area')

  end subroutine wrtlyr



  subroutine wrtlvl(pos,frmt,sfac,offs,cmpflg,msk,gridid, &
       vnm,vlngnm,vstdnm,vunits)
  !---------------------------------------------------------------
  ! Description: writes diagnostic level field to file
  ! ------------------------------------------------------------------

    ! Arguments
    real,                  intent(in) :: sfac   ! user defined scale factor to be applied
    real,                  intent(in) :: offs   ! user defined offset to be added
    integer,               intent(in) :: frmt   ! format/precision of output
    ! 0=field is not written
    ! 2=field is written as int2 with scale
    ! factor and offset
    ! 4=field is written as real4
    ! 8=field is written as real8
    integer,               intent(in) :: cmpflg ! compression flag; only wet points are
    ! written if flag is set to 1
    integer,               intent(in) :: pos    ! variable position in common buffer
    character(len=*),      intent(in) :: gridid ! grid identifier ('p','u' or 'v')
    character(len=*),      intent(in) :: vnm    ! variable name used in nc-file
    character(len=*),      intent(in) :: vlngnm ! variable long name (skipped if ' ')
    character(len=*),      intent(in) :: vstdnm ! variable standard name (skipped if ' ')
    character(len=*),      intent(in) :: vunits ! variable units (skipped if ' ')
    integer, dimension(*), intent(in) :: msk    ! ocean mask

    ! Local variables
    character(len = 100) :: dims

    ! Check whether field should be written
    if (frmt == 0) return

    ! Create dimension string
    if (cmpflg == 1) then
      dims = gridid//'comp depth time'
    else
      dims = 'x y depth time'
    end if

    ! Check output format
    if (frmt == 2) then
      if (cmpflg == 1) then
        call nccopa(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,sfac, &
             offs)
      else
        call ncpack(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk, &
             2,sfac,offs)
      end if
    else if (frmt == 4) then
      if (cmpflg == 1) then
        call nccomp(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,sfac, &
             offs,4)
      else
        call ncwrtr(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk, &
             1,sfac,offs,4)
      end if
    else if (frmt == 8) then
      if (cmpflg == 1) then
        call nccomp(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,sfac, &
             offs,8)
      else
        call ncwrtr(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk, &
             1,sfac,offs,8)
      end if
    else
      write (lp,*) 'unknown output format!'
      call xchalt('(wrtlvl)')
      stop '(wrtlvl)'
    end if

    ! Define attributes
    !      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
    !      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
    !      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
    !      call ncattr('coordinates',
    !     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
    !      call ncattr('cell_measures','area: '//gridid(1:1)//'area')

  end subroutine wrtlvl



  subroutine logh2d(pos,gridid,sfac,offs)
  ! ------------------------------------------------------------------
  ! Description: replace 2d field with log10(field)
  ! ------------------------------------------------------------------

    ! Arguments
    real,      intent(in) :: sfac   ! scale factor to be applied before log10
    real,      intent(in) :: offs   ! offset to be added before log10
    integer,   intent(in) :: pos    ! field position in layer buffer
    character, intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j,l
    real :: epsil = 1e-11

    ! Check whether field should be processed
    if (pos == 0) return

    if (gridid == 'u') then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
            if (phyh2d(i,j,pos) < epsil) then
              phyh2d(i,j,pos) = 0.
            else
              phyh2d(i,j,pos) = log10(phyh2d(i,j,pos)*sfac+offs)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'v') then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (phyh2d(i,j,pos) < epsil) then
              phyh2d(i,j,pos) = 0.
            else
              phyh2d(i,j,pos) = log10(phyh2d(i,j,pos)*sfac+offs)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'p') then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (phyh2d(i,j,pos) < epsil) then
              phyh2d(i,j,pos) = 0.
            else
              phyh2d(i,j,pos) = log10(phyh2d(i,j,pos)*sfac+offs)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    else
      write (lp,*) 'cannot identify grid '//gridid//'!'
      call xchalt('(logh2d)')
      stop '(logh2d)'
    end if

  end subroutine logh2d


  subroutine loglyr(pos,gridid,sfac,offs)
  ! ------------------------------------------------------------------
  ! Description: replace 3d layer field with log10(field)
  ! ------------------------------------------------------------------

    ! Arguments
    real,      intent(in) :: sfac   ! scale factor to be applied before log10
    real,      intent(in) :: offs   ! offset to be added before log10
    integer,   intent(in) :: pos    ! field position in layer buffer
    character, intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j,k,l
    real :: epsil = 1e-11

    ! Check whether field should be processed
    if (pos == 0) return

    if (gridid == 'u') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,kk
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (phylyr(i,j,k,pos) < epsil) then
                phylyr(i,j,k,pos) = 0.
              else if (phylyr(i,j,k,pos) /= nf90_fill_double) then
                phylyr(i,j,k,pos) = log10(phylyr(i,j,k,pos)*sfac+offs)
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'v') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,kk
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (phylyr(i,j,k,pos) < epsil) then
                phylyr(i,j,k,pos) = 0.
              else if (phylyr(i,j,k,pos) /= nf90_fill_double) then
                phylyr(i,j,k,pos) = log10(phylyr(i,j,k,pos)*sfac+offs)
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'p') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (phylyr(i,j,k,pos) < epsil) then
                phylyr(i,j,k,pos) = 0.
              else if (phylyr(i,j,k,pos) /= nf90_fill_double) then
                phylyr(i,j,k,pos) = log10(phylyr(i,j,k,pos)*sfac+offs)
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else
      write (lp,*) 'cannot identify grid '//gridid//'!'
      call xchalt('(loglyr)')
      stop '(loglyr)'
    end if

  end subroutine loglyr


  subroutine loglvl(pos,gridid,sfac,offs)
  ! ------------------------------------------------------------------
  ! Description: replace 3d level field with log10(field)
  ! ------------------------------------------------------------------

    ! Arguments
    integer,          intent(in) :: pos    ! field position in layer buffer
    real,             intent(in) :: sfac   ! scale factor to be applied before log10
    real,             intent(in) :: offs   ! offset to be added before log10
    character(len=*), intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j,k,l
    real :: epsil = 1e-11

    ! Check whether field should be processed
    if (pos == 0) return

    if (gridid == 'u') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,ddm
          do l = 1,isu(j)
            do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (phylvl(i,j,k,pos) < epsil) then
                phylvl(i,j,k,pos) = 0.
              else if (phylvl(i,j,k,pos) /= nf90_fill_double) then
                phylvl(i,j,k,pos) = log10(phylvl(i,j,k,pos)*sfac+offs)
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'v') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,ddm
          do l = 1,isv(j)
            do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (phylvl(i,j,k,pos) < epsil) then
                phylvl(i,j,k,pos) = 0.
              else if (phylvl(i,j,k,pos) /= nf90_fill_double) then
                phylvl(i,j,k,pos) = log10(phylvl(i,j,k,pos)*sfac+offs)
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'p') then
      !$omp parallel do private(k,l,i)
      do j = 1,jj
        do k = 1,ddm
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (phylvl(i,j,k,pos) < epsil) then
                phylvl(i,j,k,pos) = 0.
              else if (phylvl(i,j,k,pos) /= nf90_fill_double) then
                phylvl(i,j,k,pos) = log10(phylvl(i,j,k,pos)*sfac+offs)
              end if
            end do
          end do
        end do
      end do
      !$omp end parallel do
    else
      write (lp,*) 'cannot identify grid '//gridid//'!'
      call xchalt('(loglvl)')
      stop '(loglvl)'
    end if

  end subroutine loglvl


  subroutine msklvl(pos,gridid)
  ! ------------------------------------------------------------------
  ! Description: set sea floor points to NaN in level fields
  ! ------------------------------------------------------------------

    integer,          intent(in) :: pos    ! field position in level buffer
    character(len=*), intent(in) :: gridid ! grid identifier ('p','u' or 'v')

    ! Local variables
    integer :: i,j,k
    logical :: iniflg = .true.
    integer, dimension(idm,jdm), save :: kmaxu,kmaxv,kmaxp
    real, parameter :: mskval = nf90_fill_double

    ! Check whether field should be processed
    if (pos == 0) return

    ! Prepare index fields for masking

    if (iniflg) then
      !$omp parallel do private(i,k)
      do j = 1,jj
        do i = 1,ii
          kmaxp(i,j) = 0
          kmaxu(i,j) = 0
          kmaxv(i,j) = 0
        end do
        do k = 1,ddm
          do i = 1,ii
            if (depths(i,j) > depthslev_bnds(1,k)) kmaxp(i,j) = k
            if (min(depths(i,j),depths(i-1,j)) > depthslev_bnds(1,k)) then
              kmaxu(i,j) = k
            end if
            if (min(depths(i,j),depths(i,j-1)) > depthslev_bnds(1,k)) then
              kmaxv(i,j) = k
            end if
          end do
        end do
      end do
      !$omp end parallel do
      iniflg = .false.
    end if

    if (gridid == 'u') then
      !$omp parallel do private(i,k)
      do j = 1,jj
        do i = 1,ii
          do k = kmaxu(i,j)+1,ddm
            phylvl(i,j,k,pos) = mskval
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'v') then
      !$omp parallel do private(i,k)
      do j = 1,jj
        do i = 1,ii
          do k = kmaxv(i,j)+1,ddm
            phylvl(i,j,k,pos) = mskval
          end do
        end do
      end do
      !$omp end parallel do
    else if (gridid == 'p') then
      !$omp parallel do private(i,k)
      do j = 1,jj
        do i = 1,ii
          do k = kmaxp(i,j)+1,ddm
            phylvl(i,j,k,pos) = mskval
          end do
        end do
      end do
      !$omp end parallel do
    else
      write (lp,*) 'cannot identify grid '//gridid//'!'
      call xchalt('(msklvl)')
      stop '(msklvl)'
    end if

  end subroutine msklvl


  subroutine definevar(irec,iogrp,cmpflg,timeunits,calendar)
  ! ------------------------------------------------------------------
  ! Define output variables
  ! ------------------------------------------------------------------

    ! Arguments
    integer,           intent(in) :: irec,iogrp,cmpflg
    character(len=30), intent(in) :: timeunits
    character(len=19), intent(in) :: calendar

    ! Local variables
    integer :: isize_lyr,nt,nat,km
    character :: trcnm*80,trcnml*80,dimslyr*100

    call ncdefvar('time','time',ndouble,0)
    call ncattr('long_name','time')
    call ncattr('units',timeunits)
    call ncattr('calendar',calendar)
    if (irec == 1) then
      ! define sigma levels
      if (vcoord_tag /= vcoord_isopyc_bulkml) then
        call ncdefvar('sigma','layer',ndouble,8)
      else
        call ncdefvar('sigma','sigma',ndouble,8)
      end if
      call ncattr('long_name','Potential density')
      call ncattr('standard_name','sea_water_sigma_theta')
      call ncattr('units','kg m-3')
      call ncattr('positive','down')
      ! define zlevel
      call ncdefvar('depth','depth',ndouble,8)
      call ncattr('long_name','z level')
      call ncattr('units','m')
      call ncattr('positive','down')
      call ncattr('bounds','depth_bnds')
      call ncdefvar('depth_bnds','bounds depth',ndouble,8)
      if (MSC_MMFLXL(iogrp)+MSC_MMFLXD(iogrp)+MSC_MMFTDL(iogrp) &
           +MSC_MMFSML(iogrp)+MSC_MMFTDD(iogrp)+MSC_MMFSMD(iogrp) &
           +MSC_MHFLX (iogrp)+MSC_MHFTD (iogrp)+MSC_MHFSM (iogrp) &
           +MSC_MHFLD (iogrp)+MSC_MSFLX (iogrp)+MSC_MSFTD (iogrp) &
           +msc_msfsm (iogrp)+msc_msfld (iogrp) /= 0) then
        call ncdefvar('lat','lat',ndouble,8)
        call ncattr('long_name','Latitude')
        call ncattr('standard_name','latitude')
        call ncattr('units','degree_north')
        call ncdefvar('region','slenmax region',nchar,0)
        call ncattr('long_name','Region name')
      end if
      if (msc_voltr(iogrp) /= 0) then
        call ncdefvar('section','slenmax section',nchar,0)
        call ncattr('long_name','Section name')
      end if
    end if

    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      isize_lyr=11
    else
      isize_lyr=1
    end if

    ! define 2d fields
    call ncdefvar3d(H2D_SIGMX(iogrp),cmpflg,'p','sigmx', &
         'Mixed layer density',' ','kg m-3',0)

    call ncdefvar3d(H2D_UB(iogrp),cmpflg,'u','ubaro', &
         'Barotropic velocity x-component',' ','m s-1',0)

    call ncdefvar3d(H2D_VB(iogrp),cmpflg,'v','vbaro', &
         'Barotropic velocity y-component',' ','m s-1',0)

    call ncdefvar3d(H2D_PSRF(iogrp),cmpflg,'p','psrf', &
         'Surface pressure',' ','Pa',0)

    call ncdefvar3d(H2D_PBOT(iogrp),cmpflg,'p','pbot', &
         'Bottom pressure',' ','Pa',0)

    call ncdefvar3d(H2D_SEALV(iogrp),cmpflg,'p','sealv', &
         'Sea level',' ','m',0)

    call ncdefvar3d(H2D_SLVSQ(iogrp),cmpflg,'p','slvsq', &
         'Sea level squared',' ','m2',0)

    call ncdefvar3d(H2D_BTMSTR(iogrp),cmpflg,'p','btmstr', &
         'Barotropic mass streamfunction',' ','kg s-1',0)

    call ncdefvar3d(H2D_HICE(iogrp),cmpflg,'p','hice', &
         'Ice thickness',' ','m',0)

    call ncdefvar3d(H2D_TICE(iogrp),cmpflg,'p','tice', &
         'Ice temperature',' ','degC',0)

    call ncdefvar3d(H2D_HSNW(iogrp),cmpflg,'p','hsnw', &
         'Snow depth',' ','m',0)

    call ncdefvar3d(H2D_FICE(iogrp),cmpflg,'p','fice', &
         'Ice concentration',' ','%',0)

    call ncdefvar3d(H2D_TSRF(iogrp),cmpflg,'p','tsrf', &
         'Surface temperature',' ','degC',0)

    call ncdefvar3d(H2D_IAGE(iogrp),cmpflg,'p','iage', &
         'Ice age',' ','day',0)

    call ncdefvar3d(H2D_UICE(iogrp),cmpflg,'u','uice', &
         'Ice velocity x-component',' ','m s-1',0)

    call ncdefvar3d(H2D_VICE(iogrp),cmpflg,'v','vice', &
         'Ice velocity y-component',' ','m s-1',0)

    call ncdefvar3d(H2D_SWA(iogrp),cmpflg,'p','swa', &
         'Short-wave heat flux',' ','W m-2',0)

    call ncdefvar3d(H2D_NSF(iogrp),cmpflg,'p','nsf', &
         'Non-solar heat flux',' ','W m-2',0)

    call ncdefvar3d(H2D_HMAT(iogrp),cmpflg,'p','hmat', &
         'Heat flux due to material enthalpy flux',' ','W m-2',0)

    call ncdefvar3d(H2D_HMLTFZ(iogrp),cmpflg,'p','hmltfz', &
         'Heat flux due to melting/freezing',' ','W m-2',0)

    call ncdefvar3d(H2D_DFL(iogrp),cmpflg,'p','dfl', &
         'Non-solar heat flux derivative',' ','W m-2 K-1',0)

    call ncdefvar3d(H2D_SURFLX(iogrp),cmpflg,'p','hflx', &
         'Heat flux received by ocean',' ','W m-2',0)

    call ncdefvar3d(H2D_SURRLX(iogrp),cmpflg,'p','hrflx', &
         'Restoring heat flux received by ocean',' ','W m-2',0)

    call ncdefvar3d(H2D_LIP(iogrp),cmpflg,'p','lip', &
         'Liquid precipitation',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_SOP(iogrp),cmpflg,'p','sop', &
         'Solid precipitation',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_EVA(iogrp),cmpflg,'p','eva', &
         'Evaporation',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_FMLTFZ(iogrp),cmpflg,'p','fmltfz', &
         'Fresh water flux due to melting/freezing',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_RNFFLX(iogrp),cmpflg,'p','rnf', &
         'Liquid runoff',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_RFIFLX(iogrp),cmpflg,'p','rfi', &
         'Frozen runoff',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_SALFLX(iogrp),cmpflg,'p','sflx', &
         'Salt flux received by ocean',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_SALRLX(iogrp),cmpflg,'p','srflx', &
         'Restoring salt flux received by ocean',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_BRNFLX(iogrp),cmpflg,'p','bflx', &
         'Brine flux',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_ZTX(iogrp),cmpflg,'u','ztx', &
         'Wind stress x-component',' ','N m-2',0)

    call ncdefvar3d(H2D_MTY(iogrp),cmpflg,'v','mty', &
         'Wind stress y-component',' ','N m-2',0)

    call ncdefvar3d(H2D_TAUX(iogrp),cmpflg,'u','taux', &
         'Momentum flux received by ocean x-component',' ','N m-2',0)

    call ncdefvar3d(H2D_TAUY(iogrp),cmpflg,'v','tauy', &
         'Momentum flux received by ocean y-component',' ','N m-2',0)

    call ncdefvar3d(H2D_IDKEDT(iogrp),cmpflg,'p','idkedt', &
         'Mixed layer inertial kinetic energy tendency per unit area', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_USTAR(iogrp),cmpflg,'p','ustar', &
         'Friction velocity',' ','m s-1',0)

    call ncdefvar3d(H2D_USTAR3(iogrp),cmpflg,'p','ustar3', &
         'Friction velocity cubed',' ','m3 s-3',0)

    call ncdefvar3d(H2D_ABSWND(iogrp),cmpflg,'p','abswnd', &
         'Absolute wind speed',' ','m s-1',0)

    call ncdefvar3d(H2D_MTKEUS(iogrp),cmpflg,'p','mtkeus', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to friction velocity', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_MTKENI(iogrp),cmpflg,'p','mtkeni', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to near inertial motions', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_MTKEBF(iogrp),cmpflg,'p','mtkebf', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to buoyancy forcing', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_MTKERS(iogrp),cmpflg,'p','mtkers', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to eddy restratification', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_MTKEPE(iogrp),cmpflg,'p','mtkepe', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to potential energy change', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_MTKEKE(iogrp),cmpflg,'p','mtkeke', &
         'Mixed layer turbulent kinetic energy tendency '// &
         'per unit area related to kinetic energy change', &
         ' ','kg s-3',0)

    call ncdefvar3d(H2D_LAMULT(iogrp),cmpflg,'p','lamult', &
         'Langmuir enhancement factor',' ','1',0)

    call ncdefvar3d(H2D_LASL(iogrp),cmpflg,'p','lasl', &
         'Surface layer averaged Langmuir number',' ','1',0)

    call ncdefvar3d(H2D_USTOKES(iogrp),cmpflg,'u','ustokes', &
         'Surface Stokes drift x-component',' ','m s-1',0)

    call ncdefvar3d(H2D_VSTOKES(iogrp),cmpflg,'v','vstokes', &
         'Surface Stokes drift y-component',' ','m s-1',0)

    call ncdefvar3d(H2D_SFL(iogrp),cmpflg,'p','sfl', &
         'Salt flux',' ','kg m-2 s-1',0)

    call ncdefvar3d(H2D_ALB(iogrp),cmpflg,'p','alb', &
         'Surface albedo',' ','1',0)

    call ncdefvar3d(H2D_MLD(iogrp),cmpflg,'p','mld', &
         'Mixed layer depth',' ','m',0)

    call ncdefvar3d(H2D_MAXMLD(iogrp),cmpflg,'p','maxmld', &
         'Maximum mixed layer depth',' ','m',0)

    call ncdefvar3d(H2D_MLTS(iogrp),cmpflg,'p','mlts', &
         'Mixed layer thickness defined by sigma t',' ','m',0)

    call ncdefvar3d(H2D_MLTSMN(iogrp),cmpflg,'p','mltsmn', &
         'Minimum mixed layer thickness defined by sigma t',' ','m',0)

    call ncdefvar3d(H2D_MLTSMX(iogrp),cmpflg,'p','mltsmx', &
         'Maximum mixed layer thickness defined by sigma t',' ','m',0)

    call ncdefvar3d(H2D_MLTSSQ(iogrp),cmpflg,'p','mltssq', &
         'Maximum mixed layer thickness squared defined by sigma t',' ', &
         'm2',0)

    call ncdefvar3d(H2D_T20D(iogrp),cmpflg,'p','t20d', &
         '20C isoterm depth',' ','m',0)

    call ncdefvar3d(H2D_BRNPD(iogrp),cmpflg,'p','brnpd', &
         'Brine plume depth',' ','m',0)

    call ncdefvar3d(H2D_SSS(iogrp),cmpflg,'p','sss', &
         'Ocean surface salinity',' ','g kg-1',0)

    call ncdefvar3d(H2D_SSSSQ(iogrp),cmpflg,'p','ssssq', &
         'Ocean surface salinity squared',' ','g2 kg-2',0)

    call ncdefvar3d(H2D_SBOT(iogrp),cmpflg,'p','sbot', &
         'Bottom salinity',' ','g kg-1',0)

    call ncdefvar3d(H2D_SST(iogrp),cmpflg,'p','sst', &
         'Ocean surface temperature',' ','degC',0)

    call ncdefvar3d(H2D_SSTSQ(iogrp),cmpflg,'p','sstsq', &
         'Ocean surface temperature squared',' ','degC2',0)

    call ncdefvar3d(H2D_TBOT(iogrp),cmpflg,'p','tbot', &
         'Bottom temperature',' ','degC',0)

    ! define 3d layer fields
    call ncdefvar3d(LYR_DP(iogrp),cmpflg,'p','dp', &
         'Layer pressure thickness',' ','Pa',isize_lyr)

    call ncdefvar3d(LYR_DZ(iogrp),cmpflg,'p','dz', &
         'Layer thickness',' ','m',isize_lyr)

    call ncdefvar3d(LYR_TEMP(iogrp),cmpflg,'p','temp', &
         'Temperature','Ocean temperature','degC',isize_lyr)

    call ncdefvar3d(LYR_SALN(iogrp),cmpflg,'p','saln', &
         'Salinity','Ocean salinity','g kg-1',isize_lyr)

    call ncdefvar3d(LYR_UVEL(iogrp),cmpflg,'u','uvel', &
         'Velocity x-component',' ','m s-1',isize_lyr)

    call ncdefvar3d(LYR_VVEL(iogrp),cmpflg,'v','vvel', &
         'Velocity y-component',' ','m s-1',isize_lyr)

    call ncdefvar3d(LYR_UFLX(iogrp),cmpflg,'u','uflx', &
         'Mass flux in x-direction',' ','kg s-1',isize_lyr)

    call ncdefvar3d(LYR_VFLX(iogrp),cmpflg,'v','vflx', &
         'Mass flux in y-direction',' ','kg s-1',isize_lyr)

    call ncdefvar3d(LYR_UTFLX(iogrp),cmpflg,'u','uhflx', &
         'Heat flux in x-direction',' ','W',isize_lyr)

    call ncdefvar3d(LYR_VTFLX(iogrp),cmpflg,'v','vhflx', &
         'Heat flux in y-direction',' ','W',isize_lyr)

    call ncdefvar3d(LYR_USFLX(iogrp),cmpflg,'u','usflx', &
         'Salt flux in x-direction',' ','kg s-1',isize_lyr)

    call ncdefvar3d(LYR_VSFLX(iogrp),cmpflg,'v','vsflx', &
         'Salt flux in y-direction',' ','kg s-1',isize_lyr)

    call ncdefvar3d(LYR_UMFLTD(iogrp),cmpflg,'u','umfltd', &
         'Mass flux due to thickness diffusion in x-direction',' ', &
         'kg s-1',isize_lyr)

    call ncdefvar3d(LYR_VMFLTD(iogrp),cmpflg,'v','vmfltd', &
         'Mass flux due to thickness diffusion in y-direction',' ', &
         'kg s-1',isize_lyr)

    call ncdefvar3d(LYR_UMFLSM(iogrp),cmpflg,'u','umflsm', &
         'Mass flux due to submesoscale transport in x-direction',' ', &
         'kg s-1',isize_lyr)

    call ncdefvar3d(LYR_VMFLSM(iogrp),cmpflg,'v','vmflsm', &
         'Mass flux due to submesoscale transport in y-direction',' ', &
         'kg s-1',isize_lyr)

    call ncdefvar3d(LYR_UTFLTD(iogrp),cmpflg,'u','uhfltd', &
         'Heat flux due to thickness diffusion in x-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_VTFLTD(iogrp),cmpflg,'v','vhfltd', &
         'Heat flux due to thickness diffusion in y-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_UTFLSM(iogrp),cmpflg,'u','uhflsm', &
         'Heat flux due to submesoscale transport in x-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_VTFLSM(iogrp),cmpflg,'v','vhflsm', &
         'Heat flux due to submesoscale transport in y-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_UTFLLD(iogrp),cmpflg,'u','uhflld', &
         'Heat flux due to lateral diffusion in x-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_VTFLLD(iogrp),cmpflg,'v','vhflld', &
         'Heat flux due to lateral diffusion in y-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_USFLTD(iogrp),cmpflg,'u','usfltd', &
         'Salt flux due to thickness diffusion in x-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_VSFLTD(iogrp),cmpflg,'v','vsfltd', &
         'Salt flux due to thickness diffusion in y-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_USFLSM(iogrp),cmpflg,'u','usflsm', &
         'Salt flux due to submesoscale transport in x-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_VSFLSM(iogrp),cmpflg,'v','vsflsm', &
         'Salt flux due to submesoscale transport in y-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_USFLLD(iogrp),cmpflg,'u','usflld', &
         'Salt flux due to lateral diffusion in x-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_VSFLLD(iogrp),cmpflg,'v','vsflld', &
         'Salt flux due to lateral diffusion in y-direction',' ', &
         'W',isize_lyr)

    call ncdefvar3d(LYR_WFLX(iogrp),cmpflg,'p','wflx', &
         'Vertical mass flux',' ','kg s-1',isize_lyr)

    call ncdefvar3d(LYR_WFLX2(iogrp),cmpflg,'p','wflx2', &
         'Vertical mass flux squared',' ','kg2 s-2',isize_lyr)

    call ncdefvar3d(LYR_BFSQ(iogrp),cmpflg,'p','bfsq', &
         'Squared buoyancy frequency',' ','s-1',isize_lyr)

    call ncdefvar3d(LYR_PV(iogrp),cmpflg,'p','pv', &
         'Potential vorticity',' ','m-1 s-1',isize_lyr)

    if (lyr_difint(iogrp) == 2) then
      call ncdefvar3d(LYR_DIFINT(iogrp),cmpflg,'p','difint', &
           'Layer interface diffusivity',' ','log10(m2 s-1)',isize_lyr)
    else
      call ncdefvar3d(LYR_DIFINT(iogrp),cmpflg,'p','difint', &
           'Layer interface diffusivity',' ','m2 s-1',isize_lyr)
    end if

    if (lyr_difiso(iogrp) == 2) then
      call ncdefvar3d(LYR_DIFISO(iogrp),cmpflg,'p','difiso', &
           'Isopycnal diffusivity',' ','log10(m2 s-1)',isize_lyr)
    else
      call ncdefvar3d(LYR_DIFISO(iogrp),cmpflg,'p','difiso', &
           'Isopycnal diffusivity',' ','m2 s-1',isize_lyr)
    end if

    if (lyr_difdia(iogrp) == 2) then
      call ncdefvar3d(LYR_DIFDIA(iogrp),cmpflg,'p','difdia', &
           'Vertical diffusivity',' ','log10(m2 s-1)',isize_lyr)
    else
      call ncdefvar3d(LYR_DIFDIA(iogrp),cmpflg,'p','difdia', &
           'Vertical diffusivity',' ','m2 s-1',isize_lyr)
    end if

    if (lyr_difvmo(iogrp) == 2) then
      call ncdefvar3d(LYR_DIFVMO(iogrp),cmpflg,'p','difvmo', &
           'Vertical momentum diffusivity',' ','log10(m2 s-1)',isize_lyr)
    else
      call ncdefvar3d(LYR_DIFVMO(iogrp),cmpflg,'p','difvmo', &
           'Vertical momentum diffusivity',' ','m2 s-1',isize_lyr)
    end if

    if (lyr_difvho(iogrp) == 2) then
      call ncdefvar3d(LYR_DIFVHO(iogrp),cmpflg,'p','difvho', &
           'Vertical heat diffusivity',' ','log10(m2 s-1)',isize_lyr)
    else
      call ncdefvar3d(LYR_DIFVHO(iogrp),cmpflg,'p','difvho', &
           'Vertical heat diffusivity',' ','m2 s-1',isize_lyr)
    end if

    if (lyr_difvso(iogrp) == 2) then
      call ncdefvar3d(LYR_DIFVSO(iogrp),cmpflg,'p','difvso', &
           'Vertical salt diffusivity',' ','log10(m2 s-1)',isize_lyr)
    else
      call ncdefvar3d(LYR_DIFVSO(iogrp),cmpflg,'p','difvso', &
           'Vertical salt diffusivity',' ','m2 s-1',isize_lyr)
    end if

    if (use_TKE) then
      call ncdefvar3d(LYR_TKE(iogrp),cmpflg,'p','tke', &
           'TKE','Turbulent kinetic energy','m2 s-2',isize_lyr)

      call ncdefvar3d(LYR_GLS_PSI(iogrp),cmpflg,'p','gls_psi', &
           'GLS_PSI','Generic length scale','m2 s-3',isize_lyr)
    end if

    ! define 3d depth fields
    call ncdefvar3d(LVL_DZ(iogrp),cmpflg,'p','dzlvl', &
         'Layer thickness',' ','m',2)

    call ncdefvar3d(LVL_TEMP(iogrp),cmpflg,'p','templvl', &
         'Temperature','Ocean temperature','degC',2)

    call ncdefvar3d(LVL_SALN(iogrp),cmpflg,'p','salnlvl', &
         'Salinity','Ocean salinity','g kg-1',2)

    call ncdefvar3d(LVL_UVEL(iogrp),cmpflg,'u','uvellvl', &
         'Velocity x-component',' ','m s-1',2)

    call ncdefvar3d(LVL_VVEL(iogrp),cmpflg,'v','vvellvl', &
         'Velocity y-component',' ','m s-1',2)

    call ncdefvar3d(LVL_UFLX(iogrp),cmpflg,'u','uflxlvl', &
         'Mass flux in x-direction',' ','kg s-1',2)

    call ncdefvar3d(LVL_VFLX(iogrp),cmpflg,'v','vflxlvl', &
         'Mass flux in y-direction',' ','kg s-1',2)

    call ncdefvar3d(LVL_UTFLX(iogrp),cmpflg,'u','uhflxlvl', &
         'Heat flux in x-direction',' ','W',2)

    call ncdefvar3d(LVL_VTFLX(iogrp),cmpflg,'v','vhflxlvl', &
         'Heat flux in y-direction',' ','W',2)

    call ncdefvar3d(LVL_USFLX(iogrp),cmpflg,'u','usflxlvl', &
         'Salt flux in x-direction',' ','kg s-1',2)

    call ncdefvar3d(LVL_VSFLX(iogrp),cmpflg,'v','vsflxlvl', &
         'Salt flux in y-direction',' ','kg s-1',2)

    call ncdefvar3d(LVL_UMFLTD(iogrp),cmpflg,'u','umfltdlvl', &
         'Mass flux due to thickness diffusion in x-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_VMFLTD(iogrp),cmpflg,'v','vmfltdlvl', &
         'Mass flux due to thickness diffusion in y-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_UMFLSM(iogrp),cmpflg,'u','umflsmlvl', &
         'Mass flux due to submesoscale transport in x-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_VMFLSM(iogrp),cmpflg,'v','vmflsmlvl', &
         'Mass flux due to submesoscale transport in y-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_UTFLTD(iogrp),cmpflg,'u','uhfltdlvl', &
         'Heat flux due to thickness diffusion in x-direction',' ', &
         'W',2)

    call ncdefvar3d(LVL_VTFLTD(iogrp),cmpflg,'v','vhfltdlvl', &
         'Heat flux due to thickness diffusion in y-direction',' ', &
         'W',2)

    call ncdefvar3d(LVL_UTFLSM(iogrp),cmpflg,'u','uhflsmlvl', &
         'Heat flux due to submesoscale transport in x-direction',' ', &
         'W',2)

    call ncdefvar3d(LVL_VTFLSM(iogrp),cmpflg,'v','vhflsmlvl', &
         'Heat flux due to submesoscale transport in y-direction',' ', &
         'W',2)

    call ncdefvar3d(LVL_UTFLLD(iogrp),cmpflg,'u','uhflldlvl', &
         'Heat flux due to lateral diffusion in x-direction',' ', &
         'W',2)

    call ncdefvar3d(LVL_VTFLLD(iogrp),cmpflg,'v','vhflldlvl', &
         'Heat flux due to lateral diffusion in y-direction',' ', &
         'W',2)

    call ncdefvar3d(LVL_USFLTD(iogrp),cmpflg,'u','usfltdlvl', &
         'Salt flux due to thickness diffusion in x-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_VSFLTD(iogrp),cmpflg,'v','vsfltdlvl', &
         'Salt flux due to thickness diffusion in y-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_USFLSM(iogrp),cmpflg,'u','usflsmlvl', &
         'Salt flux due to submesoscale transport in x-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_VSFLSM(iogrp),cmpflg,'v','vsflsmlvl', &
         'Salt flux due to submesoscale transport in y-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_USFLLD(iogrp),cmpflg,'u','usflldlvl', &
         'Salt flux due to lateral diffusion in x-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_VSFLLD(iogrp),cmpflg,'v','vsflldlvl', &
         'Salt flux due to lateral diffusion in y-direction',' ', &
         'kg s-1',2)

    call ncdefvar3d(LVL_WFLX(iogrp),cmpflg,'p','wflxlvl', &
         'Vertical mass flux',' ','kg s-1',2)

    call ncdefvar3d(LVL_WFLX2(iogrp),cmpflg,'p','wflx2lvl', &
         'Vertical mass flux squared',' ','kg2 s-2',2)

    call ncdefvar3d(LVL_BFSQ(iogrp),cmpflg,'p','bfsqlvl', &
         'Squared buoyancy frequancy',' ','s-1',2)

    call ncdefvar3d(LVL_PV(iogrp),cmpflg,'p','pvlvl', &
         'Potential vorticity',' ','m-1 s-1',2)

    if (lvl_difint(iogrp) == 2) then
      call ncdefvar3d(LVL_DIFINT(iogrp),cmpflg,'p','difintlvl', &
           'Layer interface diffusivity',' ','log10(m2 s-1)',2)
    else
      call ncdefvar3d(LVL_DIFINT(iogrp),cmpflg,'p','difintlvl', &
           'Layer interface diffusivity',' ','m2 s-1',2)
    end if

    if (lvl_difiso(iogrp) == 2) then
      call ncdefvar3d(LVL_DIFISO(iogrp),cmpflg,'p','difisolvl', &
           'Isopycnal diffusivity',' ','log10(m2 s-1)',2)
    else
      call ncdefvar3d(LVL_DIFISO(iogrp),cmpflg,'p','difisolvl', &
           'Isopycnal diffusivity',' ','m2 s-1',2)
    end if

    if (lvl_difdia(iogrp) == 2) then
      call ncdefvar3d(LVL_DIFDIA(iogrp),cmpflg,'p','difdialvl', &
           'Vertical diffusivity',' ','log10(m2 s-1)',2)
    else
      call ncdefvar3d(LVL_DIFDIA(iogrp),cmpflg,'p','difdialvl', &
           'Vertical diffusivity',' ','m2 s-1',2)
    end if

    if (lvl_difvmo(iogrp) == 2) then
      call ncdefvar3d(LVL_DIFVMO(iogrp),cmpflg,'p','difvmolvl', &
           'Vertical momentum diffusivity',' ','log10(m2 s-1)',2)
    else
      call ncdefvar3d(LVL_DIFVMO(iogrp),cmpflg,'p','difvmolvl', &
           'Vertical momentum diffusivity',' ','m2 s-1',2)
    end if

    if (lvl_difvho(iogrp) == 2) then
      call ncdefvar3d(LVL_DIFVHO(iogrp),cmpflg,'p','difvholvl', &
           'Vertical heat diffusivity',' ','log10(m2 s-1)',2)
    else
      call ncdefvar3d(LVL_DIFVHO(iogrp),cmpflg,'p','difvholvl', &
           'Vertical heat diffusivity',' ','m2 s-1',2)
    end if

    if (lvl_difvso(iogrp) == 2) then
      call ncdefvar3d(LVL_DIFVSO(iogrp),cmpflg,'p','difvsolvl', &
           'Vertical salt diffusivity',' ','log10(m2 s-1)',2)
    else
      call ncdefvar3d(LVL_DIFVSO(iogrp),cmpflg,'p','difvsolvl', &
           'Vertical salt diffusivity',' ','m2 s-1',2)
    end if

    if (use_TRC .and. use_TKE) then
      call ncdefvar3d(LVL_TKE(iogrp),cmpflg,'p','tkelvl', &
           'Turbulent Kinetic Energy',' ','m2 s-2',2)

      call ncdefvar3d(LVL_GLS_PSI(iogrp),cmpflg,'p','gls_psilvl', &
           'Generic length scale',' ','m2 s-3',2)
    endif

    ! define meridional transports
    if (vcoord_tag /= vcoord_isopyc_bulkml) then
      dimslyr='lat layer region time'
    else
      dimslyr='lat sigma region time'
    end if
    if (msc_mmflxl(iogrp) /= 0) then
      call ncdefvar('mmflxl',dimslyr,ndouble,8)
      call ncattr('long_name', &
           'Overturning stream-function on layers')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmflxd(iogrp) /= 0) then
      call ncdefvar('mmflxd','lat depth region time',ndouble,8)
      call ncattr('long_name', &
           'Overturning stream-function on z-levels')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmftdl(iogrp) /= 0) then
      call ncdefvar('mmftdl',dimslyr,ndouble,8)
      call ncattr('long_name', &
           'Overturning stream-function due to thickness diffusion '// &
           'on layers')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmfsml(iogrp) /= 0) then
      call ncdefvar('mmfsml',dimslyr,ndouble,8)
      call ncattr('long_name', &
           'Overturning stream-function due to submesoscale transport '// &
           'on layers')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmftdd(iogrp) /= 0) then
      call ncdefvar('mmftdd','lat depth region time',ndouble,8)
      call ncattr('long_name', &
           'Overturning stream-function due to thickness diffusion '// &
           'on z-levels')
      call ncattr('units','kg s-1')
    end if
    if (msc_mmfsmd(iogrp) /= 0) then
      call ncdefvar('mmfsmd','lat depth region time',ndouble,8)
      call ncattr('long_name', &
           'Overturning stream-function due to submesoscale transport '// &
           'on z-levels')
      call ncattr('units','kg s-1')
    end if
    if (msc_mhflx(iogrp) /= 0) then
      call ncdefvar('mhflx','lat region time',ndouble,8)
      call ncattr('long_name','Meridional heat flux')
      call ncattr('units','W')
    end if
    if (msc_mhftd(iogrp) /= 0) then
      call ncdefvar('mhftd','lat region time',ndouble,8)
      call ncattr('long_name', &
           'Meridional heat flux due to thickness diffusion')
      call ncattr('units','W')
    end if
    if (msc_mhfsm(iogrp) /= 0) then
      call ncdefvar('mhfsm','lat region time',ndouble,8)
      call ncattr('long_name', &
           'Meridional heat flux due to submesoscale transport')
      call ncattr('units','W')
    end if
    if (msc_mhfld(iogrp) /= 0) then
      call ncdefvar('mhfld','lat region time',ndouble,8)
      call ncattr('long_name', &
           'Meridional heat flux due to lateral diffusion')
      call ncattr('units','W')
    end if
    if (msc_msflx(iogrp) /= 0) then
      call ncdefvar('msflx','lat region time',ndouble,8)
      call ncattr('long_name','Meridional salt flux')
      call ncattr('units','kg s-1')
    end if
    if (msc_msftd(iogrp) /= 0) then
      call ncdefvar('msftd','lat region time',ndouble,8)
      call ncattr('long_name', &
           'Meridional salt flux due to thickness diffusion')
      call ncattr('units','kg s-1')
    end if
    if (msc_msfsm(iogrp) /= 0) then
      call ncdefvar('msfsm','lat region time',ndouble,8)
      call ncattr('long_name', &
           'Meridional salt flux due to submesoscale transport')
      call ncattr('units','kg s-1')
    end if
    if (msc_msfld(iogrp) /= 0) then
      call ncdefvar('msfld','lat region time',ndouble,8)
      call ncattr('long_name', &
           'Meridional salt flux due to lateral diffusion')
      call ncattr('units','kg s-1')
    end if

    ! store section transports
    if (msc_voltr(iogrp) /= 0) then
      call ncdefvar('voltr','section time',ndouble,8)
      call ncattr('long_name','Section transports')
      call ncattr('units','kg s-1')
    end if

    ! store global sums and averages
    if (msc_massgs(iogrp) /= 0) then
      call ncdefvar('massgs','time',ndouble,8)
      call ncattr('long_name','Sea water mass')
      call ncattr('units','kg')
    end if
    if (msc_volgs(iogrp) /= 0) then
      call ncdefvar('volgs','time',ndouble,8)
      call ncattr('long_name','Sea water volume')
      call ncattr('units','m3')
    end if
    if (msc_salnga(iogrp) /= 0) then
      call ncdefvar('salnga','time',ndouble,8)
      call ncattr('long_name','Global average salinity')
      call ncattr('units','g kg-1')
    end if
    if (msc_tempga(iogrp) /= 0) then
      call ncdefvar('tempga','time',ndouble,8)
      call ncattr('long_name','Global average temperature')
      call ncattr('units','degC')
    end if
    if (msc_sssga(iogrp) /= 0) then
      call ncdefvar('sssga','time',ndouble,8)
      call ncattr('long_name','Global average sea surface salinity')
      call ncattr('units','g kg-1')
    end if
    if (msc_sstga(iogrp) /= 0) then
      call ncdefvar('sstga','time',ndouble,8)
      call ncattr('long_name', &
           'Global average sea surface temperature')
      call ncattr('units','degC')
    end if

    if (use_TRC) then
      if (lyr_idlage(iogrp) /= 0.or.lyr_trc(iogrp) /= 0) then
        call ncdefvar3d(max(LYR_IDLAGE(iogrp),LYR_TRC(iogrp)),cmpflg, &
             'p','dp_trc','Layer pressure thickness',' ','Pa',isize_lyr)
      end if
      ! ideal age tracer
      if (use_IDLAGE) then
        call ncdefvar3d(LYR_IDLAGE(iogrp),cmpflg,'p','idlage', &
             'Ideal age','sea_water_age_since_surface_contact','year',isize_lyr)

        if (lvl_idlage(iogrp) /= 0) then
          call ncdefvar3d(LVL_IDLAGE(iogrp),cmpflg,'p','idlagelvl', &
               'Ideal age','sea_water_age_since_surface_contact','year',2)
        end if
      end if

      ! ocean tracers
      if (lyr_trc(iogrp) > 0.and.ntrocn > 0) then
        if (use_ATRC) then
          do nt = 1,ntrocn-natr
            write (trcnm,'(a,i3.3)') 'trc',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call ncdefvar3d(LYR_TRC(iogrp),cmpflg,'p',trim(trcnm), &
                 trim(trcnml),' ',' ',isize_lyr)
          end do
          do nt = 1,natr
            nat = ntr-natr+nt
            write (trcnm,'(a,i3.3)') 'atrc',nt
            write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
            call ncdefvar3d(LYR_TRC(iogrp),cmpflg,'p',trim(trcnm), &
                 trim(trcnml),' ',' ',isize_lyr)
          end do
        else
          do nt = 1,ntrocn
            write (trcnm,'(a,i3.3)') 'trc',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call ncdefvar3d(LYR_TRC(iogrp),cmpflg,'p',trim(trcnm), &
                 trim(trcnml),' ',' ',isize_lyr)
          end do
        end if
      end if
      if (lvl_trc(iogrp) > 0.and.ntrocn > 0) then
        if (use_ATRC) then
          do nt = 1,ntrocn-natr
            write (trcnm,'(a,i3.3)') 'trclvl',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call ncdefvar3d(LVL_TRC(iogrp),cmpflg,'p',trim(trcnm), &
                 trim(trcnml),' ',' ',2)
          end do
          do nt = 1,natr
            nat = ntr-natr+nt
            write (trcnm,'(a,i3.3)') 'atrclvl',nt
            write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
            call ncdefvar3d(LVL_TRC(iogrp),cmpflg,'p',trim(trcnm), &
                 trim(trcnml),' ',' ',2)
          end do
        else
          do nt = 1,ntrocn
            write (trcnm,'(a,i3.3)') 'trclvl',nt
            write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
            call ncdefvar3d(LVL_TRC(iogrp),cmpflg,'p',trim(trcnm), &
                 trim(trcnml),' ',' ',2)
          end do
        end if
      end if
    end if

    call ncedef

  end subroutine definevar

end module mod_dia
