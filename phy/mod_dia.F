! ------------------------------------------------------------------------------
! Copyright (C) 2010-2020 Ingo Bethke, Mats Bentsen, Mehmet Ilicak,
!                         Alok Kumar Gupta, Jörg Schwinger
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
c
      use mod_xc 
      use mod_nctools
      use netcdf, only : nf90_fill_double
c
      implicit none
c
c --- ------------------------------------------------------------------
c --- common blocks related to the accumulation and averaging of
c --- diagnostic fields
c --- ------------------------------------------------------------------
c
c --- Averaging and writing frequencies for diagnostics output
      integer, save :: nphy
      integer, parameter :: nphymax=10
      real, dimension(nphymax), save :: diagfq_phy,filefq_phy
      integer, dimension(nphymax), save :: nacc_phy
      logical, dimension(nphymax), save ::
     .  diagmon_phy,filemon_phy,diagann_phy,fileann_phy
c
c --- Restart parameters
      real, save :: rstfrq
      integer, save :: rstcmp,rstfmt,iotype
      logical, save :: rstmon,rstann
c
c --- Copies of BLOM variables that are used for HAMOCC diagnostics
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  pbath,ubath,vbath 
      integer, save :: nstepinday 
 
c --- 2d and 3d diagnostic variables     
      integer, save :: nphyh2d,nphylyr,nphylvl 
      real, save, allocatable, dimension(:,:,:) :: phyh2d
      real, save, allocatable, dimension(:,:,:,:) :: phylyr,phylvl 
c
c --- Levitus levels 
#ifndef LEVITUS2X 
      integer, parameter :: ddm=35,k350=12
      real, parameter :: w350=1.
      real,  parameter, dimension(ddm) :: depthslev=(/
     .  0000.0,0010.0,0020.0,0030.0,0050.0,0075.0,0100.0,0125.0,0150.0, 
     .  0200.0,0250.0,0300.0,0400.0,0500.0,0600.0,0700.0,0800.0,0900.0, 
     .  1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1750.0,2000.0,2500.0, 
     .  3000.0,3500.0,4000.0,4500.0,5000.0,5500.0,6000.0,6500.0/)
      real,  parameter, dimension(2,ddm) :: 
     .  depthslev_bnds=reshape((/
     .  0000.0,0005.0,0005.0,0015.0,0015.0,0025.0,0025.0,0040.0,0040.0,
     .  0062.5,0062.5,0087.5,0087.5,0112.5,0112.5,0137.5,0137.5,0175.0,
     .  0175.0,0225.0,0225.0,0275.0,0275.0,0350.0,0350.0,0450.0,0450.0,
     .  0550.0,0550.0,0650.0,0650.0,0750.0,0750.0,0850.0,0850.0,0950.0,
     .  0950.0,1050.0,1050.0,1150.0,1150.0,1250.0,1250.0,1350.0,1350.0,
     .  1450.0,1450.0,1625.0,1625.0,1875.0,1875.0,2250.0,2250.0,2750.0,
     .  2750.0,3250.0,3250.0,3750.0,3750.0,4250.0,4250.0,4750.0,4750.0,
     .  5250.0,5250.0,5750.0,5750.0,6250.0,6250.0,8000.0/),(/2,ddm/))
#else 
      integer, parameter :: ddm=70,k350=25
      real, parameter :: w350=0.5
      real,  parameter, dimension(ddm) :: depthslev=(/
     .  0000.0,0005.0,0010.0,0015.0,0020.0,0025.0,0030.0,0040.0,0050.0, 
     .  0062.5,0075.0,0087.5,0100.0,0112.5,0125.0,0137.5,0150.0,0175.0, 
     .  0200.0,0225.0,0250.0,0275.0,0300.0,0350.0,0400.0,0450.0,0500.0, 
     .  0550.0,0600.0,0650.0,0700.0,0750.0,0800.0,0850.0,0900.0,0950.0, 
     .  1000.0,1050.0,1100.0,1150.0,1200.0,1250.0,1300.0,1350.0,1400.0, 
     .  1450.0,1500.0,1625.0,1750.0,1875.0,2000.0,2250.0,2500.0,2750.0, 
     .  3000.0,3250.0,3500.0,3750.0,4000.0,4250.0,4500.0,4750.0,5000.0, 
     .  5250.0,5500.0,5750.0,6000.0,6250.0,6500.0,6750.0/)
      real,  parameter, dimension(2,ddm) :: 
     .  depthslev_bnds=reshape((/
     .  0000.0,0002.5,0002.5,0007.5,0007.5,0012.5,0012.5,0017.5,0017.5, 
     .  0022.5,0022.5,0027.5,0027.5,0035.0,0035.0,0045.0,0045.0,0056.2, 
     .  0056.2,0068.8,0068.8,0081.2,0081.2,0093.8,0093.8,0106.2,0106.2, 
     .  0118.8,0118.8,0131.2,0131.2,0143.8,0143.8,0162.5,0162.5,0187.5, 
     .  0187.5,0212.5,0212.5,0237.5,0237.5,0262.5,0262.5,0287.5,0287.5, 
     .  0325.0,0325.0,0375.0,0375.0,0425.0,0425.0,0475.0,0475.0,0525.0, 
     .  0525.0,0575.0,0575.0,0625.0,0625.0,0675.0,0675.0,0725.0,0725.0, 
     .  0775.0,0775.0,0825.0,0825.0,0875.0,0875.0,0925.0,0925.0,0975.0, 
     .  0975.0,1025.0,1025.0,1075.0,1075.0,1125.0,1125.0,1175.0,1175.0, 
     .  1225.0,1225.0,1275.0,1275.0,1325.0,1325.0,1375.0,1375.0,1425.0, 
     .  1425.0,1475.0,1475.0,1562.5,1562.5,1687.5,1687.5,1812.5,1812.5, 
     .  1937.5,1937.5,2125.0,2125.0,2375.0,2375.0,2625.0,2625.0,2875.0, 
     .  2875.0,3125.0,3125.0,3375.0,3375.0,3625.0,3625.0,3875.0,3875.0, 
     .  4125.0,4125.0,4375.0,4375.0,4625.0,4625.0,4875.0,4875.0,5125.0, 
     .  5125.0,5375.0,5375.0,5625.0,5625.0,5875.0,5875.0,6125.0,6125.0, 
     .  6375.0,6375.0,6625.0,6625.0,8000.0/),(/2,ddm/))
#endif 
c
c --- Meridional overturning and flux diagnostics 
      integer, parameter ::
     .  ldm=itdm+jtdm,sdm=ldm,odm=10,slenmax=50,rflgdm=20
      character(len=slenmax), save, dimension(odm) :: mer_regnam=''
      character(len=256), save :: mer_orfile,mer_mifile
      integer, save, dimension(odm,rflgdm) :: mer_regflg=-1
      integer, save, dimension(odm) :: mer_nflg
      real, save, dimension(odm) :: mer_minlat=-90.,mer_maxlat=90.
      integer, save :: mer_nreg,lmax 
      real, save, dimension(ldm) :: mtlat
      real, save, dimension(kdm) :: sigmar1 
      real, save, allocatable, dimension(:,:,:) ::
     .  mmflxl,mmftdl,mmflxd,mmftdd
      real, save, allocatable, dimension(:,:) ::
     .  mhflx,mhftd,mhfld,msflx,msftd,msfld
c
c --- Section transports 
      character(len=256), save :: sec_sifile
      integer, save :: sec_num
      integer, parameter :: max_sec=400
      character(len=slenmax), save :: sec_name(max_sec)
      real, save, dimension(max_sec) :: voltr
c
c --- Global sums and averages
      real, save, dimension(1) :: massgs,volgs,salnga,tempga,sssga,sstga
c
c --- Pressure thickness [g cm-1 s-2] of region for bottom salinity and
c --- temperature diagnostics
      real, parameter :: dpbot=98060.
c
c --- Critical buoyancy difference [cm s-2] used in the mixed layer
c --- thickness estimation (Levitus, 1982)
      real, parameter :: dbcrit=.03
c
c --- Namelist 
      integer, dimension(nphymax), save ::
     . H2D_ABSWND ,H2D_ALB    ,H2D_BTMSTR ,H2D_BRNFLX ,H2D_BRNPD  ,
     . H2D_DFL    ,H2D_EVA    ,H2D_FICE   ,H2D_FMLTFZ ,H2D_HICE   ,
     . H2D_HMLTFZ ,H2D_HSNW   ,H2D_IAGE   ,H2D_IDKEDT ,H2D_LIP    ,
     . H2D_MAXMLD ,H2D_MLD    ,H2D_MLDU   ,H2D_MLDV   ,H2D_MLTS   ,
     . H2D_MLTSMN ,H2D_MLTSMX ,H2D_MLTSSQ ,H2D_MTKEUS ,H2D_MTKENI ,
     . H2D_MTKEBF ,H2D_MTKERS ,H2D_MTKEPE ,H2D_MTKEKE ,H2D_MTY    ,
     . H2D_MXLU   ,H2D_MXLV   ,H2D_NSF    ,H2D_PBOT   ,H2d_PSRF   ,
     . H2D_RFIFLX ,H2D_RNFFLX ,H2D_SALFLX ,H2D_SALRLX ,H2D_SBOT   ,
     . H2D_SEALV  ,H2D_SLVSQ  ,H2D_SFL    ,H2D_SOP    ,H2D_SIGMX  ,
     . H2D_SSS    ,H2D_SSSSQ  ,H2D_SST    ,H2D_SSTSQ  ,H2D_SURFLX ,
     . H2D_SURRLX ,H2D_SWA    ,H2D_T20D   ,H2D_TAUX   ,H2D_TAUY   ,
     . H2D_TBOT   ,H2D_TICE   ,H2D_TSRF   ,H2D_UB     ,H2D_UICE   ,
     . H2D_USTAR  ,H2D_USTAR3 ,H2D_VB     ,H2D_VICE   ,H2D_ZTX    ,
     . LYR_BFSQ   ,LYR_DIFDIA ,LYR_DIFINT ,LYR_DIFISO ,LYR_DP     ,
     . LYR_DPU    ,LYR_DPV    ,LYR_DZ     ,LYR_SALN   ,LYR_TEMP   ,
     . LYR_TRC    ,LYR_UFLX   ,LYR_UTFLX  ,LYR_USFLX  ,LYR_UMFLTD ,
     . LYR_UTFLTD ,LYR_UTFLLD ,LYR_USFLTD ,LYR_USFLLD ,LYR_UVEL   ,
     . LYR_VFLX   ,LYR_VTFLX  ,LYR_VSFLX  ,LYR_VMFLTD ,LYR_VTFLTD ,
     . LYR_VTFLLD ,LYR_VSFLTD ,LYR_VSFLLD ,LYR_VVEL   ,LYR_WFLX   ,
     . LYR_WFLX2  ,LYR_PV     ,LYR_TKE    ,LYR_GLS_PSI,LYR_IDLAGE ,
     . LVL_BFSQ   ,LVL_DIFDIA ,LVL_DIFINT ,LVL_DIFISO ,LVL_DZ     ,
     . LVL_SALN   ,LVL_TEMP   ,LVL_TRC    ,LVL_UFLX   ,LVL_UTFLX  ,
     . LVL_USFLX  ,LVL_UMFLTD ,LVL_UTFLTD ,LVL_UTFLLD ,LVL_USFLTD ,
     . LVL_USFLLD ,LVL_UVEL   ,LVL_VFLX   ,LVL_VTFLX  ,LVL_VSFLX  ,
     . LVL_VMFLTD ,LVL_VTFLTD ,LVL_VTFLLD ,LVL_VSFLTD ,LVL_VSFLLD ,
     . LVL_VVEL   ,LVL_WFLX   ,LVL_WFLX2  ,LVL_PV     ,LVL_TKE    ,
     . LVL_GLS_PSI,LVL_IDLAGE ,
     . MSC_MMFLXL ,MSC_MMFLXD ,MSC_MMFTDL ,MSC_MMFTDD ,MSC_MHFLX  ,
     . MSC_MHFTD  ,MSC_MHFLD  ,MSC_MSFLX  ,MSC_MSFTD  ,MSC_MSFLD  ,
     . MSC_VOLTR  ,MSC_MASSGS ,MSC_VOLGS  ,MSC_SALNGA ,MSC_TEMPGA ,
     . MSC_SSSGA  ,MSC_SSTGA  ,
     . GLB_AVEPERIO,GLB_FILEFREQ,GLB_COMPFLAG,GLB_NCFORMAT
      character(len=10), dimension(nphymax), save :: GLB_FNAMETAG
      integer, dimension(nphymax), save ::
     . ACC_ABSWND ,ACC_ALB    ,ACC_BRNFLX ,ACC_BRNPD  ,ACC_DFL    ,
     . ACC_EVA    ,ACC_FICE   ,ACC_FMLTFZ ,ACC_HICE   ,ACC_HMLTFZ ,
     . ACC_HSNW   ,ACC_IAGE   ,ACC_IDKEDT ,ACC_LIP    ,ACC_MAXMLD ,
     . ACC_MLD    ,ACC_MLDU   ,ACC_MLDV   ,ACC_MLTS   ,ACC_MLTSMN ,
     . ACC_MLTSMX ,ACC_MLTSSQ ,ACC_MTKEUS ,ACC_MTKENI ,ACC_MTKEBF ,
     . ACC_MTKERS ,ACC_MTKEPE ,ACC_MTKEKE ,ACC_MTY    ,ACC_MXLU   ,
     . ACC_MXLV   ,ACC_NSF    ,ACC_PBOT   ,ACC_PSRF   ,ACC_RFIFLX ,
     . ACC_RNFFLX ,ACC_SALFLX ,ACC_SALRLX ,ACC_SBOT   ,ACC_SEALV  ,
     . ACC_SLVSQ  ,ACC_SFL    ,ACC_SOP    ,ACC_SIGMX  ,ACC_SSS    ,
     . ACC_SSSSQ  ,ACC_SST    ,ACC_SSTSQ  ,ACC_SURFLX ,ACC_SURRLX ,
     . ACC_SWA    ,ACC_T20D   ,ACC_TAUX   ,ACC_TAUY   ,ACC_TBOT   ,
     . ACC_TICE   ,ACC_TSRF   ,ACC_UB     ,ACC_UBFLXS ,ACC_UICE   ,
     . ACC_USTAR  ,ACC_USTAR3 ,ACC_VB     ,ACC_VBFLXS ,ACC_VICE   ,
     . ACC_ZTX    ,ACC_IVOLU  ,ACC_IVOLV  ,ACC_UTILH2D,
     . ACC_BFSQ   ,ACC_DIFDIA ,ACC_DIFINT ,ACC_DIFISO ,ACC_DP     ,
     . ACC_DPU    ,ACC_DPV    ,ACC_DZ     ,ACC_SALN   ,ACC_TEMP   ,
     . ACC_UFLX   ,ACC_UTFLX  ,ACC_USFLX  ,ACC_UMFLTD ,ACC_UTFLTD ,
     . ACC_UTFLLD ,ACC_USFLTD ,ACC_USFLLD ,ACC_UVEL   ,ACC_VFLX   ,
     . ACC_VTFLX  ,ACC_VSFLX  ,ACC_VMFLTD ,ACC_VTFLTD ,ACC_VTFLLD ,
     . ACC_VSFLTD ,ACC_VSFLLD ,ACC_VVEL   ,ACC_WFLX   ,ACC_WFLX2  ,
     . ACC_AVDSG  ,ACC_DPVOR  ,ACC_TKE    ,ACC_GLS_PSI,ACC_UTILLYR,
     . ACC_BFSQLVL  ,ACC_DIFDIALVL ,ACC_DIFINTLVL,ACC_DIFISOLVL,
     . ACC_DZLVL    ,ACC_SALNLVL   ,ACC_TEMPLVL  ,ACC_UFLXLVL  ,
     . ACC_UTFLXLVL ,ACC_USFLXLVL  ,ACC_UMFLTDLVL,ACC_UTFLTDLVL,
     . ACC_UTFLLDLVL,ACC_USFLTDLVL ,ACC_USFLLDLVL,ACC_UVELLVL  ,
     . ACC_VFLXLVL  ,ACC_VTFLXLVL  ,ACC_VSFLXLVL ,ACC_VMFLTDLVL,
     . ACC_VTFLTDLVL,ACC_VTFLLDLVL ,ACC_VSFLTDLVL,ACC_VSFLLDLVL,
     . ACC_VVELLVL  ,ACC_WFLXLVL   ,ACC_WFLX2LVL ,ACC_PVLVL    ,
     . ACC_TKELVL   ,ACC_GLS_PSILVL,ACC_UFLXOLD  ,ACC_VFLXOLD  ,
     . ACC_UTILLVL  ,
     . ACC_MMFLXL,ACC_MMFLXD,ACC_MMFTDL,ACC_MMFTDD,ACC_MHFLX,ACC_MHFTD,
     . ACC_MHFLD ,ACC_MSFLX ,ACC_MSFTD ,ACC_MSFLD ,ACC_VOLTR 
      namelist /MERDIA/
     . MER_ORFILE,MER_MIFILE,MER_REGNAM,MER_REGFLG,MER_MINLAT,MER_MAXLAT
      namelist /SECDIA/
     . SEC_SIFILE
      namelist /DIAPHY/
     . H2D_ABSWND ,H2D_ALB    ,H2D_BTMSTR ,H2D_BRNFLX ,H2D_BRNPD  ,
     . H2D_DFL    ,H2D_EVA    ,H2D_FICE   ,H2D_FMLTFZ ,H2D_HICE   ,
     . H2D_HMLTFZ ,H2D_HSNW   ,H2D_IAGE   ,H2D_IDKEDT ,H2D_LIP    ,
     . H2D_MAXMLD ,H2D_MLD    ,H2D_MLDU   ,H2D_MLDV   ,H2D_MLTS   ,
     . H2D_MLTSMN ,H2D_MLTSMX ,H2D_MLTSSQ ,H2D_MTKEUS ,H2D_MTKENI ,
     . H2D_MTKEBF ,H2D_MTKERS ,H2D_MTKEPE ,H2D_MTKEKE ,H2D_MTY    ,
     . H2D_MXLU   ,H2D_MXLV   ,H2D_NSF    ,H2D_PBOT   ,H2D_PSRF   ,
     . H2D_RFIFLX ,H2D_RNFFLX ,H2D_SALFLX ,H2D_SALRLX ,H2D_SBOT   ,
     . H2D_SEALV  ,H2D_SLVSQ  ,H2D_SFL    ,H2D_SOP    ,H2D_SIGMX  ,
     . H2D_SSS    ,H2D_SSSSQ  ,H2D_SST    ,H2D_SSTSQ  ,H2D_SURFLX ,
     . H2D_SURRLX ,H2D_SWA    ,H2d_T20D   ,H2D_TAUX   ,H2D_TAUY   ,
     . H2D_TBOT   ,H2D_TICE   ,H2D_TSRF   ,H2D_UB     ,H2D_UICE   ,
     . H2D_USTAR  ,H2D_USTAR3 ,H2D_VB     ,H2D_VICE   ,H2D_ZTX    ,
     . LYR_BFSQ   ,LYR_DIFDIA ,LYR_DIFINT ,LYR_DIFISO ,LYR_DP     ,
     . LYR_DPU    ,LYR_DPV    ,LYR_DZ     ,LYR_SALN   ,LYR_TEMP   ,
     . LYR_TRC    ,LYR_UFLX   ,LYR_UTFLX  ,LYR_USFLX  ,LYR_UMFLTD ,
     . LYR_UTFLTD ,LYR_UTFLLD ,LYR_USFLTD ,LYR_USFLLD ,LYR_UVEL   ,
     . LYR_VFLX   ,LYR_VTFLX  ,LYR_VSFLX  ,LYR_VMFLTD ,LYR_VTFLTD ,
     . LYR_VTFLLD ,LYR_VSFLTD ,LYR_VSFLLD ,LYR_VVEL   ,LYR_WFLX   ,
     . LYR_WFLX2  ,LYR_PV     ,LYR_TKE    ,LYR_GLS_PSI,LYR_IDLAGE ,
     . LVL_BFSQ   ,LVL_DIFDIA ,LVL_DIFINT ,LVL_DIFISO ,LVL_DZ     ,
     . LVL_SALN   ,LVL_TEMP   ,LVL_TRC    ,LVL_UFLX   ,LVL_UTFLX  ,
     . LVL_USFLX  ,LVL_UMFLTD ,LVL_UTFLTD ,LVL_UTFLLD ,LVL_USFLTD ,
     . LVL_USFLLD ,LVL_UVEL   ,LVL_VFLX   ,LVL_VTFLX  ,LVL_VSFLX  ,
     . LVL_VMFLTD ,LVL_VTFLTD ,LVL_VTFLLD ,LVL_VSFLTD ,LVL_VSFLLD ,
     . LVL_VVEL   ,LVL_WFLX   ,LVL_WFLX2  ,LVL_PV     ,LVL_TKE    ,
     . LVL_GLS_PSI,LVL_IDLAGE ,
     . MSC_MMFLXL ,MSC_MMFLXD ,MSC_MMFTDL ,MSC_MMFTDD ,MSC_MHFLX  ,
     . MSC_MHFTD  ,MSC_MHFLD  ,MSC_MSFLX  ,MSC_MSFTD  ,MSC_MSFLD  ,
     . MSC_VOLTR  ,MSC_MASSGS ,MSC_VOLGS  ,MSC_SALNGA ,MSC_TEMPGA ,
     . MSC_SSSGA  ,MSC_SSTGA  ,
     . GLB_AVEPERIO,GLB_FILEFREQ,GLB_COMPFLAG,GLB_NCFORMAT, 
     . GLB_FNAMETAG



      contains 



      subroutine diafnm(runid,runid_len,expcnf,ctag,nstep,diagfq,diagmon
     .  ,diagann,fname)
c
c --- ------------------------------------------------------------------
c --- Description: creates file name for the diagnostic output
c ---  
c --- Arguments: 
c ---   char runid     (in)    : acronym of simulation 
c ---   int  runid_len (in)    : actual length of runid 
c ---   char expcnf    (in)    : experiment configuration
c ---   char ctag      (in)    : string used in middle of file name 
c ---   real time      (in)    : BLOM time step 
c ---   real diagfq    (in)    : diagnostic frequency     
c ---   logi diagmon   (in)    : switch to show whether diagfq=month 
c ---   logi diagann   (in)    : switch to show whether diagfq=year 
c ---   char fname     (out)   : file name 
c --- ------------------------------------------
c
      implicit none
c
#include "common_clndr.h"
c
      character runid*(*),expcnf*(*),ctag*(*),fname*(*)
      integer runid_len,nstep
      real diagfq
      logical diagmon,diagann
c
      real epsil
      parameter (epsil=1.e-11)
c
      integer ny,nm,nd,ns
c
      if (diagfq+epsil.gt.1.) then
        ny=nyear
        nm=nmonth
        nd=nday
        if   (expcnf.eq.'ben02clim') then
          call pstdat(ny,nm,nd,nint(diagfq))
          if     (mod(nint(diagfq),360).eq.0) then
            write(fname,'(4a,i4.4,a)')
     .        runid(1:runid_len),'_',ctag,'_',ny,'.nc'
          elseif (mod(nint(diagfq),30).eq.0) then
            write(fname,'(4a,i4.4,a,i2.2,a)')
     .        runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.nc'
          else
            write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)')
     .        runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.',nd,'.nc'
          endif
        elseif (expcnf.eq.'ben02syn') then
          if (diagmon) then
            call pstdat(ny,nm,nd,1)
            write(fname,'(4a,i4.4,a,i2.2,a)')
     .        runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.nc'
          elseif (diagann) then
            call pstdat(ny,nm,nd,1)
            write(fname,'(4a,i4.4,a)')
     .        runid(1:runid_len),'_',ctag,'_',ny,'.nc'
          else
            call pstdat(ny,nm,nd,nint(diagfq))
            write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)')
     .        runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.',nd,'.nc'
          endif
        elseif (expcnf.eq.'cesm') then
          if (diagmon) then
            call pstdat(ny,nm,nd,1)
            write(fname,'(4a,i4.4,a,i2.2,a)')
     .        runid(1:runid_len),'.blom.',ctag,'.',ny,'-',nm,'.nc'
          elseif (diagann) then
            call pstdat(ny,nm,nd,1)
            write(fname,'(4a,i4.4,a)')
     .        runid(1:runid_len),'.blom.',ctag,'.',ny,'.nc'
          else
            call pstdat(ny,nm,nd,nint(diagfq))
            write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)')
     .        runid(1:runid_len),'.blom.',ctag,'.',ny,'-',nm,'-',nd,
     .        '.nc'
          endif
        else
          if (mnproc.eq.1) then
            write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
          endif
          call xcstop('(inifrc)')
                 stop '(inifrc)'
        endif
      else
        ny=nyear
        nm=nmonth
        nd=nday
        if (mod(nstep,nstep_in_day).eq.0) call pstdat(ny,nm,nd,1)
        ns=nint(nint(mod(nstep+.5-diagfq,real(nstep_in_day))/diagfq)
     .          *diagfq/nstep_in_day*86400)
        if (expcnf.eq.'cesm') then
          write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)')
     .      runid(1:runid_len),'.blom.',ctag,'.',ny,'-',nm,'-',nd,'-',
     .      ns,'.nc'
        else
          write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.',nd,'.',ns,
     .      '.nc'
        endif
      endif
      end subroutine diafnm



      subroutine diaini
c
c --- ------------------------------------------------------------------
c --- initialize diagnostic variables
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
c
      integer i,j,l,n,istat,istatsum
      integer, parameter :: imn=1-nbdy,imx=idm+nbdy,jmn=imn,jmx=jdm+nbdy
      logical fexist
c
c --- Check existence of data files for meridional and section transport
c --- diagnostics 
      if (mnproc.eq.1) then
        inquire(file=mer_orfile,exist=fexist)
        if (.not.fexist) then
          write (lp,'(3a)') ' Could not find file ',trim(mer_orfile),'!'
          call xchalt('(diaini)')
                 stop '(diaini)'
        endif
        inquire(file=mer_mifile,exist=fexist)
        if (.not.fexist) then
          write (lp,'(3a)') ' Could not find file ',trim(mer_mifile),'!'
          call xchalt('(diaini)')
                 stop '(diaini)'
        endif
        inquire(file=sec_sifile,exist=fexist)
        if (.not.fexist) then
          write (lp,'(3a)') ' Could not find file ',trim(sec_sifile),'!'
          call xchalt('(diaini)')
                 stop '(diaini)'
        endif
      endif
c
      nphyh2d=0
      nphylyr=0
      nphylvl=0
c
c --- Loop over io groups
      do n=1,nphy
        nacc_phy(n)=0
c
c --- - Solve dependencies for diagnostic variables (0=skipped) 
        ACC_ABSWND(n)   = H2D_ABSWND(n)
        ACC_ALB(n)      = H2D_ALB(n)
        ACC_BRNFLX(n)   = H2D_BRNFLX(n)
        ACC_BRNPD(n)    = H2D_BRNPD(n)
        ACC_DFL(n)      = H2D_DFL(n)
        ACC_EVA(n)      = H2D_EVA(n)
        ACC_FMLTFZ(n)   = H2D_FMLTFZ(n)
        ACC_FICE(n)     = H2D_FICE(n)   + H2D_HICE(n)   + H2D_UICE(n)  +
     .                    H2D_VICE(n)   + H2D_HSNW(n)
        ACC_HICE(n)     = H2D_HICE(n)   + H2D_UICE(n)   + H2D_VICE(n)
        ACC_HMLTFZ(n)   = H2D_HMLTFZ(n)
        ACC_HSNW(n)     = H2D_HSNW(n)
        ACC_IAGE(n)     = H2D_IAGE(n)
        ACC_IDKEDT(n)   = H2D_IDKEDT(n)
        ACC_IVOLU(n)    = H2D_UICE(n)
        ACC_IVOLV(n)    = H2D_VICE(n)
        ACC_LIP(n)      = H2D_LIP(n)
        ACC_MAXMLD(n)   = H2D_MAXMLD(n)
        ACC_MLD(n)      = H2D_MLD(n)
        ACC_MLDU(n)     = H2D_MLDU(n)   + H2D_MXLU(n)
        ACC_MLDV(n)     = H2D_MLDV(n)   + H2D_MXLV(n)
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
        ACC_MXLU(n)     = H2D_MXLU(n)
        ACC_MXLV(n)     = H2D_MXLV(n)
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
        ACC_SIGMX(n)    = H2D_SIGMX(n)  + MSC_MMFLXL(n) + MSC_MMFTDL(n)
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
        ACC_VB(n)       = H2D_VB(n)
        ACC_VBFLXS(n)   = H2D_BTMSTR(n)
        ACC_VICE(n)     = H2D_VICE(n)
        ACC_ZTX(n)      = H2D_ZTX(n)
        ACC_BFSQ(n)     = LYR_BFSQ(n)
        ACC_BFSQLVL(n)  = LVL_BFSQ(n)
        ACC_DIFDIA(n)   = LYR_DIFDIA(n)
        ACC_DIFDIALVL(n)= LVL_DIFDIA(n)
        ACC_DIFINT(n)   = LYR_DIFINT(n)
        ACC_DIFINTLVL(n)= LVL_DIFINT(n)
        ACC_DIFISO(n)   = LYR_DIFISO(n)
        ACC_DIFISOLVL(n)= LVL_DIFISO(n)
        ACC_DP(n)       = LYR_DP(n)     + LYR_BFSQ(n)   +
     .                    LYR_SALN(n)   + LYR_TEMP(n)   +
     .                    LYR_DIFDIA(n) + LYR_DIFINT(n) + LYR_DIFISO(n)+
     .                    LYR_TKE(n)    + LYR_GLS_PSI(n)+
     .                    LVL_BFSQ(n)   + LVL_SALN(n)   + LVL_TEMP(n)  +
     .                    LVL_DIFDIA(n) + LVL_DIFINT(n) + LVL_DIFISO(n)+
     .                    LVL_TKE(n)    + LVL_GLS_PSI(n)+
     .                    MSC_MASSGS(n) + MSC_SALNGA(n) + MSC_TEMPGA(n)
        ACC_DPU(n)      = LYR_DPU(n)    + LYR_UVEL(n)
        ACC_DPV(n)      = LYR_DPV(n)    + LYR_VVEL(n)
        ACC_DZ(n)       = LYR_DZ(n)     + MSC_VOLGS(n)
        ACC_DZLVL(n)    = LVL_DZ(n)
        ACC_SALN(n)     = LYR_SALN(n)   + MSC_SALNGA(n)
        ACC_SALNLVL(n)  = LVL_SALN(n)
        ACC_TEMP(n)     = LYR_TEMP(n)   + MSC_TEMPGA(n)
        ACC_TEMPLVL(n)  = LVL_TEMP(n)
        ACC_UFLX(n)     = LYR_UFLX(n)   + MSC_MMFLXL(n) + LYR_WFLX(n)  +
     .                    LYR_WFLX2(n)
        ACC_UFLXLVL(n)  = LVL_UFLX(n)   + MSC_MMFLXD(n) + MSC_VOLTR(n) +
     .                    LVL_WFLX(n)   + LVL_WFLX2(n)
        ACC_UFLXOLD(n)  = LVL_WFLX(n)   + LVL_WFLX2(n)
        ACC_UTFLX(n)    = LYR_UTFLX(n)  + MSC_MHFLX(n)
        ACC_UTFLXLVL(n) = LVL_UTFLX(n)
        ACC_USFLX(n)    = LYR_USFLX(n)  + MSC_MSFLX(n)
        ACC_USFLXLVL(n) = LVL_USFLX(n)
        ACC_UMFLTD(n)   = LYR_UMFLTD(n) + MSC_MMFTDL(n)
        ACC_UMFLTDLVL(n)= LVL_UMFLTD(n) + MSC_MMFTDD(n)
        ACC_UTFLTD(n)   = LYR_UTFLTD(n) + MSC_MHFTD(n)
        ACC_UTFLTDLVL(n)= LVL_UTFLTD(n)
        ACC_UTFLLD(n)   = LYR_UTFLLD(n) + MSC_MHFLD(n)
        ACC_UTFLLDLVL(n)= LVL_UTFLLD(n)
        ACC_USFLTD(n)   = LYR_USFLTD(n) + MSC_MSFTD(n)
        ACC_USFLTDLVL(n)= LVL_USFLTD(n)
        ACC_USFLLD(n)   = LYR_USFLLD(n) + MSC_MSFLD(n)
        ACC_USFLLDLVL(n)= LVL_USFLLD(n)
        ACC_UVEL(n)     = LYR_UVEL(n)
        ACC_UVELLVL(n)  = LVL_UVEL(n)
        ACC_VFLX(n)     = LYR_VFLX(n)   + MSC_MMFLXL(n) + LYR_WFLX(n)  +
     .                    LYR_WFLX2(n)
        ACC_VFLXLVL(n)  = LVL_VFLX(n)   + MSC_MMFLXD(n) + MSC_VOLTR(n) +
     .                    LVL_WFLX(n)   + LVL_WFLX2(n)
        ACC_VFLXOLD(n)  = LVL_WFLX(n)   + LVL_WFLX2(n)
        ACC_VTFLX(n)    = LYR_VTFLX(n)  + MSC_MHFLX(n)
        ACC_VTFLXLVL(n) = LVL_VTFLX(n)
        ACC_VSFLX(n)    = LYR_VSFLX(n)  + MSC_MSFLX(n)
        ACC_VSFLXLVL(n) = LVL_VSFLX(n)
        ACC_VMFLTD(n)   = LYR_VMFLTD(n) + MSC_MMFTDL(n)
        ACC_VMFLTDLVL(n)= LVL_VMFLTD(n) + MSC_MMFTDD(n)
        ACC_VTFLTD(n)   = LYR_VTFLTD(n) + MSC_MHFTD(n)
        ACC_VTFLTDLVL(n)= LVL_VTFLTD(n)
        ACC_VTFLLD(n)   = LYR_VTFLLD(n) + MSC_MHFLD(n)
        ACC_VTFLLDLVL(n)= LVL_VTFLLD(n)
        ACC_VSFLTD(n)   = LYR_VSFLTD(n) + MSC_MSFTD(n)
        ACC_VSFLTDLVL(n)= LVL_VSFLTD(n)
        ACC_VSFLLD(n)   = LYR_VSFLLD(n) + MSC_MSFLD(n)
        ACC_VSFLLDLVL(n)= LVL_VSFLLD(n)
        ACC_VVEL(n)     = LYR_VVEL(n)
        ACC_VVELLVL(n)  = LVL_VVEL(n)
        ACC_WFLX(n)     = LYR_WFLX(n)   + LYR_WFLX2(n)  + LVL_WFLX(n)  +
     .                    LVL_WFLX2(n)
        ACC_WFLXLVL(n)  = LVL_WFLX(n)   + LVL_WFLX2(n)  + LYR_WFLX(n)  +
     .                    LYR_WFLX2(n)
        ACC_WFLX2(n)    = LYR_WFLX2(n)  + LYR_WFLX(n)   + LVL_WFLX(n)  +
     .                    LVL_WFLX2(n)
        ACC_WFLX2LVL(n) = LVL_WFLX2(n)  + LVL_WFLX(n)   + LYR_WFLX(n)  +
     .                    LYR_WFLX2(n)
        ACC_AVDSG(n)    = LYR_PV(n)
        ACC_DPVOR(n)    = LYR_PV(n)
        ACC_PVLVL(n)    = LVL_PV(n)
        ACC_TKE(n)      = LYR_TKE(n)
        ACC_TKELVL(n)   = LVL_TKE(n)
        ACC_GLS_PSI(n)  = LYR_GLS_PSI(n)
        ACC_GLS_PSILVL(n) = LVL_GLS_PSI(n)
        ACC_MMFLXL(n)   = MSC_MMFLXL(n) 
        ACC_MMFLXD(n)   = MSC_MMFLXD(n)
        ACC_MMFTDL(n)   = MSC_MMFTDL(n) 
        ACC_MMFTDD(n)   = MSC_MMFTDD(n) 
        ACC_MHFLX(n)    = MSC_MHFLX(n)
        ACC_MHFTD(n)    = MSC_MHFTD(n)
        ACC_MHFLD(n)    = MSC_MHFLD(n)
        ACC_MSFLX(n)    = MSC_MSFLX(n)
        ACC_MSFTD(n)    = MSC_MSFTD(n)
        ACC_MSFLD(n)    = MSC_MSFLD(n)
        ACC_VOLTR(n)    = MSC_VOLTR(n)
c
c --- - Determine position in buffer 
        if (ACC_ABSWND(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_ABSWND(n)=nphyh2d*min(1,ACC_ABSWND(n))
        if (ACC_ALB(n).ne.0) nphyh2d=nphyh2d+1
         ACC_ALB(n)=nphyh2d*min(1,ACC_ALB(n)) 
        if (ACC_BRNFLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_BRNFLX(n)=nphyh2d*min(1,ACC_BRNFLX(n))
        if (ACC_BRNPD(n).ne.0) nphyh2d=nphyh2d+1
         ACC_BRNPD(n)=nphyh2d*min(1,ACC_BRNPD(n))
        if (ACC_DFL(n).ne.0) nphyh2d=nphyh2d+1
         ACC_DFL(n)=nphyh2d*min(1,ACC_DFL(n)) 
        if (ACC_EVA(n).ne.0) nphyh2d=nphyh2d+1
         ACC_EVA(n)=nphyh2d*min(1,ACC_EVA(n)) 
        if (ACC_FMLTFZ(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_FMLTFZ(n)=nphyh2d*min(1,ACC_FMLTFZ(n))
        if (ACC_FICE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_FICE(n)=nphyh2d*min(1,ACC_FICE(n)) 
        if (ACC_HICE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_HICE(n)=nphyh2d*min(1,ACC_HICE(n)) 
        if (ACC_HMLTFZ(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_HMLTFZ(n)=nphyh2d*min(1,ACC_HMLTFZ(n))
        if (ACC_HSNW(n).ne.0) nphyh2d=nphyh2d+1
         ACC_HSNW(n)=nphyh2d*min(1,ACC_HSNW(n)) 
        if (ACC_IAGE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_IAGE(n)=nphyh2d*min(1,ACC_IAGE(n)) 
        if (ACC_IDKEDT(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_IDKEDT(n)=nphyh2d*min(1,ACC_IDKEDT(n))
        if (ACC_IVOLU(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_IVOLU(n)=nphyh2d*min(1,ACC_IVOLU(n))
        if (ACC_IVOLV(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_IVOLV(n)=nphyh2d*min(1,ACC_IVOLV(n))
        if (ACC_LIP(n).ne.0) nphyh2d=nphyh2d+1
         ACC_LIP(n)=nphyh2d*min(1,ACC_LIP(n)) 
        if (ACC_MAXMLD(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_MAXMLD(n)=nphyh2d*min(1,ACC_MAXMLD(n))
        if (ACC_MLD(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLD(n)=nphyh2d*min(1,ACC_MLD(n)) 
        if (ACC_MLDU(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLDU(n)=nphyh2d*min(1,ACC_MLDU(n)) 
        if (ACC_MLDV(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLDV(n)=nphyh2d*min(1,ACC_MLDV(n))
        if (ACC_MLTS(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLTS(n)=nphyh2d*min(1,ACC_MLTS(n))
        if (ACC_MLTSMN(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLTSMN(n)=nphyh2d*min(1,ACC_MLTSMN(n))
        if (ACC_MLTSMX(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLTSMX(n)=nphyh2d*min(1,ACC_MLTSMX(n))
        if (ACC_MLTSSQ(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MLTSSQ(n)=nphyh2d*min(1,ACC_MLTSSQ(n))
        if (ACC_MTKEUS(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTKEUS(n)=nphyh2d*min(1,ACC_MTKEUS(n)) 
        if (ACC_MTKENI(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTKENI(n)=nphyh2d*min(1,ACC_MTKENI(n)) 
        if (ACC_MTKEBF(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTKEBF(n)=nphyh2d*min(1,ACC_MTKEBF(n)) 
        if (ACC_MTKERS(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTKERS(n)=nphyh2d*min(1,ACC_MTKERS(n)) 
        if (ACC_MTKEPE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTKEPE(n)=nphyh2d*min(1,ACC_MTKEPE(n)) 
        if (ACC_MTKEKE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTKEKE(n)=nphyh2d*min(1,ACC_MTKEKE(n)) 
        if (ACC_MTY(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MTY(n)=nphyh2d*min(1,ACC_MTY(n)) 
        if (ACC_MXLU(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MXLU(n)=nphyh2d*min(1,ACC_MXLU(n)) 
        if (ACC_MXLV(n).ne.0) nphyh2d=nphyh2d+1
         ACC_MXLV(n)=nphyh2d*min(1,ACC_MXLV(n))
        if (ACC_NSF(n).ne.0) nphyh2d=nphyh2d+1
         ACC_NSF(n)=nphyh2d*min(1,ACC_NSF(n)) 
        if (ACC_PBOT(n).ne.0) nphyh2d=nphyh2d+1
         ACC_PBOT(n)=nphyh2d*min(1,ACC_PBOT(n)) 
        if (ACC_PSRF(n).ne.0) nphyh2d=nphyh2d+1
         ACC_PSRF(n)=nphyh2d*min(1,ACC_PSRF(n)) 
        if (ACC_RFIFLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_RFIFLX(n)=nphyh2d*min(1,ACC_RFIFLX(n))
        if (ACC_RNFFLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_RNFFLX(n)=nphyh2d*min(1,ACC_RNFFLX(n))
        if (ACC_SURFLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SURFLX(n)=nphyh2d*min(1,ACC_SURFLX(n))
        if (ACC_SURRLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SURRLX(n)=nphyh2d*min(1,ACC_SURRLX(n))
        if (ACC_SALFLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SALFLX(n)=nphyh2d*min(1,ACC_SALFLX(n))
        if (ACC_SALRLX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SALRLX(n)=nphyh2d*min(1,ACC_SALRLX(n))
        if (ACC_SBOT(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SBOT(n)=nphyh2d*min(1,ACC_SBOT(n))
        if (ACC_SEALV(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SEALV(n)=nphyh2d*min(1,ACC_SEALV(n))
        if (ACC_SLVSQ(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SLVSQ(n)=nphyh2d*min(1,ACC_SLVSQ(n))
        if (ACC_SFL(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SFL(n)=nphyh2d*min(1,ACC_SFL(n)) 
        if (ACC_SIGMX(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_SIGMX(n)=nphyh2d*min(1,ACC_SIGMX(n))
        if (ACC_SOP(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SOP(n)=nphyh2d*min(1,ACC_SOP(n)) 
        if (ACC_SSS(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SSS(n)=nphyh2d*min(1,ACC_SSS(n)) 
        if (ACC_SSSSQ(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SSSSQ(n)=nphyh2d*min(1,ACC_SSSSQ(n)) 
        if (ACC_SST(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SST(n)=nphyh2d*min(1,ACC_SST(n)) 
        if (ACC_SSTSQ(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SSTSQ(n)=nphyh2d*min(1,ACC_SSTSQ(n)) 
        if (ACC_SWA(n).ne.0) nphyh2d=nphyh2d+1
         ACC_SWA(n)=nphyh2d*min(1,ACC_SWA(n)) 
        if (ACC_T20D(n).ne.0) nphyh2d=nphyh2d+1
         ACC_T20D(n)=nphyh2d*min(1,ACC_T20D(n)) 
        if (ACC_TAUX(n).ne.0) nphyh2d=nphyh2d+1
         ACC_TAUX(n)=nphyh2d*min(1,ACC_TAUX(n)) 
        if (ACC_TAUY(n).ne.0) nphyh2d=nphyh2d+1
         ACC_TAUY(n)=nphyh2d*min(1,ACC_TAUY(n)) 
        if (ACC_TBOT(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_TBOT(n)=nphyh2d*min(1,ACC_TBOT(n))
        if (ACC_TICE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_TICE(n)=nphyh2d*min(1,ACC_TICE(n)) 
        if (ACC_TSRF(n).ne.0) nphyh2d=nphyh2d+1
         ACC_TSRF(n)=nphyh2d*min(1,ACC_TSRF(n))
        if (ACC_UB(n).ne.0) nphyh2d=nphyh2d+1
         ACC_UB(n)=nphyh2d*min(1,ACC_UB(n)) 
        if (ACC_UBFLXS(n).ne.0) nphyh2d=nphyh2d+1
         ACC_UBFLXS(n)=nphyh2d*min(1,ACC_UBFLXS(n)) 
        if (ACC_UICE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_UICE(n)=nphyh2d*min(1,ACC_UICE(n)) 
        if (ACC_USTAR(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_USTAR(n)=nphyh2d*min(1,ACC_USTAR(n))
        if (ACC_USTAR3(n).ne.0) nphyh2d=nphyh2d+1 
         ACC_USTAR3(n)=nphyh2d*min(1,ACC_USTAR3(n))
        if (ACC_VB(n).ne.0) nphyh2d=nphyh2d+1
         ACC_VB(n)=nphyh2d*min(1,ACC_VB(n)) 
        if (ACC_VBFLXS(n).ne.0) nphyh2d=nphyh2d+1
         ACC_VBFLXS(n)=nphyh2d*min(1,ACC_VBFLXS(n)) 
        if (ACC_VICE(n).ne.0) nphyh2d=nphyh2d+1
         ACC_VICE(n)=nphyh2d*min(1,ACC_VICE(n))
        if (ACC_ZTX(n).ne.0) nphyh2d=nphyh2d+1
         ACC_ZTX(n)=nphyh2d*min(1,ACC_ZTX(n))
c
        if (ACC_BFSQ(n).ne.0) nphylyr=nphylyr+1
         ACC_BFSQ(n)=nphylyr*min(1,ACC_BFSQ(n))
        if (ACC_DIFDIA(n).ne.0) nphylyr=nphylyr+1
         ACC_DIFDIA(n)=nphylyr*min(1,ACC_DIFDIA(n))
        if (ACC_DIFINT(n).ne.0) nphylyr=nphylyr+1
         ACC_DIFINT(n)=nphylyr*min(1,ACC_DIFINT(n))
        if (ACC_DIFISO(n).ne.0) nphylyr=nphylyr+1
         ACC_DIFISO(n)=nphylyr*min(1,ACC_DIFISO(n))
        if (ACC_DP(n).ne.0) nphylyr=nphylyr+1
         ACC_DP(n)=nphylyr*min(1,ACC_DP(n))
        if (ACC_DPU(n).ne.0) nphylyr=nphylyr+1
         ACC_DPU(n)=nphylyr*min(1,ACC_DPU(n))
        if (ACC_DPV(n).ne.0) nphylyr=nphylyr+1
         ACC_DPV(n)=nphylyr*min(1,ACC_DPV(n))
        if (ACC_DZ(n).ne.0) nphylyr=nphylyr+1
         ACC_DZ(n)=nphylyr*min(1,ACC_DZ(n))
        if (ACC_SALN(n).ne.0) nphylyr=nphylyr+1
         ACC_SALN(n)=nphylyr*min(1,ACC_SALN(n))
        if (ACC_TEMP(n).ne.0) nphylyr=nphylyr+1
         ACC_TEMP(n)=nphylyr*min(1,ACC_TEMP(n))
        if (ACC_UFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_UFLX(n)=nphylyr*min(1,ACC_UFLX(n))
        if (ACC_UTFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_UTFLX(n)=nphylyr*min(1,ACC_UTFLX(n))
        if (ACC_USFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_USFLX(n)=nphylyr*min(1,ACC_USFLX(n))
        if (ACC_UMFLTD(n).ne.0) nphylyr=nphylyr+1
         ACC_UMFLTD(n)=nphylyr*min(1,ACC_UMFLTD(n))
        if (ACC_UTFLTD(n).ne.0) nphylyr=nphylyr+1
         ACC_UTFLTD(n)=nphylyr*min(1,ACC_UTFLTD(n))
        if (ACC_UTFLLD(n).ne.0) nphylyr=nphylyr+1
         ACC_UTFLLD(n)=nphylyr*min(1,ACC_UTFLLD(n))
        if (ACC_USFLTD(n).ne.0) nphylyr=nphylyr+1
         ACC_USFLTD(n)=nphylyr*min(1,ACC_USFLTD(n))
        if (ACC_USFLLD(n).ne.0) nphylyr=nphylyr+1
         ACC_USFLLD(n)=nphylyr*min(1,ACC_USFLLD(n))
        if (ACC_UVEL(n).ne.0) nphylyr=nphylyr+1
         ACC_UVEL(n)=nphylyr*min(1,ACC_UVEL(n))
        if (ACC_VFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_VFLX(n)=nphylyr*min(1,ACC_VFLX(n))
        if (ACC_VTFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_VTFLX(n)=nphylyr*min(1,ACC_VTFLX(n))
        if (ACC_VSFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_VSFLX(n)=nphylyr*min(1,ACC_VSFLX(n))
        if (ACC_VMFLTD(n).ne.0) nphylyr=nphylyr+1
         ACC_VMFLTD(n)=nphylyr*min(1,ACC_VMFLTD(n))
        if (ACC_VTFLTD(n).ne.0) nphylyr=nphylyr+1
         ACC_VTFLTD(n)=nphylyr*min(1,ACC_VTFLTD(n))
        if (ACC_VTFLLD(n).ne.0) nphylyr=nphylyr+1
         ACC_VTFLLD(n)=nphylyr*min(1,ACC_VTFLLD(n))
        if (ACC_VSFLTD(n).ne.0) nphylyr=nphylyr+1
         ACC_VSFLTD(n)=nphylyr*min(1,ACC_VSFLTD(n))
        if (ACC_VSFLLD(n).ne.0) nphylyr=nphylyr+1
         ACC_VSFLLD(n)=nphylyr*min(1,ACC_VSFLLD(n))
        if (ACC_VVEL(n).ne.0) nphylyr=nphylyr+1
         ACC_VVEL(n)=nphylyr*min(1,ACC_VVEL(n))
        if (ACC_WFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_WFLX(n)=nphylyr*min(1,ACC_WFLX(n))
        if (ACC_WFLX2(n).ne.0) nphylyr=nphylyr+1
         ACC_WFLX2(n)=nphylyr*min(1,ACC_WFLX2(n))
        if (ACC_AVDSG(n).ne.0) nphylyr=nphylyr+1
         ACC_AVDSG(n)=nphylyr*min(1,ACC_AVDSG(n))
        if (ACC_DPVOR(n).ne.0) nphylyr=nphylyr+1
         ACC_DPVOR(n)=nphylyr*min(1,ACC_DPVOR(n))
        if (ACC_TKE(n).ne.0) nphylyr=nphylyr+1
         ACC_TKE(n)=nphylyr*min(1,ACC_TKE(n))
        if (ACC_GLS_PSI(n).ne.0) nphylyr=nphylyr+1
         ACC_GLS_PSI(n)=nphylyr*min(1,ACC_GLS_PSI(n))
c
        if (ACC_BFSQLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_BFSQLVL(n)=nphylvl*min(1,ACC_BFSQLVL(n))
        if (ACC_DIFDIALVL(n).ne.0) nphylvl=nphylvl+1
         ACC_DIFDIALVL(n)=nphylvl*min(1,ACC_DIFDIALVL(n))
        if (ACC_DIFINTLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_DIFINTLVL(n)=nphylvl*min(1,ACC_DIFINTLVL(n))
        if (ACC_DIFISOLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_DIFISOLVL(n)=nphylvl*min(1,ACC_DIFISOLVL(n))
        if (ACC_DZLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_DZLVL(n)=nphylvl*min(1,ACC_DZLVL(n))
        if (ACC_SALNLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_SALNLVL(n)=nphylvl*min(1,ACC_SALNLVL(n))
        if (ACC_TEMPLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_TEMPLVL(n)=nphylvl*min(1,ACC_TEMPLVL(n))
        if (ACC_UFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_UFLXLVL(n)=nphylvl*min(1,ACC_UFLXLVL(n))
        if (ACC_UFLXOLD(n).ne.0) nphylvl=nphylvl+1
         ACC_UFLXOLD(n)=nphylvl*min(1,ACC_UFLXOLD(n))
        if (ACC_UTFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_UTFLXLVL(n)=nphylvl*min(1,ACC_UTFLXLVL(n))
        if (ACC_USFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_USFLXLVL(n)=nphylvl*min(1,ACC_USFLXLVL(n))
        if (ACC_UMFLTDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_UMFLTDLVL(n)=nphylvl*min(1,ACC_UMFLTDLVL(n))
        if (ACC_UTFLTDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_UTFLTDLVL(n)=nphylvl*min(1,ACC_UTFLTDLVL(n))
        if (ACC_UTFLLDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_UTFLLDLVL(n)=nphylvl*min(1,ACC_UTFLLDLVL(n))
        if (ACC_USFLTDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_USFLTDLVL(n)=nphylvl*min(1,ACC_USFLTDLVL(n))
        if (ACC_USFLLDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_USFLLDLVL(n)=nphylvl*min(1,ACC_USFLLDLVL(n))
        if (ACC_UVELLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_UVELLVL(n)=nphylvl*min(1,ACC_UVELLVL(n))
        if (ACC_VFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VFLXLVL(n)=nphylvl*min(1,ACC_VFLXLVL(n))
        if (ACC_VFLXOLD(n).ne.0) nphylvl=nphylvl+1
         ACC_VFLXOLD(n)=nphylvl*min(1,ACC_VFLXOLD(n))
        if (ACC_VTFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VTFLXLVL(n)=nphylvl*min(1,ACC_VTFLXLVL(n))
        if (ACC_VSFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VSFLXLVL(n)=nphylvl*min(1,ACC_VSFLXLVL(n))
        if (ACC_VMFLTDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VMFLTDLVL(n)=nphylvl*min(1,ACC_VMFLTDLVL(n))
        if (ACC_VTFLTDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VTFLTDLVL(n)=nphylvl*min(1,ACC_VTFLTDLVL(n))
        if (ACC_VTFLLDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VTFLLDLVL(n)=nphylvl*min(1,ACC_VTFLLDLVL(n))
        if (ACC_VSFLTDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VSFLTDLVL(n)=nphylvl*min(1,ACC_VSFLTDLVL(n))
        if (ACC_VSFLLDLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VSFLLDLVL(n)=nphylvl*min(1,ACC_VSFLLDLVL(n))
        if (ACC_VVELLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_VVELLVL(n)=nphylvl*min(1,ACC_VVELLVL(n))
        if (ACC_WFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_WFLXLVL(n)=nphylvl*min(1,ACC_WFLXLVL(n))
        if (ACC_WFLX2LVL(n).ne.0) nphylvl=nphylvl+1
         ACC_WFLX2LVL(n)=nphylvl*min(1,ACC_WFLX2LVL(n))
        if (ACC_PVLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_PVLVL(n)=nphylvl*min(1,ACC_PVLVL(n))
        if (ACC_TKELVL(n).ne.0) nphylvl=nphylvl+1
         ACC_TKELVL(n)=nphylvl*min(1,ACC_TKELVL(n))
        if (ACC_GLS_PSILVL(n).ne.0) nphylvl=nphylvl+1
         ACC_GLS_PSILVL(n)=nphylvl*min(1,ACC_GLS_PSILVL(n))
c
c --- End loop over io groups
      enddo
c 
c --- Assign buffer positions for utility fields
      ACC_UTILH2D=0 
      nphyh2d=nphyh2d+1
      ACC_UTILH2D(1)=nphyh2d 
c
      ACC_UTILLYR=0 
      nphylyr=nphylyr+1
      ACC_UTILLYR(1)=nphylyr 
c
      ACC_UTILLVL=0 
      nphylvl=nphylvl+1
      ACC_UTILLVL(1)=nphylvl 
c 
c --- Allocate buffers
      istatsum=0
      istat=0
      if (nphyh2d.ne.0)
     .  allocate(phyh2d(imn:imx,jmn:jmx,nphyh2d),stat=istat)
      istatsum=istatsum+istat
      if (nphylyr.ne.0)
     .  allocate(phylyr(imn:imx,jmn:jmx,kdm,nphylyr),stat=istat)
      istatsum=istatsum+istat
      if (nphylvl.ne.0)
     .  allocate(phylvl(imn:imx,jmn:jmx,ddm,nphylvl),stat=istat)
      istatsum=istatsum+istat
      if (istatsum.ne.0) then
        write (lp,*) 'Cannot allocate enough memory!'
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
c
c --- initialisation of h2h, lyr and lvl fields
      do n=1,nphy
        call inifld(n)
      enddo
c
c --- Load bathymetry into module mod_dia (used for vertical 
c --- interpolation in BLOM and HAMOCC)
      nstepinday=nstep_in_day
c$OMP PARALLEL DO PRIVATE(l,i) 
      do j=1,jj+1
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pbath(i,j)=depths(i,j)
          enddo
        enddo
        do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            ubath(i,j)=min(depths(i,j),depths(i-1,j))
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vbath(i,j)=min(depths(i,j),depths(i,j-1))
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end subroutine diaini


      subroutine diasg1
c
c --- ------------------------------------------------------------------
c --- Extract reference potential density vector representative of the
c --- dominating ocean domain
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
c
      integer :: i,j,k,i1,j1
      logical lsigmar1
      real, dimension (itdm,jtdm) :: tmp2d
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          util1(i,j)=real(ipwocn(i,j))
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(tmp2d,util1,1)
      if (mnproc.eq.1) then
        lsigmar1=.false.
        do j=1,jtdm
          do i=1,itdm
            if (tmp2d(i,j).gt.0.) then
              i1=i
              j1=j 
              lsigmar1=.true.
              exit
            endif
          enddo
          if (lsigmar1) exit
        enddo
      endif
      call xcbcst(i1)
      call xcbcst(j1)
      do k=1,kk
        call xceget(sigmar1(k),sigmar(1-nbdy,1-nbdy,k),i1,j1)
        sigmar1(k)=sigmar1(k)*1.e3 ! Convert units from g cm-3 to kg m-3
      enddo
      if (mnproc.eq.1) then
        write(lp,*) 'sigma layers=',sigmar1 
      endif
c
      end subroutine diasg1



      subroutine diaacc(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- accumulate diagnostic variables
c --- ------------------------------------------------------------------
c
      use mod_thdysi, only: tsrfm,ticem
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_seaice.h"
#if defined(TRC) && defined(TKE)
#  include "param_trc.h"
#  include "common_trc.h"
#endif
c
      integer i,j,k,l,km,kup,iogrp
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2,
     .  ipsw,ipse,ipnw,ipne
c     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts,
     .  wghtsflx
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: uvel,vvel,
     .  dz,avdsg_p,dpvor_p,pv_p,dummy
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     .  dpml,sbot,tbot,dps,mlts,t20d
      real dsig,q,zup,zlo,plo,dbup,dblo,tup,tlo
c
      real p_alpha,rho
      external p_alpha,rho
c
c --- Increase counter 
      do iogrp=1,nphy
        nacc_phy(iogrp)=nacc_phy(iogrp)+1
      enddo
c
c --- Define auxillary variables 
c
      if (sum(ACC_UICE(1:nphy)+ACC_VICE(1:nphy)).ne.0) then
        call xctilr(hicem, 1,1, 1,1, halo_ps)
        call xctilr(ficem, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              util1(i,j)=hicem(i-1,j)*ficem(i-1,j)+hicem(i,j)*
     .          ficem(i,j)
            enddo
          enddo
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              util3(i,j)=hicem(i,j-1)*ficem(i,j-1)+hicem(i,j)*
     .          ficem(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_MXLU(1:nphy)+ACC_MXLV(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              util2(i,j)=(u(i,j,k1m)+ub(i,j,m))*dpu(i,j,k1m)
            enddo
          enddo
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              util4(i,j)=(v(i,j,k1m)+vb(i,j,m))*dpv(i,j,k1m)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_MLD(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              dpml(i,j)=dp(i,j,1+mm)+dp(i,j,2+mm)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_UVEL(1:nphy)+ACC_UVELLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
        do j=1,jj
          do k=1,kk
            km=k+mm
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                uvel(i,j,k)=u(i,j,km)+ub(i,j,m)
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_VVEL(1:nphy)+ACC_VVELLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
        do j=1,jj+1
          do k=1,kk
            km=k+mm
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vvel(i,j,k)=v(i,j,km)+vb(i,j,m)
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_DZ(1:nphy)+ACC_DZLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            z(i,j,kk+1)=-phi(i,j,kk+1)/g
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
       do j=1,jj
         do k=kk,1,-1
           km=k+mm
           do l=1,isp(j)
           do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
             if (dp(i,j,km).lt.epsil) then
               z(i,j,k)=z(i,j,k+1)
             else
               z(i,j,k)=z(i,j,k+1)+p_alpha(p(i,j,k+1),p(i,j,k),
     .                                     temp(i,j,km),saln(i,j,km))/g
             endif
             dz(i,j,k)=z(i,j,k+1)-z(i,j,k)
           enddo
           enddo
         enddo
       enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_AVDSG(1:nphy)+ACC_PVLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i,k,km,dsig)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              ipsw(i,j)=min(1,iu(i  ,j)+iv(i,j  ))
              ipse(i,j)=min(1,iu(i+1,j)+iv(i,j  ))
              ipnw(i,j)=min(1,iu(i  ,j)+iv(i,j+1))
              ipne(i,j)=min(1,iu(i+1,j)+iv(i,j+1))
            enddo
          enddo
          do k=1,kk
            km=k+mm
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (k.le.2) then
                  dsig=max(0.,sigma(i,j,2+mm)-sigma(i,j,1+mm))
                elseif (k.eq.kk) then
                  dsig=max(0.,sigma(i,j,km)
     .                       -max(sigma(i,j,km-1),sigma(i,j,2+mm)))
                else
                  dsig=.5*max(0.,sigma(i,j,km+1)
     .                          -max(sigma(i,j,km-1),sigma(i,j,2+mm)))
                endif
                avdsg_p(i,j,k)=(absvor(i  ,j  ,k)*ipsw(i,j)
     .                         +absvor(i+1,j  ,k)*ipse(i,j)
     .                         +absvor(i  ,j+1,k)*ipnw(i,j)
     .                         +absvor(i+1,j+1,k)*ipne(i,j))*dsig
                dpvor_p(i,j,k)=dpvor(i  ,j  ,k)*ipsw(i,j)
     .                        +dpvor(i+1,j  ,k)*ipse(i,j)
     .                        +dpvor(i  ,j+1,k)*ipnw(i,j)
     .                        +dpvor(i+1,j+1,k)*ipne(i,j)
                pv_p(i,j,k)=avdsg_p(i,j,k)/dpvor_p(i,j,k)
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_SBOT(1:nphy)+ACC_TBOT(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              sbot(i,j)=0.
              tbot(i,j)=0.
              dps(i,j)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,km,l,i,q)
        do j=1,jj
          do k=1,kk
            km=k+mm
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                q=max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpbot,p(i,j,k)))
                sbot(i,j)=sbot(i,j)+saln(i,j,km)*q
                tbot(i,j)=tbot(i,j)+temp(i,j,km)*q
                dps(i,j)=dps(i,j)+q
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (dps(i,j).gt.onemm) then
                q=1./dps(i,j)
                sbot(i,j)=sbot(i,j)*q
                tbot(i,j)=tbot(i,j)*q
              else
                sbot(i,j)=saln(i,j,1+mm)
                tbot(i,j)=temp(i,j,1+mm)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_MLTS(1:nphy)+ACC_MLTSMN(1:nphy)
     .       +ACC_MLTSMX(1:nphy)+ACC_MLTSSQ(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i,k,km,zup,dbup,plo,zlo,dblo)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              k=2
              km=k+mm
              zup=z(i,j,1)+.5*dz(i,j,1) ! .5*(z(i,j,1)+z(i,j,2))
              dbup=0.
              do
                if (dp(i,j,km).gt.onecm) then
                  plo=p(i,j,k)+.5*dp(i,j,km) ! .5*(p(i,j,k)+p(i,j,k+1))
                  zlo=z(i,j,k)+.5*dz(i,j,k ) ! .5*(z(i,j,k)+z(i,j,k+1))
                  dblo=g*(1.-rho(plo,temp(i,j,k1m),saln(i,j,k1m))
     .                      /rho(plo,temp(i,j,km ),saln(i,j,km )))
                  if (dblo.le.dbcrit) then
                    zup=zlo
                    dbup=dblo
                  else
                    dbup=min(dbup,dbcrit-epsil)
                    mlts(i,j)=(zup*(dblo-dbcrit)
     .                        +zlo*(dbcrit-dbup))/(dblo-dbup)-z(i,j,1)
                    exit
                  endif
                endif
                k=k+1
                if (k.gt.kk) then
                  mlts(i,j)=z(i,j,kk+1)-z(i,j,1)
                  exit
                endif
                km=k+mm
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_T20D(1:nphy)).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i,k,km,kup,zup,zlo,tup,tlo)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              k=1
              km=k+mm
              do 
                if (dp(i,j,km).gt.onecm) then
                  if (temp(i,j,km).gt.20.) then
                    kup=k
                  else
                    exit
                  endif
                endif
                k=k+1
                km=k+mm
                if (k.gt.kk) exit
              enddo
              if     (k.eq.1) then
                t20d(i,j)=0.
              elseif (k.gt.kk) then
                t20d(i,j)=z(i,j,kk+1)
              else
                zup=z(i,j,kup)+.5*dz(i,j,kup)
                zlo=z(i,j,k  )+.5*dz(i,j,k  )
                tup=temp(i,j,kup+mm)
                tlo=min(temp(i,j,km),tup-epsil)
                t20d(i,j)=(zup*(tlo-20.)+zlo*(20.-tup))/(tlo-tup)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- ------------------------------------------------------------------
c --- accumulate 2d diagnostic variables
c --- ------------------------------------------------------------------
c
c --- u-component of barotropic velocity [cm/s]
      call acch2d(ACC_UB,ub(1-nbdy,1-nbdy,m),dummy,0,'u')
c
c --- u-component of barotropic mass flux [g*cm/s^3]
      call acch2d(ACC_UBFLXS,ubflxs(1-nbdy,1-nbdy,n),dummy,0,'u')
c
c --- u-component of wind stress [N/m^2]
      call acch2d(ACC_ZTX,ztx,dummy,0,'u')
c
c --- u-component of momentum flux received by the ocean [dyn/cm^2]
      call acch2d(ACC_TAUX,taux,dummy,0,'u')
c
c --- weighted u-component of ice velocity [m^2/s]
      call acch2d(ACC_UICE,uicem,util1,1,'u')
c
c --- weighted u-component of total velocity [g/s^3]
      call acch2d(ACC_MXLU,util2,dummy,0,'u')
c
c --- mixed layer pressure thickness at u-point [g/cm/s^2]
      call acch2d(ACC_MLDU,dpu(1-nbdy,1-nbdy,k1m),dummy,0,'u')
c
c --- v-component of barotropic velocity [cm/s]
      call acch2d(ACC_VB,vb(1-nbdy,1-nbdy,m),dummy,0,'v')
c
c --- v-component of barotropic mass flux [g*cm/s^3]
      call acch2d(ACC_VBFLXS,vbflxs(1-nbdy,1-nbdy,n),dummy,0,'v')
c
c --- v-component of wind stress [N/m^2]
      call acch2d(ACC_MTY,mty,dummy,0,'v')
c
c --- v-component of momentum flux received by the ocean [dyn/cm^2]
      call acch2d(ACC_TAUY,tauy,dummy,0,'v')
c
c --- weighted v-component of ice velocity [m^2/s]
      call acch2d(ACC_VICE,vicem,util3,1,'v')
c
c --- weighted v-component of total velocity [g/s^3]
      call acch2d(ACC_MXLV,util4,dummy,0,'v')
c
c --- mixed layer pressure thickness at v-point [g/cm/s^2]
      call acch2d(ACC_MLDV,dpv(1-nbdy,1-nbdy,k1m),dummy,0,'v')
c
c --- surface pressure [g/cm/s^2]
      call acch2d(ACC_PSRF,p(1-nbdy,1-nbdy,1),dummy,0,'p')
c
c --- bottom pressure [g/cm/s^2]
      call acch2d(ACC_PBOT,p(1-nbdy,1-nbdy,kk+1),dummy,0,'p')
c
c --- sea level height [cm]
      call acch2d(ACC_SEALV,z(1-nbdy,1-nbdy,1),dummy,0,'p')
c
c --- mixed layer density (sigma units)
      call acch2d(ACC_SIGMX,sigma(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- weighted ice thickness [m]
      call acch2d(ACC_HICE,hicem,ficem,1,'p')
c
c --- weighted snow thickness [m]
      call acch2d(ACC_HSNW,hsnwm,ficem,1,'p')
c
c --- fractional ice cover
      call acch2d(ACC_FICE,ficem,dummy,0,'p')
c 
c --- ice volume in u-points[m] 
      call acch2d(ACC_IVOLU,util1,dummy,0,'u')
c 
c --- ice volume in v-points[m] 
      call acch2d(ACC_IVOLV,util3,dummy,0,'v')
c
c --- surface temperature [K]
      call acch2d(ACC_TSRF,tsrfm,dummy,0,'p')
c
c --- ice temperature [K]
      call acch2d(ACC_TICE,ticem,dummy,0,'p')
c
c --- short wave heat flux [W/m^2]
      call acch2d(ACC_SWA,swa,dummy,0,'p')
c
c --- non-solar heat flux [W/m^2]
      call acch2d(ACC_NSF,nsf,dummy,0,'p')
c
c --- heat flux due to melting/freezing [W/m^2]
      call acch2d(ACC_HMLTFZ,hmltfz,dummy,0,'p')
c
c --- derivative of non-solar heat flux by surface temperature [W/m^2/K]
      call acch2d(ACC_DFL,dfl,dummy,0,'p')
c
c --- liquid precipitation [mm/day]
      call acch2d(ACC_LIP,lip,dummy,0,'p')
c
c --- solid precipitation [mm/day]
      call acch2d(ACC_SOP,sop,dummy,0,'p')
c
c --- evaporation [mm/day]
      call acch2d(ACC_EVA,eva,dummy,0,'p')
c
c --- fresh water flux due to melting/freezing [kg/m^2/s]
      call acch2d(ACC_FMLTFZ,fmltfz,dummy,0,'p')
c
c --- salt flux [kg/m^2/s]
      call acch2d(ACC_SFL,sfl,dummy,0,'p')
c
c --- albedo
      call acch2d(ACC_ALB,alb,dummy,0,'p')
c
c --- liquid runoff [kg m-2 s-1]
      call acch2d(ACC_RNFFLX,rnfflx,dummy,0,'p')
c
c --- frozen runoff [kg m-2 s-1]
      call acch2d(ACC_RFIFLX,rfiflx,dummy,0,'p')
c
c --- friction velocity [cm s-1]
      call acch2d(ACC_USTAR,ustar,dummy,0,'p')
c
c --- friction velocity cubed [cm3 s-3]
      call acch2d(ACC_USTAR3,ustar3,dummy,0,'p')
c
c --- mixed layer integrated inertial kinetic energy tendency [cm3 s-3]
      call acch2d(ACC_IDKEDT,idkedt,dummy,0,'p')
c
c --- absolute wind speed [m s-1]
      call acch2d(ACC_ABSWND,abswnd,dummy,0,'p')
c
c --- mixed layer TKE tendency related to friction velocity [cm3 s-3]
      call acch2d(ACC_MTKEUS,mtkeus,dummy,0,'p')
c
c --- mixed layer TKE tendency related to near inertial motions [cm3 s-3]
      call acch2d(ACC_MTKENI,mtkeni,dummy,0,'p')
c
c --- mixed layer TKE tendency related to buoyancy forcing [cm3 s-3]
      call acch2d(ACC_MTKEBF,mtkebf,dummy,0,'p')
c
c --- mixed layer TKE tendency related to eddy restratification [cm3 s-3]
      call acch2d(ACC_MTKERS,mtkers,dummy,0,'p')
c
c --- mixed layer TKE tendency related to potential energy change [cm3 s-3]
      call acch2d(ACC_MTKEPE,mtkepe,dummy,0,'p')
c
c --- mixed layer TKE tendency related to kinetic energy change [cm3 s-3]
      call acch2d(ACC_MTKEKE,mtkeke,dummy,0,'p')
c
c --- sea surface salinity [g kg-1]
      call acch2d(ACC_SSS,saln(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- sea surface temperature [degC]
      call acch2d(ACC_SST,temp(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- bottom salinity [g kg-1]
      call acch2d(ACC_SBOT,sbot,dummy,0,'p')
c
c --- bottom temperature [g kg-1]
      call acch2d(ACC_TBOT,tbot,dummy,0,'p')
c
c --- mixed layer pressure thickness [g/cm/s^2]
      call acch2d(ACC_MLD,dpml,dummy,0,'p')
c
c --- mixed layer thickness using "sigma-t" criterion [cm]
      call acch2d(ACC_MLTS,mlts,dummy,0,'p')
c
c --- 20C isoterm depth [cm]
      call acch2d(ACC_T20D,t20d,dummy,0,'p')
c
c --- heat flux given by the ocean [W/m^2]
      call acch2d(ACC_SURFLX,surflx,dummy,0,'p')
c
c --- salt flux given by the ocean [10^-3 g/cm^2/s]
      call acch2d(ACC_SALFLX,salflx,dummy,0,'p')
c
c --- restoring heat flux received by the ocean [W/m^2]
      call acch2d(ACC_SURRLX,surrlx,dummy,0,'p')
c
c --- restoring salt flux received by the ocean [10^-3 g/cm^2/s]
      call acch2d(ACC_SALRLX,salrlx,dummy,0,'p')
c
c --- brine flux received by the ocean [10^-3 g/cm^2/s]
      call acch2d(ACC_BRNFLX,brnflx,dummy,0,'p')
c
c --- brine plume pressure depth [g/cm/s^2]
      call acch2d(ACC_BRNPD,pbrnda,dummy,0,'p')
c
c --- ------------------------------------------------------------------
c --- store minimum or maximum of 2d diagnostic variables
c --- ------------------------------------------------------------------
c
c --- maximum mixed layer pressure thickness [g/cm/s^2]
      call maxh2d(ACC_MAXMLD,dpml,'p')
c
c --- minimum mixed layer thickness using "sigma-t" criterion [cm]
      call minh2d(ACC_MLTSMN,mlts,'p')
c
c --- maximum mixed layer thickness using "sigma-t" criterion [cm]
      call maxh2d(ACC_MLTSMX,mlts,'p')
c
c --- ------------------------------------------------------------------
c --- store squared of 2d diagnostic variables
c --- ------------------------------------------------------------------
c
c --- mixed layer thickness squared using "sigma-t" criterion [cm^2]
      call sqh2d(ACC_MLTSSQ,mlts,'p')
c
c --- sea level height squared [cm^2]
      call sqh2d(ACC_SLVSQ,z(1-nbdy,1-nbdy,1),'p')
c
c --- sea surface salinity squared [g2 kg-2]
      call sqh2d(ACC_SSSSQ,saln(1-nbdy,1-nbdy,k1m),'p')
c
c --- sea surface temperature squared [degC2]
      call sqh2d(ACC_SSTSQ,temp(1-nbdy,1-nbdy,k1m),'p')
c
c --- ------------------------------------------------------------------
c --- accumulate 3d diagnostic variables
c --- ------------------------------------------------------------------
c
c --- weighted u-component of total velocity [g/s^3]
      call acclyr(ACC_UVEL,uvel,dpu(1-nbdy,1-nbdy,k1m),1,'u')
c
c --- layer pressure thickness at u-point [g/cm/s^2]
      call acclyr(ACC_DPU,dpu(1-nbdy,1-nbdy,k1m),dummy,0,'u')
c
c --- u-component of mass flux [g*cm/s^2]
      call acclyr(ACC_UFLX,uflx(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of heat flux [K*g*cm/s^2]
      call acclyr(ACC_UTFLX,utflx(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of salt flux [g^2*cm/kg/s^2]
      call acclyr(ACC_USFLX,usflx(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of mass flux due to thickness diffusion [g*cm/s^2]
      call acclyr(ACC_UMFLTD,umfltd(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of heat flux due to thickness diffusion [K*g*cm/s^2]
      call acclyr(ACC_UTFLTD,utfltd(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of heat flux due to lateral diffusion [K*g*cm/s^2]
      call acclyr(ACC_UTFLLD,utflld(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of salt flux due to thickness diffusion [g^2*cm/kg/s^2]
      call acclyr(ACC_USFLTD,usfltd(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- u-component of salt flux due to lateral diffusion [g^2*cm/kg/s^2]
      call acclyr(ACC_USFLLD,usflld(1-nbdy,1-nbdy,k1n),dummy,0,'u')
c
c --- weighted v-component of total velocity [g/s^3]
      call acclyr(ACC_VVEL,vvel,dpv(1-nbdy,1-nbdy,k1m),1,'v')
c
c --- layer pressure thickness at v-point [g/cm/s^2]
      call acclyr(ACC_DPV,dpv(1-nbdy,1-nbdy,k1m),dummy,0,'v')
c
c --- v-component of mass flux [g*cm/s^2]
      call acclyr(ACC_VFLX,vflx(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of heat flux [K*g*cm/s^2]
      call acclyr(ACC_VTFLX,vtflx(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of salt flux [g^2*cm/kg/s^2]
      call acclyr(ACC_VSFLX,vsflx(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of mass flux due to thickness diffusion [g*cm/s^2]
      call acclyr(ACC_VMFLTD,vmfltd(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of heat flux due to thickness diffusion [K*g*cm/s^2]
      call acclyr(ACC_VTFLTD,vtfltd(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of heat flux due to lateral diffusion [K*g*cm/s^2]
      call acclyr(ACC_VTFLLD,vtflld(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of salt flux due to thickness diffusion [g^2*cm/kg/s^2]
      call acclyr(ACC_VSFLTD,vsfltd(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- v-component of salt flux due to lateral diffusion [g^2*cm/kg/s^2]
      call acclyr(ACC_VSFLLD,vsflld(1-nbdy,1-nbdy,k1n),dummy,0,'v')
c
c --- weighted salinity [g^2/kg/cm/s^2]
      call acclyr(ACC_SALN,saln(1-nbdy,1-nbdy,k1m),
     .            dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- weighted temperature [degC*g/cm/s^2]
      call acclyr(ACC_TEMP,temp(1-nbdy,1-nbdy,k1m),
     .            dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- layer pressure thickness [g/cm/s^2]
      call acclyr(ACC_DP,dp(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- layer geopotential thickness [cm^2/s^2]
      call acclyr(ACC_DZ,dz,dummy,0,'p')
c
c --- buoyancy frequency squared [1/s]
      call acclyr(ACC_BFSQ,bfsql,dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- layer interface diffusivity [cm^2/s]
      call acclyr(ACC_DIFINT,difint,dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- isopycnal diffusivity [cm^2/s]
      call acclyr(ACC_DIFISO,difiso,dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- diapycnal diffusivity [cm^2/s]
      call acclyr(ACC_DIFDIA,difdia,dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- absolute vorticity multiplied with potential density difference
c --- over layer [g/cm^3/s]
      call acclyr(ACC_AVDSG,avdsg_p,dummy,0,'p')
c
c --- layer pressure thickness used in vorticity computation [g/cm/s^2]
      call acclyr(ACC_DPVOR,dpvor_p,dummy,0,'p')
c
#if defined(TRC) && defined(TKE)
c --- weighted tke [cm2/s2*g/cm/s^2]
      call acclyr(ACC_TKE,trc(1-nbdy,1-nbdy,k1m,itrtke),
     .            dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- weighted gls_psi [cm2/s3*g/cm/s^2]
      call acclyr(ACC_GLS_PSI,trc(1-nbdy,1-nbdy,k1m,itrgls),
     .            dp(1-nbdy,1-nbdy,k1m),1,'p')
c
#endif
c --- ------------------------------------------------------------------
c --- accumulate 3d diagnostic variables on Levitus levels 
c --- ------------------------------------------------------------------
c
      do iogrp=1,nphy 
        if (ACC_WFLXLVL(iogrp)+ACC_WFLX2LVL(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(k,l,i)
          do j=1,jj+1
            do k=1,ddm
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l)) 
                  phylvl(i,j,k,ACC_UFLXOLD(iogrp))=
     .              phylvl(i,j,k,ACC_UFLXLVL(iogrp))
                enddo
              enddo
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l)) 
                  phylvl(i,j,k,ACC_VFLXOLD(iogrp))=
     .              phylvl(i,j,k,ACC_VFLXLVL(iogrp))
                enddo
              enddo
            enddo  
          enddo 
c$OMP END PARALLEL DO
        endif
      enddo
c
      if (sum(ACC_UVELLVL(1:nphy)+ACC_UFLXLVL(1:nphy)+
     .  ACC_UTFLXLVL(1:nphy)+ACC_USFLXLVL(1:nphy)+
     .  ACC_UMFLTDLVL(1:nphy)+ACC_UTFLTDLVL(1:nphy)+
     .  ACC_UTFLLDLVL(1:nphy)+ACC_USFLTDLVL(1:nphy)+
     .  ACC_USFLLDLVL(1:nphy)).ne.0) then
        do k=1,kk
          call diazlv('u',k,mm,nn,ind1,ind2,wghts,wghtsflx)
c
c --- --- weighted u-component of total velocity [g/s^3]
          call acclvl(ACC_UVELLVL,uvel,'u',k,ind1,ind2,wghts)
c
c --- --- u-component of mass flux [g*cm/s^2]
          call acclvl(ACC_UFLXLVL,uflx(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of heat flux [K*g*cm/s^2]
          call acclvl(ACC_UTFLXLVL,utflx(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of salt flux [g^2*cm/kg/s^2]
          call acclvl(ACC_USFLXLVL,usflx(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of mass flux due to thickness diffusion [g*cm/s^2]
          call acclvl(ACC_UMFLTDLVL,umfltd(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of heat flux due to thickness diffusion [K*g*cm/s^2]
          call acclvl(ACC_UTFLTDLVL,utfltd(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of heat flux due to lateral diffusion [K*g*cm/s^2]
          call acclvl(ACC_UTFLLDLVL,utflld(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of salt flux due to thickness diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_USFLTDLVL,usfltd(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
c
c --- --- u-component of salt flux due to lateral diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_USFLLDLVL,usflld(1-nbdy,1-nbdy,k1n),
     .                'u',k,ind1,ind2,wghtsflx)
        enddo
      endif
c
      if (sum(ACC_VVELLVL(1:nphy)+ACC_VFLXLVL(1:nphy)+
     .  ACC_VTFLXLVL(1:nphy)+ACC_VSFLXLVL(1:nphy)+
     .  ACC_VMFLTDLVL(1:nphy)+ACC_VTFLTDLVL(1:nphy)+
     .  ACC_VTFLLDLVL(1:nphy)+ACC_VSFLTDLVL(1:nphy)+
     .  ACC_VSFLLDLVL(1:nphy)).ne.0) then
        do k=1,kk
          call diazlv('v',k,mm,nn,ind1,ind2,wghts,wghtsflx)
c
c --- --- weighted v-component of total velocity [g/s^3]
          call acclvl(ACC_VVELLVL,vvel,'v',k,ind1,ind2,wghts)
c
c --- --- v-component of mass flux [g*cm/s^2]
          call acclvl(ACC_VFLXLVL,vflx(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of heat flux [K*g*cm/s^2]
          call acclvl(ACC_VTFLXLVL,vtflx(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of salt flux [g^2*cm/kg/s^2]
          call acclvl(ACC_VSFLXLVL,vsflx(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of mass flux due to thickness diffusion [g*cm/s^2]
          call acclvl(ACC_VMFLTDLVL,vmfltd(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of heat flux due to thickness diffusion [K*g*cm/s^2]
          call acclvl(ACC_VTFLTDLVL,vtfltd(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of heat flux due to lateral diffusion [K*g*cm/s^2]
          call acclvl(ACC_VTFLLDLVL,vtflld(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of salt flux due to thickness diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_VSFLTDLVL,vsfltd(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
c
c --- --- v-component of salt flux due to lateral diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_VSFLLDLVL,vsflld(1-nbdy,1-nbdy,k1n),
     .                'v',k,ind1,ind2,wghtsflx)
        enddo
      endif
c
      if (sum(ACC_SALNLVL(1:nphy)+ACC_TEMPLVL(1:nphy)+
     .  ACC_BFSQLVL(1:nphy)+ACC_DIFDIALVL(1:nphy)+
     .  ACC_DIFINTLVL(1:nphy)+ACC_DIFISOLVL(1:nphy)+
     .  ACC_TKELVL(1:nphy)+ACC_GLS_PSILVL(1:nphy)+
     .  ACC_PVLVL(1:nphy)+ACC_DZLVL(1:nphy)).ne.0) then
        do k=1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
c
c --- --- salinity [g^2/kg/cm/s^2]
          call acclvl(ACC_SALNLVL,saln(1-nbdy,1-nbdy,k1m),
     .                'p',k,ind1,ind2,wghts)
c
c --- --- temperature [degC*g/cm/s^2]
          call acclvl(ACC_TEMPLVL,temp(1-nbdy,1-nbdy,k1m),
     .                'p',k,ind1,ind2,wghts)
c
c --- --- buoyancy frequency squared [1/s]
          call acclvl(ACC_BFSQLVL,bfsql,'p',k,ind1,ind2,wghts)
c
c --- --- layer interface diffusivity [cm^2/s]
          call acclvl(ACC_DIFINTLVL,difint,'p',k,ind1,ind2,wghts)
c
c --- --- isopycnal diffusivity [cm^2/s]
          call acclvl(ACC_DIFISOLVL,difiso,'p',k,ind1,ind2,wghts)
c
c --- --- diapycnal diffusivity [cm^2/s]
          call acclvl(ACC_DIFDIALVL,difdia,'p',k,ind1,ind2,wghts)
c
c --- --- potential vorticity [s m-2]
          call acclvl(ACC_PVLVL,pv_p,'p',k,ind1,ind2,wghts)
c
#if defined(TRC) && defined(TKE)
c --- --- tke [cm2/s2*g/cm/s^2]
          call acclvl(ACC_TKELVL,trc(1-nbdy,1-nbdy,k1m,itrtke),
     .                'p',k,ind1,ind2,wghts)
c
c --- --- gls_psi [cm2/s3*g/cm/s^2]
          call acclvl(ACC_GLS_PSILVL,trc(1-nbdy,1-nbdy,k1m,itrgls),
     .                'p',k,ind1,ind2,wghts)
c
#endif
c --- --- layer geopotential thickness [cm^2/s^2]
          call acclvl(ACC_DZLVL,dz,'p',k,ind1,ind2,wghts)
c
        enddo
      endif
c
c --- Accumulate vertical velocity
      do iogrp=1,nphy
        if (ACC_WFLX(iogrp)+ACC_WFLX2(iogrp)+ACC_WFLXLVL(iogrp)+
     .    ACC_WFLX2LVL(iogrp).ne.0)
     .    call diavfl(iogrp,m,n,mm,nn,k1m,k1n)
      enddo
c
      end subroutine diaacc



      subroutine diaout(iogrp,m,n,mm,nn,k1m,k1n)
c
c --- -------------------------------------
c --- Write high frequency diagnostic fields
c --- -------------------------------------
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_seaice.h"
c
      integer iogrp
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv
      logical iniflg,append2file(nphymax)
      data iniflg /.true./ , append2file /nphymax*.false./
      save iniflg,iuu,ivv,append2file
c
      integer i,j,k,l,cmpflg,irec(nphymax)
      character startdate*20,timeunits*30,fname(nphymax)*256
      real datenum,rnacc
      save fname,irec
c
      real, dimension(itdm,jtdm) :: bflxg,strg
c
#if defined(TRC)
#  include "param_trc.h"
#  include "common_trc.h"
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts,
     .   wghtsflx
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: tmp3d
      integer nt,nat,km
      character trcnm*80,trcnml*80
      real treps
      parameter (treps=1.e-14)
#endif
c
c --- prepare output fields
      if (mnproc.eq.1)
     .  write (lp,'(a,f6.2,a)') 'diaout: fields averaged over ',
     .  real(nacc_phy(iogrp))/nstep_in_day,' days'
c
      rnacc=1./real(nacc_phy(iogrp))
      cmpflg=GLB_COMPFLAG(iogrp)
c 
c --- compute meridional transports and transports through sections
      if (ACC_MMFLXL(iogrp)+ACC_MMFLXD(iogrp)+ACC_MMFTDL(iogrp)
     .   +ACC_MMFTDD(iogrp)+ACC_MHFLX(iogrp)+ACC_MHFTD(iogrp)
     .   +ACC_MHFLD(iogrp)+ACC_MSFLX(iogrp)+ACC_MSFTD(iogrp)
     .   +ACC_MSFLD(iogrp).ne.0)
     .  call diamer(iogrp)
      if (ACC_VOLTR(iogrp).ne.0) call diasec(iogrp)
c
c --- compute barotropic mass streamfunction
      if (H2D_BTMSTR(iogrp).ne.0) then
        if     (nreg.le.2) then
          call xcaget(bflxg,phyh2d(1-nbdy,1-nbdy,ACC_UBFLXS(iogrp)),1)
          if (mnproc.eq.1) then
            do i=1,itdm
              strg(i,1)=0.
            enddo
            do j=1,jtdm-1
              do i=1,itdm
                strg(i,j+1)=strg(i,j)-bflxg(i,j)
              enddo
            enddo
          endif
        elseif (nreg.eq.4) then
          call xcaget(bflxg,phyh2d(1-nbdy,1-nbdy,ACC_VBFLXS(iogrp)),1)
          if (mnproc.eq.1) then
            do j=1,jtdm
              strg(itdm,j)=0.
            enddo
            do j=1,jtdm
              do i=itdm-1,1,-1
                strg(i,j)=strg(i+1,j)-bflxg(i,j)
              enddo
            enddo
          endif
        else
          if (mnproc.eq.1) then
            write (lp,'(a,i2,a)')
     .        'mod_dia: cannot compute streamfunction for nreg =',nreg,
     .        '!'
          endif
          call xcstop('(mod_dia)')
                 stop '(mod_dia)'
        endif
        call xcaput(strg,util1,1)
        call xctilr(util1,1,1, 1,1, halo_qs)
        call inih2d(ACC_UTILH2D(1),'p',0.)
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              phyh2d(i,j,ACC_UTILH2D(1))=
     .          .25*(util1(i  ,j  )+util1(i+1,j  )
     .              +util1(i  ,j+1)+util1(i+1,j+1))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- compute global sums and averages
      if (MSC_MASSGS(iogrp)+MSC_SALNGA(iogrp)
     .   +MSC_TEMPGA(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=phylyr(i,j,1,ACC_DP(iogrp))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=util1(i,j)+phylyr(i,j,k,ACC_DP(iogrp))
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=util1(i,j)*scp2(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(massgs(1),util1,ips)
      endif
      if (MSC_VOLGS(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=phylyr(i,j,1,ACC_DZ(iogrp))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=util1(i,j)+phylyr(i,j,k,ACC_DZ(iogrp))
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=util1(i,j)*scp2(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(volgs(1),util1,ips)
        volgs(1)=rnacc*1e-6*volgs(1)/g
      endif
      if (MSC_SALNGA(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=phylyr(i,j,1,ACC_SALN(iogrp))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=util1(i,j)+phylyr(i,j,k,ACC_SALN(iogrp))
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=util1(i,j)*scp2(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(salnga(1),util1,ips)
        salnga(1)=salnga(1)/massgs(1)
      endif
      if (MSC_TEMPGA(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=phylyr(i,j,1,ACC_TEMP(iogrp))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=util1(i,j)+phylyr(i,j,k,ACC_TEMP(iogrp))
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=util1(i,j)*scp2(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(tempga(1),util1,ips)
        tempga(1)=tempga(1)/massgs(1)
      endif
      if (MSC_MASSGS(iogrp).ne.0) then
        massgs(1)=rnacc*1e-3*massgs(1)/g
      endif
      if (MSC_SSSGA(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=phyh2d(i,j,ACC_SSS(iogrp))*scp2(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(sssga(1),util1,ips)
        sssga(1)=rnacc*sssga(1)/area
      endif
      if (MSC_SSTGA(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=phyh2d(i,j,ACC_SST(iogrp))*scp2(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcsum(sstga(1),util1,ips)
        sstga(1)=rnacc*sstga(1)/area
      endif
c
c --- finalize accumulation of 2d fields
      call finh2d(ACC_HICE(iogrp),ACC_FICE(iogrp),'p')
      call finh2d(ACC_HSNW(iogrp),ACC_FICE(iogrp),'p')
      call finh2d(ACC_MXLU(iogrp),ACC_MLDU(iogrp),'u')
      call finh2d(ACC_MXLV(iogrp),ACC_MLDV(iogrp),'v')
      call finh2d(ACC_UICE(iogrp),ACC_IVOLU(iogrp),'u')
      call finh2d(ACC_VICE(iogrp),ACC_IVOLV(iogrp),'v')
c
c --- finalize accumulation of layer fields
      call finlyr(ACC_UVEL(iogrp),ACC_DPU(iogrp),'u')
      call finlyr(ACC_VVEL(iogrp),ACC_DPV(iogrp),'v')
      call finlyr(ACC_SALN(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_TEMP(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_BFSQ(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_DIFDIA(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_DIFINT(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_DIFISO(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_AVDSG(iogrp),ACC_DPVOR(iogrp),'p')
#if defined(TRC) && defined(TKE)
      call finlyr(ACC_TKE(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_GLS_PSI(iogrp),ACC_DP(iogrp),'p')
#endif
c 
c --- compute log10 of diffusivities
      if (LYR_DIFDIA(iogrp).eq.2)
     .  call loglyr(ACC_DIFDIA(iogrp),'p',1e-4,0.)
      if (LYR_DIFINT(iogrp).eq.2)
     .  call loglyr(ACC_DIFINT(iogrp),'p',1e-4,0.)
      if (LYR_DIFISO(iogrp).eq.2)
     .  call loglyr(ACC_DIFISO(iogrp),'p',1e-4,0.)
c
      if (LVL_DIFDIA(iogrp).eq.2)
     .  call loglvl(ACC_DIFDIALVL(iogrp),'p',1e-4*rnacc,0.)
      if (LVL_DIFINT(iogrp).eq.2)
     .  call loglvl(ACC_DIFINTLVL(iogrp),'p',1e-4*rnacc,0.)
      if (LVL_DIFISO(iogrp).eq.2)
     .  call loglvl(ACC_DIFISOLVL(iogrp),'p',1e-4*rnacc,0.)
c
c --- mask sea floor of level fields
      call msklvl(ACC_BFSQLVL(iogrp),'p')
      call msklvl(ACC_DIFDIALVL(iogrp),'p')
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
      call msklvl(ACC_UTFLTDLVL(iogrp),'u')
      call msklvl(ACC_VTFLTDLVL(iogrp),'v')
      call msklvl(ACC_UTFLLDLVL(iogrp),'u')
      call msklvl(ACC_VTFLLDLVL(iogrp),'v')
      call msklvl(ACC_USFLTDLVL(iogrp),'u')
      call msklvl(ACC_VSFLTDLVL(iogrp),'v')
      call msklvl(ACC_USFLLDLVL(iogrp),'u')
      call msklvl(ACC_VSFLLDLVL(iogrp),'v')
      call msklvl(ACC_SALNLVL(iogrp),'p')
      call msklvl(ACC_TEMPLVL(iogrp),'p')
      call msklvl(ACC_WFLXLVL(iogrp),'p')
      call msklvl(ACC_WFLX2LVL(iogrp),'p')
      call msklvl(ACC_PVLVL(iogrp),'p')
#if defined(TRC) && defined(TKE)
      call msklvl(ACC_TKELVL(iogrp),'p')
      call msklvl(ACC_GLS_PSILVL(iogrp),'p')
#endif
c 
c --- get instantaneous values for ice age 
      if (ACC_IAGE(iogrp).ne.0) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              phyh2d(i,j,ACC_IAGE(iogrp))=iagem(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- set time information
      timeunits=' '
      startdate=' '
      write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)')
     . 'days since ',min(1800,nyear0),'-',1,'-',1,' 00:00'
      write(startdate,'(i4.4,a1,i2.2,a1,i2.2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=time-time0-0.5*diagfq_phy(iogrp)/nstep_in_day
c
c --- create file name 
      if (.not.append2file(iogrp)) then
        call diafnm(runid,runid_len,expcnf,trim(GLB_FNAMETAG(IOGRP)),
     .  nstep,filefq_phy(iogrp)/real(nstep_in_day),filemon_phy(iogrp),
     .  fileann_phy(iogrp),fname(iogrp))
        append2file(iogrp)=.true.
        irec(iogrp)=1
      else
        irec(iogrp)=irec(iogrp)+1
      endif
      if (((fileann_phy(iogrp).and.nday_of_year.eq.1.or.
     .  filemon_phy(iogrp).and.nday.eq.1).and.mod(nstep,nstep_in_day)
     .  .eq.0).or..not.(fileann_phy(iogrp).or.filemon_phy(iogrp)).and.
     .  mod(nstep+.5,filefq_phy(iogrp)).lt.1.) then
        append2file(iogrp)=.false.
      endif
c
c --- open output file
      if (mnproc.eq.1)
     .   write (lp,'(2a)') 'Writing physical diagnostics to file: ',
     .                     trim(fname(iogrp))
      if (GLB_NCFORMAT(iogrp).eq.1) then
        call ncfopn(fname(iogrp),'w','6',irec(iogrp),iotype)
      elseif (GLB_NCFORMAT(iogrp).eq.2) then
        call ncfopn(fname(iogrp),'w','h',irec(iogrp),iotype)
      else
        call ncfopn(fname(iogrp),'w','c',irec(iogrp),iotype)
      endif
c
c --- compute extended ocean masks 
      if (iniflg) then
         iniflg=.false.
c$OMP PARALLEL DO PRIVATE(i)
         do j=1,jj
            do i=1,ii
               if((ip(i,j)+ip(i-1,j)).ge.1) then
                  iuu(i,j)=1
               else
                  iuu(i,j)=0
               endif
               if((ip(i,j)+ip(i,j-1)).ge.1) then
                  ivv(i,j)=1
               else
                  ivv(i,j)=0
               endif
            enddo
         enddo
c$OMP END PARALLEL DO
      endif
c
c --- define output dimensions  
      if (cmpflg.ne.0) then
        call ncdimc('pcomp',ip,0)
        call ncdimc('ucomp',iuu,0)
        call ncdimc('vcomp',ivv,0)
      else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
      endif
      call ncdims('sigma',kdm)
      call ncdims('depth',ddm)
      call ncdims('region',mer_nreg)
      call ncdims('bounds',2) 
      call ncdims('slenmax',slenmax)
      call ncdims('time',0)
c
      if (ACC_VOLTR(iogrp).ne.0) then
        if ((sec_num.gt.0.and.sec_num.le.max_sec)) then
          call ncdims('section',sec_num)
        else
          write (lp,*) 'Illegal number of sections: sec_num=',sec_num,
     .                 ' max_sec=',max_sec
          call xchalt('(diaout)')
                 stop '(diaout)'
        endif
      endif
c
      if (ACC_MMFLXL(iogrp)+ACC_MMFLXD(iogrp)+ACC_MMFTDL(iogrp)
     .   +ACC_MMFTDD(iogrp)+ACC_MHFLX(iogrp)+ACC_MHFTD(iogrp)
     .   +ACC_MHFLD(iogrp)+ACC_MSFLX(iogrp)+ACC_MSFTD(iogrp)
     .   +ACC_MSFLD(iogrp).ne.0) then
        if ((lmax.gt.0.and.lmax.le.ldm)) then
          call ncdims('lat',lmax)
        else
          write (lp,*)
     .      'Illegal dimension of meridional diagnostics: lmax=',lmax,
     .      ' ldm=',ldm
          call xchalt('(diaout)')
                 stop '(diaout)'
        endif
      endif
c
       call definevar(irec(iogrp),iogrp,cmpflg,timeunits,calendar)

       call nctime(datenum,calendar,timeunits,startdate)
c
c --- write auxillary dimension information 
      if (irec(iogrp).eq.1) then
c --- sigma levels
        call ncwrt1('sigma','sigma',sigmar1)
        call ncattr('long_name','Potential density')
        call ncattr('standard_name','sea_water_sigma_theta')
        call ncattr('units','kg m-3')
        call ncattr('positive','down')
c --- zlevel
        call ncwrt1('depth','depth',depthslev)
        call ncattr('long_name','z level')
        call ncattr('units','m')
        call ncattr('positive','down')
        call ncattr('bounds','depth_bnds')
        call ncwrt1('depth_bnds','bounds depth',depthslev_bnds)
        if (MSC_MMFLXL(iogrp)+MSC_MMFLXD(iogrp)+MSC_MMFTDL(iogrp)
     .     +MSC_MMFTDD(iogrp)+MSC_MHFLX(iogrp)+MSC_MHFTD(iogrp)
     .     +MSC_MHFLD(iogrp)+MSC_MSFLX(iogrp)+MSC_MSFTD(iogrp)
     .     +MSC_MSFLD(iogrp)+MSC_VOLTR(iogrp).ne.0) then 
          call ncwrt1('lat','lat',mtlat)  
          call ncattr('long_name','Latitude')
          call ncattr('standard_name','latitude')
          call ncattr('units','degree_north')
        endif 
        if (MSC_MMFLXL(iogrp)+MSC_MMFLXD(iogrp)+MSC_MMFTDL(iogrp)
     .     +MSC_MMFTDD(iogrp).ne.0) then
          call ncwrtc('region','slenmax region',mer_regnam)  
          call ncattr('long_name','Region name') 
        endif
        if (MSC_MHFLX(iogrp)+MSC_MHFTD(iogrp)+MSC_MHFLD(iogrp)
     .     +MSC_MSFLX(iogrp)+MSC_MSFTD(iogrp)+MSC_MSFLD(iogrp).ne.0)
     .  then 
          call ncwrtc('section','slenmax section',sec_name)
          call ncattr('long_name','Section name') 
        endif 
      endif
c
c --- write 2d fields
      call wrth2d(ACC_SIGMX(iogrp),H2D_SIGMX(iogrp),rnacc*1e3,
     ,  0.,cmpflg,ip,'p','sigmx','Mixed layer density',' ','kg m-3')
c
      call wrth2d(ACC_UB(iogrp),H2D_UB(iogrp),rnacc*1e-2,
     .  0.,cmpflg,iuu,'u','ubaro','Barotropic velocity x-component',
     .  ' ','m s-1')
c
      call wrth2d(ACC_VB(iogrp),H2D_VB(iogrp),rnacc*1e-2,
     .  0.,cmpflg,ivv,'v','vbaro','Barotropic velocity y-component',
     .  ' ','m s-1')
c
      call wrth2d(ACC_PSRF(iogrp),H2D_PSRF(iogrp),
     .  rnacc*.1,0.,cmpflg,ip,'p','psrf','Surface pressure',' ','Pa')
c
      call wrth2d(ACC_PBOT(iogrp),H2D_PBOT(iogrp),
     .  rnacc*.1,0.,cmpflg,ip,'p','pbot','Bottom pressure',' ','Pa')
c
      call wrth2d(ACC_SEALV(iogrp),H2D_SEALV(iogrp),
     .  -rnacc*1e-2,0.,cmpflg,ip,'p','sealv','Sea level',' ','m')
c
      call wrth2d(ACC_SLVSQ(iogrp),H2D_SLVSQ(iogrp),
     .  rnacc*1e-4,0.,cmpflg,ip,'p','slvsq','Sea level squared',' ',
     .  'm2')
c
      call wrth2d(ACC_UTILH2D(1),H2D_BTMSTR(iogrp),
     .  rnacc*0.5e-3*dlt/(g*baclin),0.,cmpflg,ip,'p','btmstr',
     .  'Barotropic mass streamfunction',' ','kg s-1')
c
      call wrth2d(ACC_HICE(iogrp),H2D_HICE(iogrp),1.,0.,
     .  cmpflg,ip,'p','hice','Ice thickness',' ','m')
c
      call wrth2d(ACC_TICE(iogrp),H2D_TICE(iogrp),rnacc,
     .  -t0deg,cmpflg,ip,'p','tice','Ice temperature',' ','degC')
c
      call wrth2d(ACC_HSNW(iogrp),H2D_HSNW(iogrp),1.,0.,
     .  cmpflg,ip,'p','hsnw','Snow depth',' ','m')
c
      call wrth2d(ACC_FICE(iogrp),H2D_FICE(iogrp),
     .  rnacc*1e2,0.,cmpflg,ip,'p','fice','Ice concentration',' ','%')
c
      call wrth2d(ACC_TSRF(iogrp),H2D_TSRF(iogrp),rnacc,
     .  -t0deg,cmpflg,ip,'p','tsrf','Surface temperature',' ','degC')
c
      call wrth2d(ACC_IAGE(iogrp),H2D_IAGE(iogrp),1.,0.,
     .  cmpflg,ip,'p','iage','Ice age',' ','day')
c
      call wrth2d(ACC_UICE(iogrp),H2D_UICE(iogrp),1e-2,0.,
     .  cmpflg,iuu,'u','uice','Ice velocity x-component',' ','m s-1')
c
      call wrth2d(ACC_VICE(iogrp),H2D_VICE(iogrp),1e-2,0.,
     .  cmpflg,ivv,'v','vice','Ice velocity y-component',' ','m s-1')
c
      call wrth2d(ACC_SWA(iogrp),H2D_SWA(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','swa','Short-wave heat flux',' ','W m-2')
c
      call wrth2d(ACC_NSF(iogrp),H2D_NSF(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','nsf','Non-solar heat flux',' ','W m-2')
c
      call wrth2d(ACC_HMLTFZ(iogrp),H2D_HMLTFZ(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','hmltfz',
     .  'Heat flux due to melting/freezing',' ','W m-2')
c
      call wrth2d(ACC_DFL(iogrp),H2D_DFL(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','dfl','Non-solar heat flux derivative',' ',
     .  'W m-2 K-1')
c
      call wrth2d(ACC_SURFLX(iogrp),H2D_SURFLX(iogrp),
     .  -rnacc*1e4,0.,cmpflg,ip,'p','hflx',
     .  'Heat flux received by ocean',' ','W m-2')
c
      call wrth2d(ACC_SURRLX(iogrp),H2D_SURRLX(iogrp),
     .  -rnacc*1e4,0.,cmpflg,ip,'p','hrflx',
     .  'Restoring heat flux received by ocean',' ','W m-2')
c
      call wrth2d(ACC_LIP(iogrp),H2D_LIP(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','lip','Liquid precipitation',' ','kg m-2 s-1')
c
      call wrth2d(ACC_SOP(iogrp),H2D_SOP(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sop','Solid precipitation',' ','kg m-2 s-1')
c
      call wrth2d(ACC_EVA(iogrp),H2D_EVA(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','eva','Evaporation',' ','kg m-2 s-1')
c
      call wrth2d(ACC_FMLTFZ(iogrp),H2D_FMLTFZ(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','fmltfz',
     .  'Fresh water flux due to melting/freezing',' ','kg m-2 s-1')
c
      call wrth2d(ACC_RNFFLX(iogrp),H2D_RNFFLX(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','rnf','Liquid runoff',' ','kg m-2 s-1')
c
      call wrth2d(ACC_RFIFLX(iogrp),H2D_RFIFLX(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','rfi','Frozen runoff',' ','kg m-2 s-1')
c
      call wrth2d(ACC_SALFLX(iogrp),H2D_SALFLX(iogrp),
     .  -rnacc*1e-2,0.,cmpflg,ip,'p','sflx',
     .  'Salt flux received by ocean',' ','kg m-2 s-1')
c
      call wrth2d(ACC_SALRLX(iogrp),H2D_SALRLX(iogrp),
     .  -rnacc*1e-2,0.,cmpflg,ip,'p','srflx',
     .  'Restoring salt flux received by ocean',' ','kg m-2 s-1')
c
      call wrth2d(ACC_BRNFLX(iogrp),H2D_BRNFLX(iogrp),
     .  rnacc*(-1e-2),0.,cmpflg,ip,'p','bflx','Brine flux',' ',
     .  'kg m-2 s-1')
c
      call wrth2d(ACC_ZTX(iogrp),H2D_ZTX(iogrp),rnacc,0.,
     .  cmpflg,iuu,'u','ztx','Wind stress x-component',' ','N m-2')
c
      call wrth2d(ACC_MTY(iogrp),H2D_MTY(iogrp),rnacc,0.,
     .  cmpflg,ivv,'v','mty','Wind stress y-component',' ','N m-2')
c
      call wrth2d(ACC_TAUX(iogrp),H2D_TAUX(iogrp),rnacc*.1,
     .  0.,cmpflg,iuu,'u','taux',
     .  'Momentum flux received by ocean x-component',' ','N m-2')
c
      call wrth2d(ACC_TAUY(iogrp),H2D_TAUY(iogrp),rnacc*.1,
     .  0.,cmpflg,ivv,'v','tauy', 
     .  'Momentum flux received by ocean y-component',' ','N m-2')
c
      call wrth2d(ACC_IDKEDT(iogrp),H2D_IDKEDT(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','idkedt',
     .  'Mixed layer inertial kinetic energy tendency per unit area',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_USTAR(iogrp),H2D_USTAR(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','ustar','Friction velocity',' ',
     .  'm s-1')
c
      call wrth2d(ACC_USTAR3(iogrp),H2D_USTAR3(iogrp),
     .  rnacc*1.e-6,0.,cmpflg,ip,'p','ustar3',
     .  'Friction velocity cubed',' ','m3 s-3')
c
      call wrth2d(ACC_ABSWND(iogrp),H2D_ABSWND(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','abswnd','Absolute wind speed',' ',
     .  'm s-1')
c
      call wrth2d(ACC_MTKEUS(iogrp),H2D_MTKEUS(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','mtkeus',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to friction velocity',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_MTKENI(iogrp),H2D_MTKENI(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','mtkeni',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to near inertial motions',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_MTKEBF(iogrp),H2D_MTKEBF(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','mtkebf',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to buoyancy forcing',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_MTKERS(iogrp),H2D_MTKERS(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','mtkers',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to eddy restratification',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_MTKEPE(iogrp),H2D_MTKEPE(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','mtkepe',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to potential energy change',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_MTKEKE(iogrp),H2D_MTKEKE(iogrp),
     .  rnacc*1.e-3/alpha0,0.,cmpflg,ip,'p','mtkeke',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to kinetic energy change',
     .  ' ','kg s-3')
c
      call wrth2d(ACC_SFL(iogrp),H2D_SFL(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sfl','Salt flux',' ','kg m-2 s-1')
c
      call wrth2d(ACC_ALB(iogrp),H2D_ALB(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','alb','Surface albedo',' ','1')
c
      call wrth2d(ACC_MLD(iogrp),H2D_MLD(iogrp),rnacc/onem,
     .  0.,cmpflg,ip,'p','mld','Mixed layer depth',' ','m')
c
      call wrth2d(ACC_MLDU(iogrp),H2D_MLDU(iogrp),
     .  rnacc/onem,0.,cmpflg,iuu,'u','mldu',
     .  'Mixed layer depth at u-point',' ','m')
c
      call wrth2d(ACC_MLDV(iogrp),H2D_MLDV(iogrp),
     .  rnacc/onem,0.,cmpflg,ivv,'v','mldv',
     .  'Mixed layer depth at v-point',' ','m')
c
      call wrth2d(ACC_MAXMLD(iogrp),H2D_MAXMLD(iogrp),
     .  1./onem,0.,cmpflg,ip,'p','maxmld','Maximum mixed layer depth',
     .  ' ','m')
c
      call wrth2d(ACC_MLTS(iogrp),H2D_MLTS(iogrp),rnacc*1e-2,
     .  0.,cmpflg,ip,'p','mlts',
     . 'Mixed layer thickness defined by sigma t',' ','m')
c
      call wrth2d(ACC_MLTSMN(iogrp),H2D_MLTSMN(iogrp),1e-2,
     .  0.,cmpflg,ip,'p','mltsmn',
     . 'Minimum mixed layer thickness defined by sigma t',' ','m')
c
      call wrth2d(ACC_MLTSMX(iogrp),H2D_MLTSMX(iogrp),1e-2,
     .  0.,cmpflg,ip,'p','mltsmx',
     . 'Maximum mixed layer thickness defined by sigma t',' ','m')
c
      call wrth2d(ACC_MLTSSQ(iogrp),H2D_MLTSSQ(iogrp),rnacc*1e-4,
     .  0.,cmpflg,ip,'p','mltssq',
     . 'Mixed layer thickness squared defined by sigma t',' ','m')
c
      call wrth2d(ACC_T20D(iogrp),H2D_T20D(iogrp),rnacc*1e-2,
     .  0.,cmpflg,ip,'p','t20d','20C isoterm depth',' ','m')
c
      call wrth2d(ACC_BRNPD(iogrp),H2D_BRNPD(iogrp),rnacc/onem,
     .  0.,cmpflg,ip,'p','brnpd','Brine plume depth',' ','m')
c
      call wrth2d(ACC_SSS(iogrp),H2D_SSS(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sss','Ocean surface salinity',' ','g kg-1')
c
      call wrth2d(ACC_SSSSQ(iogrp),H2D_SSSSQ(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','ssssq','Ocean surface salinity squared',' ',
     .  'g2 kg-2')
c
      call wrth2d(ACC_SBOT(iogrp),H2D_SBOT(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sbot','Bottom salinity',' ','g kg-1')
c
      call wrth2d(ACC_SST(iogrp),H2D_SST(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sst','Ocean surface temperature',' ','degC')
c
      call wrth2d(ACC_SSTSQ(iogrp),H2D_SSTSQ(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sstsq','Ocean surface temperature squared',' ',
     . 'degC2')
c
      call wrth2d(ACC_TBOT(iogrp),H2D_TBOT(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','tbot','Bottom temperature',' ','degC')
c
      call wrth2d(ACC_MXLU(iogrp),H2D_MXLU(iogrp),1e-2,0.,
     .  cmpflg,iuu,'u','mxlu','Mixed layer velocity x-component',' ',
     .  'm s-1')
c
      call wrth2d(ACC_MXLV(iogrp),H2D_MXLV(iogrp),1e-2,0.,
     .  cmpflg,ivv,'v','mxlv','Mixed layer velocity y-component',' ',
     .  'm s-1')
c
c --- write 3d layer fields
      call wrtlyr(ACC_DP(iogrp),LYR_DP(iogrp),rnacc*.1,0.
     .  ,cmpflg,ip,'p','dp','Layer pressure thickness',' ','Pa')
c
      call wrtlyr(ACC_DZ(iogrp),LYR_DZ(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','dz','Layer thickness',' ','m')
c
      call wrtlyr(ACC_TEMP(iogrp),LYR_TEMP(iogrp),1.,0.,
     .  cmpflg,ip,'p','temp','Temperature','Ocean temperature',
     .  'degC')
c
      call wrtlyr(ACC_SALN(iogrp),LYR_SALN(iogrp),1.,0.,
     .  cmpflg,ip,'p','saln','Salinity','Ocean salinity','g kg-1')
c
      call wrtlyr(ACC_UVEL(iogrp),LYR_UVEL(iogrp),1e-2,
     .  0.,cmpflg,iuu,'u','uvel','Velocity x-component',' ','m s-1')
c
      call wrtlyr(ACC_VVEL(iogrp),LYR_VVEL(iogrp),1e-2,
     .  0.,cmpflg,ivv,'v','vvel','Velocity y-component',' ','m s-1')
c
      call wrtlyr(ACC_UFLX(iogrp),LYR_UFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,iuu,'u','uflx',
     .  'Mass flux in x-direction',' ','kg s-1')
c
      call wrtlyr(ACC_VFLX(iogrp),LYR_VFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ivv,'v','vflx',
     .  'Mass flux in y-direction',' ','kg s-1')
c
      call wrtlyr(ACC_UTFLX(iogrp),LYR_UTFLX(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,iuu,'u','uhflx',
     .  'Heat flux in x-direction',' ','W')
c
      call wrtlyr(ACC_VTFLX(iogrp),LYR_VTFLX(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,ivv,'v','vhflx',
     .  'Heat flux in y-direction',' ','W')
c
      call wrtlyr(ACC_USFLX(iogrp),LYR_USFLX(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usflx',
     .  'Salt flux in x-direction',' ','kg s-1')
c
      call wrtlyr(ACC_VSFLX(iogrp),LYR_VSFLX(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsflx',
     .  'Salt flux in y-direction',' ','kg s-1')
c
      call wrtlyr(ACC_UMFLTD(iogrp),LYR_UMFLTD(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,iuu,'u','umfltd',
     .  'Mass flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1')
c
      call wrtlyr(ACC_VMFLTD(iogrp),LYR_VMFLTD(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ivv,'v','vmfltd',
     .  'Mass flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1')
c
      call wrtlyr(ACC_UTFLTD(iogrp),LYR_UTFLTD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,iuu,'u','uhfltd',
     .  'Heat flux due to thickness diffusion in x-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_VTFLTD(iogrp),LYR_VTFLTD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,ivv,'v','vhfltd',
     .  'Heat flux due to thickness diffusion in y-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_UTFLLD(iogrp),LYR_UTFLLD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,iuu,'u','uhflld',
     .  'Heat flux due to lateral diffusion in x-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_VTFLLD(iogrp),LYR_VTFLLD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,ivv,'v','vhflld',
     .  'Heat flux due to lateral diffusion in y-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_USFLTD(iogrp),LYR_USFLTD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usfltd',
     .  'Salt flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1')
c
      call wrtlyr(ACC_VSFLTD(iogrp),LYR_VSFLTD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsfltd',
     .  'Salt flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1')
c
      call wrtlyr(ACC_USFLLD(iogrp),LYR_USFLLD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usflld',
     .  'Salt flux due to lateral diffusion in x-direction',' ',
     .  'kg s-1')
c
      call wrtlyr(ACC_VSFLLD(iogrp),LYR_VSFLLD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsflld',
     .  'Salt flux due to lateral diffusion in y-direction',' ',
     .  'kg s-1')
c
      call wrtlyr(ACC_WFLX(iogrp),LYR_WFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ip,'p','wflx',
     .  'Vertical mass flux',' ','kg s-1')
c
      call wrtlyr(ACC_WFLX2(iogrp),LYR_WFLX2(iogrp),
     .  rnacc*(0.5e-3/(g*baclin))**2,0.,cmpflg,ip,'p','wflx2',
     .  'Vertical mass flux squared',' ','kg2 s-2')
c
      call wrtlyr(ACC_BFSQ(iogrp),LYR_BFSQ(iogrp),1.,0.,
     .  cmpflg,ip,'p','bfsq','Squared buoyancy frequency',' ',
     .  's-1')
c
      call wrtlyr(ACC_AVDSG(iogrp),LYR_PV(iogrp),
     .  1.e2*g,0.,cmpflg,ip,'p','pv','Potential vorticity',' ',
     .  'm-1 s-1')
c
      if (LYR_DIFINT(iogrp).eq.2) then
        call wrtlyr(ACC_DIFINT(iogrp),LYR_DIFINT(iogrp),1.,
     .    0.,cmpflg,ip,'p','difint','Layer interface diffusivity',
     .    ' ','log10(m2 s-1)')
      else
        call wrtlyr(ACC_DIFINT(iogrp),LYR_DIFINT(iogrp),1e-4,
     .    0.,cmpflg,ip,'p','difint','Layer interface diffusivity',
     .    ' ','m2 s-1')
      endif
c
      if (LYR_DIFISO(iogrp).eq.2) then
        call wrtlyr(ACC_DIFISO(iogrp),LYR_DIFISO(iogrp),1.,
     .    0.,cmpflg,ip,'p','difiso','Isopycnal diffusivity',' ',
     .    'log10(m2 s-1)')
      else
        call wrtlyr(ACC_DIFISO(iogrp),LYR_DIFISO(iogrp),1e-4,
     .    0.,cmpflg,ip,'p','difiso','Isopycnal diffusivity',' ',
     .    'm2 s-1')
      endif
c
      if (LYR_DIFDIA(iogrp).eq.2) then
        call wrtlyr(ACC_DIFDIA(iogrp),LYR_DIFDIA(iogrp),1.,
     .    0.,cmpflg,ip,'p','difdia','Diapycnal diffusivity',' ',
     .    'log10(m2 s-1)')
      else
        call wrtlyr(ACC_DIFDIA(iogrp),LYR_DIFDIA(iogrp),1e-4,
     .    0.,cmpflg,ip,'p','difdia','Diapycnal diffusivity',' ',
     .    'm2 s-1')
      endif
c
#if defined(TRC) && defined(TKE)
      call wrtlyr(ACC_TKE(iogrp),LYR_TKE(iogrp),1e-4,0.,
     .  cmpflg,ip,'p','tke','TKE','Turbulent kinetic energy',
     .  'm2 s-2')
c
      call wrtlyr(ACC_GLS_PSI(iogrp),LYR_GLS_PSI(iogrp),1.e-4,0.,
     .  cmpflg,ip,'p','gls_psi','GLS_PSI','Generic length scale',
     .  'm2 s-3')
c
#endif
c --- Write 3d depth fields
      call wrtlvl(ACC_DZLVL(iogrp),LVL_DZ(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','dzlvl','Layer thickness',' ','m')
c
      call wrtlvl(ACC_TEMPLVL(iogrp),LVL_TEMP(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','templvl','Temperature',
     .  'Ocean temperature','degC')
c
      call wrtlvl(ACC_SALNLVL(iogrp),LVL_SALN(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','salnlvl','Salinity',
     .  'Ocean salinity','g kg-1')
c
      call wrtlvl(ACC_UVELLVL(iogrp),LVL_UVEL(iogrp),
     .  rnacc*1e-2,0.,cmpflg,iuu,'u','uvellvl',
     .  'Velocity x-component',' ','m s-1')
c
      call wrtlvl(ACC_VVELLVL(iogrp),LVL_VVEL(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ivv,'v','vvellvl',
     .  'Velocity y-component',' ','m s-1')
c 
      call wrtlvl(ACC_UFLXLVL(iogrp),LVL_UFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,iuu,'u','uflxlvl',
     .  'Mass flux in x-direction',' ','kg s-1')
c
      call wrtlvl(ACC_VFLXLVL(iogrp),LVL_VFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ivv,'v','vflxlvl',
     .  'Mass flux in y-direction',' ','kg s-1')
c
      call wrtlvl(ACC_UTFLXLVL(iogrp),LVL_UTFLX(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,iuu,'u','uhflxlvl',
     .  'Heat flux in x-direction',' ','W')
c
      call wrtlvl(ACC_VTFLXLVL(iogrp),LVL_VTFLX(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,ivv,'v','vhflxlvl',
     .  'Heat flux in y-direction',' ','W')
c
      call wrtlvl(ACC_USFLXLVL(iogrp),LVL_USFLX(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usflxlvl',
     .  'Salt flux in x-direction',' ','kg s-1')
c
      call wrtlvl(ACC_VSFLXLVL(iogrp),LVL_VSFLX(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsflxlvl',
     .  'Salt flux in y-direction',' ','kg s-1')
c
      call wrtlvl(ACC_UMFLTDLVL(iogrp),LVL_UMFLTD(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,iuu,'u','umfltdlvl',
     .  'Mass flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1')
c
      call wrtlvl(ACC_VMFLTDLVL(iogrp),LVL_VMFLTD(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ivv,'v','vmfltdlvl',
     .  'Mass flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1')
c
      call wrtlvl(ACC_UTFLTDLVL(iogrp),LVL_UTFLTD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,iuu,'u','uhfltdlvl',
     .  'Heat flux due to thickness diffusion in x-direction',' ',
     .  'W')
c
      call wrtlvl(ACC_VTFLTDLVL(iogrp),LVL_VTFLTD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,ivv,'v','vhfltdlvl',
     .  'Heat flux due to thickness diffusion in y-direction',' ',
     .  'W')
c
      call wrtlvl(ACC_UTFLLDLVL(iogrp),LVL_UTFLLD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,iuu,'u','uhflldlvl',
     .  'Heat flux due to lateral diffusion in x-direction',' ',
     .  'W')
c
      call wrtlvl(ACC_VTFLLDLVL(iogrp),LVL_VTFLLD(iogrp),
     .  rnacc*spcifh*0.5/(g*baclin),0.,cmpflg,ivv,'v','vhflldlvl',
     .  'Heat flux due to lateral diffusion in y-direction',' ',
     .  'W')
c
      call wrtlvl(ACC_USFLTDLVL(iogrp),LVL_USFLTD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usfltdlvl',
     .  'Salt flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1')
c
      call wrtlvl(ACC_VSFLTDLVL(iogrp),LVL_VSFLTD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsfltdlvl',
     .  'Salt flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1')
c
      call wrtlvl(ACC_USFLLDLVL(iogrp),LVL_USFLLD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usflldlvl',
     .  'Salt flux due to lateral diffusion in x-direction',' ',
     .  'kg s-1')
c
      call wrtlvl(ACC_VSFLLDLVL(iogrp),LVL_VSFLLD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsflldlvl',
     .  'Salt flux due to lateral diffusion in y-direction',' ',
     .  'kg s-1')
c
      call wrtlvl(ACC_WFLXLVL(iogrp),LVL_WFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ip,'p','wflxlvl',
     .  'Vertical mass flux',' ','kg s-1')
c
      call wrtlvl(ACC_WFLX2LVL(iogrp),LVL_WFLX2(iogrp),
     .  rnacc*(0.5e-3/(g*baclin))**2,0.,cmpflg,ip,'p','wflx2lvl',
     .  'Vertical mass flux squared',' ','kg2 s-2')
c
      call wrtlvl(ACC_BFSQLVL(iogrp),LVL_BFSQ(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','bfsqlvl','Squared buoyancy frequency',
     .  ' ','s-1')
c
      call wrtlvl(ACC_PVLVL(iogrp),LVL_PV(iogrp),
     .  rnacc*1.e2*g,0.,cmpflg,ip,'p','pvlvl','Potential vorticity',
     .  ' ','m-1 s-1')
c
      if (LVL_DIFINT(iogrp).eq.2) then
        call wrtlvl(ACC_DIFINTLVL(iogrp),LVL_DIFINT(iogrp),1.,
     .    0.,cmpflg,ip,'p','difintlvl','Layer interface diffusivity',
     .    ' ','log10(m2 s-1)')
      else
        call wrtlvl(ACC_DIFINTLVL(iogrp),LVL_DIFINT(iogrp),1e-4*rnacc,
     .    0.,cmpflg,ip,'p','difintlvl','Layer interface diffusivity',
     .    ' ','m2 s-1')
      endif
c
      if (LVL_DIFISO(iogrp).eq.2) then
        call wrtlvl(ACC_DIFISOLVL(iogrp),LVL_DIFISO(iogrp),1.,
     .    0.,cmpflg,ip,'p','difisolvl','Isopycnal diffusivity',' ',
     .    'log10(m2 s-1)')
      else
        call wrtlvl(ACC_DIFISOLVL(iogrp),LVL_DIFISO(iogrp),1e-4*rnacc,
     .    0.,cmpflg,ip,'p','difisolvl','Isopycnal diffusivity',' ',
     .    'm2 s-1')
      endif
c
      if (LVL_DIFDIA(iogrp).eq.2) then
        call wrtlvl(ACC_DIFDIALVL(iogrp),LVL_DIFDIA(iogrp),1.,
     .    0.,cmpflg,ip,'p','difdialvl','Diapycnal diffusivity',' ',
     .    'log10(m2 s-1)')
      else
        call wrtlvl(ACC_DIFDIALVL(iogrp),LVL_DIFDIA(iogrp),1e-4*rnacc,
     .    0.,cmpflg,ip,'p','difdialvl','Diapycnal diffusivity',' ',
     .    'm2 s-1')
      endif
c
#if defined(TRC) && defined(TKE)
      call wrtlvl(ACC_TKELVL(iogrp),LVL_TKE(iogrp),rnacc*1.e-4,
     .  0.,cmpflg,ip,'p','tkelvl','Turbulent kinetic energy',' ',
     .  'm2 s-2')
c
      call wrtlvl(ACC_GLS_PSILVL(iogrp),LVL_GLS_PSI(iogrp),rnacc*1.e-4,
     .  0.,cmpflg,ip,'p','gls_psilvl','Generic length scale',' ',
     .  'm2 s-3')
c
#endif
c
c --- store meridional transports 
      if (MSC_MMFLXL(iogrp).ne.0) then 
        call ncwrt1('mmflxl','lat sigma region time',mmflxl)
        call ncattr('long_name',
     .    'Overturning stream-function on isopycnic layers') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MMFLXD(iogrp).ne.0) then 
        call ncwrt1('mmflxd','lat depth region time',mmflxd)
        call ncattr('long_name',
     .    'Overturning stream-function on z-levels') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MMFTDL(iogrp).ne.0) then 
        call ncwrt1('mmftdl','lat sigma region time',mmftdl)
        call ncattr('long_name',
     .    'Overturning stream-function due to thickness diffusion '//
     .    'on isopycnic layers') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MMFTDD(iogrp).ne.0) then 
        call ncwrt1('mmftdd','lat depth region time',mmftdd)
        call ncattr('long_name',
     .    'Overturning stream-function due to thickness diffusion '//
     .    'on z-levels') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MHFLX(iogrp).ne.0) then 
        call ncwrt1('mhflx','lat region time',mhflx)
        call ncattr('long_name','Meridional heat flux') 
        call ncattr('units','W') 
      endif
      if (MSC_MHFTD(iogrp).ne.0) then 
        call ncwrt1('mhftd','lat region time',mhftd)
        call ncattr('long_name',
     .    'Meridional heat flux due to thickness diffusion') 
        call ncattr('units','W') 
      endif
      if (MSC_MHFLD(iogrp).ne.0) then 
        call ncwrt1('mhfld','lat region time',mhfld)
        call ncattr('long_name',
     .    'Meridional heat flux due to lateral diffusion') 
        call ncattr('units','W') 
      endif
      if (MSC_MSFLX(iogrp).ne.0) then 
        call ncwrt1('msflx','lat region time',msflx)
        call ncattr('long_name','Meridional salt flux') 
        call ncattr('units','kg s-1') 
      endif
      if (MSC_MSFTD(iogrp).ne.0) then 
        call ncwrt1('msftd','lat region time',msftd)
        call ncattr('long_name',
     .    'Meridional salt flux due to thickness diffusion') 
        call ncattr('units','kg s-1') 
      endif
      if (MSC_MSFLD(iogrp).ne.0) then 
        call ncwrt1('msfld','lat region time',msfld)
        call ncattr('long_name',
     .    'Meridional salt flux due to lateral diffusion') 
        call ncattr('units','kg s-1') 
      endif
c
c --- store section transports 
      if (MSC_VOLTR(iogrp).ne.0) then       
        call ncwrt1('voltr','section time',voltr)
        call ncattr('long_name','Section transports')
        call ncattr('units','kg s-1')
      endif
c
c --- store global sums and averages
      if (MSC_MASSGS(iogrp).ne.0) then
        call ncwrt1('massgs','time',massgs)
        call ncattr('long_name','Sea water mass')
        call ncattr('units','kg')
      endif
      if (MSC_VOLGS(iogrp).ne.0) then
        call ncwrt1('volgs','time',volgs)
        call ncattr('long_name','Sea water volume')
        call ncattr('units','m3')
      endif
      if (MSC_SALNGA(iogrp).ne.0) then
        call ncwrt1('salnga','time',salnga)
        call ncattr('long_name','Global average salinity')
        call ncattr('units','g kg-1')
      endif
      if (MSC_TEMPGA(iogrp).ne.0) then
        call ncwrt1('tempga','time',tempga)
        call ncattr('long_name','Global average temperature')
        call ncattr('units','degC')
      endif
      if (MSC_SSSGA(iogrp).ne.0) then
        call ncwrt1('sssga','time',sssga)
        call ncattr('long_name','Global average sea surface salinity')
        call ncattr('units','g kg-1')
      endif
      if (MSC_SSTGA(iogrp).ne.0) then
        call ncwrt1('sstga','time',sstga)
        call ncattr('long_name',
     .              'Global average sea surface temperature')
        call ncattr('units','degC')
      endif
c
#ifdef TRC
      if (LYR_IDLAGE(iogrp).ne.0.or.LYR_TRC(iogrp).ne.0) then
        call inilyr(ACC_UTILLYR(1),'p',0.)
        call acclyr(ACC_UTILLYR,dp(1-nbdy,1-nbdy,k1m),tmp3d,0,'p')
        call wrtlyr(ACC_UTILLYR(1),
     .    max(LYR_IDLAGE(iogrp),LYR_TRC(iogrp)),1.,0.,cmpflg,ip,'p',
     .    'dp_trc','Layer pressure thickness',' ','Pa')
      endif
#  ifdef IDLAGE
c
c --- ideal age tracer 
      if (LYR_IDLAGE(iogrp).ne.0) then
        call inilyr(ACC_UTILLYR(1),'p',0.)
        call acclyr(ACC_UTILLYR,trc(1-nbdy,1-nbdy,k1m,itriag),tmp3d,0,
     .    'p')
        call wrtlyr(ACC_UTILLYR(1),LYR_IDLAGE(iogrp),1.,0.,
     .    cmpflg,ip,'p','idlage','Ideal age',
     .    'sea_water_age_since_surface_contact','year')
      endif
      if (LVL_IDLAGE(iogrp).ne.0) then
        call inilvl(ACC_UTILLVL(1),'p',0.)
        do k=1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
          call acclvl(ACC_UTILLVL,trc(1-nbdy,1-nbdy,k1m,itriag),'p',
     .      k,ind1,ind2,wghts)
        enddo
        call msklvl(ACC_UTILLVL(1),'p')
        call wrtlvl(ACC_UTILLVL(1),LVL_IDLAGE(iogrp),
     .  1.,0.,cmpflg,ip,'p','idlagelvl','Ideal age',
     .   'sea_water_age_since_surface_contact','year')
      endif
#  endif
c
c --- ocean tracers
      if (LYR_TRC(iogrp).gt.0.and.ntrocn.gt.0) then
#  ifdef ATRC
        do nt=1,ntrocn-natr
          call inilyr(ACC_UTILLYR(1),'p',0.)
          call acclyr(ACC_UTILLYR,trc(1-nbdy,1-nbdy,k1m,nt),tmp3d,0,'p')
          write (trcnm,'(a,i3.3)') 'trc',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call wrtlyr(ACC_UTILLYR(1),LYR_TRC(iogrp),1.,0.,
     .      cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
        enddo
        do nt=1,natr
          nat=ntr-natr+nt
          call inilyr(ACC_UTILLYR(1),'p',0.)
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
          do j=1,jj
            do k=1,kk
              km=k+mm
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                tmp3d(i,j,k)=trc(i,j,km,nat)/max(trc(i,j,km,nt),treps)
              enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
          call acclyr(ACC_UTILLYR,tmp3d,tmp3d,0,'p')
          write (trcnm,'(a,i3.3)') 'atrc',nt
          write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
          call wrtlyr(ACC_UTILLYR(1),LYR_TRC(iogrp),1.,0.,
     .      cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
        enddo
#  else
        do nt=1,ntrocn
          call inilyr(ACC_UTILLYR(1),'p',0.)
          call acclyr(ACC_UTILLYR,trc(1-nbdy,1-nbdy,k1m,nt),tmp3d,0,'p')
          write (trcnm,'(a,i3.3)') 'trc',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call wrtlyr(ACC_UTILLYR(1),LYR_TRC(iogrp),1.,0.,
     .      cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
        enddo
#  endif
      endif
      if (LVL_TRC(iogrp).gt.0.and.ntrocn.gt.0) then
#  ifdef ATRC
        do nt=1,ntrocn-natr
          call inilvl(ACC_UTILLVL(1),'p',0.)
          do k=1,kk
            call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
            call acclvl(ACC_UTILLVL,trc(1-nbdy,1-nbdy,k1m,nt),'p',
     .        k,ind1,ind2,wghts)
          enddo
          call msklvl(ACC_UTILLVL(1),'p')
          write (trcnm,'(a,i3.3)') 'trclvl',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call wrtlvl(ACC_UTILLVL(1),LVL_TRC(iogrp),
     .    1.,0.,cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
        enddo
        do nt=1,natr
          nat=ntr-natr+nt
          call inilvl(ACC_UTILLVL(1),'p',0.)
c$OMP PARALLEL DO PRIVATE(k,km,l,i)
          do j=1,jj
            do k=1,kk
              km=k+mm
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                tmp3d(i,j,k)=trc(i,j,km,nat)/max(trc(i,j,km,nt),treps)
              enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
          do k=1,kk
            call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
            call acclvl(ACC_UTILLVL,tmp3d,'p',k,ind1,ind2,wghts)
          enddo
          call msklvl(ACC_UTILLVL(1),'p')
          write (trcnm,'(a,i3.3)') 'atrclvl',nt
          write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
          call wrtlvl(ACC_UTILLVL(1),LVL_TRC(iogrp),
     .    1.,0.,cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
        enddo
#  else
        do nt=1,ntrocn
          call inilvl(ACC_UTILLVL(1),'p',0.)
          do k=1,kk
            call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
            call acclvl(ACC_UTILLVL,trc(1-nbdy,1-nbdy,k1m,nt),'p',
     .        k,ind1,ind2,wghts)
          enddo
          call msklvl(ACC_UTILLVL(1),'p')
          write (trcnm,'(a,i3.3)') 'trclvl',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call wrtlvl(ACC_UTILLVL(1),LVL_TRC(iogrp),
     .    1.,0.,cmpflg,ip,'p',trim(trcnm),trim(trcnml),' ',' ')
        enddo
#  endif
      endif
#endif
c
c --- close netcdf file
      call ncfcls
c
c --- initialisation of fields
      call inifld(iogrp)
c
c --- reset accumulation counter 
      nacc_phy(iogrp)=0
c
      end subroutine diaout



      subroutine diasec(iogrp)
c
      implicit none
c
#include "common_blocks.h"

      integer :: iogrp
c
      integer :: n,i,j,k,s,l
      integer :: iostatus
      integer :: nsi(max_sec),isi(max_sec,sdm),jsi(max_sec,sdm),
     .           usi(max_sec,sdm),vsi(max_sec,sdm),equat_sec
      character(len=120) :: char120
      logical :: iniflg=.true.
      save nsi,isi,jsi,usi,vsi,iniflg,equat_sec
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .       uflx_cum,vflx_cum,uflx_cum350,vflx_cum350
      real, dimension(itdm,jtdm) :: uflx_cumt,vflx_cumt
     .       ,uflx_cum350t,vflx_cum350t
      real*8 :: volu,volv
c
c --- ------------------------------------------------------------------
c --- read section information
c --- ------------------------------------------------------------------
      if (iniflg) then
        if (mnproc.eq.1) then
          equat_sec=-1
          open(nfu,file=sec_sifile,status='old')
          sec_num=0
          do
            read(nfu,'(a120)',iostat=iostatus) char120
            if (iostatus.ne.0) exit
            if (char120(1:4).eq.'Name') then
              sec_num=sec_num+1
              sec_name(sec_num)=char120(7:120)
              if (index(char120,'equatorial_undercurrent').gt.0) then
                equat_sec=sec_num
              endif
              nsi(sec_num)=0
            else
              nsi(sec_num)=nsi(sec_num)+1
              read(char120,*) isi(sec_num,nsi(sec_num)),
     .                        jsi(sec_num,nsi(sec_num)),
     .                        usi(sec_num,nsi(sec_num)),
     .                        vsi(sec_num,nsi(sec_num))
            endif
          enddo
          close(nfu)
          write(lp,*) 'number of sections=',sec_num
        endif
        call xcbcst(sec_num)
        iniflg=.false.
      endif
c      
c --- Prepare 2d field   
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          uflx_cum(i,j)=0.
          vflx_cum(i,j)=0.
          uflx_cum350(i,j)=0.
          vflx_cum350(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Compute accumulated transports
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1,jj
        do k=1,ddm
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              uflx_cum(i,j)=uflx_cum(i,j)+
     .          phylvl(i,j,k,ACC_UFLXLVL(iogrp))
     .          *0.5e-3/(g*baclin*nacc_phy(iogrp))
            enddo
          enddo
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vflx_cum(i,j)=vflx_cum(i,j)+
     .          phylvl(i,j,k,ACC_VFLXLVL(iogrp))
     .          *0.5e-3/(g*baclin*nacc_phy(iogrp))
            enddo
          enddo
c
c --- -- the upper 350 m  for equatorial_undercurrent
          if (k.eq.k350-1) then
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                uflx_cum350(i,j)=uflx_cum(i,j)
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vflx_cum350(i,j)=vflx_cum(i,j)
              enddo
            enddo
          elseif (k.eq.k350) then 
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                uflx_cum350(i,j)=uflx_cum350(i,j)+
     .            w350*(uflx_cum(i,j)-uflx_cum350(i,j))
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vflx_cum350(i,j)=vflx_cum350(i,j)+
     .            w350*(vflx_cum(i,j)-vflx_cum350(i,j))
              enddo
            enddo
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Collect data on master node
      call xcaget(uflx_cumt,uflx_cum,1)
      call xcaget(vflx_cumt,vflx_cum,1)
      call xcaget(uflx_cum350t,uflx_cum350,1)
      call xcaget(vflx_cum350t,vflx_cum350,1)
c
c --- Compute section transports
      if (mnproc.eq.1) then
        do s=1,sec_num
          voltr(s)=0.
          if (s.eq.equat_sec) then
            do n=1,nsi(s)
              i=isi(s,n)
              j=jsi(s,n)
              volu=uflx_cum350t(i,j)*real(usi(s,n))
              volv=vflx_cum350t(i,j)*real(vsi(s,n))
              voltr(s)=voltr(s)+volu+volv
            enddo
          else
            do n=1,nsi(s)
              i=isi(s,n)
              j=jsi(s,n)
              volu=uflx_cumt(i,j)*real(usi(s,n))
              volv=vflx_cumt(i,j)*real(vsi(s,n))
              voltr(s)=voltr(s)+volu+volv
            enddo
          endif
        enddo
      endif
c
      end subroutine diasec



      subroutine diamer(iogrp)
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c
      integer :: iogrp
c 
      integer :: ncid,dimid,varid,i,j,k,l,m,n,o,s,ocn_nreg,iostatus,
     .  istat,iind1,jind1,uflg1,vflg1,nind1,
     .  nfld,ACC_UIND,ACC_VIND,nind(ldm),iind(sdm,ldm),jind(sdm,ldm),
     .  kmxl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer(kind=2) :: uflg(sdm,ldm),vflg(sdm,ldm),oflg(itdm,jtdm)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ucum,vcum
      real, dimension(itdm,jtdm) :: depthst,ucumg,vcumg
      real, save, allocatable, dimension(:,:) ::
     .  mflx_or,mflx_mr,mflx_last_mr
      integer, save, allocatable, dimension(:,:) :: mcnt_or,mcnt_mr,kmax
      real :: r
c
      character :: c20*20
c
      logical :: iniflg=.true.
c
      save nind,iind,jind,oflg,uflg,vflg,depthst,iniflg,ocn_nreg
c
      if (iniflg) then
c
        if (mnproc.eq.1) then
c
c --- --- Read ocean region flags from ocean_regions.nc
          call ncerro(nf90_open(mer_orfile,nf90_nowrite,ncid))
          call ncerro(nf90_inq_dimid(ncid,'x',dimid))
          call ncerro(nf90_inquire_dimension(ncid,dimid,len=i))
          call ncerro(nf90_inq_dimid(ncid,'y',dimid))
          call ncerro(nf90_inquire_dimension(ncid,dimid,len=j))
          if (i.ne.itdm.or.j.ne.jtdm) then
            write (lp,'(2a)') ' wrong dimensions in',trim(mer_orfile)
            call xchalt('(diamer)')
                   stop '(diamer)'
          endif
          call ncerro(nf90_inq_dimid(ncid,'regions',dimid))
          call ncerro(nf90_inquire_dimension(ncid,dimid,len=ocn_nreg))
          call ncerro(nf90_inq_varid(ncid,'region',varid))
          call ncerro(nf90_get_var(ncid,varid,oflg))
          call ncerro(nf90_close(ncid))
c
c --- --- Read section file metra_index.dat
          open(nfu,file=mer_mifile,status='old')
          lmax=0
          do l=1,ldm
            c20=' '
            read(nfu,'(a)',iostat=iostatus) c20
            if (iostatus.ne.0) exit
            if (c20(1:7).eq.'Section') then
              lmax=l
              read(c20(9:15),'(f7.3)') mtlat(l)
              read(c20(16:20),'(i7)') nind1
            else
              write(lp,'(2a)') ' problem in ',trim(mer_mifile)
              call xchalt('(diamer)')
                     stop '(diamer)'
            endif
            nind(l)=0
            do s=1,nind1
              read(nfu,'(a)',iostat=iostatus) c20
              if (iostatus.ne.0) then
                write(lp,*) 'section too short?'
                call xchalt('(diamer)')
                       stop '(diamer)'
              endif
              read(c20,*) iind1,jind1,uflg1,vflg1
              if (oflg(iind1,jind1).gt.0) then
                nind(l)=nind(l)+1
                iind(nind(l),l)=iind1
                jind(nind(l),l)=jind1
                uflg(nind(l),l)=uflg1
                vflg(nind(l),l)=vflg1
                if (iind(nind(l),l).lt.1.or.
     .              iind(nind(l),l).gt.itdm.or.
     .              jind(nind(l),l).lt.1.or.
     .              jind(nind(l),l).gt.jtdm) then
                  write(lp,*) 'iind=',iind(nind(l),l),' itdm=',itdm
                  write(lp,*) 'jind=',jind(nind(l),l),' jtdm=',jtdm
                  call flush(lp)
                  write(lp,*) 'indices out of range!'
                  call xchalt('(diamer)')
                         stop '(diamer)'
                endif
              endif
            enddo
            if (lmax.gt.ldm) then
              write(lp,*) 'too many or too long sections!'
              call xchalt('(diamer)')
                     stop '(diamer)'
            endif
          enddo
          close(nfu)
c
c --- --- Allocate utility arrays for meridional fluxes
          allocate(mflx_or(lmax,ocn_nreg),mflx_mr(lmax,mer_nreg),
     .             mflx_last_mr(lmax,mer_nreg),
     .             mcnt_or(lmax,ocn_nreg),mcnt_mr(lmax,mer_nreg),
     .             kmax(lmax,mer_nreg),
     .             stat=istat)
          if (istat.ne.0) then
             write (lp,*) 'Cannot allocate enough memory!'
             call xchalt('(diamer)')
                    stop '(diamer)'
          endif
c
        endif
c
        call xcbcst(lmax)
c
c --- - Allocate arrays for meridional fluxes
        allocate(mmflxl(lmax,kdm,mer_nreg),mmftdl(lmax,kdm,mer_nreg),
     .           mmflxd(lmax,ddm,mer_nreg),mmftdd(lmax,ddm,mer_nreg),
     .           mhflx(lmax,mer_nreg),mhftd(lmax,mer_nreg),
     .           mhfld(lmax,mer_nreg),msflx(lmax,mer_nreg),
     .           msftd(lmax,mer_nreg),msfld(lmax,mer_nreg),
     .           stat=istat)
        if (istat.ne.0) then
           write (lp,*) 'Cannot allocate enough memory!'
           call xchalt('(diamer)')
                  stop '(diamer)'
        endif
c
      endif
c
c --- Compute vertical integrated heat and salt transports 
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          ucum(i,j)=0.
          vcum(i,j)=0.
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do nfld=1,6
c
        if     (nfld.eq.1) then
          if (ACC_MHFLX(iogrp).eq.0) exit
          ACC_UIND=ACC_UTFLX(iogrp)
          ACC_VIND=ACC_VTFLX(iogrp)
          r=spcifh*0.5/(g*baclin*nacc_phy(iogrp))
        elseif (nfld.eq.2) then
          if (ACC_MHFTD(iogrp).eq.0) exit
          ACC_UIND=ACC_UTFLTD(iogrp)
          ACC_VIND=ACC_VTFLTD(iogrp)
          r=spcifh*0.5/(g*baclin*nacc_phy(iogrp))
        elseif (nfld.eq.3) then
          if (ACC_MHFLD(iogrp).eq.0) exit
          ACC_UIND=ACC_UTFLLD(iogrp)
          ACC_VIND=ACC_VTFLLD(iogrp)
          r=spcifh*0.5/(g*baclin*nacc_phy(iogrp))
        elseif (nfld.eq.4) then
          if (ACC_MSFLX(iogrp).eq.0) exit
          ACC_UIND=ACC_USFLX(iogrp)
          ACC_VIND=ACC_VSFLX(iogrp)
          r=0.5e-6/(g*baclin*nacc_phy(iogrp))
        elseif (nfld.eq.5) then
          if (ACC_MSFTD(iogrp).eq.0) exit
          ACC_UIND=ACC_USFLTD(iogrp)
          ACC_VIND=ACC_VSFLTD(iogrp)
          r=0.5e-6/(g*baclin*nacc_phy(iogrp))
        elseif (nfld.eq.6) then
          if (ACC_MSFLD(iogrp).eq.0) exit
          ACC_UIND=ACC_USFLLD(iogrp)
          ACC_VIND=ACC_VSFLLD(iogrp)
          r=0.5e-6/(g*baclin*nacc_phy(iogrp))
        else
          write(lp,*) 'field index out of range'
          call xchalt('(diamer)')
                 stop '(diamer)'
        endif
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              ucum(i,j)=phylyr(i,j,1,ACC_UIND)*r
            enddo
          enddo
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vcum(i,j)=phylyr(i,j,1,ACC_VIND)*r
            enddo
          enddo
          do k=2,kk
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                ucum(i,j)=ucum(i,j)+phylyr(i,j,k,ACC_UIND)*r
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vcum(i,j)=vcum(i,j)+phylyr(i,j,k,ACC_VIND)*r
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        call xcaget(ucumg,ucum,1)
        call xcaget(vcumg,vcum,1)
        if (mnproc.eq.1) then
          do l=1,lmax
c ---- ---- Accumulate meridional fluxes in seperate ocean regions
            do o=1,ocn_nreg
              mflx_or(l,o)=0.
              mcnt_or(l,o)=0
            enddo 
            do s=1,nind(l) 
              i=iind(s,l)
              j=jind(s,l)
              o=oflg(i,j)
              mflx_or(l,o)=mflx_or(l,o)
     .                    +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
              mcnt_or(l,o)=mcnt_or(l,o)+1
            enddo
c --- ----- Add together the ocean regions that belong to the regions of
c --- ----- meridional flux diagnostics and apply a fill value outside
c --- ----- the regions latitude bounds or when no values have been
c --- ----- accumulated
            do m=1,mer_nreg
              if (mtlat(l).lt.mer_minlat(m).or.
     .            mtlat(l).gt.mer_maxlat(m)) then
                mflx_mr(l,m)=nf90_fill_double
              else
                mflx_mr(l,m)=0.
                mcnt_mr(l,m)=0
                if (mer_regflg(m,1).eq.0) then
                  do o=1,ocn_nreg
                    mflx_mr(l,m)=mflx_mr(l,m)+mflx_or(l,o)
                    mcnt_mr(l,m)=mcnt_mr(l,m)+mcnt_or(l,o)
                  enddo
                else
                  do n=1,mer_nflg(m)
                    o=mer_regflg(m,n)
                    mflx_mr(l,m)=mflx_mr(l,m)+mflx_or(l,o)
                    mcnt_mr(l,m)=mcnt_mr(l,m)+mcnt_or(l,o)
                  enddo
                endif
                if (mcnt_mr(l,m).eq.0) mflx_mr(l,m)=nf90_fill_double
              endif
            enddo
          enddo
c
          if     (nfld.eq.1) then
            do l=1,lmax
              do m=1,mer_nreg
                mhflx(l,m)=mflx_mr(l,m)
              enddo
            enddo
          elseif (nfld.eq.2) then
            do l=1,lmax
              do m=1,mer_nreg
                mhftd(l,m)=mflx_mr(l,m)
              enddo
            enddo
          elseif (nfld.eq.3) then
            do l=1,lmax
              do m=1,mer_nreg
                mhfld(l,m)=mflx_mr(l,m)
              enddo
            enddo
          elseif (nfld.eq.4) then
            do l=1,lmax
              do m=1,mer_nreg
                msflx(l,m)=mflx_mr(l,m)
              enddo
            enddo
          elseif (nfld.eq.5) then
            do l=1,lmax
              do m=1,mer_nreg
                msftd(l,m)=mflx_mr(l,m)
              enddo
            enddo
          elseif (nfld.eq.6) then
            do l=1,lmax
              do m=1,mer_nreg
                msfld(l,m)=mflx_mr(l,m)
              enddo
            enddo
          else
            write(lp,*) 'field index out of range'
            call xchalt('(diamer)')
                   stop '(diamer)'
          endif
c
        endif
      enddo
c
c --- Compute overturning stream function at isopycnic layer interfaces
c
      r=1./real(nacc_phy(iogrp))
      kmxl=0
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          k=1
          do while (phyh2d(i,j,ACC_SIGMX(IOGRP))*r.ge.sigmar(i,j,k))
            k=k+1
            if (k.eq.kk) exit
          enddo
          kmxl(i,j)=k
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kmxl(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmxl(i,j)=nint(util1(i,j))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      r=0.5e-3/(g*baclin*nacc_phy(iogrp))
c
      do nfld=1,2
c
        if     (nfld.eq.1) then
          if (ACC_MMFLXL(iogrp).eq.0) exit
          ACC_UIND=ACC_UFLX(iogrp)
          ACC_VIND=ACC_VFLX(iogrp)
        elseif (nfld.eq.2) then
          if (ACC_MMFTDL(iogrp).eq.0) exit
          ACC_UIND=ACC_UMFLTD(iogrp)
          ACC_VIND=ACC_VMFLTD(iogrp)
        else
          write(lp,*) 'field index out of range'
          call xchalt('(diamer)')
                 stop '(diamer)'
        endif
c
        if (mnproc.eq.1) then
          do l=1,lmax
            do m=1,mer_nreg
              mflx_last_mr(l,m)=0.
            enddo
          enddo
        endif
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              ucum(i,j)=0.
            enddo
          enddo
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vcum(i,j)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        do k=1,kk
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                if (k.gt.2) then
                  ucum(i,j)=ucum(i,j)+r*phylyr(i,j,k,ACC_UIND)
                endif
                if (k.eq.min(kmxl(i-1,j),kmxl(i,j))) then
                  ucum(i,j)=ucum(i,j)+r*(phylyr(i,j,1,ACC_UIND)
     .                                  +phylyr(i,j,2,ACC_UIND))
                endif
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                if (k.gt.2) then
                  vcum(i,j)=vcum(i,j)+r*phylyr(i,j,k,ACC_VIND)
                endif
                if (k.eq.min(kmxl(i,j-1),kmxl(i,j))) then
                  vcum(i,j)=vcum(i,j)+r*(phylyr(i,j,1,ACC_VIND)
     .                                  +phylyr(i,j,2,ACC_VIND))
                endif
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
          call xcaget(ucumg,ucum,1)
          call xcaget(vcumg,vcum,1)
          if (mnproc.eq.1) then
            do l=1,lmax
c ---- ------ Accumulate meridional fluxes in seperate ocean regions
              do o=1,ocn_nreg 
                mflx_or(l,o)=0.
              enddo 
              do s=1,nind(l) 
                i=iind(s,l)
                j=jind(s,l)
                o=oflg(i,j)
                mflx_or(l,o)=mflx_or(l,o)
     .                      +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
              enddo
c --- ------- Add together the ocean regions that belong to the regions
c --- ------- of meridional flux diagnostics and apply a fill value
c --- ------- outside the regions latitude bounds or when no values have
c --- ------- been accumulated
              do m=1,mer_nreg
                if (mtlat(l).lt.mer_minlat(m).or.
     .              mtlat(l).gt.mer_maxlat(m)) then
                  mflx_mr(l,m)=nf90_fill_double
                else
                  mflx_mr(l,m)=0.
                  if (mer_regflg(m,1).eq.0) then
                    do o=1,ocn_nreg
                      mflx_mr(l,m)=mflx_mr(l,m)+mflx_or(l,o)
                    enddo
                  else
                    do n=1,mer_nflg(m)
                      o=mer_regflg(m,n)
                      mflx_mr(l,m)=mflx_mr(l,m)+mflx_or(l,o)
                    enddo
                  endif
                  if (abs(mflx_mr(l,m)-mflx_last_mr(l,m)).lt.
     .                1.e5*epsil) then
                    mflx_last_mr(l,m)=mflx_mr(l,m)
                    mflx_mr(l,m)=nf90_fill_double
                  else
                    mflx_last_mr(l,m)=mflx_mr(l,m)
                  endif
                endif
              enddo
            enddo
c
            if     (nfld.eq.1) then
              do l=1,lmax
                do m=1,mer_nreg
                  mmflxl(l,k,m)=mflx_mr(l,m)
                enddo
              enddo
            elseif (nfld.eq.2) then
              do l=1,lmax
                do m=1,mer_nreg
                  mmftdl(l,k,m)=mflx_mr(l,m)
                enddo
              enddo
            else
              write(lp,*) 'field index out of range'
              call xchalt('(diamer)')
                     stop '(diamer)'
            endif
c
          endif
        enddo
      enddo
c
c --- Compute overturning stream function at levitus level interfaces
c --- Prepare depth mask 
c
      if (iniflg) call xcaget(depthst,depths,1)
      if (iniflg.and.mnproc.eq.1) then
        do l=1,lmax
          do m=1,mer_nreg
            kmax(l,m)=0
          enddo
        enddo 
        do k=1,ddm
          do l=1,lmax
            do s=1,nind(l)
              i=iind(s,l)
              j=jind(s,l)
              if (depthslev_bnds(1,k).lt.depthst(i,j)) then
                do m=1,mer_nreg
                  if (mer_regflg(m,1).eq.0) then
                    kmax(l,m)=k
                  else
                    do n=1,mer_nflg(m)
                      if (mer_regflg(m,n).eq.oflg(i,j)) kmax(l,m)=k
                    enddo
                  endif
                enddo
              endif
            enddo 
          enddo
        enddo
      endif
c
      r=0.5e-3/(g*baclin*nacc_phy(iogrp))
c
      do nfld=1,2
c
        if     (nfld.eq.1) then
          if (ACC_MMFLXD(iogrp).eq.0) exit
          ACC_UIND=ACC_UFLXLVL(iogrp)
          ACC_VIND=ACC_VFLXLVL(iogrp)
        elseif (nfld.eq.2) then
          if (ACC_MMFTDD(iogrp).eq.0) exit
          ACC_UIND=ACC_UMFLTDLVL(iogrp)
          ACC_VIND=ACC_VMFLTDLVL(iogrp)
        else
          write(lp,*) 'field index out of range'
          call xchalt('(diamer)')
                 stop '(diamer)'
        endif
c
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              ucum(i,j)=0.
            enddo
          enddo
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              vcum(i,j)=0.
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        do k=1,ddm
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                ucum(i,j)=ucum(i,j)+r*phylvl(i,j,k,ACC_UIND)
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vcum(i,j)=vcum(i,j)+r*phylvl(i,j,k,ACC_VIND)
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
          call xcaget(ucumg,ucum,1)
          call xcaget(vcumg,vcum,1)
          if (mnproc.eq.1) then
            do l=1,lmax
c ---- ------ Accumulate meridional fluxes in seperate ocean regions
              do o=1,ocn_nreg 
                mflx_or(l,o)=0.
              enddo 
              do s=1,nind(l) 
                i=iind(s,l)
                j=jind(s,l)
                o=oflg(i,j)
                mflx_or(l,o)=mflx_or(l,o)
     .                      +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
              enddo
c --- ------- Add together the ocean regions that belong to the regions
c --- ------- of meridional flux diagnostics and apply a fill value
c --- ------- outside the regions latitude and depth bounds
              do m=1,mer_nreg
                if (mtlat(l).lt.mer_minlat(m).or.
     .              mtlat(l).gt.mer_maxlat(m).or.
     .              kmax(l,m).lt.k) then
                  mflx_mr(l,m)=nf90_fill_double
                else
                  mflx_mr(l,m)=0.
                  if (mer_regflg(m,1).eq.0) then
                    do o=1,ocn_nreg
                      mflx_mr(l,m)=mflx_mr(l,m)+mflx_or(l,o)
                    enddo
                  else
                    do n=1,mer_nflg(m)
                      o=mer_regflg(m,n)
                      mflx_mr(l,m)=mflx_mr(l,m)+mflx_or(l,o)
                    enddo
                  endif
                endif
              enddo
            enddo
c
            if     (nfld.eq.1) then
              do l=1,lmax
                do m=1,mer_nreg
                  mmflxd(l,k,m)=mflx_mr(l,m)
                enddo
              enddo
            elseif (nfld.eq.2) then
              do l=1,lmax
                do m=1,mer_nreg
                  mmftdd(l,k,m)=mflx_mr(l,m)
                enddo
              enddo
            else
              write(lp,*) 'field index out of range'
              call xchalt('(diamer)')
                     stop '(diamer)'
            endif
c
          endif
        enddo
      enddo
c
      if (iniflg) iniflg=.false.
c
      end subroutine diamer



      subroutine diavfl(iogrp,m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- computation of vertical mass flux at isopycnic layer interfaces  
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer iogrp,m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: wflx,ucum,vcum
      integer :: i,j,k,km,kn,l
      real :: q
c
c ---
c --- Compute vertical mass flux at isopycnic layer interfaces 
c ---
      if (ACC_WFLX(iogrp)+ACC_WFLX2(IOGRP).ne.0) then
c
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            wflx(i,j)=0.
          enddo
        enddo  
c$OMP END PARALLEL DO
        do k=kk,1,-1
          km=k+mm
          kn=k+nn
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj  
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l)) 
                wflx(i,j)=wflx(i,j)
     .                   -(uflx(i+1,j,kn)-uflx(i,j,kn)
     .                    +vflx(i,j+1,kn)-vflx(i,j,kn))
     .                   -(dp(i,j,km)-dpold(i,j,km))*scp2(i,j)
                phylyr(i,j,k,ACC_WFLX(iogrp))=
     .            phylyr(i,j,k,ACC_WFLX(iogrp))+wflx(i,j)
                phylyr(i,j,k,ACC_WFLX2(iogrp))=
     .            phylyr(i,j,k,ACC_WFLX2(iogrp))+wflx(i,j)**2
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        enddo
      endif
c 
c --- Computation of vertical mass flux at levitus layer interfaces 
      if (ACC_WFLXLVL(iogrp)+ACC_WFLX2LVL(iogrp).ne.0) then
c
        call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_UFLXLVL(iogrp)),
     .              1,ddm, 1,1, halo_uv)
        call xctilr(phylvl(1-nbdy,1-nbdy,1,ACC_VFLXLVL(iogrp)),
     .              1,ddm, 1,1, halo_vv)
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj+1
          do i=1,ii+1
            ucum(i,j)=0.
            vcum(i,j)=0.
          enddo
        enddo  
c$OMP END PARALLEL DO
        do k=ddm,1,-1
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l)) 
                ucum(i,j)=ucum(i,j)
     .                   +phylvl(i,j,k,ACC_UFLXLVL(iogrp))
     .                   -phylvl(i,j,k,ACC_UFLXOLD(iogrp))
              enddo
            enddo
          enddo  
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj+1 
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l)) 
                vcum(i,j)=vcum(i,j)
     .                   +phylvl(i,j,k,ACC_VFLXLVL(iogrp))
     .                   -phylvl(i,j,k,ACC_VFLXOLD(iogrp))
              enddo
            enddo
          enddo  
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i,q)
          do j=1,jj  
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l)) 
                q=-(ucum(i+1,j)-ucum(i,j)
     .             +vcum(i,j+1)-vcum(i,j))
                phylvl(i,j,k,ACC_WFLXLVL(iogrp))=
     .            phylvl(i,j,k,ACC_WFLXLVL(iogrp))+q
                phylvl(i,j,k,ACC_WFLX2LVL(iogrp))=
     .            phylvl(i,j,k,ACC_WFLX2LVL(iogrp))+q**2
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        enddo
      endif
c
      end subroutine diavfl



      subroutine diazlv(gridid,k,mm,nn,ind1,ind2,weights,weightsflx)
c-----------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
c
      integer :: d,i,j,k,mm,nn,l,kl,km,kn,kml,k1m
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)  :: ind1,ind2
c 
      real :: r,dzeps,dpeps,flxeps
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk)  :: ztop,zbot
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: weights,
     .  weightsflx,dlevp,dlevu,dlevv
c
      character :: gridid
c
      logical :: iniflg=.true.
c
      save ztop,zbot,dlevp,dlevu,dlevv,iniflg
c
c --- Define thresholds 
      dzeps=1e1*epsil
      dpeps=1e5*epsil
      flxeps=1e5*epsil
c
c --- Sort out stuff related to time stepping 
      km=k+mm
      kn=k+nn
      k1m=1+mm
c
c --- Adjust bounds of levitus levels according to model bathymetry
      if (iniflg) then
c$OMP PARALLEL DO PRIVATE(d,l,i)
        do j=1,jj+1
          do d=1,ddm
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                 dlevp(i,j,d)=max(dzeps,min(pbath(i,j),
     .             depthslev_bnds(2,d))-depthslev_bnds(1,d))
              enddo
            enddo
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                 dlevu(i,j,d)=max(dzeps,min(ubath(i,j),
     .             depthslev_bnds(2,d))-depthslev_bnds(1,d))
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                 dlevv(i,j,d)=max(dzeps,min(vbath(i,j),
     .             depthslev_bnds(2,d))-depthslev_bnds(1,d))
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        iniflg=.false.
      endif
c
c
c --- Compute top and bottom depths of density layers 
      if (k.eq.1) then
        if (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i,kl,kml)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,1)=dp(i,j,k1m)
              enddo
            enddo
            do kl=2,kk
              kml=kl+mm
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  zbot(i,j,kl)=zbot(i,j,kl-1)+dp(i,j,kml)
                enddo
              enddo
            enddo
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,1)=zbot(i,j,1)*pbath(i,j)/zbot(i,j,kk)
                ztop(i,j,1)=0.
                ind1(i,j)=1
              enddo
            enddo
            do kl=2,kk
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  zbot(i,j,kl)=zbot(i,j,kl)*pbath(i,j)/zbot(i,j,kk)
                  ztop(i,j,kl)=zbot(i,j,kl-1)
                enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i,kl,kml)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                zbot(i,j,1)=max(dpu(i,j,k1m),dpeps)
              enddo
            enddo
            do kl=2,kk
              kml=kl+mm
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                  zbot(i,j,kl)=zbot(i,j,kl-1)+max(dpu(i,j,kml),dpeps)
                enddo
              enddo
            enddo
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                zbot(i,j,1)=zbot(i,j,1)*ubath(i,j)/zbot(i,j,kk)
                ztop(i,j,1)=0.
                ind1(i,j)=1
              enddo
            enddo
            do kl=2,kk
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                  zbot(i,j,kl)=zbot(i,j,kl)*ubath(i,j)/zbot(i,j,kk)
                  ztop(i,j,kl)=zbot(i,j,kl-1)
                enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i,kl,kml)
          do j=1,jj+1
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                zbot(i,j,1)=max(dpv(i,j,k1m),dpeps)
              enddo
            enddo
            do kl=2,kk
              kml=kl+mm
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                  zbot(i,j,kl)=zbot(i,j,kl-1)+max(dpv(i,j,kml),dpeps)
                enddo
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                zbot(i,j,1)=zbot(i,j,1)*vbath(i,j)/zbot(i,j,kk)
                ztop(i,j,1)=0.
                ind1(i,j)=1
              enddo
            enddo
            do kl=2,kk
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                  zbot(i,j,kl)=zbot(i,j,kl)*vbath(i,j)/zbot(i,j,kk)
                  ztop(i,j,kl)=zbot(i,j,kl-1)
                enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          write (lp,'(a)') 'cannot identify grid'
          flush(lp)
          call xchalt('(diazlv)')
                 stop '(diazlv)'
        endif
      endif
c
c --- Compute interpolation weights 
      if (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i,d)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              ind2(i,j)=0
              if (dp(i,j,km).gt.dpeps) then
                do d=ind1(i,j),ddm
                  if (depthslev_bnds(2,d).le.ztop(i,j,k)) then
                    ind1(i,j)=d+1
                    cycle
                  elseif (depthslev_bnds(1,d).ge.zbot(i,j,k)) then
                    exit
                  endif
                  ind2(i,j)=d
                  weights(i,j,d)=(min(zbot(i,j,k),
     .              depthslev_bnds(2,d))-max(ztop(i,j,k),
     .              depthslev_bnds(1,d)))/dlevp(i,j,d)
                enddo
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i,d,r) 
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
              ind2(i,j)=0
              if (abs(uflx(i,j,kn)).gt.flxeps) then
                do d=ind1(i,j),ddm
                  if (depthslev_bnds(2,d).lt.ztop(i,j,k)) then
                    ind1(i,j)=d+1
                    cycle
                  elseif (depthslev_bnds(1,d).gt.zbot(i,j,k)) then
                    exit
                  endif
                  ind2(i,j)=d
                  r=(min(zbot(i,j,k),depthslev_bnds(2,d))-
     .              max(ztop(i,j,k),depthslev_bnds(1,d)))
                  weights(i,j,d)=r/dlevu(i,j,d)
                  weightsflx(i,j,d)=r/(zbot(i,j,k)-ztop(i,j,k))
                enddo
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i,d,r)
        do j=1,jj+1
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              ind2(i,j)=0
              if (abs(vflx(i,j,kn)).gt.flxeps) then
                do d=ind1(i,j),ddm
                  if (depthslev_bnds(2,d).lt.ztop(i,j,k)) then
                    ind1(i,j)=d+1
                    cycle
                  elseif (depthslev_bnds(1,d).gt.zbot(i,j,k)) then
                    exit
                  endif
                  ind2(i,j)=d
                  r=(min(zbot(i,j,k),depthslev_bnds(2,d))-
     .               max(ztop(i,j,k),depthslev_bnds(1,d)))
                  weights(i,j,d)=r/dlevv(i,j,d)
                  weightsflx(i,j,d)=r/(zbot(i,j,k)-ztop(i,j,k))
                enddo
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      else
        write (lp,'(a)') 'cannot identify grid'
        flush(lp)
        call xchalt('(diazlv)')
               stop '(diazlv)'
      endif
c
      end subroutine diazlv


c --- ------------------------------------------------------------------
c --- ------------------------------------------------------------------
c --- ------------------------------------------------------------------


      subroutine inih2d(pos,gridid,inival)
c
c --- ------------------------------------------------------------------
c --- Description: initialise 2d diagnostic field
c ---   
c --- Arguments:
c ---   int  pos      (in)     : position in common buffer  
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   real inival   (in)     : value used for initalisation
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos
      real :: inival
      character :: gridid
c 
      integer :: i,j
c
c --- Check whether field should be initialised
      if (pos.eq.0) return
c
      if (gridid(1:1).eq.'u') then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            phyh2d(i,j,pos)=inival*iu(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            phyh2d(i,j,pos)=inival*iv(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            phyh2d(i,j,pos)=inival*ip(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            phyh2d(i,j,pos)=inival
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      end subroutine inih2d



      subroutine inilyr(pos,gridid,inival)
c
c --- ------------------------------------------------------------------
c --- Description: initialise layer diagnostic field
c ---   
c --- Arguments:
c ---   int  pos      (in)     : position in common buffer  
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   real inival   (in)     : value used for initalisation
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos
      real :: inival
      character :: gridid
c 
      integer :: i,j,k
c
c --- Check whether field should be initialised
      if (pos.eq.0) return
c
      if (gridid(1:1).eq.'u') then
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,kk 
            do i=1-nbdy,ii+nbdy
              phylyr(i,j,k,pos)=inival*iu(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,kk 
            do i=1-nbdy,ii+nbdy
              phylyr(i,j,k,pos)=inival*iv(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,kk 
            do i=1-nbdy,ii+nbdy
              phylyr(i,j,k,pos)=inival*ip(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,kk 
            do i=1-nbdy,ii+nbdy
              phylyr(i,j,k,pos)=inival
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      end subroutine inilyr



      subroutine inilvl(pos,gridid,inival)
c
c --- ------------------------------------------------------------------
c --- Description: initialise level diagnostic field
c ---   
c --- Arguments:
c ---   int  pos      (in)     : position in common buffer  
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   real inival   (in)     : value used for initalisation
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos
      real :: inival
      character :: gridid
c 
      integer :: i,j,k
c
c --- Check whether field should be initialised
      if (pos.eq.0) return
c
      if (gridid(1:1).eq.'u') then
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,ddm
            do i=1-nbdy,ii+nbdy
              phylvl(i,j,k,pos)=inival*iu(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,ddm
            do i=1-nbdy,ii+nbdy
              phylvl(i,j,k,pos)=inival*iv(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,ddm
            do i=1-nbdy,ii+nbdy
              phylvl(i,j,k,pos)=inival*ip(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO PRIVATE(k,i)
        do j=1-nbdy,jj+nbdy
          do k=1,ddm
            do i=1-nbdy,ii+nbdy
              phylvl(i,j,k,pos)=inival
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      end subroutine inilvl



      subroutine acch2d(pos,fld,wghts,wghtsflg,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: accumulate 2d fields 
c ---  
c --- Arguments: 
c ---   int  pos      (in)     : position in 2d buffer  
c ---   real fld      (in)     : input data used for accumulation
c ---   real wghts    (in)     : weights used for accumulation
c ---   int  wghtsflg (in)     : weights flag (0=no weighting)  
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos(nphymax),wghtsflg
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld,wghts
      character :: gridid
c 
      integer :: i,j,l,o
c
c --- Check whether field should be accumulated
      do o=1,nphy
        if (pos(o).eq.0) cycle
c
        if (gridid.eq.'u') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                  phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                  phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)*
     .              wghts(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        elseif (gridid.eq.'v') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                  phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                  phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)*
     .              wghts(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        elseif (gridid.eq.'p') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)*
     .              wghts(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(acch2d)')
                 stop '(acch2d)'
        endif
c
      enddo
c   
      end subroutine acch2d



      subroutine maxh2d(pos,fld,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: store maximum of 2d fields 
c ---  
c --- Arguments: 
c ---   int  pos      (in)     : position in 2d buffer  
c ---   real fld      (in)     : input data used for finding maximum
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos(nphymax)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
      character :: gridid
c 
      integer :: i,j,l,o
c
c --- Check whether maximum of field should be stored
      do o=1,nphy
        if (pos(o).eq.0) cycle
c
        if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                phyh2d(i,j,pos(o))=max(phyh2d(i,j,pos(o)),fld(i,j))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                phyh2d(i,j,pos(o))=max(phyh2d(i,j,pos(o)),fld(i,j))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                phyh2d(i,j,pos(o))=max(phyh2d(i,j,pos(o)),fld(i,j))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(maxh2d)')
                 stop '(maxh2d)'
        endif
c
      enddo
c   
      end subroutine maxh2d



      subroutine minh2d(pos,fld,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: store minimum of 2d fields 
c ---  
c --- Arguments: 
c ---   int  pos      (in)     : position in 2d buffer  
c ---   real fld      (in)     : input data used for finding minimum
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos(nphymax)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
      character :: gridid
c 
      integer :: i,j,l,o
c
c --- Check whether minimum of field should be stored
      do o=1,nphy
        if (pos(o).eq.0) cycle
c
        if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                phyh2d(i,j,pos(o))=min(phyh2d(i,j,pos(o)),fld(i,j))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                phyh2d(i,j,pos(o))=min(phyh2d(i,j,pos(o)),fld(i,j))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                phyh2d(i,j,pos(o))=min(phyh2d(i,j,pos(o)),fld(i,j))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(minh2d)')
                 stop '(minh2d)'
        endif
c
      enddo
c   
      end subroutine minh2d



      subroutine sqh2d(pos,fld,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: accumulate square of 2d fields 
c ---  
c --- Arguments: 
c ---   int  pos      (in)     : position in 2d buffer  
c ---   real fld      (in)     : input data used for accumulation
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos(nphymax)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
      character :: gridid
c 
      integer :: i,j,l,o
c
c --- Check whether field should be accumulated
      do o=1,nphy
        if (pos(o).eq.0) cycle
c
        if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)*fld(i,j)
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)*fld(i,j)
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                phyh2d(i,j,pos(o))=phyh2d(i,j,pos(o))+fld(i,j)*fld(i,j)
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(minh2d)')
                 stop '(minh2d)'
        endif
c
      enddo
c   
      end subroutine sqh2d



      subroutine acclyr(pos,fld,wghts,wghtsflg,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: accumulate layer fields 
c ---  
c --- Arguments: 
c ---   int  pos      (in)     : position in 3d layer buffer  
c ---   real fld      (in)     : input data used for accumulation
c ---   real wghts    (in)     : weights used for accumulation
c ---   int  wghtsflg (in)     : weights flag (0=no weighting)  
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos(nphymax),wghtsflg
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: fld,wghts
      character :: gridid
c 
      integer :: i,j,k,l,o
c
c --- Check whether field should be accumulated
      do o=1,nphy
        if (pos(o).eq.0) cycle
c
        if (gridid.eq.'u') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO PRIVATE(k,l,i)
            do j=1,jj
              do k=1,kk
                do l=1,isu(j)
                  do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(k,l,i)
            do j=1,jj
              do k=1,kk
                do l=1,isu(j)
                  do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
     .                *wghts(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        elseif (gridid.eq.'v') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO PRIVATE(k,l,i)
            do j=1,jj
              do k=1,kk
                do l=1,isv(j)
                  do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(k,l,i)
            do j=1,jj
              do k=1,kk
                do l=1,isv(j)
                  do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
     .                *wghts(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        elseif (gridid.eq.'p') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO PRIVATE(k,l,i)
            do j=1,jj
              do k=1,kk
                do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO PRIVATE(k,l,i)
            do j=1,jj
              do k=1,kk
                do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
     .                *wghts(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(acclyr)')
                 stop '(acclyr)'
        endif
c
      enddo
c   
      end subroutine acclyr



      subroutine acclvl(pos,fld,gridid,k,ind1,ind2,wghts)
c
c --- ------------------------------------------------------------------
c --- Description: accumulate 3d level fields
c ---  
c --- Arguments: 
c ---   int  pos      (in)     : position in buffer  
c ---   real fld      (in)     : input data used for accumulation
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   int  k        (in)     : layer index of fld  
c ---   int  ind1     (in)     : index field for first accumulated level 
c ---   int  ind2     (in)     : index field for last accumulated level 
c ---   real wghts    (in)     : weights used for accumulation
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: pos(nphymax),k
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: fld
      character :: gridid
c 
      integer :: d,i,j,l,o
c
c --- Check whether field should be accumulated
      do o=1,nphy     
        if (pos(o).eq.0) cycle
c
        if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i,d)
          do j=1,jj
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                do d=ind1(i,j),ind2(i,j)
                  phylvl(i,j,d,pos(o))=phylvl(i,j,d,pos(o))+
     .              fld(i,j,k)*wghts(i,j,d)
                enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i,d)
          do j=1,jj+1
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                do d=ind1(i,j),ind2(i,j)
                  phylvl(i,j,d,pos(o))=phylvl(i,j,d,pos(o))+
     .              fld(i,j,k)*wghts(i,j,d)
                enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i,d)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                do d=ind1(i,j),ind2(i,j)
                  phylvl(i,j,d,pos(o))=phylvl(i,j,d,pos(o))+
     .              fld(i,j,k)*wghts(i,j,d)
                enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(acclvl)')
                 stop '(acclvl)'
        endif
      enddo
c     
      end subroutine acclvl



      subroutine inifld(iogrp)
c
      implicit none
c
#include "common_blocks.h"

      integer :: iogrp
c
c --- initialisation of 2d fields
      call inih2d(ACC_UB(iogrp),'u',0.)
      call inih2d(ACC_UBFLXS(iogrp),'u',0.)
      call inih2d(ACC_ZTX(iogrp),'u',0.)
      call inih2d(ACC_TAUX(iogrp),'u',0.)
      call inih2d(ACC_MXLU(iogrp),'u',0.)
      call inih2d(ACC_MLDU(iogrp),'u',0.)
      call inih2d(ACC_UICE(iogrp),'u',0.)
      call inih2d(ACC_IVOLU(iogrp),'u',0.)
c
      call inih2d(ACC_VB(iogrp),'v',0.)
      call inih2d(ACC_VBFLXS(iogrp),'v',0.)
      call inih2d(ACC_MTY(iogrp),'v',0.)
      call inih2d(ACC_TAUY(iogrp),'v',0.)
      call inih2d(ACC_MXLV(iogrp),'v',0.)
      call inih2d(ACC_MLDV(iogrp),'v',0.)
      call inih2d(ACC_VICE(iogrp),'v',0.)
      call inih2d(ACC_IVOLV(iogrp),'v',0.)
c
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
      call inih2d(ACC_FMLTFZ(iogrp),'p',0.)
      call inih2d(ACC_HMLTFZ(iogrp),'p',0.)
      call inih2d(ACC_HICE(iogrp),'p',0.)
      call inih2d(ACC_HSNW(iogrp),'p',0.)
      call inih2d(ACC_FICE(iogrp),'p',0.)
      call inih2d(ACC_TSRF(iogrp),'p',0.)
      call inih2d(ACC_TICE(iogrp),'p',0.)
c 
c --- initialisation of 3d layer fields
      call inilyr(ACC_UVEL(iogrp),'u',0.)
      call inilyr(ACC_DPU(iogrp),'u',0.)
      call inilyr(ACC_UFLX(iogrp),'u',0.)
      call inilyr(ACC_UTFLX(iogrp),'u',0.)
      call inilyr(ACC_USFLX(iogrp),'u',0.)
      call inilyr(ACC_UMFLTD(iogrp),'u',0.)
      call inilyr(ACC_UTFLTD(iogrp),'u',0.)
      call inilyr(ACC_UTFLLD(iogrp),'u',0.)
      call inilyr(ACC_USFLTD(iogrp),'u',0.)
      call inilyr(ACC_USFLLD(iogrp),'u',0.)
c
      call inilyr(ACC_VVEL(iogrp),'v',0.)
      call inilyr(ACC_DPV(iogrp),'v',0.)
      call inilyr(ACC_VFLX(iogrp),'v',0.)
      call inilyr(ACC_VTFLX(iogrp),'v',0.)
      call inilyr(ACC_VSFLX(iogrp),'v',0.)
      call inilyr(ACC_VMFLTD(iogrp),'v',0.)
      call inilyr(ACC_VTFLTD(iogrp),'v',0.)
      call inilyr(ACC_VTFLLD(iogrp),'v',0.)
      call inilyr(ACC_VSFLTD(iogrp),'v',0.)
      call inilyr(ACC_VSFLLD(iogrp),'v',0.)
c
      call inilyr(ACC_SALN(iogrp),'p',0.)
      call inilyr(ACC_TEMP(iogrp),'p',0.)
      call inilyr(ACC_DP(iogrp),'p',0.)
      call inilyr(ACC_DZ(iogrp),'p',0.)
      call inilyr(ACC_BFSQ(iogrp),'p',0.)
      call inilyr(ACC_DIFDIA(iogrp),'p',0.)
      call inilyr(ACC_DIFINT(iogrp),'p',0.)
      call inilyr(ACC_DIFISO(iogrp),'p',0.)
      call inilyr(ACC_WFLX(iogrp),'p',0.)
      call inilyr(ACC_WFLX2(iogrp),'p',0.)
      call inilyr(ACC_AVDSG(iogrp),'p',0.)
      call inilyr(ACC_DPVOR(iogrp),'p',0.)
#if defined(TRC) && defined(TKE)
      call inilyr(ACC_TKE(iogrp),'p',0.)
      call inilyr(ACC_GLS_PSI(iogrp),'p',0.)
#endif
c
c --- initialsation of 3d level fields
      call inilvl(ACC_UVELLVL(iogrp),'u',0.)
      call inilvl(ACC_UFLXLVL(iogrp),'u',0.)
      call inilvl(ACC_UTFLXLVL(iogrp),'u',0.)
      call inilvl(ACC_USFLXLVL(iogrp),'u',0.)
      call inilvl(ACC_UMFLTDLVL(iogrp),'u',0.)
      call inilvl(ACC_UTFLTDLVL(iogrp),'u',0.)
      call inilvl(ACC_UTFLLDLVL(iogrp),'u',0.)
      call inilvl(ACC_USFLTDLVL(iogrp),'u',0.)
      call inilvl(ACC_USFLLDLVL(iogrp),'u',0.)
c
      call inilvl(ACC_VVELLVL(iogrp),'v',0.)
      call inilvl(ACC_VFLXLVL(iogrp),'v',0.)
      call inilvl(ACC_VTFLXLVL(iogrp),'v',0.)
      call inilvl(ACC_VSFLXLVL(iogrp),'v',0.)
      call inilvl(ACC_VMFLTDLVL(iogrp),'v',0.)
      call inilvl(ACC_VTFLTDLVL(iogrp),'v',0.)
      call inilvl(ACC_VTFLLDLVL(iogrp),'v',0.)
      call inilvl(ACC_VSFLTDLVL(iogrp),'v',0.)
      call inilvl(ACC_VSFLLDLVL(iogrp),'v',0.)
c
      call inilvl(ACC_BFSQLVL(iogrp),'p',0.)
      call inilvl(ACC_DIFDIALVL(iogrp),'p',0.)
      call inilvl(ACC_DIFINTLVL(iogrp),'p',0.)
      call inilvl(ACC_DIFISOLVL(iogrp),'p',0.)
      call inilvl(ACC_DZLVL(iogrp),'p',0.)
      call inilvl(ACC_SALNLVL(iogrp),'p',0.)
      call inilvl(ACC_TEMPLVL(iogrp),'p',0.)
      call inilvl(ACC_WFLXLVL(iogrp),'p',0.)
      call inilvl(ACC_WFLX2LVL(iogrp),'p',0.)
      call inilvl(ACC_PVLVL(iogrp),'p',0.)
#if defined(TRC) && defined(TKE)
      call inilvl(ACC_TKELVL(iogrp),'p',0.)
      call inilvl(ACC_GLS_PSILVL(iogrp),'p',0.)
#endif
c
      end subroutine inifld
      


      subroutine finh2d(posacc,poswgt,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: finalise accumulation of weighted 2d fields 
c ---   
c --- Arguments:
c ---   real posacc   (in)     : position of accumulated field in buffer
c ---   real poswgt   (in)     : position of accumulated weights 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: posacc,poswgt
      character :: gridid
c 
      integer :: i,j,l
      real, parameter :: epsil=1e-11
c
c --- Check whether field should be initialised
      if (posacc.eq.0) return
c
      if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              phyh2d(i,j,posacc)=phyh2d(i,j,posacc)/
     .          max(epsil,phyh2d(i,j,poswgt))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              phyh2d(i,j,posacc)=phyh2d(i,j,posacc)/
     .          max(epsil,phyh2d(i,j,poswgt))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              phyh2d(i,j,posacc)=phyh2d(i,j,posacc)/
     .          max(epsil,phyh2d(i,j,poswgt))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(finh2d)')
               stop '(finh2d)'
      endif
c     
      end subroutine finh2d
     


      subroutine finlyr(posacc,poswgt,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: finalise accumulation of weighted 3d layer fields 
c ---   
c --- Arguments:
c ---   real posacc   (in)     : position of accumulated field in buffer
c ---   real poswgt   (in)     : position of accumulated weights 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c 
      integer :: posacc,poswgt
      character :: gridid
c 
      integer :: i,j,k,l
      real, parameter :: epsil=1e-11
c
c --- Check whether field should be initialised
      if (posacc.eq.0) return
c
      if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                if (phylyr(i,j,k,poswgt).gt.epsil) then
                  phylyr(i,j,k,posacc)=phylyr(i,j,k,posacc)/
     .              phylyr(i,j,k,poswgt) 
                else 
                  phylyr(i,j,k,posacc)=nf90_fill_double
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                if (phylyr(i,j,k,poswgt).gt.epsil) then
                  phylyr(i,j,k,posacc)=phylyr(i,j,k,posacc)/
     .              phylyr(i,j,k,poswgt) 
                else 
                  phylyr(i,j,k,posacc)=nf90_fill_double
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (phylyr(i,j,k,poswgt).gt.epsil) then
                  phylyr(i,j,k,posacc)=phylyr(i,j,k,posacc)/
     .              phylyr(i,j,k,poswgt) 
                else 
                  phylyr(i,j,k,posacc)=nf90_fill_double
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(finlyr)')
               stop '(finlyr)'
      endif
c     
      end subroutine finlyr



      subroutine wrth2d(pos,frmt,sfac,offs,cmpflg,msk,gridid,
     .  vnm,vlngnm,vstdnm,vunits)
c
c --- ------------------------------------------------------------------
c --- Description: writes diagnostic 2d field to file  
c ---   
c --- Arguments:
c ---   int  pos      (in)     : variable position in common buffer
c ---   int  frmt     (in)     : format/precision of output 
c ---                            0=field is not written  
c ---                            2=field is written as int2 with scale 
c ---                              factor and offset 
c ---                            4=field is written as real4
c ---                            8=field is written as real8
c ---   real sfac     (in)     : user defined scale factor to be applied   
c ---   real offs     (in)     : user defined offset to be added 
c ---   int  cmpflg   (in)     : compression flag; only wet points are 
c ---                            written if flag is set to 1 
c ---   int  msk      (in)       ocean mask 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   char vnm      (in)     : variable name used in nc-file 
c ---   char vlngnm   (in)     : variable long name (skipped if ' ') 
c ---   char vstdnm   (in)     : variable standard name (skipped if ' ') 
c ---   char vunits   (in)     : variable units (skipped if ' ') 
c --- ------------------------------------------------------------------
c
      implicit none
c 
      real :: sfac,offs
      integer :: frmt,cmpflg,pos
      character(len=*) :: gridid,vnm,vlngnm,vstdnm,vunits
      integer, dimension(*) :: msk
c
      character(len=100) :: dims
c
c --- Check whether field should be written
      if (frmt.eq.0) return
c
c --- Create dimension string 
      if (cmpflg.eq.1) then
        dims=gridid(1:1)//'comp time'
      else
        dims='x y time'
      endif
c
c --- Check output format
      if (frmt.eq.2) then
        if (cmpflg.eq.1) then
          call nccopa(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,sfac,
     .      offs)
        else
          call ncpack(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,
     .      1,sfac,offs)
        endif
      elseif (frmt.eq.4) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,sfac,
     .      offs,4)
        else
          call ncwrtr(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,
     .      1,sfac,offs,4)
        endif
      elseif (frmt.eq.8) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,sfac,
     .      offs,8)
        else
          call ncwrtr(vnm,dims,phyh2d(1-nbdy,1-nbdy,pos),msk,
     .      1,sfac,offs,8)
        endif
      else
        write (lp,*) 'unknown output format!'
        call xchalt('(wrth2d)')
               stop '(wrth2d)'
      endif
c
c --- Define attributes
!      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
!      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
!      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
!      call ncattr('coordinates',
!     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
!      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
c
      end subroutine wrth2d



      subroutine wrtlyr(pos,frmt,sfac,offs,cmpflg,msk,gridid,
     .  vnm,vlngnm,vstdnm,vunits)
c
c --- ------------------------------------------------------------------
c --- Description: writes diagnostic layer field to file  
c ---   
c --- Arguments:
c ---   int  pos      (in)     : variable position in common buffer
c ---   int  frmt     (in)     : format/precision of output 
c ---                            0=field is not written  
c ---                            2=field is written as int2 with scale 
c ---                              factor and offset 
c ---                            4=field is written as real4
c ---                            8=field is written as real8
c ---   real sfac     (in)     : user defined scale factor to be applied   
c ---   real offs     (in)     : user defined offset to be added 
c ---   int  cmpflg   (in)     : compression flag; only wet points are 
c ---                            written if flag is set to 1 
c ---   int  msk      (in)       ocean mask 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   char vnm      (in)     : variable name used in nc-file 
c ---   char vlngnm   (in)     : variable long name (skipped if ' ') 
c ---   char vstdnm   (in)     : variable standard name (skipped if ' ') 
c ---   char vunits   (in)     : variable units (skipped if ' ') 
c --- ------------------------------------------------------------------
c
      implicit none
c 
      real :: sfac,offs
      integer :: frmt,cmpflg,pos
      character(len=*) :: gridid,vnm,vlngnm,vstdnm,vunits
      integer, dimension(*) :: msk
c
      character(len=100) :: dims
c
c --- Check whether field should be written
      if (frmt.eq.0) return
c
c --- Create dimension string 
      if (cmpflg.eq.1) then
        dims=gridid(1:1)//'comp sigma time'
      else
        dims='x y sigma time'
      endif
c
c --- Check output format
      if (frmt.eq.2) then
        if (cmpflg.eq.1) then
          call nccopa(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,sfac,
     .      offs)
        else
          call ncpack(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,
     .      2,sfac,offs)
        endif
      elseif (frmt.eq.4) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,sfac,
     .      offs,4)
        else
          call ncwrtr(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,
     .      1,sfac,offs,4)
        endif
      elseif (frmt.eq.8) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,sfac,
     .      offs,8)
        else
          call ncwrtr(vnm,dims,phylyr(1-nbdy,1-nbdy,1,pos),msk,
     .      1,sfac,offs,8)
        endif
      else
        write (lp,*) 'unknown output format!'
        call xchalt('(wrtlyr)')
               stop '(wrtlyr)'
      endif
c
c --- Define attributes
!      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
!      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
!      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
!      call ncattr('coordinates',
!     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
!      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
c
      end subroutine wrtlyr



      subroutine wrtlvl(pos,frmt,sfac,offs,cmpflg,msk,gridid,
     .  vnm,vlngnm,vstdnm,vunits)
c
c --- ------------------------------------------------------------------
c --- Description: writes diagnostic level field to file  
c ---   
c --- Arguments:
c ---   int  pos      (in)     : variable position in common buffer
c ---   int  frmt     (in)     : format/precision of output 
c ---                            0=field is not written  
c ---                            2=field is written as int2 with scale 
c ---                              factor and offset 
c ---                            4=field is written as real4
c ---                            8=field is written as real8
c ---   real sfac     (in)     : user defined scale factor to be applied   
c ---   real offs     (in)     : user defined offset to be added 
c ---   int  cmpflg   (in)     : compression flag; only wet points are 
c ---                            written if flag is set to 1 
c ---   int  msk      (in)       ocean mask 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   char vnm      (in)     : variable name used in nc-file 
c ---   char vlngnm   (in)     : variable long name (skipped if ' ') 
c ---   char vstdnm   (in)     : variable standard name (skipped if ' ') 
c ---   char vunits   (in)     : variable units (skipped if ' ') 
c --- ------------------------------------------------------------------
c
      implicit none
c 
      real :: sfac,offs
      integer :: frmt,cmpflg,pos
      character(len=*) :: gridid,vnm,vlngnm,vstdnm,vunits
      integer, dimension(*) :: msk
c
      character(len=100) :: dims
c
c --- Check whether field should be written
      if (frmt.eq.0) return
c
c --- Create dimension string 
      if (cmpflg.eq.1) then
        dims=gridid//'comp depth time'
      else
        dims='x y depth time'
      endif
c
c --- Check output format
      if (frmt.eq.2) then
        if (cmpflg.eq.1) then
          call nccopa(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,sfac,
     .      offs)
        else
          call ncpack(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,
     .      2,sfac,offs)
        endif
      elseif (frmt.eq.4) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,sfac,
     .      offs,4)
        else
          call ncwrtr(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,
     .      1,sfac,offs,4)
        endif
      elseif (frmt.eq.8) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,sfac,
     .      offs,8)
        else
          call ncwrtr(vnm,dims,phylvl(1-nbdy,1-nbdy,1,pos),msk,
     .      1,sfac,offs,8)
        endif
      else
        write (lp,*) 'unknown output format!'
        call xchalt('(wrtlvl)')
               stop '(wrtlvl)'
      endif
c
c --- Define attributes
!      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
!      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
!      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
!      call ncattr('coordinates',
!     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
!      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
c
      end subroutine wrtlvl



      subroutine logh2d(pos,gridid,sfac,offs)
c
c --- ------------------------------------------------------------------
c --- Description: replace 2d field with log10(field) 
c ---   
c --- Arguments:
c ---   int  pos      (in)     : field position in layer buffer 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   real sfac     (in)     : scale factor to be applied before log10   
c ---   real offs     (in)     : offset to be added before log10   
c --- ------------------------------------------------------------------
c
      implicit none
c 
      real :: sfac,offs
      integer :: pos
      character :: gridid
c 
      integer :: i,j,l
      real :: epsil=1e-11
c
c --- Check whether field should be processed
      if (pos.eq.0) return
c
      if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (phyh2d(i,j,pos).lt.epsil) then
                phyh2d(i,j,pos)=0.
              else
                phyh2d(i,j,pos)=log10(phyh2d(i,j,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (phyh2d(i,j,pos).lt.epsil) then
                phyh2d(i,j,pos)=0.
              else
                phyh2d(i,j,pos)=log10(phyh2d(i,j,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (phyh2d(i,j,pos).lt.epsil) then
                phyh2d(i,j,pos)=0.
              else
                phyh2d(i,j,pos)=log10(phyh2d(i,j,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(logh2d)')
               stop '(logh2d)'
      endif
c
      end subroutine logh2d


      subroutine loglyr(pos,gridid,sfac,offs)
c
c --- ------------------------------------------------------------------
c --- Description: replace 3d layer field with log10(field) 
c ---   
c --- Arguments:
c ---   int  pos      (in)     : field position in layer buffer 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   real sfac     (in)     : scale factor to be applied before log10   
c ---   real offs     (in)     : offset to be added before log10   
c --- ------------------------------------------------------------------
c
      implicit none
c 
      real :: sfac,offs
      integer :: pos
      character :: gridid
c 
      integer :: i,j,k,l
      real :: epsil=1e-11
c
c --- Check whether field should be processed
      if (pos.eq.0) return
c
      if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                if (phylyr(i,j,k,pos).lt.epsil) then
                  phylyr(i,j,k,pos)=0.
                elseif (phylyr(i,j,k,pos).ne.nf90_fill_double) then
                  phylyr(i,j,k,pos)=log10(phylyr(i,j,k,pos)*sfac+offs)
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                if (phylyr(i,j,k,pos).lt.epsil) then
                  phylyr(i,j,k,pos)=0.
                elseif (phylyr(i,j,k,pos).ne.nf90_fill_double) then
                  phylyr(i,j,k,pos)=log10(phylyr(i,j,k,pos)*sfac+offs)
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (phylyr(i,j,k,pos).lt.epsil) then
                  phylyr(i,j,k,pos)=0.
                elseif (phylyr(i,j,k,pos).ne.nf90_fill_double) then
                  phylyr(i,j,k,pos)=log10(phylyr(i,j,k,pos)*sfac+offs)
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(loglyr)')
               stop '(loglyr)'
      endif
c
      end subroutine loglyr


      subroutine loglvl(pos,gridid,sfac,offs)
c
c --- ------------------------------------------------------------------
c --- Description: replace 3d level field with log10(field) 
c ---   
c --- Arguments:
c ---   int  pos      (in)     : field position in layer buffer 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c ---   real sfac     (in)     : scale factor to be applied before log10   
c ---   real offs     (in)     : offset to be added before log10   
c --- ------------------------------------------------------------------
c
      implicit none
c 
      real :: sfac,offs
      integer :: pos
      character :: gridid
c 
      integer :: i,j,k,l
      real :: epsil=1e-11
c
c --- Check whether field should be processed
      if (pos.eq.0) return
c
      if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,ddm
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                if (phylvl(i,j,k,pos).lt.epsil) then
                  phylvl(i,j,k,pos)=0.
                elseif (phylvl(i,j,k,pos).ne.nf90_fill_double) then
                  phylvl(i,j,k,pos)=log10(phylvl(i,j,k,pos)*sfac+offs)
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,ddm
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                if (phylvl(i,j,k,pos).lt.epsil) then
                  phylvl(i,j,k,pos)=0.
                elseif (phylvl(i,j,k,pos).ne.nf90_fill_double) then
                  phylvl(i,j,k,pos)=log10(phylvl(i,j,k,pos)*sfac+offs)
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,ddm
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (phylvl(i,j,k,pos).lt.epsil) then
                  phylvl(i,j,k,pos)=0.
                elseif (phylvl(i,j,k,pos).ne.nf90_fill_double) then
                  phylvl(i,j,k,pos)=log10(phylvl(i,j,k,pos)*sfac+offs)
                endif
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(loglvl)')
               stop '(loglvl)'
      endif
c
      end subroutine loglvl


      subroutine msklvl(pos,gridid)
c
c --- ------------------------------------------------------------------
c --- Description: set sea floor points to NaN in level fields 
c ---   
c --- Arguments:
c ---   int  pos      (in)     : field position in level buffer 
c ---   char gridid   (in)     : grid identifier ('p','u' or 'v')
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
c 
      integer :: pos
      character :: gridid
c 
      integer :: i,j,k
      logical :: iniflg=.true.
      integer, dimension(idm,jdm) :: kmaxu,kmaxv,kmaxp
      real, parameter :: mskval=nf90_fill_double
c
      save iniflg,kmaxu,kmaxv,kmaxp
c     
c --- Check whether field should be processed
      if (pos.eq.0) return
c
c --- Prepare index fields for masking

      if (iniflg) then
c$OMP PARALLEL DO PRIVATE(i,k)
        do j=1,jj
          do i=1,ii
            kmaxp(i,j)=0
            kmaxu(i,j)=0
            kmaxv(i,j)=0
          enddo
          do k=1,ddm
            do i=1,ii
              if (depths(i,j).gt.depthslev_bnds(1,k)) kmaxp(i,j)=k
              if (min(depths(i,j),depths(i-1,j)).gt.depthslev_bnds(1,k))
     .          kmaxu(i,j)=k
              if (min(depths(i,j),depths(i,j-1)).gt.depthslev_bnds(1,k))
     .          kmaxv(i,j)=k
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        iniflg=.false.
      endif
c
      if (gridid.eq.'u') then
c$OMP PARALLEL DO PRIVATE(i,k)
        do j=1,jj
          do i=1,ii
            do k=kmaxu(i,j)+1,ddm
              phylvl(i,j,k,pos)=mskval
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO PRIVATE(i,k)
        do j=1,jj
          do i=1,ii
            do k=kmaxv(i,j)+1,ddm
              phylvl(i,j,k,pos)=mskval
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(i,k)
        do j=1,jj
          do i=1,ii
            do k=kmaxp(i,j)+1,ddm
              phylvl(i,j,k,pos)=mskval
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(msklvl)')
               stop '(msklvl)'
      endif
c
      end subroutine msklvl

!define output variables subroutine definevar
      subroutine definevar(irec,iogrp,cmpflg,timeunits,calendar)
      integer irec,iogrp,cmpflg
      character timeunits*30,calendar*19
#if defined(TRC)
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt,nat,km
      character trcnm*80,trcnml*80
#endif
      call ncdefvar('time','time',ndouble,0)
      call ncattr('long_name','time')
      call ncattr('units',timeunits)
      call ncattr('calendar',calendar)
      if (irec.eq.1) then
c --- define sigma levels
       call ncdefvar('sigma','sigma',ndouble,8)
       call ncattr('long_name','Potential density')
       call ncattr('standard_name','sea_water_sigma_theta')
       call ncattr('units','kg m-3')
       call ncattr('positive','down')
c --- define zlevel
       call ncdefvar('depth','depth',ndouble,8)
       call ncattr('long_name','z level')
       call ncattr('units','m')
       call ncattr('positive','down')
       call ncattr('bounds','depth_bnds')
       call ncdefvar('depth_bnds','bounds depth',ndouble,8)
        if (MSC_MMFLXL(iogrp)+MSC_MMFLXD(iogrp)+MSC_MMFTDL(iogrp)
     .     +MSC_MMFTDD(iogrp)+MSC_MHFLX(iogrp)+MSC_MHFTD(iogrp)
     .     +MSC_MHFLD(iogrp)+MSC_MSFLX(iogrp)+MSC_MSFTD(iogrp)
     .     +MSC_MSFLD(iogrp)+MSC_VOLTR(iogrp).ne.0) then
        call ncdefvar('lat','lat',ndouble,8)
        call ncattr('long_name','Latitude')
        call ncattr('standard_name','latitude')
        call ncattr('units','degree_north')
        endif
        if (MSC_MMFLXL(iogrp)+MSC_MMFLXD(iogrp)+MSC_MMFTDL(iogrp)
     .     +MSC_MMFTDD(iogrp).ne.0) then
        call ncdefvar('region','slenmax region',nchar,0)
        call ncattr('long_name','Region name')
        endif
        if (MSC_MHFLX(iogrp)+MSC_MHFTD(iogrp)+MSC_MHFLD(iogrp)
     .     +MSC_MSFLX(iogrp)+MSC_MSFTD(iogrp)+MSC_MSFLD(iogrp).ne.0)
     .  then
        call ncdefvar('section','slenmax section',nchar,0)
        call ncattr('long_name','Section name')
        endif
       endif
c --- define 2d fields
      call ncdefvar3d(H2D_SIGMX(iogrp),cmpflg,'p','sigmx',
     . 'Mixed layer density',' ','kg m-3',0)
c
      call ncdefvar3d(H2D_UB(iogrp),cmpflg,'u','ubaro',
     . 'Barotropic velocity x-component',' ','m s-1',0)
c
      call ncdefvar3d(H2D_VB(iogrp),cmpflg,'v','vbaro',
     . 'Barotropic velocity y-component',' ','m s-1',0)
c
      call ncdefvar3d(H2D_PSRF(iogrp),cmpflg,'p','psrf',
     . 'Surface pressure',' ','Pa',0)
c
      call ncdefvar3d(H2D_PBOT(iogrp),cmpflg,'p','pbot',
     . 'Bottom pressure',' ','Pa',0)
c
      call ncdefvar3d(H2D_SEALV(iogrp),cmpflg,'p','sealv',
     . 'Sea level',' ','m',0)
c
      call ncdefvar3d(H2D_SLVSQ(iogrp),cmpflg,'p','slvsq',
     . 'Sea level squared',' ','m2',0)
c
      call ncdefvar3d(H2D_BTMSTR(iogrp),cmpflg,'p','btmstr',
     . 'Barotropic mass streamfunction',' ','kg s-1',0)
c
      call ncdefvar3d(H2D_HICE(iogrp),cmpflg,'p','hice',
     . 'Ice thickness',' ','m',0)
c
      call ncdefvar3d(H2D_TICE(iogrp),cmpflg,'p','tice',
     . 'Ice temperature',' ','degC',0)
c
      call ncdefvar3d(H2D_HSNW(iogrp),cmpflg,'p','hsnw',
     . 'Snow depth',' ','m',0)
c
      call ncdefvar3d(H2D_FICE(iogrp),cmpflg,'p','fice',
     . 'Ice concentration',' ','%',0)
c
      call ncdefvar3d(H2D_TSRF(iogrp),cmpflg,'p','tsrf',
     . 'Surface temperature',' ','degC',0)
c
      call ncdefvar3d(H2D_IAGE(iogrp),cmpflg,'p','iage',
     . 'Ice age',' ','day',0)
c
      call ncdefvar3d(H2D_UICE(iogrp),cmpflg,'u','uice',
     . 'Ice velocity x-component',' ','m s-1',0)
c
      call ncdefvar3d(H2D_VICE(iogrp),cmpflg,'v','vice',
     . 'Ice velocity y-component',' ','m s-1',0)
c
      call ncdefvar3d(H2D_SWA(iogrp),cmpflg,'p','swa',
     . 'Short-wave heat flux',' ','W m-2',0)
c
      call ncdefvar3d(H2D_NSF(iogrp),cmpflg,'p','nsf',
     . 'Non-solar heat flux',' ','W m-2',0)
c
      call ncdefvar3d(H2D_HMLTFZ(iogrp),cmpflg,'p','hmltfz',
     .  'Heat flux due to melting/freezing',' ','W m-2',0)
c
      call ncdefvar3d(H2D_DFL(iogrp),cmpflg,'p','dfl',
     . 'Non-solar heat flux derivative',' ','W m-2 K-1',0)
c
      call ncdefvar3d(H2D_SURFLX(iogrp),cmpflg,'p','hflx',
     .  'Heat flux received by ocean',' ','W m-2',0)
c
      call ncdefvar3d(H2D_SURRLX(iogrp),cmpflg,'p','hrflx',
     .  'Restoring heat flux received by ocean',' ','W m-2',0)
c
      call ncdefvar3d(H2D_LIP(iogrp),cmpflg,'p','lip',
     . 'Liquid precipitation',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_SOP(iogrp),cmpflg,'p','sop',
     . 'Solid precipitation',' ','kg m-2 s-1',0)

      call ncdefvar3d(H2D_EVA(iogrp),cmpflg,'p','eva',
     . 'Evaporation',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_FMLTFZ(iogrp),cmpflg,'p','fmltfz',
     .  'Fresh water flux due to melting/freezing',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_RNFFLX(iogrp),cmpflg,'p','rnf',
     . 'Liquid runoff',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_RFIFLX(iogrp),cmpflg,'p','rfi',
     . 'Frozen runoff',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_SALFLX(iogrp),cmpflg,'p','sflx',
     .  'Salt flux received by ocean',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_SALRLX(iogrp),cmpflg,'p','srflx',
     .  'Restoring salt flux received by ocean',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_BRNFLX(iogrp),cmpflg,'p','bflx',
     . 'Brine flux',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_ZTX(iogrp),cmpflg,'u','ztx',
     . 'Wind stress x-component',' ','N m-2',0)
c
      call ncdefvar3d(H2D_MTY(iogrp),cmpflg,'v','mty',
     . 'Wind stress y-component',' ','N m-2',0)
c
      call ncdefvar3d(H2D_TAUX(iogrp),cmpflg,'u','taux',
     .  'Momentum flux received by ocean x-component',' ','N m-2',0)
c
      call ncdefvar3d(H2D_TAUY(iogrp),cmpflg,'v','tauy',
     .  'Momentum flux received by ocean y-component',' ','N m-2',0)
c
      call ncdefvar3d(H2D_IDKEDT(iogrp),cmpflg,'p','idkedt',
     .  'Mixed layer inertial kinetic energy tendency per unit area',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_USTAR(iogrp),cmpflg,'p','ustar',
     . 'Friction velocity',' ','m s-1',0)
c
      call ncdefvar3d(H2D_USTAR3(iogrp),cmpflg,'p','ustar3',
     .  'Friction velocity cubed',' ','m3 s-3',0)
c
      call ncdefvar3d(H2D_ABSWND(iogrp),cmpflg,'p','abswnd',
     . 'Absolute wind speed',' ','m s-1',0)
c
      call ncdefvar3d(H2D_MTKEUS(iogrp),cmpflg,'p','mtkeus',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to friction velocity',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_MTKENI(iogrp),cmpflg,'p','mtkeni',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to near inertial motions',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_MTKEBF(iogrp),cmpflg,'p','mtkebf',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to buoyancy forcing',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_MTKERS(iogrp),cmpflg,'p','mtkers',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to eddy restratification',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_MTKEPE(iogrp),cmpflg,'p','mtkepe',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to potential energy change',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_MTKEKE(iogrp),cmpflg,'p','mtkeke',
     .  'Mixed layer turbulent kinetic energy tendency '//
     .  'per unit area related to kinetic energy change',
     .  ' ','kg s-3',0)
c
      call ncdefvar3d(H2D_SFL(iogrp),cmpflg,'p','sfl',
     . 'Salt flux',' ','kg m-2 s-1',0)
c
      call ncdefvar3d(H2D_ALB(iogrp),cmpflg,'p','alb',
     . 'Surface albedo',' ','1',0)
c
      call ncdefvar3d(H2D_MLD(iogrp),cmpflg,'p','mld',
     . 'Mixed layer depth',' ','m',0)
c
      call ncdefvar3d(H2D_MLDU(iogrp),cmpflg,'u','mldu',
     .  'Mixed layer depth at u-point',' ','m',0)
c
      call ncdefvar3d(H2D_MLDV(iogrp),cmpflg,'v','mldv',
     .  'Mixed layer depth at v-point',' ','m',0)
c
      call ncdefvar3d(H2D_MAXMLD(iogrp),cmpflg,'p','maxmld',
     . 'Maximum mixed layer depth',' ','m',0)
c
      call ncdefvar3d(H2D_MLTS(iogrp),cmpflg,'p','mlts',
     . 'Mixed layer thickness defined by sigma t',' ','m',0)
c
      call ncdefvar3d(H2D_MLTSMN(iogrp),cmpflg,'p','mltsmn',
     . 'Minimum mixed layer thickness defined by sigma t',' ','m',0)
c
      call ncdefvar3d(H2D_MLTSMX(iogrp),cmpflg,'p','mltsmx',
     . 'Maximum mixed layer thickness defined by sigma t',' ','m',0)
c
      call ncdefvar3d(H2D_MLTSSQ(iogrp),cmpflg,'p','mltssq',
     . 'Maximum mixed layer thickness squared defined by sigma t',' ',
     . 'm2',0)
c
      call ncdefvar3d(H2D_T20D(iogrp),cmpflg,'p','t20d',
     . '20C isoterm depth',' ','m',0)
c
      call ncdefvar3d(H2D_BRNPD(iogrp),cmpflg,'p','brnpd',
     . 'Brine plume depth',' ','m',0)
c
      call ncdefvar3d(H2D_SSS(iogrp),cmpflg,'p','sss',
     . 'Ocean surface salinity',' ','g kg-1',0)
c
      call ncdefvar3d(H2D_SSSSQ(iogrp),cmpflg,'p','ssssq',
     . 'Ocean surface salinity squared',' ','g2 kg-2',0)
c
      call ncdefvar3d(H2D_SBOT(iogrp),cmpflg,'p','sbot',
     . 'Bottom salinity',' ','g kg-1',0)
c
      call ncdefvar3d(H2D_SST(iogrp),cmpflg,'p','sst',
     . 'Ocean surface temperature',' ','degC',0)
c
      call ncdefvar3d(H2D_SSTSQ(iogrp),cmpflg,'p','sstsq',
     . 'Ocean surface temperature squared',' ','degC2',0)
c
      call ncdefvar3d(H2D_TBOT(iogrp),cmpflg,'p','tbot',
     . 'Bottom temperature',' ','degC',0)
c
      call ncdefvar3d(H2D_MXLU(iogrp),cmpflg,'u','mxlu',
     . 'Mixed layer velocity x-component',' ','m s-1',0)
c
      call ncdefvar3d(H2D_MXLV(iogrp),cmpflg,'v','mxlv',
     . 'Mixed layer velocity y-component',' ','m s-1',0)
c
c --- define 3d layer fields
      call ncdefvar3d(LYR_DP(iogrp),cmpflg,'p','dp',
     . 'Layer pressure thickness',' ','Pa',1)
c
      call ncdefvar3d(LYR_DZ(iogrp),cmpflg,'p','dz',
     . 'Layer thickness',' ','m',1)
c
      call ncdefvar3d(LYR_TEMP(iogrp),cmpflg,'p','temp',
     . 'Temperature','Ocean temperature','degC',1)
c
      call ncdefvar3d(LYR_SALN(iogrp),cmpflg,'p','saln',
     . 'Salinity','Ocean salinity','g kg-1',1)
c
      call ncdefvar3d(LYR_UVEL(iogrp),cmpflg,'u','uvel',
     . 'Velocity x-component',' ','m s-1',1)
c
      call ncdefvar3d(LYR_VVEL(iogrp),cmpflg,'v','vvel',
     . 'Velocity y-component',' ','m s-1',1)
c
      call ncdefvar3d(LYR_UFLX(iogrp),cmpflg,'u','uflx',
     .  'Mass flux in x-direction',' ','kg s-1',1)
c
      call ncdefvar3d(LYR_VFLX(iogrp),cmpflg,'v','vflx',
     .  'Mass flux in y-direction',' ','kg s-1',1)
c
      call ncdefvar3d(LYR_UTFLX(iogrp),cmpflg,'u','uhflx',
     .  'Heat flux in x-direction',' ','W',1)
c
      call ncdefvar3d(LYR_VTFLX(iogrp),cmpflg,'v','vhflx',
     .  'Heat flux in y-direction',' ','W',1)
c
      call ncdefvar3d(LYR_USFLX(iogrp),cmpflg,'u','usflx',
     .  'Salt flux in x-direction',' ','kg s-1',1)
c
      call ncdefvar3d(LYR_VSFLX(iogrp),cmpflg,'v','vsflx',
     .  'Salt flux in y-direction',' ','kg s-1',1)
c
      call ncdefvar3d(LYR_UMFLTD(iogrp),cmpflg,'u','umfltd',
     .  'Mass flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1',1)
c
      call ncdefvar3d(LYR_VMFLTD(iogrp),cmpflg,'v','vmfltd',
     .  'Mass flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1',1)
c
      call ncdefvar3d(LYR_UTFLTD(iogrp),cmpflg,'u','uhfltd',
     .  'Heat flux due to thickness diffusion in x-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_VTFLTD(iogrp),cmpflg,'v','vhfltd',
     .  'Heat flux due to thickness diffusion in y-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_UTFLLD(iogrp),cmpflg,'u','uhflld',
     .  'Heat flux due to lateral diffusion in x-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_VTFLLD(iogrp),cmpflg,'v','vhflld',
     .  'Heat flux due to lateral diffusion in y-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_USFLTD(iogrp),cmpflg,'u','usfltd',
     .  'Salt flux due to thickness diffusion in x-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_VSFLTD(iogrp),cmpflg,'v','vsfltd',
     .  'Salt flux due to thickness diffusion in y-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_USFLLD(iogrp),cmpflg,'u','usflld',
     .  'Salt flux due to lateral diffusion in x-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_VSFLLD(iogrp),cmpflg,'v','vsflld',
     .  'Salt flux due to lateral diffusion in y-direction',' ',
     .  'W',1)
c
      call ncdefvar3d(LYR_WFLX(iogrp),cmpflg,'p','wflx',
     .  'Vertical mass flux',' ','kg s-1',1)
c
      call ncdefvar3d(LYR_WFLX2(iogrp),cmpflg,'p','wflx2',
     .  'Vertical mass flux squared',' ','kg2 s-2',1)
c
      call ncdefvar3d(LYR_BFSQ(iogrp),cmpflg,'p','bfsq',
     .  'Squared buoyancy frequency',' ','s-1',1)
c
      call ncdefvar3d(LYR_PV(iogrp),cmpflg,'p','pv',
     .  'Potential vorticity',' ','m-1 s-1',1)
c
      if (LYR_DIFINT(iogrp).eq.2) then
        call ncdefvar3d(LYR_DIFINT(iogrp),cmpflg,'p','difint',
     . 'Layer interface diffusivity',' ','log10(m2 s-1)',1)
      else
        call ncdefvar3d(LYR_DIFINT(iogrp),cmpflg,'p','difint',
     . 'Layer interface diffusivity',' ','m2 s-1',1)
      endif
c
      if (LYR_DIFISO(iogrp).eq.2) then
        call ncdefvar3d(LYR_DIFISO(iogrp),cmpflg,'p','difiso',
     . 'Isopycnal diffusivity',' ','log10(m2 s-1)',1)
      else
        call ncdefvar3d(LYR_DIFISO(iogrp),cmpflg,'p','difiso',
     . 'Isopycnal diffusivity',' ','m2 s-1',1)
      endif
c
      if (LYR_DIFDIA(iogrp).eq.2) then
        call ncdefvar3d(LYR_DIFDIA(iogrp),cmpflg,'p','difdia',
     . 'Diapycnal diffusivity',' ','log10(m2 s-1)',1)
      else
        call ncdefvar3d(LYR_DIFDIA(iogrp),cmpflg,'p','difdia',
     . 'Diapycnal diffusivity',' ','m2 s-1',1)
      endif
c
#if defined TKE 

      call ncdefvar3d(LYR_TKE(iogrp),cmpflg,'p','tke',
     . 'TKE','Turbulent kinetic energy','m2 s-2',1)
c
      call ncdefvar3d(LYR_GLS_PSI(iogrp),cmpflg,'p','gls_psi',
     . 'GLS_PSI','Generic length scale','m2 s-3',1)
c
#endif
c --- define 3d depth fields
      call ncdefvar3d(LVL_DZ(iogrp),cmpflg,'p','dzlvl',
     . 'Layer thickness',' ','m',2)
c
      call ncdefvar3d(LVL_TEMP(iogrp),cmpflg,'p','templvl',
     . 'Temperature','Ocean temperature','degC',2)
c
      call ncdefvar3d(LVL_SALN(iogrp),cmpflg,'p','salnlvl',
     . 'Salinity','Ocean salinity','g kg-1',2)
c
      call ncdefvar3d(LVL_UVEL(iogrp),cmpflg,'u','uvellvl',
     .  'Velocity x-component',' ','m s-1',2)
c
      call ncdefvar3d(LVL_VVEL(iogrp),cmpflg,'v','vvellvl',
     .  'Velocity y-component',' ','m s-1',2)
c 
      call ncdefvar3d(LVL_UFLX(iogrp),cmpflg,'u','uflxlvl',
     .  'Mass flux in x-direction',' ','kg s-1',2)
c
      call ncdefvar3d(LVL_VFLX(iogrp),cmpflg,'v','vflxlvl',
     .  'Mass flux in y-direction',' ','kg s-1',2)
c
      call ncdefvar3d(LVL_UTFLX(iogrp),cmpflg,'u','uhflxlvl',
     .  'Heat flux in x-direction',' ','W',2)
c
      call ncdefvar3d(LVL_VTFLX(iogrp),cmpflg,'v','vhflxlvl',
     .  'Heat flux in y-direction',' ','W',2)
c
      call ncdefvar3d(LVL_USFLX(iogrp),cmpflg,'u','usflxlvl',
     .  'Salt flux in x-direction',' ','kg s-1',2)
c
      call ncdefvar3d(LVL_VSFLX(iogrp),cmpflg,'v','vsflxlvl',
     .  'Salt flux in y-direction',' ','kg s-1',2)
c
      call ncdefvar3d(LVL_UMFLTD(iogrp),cmpflg,'u','umfltdlvl',
     .  'Mass flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1',2)
c
      call ncdefvar3d(LVL_VMFLTD(iogrp),cmpflg,'v','vmfltdlvl',
     .  'Mass flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1',2)
c
      call ncdefvar3d(LVL_UTFLTD(iogrp),cmpflg,'u','uhfltdlvl',
     .  'Heat flux due to thickness diffusion in x-direction',' ',
     .  'W',2)
c
      call ncdefvar3d(LVL_VTFLTD(iogrp),cmpflg,'v','vhfltdlvl',
     .  'Heat flux due to thickness diffusion in y-direction',' ',
     .  'W',2)
c
      call ncdefvar3d(LVL_UTFLLD(iogrp),cmpflg,'u','uhflldlvl',
     .  'Heat flux due to lateral diffusion in x-direction',' ',
     .  'W',2)
c
      call ncdefvar3d(LVL_VTFLLD(iogrp),cmpflg,'v','vhflldlvl',
     .  'Heat flux due to lateral diffusion in y-direction',' ',
     .  'W',2)
c
      call ncdefvar3d(LVL_USFLTD(iogrp),cmpflg,'u','usfltdlvl',
     .  'Salt flux due to thickness diffusion in x-direction',' ',
     .  'kg s-1',2)
c
      call ncdefvar3d(LVL_VSFLTD(iogrp),cmpflg,'v','vsfltdlvl',
     .  'Salt flux due to thickness diffusion in y-direction',' ',
     .  'kg s-1',2)
c
      call ncdefvar3d(LVL_USFLLD(iogrp),cmpflg,'u','usflldlvl',
     .  'Salt flux due to lateral diffusion in x-direction',' ',
     .  'kg s-1',2)
c
      call ncdefvar3d(LVL_VSFLLD(iogrp),cmpflg,'v','vsflldlvl',
     .  'Salt flux due to lateral diffusion in y-direction',' ',
     .  'kg s-1',2)
c
      call ncdefvar3d(LVL_WFLX(iogrp),cmpflg,'p','wflxlvl',
     .  'Vertical mass flux',' ','kg s-1',2)
c
      call ncdefvar3d(LVL_WFLX2(iogrp),cmpflg,'p','wflx2lvl',
     .  'Vertical mass flux squared',' ','kg2 s-2',2)
c
      call ncdefvar3d(LVL_BFSQ(iogrp),cmpflg,'p','bfsqlvl',
     .  'Squared buoyancy frequancy',' ','s-1',2)
c
      call ncdefvar3d(LVL_PV(iogrp),cmpflg,'p','pvlvl',
     .  'Potential vorticity',' ','m-1 s-1',2)
c
      if (LVL_DIFINT(iogrp).eq.2) then
        call ncdefvar3d(LVL_DIFINT(iogrp),cmpflg,'p','difintlvl',
     . 'Layer interface diffusivity',' ','log10(m2 s-1)',2)
      else
        call ncdefvar3d(LVL_DIFINT(iogrp),cmpflg,'p','difintlvl',
     . 'Layer interface diffusivity',' ','m2 s-1',2)
      endif
c
      if (LVL_DIFISO(iogrp).eq.2) then
        call ncdefvar3d(LVL_DIFISO(iogrp),cmpflg,'p','difisolvl',
     . 'Isopycnal diffusivity',' ','log10(m2 s-1)',2)
      else
        call ncdefvar3d(LVL_DIFISO(iogrp),cmpflg,'p','difisolvl',
     . 'Isopycnal diffusivity',' ','m2 s-1',2)
      endif
c
      if (LVL_DIFDIA(iogrp).eq.2) then
        call ncdefvar3d(LVL_DIFDIA(iogrp),cmpflg,'p','difdialvl',
     . 'Diapycnal diffusivity',' ','log10(m2 s-1)',2)
      else
        call ncdefvar3d(LVL_DIFDIA(iogrp),cmpflg,'p','difdialvl',
     . 'Diapycnal diffusivity',' ','m2 s-1',2)
      endif
c
#if defined(TRC) && defined(TKE)
      call ncdefvar3d(LVL_TKE(iogrp),cmpflg,'p','tkelvl',
     . 'Turbulent Kinetic Energy',' ','m2 s-2',2)
c
      call ncdefvar3d(LVL_GLS_PSI(iogrp),cmpflg,'p','gls_psilvl',
     . 'Generic length scale',' ','m2 s-3',2)
c
#endif
c
c --- define meridional transports 
      if (MSC_MMFLXL(iogrp).ne.0) then
        call ncdefvar('mmflxl','lat sigma region time',ndouble,8)
        call ncattr('long_name',
     .    'Overturning stream-function on isopycnic layers')
        call ncattr('units','kg s-1')
      endif
      if (MSC_MMFLXD(iogrp).ne.0) then
        call ncdefvar('mmflxd','lat depth region time',ndouble,8)
        call ncattr('long_name',
     .    'Overturning stream-function on z-levels')
        call ncattr('units','kg s-1')
      endif
      if (MSC_MMFTDL(iogrp).ne.0) then
        call ncdefvar('mmftdl','lat sigma region time',ndouble,8)
        call ncattr('long_name',
     .    'Overturning stream-function due to thickness diffusion '//
     .    'on isopycnic layers')
        call ncattr('units','kg s-1')
      endif
      if (MSC_MMFTDD(iogrp).ne.0) then
        call ncdefvar('mmftdd','lat depth region time',ndouble,8)
        call ncattr('long_name',
     .    'Overturning stream-function due to thickness diffusion '//
     .    'on z-levels')
        call ncattr('units','kg s-1')
      endif
      if (MSC_MHFLX(iogrp).ne.0) then
        call ncdefvar('mhflx','lat region time',ndouble,8)
        call ncattr('long_name','Meridional heat flux')
        call ncattr('units','W')
      endif
      if (MSC_MHFTD(iogrp).ne.0) then
        call ncdefvar('mhftd','lat region time',ndouble,8)
        call ncattr('long_name',
     .    'Meridional heat flux due to thickness diffusion')
        call ncattr('units','W')
      endif
      if (MSC_MHFLD(iogrp).ne.0) then
        call ncdefvar('mhfld','lat region time',ndouble,8)
        call ncattr('long_name',
     .    'Meridional heat flux due to lateral diffusion')
        call ncattr('units','W')
      endif
      if (MSC_MSFLX(iogrp).ne.0) then
        call ncdefvar('msflx','lat region time',ndouble,8)
        call ncattr('long_name','Meridional salt flux')
        call ncattr('units','kg s-1')
      endif
      if (MSC_MSFTD(iogrp).ne.0) then
        call ncdefvar('msftd','lat region time',ndouble,8)
        call ncattr('long_name',
     .    'Meridional salt flux due to thickness diffusion')
        call ncattr('units','kg s-1')
      endif
      if (MSC_MSFLD(iogrp).ne.0) then
        call ncdefvar('msfld','lat region time',ndouble,8)
        call ncattr('long_name',
     .    'Meridional salt flux due to lateral diffusion')
        call ncattr('units','kg s-1')
      endif
c
c --- store section transports 
      if (MSC_VOLTR(iogrp).ne.0) then
        call ncdefvar('voltr','section time',ndouble,8)
        call ncattr('long_name','Section transports')
        call ncattr('units','kg s-1')
      endif
c
c --- store global sums and averages
      if (MSC_MASSGS(iogrp).ne.0) then
        call ncdefvar('massgs','time',ndouble,8)
        call ncattr('long_name','Sea water mass')
        call ncattr('units','kg')
      endif
      if (MSC_VOLGS(iogrp).ne.0) then
        call ncdefvar('volgs','time',ndouble,8)
        call ncattr('long_name','Sea water volume')
        call ncattr('units','m3')
      endif
      if (MSC_SALNGA(iogrp).ne.0) then
        call ncdefvar('salnga','time',ndouble,8)
        call ncattr('long_name','Global average salinity')
        call ncattr('units','g kg-1')
      endif
      if (MSC_TEMPGA(iogrp).ne.0) then
        call ncdefvar('tempga','time',ndouble,8)
        call ncattr('long_name','Global average temperature')
        call ncattr('units','degC')
      endif
      if (MSC_SSSGA(iogrp).ne.0) then
        call ncdefvar('sssga','time',ndouble,8)
        call ncattr('long_name','Global average sea surface salinity')
        call ncattr('units','g kg-1')
      endif
      if (MSC_SSTGA(iogrp).ne.0) then
        call ncdefvar('sstga','time',ndouble,8)
        call ncattr('long_name',
     .              'Global average sea surface temperature')
        call ncattr('units','degC')
      endif
c
#ifdef TRC
      if (LYR_IDLAGE(iogrp).ne.0.or.LYR_TRC(iogrp).ne.0) then
        call ncdefvar3d(max(LYR_IDLAGE(iogrp),LYR_TRC(iogrp)),cmpflg,
     . 'p','dp_trc','Layer pressure thickness',' ','Pa',1)
      endif
c --- ideal age tracer 
#if IDLAGE
      call ncdefvar3d(LYR_IDLAGE(iogrp),cmpflg,'p','idlage',
     . 'Ideal age','sea_water_age_since_surface_contact','year',1)

      if (LVL_IDLAGE(iogrp).ne.0) then
        call ncdefvar3d(LVL_IDLAGE(iogrp),cmpflg,'p','idlagelvl',
     . 'Ideal age','sea_water_age_since_surface_contact','year',2)
      endif
#endif
c
c --- ocean tracers
      if (LYR_TRC(iogrp).gt.0.and.ntrocn.gt.0) then
#  ifdef ATRC
      do nt=1,ntrocn-natr
      write (trcnm,'(a,i3.3)') 'trc',nt
      write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
      call ncdefvar3d(LYR_TRC(iogrp),cmpflg,'p',trim(trcnm),
     .     trim(trcnml),' ',' ',1)
      enddo
      do nt=1,natr
          nat=ntr-natr+nt
      write (trcnm,'(a,i3.3)') 'atrc',nt
          write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
          call ncdefvar3d(LYR_TRC(iogrp),cmpflg,'p',trim(trcnm),
     .     trim(trcnml),' ',' ',1)
        enddo
#  else
      do nt=1,ntrocn
      write (trcnm,'(a,i3.3)') 'trc',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call ncdefvar3d(LYR_TRC(iogrp),cmpflg,'p',trim(trcnm),
     .         trim(trcnml),' ',' ',1)
        enddo
#  endif
      endif
      if (LVL_TRC(iogrp).gt.0.and.ntrocn.gt.0) then
#  ifdef ATRC
        do nt=1,ntrocn-natr
          write (trcnm,'(a,i3.3)') 'trclvl',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call ncdefvar3d(LVL_TRC(iogrp),cmpflg,'p',trim(trcnm),
     .         trim(trcnml),' ',' ',2)
        enddo
        do nt=1,natr
          nat=ntr-natr+nt
          write (trcnm,'(a,i3.3)') 'atrclvl',nt
          write (trcnml,'(a,i3.3)') 'Ocean age tracer ',nt
          call ncdefvar3d(LVL_TRC(iogrp),cmpflg,'p',trim(trcnm),
     .      trim(trcnml),' ',' ',2)
        enddo
#  else
        do nt=1,ntrocn
          write (trcnm,'(a,i3.3)') 'trclvl',nt
          write (trcnml,'(a,i3.3)') 'Ocean tracer ',nt
          call ncdefvar3d(LVL_TRC(iogrp),cmpflg,'p',trim(trcnm),
     .        trim(trcnml),' ',' ',2)
        enddo
#  endif
      endif
#endif
      call ncedef
      end subroutine definevar
      end module mod_dia
