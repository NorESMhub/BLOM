! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen, Mehmet Ilicak, Ingo Bethke
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

      subroutine rdlim
c
c --- ------------------------------------------------------------------
c --- Read limits file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_swabs, only: swamth, jwtype, chlopt, ccfile
      use mod_dia
      use mod_ben02, only: atm_path, atm_path_len
      use data_mct, only: runid_mct, ocn_cpl_dt_mct
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
c
      integer :: m,n,nd0rst,nm0rst,ny0rst
      character*256 :: runtyp,rstfnm
      logical :: fexist
      integer :: idate,idate0,nstep_in_cpl
      namelist /limits/ nday1,nday2,idate,idate0,runid,expcnf,runtyp,
     .  grfile,icfile,baclin,batrop,mdv2hi,mdv2lo,mdv4hi,mdv4lo,mdc2hi,
     .  mdc2lo,vsc2hi,vsc2lo,vsc4hi,vsc4lo,cbar,cb,cwbdts,cwbdls,mommth,
     .  eitmth,edritp,bmcmth,rmpmth,edwmth,mlrttp,edsprs,egc,eggam,
     .  eglsmn,egmndf,egmxdf,egidfq,ri0,rm0,rm5,ce,bdmtyp,bdmc1,bdmc2,
     .  tdfile,tkepf,niwgf,niwbf,niwlf,swamth,jwtype,chlopt,ccfile,
     .  trxday,srxday,trxdpt,srxdpt,trxlim,srxlim,aptflx,apsflx,ditflx,
     .  disflx,srxbal,scfile,smtfrc,sprfac,atm_path,rstfrq,rstfmt,
     .  rstcmp,iotype
c
c --- read namelists
      if (mnproc.eq.1) then
c
        GLB_AVEPERIO=-999
        inquire(file='ocn_in',exist=fexist)
        if (fexist) then
          open (unit=nfu,file='ocn_in',status='old',action='read',
     .          recl=80)
        else
          inquire(file='limits',exist=fexist)
          if (fexist) then
            open (unit=nfu,file='limits',status='old',action='read',
     .            recl=80)
          else
            write (lp,*) 'rdlim: could not find namelist file!'
            call xchalt('(rdlim)')
                   stop '(rdlim)'
          endif
        endif
        read (unit=nfu,nml=LIMITS)
        read (unit=nfu,nml=MERDIA)
        read (unit=nfu,nml=SECDIA)
        read (unit=nfu,nml=DIAPHY)
        close (unit=nfu)
c
c --- - print namelist on stdout 
        write (lp,*)
        write (lp,*) 'rdlim: BLOM LIMITS NAMELIST GROUP:' 
        write (lp,*) 'NDAY1',NDAY1
        write (lp,*) 'NDAY2',NDAY2
        write (lp,*) 'IDATE',IDATE
        write (lp,*) 'IDATE0',IDATE0
        write (lp,*) 'RUNID ',trim(RUNID)
        write (lp,*) 'EXPCNF ',trim(EXPCNF)
        write (lp,*) 'RUNTYP ',trim(RUNTYP)
        write (lp,*) 'GRFILE ',trim(GRFILE)
        write (lp,*) 'ICFILE ',trim(ICFILE)
        write (lp,*) 'BACLIN',BACLIN
        write (lp,*) 'BATROP',BATROP
        write (lp,*) 'MDV2HI',MDV2HI
        write (lp,*) 'MDV2LO',MDV2LO
        write (lp,*) 'MDV4HI',MDV4HI
        write (lp,*) 'MDV4LO',MDV4LO
        write (lp,*) 'MDC2HI',MDC2HI
        write (lp,*) 'MDC2LO',MDC2LO
        write (lp,*) 'VSC2HI',VSC2HI
        write (lp,*) 'VSC2LO',VSC2LO
        write (lp,*) 'VSC4HI',VSC4HI
        write (lp,*) 'VSC4LO',VSC4LO
        write (lp,*) 'CBAR',CBAR
        write (lp,*) 'CB',CB
        write (lp,*) 'CWBDTS',CWBDTS
        write (lp,*) 'CWBDLS',CWBDLS
        write (lp,*) 'MOMMTH ',trim(MOMMTH)
        write (lp,*) 'EITMTH ',trim(EITMTH)
        write (lp,*) 'EDRITP ',trim(EDRITP)
        write (lp,*) 'BMCMTH ',trim(BMCMTH)
        write (lp,*) 'RMPMTH ',trim(RMPMTH)
        write (lp,*) 'EDWMTH ',trim(EDWMTH)
        write (lp,*) 'EDSPRS ',EDSPRS
        write (lp,*) 'EGC',EGC
        write (lp,*) 'EGGAM',EGGAM
        write (lp,*) 'EGLSMN',EGLSMN
        write (lp,*) 'EGMNDF',EGMNDF
        write (lp,*) 'EGMXDF',EGMXDF
        write (lp,*) 'EGIDFQ',EGIDFQ
        write (lp,*) 'RI0',RI0
        write (lp,*) 'RM0',RM0
        write (lp,*) 'RM5',RM5
        write (lp,*) 'CE',CE
        write (lp,*) 'BDMTYP',BDMTYP
        write (lp,*) 'BDMC1',BDMC1
        write (lp,*) 'BDMC2',BDMC2
        write (lp,*) 'TDFILE',trim(TDFILE)
        write (lp,*) 'TKEPF',TKEPF
        write (lp,*) 'NIWGF',NIWGF
        write (lp,*) 'NIWBF',NIWBF
        write (lp,*) 'NIWLF',NIWLF
        write (lp,*) 'SWAMTH ',trim(SWAMTH)
        write (lp,*) 'JWTYPE',JWTYPE
        write (lp,*) 'CHLOPT ',trim(CHLOPT)
        write (lp,*) 'CCFILE ',trim(CCFILE)
        write (lp,*) 'TRXDAY',TRXDAY
        write (lp,*) 'SRXDAY',SRXDAY
        write (lp,*) 'TRXDPT',TRXDPT
        write (lp,*) 'SRXDPT',SRXDPT
        write (lp,*) 'TRXLIM',TRXLIM
        write (lp,*) 'SRXLIM',SRXLIM
        write (lp,*) 'APTFLX',APTFLX
        write (lp,*) 'APSFLX',APSFLX
        write (lp,*) 'DITFLX',DITFLX
        write (lp,*) 'DISFLX',DISFLX
        write (lp,*) 'SRXBAL',SRXBAL
        write (lp,*) 'SCFILE ',trim(SCFILE)
        write (lp,*) 'SMTFRC',SMTFRC
        write (lp,*) 'SPRFAC',SPRFAC
        write (lp,*) 'ATM_PATH ',trim(ATM_PATH)
        write (lp,*) 'RSTFRQ',RSTFRQ
        write (lp,*) 'RSTFMT',RSTFMT
        write (lp,*) 'RSTCMP',RSTCMP
        write (lp,*) 'IOTYPE',IOTYPE
c
c --- - determine number of regions for meridional overturning and flux
c --- - diagnostics and print namelist
        mer_nreg=1
        do while (mer_regnam(mer_nreg).ne.'')
          mer_nreg=mer_nreg+1
          if (mer_nreg.gt.odm) exit
        enddo
        mer_nreg=mer_nreg-1
c
        do m=1,mer_nreg
          mer_nflg(m)=1
          do while (mer_regflg(m,mer_nflg(m)).ne.-1)
            mer_nflg(m)=mer_nflg(m)+1
            if (mer_nflg(m).gt.rflgdm) exit
          enddo
          mer_nflg(m)=mer_nflg(m)-1
        enddo
c
        write (lp,*)
        write (lp,*) 'rdlim: BLOM MERDIA NAMELIST GROUP:'
        write (lp,*) 'MER_ORFILE ',trim(mer_orfile)
        write (lp,*) 'MER_MIFILE ',trim(mer_mifile)
        do m=1,mer_nreg
          write (lp,*)
          write (lp,*) 'MER_REGNAM ',trim(mer_regnam(m))
          write (lp,*) 'MER_MINLAT ',mer_minlat(m)
          write (lp,*) 'MER_MAXLAT ',mer_maxlat(m)
          write (lp,*) 'MER_REGFLG ',(mer_regflg(m,n),n=1,mer_nflg(m))
        enddo
c
        write (lp,*)
        write (lp,*) 'rdlim: BLOM SECDIA NAMELIST GROUP:'
        write (lp,*) 'SEC_SIFILE ',trim(sec_sifile)
c
c --- - determine number of io groups and print namelist
        nphy=0  
        do n=1,nphymax 
          if (GLB_AVEPERIO(n).ne.-999) nphy=nphy+1 
        enddo
        if (trxday.eq.0.) then
          H2D_SURRLX(1:nphy)=0
        endif
        if (srxday.eq.0.) then
          H2D_SALRLX(1:nphy)=0
        endif
        write (lp,*)
        write (lp,*) 'rdlim: BLOM DIAPHY NAMELIST GROUP:' 
        write (lp,*) 'GLB_FNAMETAG',GLB_FNAMETAG(1:nphy)
        write (lp,*) 'GLB_AVEPERIO',GLB_AVEPERIO(1:nphy)
        write (lp,*) 'GLB_FILEFREQ',GLB_FILEFREQ(1:nphy)
        write (lp,*) 'GLB_COMPFLAG',GLB_COMPFLAG(1:nphy)
        write (lp,*) 'GLB_NCFORMAT',GLB_NCFORMAT(1:nphy)
        write (lp,*) 'H2D_ABSWND  ',H2D_ABSWND(1:nphy)
        write (lp,*) 'H2D_ALB     ',H2D_ALB(1:nphy)
        write (lp,*) 'H2D_BTMSTR  ',H2D_BTMSTR(1:nphy)
        write (lp,*) 'H2D_BRNFLX  ',H2D_BRNFLX(1:nphy)
        write (lp,*) 'H2D_BRNPD   ',H2D_BRNPD(1:nphy)
        write (lp,*) 'H2D_DFL     ',H2D_DFL(1:nphy)
        write (lp,*) 'H2D_EVA     ',H2D_EVA(1:nphy)
        write (lp,*) 'H2D_FMLTFZ  ',H2D_FMLTFZ(1:nphy)
        write (lp,*) 'H2D_FICE    ',H2D_FICE(1:nphy)
        write (lp,*) 'H2D_HICE    ',H2D_HICE(1:nphy)
        write (lp,*) 'H2D_HMLTFZ  ',H2D_HMLTFZ(1:nphy)
        write (lp,*) 'H2D_HSNW    ',H2D_HSNW(1:nphy)
        write (lp,*) 'H2D_IAGE    ',H2D_IAGE(1:nphy)
        write (lp,*) 'H2D_IDKEDT  ',H2D_IDKEDT(1:nphy)
        write (lp,*) 'H2D_LIP     ',H2D_LIP(1:nphy)
        write (lp,*) 'H2D_MAXMLD  ',H2D_MAXMLD(1:nphy)
        write (lp,*) 'H2D_MLD     ',H2D_MLD(1:nphy)
        write (lp,*) 'H2D_MLDU    ',H2D_MLDU(1:nphy)
        write (lp,*) 'H2D_MLDV    ',H2D_MLDV(1:nphy)
        write (lp,*) 'H2D_MLTS    ',H2D_MLTS(1:nphy)
        write (lp,*) 'H2D_MLTSMN  ',H2D_MLTSMN(1:nphy)
        write (lp,*) 'H2D_MLTSMX  ',H2D_MLTSMX(1:nphy)
        write (lp,*) 'H2D_MLTSSQ  ',H2D_MLTSSQ(1:nphy)
        write (lp,*) 'H2D_MTKEUS  ',H2D_MTKEUS(1:nphy)
        write (lp,*) 'H2D_MTKENI  ',H2D_MTKENI(1:nphy)
        write (lp,*) 'H2D_MTKEBF  ',H2D_MTKEBF(1:nphy)
        write (lp,*) 'H2D_MTKERS  ',H2D_MTKERS(1:nphy)
        write (lp,*) 'H2D_MTKEPE  ',H2D_MTKEPE(1:nphy)
        write (lp,*) 'H2D_MTKEKE  ',H2D_MTKEKE(1:nphy)
        write (lp,*) 'H2D_MTY     ',H2D_MTY(1:nphy)
        write (lp,*) 'H2D_MXLU    ',H2D_MXLU(1:nphy)
        write (lp,*) 'H2D_MXLV    ',H2D_MXLV(1:nphy)
        write (lp,*) 'H2D_NSF     ',H2D_NSF(1:nphy)
        write (lp,*) 'H2D_PBOT    ',H2D_PBOT(1:nphy)
        write (lp,*) 'H2D_PSRF    ',H2D_PSRF(1:nphy)
        write (lp,*) 'H2D_RFIFLX  ',H2D_RFIFLX(1:nphy)
        write (lp,*) 'H2D_RNFFLX  ',H2D_RNFFLX(1:nphy)
        write (lp,*) 'H2D_SALFLX  ',H2D_SALFLX(1:nphy)
        write (lp,*) 'H2D_SALRLX  ',H2D_SALRLX(1:nphy)
        write (lp,*) 'H2D_SBOT    ',H2D_SBOT(1:nphy)
        write (lp,*) 'H2D_SEALV   ',H2D_SEALV(1:nphy)
        write (lp,*) 'H2D_SLVSQ   ',H2D_SLVSQ(1:nphy)
        write (lp,*) 'H2D_SFL     ',H2D_SFL(1:nphy)
        write (lp,*) 'H2D_SIGMX   ',H2D_SIGMX(1:nphy)
        write (lp,*) 'H2D_SOP     ',H2D_SOP(1:nphy)
        write (lp,*) 'H2D_SSS     ',H2D_SSS(1:nphy)
        write (lp,*) 'H2D_SSSSQ   ',H2D_SSSSQ(1:nphy)
        write (lp,*) 'H2D_SST     ',H2D_SST(1:nphy)
        write (lp,*) 'H2D_SSTSQ   ',H2D_SSTSQ(1:nphy)
        write (lp,*) 'H2D_SURFLX  ',H2D_SURFLX(1:nphy)
        write (lp,*) 'H2D_SURRLX  ',H2D_SURRLX(1:nphy)
        write (lp,*) 'H2D_SWA     ',H2D_SWA(1:nphy)
        write (lp,*) 'H2D_T20D    ',H2D_T20D(1:nphy)
        write (lp,*) 'H2D_TAUX    ',H2D_TAUX(1:nphy)
        write (lp,*) 'H2D_TAUY    ',H2D_TAUY(1:nphy)
        write (lp,*) 'H2D_TBOT    ',H2D_TBOT(1:nphy)
        write (lp,*) 'H2D_TICE    ',H2D_TICE(1:nphy)
        write (lp,*) 'H2D_TSRF    ',H2D_TSRF(1:nphy)
        write (lp,*) 'H2D_UB      ',H2D_UB(1:nphy)
        write (lp,*) 'H2D_UICE    ',H2D_UICE(1:nphy)
        write (lp,*) 'H2D_USTAR   ',H2D_USTAR(1:nphy)
        write (lp,*) 'H2D_USTAR3  ',H2D_USTAR3(1:nphy)
        write (lp,*) 'H2D_VB      ',H2D_VB(1:nphy)
        write (lp,*) 'H2D_VICE    ',H2D_VICE(1:nphy)
        write (lp,*) 'H2D_ZTX     ',H2D_ZTX(1:nphy)
        write (lp,*) 'LYR_BFSQ    ',LYR_BFSQ(1:nphy)
        write (lp,*) 'LYR_DIFDIA  ',LYR_DIFDIA(1:nphy)
        write (lp,*) 'LYR_DIFINT  ',LYR_DIFINT(1:nphy)
        write (lp,*) 'LYR_DIFISO  ',LYR_DIFISO(1:nphy)
        write (lp,*) 'LYR_DP      ',LYR_DP(1:nphy)
        write (lp,*) 'LYR_DZ      ',LYR_DZ(1:nphy)
        write (lp,*) 'LYR_SALN    ',LYR_SALN(1:nphy)
        write (lp,*) 'LYR_TEMP    ',LYR_TEMP(1:nphy)
        write (lp,*) 'LYR_TRC     ',LYR_TRC(1:nphy)
        write (lp,*) 'LYR_UFLX    ',LYR_UFLX(1:nphy)
        write (lp,*) 'LYR_UTFLX   ',LYR_UTFLX(1:nphy)
        write (lp,*) 'LYR_USFLX   ',LYR_USFLX(1:nphy)
        write (lp,*) 'LYR_UMFLTD  ',LYR_UMFLTD(1:nphy)
        write (lp,*) 'LYR_UTFLTD  ',LYR_UTFLTD(1:nphy)
        write (lp,*) 'LYR_UTFLLD  ',LYR_UTFLLD(1:nphy)
        write (lp,*) 'LYR_USFLTD  ',LYR_USFLTD(1:nphy)
        write (lp,*) 'LYR_USFLLD  ',LYR_USFLLD(1:nphy)
        write (lp,*) 'LYR_UVEL    ',LYR_UVEL(1:nphy)
        write (lp,*) 'LYR_VFLX    ',LYR_VFLX(1:nphy)
        write (lp,*) 'LYR_VTFLX   ',LYR_VTFLX(1:nphy)
        write (lp,*) 'LYR_VSFLX   ',LYR_VSFLX(1:nphy)
        write (lp,*) 'LYR_VMFLTD  ',LYR_VMFLTD(1:nphy)
        write (lp,*) 'LYR_VTFLTD  ',LYR_VTFLTD(1:nphy)
        write (lp,*) 'LYR_VTFLLD  ',LYR_VTFLLD(1:nphy)
        write (lp,*) 'LYR_VSFLTD  ',LYR_VSFLTD(1:nphy)
        write (lp,*) 'LYR_VSFLLD  ',LYR_VSFLLD(1:nphy)
        write (lp,*) 'LYR_VVEL    ',LYR_VVEL(1:nphy)
        write (lp,*) 'LYR_WFLX    ',LYR_WFLX(1:nphy)
        write (lp,*) 'LYR_WFLX2   ',LYR_WFLX2(1:nphy)
        write (lp,*) 'LYR_PV      ',LYR_PV(1:nphy)
        write (lp,*) 'LYR_TKE     ',LYR_TKE(1:nphy)
        write (lp,*) 'LYR_GLS_PSI ',LYR_GLS_PSI(1:nphy)
        write (lp,*) 'LYR_IDLAGE  ',LYR_IDLAGE(1:nphy)
        write (lp,*) 'LVL_BFSQ    ',LVL_BFSQ(1:nphy)
        write (lp,*) 'LVL_SALN    ',LVL_SALN(1:nphy)
        write (lp,*) 'LVL_TEMP    ',LVL_TEMP(1:nphy)
        write (lp,*) 'LVL_TRC     ',LVL_TRC(1:nphy)
        write (lp,*) 'LVL_UFLX    ',LVL_UFLX(1:nphy)
        write (lp,*) 'LVL_UTFLX   ',LVL_UTFLX(1:nphy)
        write (lp,*) 'LVL_USFLX   ',LVL_USFLX(1:nphy)
        write (lp,*) 'LVL_UMFLTD  ',LVL_UMFLTD(1:nphy)
        write (lp,*) 'LVL_UTFLTD  ',LVL_UTFLTD(1:nphy)
        write (lp,*) 'LVL_UTFLLD  ',LVL_UTFLLD(1:nphy)
        write (lp,*) 'LVL_USFLTD  ',LVL_USFLTD(1:nphy)
        write (lp,*) 'LVL_USFLLD  ',LVL_USFLLD(1:nphy)
        write (lp,*) 'LVL_UVEL    ',LVL_UVEL(1:nphy)
        write (lp,*) 'LVL_VFLX    ',LVL_VFLX(1:nphy)
        write (lp,*) 'LVL_VTFLX   ',LVL_VTFLX(1:nphy)
        write (lp,*) 'LVL_VSFLX   ',LVL_VSFLX(1:nphy)
        write (lp,*) 'LVL_VMFLTD  ',LVL_VMFLTD(1:nphy)
        write (lp,*) 'LVL_VTFLTD  ',LVL_VTFLTD(1:nphy)
        write (lp,*) 'LVL_VTFLLD  ',LVL_VTFLLD(1:nphy)
        write (lp,*) 'LVL_VSFLTD  ',LVL_VSFLTD(1:nphy)
        write (lp,*) 'LVL_VSFLLD  ',LVL_VSFLLD(1:nphy)
        write (lp,*) 'LVL_VVEL    ',LVL_VVEL(1:nphy)
        write (lp,*) 'LVL_WFLX    ',LVL_WFLX(1:nphy)
        write (lp,*) 'LVL_WFLX2   ',LVL_WFLX2(1:nphy)
        write (lp,*) 'LVL_PV      ',LVL_PV(1:nphy)
        write (lp,*) 'LVL_TKE     ',LVL_TKE(1:nphy)
        write (lp,*) 'LVL_GLS_PSI ',LVL_GLS_PSI(1:nphy)
        write (lp,*) 'LVL_IDLAGE  ',LVL_IDLAGE(1:nphy)
        write (lp,*) 'MSC_MMFLXL  ',MSC_MMFLXL(1:nphy)
        write (lp,*) 'MSC_MMFLXD  ',MSC_MMFLXD(1:nphy)
        write (lp,*) 'MSC_MMFTDL  ',MSC_MMFTDL(1:nphy)
        write (lp,*) 'MSC_MMFTDD  ',MSC_MMFTDD(1:nphy)
        write (lp,*) 'MSC_MHFLX   ',MSC_MHFLX(1:nphy)
        write (lp,*) 'MSC_MHFTD   ',MSC_MHFTD(1:nphy)
        write (lp,*) 'MSC_MHFLD   ',MSC_MHFLD(1:nphy)
        write (lp,*) 'MSC_MSFLX   ',MSC_MSFLX(1:nphy)
        write (lp,*) 'MSC_MSFTD   ',MSC_MSFTD(1:nphy)
        write (lp,*) 'MSC_MSFLD   ',MSC_MSFLD(1:nphy)
        write (lp,*) 'MSC_VOLTR   ',MSC_VOLTR(1:nphy)
        write (lp,*) 'MSC_MASSGS  ',MSC_MASSGS(1:nphy)
        write (lp,*) 'MSC_VOLGS   ',MSC_VOLGS(1:nphy)
        write (lp,*) 'MSC_SALNGA  ',MSC_SALNGA(1:nphy)
        write (lp,*) 'MSC_TEMPGA  ',MSC_TEMPGA(1:nphy)
        write (lp,*) 'MSC_SSSGA   ',MSC_SSSGA(1:nphy)
        write (lp,*) 'MSC_SSTGA   ',MSC_SSTGA(1:nphy)
        write (lp,*)
c
      endif
c
c --- broadcast variables set by namelists
c
      call xcbcst(nday1)
      call xcbcst(nday2)
      call xcbcst(idate)
      call xcbcst(idate0)
      call xcbcst(runid)
      call xcbcst(expcnf)
      call xcbcst(runtyp)
      call xcbcst(grfile)
      call xcbcst(icfile)
      call xcbcst(baclin)
      call xcbcst(batrop)
      call xcbcst(mdv2hi)
      call xcbcst(mdv2lo)
      call xcbcst(mdv4hi)
      call xcbcst(mdv4lo)
      call xcbcst(mdc2hi)
      call xcbcst(mdc2lo)
      call xcbcst(vsc2hi)
      call xcbcst(vsc2lo)
      call xcbcst(vsc4hi)
      call xcbcst(vsc4lo)
      call xcbcst(cbar)
      call xcbcst(cb)
      call xcbcst(cwbdts)
      call xcbcst(cwbdls)
      call xcbcst(mommth)
      call xcbcst(eitmth)
      call xcbcst(edritp)
      call xcbcst(bmcmth)
      call xcbcst(rmpmth)
      call xcbcst(edwmth)
      call xcbcst(mlrttp)
      call xcbcst(edsprs)
      call xcbcst(egc)
      call xcbcst(eggam)
      call xcbcst(eglsmn)
      call xcbcst(egmndf)
      call xcbcst(egmxdf)
      call xcbcst(egidfq)
      call xcbcst(ri0)
      call xcbcst(rm0)
      call xcbcst(rm5)
      call xcbcst(ce)
      call xcbcst(bdmtyp)
      call xcbcst(bdmc1)
      call xcbcst(bdmc2)
      call xcbcst(tdfile)
      call xcbcst(tkepf)
      call xcbcst(niwgf)
      call xcbcst(niwbf)
      call xcbcst(niwlf)
      call xcbcst(swamth)
      call xcbcst(jwtype)
      call xcbcst(chlopt)
      call xcbcst(ccfile)
      call xcbcst(trxday)
      call xcbcst(srxday)
      call xcbcst(trxdpt)
      call xcbcst(srxdpt)
      call xcbcst(trxlim)
      call xcbcst(srxlim)
      call xcbcst(aptflx)
      call xcbcst(apsflx)
      call xcbcst(ditflx)
      call xcbcst(disflx)
      call xcbcst(srxbal)
      call xcbcst(scfile)
      call xcbcst(smtfrc)
      call xcbcst(sprfac)
      call xcbcst(atm_path)
      call xcbcst(rstfrq)
      call xcbcst(rstfmt)
      call xcbcst(rstcmp)
      call xcbcst(iotype)
c
      call xcbcst(mer_nreg)
      do m=1,mer_nreg
        call xcbcst(mer_regnam(m))
      enddo
c
      call xcbcst(H2D_ABSWND)
      call xcbcst(H2D_ALB)
      call xcbcst(H2D_BTMSTR)
      call xcbcst(H2D_BRNFLX)
      call xcbcst(H2D_BRNPD)
      call xcbcst(H2D_DFL)
      call xcbcst(H2D_EVA)
      call xcbcst(H2D_FICE)
      call xcbcst(H2D_FMLTFZ)
      call xcbcst(H2D_HICE)
      call xcbcst(H2D_HMLTFZ)
      call xcbcst(H2D_HSNW)
      call xcbcst(H2D_IAGE)
      call xcbcst(H2D_IDKEDT)
      call xcbcst(H2D_LIP)
      call xcbcst(H2D_MAXMLD)
      call xcbcst(H2D_MLD)
      call xcbcst(H2D_MLDU)
      call xcbcst(H2D_MLDV)
      call xcbcst(H2D_MLTS)
      call xcbcst(H2D_MLTSMN)
      call xcbcst(H2D_MLTSMX)
      call xcbcst(H2D_MLTSSQ)
      call xcbcst(H2D_MTKEUS)
      call xcbcst(H2D_MTKENI)
      call xcbcst(H2D_MTKEBF)
      call xcbcst(H2D_MTKERS)
      call xcbcst(H2D_MTKEPE)
      call xcbcst(H2D_MTKEKE)
      call xcbcst(H2D_MTY)
      call xcbcst(H2D_MXLU)
      call xcbcst(H2D_MXLV)
      call xcbcst(H2D_NSF)
      call xcbcst(H2D_PBOT)
      call xcbcst(H2D_PSRF)
      call xcbcst(H2D_RFIFLX)
      call xcbcst(H2D_RNFFLX)
      call xcbcst(H2D_SALFLX)
      call xcbcst(H2D_SALRLX)
      call xcbcst(H2D_SBOT)
      call xcbcst(H2D_SEALV)
      call xcbcst(H2D_SLVSQ)
      call xcbcst(H2D_SFL)
      call xcbcst(H2D_SOP)
      call xcbcst(H2D_SIGMX)
      call xcbcst(H2D_SSS)
      call xcbcst(H2D_SSSSQ)
      call xcbcst(H2D_SST)
      call xcbcst(H2D_SSTSQ)
      call xcbcst(H2D_SURFLX)
      call xcbcst(H2D_SURRLX)
      call xcbcst(H2D_SWA)
      call xcbcst(H2D_T20D)
      call xcbcst(H2D_TAUX)
      call xcbcst(H2D_TAUY)
      call xcbcst(H2D_TBOT)
      call xcbcst(H2D_TICE)
      call xcbcst(H2D_TSRF)
      call xcbcst(H2D_UB)
      call xcbcst(H2D_UICE)
      call xcbcst(H2D_USTAR)
      call xcbcst(H2D_USTAR3)
      call xcbcst(H2D_VB)
      call xcbcst(H2D_VICE)
      call xcbcst(H2D_ZTX)
      call xcbcst(LYR_BFSQ)
      call xcbcst(LYR_DIFDIA)
      call xcbcst(LYR_DIFINT)
      call xcbcst(LYR_DIFISO)
      call xcbcst(LYR_DP)
      call xcbcst(LYR_DPU)
      call xcbcst(LYR_DPV)
      call xcbcst(LYR_DZ)
      call xcbcst(LYR_SALN)
      call xcbcst(LYR_TEMP)
      call xcbcst(LYR_TRC)
      call xcbcst(LYR_UFLX)
      call xcbcst(LYR_UTFLX)
      call xcbcst(LYR_USFLX)
      call xcbcst(LYR_UMFLTD)
      call xcbcst(LYR_UTFLTD)
      call xcbcst(LYR_UTFLLD)
      call xcbcst(LYR_USFLTD)
      call xcbcst(LYR_USFLLD)
      call xcbcst(LYR_UVEL)
      call xcbcst(LYR_VFLX)
      call xcbcst(LYR_VTFLX)
      call xcbcst(LYR_VSFLX)
      call xcbcst(LYR_VMFLTD)
      call xcbcst(LYR_VTFLTD)
      call xcbcst(LYR_VTFLLD)
      call xcbcst(LYR_VSFLTD)
      call xcbcst(LYR_VSFLLD)
      call xcbcst(LYR_VVEL)
      call xcbcst(LYR_WFLX)
      call xcbcst(LYR_WFLX2)
      call xcbcst(LYR_PV)
      call xcbcst(LYR_TKE)
      call xcbcst(LYR_GLS_PSI)
      call xcbcst(LYR_IDLAGE)
      call xcbcst(LVL_BFSQ)
      call xcbcst(LVL_DIFDIA)
      call xcbcst(LVL_DIFINT)
      call xcbcst(LVL_DIFISO)
      call xcbcst(LVL_DZ)
      call xcbcst(LVL_SALN)
      call xcbcst(LVL_TEMP)
      call xcbcst(LVL_TRC)
      call xcbcst(LVL_UFLX)
      call xcbcst(LVL_UTFLX)
      call xcbcst(LVL_USFLX)
      call xcbcst(LVL_UMFLTD)
      call xcbcst(LVL_UTFLTD)
      call xcbcst(LVL_UTFLLD)
      call xcbcst(LVL_USFLTD)
      call xcbcst(LVL_USFLLD)
      call xcbcst(LVL_UVEL)
      call xcbcst(LVL_VFLX)
      call xcbcst(LVL_VTFLX)
      call xcbcst(LVL_VSFLX)
      call xcbcst(LVL_VMFLTD)
      call xcbcst(LVL_VTFLTD)
      call xcbcst(LVL_VTFLLD)
      call xcbcst(LVL_VSFLTD)
      call xcbcst(LVL_VSFLLD)
      call xcbcst(LVL_VVEL)
      call xcbcst(LVL_WFLX)
      call xcbcst(LVL_WFLX2)
      call xcbcst(LVL_PV)
      call xcbcst(LVL_TKE)
      call xcbcst(LVL_GLS_PSI)
      call xcbcst(LVL_IDLAGE)
      call xcbcst(MSC_MMFLXL)
      call xcbcst(MSC_MMFLXD)
      call xcbcst(MSC_MMFTDL)
      call xcbcst(MSC_MMFTDD)
      call xcbcst(MSC_MHFLX)
      call xcbcst(MSC_MHFTD)
      call xcbcst(MSC_MHFLD)
      call xcbcst(MSC_MSFLX)
      call xcbcst(MSC_MSFTD)
      call xcbcst(MSC_MSFLD)
      call xcbcst(MSC_VOLTR)
      call xcbcst(MSC_MASSGS)
      call xcbcst(MSC_VOLGS)
      call xcbcst(MSC_SALNGA)
      call xcbcst(MSC_TEMPGA)
      call xcbcst(MSC_SSSGA)
      call xcbcst(MSC_SSTGA)
      call xcbcst(GLB_AVEPERIO)
      call xcbcst(GLB_FILEFREQ)
      call xcbcst(GLB_COMPFLAG)
      call xcbcst(GLB_NCFORMAT)
      do n=1,nphymax 
        call xcbcst(GLB_FNAMETAG(n))
      enddo
c
      call xcbcst(nphy)
c
c --- convert integer dates 
      nyear=sign(abs(idate)/10000,idate)
      nmonth=abs(idate)/100-abs(nyear)*100
      nday=abs(idate)-abs(nyear)*10000-nmonth*100
      nyear0=sign(abs(idate0)/10000,idate0)
      nmonth0=abs(idate0)/100-abs(nyear0)*100
      nday0=abs(idate0)-abs(nyear0)*10000-nmonth0*100
c
c --- set atm_path length
      if (expcnf.eq.'ben02syn') then
        atm_path_len=1
        do while (atm_path_len.lt.80.and.
     .            atm_path(atm_path_len:atm_path_len).ne.' ')
          atm_path_len=atm_path_len+1
        enddo
        atm_path_len=atm_path_len-1
      endif
c
c --- verify integer number of baroclinic time steps per day
      if (mod(86400.+epsil,baclin).gt.2.*epsil) then
        if (mnproc.eq.1) then
          write (lp,*) 'rdlim: must have an integer number of '//
     .                 'baroclinic time steps pr. day!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
c --- get baroclinic time steps in a day
      nstep_in_day=nint(86400./baclin)
      if (mnproc.eq.1) then
        write (lp,*) 'rdlim: number of baroclinic time steps per '//
     .               'day:              ',nstep_in_day
      endif
c
c --- Initialize the calendar
      call inicln
c
      if (expcnf.eq.'cesm') then
c
c --- - override namelist experiment id with the one received from
c --- - coupler
        runid=runid_mct
c
c --- - verify integer number of baroclinic time steps per coupling
c --- - interval
        if (mod(ocn_cpl_dt_mct+epsil,baclin).gt.2.*epsil) then
          if (mnproc.eq.1) then
            write (lp,*) 'rdlim: must have an integer number of '//
     .                   'baroclinic time steps pr. coupling'
            write (lp,*) '       interval!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
c
c --- - get time step and correct model date after first coupling
c --- - interval
        nstep_in_cpl=nint(ocn_cpl_dt_mct/baclin)
        if (mnproc.eq.1) then
          write (lp,*) 'rdlim: number of baroclinic time steps per '//
     .                 'coupling interval:',nstep_in_cpl
        endif
        nstep0=0
        do
          nstep0=nstep0+1
          call updcln(nstep0)
          if (nstep0.eq.nstep_in_cpl) exit
        enddo
c
        if (runtyp.eq.'startup') then
c
c --- --- when runtyp equal 'startup' the ocean integration start after
c --- --- first coupling interval
c
          time0=0.
          nday1=0
          nstep1=nstep0
          nstep2=nstep1
c
        else
c
c --- --- for runtyp equal 'hybrid', 'branch' or 'continue' a
c --- --- 'rpointer.ocn' file containing the path to a valid restart
c --- --- file is expected and integration time is retrieved from
c --- --- restart file
c
          if (mnproc.eq.1) inquire(file='rpointer.ocn',exist=fexist)
          call xcbcst(fexist)
          if (.not.fexist) then
            if (mnproc.eq.1) then
              write (lp,*) 'rdlim: could not find rpointer.ocn file!'
            endif
            call xcstop('(rdlim)')
                   stop '(rdlim)'
          endif
          if (mnproc.eq.1) then
            open (unit=nfu,file='rpointer.ocn')
            read (nfu,'(a)') rstfnm
            close (unit=nfu)
            inquire(file=rstfnm,exist=fexist)
          endif
          call xcbcst(fexist)
          if (.not.fexist) then
            if (mnproc.eq.1) then
              write (lp,*) 'rdlim: could not find restart file!'
            endif
            call xcstop('(rdlim)')
                   stop '(rdlim)'
          endif
          call xcbcst(rstfnm)
          call ncfopn(rstfnm,'r',' ',1,iotype)
          call ncgeti('nday0',nd0rst)
          call ncgeti('nmonth0',nm0rst)
          call ncgeti('nyear0',ny0rst)
          call ncgetr('time0',time0)
          call ncgetr('time',time)
          nday1=nint(time-time0)
c
          if (runtyp.eq.'hybrid') then
c 
c --- ----- When runtyp equal 'hybrid' the ocean integration starts
c --- ----- after first coupling interval.
            n=1
            do while (n.lt.256.and.rstfnm(n:n).ne.' ')
              n=n+1
            enddo
            n=n-1
            read (rstfnm(n-18:n-15),'(i4)') nyear
            read (rstfnm(n-13:n-12),'(i2)') nmonth
            read (rstfnm(n-10:n-9 ),'(i2)') nday
            if (nyear.eq.nyear0.and.
     .          nmonth.eq.nmonth0.and.
     .          nday.eq.nday0) then
              nday0=nd0rst
              nmonth0=nm0rst
              nyear0=ny0rst
              nstep1=nday1*nstep_in_day+nstep0
            else
              time0=time
              nday1=0
              nstep1=nstep0
              nstep2=nstep1
              nyear=nyear0
              nmonth=nmonth0
              nday=nday0
            endif
c
          else
c
            nday0=nd0rst
            nmonth0=nm0rst
            nyear0=ny0rst
            nstep1=nday1*nstep_in_day
c
c --- ----- For runtypes 'branch' or 'continue' override namelist date
c --- ----- with date extracted from restart file name
c
            n=1
            do while (n.lt.256.and.rstfnm(n:n).ne.' ')
              n=n+1
            enddo
            n=n-1
            read (rstfnm(n-18:n-15),'(i4)') nyear
            read (rstfnm(n-13:n-12),'(i2)') nmonth
            read (rstfnm(n-10:n-9 ),'(i2)') nday
            call inicln
c
          endif
c
          call ncfcls
c
          nstep2=nstep1
c
        endif
c
      else
c
        if (nday1.lt.0.or.nday2.lt.0) then
          if (mnproc.eq.1) then
            write (lp,*) 'rdlim: integrations days must be positive!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
c
        if (nday2.lt.nday1) then
          if (mnproc.eq.1) then
            write (lp,*) 'rdlim: first day of integration must be '//
     .                   'less than or equal to last day!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
c
        if (nday1.eq.0.and.
     .      (nyear.ne.nyear0.or.nmonth.ne.nmonth0.or.
     .       nday.ne.nday0)) then
          if (mnproc.eq.1) then
            write (lp,*) 'rdlim: when first integration day is zero, '//
     .                   'model date and initial '
            write (lp,*) '       experiment date must be equal!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
c
        if (nyear.lt.nyear0.or.
     .      (nyear.eq.nyear0.and.nmonth.lt.nmonth0).or.
     .      (nyear.eq.nyear0.and.nmonth.eq.nmonth0.and.
     .       nday.lt.nday0)) then
          if (mnproc.eq.1) then
            write (lp,*) 'rdlim: model date must be greater or equal '//
     .                   'to initial experiment date!'
          endif
          call xcstop('(rdlim)')
                 stop '(rdlim)'
        endif
c
c --- - model is to be integrated from time step 'nstep1' to 'nstep2'
        time0=0.
        nstep1=nday1*nstep_in_day
        nstep2=nday2*nstep_in_day
c
        if (csdiag) then
          nstep2=nstep1+2
        endif
c
      endif
c
c --- set experiment id length
      runid_len=1
      do while (runid_len.lt.256.and.runid(runid_len:runid_len).ne.' ')
        runid_len=runid_len+1
      enddo
      runid_len=runid_len-1
c
      if (trxday.eq.0.and.ditflx) then
        if (mnproc.eq.1) then
          write (lp,*)
     .      'rdlim: trxday=0. and ditflx=.true.. Inconsistent!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
      if (srxday.eq.0.and.disflx) then
        if (mnproc.eq.1) then
          write (lp,*)
     .      'rdlim: srxday=0. and disflx=.true.. Inconsistent!'
        endif
        call xcstop('(rdlim)')
               stop '(rdlim)'
      endif
c
c --- represent time between restarts in time steps
      rstmon=.false.
      rstann=.false.
      if (nint(rstfrq).eq.30) then
        rstmon=.true.
      elseif (nint(rstfrq).eq.365) then
        rstann=.true.
      endif
      rstfrq=nstep_in_day*max(1.,rstfrq)
c
c --- represent time between diagnostics in time steps
      do n=1,nphy 
        GLB_FILEFREQ(n)=max(GLB_AVEPERIO(n),GLB_FILEFREQ(n))
c
        if (GLB_AVEPERIO(n).lt.0) then
          diagfq_phy(n)=-real(nstep_in_day)/GLB_AVEPERIO(n)
        else
          diagfq_phy(n)=nstep_in_day*max(1,GLB_AVEPERIO(n))
        endif
        diagmon_phy(n)=.false.
        diagann_phy(n)=.false.
        if (GLB_AVEPERIO(n).eq.30) then
          diagmon_phy(n)=.true.
        elseif (GLB_AVEPERIO(n).eq.365) then
          diagann_phy(n)=.true.
        endif
c
        if (GLB_FILEFREQ(n).lt.0) then
          filefq_phy(n)=-real(nstep_in_day)/GLB_FILEFREQ(n)
        else
          filefq_phy(n)=nstep_in_day*max(1,GLB_FILEFREQ(n))
        endif
        filemon_phy(n)=.false.
        fileann_phy(n)=.false.
        if (GLB_FILEFREQ(n).eq.30) then
          filemon_phy(n)=.true.
        elseif (GLB_FILEFREQ(n).eq.365) then
          fileann_phy(n)=.true.
        endif
      enddo
c
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
      lstep=2*ceiling(.5*baclin/batrop)
      dlt=baclin/lstep
c
      if (mnproc.eq.1) then
        write (lp,*)
        write (lp,*) 'rdlim: resolved parameters:'
        write (lp,'(a,i10)') ' nday1:                          ',nday1
        write (lp,'(a,i10)') ' nday2:                          ',nday2
        write (lp,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5)')
     .    ' nyear ,nmonth ,nday ,nsec :     ',
     .     nyear ,' ,',nmonth ,' ,',nday ,' ,',
     .     nint(mod(nstep1,nstep_in_day)*baclin)
        write (lp,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5)')
     .    ' nyear0,nmonth0,nday0,nsec0:     ',
     .     nyear0,' ,',nmonth0,' ,',nday0,' ,',0
        write (lp,'(a,i10)') ' barotr. steps per barocl. step: ',lstep
        write (lp,'(2a)') ' runid:                          ',
     .    trim(runid)
        write (lp,'(a,f10.4)') ' rstfrq:                         ',
     .    rstfrq
        write (lp,'(a,l10)') ' rstmon:                         ',rstmon
        write (lp,'(a,l10)') ' rstann:                         ',rstann
        write (lp,*)
        call flush(lp)
      endif
c
      return
      end
