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
      integer, save :: rstcmp,rstfmt
      logical, save :: rstmon,rstann
c
c --- Copies of MICOM variables that are used for HAMOCC diagnostics
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  pbath,ubath,vbath 
      integer, save :: nstepinday 
 
c --- 2d and 3d diagnostic variables     
      integer, save :: nphysrf,nphylyr,nphylvl 
      real, save, allocatable, dimension(:,:,:) :: physrf
      real, save, allocatable, dimension(:,:,:,:) :: phylyr,phylvl 
c
c --- Levitus levels 
#ifndef LEVITUS2X 
      integer, parameter :: ddm=35,k350=12,w350=1.  
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
      integer, parameter :: ddm=70,k350=25,w350=0.5
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
c --- Zonally integrated diagnostics 
      integer, save :: lmax 
      integer, parameter :: ldm=itdm+jtdm,sdm=ldm,odm=3,slenmax=50 
      real, save, dimension(ldm) :: mtlat
      real, save, dimension(kdm) :: sigmar1 
      real, save, dimension(ldm,kdm,odm) :: mmflxl,mmftdl
      real, save, dimension(ldm,ddm,odm) :: mmflxd,mmftdd
      real, save, dimension(ldm,odm) ::
     .  mhflx,mhftd,mhfld,msflx,msftd,msfld
      character(len=slenmax), save, dimension(3) :: region_name
c
c --- Section transports 
      integer, save :: sec_num
      integer, parameter :: max_sec=400
      character(len=slenmax), save :: sec_name(max_sec)
      real, save, dimension(max_sec) :: voltr
c
c --- Namelist 
      integer, dimension(nphymax), save ::
     . SRF_ABSWND ,SRF_ALB    ,SRF_BRNFLX ,SRF_BRNPD  ,SRF_DFL    ,
     . SRF_EVA    ,SRF_FICE   ,SRF_FMLTFZ ,SRF_HICE   ,SRF_HMLTFZ ,
     . SRF_HSNW   ,SRF_IAGE   ,SRF_LIP    ,SRF_MAXMLD ,SRF_MLD    ,
     . SRF_MLDU   ,SRF_MLDV   ,SRF_MTY    ,SRF_MXLU   ,SRF_MXLV   ,
     . SRF_NSF    ,SRF_RFIFLX ,SRF_RNFFLX ,SRF_SALFLX ,SRF_SALRLX ,
     . SRF_SEALV  ,SRF_SFL    ,SRF_SOP    ,SRF_SIGMX  ,SRF_SSS    ,
     . SRF_SST    ,SRF_SURFLX ,SRF_SURRLX ,SRF_SWA    ,SRF_TAUX   ,
     . SRF_TAUY   ,SRF_TICE   ,SRF_TSRF   ,SRF_UB     ,SRF_UICE   ,
     . SRF_USTAR  ,SRF_VB     ,SRF_VICE   ,SRF_ZTX    ,
     . LYR_DIAFLX ,LYR_DIFDIA ,LYR_DIFINT ,LYR_DIFISO ,LYR_DP     ,
     . LYR_DPU    ,LYR_DPV    ,LYR_DZ     ,LYR_SALN   ,LYR_TEMP   ,
     . LYR_TRC    ,LYR_UFLX   ,LYR_UTFLX  ,LYR_USFLX  ,LYR_UMFLTD ,
     . LYR_UTFLTD ,LYR_UTFLLD ,LYR_USFLTD ,LYR_USFLLD ,LYR_UVEL   ,
     . LYR_VFLX   ,LYR_VTFLX  ,LYR_VSFLX  ,LYR_VMFLTD ,LYR_VTFLTD ,
     . LYR_VTFLLD ,LYR_VSFLTD ,LYR_VSFLLD ,LYR_VVEL   ,LYR_WFLX   ,
     . LYR_WFLX2  ,LYR_TKE    ,LYR_GLS_PSI,LYR_IDLAGE ,
     . LVL_DIAFLX ,LVL_DIFDIA ,LVL_DIFINT ,LVL_DIFISO ,LVL_DZ     ,
     . LVL_SALN   ,LVL_TEMP   ,LVL_TRC    ,LVL_UFLX   ,LVL_UTFLX  ,
     . LVL_USFLX  ,LVL_UMFLTD ,LVL_UTFLTD ,LVL_UTFLLD ,LVL_USFLTD ,
     . LVL_USFLLD ,LVL_UVEL   ,LVL_VFLX   ,LVL_VTFLX  ,LVL_VSFLX  ,
     . LVL_VMFLTD ,LVL_VTFLTD ,LVL_VTFLLD ,LVL_VSFLTD ,LVL_VSFLLD ,
     . LVL_VVEL   ,LVL_WFLX   ,LVL_WFLX2  ,LVL_TKE    ,LVL_GLS_PSI,
     . LVL_IDLAGE ,
     . MSC_MMFLXL ,MSC_MMFLXD ,MSC_MMFTDL ,MSC_MMFTDD ,MSC_MHFLX  ,
     . MSC_MHFTD  ,MSC_MHFLD  ,MSC_MSFLX  ,MSC_MSFTD  ,MSC_MSFLD  ,
     . MSC_VOLTR  , 
     . GLB_AVEPERIO,GLB_FILEFREQ,GLB_COMPFLAG,GLB_NCFORMAT
      character(len=10), dimension(nphymax), save :: GLB_FNAMETAG
      integer, dimension(nphymax), save ::
     . ACC_ABSWND ,ACC_ALB    ,ACC_BRNFLX ,ACC_BRNPD  ,ACC_DFL    ,
     . ACC_EVA    ,ACC_FICE   ,ACC_FMLTFZ ,ACC_HICE   ,ACC_HMLTFZ ,
     . ACC_HSNW   ,ACC_IAGE   ,ACC_LIP    ,ACC_MAXMLD ,ACC_MLD    ,
     . ACC_MLDU   ,ACC_MLDV   ,ACC_MTY    ,ACC_MXLU   ,ACC_MXLV   ,
     . ACC_NSF    ,ACC_RFIFLX ,ACC_RNFFLX ,ACC_SALFLX ,ACC_SALRLX ,
     . ACC_SEALV  ,ACC_SFL    ,ACC_SOP    ,ACC_SIGMX  ,ACC_SSS    ,
     . ACC_SST    ,ACC_SURFLX ,ACC_SURRLX ,ACC_SWA    ,ACC_TAUX   ,
     . ACC_TAUY   ,ACC_TICE   ,ACC_TSRF   ,ACC_UB     ,ACC_UICE   ,
     . ACC_USTAR  ,ACC_VB     ,ACC_VICE   ,ACC_ZTX    ,
     . ACC_IVOLU  ,ACC_IVOLV  ,ACC_UTIL   ,
     . ACC_DIAFLX ,ACC_DIFDIA ,ACC_DIFINT ,ACC_DIFISO ,ACC_DP     ,
     . ACC_DPU    ,ACC_DPV    ,ACC_DZ     ,ACC_SALN   ,ACC_TEMP   ,
     . ACC_UFLX   ,ACC_UTFLX  ,ACC_USFLX  ,ACC_UMFLTD ,ACC_UTFLTD ,
     . ACC_UTFLLD ,ACC_USFLTD ,ACC_USFLLD ,ACC_UVEL   ,ACC_VFLX   ,
     . ACC_VTFLX  ,ACC_VSFLX  ,ACC_VMFLTD ,ACC_VTFLTD ,ACC_VTFLLD ,
     . ACC_VSFLTD ,ACC_VSFLLD ,ACC_VVEL   ,ACC_WFLX   ,ACC_WFLX2  ,
     . ACC_TKE    ,ACC_GLS_PSI,
     . ACC_DIAFLXLVL,ACC_DIFDIALVL,ACC_DIFINTLVL,ACC_DIFISOLVL,
     . ACC_DZLVL    ,ACC_SALNLVL  ,ACC_TEMPLVL  ,ACC_UFLXLVL  ,
     . ACC_UTFLXLVL ,ACC_USFLXLVL ,ACC_UMFLTDLVL,ACC_UTFLTDLVL,
     . ACC_UTFLLDLVL,ACC_USFLTDLVL,ACC_USFLLDLVL,ACC_UVELLVL  ,
     . ACC_VFLXLVL  ,ACC_VTFLXLVL ,ACC_VSFLXLVL ,ACC_VMFLTDLVL,
     . ACC_VTFLTDLVL,ACC_VTFLLDLVL,ACC_VSFLTDLVL,ACC_VSFLLDLVL,
     . ACC_VVELLVL  ,ACC_WFLXLVL  ,ACC_WFLX2LVL ,ACC_TKELVL   ,
     . ACC_GLS_PSILVL,ACC_UFLXOLD ,ACC_VFLXOLD  ,ACC_UTILLVL  ,     
     . ACC_MMFLXL ,ACC_MMFLXD ,ACC_MMFTDL ,ACC_MMFTDD ,ACC_MHFLX  ,
     . ACC_MHFTD  ,ACC_MHFLD  ,ACC_MSFLX  ,ACC_MSFTD  ,ACC_MSFLD  ,
     . ACC_VOLTR 
      namelist /DIAPHY/
     . SRF_ABSWND ,SRF_ALB    ,SRF_BRNFLX ,SRF_BRNPD  ,SRF_DFL    ,
     . SRF_EVA    ,SRF_FICE   ,SRF_FMLTFZ ,SRF_HICE   ,SRF_HMLTFZ ,
     . SRF_HSNW   ,SRF_IAGE   ,SRF_LIP    ,SRF_MAXMLD ,SRF_MLD    ,
     . SRF_MLDU   ,SRF_MLDV   ,SRF_MTY    ,SRF_MXLU   ,SRF_MXLV   ,
     . SRF_NSF    ,SRF_RFIFLX ,SRF_RNFFLX ,SRF_SALFLX ,SRF_SALRLX ,
     . SRF_SEALV  ,SRF_SFL    ,SRF_SOP    ,SRF_SIGMX  ,SRF_SSS    ,
     . SRF_SST    ,SRF_SURFLX ,SRF_SURRLX ,SRF_SWA    ,SRF_TAUX   ,
     . SRF_TAUY   ,SRF_TICE   ,SRF_TSRF   ,SRF_UB     ,SRF_UICE   ,
     . SRF_USTAR  ,SRF_VB     ,SRF_VICE   ,SRF_ZTX    ,
     . LYR_DIAFLX ,LYR_DIFDIA ,LYR_DIFINT ,LYR_DIFISO ,LYR_DP     ,
     . LYR_DPU    ,LYR_DPV    ,LYR_DZ     ,LYR_SALN   ,LYR_TEMP   ,
     . LYR_TRC    ,LYR_UFLX   ,LYR_UTFLX  ,LYR_USFLX  ,LYR_UMFLTD ,
     . LYR_UTFLTD ,LYR_UTFLLD ,LYR_USFLTD ,LYR_USFLLD ,LYR_UVEL   ,
     . LYR_VFLX   ,LYR_VTFLX  ,LYR_VSFLX  ,LYR_VMFLTD ,LYR_VTFLTD ,
     . LYR_VTFLLD ,LYR_VSFLTD ,LYR_VSFLLD ,LYR_VVEL   ,LYR_WFLX   ,
     . LYR_WFLX2  ,LYR_TKE    ,LYR_GLS_PSI,LYR_IDLAGE ,
     . LVL_DIAFLX ,LVL_DIFDIA ,LVL_DIFINT ,LVL_DIFISO ,LVL_DZ     ,
     . LVL_SALN   ,LVL_TEMP   ,LVL_TRC    ,LVL_UFLX   ,LVL_UTFLX  ,
     . LVL_USFLX  ,LVL_UMFLTD ,LVL_UTFLTD ,LVL_UTFLLD ,LVL_USFLTD ,
     . LVL_USFLLD ,LVL_UVEL   ,LVL_VFLX   ,LVL_VTFLX  ,LVL_VSFLX  ,
     . LVL_VMFLTD ,LVL_VTFLTD ,LVL_VTFLLD ,LVL_VSFLTD ,LVL_VSFLLD ,
     . LVL_VVEL   ,LVL_WFLX   ,LVL_WFLX2  ,LVL_TKE    ,LVL_GLS_PSI,
     . LVL_IDLAGE ,
     . MSC_MMFLXL ,MSC_MMFLXD ,MSC_MMFTDL ,MSC_MMFTDD ,MSC_MHFLX  ,
     . MSC_MHFTD  ,MSC_MHFLD  ,MSC_MSFLX  ,MSC_MSFTD  ,MSC_MSFLD  ,
     . MSC_VOLTR  ,
     . GLB_AVEPERIO,GLB_FILEFREQ,GLB_COMPFLAG,GLB_NCFORMAT, 
     . GLB_FNAMETAG



      contains 



      subroutine diafnm(runid,runid_len,ctag,nstep,diagfq,diagmon
     .  ,diagann,fname)
c
c --- ------------------------------------------------------------------
c --- Description: creates file name for the diagnostic output
c ---  
c --- Arguments: 
c ---   char runid     (in)    : acronym of simulation 
c ---   int  runid_len (in)    : actual length of runid 
c ---   char ctag      (in)    : string used in middle of file name 
c ---   real time      (in)    : MICOM time step 
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
      character runid*(*),ctag*(*),fname*(*)
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
#if   defined(CLIM)
        call pstdat(ny,nm,nd,nint(diagfq))
        if     (mod(nint(diagfq),360).eq.0) then
          write(fname,'(4a,i4.4,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.nc'
        elseif (mod(nint(diagfq),30).eq.0) then
          write(fname,'(4a,i4.4,a,i2.2,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.nc'
        else
          write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.',nd,'.nc'
        endif
#elif defined(SYN)
        if (diagmon) then
          call pstdat(ny,nm,nd,1)
          write(fname,'(4a,i4.4,a,i2.2,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.nc'
        elseif (diagann) then
          call pstdat(ny,nm,nd,1)
          write(fname,'(4a,i4.4,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.nc'
        else
          call pstdat(ny,nm,nd,nint(diagfq))
          write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)')
     .      runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.',nd,'.nc'
        endif
#elif defined(CCSMCOUPLED)
        if (diagmon) then
          call pstdat(ny,nm,nd,1)
          write(fname,'(4a,i4.4,a,i2.2,a)')
     .      runid(1:runid_len),'.micom.',ctag,'.',ny,'-',nm,'.nc'
        elseif (diagann) then
          call pstdat(ny,nm,nd,1)
          write(fname,'(4a,i4.4,a)')
     .      runid(1:runid_len),'.micom.',ctag,'.',ny,'.nc'
        else
          call pstdat(ny,nm,nd,nint(diagfq))
          write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a)')
     .      runid(1:runid_len),'.micom.',ctag,'.',ny,'-',nm,'-',nd,'.nc'
        endif
#endif
      else
        ny=nyear
        nm=nmonth
        nd=nday
        if (mod(nstep,nstep_in_day).eq.0) call pstdat(ny,nm,nd,1)
        ns=nint(nint(mod(nstep+.5-diagfq,real(nstep_in_day))/diagfq)
     .          *diagfq/nstep_in_day*86400)
#ifdef CCSMCOUPLED
        write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)')
     .    runid(1:runid_len),'.micom.',ctag,'.',ny,'-',nm,'-',nd,'-',ns,
     .    '.nc'
#else
        write(fname,'(4a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)')
     .    runid(1:runid_len),'_',ctag,'_',ny,'.',nm,'.',nd,'.',ns,'.nc'
#endif
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
      integer :: i,j,k,l,n,istat,istatsum,i1,j1,ncid,rhid
      integer, parameter :: imn=1-nbdy,imx=idm+nbdy,jmn=imn,jmx=jdm+nbdy
      logical fexist,lsigmar1
      real, dimension (itdm,jtdm) :: tmp2d
c
c --- Check existence of data files for meridional and section transport
c --- diagnostics 
      inquire(file=path(1:path_len)//'mertraoceans.dat',exist=fexist)
      if (.not.fexist) then
        write (lp,*) 'Could not find file ',
     .                path(1:path_len)//'mertraoceans.dat','!'
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
      inquire(file=path(1:path_len)//'mertraindex.dat',exist=fexist)
      if (.not.fexist) then
        write (lp,*) 'Could not find file ',
     .                path(1:path_len)//'mertraindex.dat','!'
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
      inquire(file=path(1:path_len)//'secindex.dat',exist=fexist)
      if (.not.fexist) then
        write (lp,*) 'Could not find file ',
     .                path(1:path_len)//'secindex.dat','!'
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
c
      nphysrf=0
      nphylyr=0
      nphylvl=0
c 
c --- Loop over io groups
      do n=1,nphy
        nacc_phy(n)=0
c
c --- - Solve dependencies for diagnostic variables (0=skipped) 
        ACC_ABSWND(n)   = SRF_ABSWND(n)
        ACC_ALB(n)      = SRF_ALB(n)
        ACC_BRNFLX(n)   = SRF_BRNFLX(n)
        ACC_BRNPD(n)    = SRF_BRNPD(n)
        ACC_DFL(n)      = SRF_DFL(n)
        ACC_EVA(n)      = SRF_EVA(n)
        ACC_FMLTFZ(n)   = SRF_FMLTFZ(n)
        ACC_FICE(n)     = SRF_FICE(n)   + SRF_HICE(n)   + SRF_UICE(n)  +
     .                    SRF_VICE(n)   + SRF_HSNW(n)
        ACC_HICE(n)     = SRF_HICE(n)   + SRF_UICE(n)   + SRF_VICE(n)
        ACC_HMLTFZ(n)   = SRF_HMLTFZ(n)
        ACC_HSNW(n)     = SRF_HSNW(n)
        ACC_IAGE(n)     = SRF_IAGE(n)
        ACC_IVOLU(n)    = SRF_UICE(n)
        ACC_IVOLV(n)    = SRF_VICE(n)
        ACC_LIP(n)      = SRF_LIP(n)
        ACC_MAXMLD(n)   = SRF_MAXMLD(n)
        ACC_MLD(n)      = SRF_MLD(n)
        ACC_MLDU(n)     = SRF_MLDU(n)   + SRF_MXLU(n)
        ACC_MLDV(n)     = SRF_MLDV(n)   + SRF_MXLV(n)
        ACC_MTY(n)      = SRF_MTY(n)
        ACC_MXLU(n)     = SRF_MXLU(n)
        ACC_MXLV(n)     = SRF_MXLV(n)
        ACC_NSF(n)      = SRF_NSF(n)
        ACC_RFIFLX(n)   = SRF_RFIFLX(n)
        ACC_RNFFLX(n)   = SRF_RNFFLX(n)
        ACC_SURFLX(n)   = SRF_SURFLX(n)
        ACC_SURRLX(n)   = SRF_SURRLX(n)
        ACC_SALFLX(n)   = SRF_SALFLX(n)
        ACC_SALRLX(n)   = SRF_SALRLX(n)
        ACC_SEALV(n)    = SRF_SEALV(n)
        ACC_SFL(n)      = SRF_SFL(n)
        ACC_SIGMX(n)    = SRF_SIGMX(n)  + MSC_MMFLXL(n) + MSC_MMFTDL(n)
        ACC_SOP(n)      = SRF_SOP(n)
        ACC_SSS(n)      = SRF_SSS(n)
        ACC_SST(n)      = SRF_SST(n)
        ACC_SWA(n)      = SRF_SWA(n)
        ACC_TAUX(n)     = SRF_TAUX(n)
        ACC_TAUY(n)     = SRF_TAUY(n)
        ACC_TICE(n)     = SRF_TICE(n)
        ACC_TSRF(n)     = SRF_TSRF(n)
        ACC_UB(n)       = SRF_UB(n)
        ACC_UICE(n)     = SRF_UICE(n)
        ACC_USTAR(n)    = SRF_USTAR(n)
        ACC_VB(n)       = SRF_VB(n)
        ACC_VICE(n)     = SRF_VICE(n)
        ACC_ZTX(n)      = SRF_ZTX(n)
        ACC_DIAFLX(n)   = LYR_DIAFLX(n)
        ACC_DIAFLXLVL(n)= LVL_DIAFLX(n)
        ACC_DIFDIA(n)   = LYR_DIFDIA(n)
        ACC_DIFDIALVL(n)= LVL_DIFDIA(n)
        ACC_DIFINT(n)   = LYR_DIFINT(n)
        ACC_DIFINTLVL(n)= LVL_DIFINT(n)
        ACC_DIFISO(n)   = LYR_DIFISO(n)
        ACC_DIFISOLVL(n)= LVL_DIFISO(n)
        ACC_DP(n)       = LYR_DP(n)     +
     .                    LYR_SALN(n)   + LYR_TEMP(n)   + LYR_TRC(n)   +
     .                    LYR_DIFDIA(n) + LYR_DIFINT(n) + LYR_DIFISO(n)+
     .                    LYR_TKE(n)    + LYR_GLS_PSI(n)+
     .                    LVL_SALN(n)   + LVL_TEMP(n)   + LVL_TRC(n)   +
     .                    LVL_DIFDIA(n) + LVL_DIFINT(n) + LVL_DIFISO(n)+
     .                    LVL_TKE(n)    + LVL_GLS_PSI(n)
        ACC_DPU(n)      = LYR_DPU(n)    + LYR_UVEL(n)
        ACC_DPV(n)      = LYR_DPV(n)    + LYR_VVEL(n)
        ACC_DZ(n)       = LYR_DZ(n)
        ACC_DZLVL(n)    = LVL_DZ(n)
        ACC_SALN(n)     = LYR_SALN(n)
        ACC_SALNLVL(n)  = LVL_SALN(n)
        ACC_TEMP(n)     = LYR_TEMP(n)
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
        if (ACC_ABSWND(n).ne.0) nphysrf=nphysrf+1 
         ACC_ABSWND(n)=nphysrf*min(1,ACC_ABSWND(n))
        if (ACC_ALB(n).ne.0) nphysrf=nphysrf+1
         ACC_ALB(n)=nphysrf*min(1,ACC_ALB(n)) 
        if (ACC_BRNFLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_BRNFLX(n)=nphysrf*min(1,ACC_BRNFLX(n))
        if (ACC_BRNPD(n).ne.0) nphysrf=nphysrf+1
         ACC_BRNPD(n)=nphysrf*min(1,ACC_BRNPD(n))
        if (ACC_DFL(n).ne.0) nphysrf=nphysrf+1
         ACC_DFL(n)=nphysrf*min(1,ACC_DFL(n)) 
        if (ACC_EVA(n).ne.0) nphysrf=nphysrf+1
         ACC_EVA(n)=nphysrf*min(1,ACC_EVA(n)) 
        if (ACC_FMLTFZ(n).ne.0) nphysrf=nphysrf+1 
         ACC_FMLTFZ(n)=nphysrf*min(1,ACC_FMLTFZ(n))
        if (ACC_FICE(n).ne.0) nphysrf=nphysrf+1
         ACC_FICE(n)=nphysrf*min(1,ACC_FICE(n)) 
        if (ACC_HICE(n).ne.0) nphysrf=nphysrf+1
         ACC_HICE(n)=nphysrf*min(1,ACC_HICE(n)) 
        if (ACC_HMLTFZ(n).ne.0) nphysrf=nphysrf+1 
         ACC_HMLTFZ(n)=nphysrf*min(1,ACC_HMLTFZ(n))
        if (ACC_HSNW(n).ne.0) nphysrf=nphysrf+1
         ACC_HSNW(n)=nphysrf*min(1,ACC_HSNW(n)) 
        if (ACC_IAGE(n).ne.0) nphysrf=nphysrf+1
         ACC_IAGE(n)=nphysrf*min(1,ACC_IAGE(n)) 
        if (ACC_IVOLU(n).ne.0) nphysrf=nphysrf+1 
         ACC_IVOLU(n)=nphysrf*min(1,ACC_IVOLU(n))
        if (ACC_IVOLV(n).ne.0) nphysrf=nphysrf+1 
         ACC_IVOLV(n)=nphysrf*min(1,ACC_IVOLV(n))
        if (ACC_LIP(n).ne.0) nphysrf=nphysrf+1
         ACC_LIP(n)=nphysrf*min(1,ACC_LIP(n)) 
        if (ACC_MAXMLD(n).ne.0) nphysrf=nphysrf+1 
         ACC_MAXMLD(n)=nphysrf*min(1,ACC_MAXMLD(n))
        if (ACC_MLD(n).ne.0) nphysrf=nphysrf+1
         ACC_MLD(n)=nphysrf*min(1,ACC_MLD(n)) 
        if (ACC_MLDU(n).ne.0) nphysrf=nphysrf+1
         ACC_MLDU(n)=nphysrf*min(1,ACC_MLDU(n)) 
        if (ACC_MLDV(n).ne.0) nphysrf=nphysrf+1
         ACC_MLDV(n)=nphysrf*min(1,ACC_MLDV(n))
        if (ACC_MTY(n).ne.0) nphysrf=nphysrf+1
         ACC_MTY(n)=nphysrf*min(1,ACC_MTY(n)) 
        if (ACC_MXLU(n).ne.0) nphysrf=nphysrf+1
         ACC_MXLU(n)=nphysrf*min(1,ACC_MXLU(n)) 
        if (ACC_MXLV(n).ne.0) nphysrf=nphysrf+1
         ACC_MXLV(n)=nphysrf*min(1,ACC_MXLV(n))
        if (ACC_NSF(n).ne.0) nphysrf=nphysrf+1
         ACC_NSF(n)=nphysrf*min(1,ACC_NSF(n)) 
        if (ACC_RFIFLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_RFIFLX(n)=nphysrf*min(1,ACC_RFIFLX(n))
        if (ACC_RNFFLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_RNFFLX(n)=nphysrf*min(1,ACC_RNFFLX(n))
        if (ACC_SURFLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_SURFLX(n)=nphysrf*min(1,ACC_SURFLX(n))
        if (ACC_SURRLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_SURRLX(n)=nphysrf*min(1,ACC_SURRLX(n))
        if (ACC_SALFLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_SALFLX(n)=nphysrf*min(1,ACC_SALFLX(n))
        if (ACC_SALRLX(n).ne.0) nphysrf=nphysrf+1 
         ACC_SALRLX(n)=nphysrf*min(1,ACC_SALRLX(n))
        if (ACC_SEALV(n).ne.0) nphysrf=nphysrf+1 
         ACC_SEALV(n)=nphysrf*min(1,ACC_SEALV(n))
        if (ACC_SFL(n).ne.0) nphysrf=nphysrf+1
         ACC_SFL(n)=nphysrf*min(1,ACC_SFL(n)) 
        if (ACC_SIGMX(n).ne.0) nphysrf=nphysrf+1 
         ACC_SIGMX(n)=nphysrf*min(1,ACC_SIGMX(n))
        if (ACC_SOP(n).ne.0) nphysrf=nphysrf+1
         ACC_SOP(n)=nphysrf*min(1,ACC_SOP(n)) 
        if (ACC_SSS(n).ne.0) nphysrf=nphysrf+1
         ACC_SSS(n)=nphysrf*min(1,ACC_SSS(n)) 
        if (ACC_SST(n).ne.0) nphysrf=nphysrf+1
         ACC_SST(n)=nphysrf*min(1,ACC_SST(n)) 
        if (ACC_SWA(n).ne.0) nphysrf=nphysrf+1
         ACC_SWA(n)=nphysrf*min(1,ACC_SWA(n)) 
        if (ACC_TAUX(n).ne.0) nphysrf=nphysrf+1
         ACC_TAUX(n)=nphysrf*min(1,ACC_TAUX(n)) 
        if (ACC_TAUY(n).ne.0) nphysrf=nphysrf+1
         ACC_TAUY(n)=nphysrf*min(1,ACC_TAUY(n)) 
        if (ACC_TICE(n).ne.0) nphysrf=nphysrf+1
         ACC_TICE(n)=nphysrf*min(1,ACC_TICE(n)) 
        if (ACC_TSRF(n).ne.0) nphysrf=nphysrf+1
         ACC_TSRF(n)=nphysrf*min(1,ACC_TSRF(n))
        if (ACC_UB(n).ne.0) nphysrf=nphysrf+1
         ACC_UB(n)=nphysrf*min(1,ACC_UB(n)) 
        if (ACC_UICE(n).ne.0) nphysrf=nphysrf+1
         ACC_UICE(n)=nphysrf*min(1,ACC_UICE(n)) 
        if (ACC_USTAR(n).ne.0) nphysrf=nphysrf+1 
         ACC_USTAR(n)=nphysrf*min(1,ACC_USTAR(n))
        if (ACC_VB(n).ne.0) nphysrf=nphysrf+1
         ACC_VB(n)=nphysrf*min(1,ACC_VB(n)) 
        if (ACC_VICE(n).ne.0) nphysrf=nphysrf+1
         ACC_VICE(n)=nphysrf*min(1,ACC_VICE(n))
        if (ACC_ZTX(n).ne.0) nphysrf=nphysrf+1
         ACC_ZTX(n)=nphysrf*min(1,ACC_ZTX(n))
c
        if (ACC_DIAFLX(n).ne.0) nphylyr=nphylyr+1
         ACC_DIAFLX(n)=nphylyr*min(1,ACC_DIAFLX(n))
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
        if (ACC_TKE(n).ne.0) nphylyr=nphylyr+1
         ACC_TKE(n)=nphylyr*min(1,ACC_TKE(n))
        if (ACC_GLS_PSI(n).ne.0) nphylyr=nphylyr+1
         ACC_GLS_PSI(n)=nphylyr*min(1,ACC_GLS_PSI(n))
c
        if (ACC_DIAFLXLVL(n).ne.0) nphylvl=nphylvl+1
         ACC_DIAFLXLVL(n)=nphylvl*min(1,ACC_DIAFLXLVL(n))
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
        if (ACC_TKELVL(n).ne.0) nphylvl=nphylvl+1
         ACC_TKELVL(n)=nphylvl*min(1,ACC_TKELVL(n))
        if (ACC_GLS_PSILVL(n).ne.0) nphylvl=nphylvl+1
         ACC_GLS_PSILVL(n)=nphylvl*min(1,ACC_GLS_PSILVL(n))
c
c --- End loop over io groups
      enddo
c 
c --- Assign buffer positions for utility fields
      ACC_UTIL=0 
      nphylyr=nphylyr+1
      ACC_UTIL(1)=nphylyr 
c
      ACC_UTILLVL=0 
      nphylvl=nphylvl+1
      ACC_UTILLVL(1)=nphylvl 
c 
c --- Allocate buffers
      istatsum=0
      istat=0
      if (nphysrf.ne.0)
     .  allocate(physrf(imn:imx,jmn:jmx,nphysrf),stat=istat)
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
c --- Set buffers to zero (ocean+land+halo) 
c$OMP PARALLEL DO
      do j=jmn,jmx
        do n=1,nphysrf
          do i=imn,imx
            physrf(i,j,n)=0.
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=jmn,jmx
        do n=1,nphylyr
          do k=1,kdm
            do i=imn,imx
              phylyr(i,j,k,n)=0.
            enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=jmn,jmx
        do n=1,nphylvl
          do k=1,ddm
            do i=imn,imx
              phylvl(i,j,k,n)=0.
            enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=jmn,jmx
        do k=1,kdm
          do i=imn,imx
            diaflx(i,j,k)=0.
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Load bathymetry into module mod_dia (used for vertical 
c --- interpolation in MICOM and HAMOCC)
      nstepinday=nstep_in_day
c$OMP PARALLEL DO  
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
c --- Extract reference density and convert units g cm-3 -> kg m-3
      sigmar1=-99.
      if (mnproc.eq.1) then
        call ncerro(nf90_open(path(1:path_len)//'inicon.nc',
     .    nf90_nowrite,ncid))
        call ncerro(nf90_inq_varid(ncid,'sigma',rhid))
c --- - Read surface density and find first ocean grid cell
        call ncerro(nf90_get_var(ncid,rhid,tmp2d,start=(/1,1,1/),
     .    count=(/itdm,jtdm,1/)))
        lsigmar1=.false.
        do j=1,jtdm
          do i=1,itdm
            if (tmp2d(i,j).gt.0) then
              i1=i
              j1=j 
              lsigmar1=.true.
              exit
            endif
          enddo
          if (lsigmar1) exit
        enddo
c --- - Read density vector and close ini-file
        call ncerro(nf90_get_var(ncid,rhid,sigmar1,start=(/i1,j1,1/),
     .    count=(/1,1,kk/)))
        call ncerro(nf90_close(ncid))
        write(lp,*) 'sigma layers=',sigmar1 
      endif
c
      end subroutine diaini



      subroutine diaacc(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- accumulate diagnostic variables
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_ice.h"
#if defined TKE
#  include "common_tke.h"
#endif
c
      integer d,i,j,k,l,km,iogrp
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
c     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts,
     .  wghtsflx
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: uvel,vvel,
     .  dz,dummy
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: util5,util6,
     .  util7,util8,util9
c
c --- Increase counter 
      do iogrp=1,nphy
        nacc_phy(iogrp)=nacc_phy(iogrp)+1
      enddo
c
c --- Define auxillary variables 
      if (sum(ACC_UICE(1:nphy)+ACC_VICE(1:nphy)).ne.0) then
        call xctilr(hicem, 1,1, 1,1, halo_ps)
        call xctilr(ficem, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO
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
      if (sum(ACC_UB(1:nphy)+ACC_VB(1:nphy)).ne.0) then
c$OMP PARALLEL DO
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
      if (sum(ACC_SURFLX(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util5(i,j)=-surflx(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_SALFLX(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util6(i,j)=-salflx(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_SURRLX(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util7(i,j)=-surrlx(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_SALRLX(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util8(i,j)=-salrlx(i,j)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      if (sum(ACC_MLD(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util9(i,j)=dp(i,j,1+mm)+dp(i,j,2+mm)
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- maximum mixed layer pressure thickness [g/cm/s^2]
      do iogrp=1,nphy
        if (ACC_MAXMLD(iogrp).ne.0) then
c$OMP PARALLEL DO
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                physrf(i,j,ACC_MAXMLD(iogrp))=
     .            max(physrf(i,j,ACC_MAXMLD(iogrp)),
     .                dp(i,j,1+mm)+dp(i,j,2+mm))
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
      enddo
c
      if (sum(ACC_UVEL(1:nphy)+ACC_UVELLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kk
            km=k+mm
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                uvel(i,j,k)=u(i,j,km)+ub(i,j,m)
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      if (sum(ACC_VVEL(1:nphy)+ACC_VVELLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
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
      if (sum(ACC_DZ(1:nphy)+ACC_DZLVL(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                dz(i,j,k)=phi(i,j,k)-phi(i,j,k+1)
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- ------------------------------------------------------------------
c --- accumulate 2d diagnostic variables
c --- ------------------------------------------------------------------
c --- u-component of barotropic velocity [cm/s]
      call accsrf(ACC_UB,ub(1-nbdy,1-nbdy,m),dummy,0,'u')
c
c --- u-component of wind stress [N/m^2]
      call accsrf(ACC_ZTX,ztx,dummy,0,'u')
c
c --- u-component of momentum flux received by the ocean [dyn/cm^2]
      call accsrf(ACC_TAUX,taux,dummy,0,'u')
c
c --- weighted u-component of ice velocity [m^2/s]
      call accsrf(ACC_UICE,uicem,util1,1,'u')
c
c --- weighted u-component of total velocity [g/s^3]
      call accsrf(ACC_MXLU,util2,dummy,0,'u')
c
c --- mixed layer pressure thickness at u-point [g/cm/s^2]
      call accsrf(ACC_MLDU,dpu(1-nbdy,1-nbdy,k1m),dummy,0,'u')
c
c --- v-component of barotropic velocity [cm/s]
      call accsrf(ACC_VB,vb(1-nbdy,1-nbdy,m),dummy,0,'v')
c
c --- v-component of wind stress [N/m^2]
      call accsrf(ACC_MTY,mty,dummy,0,'v')
c
c --- v-component of momentum flux received by the ocean [dyn/cm^2]
      call accsrf(ACC_TAUY,tauy,dummy,0,'v')
c
c --- weighted v-component of ice velocity [m^2/s]
      call accsrf(ACC_VICE,vicem,util3,1,'v')
c
c --- weighted v-component of total velocity [g/s^3]
      call accsrf(ACC_MXLV,util4,dummy,0,'v')
c
c --- mixed layer pressure thickness at v-point [g/cm/s^2]
      call accsrf(ACC_MLDV,dpv(1-nbdy,1-nbdy,k1m),dummy,0,'v')
c
c --- sea level height [cm]
      call accsrf(ACC_SEALV,sealv,dummy,0,'p')
c
c --- mixed layer density (sigma units)
      call accsrf(ACC_SIGMX,sigma(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- weighted ice thickness [m]
      call accsrf(ACC_HICE,hicem,ficem,1,'p')
c
c --- weighted snow thickness [m]
      call accsrf(ACC_HSNW,hsnwm,ficem,1,'p')
c
c --- fractional ice cover
      call accsrf(ACC_FICE,ficem,dummy,0,'p')
c 
c --- ice volume in u-points[m] 
      call accsrf(ACC_IVOLU,util1,dummy,0,'u')
c 
c --- ice volume in v-points[m] 
      call accsrf(ACC_IVOLV,util3,dummy,0,'v')
c
c --- surface temperature [K]
      call accsrf(ACC_TSRF,tsrfm,dummy,0,'p')
c
c --- ice temperature [K]
      call accsrf(ACC_TICE,ticem,dummy,0,'p')
c
c --- short wave heat flux [W/m^2]
      call accsrf(ACC_SWA,swa,dummy,0,'p')
c
c --- non-solar heat flux [W/m^2]
      call accsrf(ACC_NSF,nsf,dummy,0,'p')
c
c --- heat flux due to melting/freezing [W/m^2]
      call accsrf(ACC_HMLTFZ,hmltfz,dummy,0,'p')
c
c --- derivative of non-solar heat flux by surface temperature [W/m^2/K]
      call accsrf(ACC_DFL,dfl,dummy,0,'p')
c
c --- liquid precipitation [mm/day]
      call accsrf(ACC_LIP,lip,dummy,0,'p')
c
c --- solid precipitation [mm/day]
      call accsrf(ACC_SOP,sop,dummy,0,'p')
c
c --- evaporation [mm/day]
      call accsrf(ACC_EVA,eva,dummy,0,'p')
c
c --- fresh water flux due to melting/freezing [kg/m^2/s]
      call accsrf(ACC_FMLTFZ,fmltfz,dummy,0,'p')
c
c --- salt flux [kg/m^2/s]
      call accsrf(ACC_SFL,sfl,dummy,0,'p')
c
c --- albedo
      call accsrf(ACC_ALB,alb,dummy,0,'p')
c
c --- liquid runoff [kg m-2 s-1]
      call accsrf(ACC_RNFFLX,rnfflx,dummy,0,'p')
c
c --- frozen runoff [kg m-2 s-1]
      call accsrf(ACC_RFIFLX,rfiflx,dummy,0,'p')
c
c --- friction velocity [cm s-1]
      call accsrf(ACC_USTAR,ustar,dummy,0,'p')
c
c --- absolute wind speed [m s-1]
      call accsrf(ACC_ABSWND,abswnd,dummy,0,'p')
c
c --- sea surface salinity [g kg-1]
      call accsrf(ACC_SSS,saln(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- sea surface temperature [degC]
      call accsrf(ACC_SST,temp(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- mixed layer pressure thickness [g/cm/s^2]
      call accsrf(ACC_MLD,util9,dummy,0,'p')
c
c --- heat flux received by the ocean [W/m^2]
      call accsrf(ACC_SURFLX,util5,dummy,0,'p')
c
c --- salt flux received by the ocean [10^-3 g/cm^2/s]
      call accsrf(ACC_SALFLX,util6,dummy,0,'p')
c
c --- restoring heat flux received by the ocean [W/m^2]
      call accsrf(ACC_SURRLX,util7,dummy,0,'p')
c
c --- restoring salt flux received by the ocean [10^-3 g/cm^2/s]
      call accsrf(ACC_SALRLX,util8,dummy,0,'p')
c
c --- brine flux received by the ocean [10^-3 g/cm^2/s]
      call accsrf(ACC_BRNFLX,brnflx,dummy,0,'p')
c
c --- brine plume pressure depth [g/cm/s^2]
      call accsrf(ACC_BRNPD,pbrnda,dummy,0,'p')
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
     .  dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- weighted temperature [degC*g/cm/s^2]
      call acclyr(ACC_TEMP,temp(1-nbdy,1-nbdy,k1m),
     .  dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- layer pressure thickness [g/cm/s^2]
      call acclyr(ACC_DP,dp(1-nbdy,1-nbdy,k1m),dummy,0,'p')
c
c --- layer geopotential thickness [cm^2/s^2]
      call acclyr(ACC_DZ,dz,dummy,0,'p')
c
c --- diapycnal mass flux [g*cm/s^2]
      call acclyr(ACC_DIAFLX,diaflx,dummy,0,'p')
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
#if defined TKE
c --- weighted tke [cm2/s2*g/cm/s^2]
      call acclyr(ACC_TKE,tke(1-nbdy,1-nbdy,k1m),
     .  dp(1-nbdy,1-nbdy,k1m),1,'p')
c
c --- weighted gls_psi [cm2/s3*g/cm/s^2]
      call acclyr(ACC_GLS_PSI,gls_psi(1-nbdy,1-nbdy,k1m),
     .  dp(1-nbdy,1-nbdy,k1m),1,'p')
c
#endif
c --- ------------------------------------------------------------------
c --- accumulate 3d diagnostic variables on Levitus levels 
c --- ------------------------------------------------------------------
c
      do iogrp=1,nphy 
        if (ACC_WFLXLVL(iogrp)+ACC_WFLX2LVL(iogrp).ne.0) then
c$OMP PARALLEL DO
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
          call acclvl(ACC_UFLXLVL,uflx(1-nbdy,1-nbdy,k1n),'u',k,ind1,
     .      ind2,wghtsflx)
c
c --- --- u-component of heat flux [K*g*cm/s^2]
          call acclvl(ACC_UTFLXLVL,utflx(1-nbdy,1-nbdy,k1n),'u',k,ind1,
     .      ind2,wghtsflx)
c
c --- --- u-component of salt flux [g^2*cm/kg/s^2]
          call acclvl(ACC_USFLXLVL,usflx(1-nbdy,1-nbdy,k1n),'u',k,ind1,
     .      ind2,wghtsflx)
c
c --- --- u-component of mass flux due to thickness diffusion [g*cm/s^2]
          call acclvl(ACC_UMFLTDLVL,umfltd(1-nbdy,1-nbdy,k1n),'u',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- u-component of heat flux due to thickness diffusion [K*g*cm/s^2]
          call acclvl(ACC_UTFLTDLVL,utfltd(1-nbdy,1-nbdy,k1n),'u',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- u-component of heat flux due to lateral diffusion [K*g*cm/s^2]
          call acclvl(ACC_UTFLLDLVL,utflld(1-nbdy,1-nbdy,k1n),'u',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- u-component of salt flux due to thickness diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_USFLTDLVL,usfltd(1-nbdy,1-nbdy,k1n),'u',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- u-component of salt flux due to lateral diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_USFLLDLVL,usflld(1-nbdy,1-nbdy,k1n),'u',k,
     .      ind1,ind2,wghtsflx)
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
          call acclvl(ACC_VFLXLVL,vflx(1-nbdy,1-nbdy,k1n),'v',k,ind1,
     .      ind2,wghtsflx)
c
c --- --- v-component of heat flux [K*g*cm/s^2]
          call acclvl(ACC_VTFLXLVL,vtflx(1-nbdy,1-nbdy,k1n),'v',k,ind1,
     .      ind2,wghtsflx)
c
c --- --- v-component of salt flux [g^2*cm/kg/s^2]
          call acclvl(ACC_VSFLXLVL,vsflx(1-nbdy,1-nbdy,k1n),'v',k,ind1,
     .      ind2,wghtsflx)
c
c --- --- v-component of mass flux due to thickness diffusion [g*cm/s^2]
          call acclvl(ACC_VMFLTDLVL,vmfltd(1-nbdy,1-nbdy,k1n),'v',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- v-component of heat flux due to thickness diffusion [K*g*cm/s^2]
          call acclvl(ACC_VTFLTDLVL,vtfltd(1-nbdy,1-nbdy,k1n),'v',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- v-component of heat flux due to lateral diffusion [K*g*cm/s^2]
          call acclvl(ACC_VTFLLDLVL,vtflld(1-nbdy,1-nbdy,k1n),'v',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- v-component of salt flux due to thickness diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_VSFLTDLVL,vsfltd(1-nbdy,1-nbdy,k1n),'v',k,
     .      ind1,ind2,wghtsflx)
c
c --- --- v-component of salt flux due to lateral diffusion [g^2*cm/kg/s^2]
          call acclvl(ACC_VSFLLDLVL,vsflld(1-nbdy,1-nbdy,k1n),'v',k,
     .      ind1,ind2,wghtsflx)
        enddo
      endif
c
      if (sum(ACC_SALNLVL(1:nphy)+ACC_TEMPLVL(1:nphy)+
     .  ACC_DIAFLXLVL(1:nphy)+ACC_DIFDIALVL(1:nphy)+
     .  ACC_DIFINTLVL(1:nphy)+ACC_DIFISOLVL(1:nphy)+
     .  ACC_TKELVL(1:nphy)+ACC_GLS_PSILVL(1:nphy)+
     .  ACC_DZLVL(1:nphy)).ne.0) then
        do k=1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
c
c --- --- salinity [g^2/kg/cm/s^2]
          call acclvl(ACC_SALNLVL,saln(1-nbdy,1-nbdy,k1m),'p',k,ind1
     .      ,ind2,wghts)
c
c --- --- temperature [degC*g/cm/s^2]
          call acclvl(ACC_TEMPLVL,temp(1-nbdy,1-nbdy,k1m),'p',k,ind1
     .      ,ind2,wghts)
c
c --- --- diapycnal mass flux [g*cm/s^2]
          call acclvl(ACC_DIAFLXLVL,diaflx,'p',k,ind1,ind2,wghts)
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
#if defined TKE
c --- --- tke [cm2/s2*g/cm/s^2]
          call acclvl(ACC_TKELVL,tke(1-nbdy,1-nbdy,k1m),'p',k,ind1
     .      ,ind2,wghts)
c
c --- --- gls_psi [cm2/s3*g/cm/s^2]
          call acclvl(ACC_GLS_PSILVL,gls_psi(1-nbdy,1-nbdy,k1m),'p',
     .       k,ind1,ind2,wghts)
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
c --- Reset diaflx 
      if (sum(ACC_DIAFLX(1:nphy)).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kdm
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                diaflx(i,j,k)=0.
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
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
#include "common_ice.h"
c
      integer iogrp
c
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iuu,ivv
      logical iniflg,append2file(nphymax)
      data iniflg /.true./ , append2file /nphymax*.false./
      save iniflg,iuu,ivv,append2file
c
      integer i,j,k,l,cmpflg,irec(nphymax)
      character startdate*20,timeunits*30,fname(nphymax)*80
      real datenum,rnacc
      save fname,irec
c
#if defined(TRC) && defined(IDLAGE)
#  include "param_trc.h"
#  include "common_trc.h"
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts,
     .   wghtsflx
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
c --- finalize accumulation of 2d fields
      call finsrf(ACC_HICE(iogrp),ACC_FICE(iogrp),'p')
      call finsrf(ACC_HSNW(iogrp),ACC_FICE(iogrp),'p')
      call finsrf(ACC_MXLU(iogrp),ACC_MLDU(iogrp),'u')
      call finsrf(ACC_MXLV(iogrp),ACC_MLDV(iogrp),'v')
      call finsrf(ACC_UICE(iogrp),ACC_IVOLU(iogrp),'u')
      call finsrf(ACC_VICE(iogrp),ACC_IVOLV(iogrp),'v')
c
c --- finalize accumulation of layer fields
      call finlyr(ACC_UVEL(iogrp),ACC_DPU(iogrp),'u')
      call finlyr(ACC_VVEL(iogrp),ACC_DPV(iogrp),'v')
      call finlyr(ACC_SALN(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_TEMP(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_DIFDIA(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_DIFINT(iogrp),ACC_DP(iogrp),'p')
      call finlyr(ACC_DIFISO(iogrp),ACC_DP(iogrp),'p')
#if defined TKE
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
      call msklvl(ACC_DIAFLXLVL(iogrp),'p')
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
#if defined TKE
      call msklvl(ACC_TKELVL(iogrp),'p')
      call msklvl(ACC_GLS_PSILVL(iogrp),'p')
#endif
c 
c --- get instantaneous values for ice age 
      if (ACC_IAGE(iogrp).ne.0) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              physrf(i,j,ACC_IAGE(iogrp))=iagem(i,j)
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
        call diafnm(runid,runid_len,trim(GLB_FNAMETAG(IOGRP)),nstep,
     .  filefq_phy(iogrp)/real(nstep_in_day),filemon_phy(iogrp),
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
     .   write (lp,'(2a)')
     .    'Writing physical diagnostics to file: ',
     .    path1(1:path1_len)//fname(iogrp)
      if (GLB_NCFORMAT(iogrp).eq.1) then
        call ncfopn(path1(1:path1_len)//fname(iogrp),'w','6',
     .    irec(iogrp))
      elseif (GLB_NCFORMAT(iogrp).eq.2) then
        call ncfopn(path1(1:path1_len)//fname(iogrp),'w','h',
     .    irec(iogrp))
      else
        call ncfopn(path1(1:path1_len)//fname(iogrp),'w','c',
     .    irec(iogrp))
      endif
c
c --- compute extended ocean masks 
      if (iniflg) then
         iniflg=.false.
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
      call ncdims('region',odm)
      call ncdims('bounds',2) 
      call ncdims('slenmax',slenmax)
c
      if (ACC_VOLTR(iogrp).ne.0) then
        if ((sec_num.gt.0.and.sec_num.le.max_sec).or.mnproc.ne.1) then
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
        if ((lmax.gt.0.and.lmax.le.ldm).or.mnproc.gt.1) then
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
          region_name(1)='atlantic_arctic_ocean'
          region_name(2)='indian_pacific_ocean'
          region_name(3)='global_ocean'
          call ncwrtc('region','slenmax region',region_name)  
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
      call wrtsrf(ACC_SIGMX(iogrp),SRF_SIGMX(iogrp),rnacc*1e3,
     ,  0.,cmpflg,ip,'p','sigmx','Mixed layer density',' ','kg m-3')
c
      call wrtsrf(ACC_UB(iogrp),SRF_UB(iogrp),rnacc*1e-2,
     .  0.,cmpflg,iuu,'u','ubaro','Barotropic velocity x-component',
     .  ' ','m s-1')
c
      call wrtsrf(ACC_VB(iogrp),SRF_VB(iogrp),rnacc*1e-2,
     .  0.,cmpflg,ivv,'v','vbaro','Barotropic velocity y-component',
     .  ' ','m s-1')
c
      call wrtsrf(ACC_SEALV(iogrp),SRF_SEALV(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','sealv','Sea level',' ','m')
c
      call wrtsrf(ACC_HICE(iogrp),SRF_HICE(iogrp),1.,0.,
     .  cmpflg,ip,'p','hice','Ice thickness',' ','m')
c
      call wrtsrf(ACC_TICE(iogrp),SRF_TICE(iogrp),rnacc,
     .  -t0deg,cmpflg,ip,'p','tice','Ice temperature',' ','degC')
c
      call wrtsrf(ACC_HSNW(iogrp),SRF_HSNW(iogrp),1.,0.,
     .  cmpflg,ip,'p','hsnw','Snow depth',' ','m')
c
      call wrtsrf(ACC_FICE(iogrp),SRF_FICE(iogrp),
     .  rnacc*1e2,0.,cmpflg,ip,'p','fice','Ice concentration',' ','%')
c
      call wrtsrf(ACC_TSRF(iogrp),SRF_TSRF(iogrp),rnacc,
     .  -t0deg,cmpflg,ip,'p','tsrf','Surface temperature',' ','degC')
c
      call wrtsrf(ACC_IAGE(iogrp),SRF_IAGE(iogrp),1.,0.,
     .  cmpflg,ip,'p','iage','Ice age',' ','day')
c
      call wrtsrf(ACC_UICE(iogrp),SRF_UICE(iogrp),1e-2,0.,
     .  cmpflg,iuu,'u','uice','Ice velocity x-component',' ','m s-1')
c
      call wrtsrf(ACC_VICE(iogrp),SRF_VICE(iogrp),1e-2,0.,
     .  cmpflg,ivv,'v','vice','Ice velocity y-component',' ','m s-1')
c
      call wrtsrf(ACC_SWA(iogrp),SRF_SWA(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','swa','Short-wave heat flux',' ','W m-2')
c
      call wrtsrf(ACC_NSF(iogrp),SRF_NSF(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','nsf','Non-solar heat flux',' ','W m-2')
c
      call wrtsrf(ACC_HMLTFZ(iogrp),SRF_HMLTFZ(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','hmltfz',
     .  'Heat flux due to melting/freezing',' ','W m-2')
c
      call wrtsrf(ACC_DFL(iogrp),SRF_DFL(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','dfl','Non-solar heat flux derivative',' ',
     .  'W m-2 K-1')
c
      call wrtsrf(ACC_SURFLX(iogrp),SRF_SURFLX(iogrp),
     .  rnacc*1e4,0.,cmpflg,ip,'p','hflx',
     .  'Heat flux received by ocean',' ','W m-2')
c
      call wrtsrf(ACC_SURRLX(iogrp),SRF_SURRLX(iogrp),
     .  rnacc*1e4,0.,cmpflg,ip,'p','hrflx',
     .  'Restoring heat flux received by ocean',' ','W m-2')
c
      call wrtsrf(ACC_LIP(iogrp),SRF_LIP(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','lip','Liquid precipitation',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_SOP(iogrp),SRF_SOP(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sop','Solid precipitation',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_EVA(iogrp),SRF_EVA(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','eva','Evaporation',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_FMLTFZ(iogrp),SRF_FMLTFZ(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','fmltfz',
     .  'Fresh water flux due to melting/freezing',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_RNFFLX(iogrp),SRF_RNFFLX(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','rnf','Liquid runoff',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_RFIFLX(iogrp),SRF_RFIFLX(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','rfi','Frozen runoff',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_SALFLX(iogrp),SRF_SALFLX(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','sflx',
     .  'Salt flux received by ocean',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_SALRLX(iogrp),SRF_SALRLX(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','srflx',
     .  'Restoring salt flux received by ocean',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_BRNFLX(iogrp),SRF_BRNFLX(iogrp),
     .  rnacc*(-1e-2),0.,cmpflg,ip,'p','bflx','Brine flux',' ',
     .  'kg m-2 s-1')
c
      call wrtsrf(ACC_ZTX(iogrp),SRF_ZTX(iogrp),rnacc,0.,
     .  cmpflg,iuu,'u','ztx','Wind stress x-component',' ','N m-2')
c
      call wrtsrf(ACC_MTY(iogrp),SRF_MTY(iogrp),rnacc,0.,
     .  cmpflg,ivv,'v','mty','Wind stress y-component',' ','N m-2')
c
      call wrtsrf(ACC_TAUX(iogrp),SRF_TAUX(iogrp),rnacc*.1,
     .  0.,cmpflg,iuu,'u','taux',
     .  'Momentum flux received by ocean x-component',' ','N m-2')
c
      call wrtsrf(ACC_TAUY(iogrp),SRF_TAUY(iogrp),rnacc*.1,
     .  0.,cmpflg,ivv,'v','tauy', 
     .  'Momentum flux received by ocean y-component',' ','N m-2')
c
      call wrtsrf(ACC_USTAR(iogrp),SRF_USTAR(iogrp),
     .  rnacc*1e-2,0.,cmpflg,ip,'p','ustar','Friction velocity',' ',
     .  'm s-1')
c
      call wrtsrf(ACC_ABSWND(iogrp),SRF_ABSWND(iogrp),
     .  rnacc,0.,cmpflg,ip,'p','abswnd','Absolute wind speed',' ',
     .  'm s-1')
c
      call wrtsrf(ACC_SFL(iogrp),SRF_SFL(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sfl','Salt flux',' ','kg m-2 s-1')
c
      call wrtsrf(ACC_ALB(iogrp),SRF_ALB(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','alb','Surface albedo',' ','1')
c
      call wrtsrf(ACC_MLD(iogrp),SRF_MLD(iogrp),rnacc/onem,
     .  0.,cmpflg,ip,'p','mld','Mixed layer depth',' ','m')
c
      call wrtsrf(ACC_MLDU(iogrp),SRF_MLDU(iogrp),
     .  rnacc/onem,0.,cmpflg,iuu,'u','mldu',
     .  'Mixed layer depth at u-point',' ','m')
c
      call wrtsrf(ACC_MLDV(iogrp),SRF_MLDV(iogrp),
     .  rnacc/onem,0.,cmpflg,ivv,'v','mldv',
     .  'Mixed layer depth at v-point',' ','m')
c
      call wrtsrf(ACC_MAXMLD(iogrp),SRF_MAXMLD(iogrp),
     .  1./onem,0.,cmpflg,ip,'p','maxmld','Maximum mixed layer depth',
     .  ' ','m')
c
      call wrtsrf(ACC_BRNPD(iogrp),SRF_BRNPD(iogrp),rnacc/onem,
     .  0.,cmpflg,ip,'p','brnpd','Brine plume depth',' ','m')
c
      call wrtsrf(ACC_SST(iogrp),SRF_SST(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sst','Ocean surface temperature',' ','degC')
c
      call wrtsrf(ACC_SSS(iogrp),SRF_SSS(iogrp),rnacc,0.,
     .  cmpflg,ip,'p','sss','Ocean surface salinity',' ','g kg-1')
c
      call wrtsrf(ACC_MXLU(iogrp),SRF_MXLU(iogrp),1e-2,0.,
     .  cmpflg,iuu,'u','mxlu','Mixed layer velocity x-component',' ',
     .  'm s-1')
c
      call wrtsrf(ACC_MXLV(iogrp),SRF_MXLV(iogrp),1e-2,0.,
     .  cmpflg,ivv,'v','mxlv','Mixed layer velocity y-component',' ',
     .  'm s-1')
c
c --- write 3d layer fields
      call wrtlyr(ACC_DP(iogrp),LYR_DP(iogrp),rnacc*.1,0.
     .  ,cmpflg,ip,'p','dp','Layer pressure thickness',' ','Pa')
c
      call wrtlyr(ACC_DZ(iogrp),LYR_DZ(iogrp),
     .  rnacc*1e-2/g,0.,cmpflg,ip,'p','dz','Layer thickness',' ',
     .  'm')
c
      call wrtlyr(ACC_TEMP(iogrp),LYR_TEMP(iogrp),1.,0.,
     .  cmpflg,ip,'p','temp','Temperature','Ocean temperature',
     .  'degC')
c
      call wrtlyr(ACC_SALN(iogrp),LYR_SALN(iogrp),1.,0.,
     .  cmpflg,ip,'p','saln','Salinity','Ocean salinity',
     .  'g kg-1')
c
      call wrtlyr(ACC_UVEL(iogrp),LYR_UVEL(iogrp),1e-2,
     .  0.,cmpflg,iuu,'u','uvel','Velocity x-component',' ',
     .  'm s-1')
c
      call wrtlyr(ACC_VVEL(iogrp),LYR_VVEL(iogrp),1e-2,
     .  0.,cmpflg,ivv,'v','vvel','Velocity y-component',' ',
     .  'm s-1')
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
     .  'W')
c
      call wrtlyr(ACC_VSFLTD(iogrp),LYR_VSFLTD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsfltd',
     .  'Salt flux due to thickness diffusion in y-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_USFLLD(iogrp),LYR_USFLLD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,iuu,'u','usflld',
     .  'Salt flux due to lateral diffusion in x-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_VSFLLD(iogrp),LYR_VSFLLD(iogrp),
     .  rnacc*0.5e-6/(g*baclin),0.,cmpflg,ivv,'v','vsflld',
     .  'Salt flux due to lateral diffusion in y-direction',' ',
     .  'W')
c
      call wrtlyr(ACC_WFLX(iogrp),LYR_WFLX(iogrp),
     .  rnacc*0.5e-3/(g*baclin),0.,cmpflg,ip,'p','wflx',
     .  'Vertical mass flux',' ','kg s-1')
c
      call wrtlyr(ACC_WFLX2(iogrp),LYR_WFLX2(iogrp),
     .  rnacc*(0.5e-3/(g*baclin))**2,0.,cmpflg,ip,'p','wflx2',
     .  'Vertical mass flux squared',' ','kg2 s-2')
c
      call wrtlyr(ACC_DIAFLX(iogrp),LYR_DIAFLX(iogrp),
     .  rnacc/(onem*baclin),0.,cmpflg,ip,'p','diaflx',
     .  'Diapycnal volume flux',' ','m s-1')
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
#if defined TKE 
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
     .  rnacc*1e-2/g,0.,cmpflg,ip,'p','dzlvl','Layer thickness',' ',
     .  'm')
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
      call wrtlvl(ACC_DIAFLXLVL(iogrp),LVL_DIAFLX(iogrp),
     .  rnacc/(onem*baclin),0.,cmpflg,ip,'p','diaflxlvl',
     .  'Diapycnal volume flux',' ','m s-1')
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
#if defined TKE 
      call wrtlvl(ACC_TKELVL(iogrp),LVL_TKE(iogrp),
     .  1e-4,0.,cmpflg,ip,'p','tkelvl','TKE',
     .  'Turbulent Kinetic Energy','m2 s-2')
c
      call wrtlvl(ACC_GLS_PSILVL(iogrp),LVL_GLS_PSI(iogrp),
     .  1e-4,0.,cmpflg,ip,'p','gls_psilvl','GLS_PSI',
     .  'Generic length scale','m2 s-3')
c
#endif
c
c --- store meridional transports 
      if (MSC_MMFLXL(iogrp).ne.0) then 
        call ncwrt1('mmflxl','lat sigma region time',
     .    mmflxl(1:lmax,1:kdm,1:odm))
        call ncattr('long_name',
     .    'Overturning stream-function on isopycnic layers') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MMFLXD(iogrp).ne.0) then 
        call ncwrt1('mmflxd','lat depth region time',
     .    mmflxd(1:lmax,1:ddm,1:odm))
        call ncattr('long_name',
     .    'Overturning stream-function on z-levels') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MMFTDL(iogrp).ne.0) then 
        call ncwrt1('mmftdl','lat sigma region time',
     .    mmftdl(1:lmax,1:kdm,1:odm))
        call ncattr('long_name',
     .    'Overturning stream-function due to thickness diffusion '//
     .    'on isopycnic layers') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MMFTDD(iogrp).ne.0) then 
        call ncwrt1('mmftdd','lat depth region time',
     .    mmftdd(1:lmax,1:ddm,1:odm))
        call ncattr('long_name',
     .    'Overturning stream-function due to thickness diffusion '//
     .    'on z-levels') 
        call ncattr('units','kg s-1')  
      endif
      if (MSC_MHFLX(iogrp).ne.0) then 
        call ncwrt1('mhflx','lat region time',mhflx(1:lmax,1:odm))
        call ncattr('long_name','Meridional heat flux') 
        call ncattr('units','W') 
      endif
      if (MSC_MHFTD(iogrp).ne.0) then 
        call ncwrt1('mhftd','lat region time',mhftd(1:lmax,1:odm))
        call ncattr('long_name',
     .    'Meridional heat flux due to thickness diffusion') 
        call ncattr('units','W') 
      endif
      if (MSC_MHFLD(iogrp).ne.0) then 
        call ncwrt1('mhfld','lat region time',mhfld(1:lmax,1:odm))
        call ncattr('long_name',
     .    'Meridional heat flux due to lateral diffusion') 
        call ncattr('units','W') 
      endif
      if (MSC_MSFLX(iogrp).ne.0) then 
        call ncwrt1('msflx','lat region time',msflx(1:lmax,1:odm))
        call ncattr('long_name','Meridional salt flux') 
        call ncattr('units','kg s-1') 
      endif
      if (MSC_MSFTD(iogrp).ne.0) then 
        call ncwrt1('msftd','lat region time',msftd(1:lmax,1:odm))
        call ncattr('long_name',
     .    'Meridional salt flux due to thickness diffusion') 
        call ncattr('units','kg s-1') 
      endif
      if (MSC_MSFLD(iogrp).ne.0) then 
        call ncwrt1('msfld','lat region time',msfld(1:lmax,1:odm))
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
c --- ideal age tracer 
#if defined(TRC) && defined(IDLAGE)
      call inilyr(ACC_UTIL(1),'p',0.)
      call acclyr(ACC_UTIL,trc(1-nbdy,1-nbdy,1,itriag),
     .  dp(1-nbdy,1-nbdy,k1m),0,'p')
      call wrtlyr(ACC_UTIL(1),LYR_IDLAGE(iogrp),1.,0.,
     .  cmpflg,ip,'p','idlage','Ideal age',
     .  'sea_water_age_since_surface_contact','year')
c
      call inilyr(ACC_UTIL(1),'p',0.)
      call acclyr(ACC_UTIL,dp(1-nbdy,1-nbdy,k1m),
     .  dp(1-nbdy,1-nbdy,k1m),0,'p')
      call wrtlyr(ACC_UTIL(1),LYR_IDLAGE(iogrp),1.,0.,
     .  cmpflg,ip,'p','dp_trc','Layer pressure thickness',' ','Pa')
c
      if (LVL_IDLAGE(iogrp).ne.0) then
        call inilvl(ACC_UTILLVL(1),'p',0.)
        do k=1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
          call acclvl(ACC_UTILLVL,trc(1-nbdy,1-nbdy,1,itriag),'p',
     .      k,ind1,ind2,wghts)
        enddo
        call msklvl(ACC_UTILLVL(1),'p')
        call wrtlvl(ACC_UTILLVL(1),LVL_IDLAGE(iogrp),
     .  1.,0.,cmpflg,ip,'p','idlagelvl','Ideal age',
     .   'sea_water_age_since_surface_contact','year')
      endif
#endif
c
c --- close netcdf file
      call ncfcls
c
c --- initialisation of 2d fields
      call inisrf(ACC_UB(iogrp),'u',0.)
      call inisrf(ACC_ZTX(iogrp),'u',0.)
      call inisrf(ACC_TAUX(iogrp),'u',0.)
      call inisrf(ACC_UICE(iogrp),'u',0.)
      call inisrf(ACC_MXLU(iogrp),'u',0.)
      call inisrf(ACC_MLDU(iogrp),'u',0.)
      call inisrf(ACC_IVOLU(iogrp),'u',0.)
c
      call inisrf(ACC_VB(iogrp),'v',0.)
      call inisrf(ACC_MTY(iogrp),'v',0.)
      call inisrf(ACC_TAUY(iogrp),'v',0.)
      call inisrf(ACC_VICE(iogrp),'v',0.)
      call inisrf(ACC_MXLV(iogrp),'v',0.)
      call inisrf(ACC_MLDV(iogrp),'v',0.)
      call inisrf(ACC_IVOLV(iogrp),'v',0.)
c
      call inisrf(ACC_SEALV(iogrp),'p',0.)
      call inisrf(ACC_SIGMX(iogrp),'p',0.)
      call inisrf(ACC_HICE(iogrp),'p',0.)
      call inisrf(ACC_HSNW(iogrp),'p',0.)
      call inisrf(ACC_FICE(iogrp),'p',0.)
      call inisrf(ACC_TSRF(iogrp),'p',0.)
      call inisrf(ACC_TICE(iogrp),'p',0.)
      call inisrf(ACC_SWA(iogrp),'p',0.)
      call inisrf(ACC_NSF(iogrp),'p',0.)
      call inisrf(ACC_DFL(iogrp),'p',0.)
      call inisrf(ACC_LIP(iogrp),'p',0.)
      call inisrf(ACC_SOP(iogrp),'p',0.)
      call inisrf(ACC_EVA(iogrp),'p',0.)
      call inisrf(ACC_SFL(iogrp),'p',0.)
      call inisrf(ACC_ALB(iogrp),'p',0.)
      call inisrf(ACC_RNFFLX(iogrp),'p',0.)
      call inisrf(ACC_USTAR(iogrp),'p',0.)
      call inisrf(ACC_SSS(iogrp),'p',0.)
      call inisrf(ACC_SST(iogrp),'p',0.)
      call inisrf(ACC_MLD(iogrp),'p',0.)
      call inisrf(ACC_MAXMLD(iogrp),'p',0.)
      call inisrf(ACC_SURFLX(iogrp),'p',0.)
      call inisrf(ACC_SURRLX(iogrp),'p',0.)
      call inisrf(ACC_SALFLX(iogrp),'p',0.)
      call inisrf(ACC_SALRLX(iogrp),'p',0.)
      call inisrf(ACC_ABSWND(iogrp),'p',0.)
      call inisrf(ACC_BRNFLX(iogrp),'p',0.)
      call inisrf(ACC_BRNPD(iogrp),'p',0.)
      call inisrf(ACC_FMLTFZ(iogrp),'p',0.)
      call inisrf(ACC_HMLTFZ(iogrp),'p',0.)
      call inisrf(ACC_RFIFLX(iogrp),'p',0.)
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
      call inilyr(ACC_DIAFLX(iogrp),'p',0.)
      call inilyr(ACC_DIFDIA(iogrp),'p',0.)
      call inilyr(ACC_DIFINT(iogrp),'p',0.)
      call inilyr(ACC_DIFISO(iogrp),'p',0.)
      call inilyr(ACC_WFLX(iogrp),'p',0.)
      call inilyr(ACC_WFLX2(iogrp),'p',0.)
#if defined TKE
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
      call inilvl(ACC_DIAFLXLVL(iogrp),'p',0.)
      call inilvl(ACC_DIFDIALVL(iogrp),'p',0.)
      call inilvl(ACC_DIFINTLVL(iogrp),'p',0.)
      call inilvl(ACC_DIFISOLVL(iogrp),'p',0.)
      call inilvl(ACC_DZLVL(iogrp),'p',0.)
      call inilvl(ACC_SALNLVL(iogrp),'p',0.)
      call inilvl(ACC_TEMPLVL(iogrp),'p',0.)
      call inilvl(ACC_WFLXLVL(iogrp),'p',0.)
      call inilvl(ACC_WFLX2LVL(iogrp),'p',0.)
#if defined TKE
      call inilyr(ACC_TKELVL(iogrp),'p',0.)
      call inilyr(ACC_GLS_PSILVL(iogrp),'p',0.)
#endif
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
      if (mnproc.eq.1.and.iniflg) then
        equat_sec=-1
        open(nfu,file=path(1:path_len)//'secindex.dat',status='old')
        sec_num=0
        do
          read(nfu,'(a120)',iostat=iostatus) char120
          if (iostatus.ne.0) exit
          if (char120(1:4).eq.'Name') then
            sec_num=sec_num+1
            sec_name(sec_num)=char120(7:120)
            if (char120(7:29).eq.'equatorial_undercurrent') then
              equat_sec = sec_num
            endif
            nsi(sec_num)=0
          else
            nsi(sec_num)=nsi(sec_num)+1
            read(char120,*) isi(sec_num,nsi(sec_num)),
     .                      jsi(sec_num,nsi(sec_num)),
     .                      usi(sec_num,nsi(sec_num)),
     .                      vsi(sec_num,nsi(sec_num))
          endif
        enddo
        close(nfu)
        iniflg=.false.
        write(lp,*) 'number of sections=',sec_num
      endif
c      
c --- Prepare 2d field   
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
      integer :: i,j,k,l,o,s,iostatus,iind1,jind1,uflg1,vflg1,nind1,
     .  nfld,ACC_UIND,ACC_VIND,nind(ldm),iind(sdm,ldm),jind(sdm,ldm),
     .  kmax(ldm,odm),kmxl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer(kind=2) :: uflg(sdm,ldm),vflg(sdm,ldm),oflg(itdm,jtdm)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ucum,vcum
      real, dimension(itdm,jtdm) :: depthst,ucumg,vcumg
      real, dimension(ldm,odm) :: mflx_lo,mflx_last_lo
      real :: r
c
      character :: c20*20
c
      logical :: iniflg=.true.
c
      save nind,iind,jind,oflg,uflg,vflg,kmax,depthst,iniflg
c
      if (iniflg.and.mnproc.eq.1) then
c
c --- - Read ocean masks from metraoceans.dat 
        open(nfu,file=path(1:path_len)//'mertraoceans.dat',
     .       status='old')
        read(nfu,'(2i6)') i,j
        if (i.ne.itdm.or.j.ne.jtdm) then
          write(lp,*) 'wrong dimensions in mertraoceans.dat'
          call xchalt('(diamer)')
                 stop '(diamer)'
        endif
        write(c20,*) '(',jtdm,'i1)'
        read(nfu,c20) ((oflg(i,j),j=1,jtdm),i=1,itdm)
        close(nfu)
c 
c --- - Merge ocean regions (1=Atlantic+Arctic,2=Indian+Pacific,3=Global Ocean)
        do j=1,jtdm
          do i=1,itdm
            if     (oflg(i,j).eq.1) then 
              oflg(i,j)=3
            elseif (oflg(i,j).eq.2) then 
              oflg(i,j)=1
            elseif (oflg(i,j).eq.3) then 
              oflg(i,j)=2
            elseif (oflg(i,j).eq.4) then 
              oflg(i,j)=2
            elseif (oflg(i,j).ne.0) then 
              write(lp,*) 'unknown region index in mertraoceans.dat'
              call xchalt('(diamer)')
                     stop '(diamer)'
            endif 
          enddo 
        enddo
c
c --- - Read section file metraindex.dat
        open(nfu,file=path(1:path_len)//'mertraindex.dat',status='old')
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
            write(lp,*) 'problem in mertraindex.dat'
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
            if (oflg(iind1,jind1).ne.0) then
              nind(l)=nind(l)+1
              iind(nind(l),l)=iind1
              jind(nind(l),l)=jind1
              uflg(nind(l),l)=uflg1
              vflg(nind(l),l)=vflg1
              if (iind(nind(l),l).lt.1.or.iind(nind(l),l).gt.itdm.or.
     .            jind(nind(l),l).lt.1.or.jind(nind(l),l).gt.jtdm) then
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
      end if
c
c --- Compute vertical integrated heat and salt transports 
c
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
            do o=1,odm 
              mflx_lo(l,o)=0.
            enddo 
            do s=1,nind(l) 
              i=iind(s,l)
              j=jind(s,l)
              o=oflg(i,j)
              mflx_lo(l,o)=mflx_lo(l,o)
     .                    +uflg(s,l)*ucumg(i,j)+vflg(s,l)*vcumg(i,j)
            enddo
            do o=1,odm-1
              mflx_lo(l,odm)=mflx_lo(l,odm)+mflx_lo(l,o)
            enddo 
c --- ----- Set land to fill value
            do o=1,odm 
              if (abs(mflx_lo(l,o)).lt.epsil) then
                mflx_lo(l,o)=nf90_fill_double
              endif
            enddo
c --- ----- Mask Southern Ocean for Atlantic region 
            if (mtlat(l).lt.-34.) then 
              mflx_lo(l,1)=nf90_fill_double
              mflx_lo(l,2)=nf90_fill_double
            endif
          enddo
c
          if     (nfld.eq.1) then
            do l=1,lmax
              do o=1,odm
                mhflx(l,o)=mflx_lo(l,o)
              enddo
            enddo
          elseif (nfld.eq.2) then
            do l=1,lmax
              do o=1,odm
                mhftd(l,o)=mflx_lo(l,o)
              enddo
            enddo
          elseif (nfld.eq.3) then
            do l=1,lmax
              do o=1,odm
                mhfld(l,o)=mflx_lo(l,o)
              enddo
            enddo
          elseif (nfld.eq.4) then
            do l=1,lmax
              do o=1,odm
                msflx(l,o)=mflx_lo(l,o)
              enddo
            enddo
          elseif (nfld.eq.5) then
            do l=1,lmax
              do o=1,odm 
                msftd(l,o)=mflx_lo(l,o)
              enddo
            enddo
          elseif (nfld.eq.6) then
            do l=1,lmax
              do o=1,odm 
                msfld(l,o)=mflx_lo(l,o)
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
c$OMP PARALLEL DO PRIVATE(k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          k=1
          do while (physrf(i,j,ACC_SIGMX(IOGRP))*r.ge.sigmar(i,j,k))
            k=k+1
            if (k.eq.kk) exit
          enddo
          kmxl(i,j)=k
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=kmxl(i,j)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xctilr(util1, 1,1, 1,1, halo_ps)
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
        do o=1,odm
          do l=1,lmax
            mflx_last_lo(l,o)=0.
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
              do o=1,odm
                mflx_lo(l,o)=0.
              enddo
              do s=1,nind(l)
                i=iind(s,l)
                j=jind(s,l)
                o=oflg(i,j)
                mflx_lo(l,o)=mflx_lo(l,o)+uflg(s,l)*ucumg(i,j)
     .                                   +vflg(s,l)*vcumg(i,j)
              enddo
              do o=1,odm-1
                mflx_lo(l,odm)=mflx_lo(l,odm)+mflx_lo(l,o)
              enddo
c --- ------- Set land to fill value
              do o=1,odm
                if (abs(mflx_lo(l,o)-mflx_last_lo(l,o)).lt.1.e5*epsil)
     .          then
                  mflx_last_lo(l,o)=mflx_lo(l,o)
                  mflx_lo(l,o)=nf90_fill_double
                else
                  mflx_last_lo(l,o)=mflx_lo(l,o)
                endif
              enddo
c
c --- ------- Mask Southern Ocean  
              if (mtlat(l).lt.-34.) then
                mflx_lo(l,1)=nf90_fill_double 
                mflx_lo(l,2)=nf90_fill_double 
              endif
            enddo
c
            if     (nfld.eq.1) then
              do l=1,lmax
                do o=1,odm
                  mmflxl(l,k,o)=mflx_lo(l,o)
                enddo
              enddo
            elseif (nfld.eq.2) then
              do l=1,lmax
                do o=1,odm
                  mmftdl(l,k,o)=mflx_lo(l,o)
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
          do o=1,odm
            kmax(l,o)=0
          enddo
        enddo 
        do k=1,ddm
          do l=1,lmax
            do s=1,nind(l)
              i=iind(s,l)
              j=jind(s,l)
              o=oflg(i,j)
              if (depthslev_bnds(1,k).lt.depthst(i,j)) kmax(l,o)=k
            enddo 
            do o=1,odm-1
              kmax(l,odm)=max(kmax(l,odm),kmax(l,o))
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
              do o=1,odm
                mflx_lo(l,o)=0.
              enddo
              do s=1,nind(l)
                i=iind(s,l)
                j=jind(s,l)
                o=oflg(i,j)
                mflx_lo(l,o)=mflx_lo(l,o)+uflg(s,l)*ucumg(i,j)
     .                      +vflg(s,l)*vcumg(i,j)
              enddo
              do o=1,odm-1
                mflx_lo(l,odm)=mflx_lo(l,odm)+mflx_lo(l,o)
              enddo
c --- ------- Set land to fill value
              do o=1,odm 
                if (kmax(l,o).lt.k) then
                  mflx_lo(l,o)=nf90_fill_double
                endif
              enddo
c --- ------- Mask Southern Ocean for Atlantic region 
              if (mtlat(l).lt.-34.) then 
                mflx_lo(l,1)=nf90_fill_double 
                mflx_lo(l,2)=nf90_fill_double 
              endif
            enddo
c
            if     (nfld.eq.1) then
              do l=1,lmax
                do o=1,odm
                  mmflxd(l,k,o)=mflx_lo(l,o)
                enddo
              enddo
            elseif (nfld.eq.2) then
              do l=1,lmax
                do o=1,odm
                  mmftdd(l,k,o)=mflx_lo(l,o)
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
c$OMP PARALLEL DO
        do j=1,jj
          do i=1,ii
            wflx(i,j)=0.
          enddo
        enddo  
c$OMP END PARALLEL DO
        do k=kk,1,-1
          km=k+mm
          kn=k+nn
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
        do j=1,jj+1
          do i=1,ii+1
            ucum(i,j)=0.
            vcum(i,j)=0.
          enddo
        enddo  
c$OMP END PARALLEL DO
        do k=ddm,1,-1
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
      integer :: d,i,j,k,mm,nn,l,kl,m,km,kn,kml,k1m
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
        do j=1,jj+1
c$OMP PARALLEL DO  
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
c$OMP END PARALLEL DO
        enddo
        iniflg=.false.
      endif
c
c
c --- Compute top and bottom depths of density layers 
      if (k.eq.1) then
        if (gridid.eq.'p') then
c$OMP PARALLEL DO PRIVATE(kml)
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
c$OMP PARALLEL DO PRIVATE(kml)
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
c$OMP PARALLEL DO PRIVATE(kml)
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
c$OMP PARALLEL DO 
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
c$OMP PARALLEL DO PRIVATE(r) 
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
c$OMP PARALLEL DO PRIVATE(r)
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


      subroutine inisrf(pos,gridid,inival)
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
      integer :: i,j,l
c
c --- Check whether field should be initialised
      if (pos.eq.0) return
c
      if (gridid(1:1).eq.'u') then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              physrf(i,j,pos)=inival
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              physrf(i,j,pos)=inival
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              physrf(i,j,pos)=inival
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO
        do j=1,jj
          do i=1,ii
            physrf(i,j,pos)=inival
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      end subroutine inisrf



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
      integer :: i,j,k,l
c
c --- Check whether field should be initialised
      if (pos.eq.0) return
c
      if (gridid(1:1).eq.'u') then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kk 
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                phylyr(i,j,k,pos)=inival
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kk 
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                phylyr(i,j,k,pos)=inival
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kk 
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                phylyr(i,j,k,pos)=inival
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kk 
            do i=1,ii
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
      integer :: i,j,k,l
c
c --- Check whether field should be initialised
      if (pos.eq.0) return
c
      if (gridid(1:1).eq.'u') then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,ddm
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
                phylvl(i,j,k,pos)=inival
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO
        do j=1,jj+1
          do k=1,ddm
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                phylvl(i,j,k,pos)=inival
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'p') then
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,ddm
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                phylvl(i,j,k,pos)=inival
              enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,ddm
            do i=1,ii
              phylvl(i,j,k,pos)=inival
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      end subroutine inilvl



      subroutine accsrf(pos,fld,wghts,wghtsflg,gridid)
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
c$OMP PARALLEL DO 
            do j=1,jj
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                  physrf(i,j,pos(o))=physrf(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO 
            do j=1,jj
              do l=1,isu(j)
                do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                  physrf(i,j,pos(o))=physrf(i,j,pos(o))+fld(i,j)*
     .              wghts(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        elseif (gridid.eq.'v') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO 
            do j=1,jj
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                  physrf(i,j,pos(o))=physrf(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO 
            do j=1,jj
              do l=1,isv(j)
                do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                  physrf(i,j,pos(o))=physrf(i,j,pos(o))+fld(i,j)*
     .              wghts(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        elseif (gridid.eq.'p') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO 
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  physrf(i,j,pos(o))=physrf(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO 
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  physrf(i,j,pos(o))=physrf(i,j,pos(o))+fld(i,j)*
     .              wghts(i,j)
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          endif 
        else
          write (lp,*) 'cannot identify grid '//gridid//'!'
          call xchalt('(diaini)')
                 stop '(diaini)'
        endif
c
      enddo
c   
      end subroutine accsrf



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
c$OMP PARALLEL DO 
            do j=1,jj
              do k=1,kdm
                do l=1,isu(j)
                  do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO 
            do j=1,jj
              do k=1,kdm
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
c$OMP PARALLEL DO 
            do j=1,jj
              do k=1,kdm
                do l=1,isv(j)
                  do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
            do k=1,kdm
c$OMP PARALLEL DO 
              do j=1,jj
                do l=1,isv(j)
                  do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
     .                *wghts(i,j,k)
                  enddo
                enddo
              enddo
c$OMP END PARALLEL DO
            enddo
          endif 
        elseif (gridid.eq.'p') then
          if (wghtsflg.eq.0) then 
c$OMP PARALLEL DO 
            do j=1,jj
              do k=1,kdm
                do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    phylyr(i,j,k,pos(o))=phylyr(i,j,k,pos(o))+fld(i,j,k)
                  enddo
                enddo
              enddo
            enddo
c$OMP END PARALLEL DO
          else
c$OMP PARALLEL DO 
            do j=1,jj
              do k=1,kdm
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
          call xchalt('(diaini)')
                 stop '(diaini)'
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
          call xchalt('(diaini)')
                 stop '(diaini)'
        endif
      enddo
c     
      end subroutine acclvl


      
      subroutine finsrf(posacc,poswgt,gridid)
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
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              physrf(i,j,posacc)=physrf(i,j,posacc)/
     .          max(epsil,physrf(i,j,poswgt))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid(1:1).eq.'v') then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              physrf(i,j,posacc)=physrf(i,j,posacc)/
     .          max(epsil,physrf(i,j,poswgt))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              physrf(i,j,posacc)=physrf(i,j,posacc)/
     .          max(epsil,physrf(i,j,poswgt))
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
c     
      end subroutine finsrf
     


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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
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
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
c     
      end subroutine finlyr



      subroutine wrtsrf(pos,frmt,sfac,offs,cmpflg,msk,gridid,
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
      integer :: frmt,cmpflg,pos,n
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
          call nccopa(vnm,dims,physrf(1-nbdy,1-nbdy,pos),msk,sfac,
     .      offs)
        else
          call ncpack(vnm,dims,physrf(1-nbdy,1-nbdy,pos),msk,
     .      1,sfac,offs)
        endif
      elseif (frmt.eq.4) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,physrf(1-nbdy,1-nbdy,pos),msk,sfac,
     .      offs,4)
        else
          call ncwrtr(vnm,dims,physrf(1-nbdy,1-nbdy,pos),msk,
     .      1,sfac,offs,4)
        endif
      elseif (frmt.eq.8) then
        if (cmpflg.eq.1) then
          call nccomp(vnm,dims,physrf(1-nbdy,1-nbdy,pos),msk,sfac,
     .      offs,8)
        else
          call ncwrtr(vnm,dims,physrf(1-nbdy,1-nbdy,pos),msk,
     .      1,sfac,offs,8)
        endif
      else
        write (lp,*) 'unknown output format!'
        call xchalt('(wrtsrf)')
               stop '(wrtsrf)'
      endif
c
c --- Define attributes
      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
      call ncattr('coordinates',
     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
c
      end subroutine wrtsrf



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
      integer :: frmt,cmpflg,pos,n
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
      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
      call ncattr('coordinates',
     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
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
      integer :: frmt,cmpflg,pos,n
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
      if (len(trim(vunits)).ne.0) call ncattr('units',vunits)
      if (len(trim(vlngnm)).ne.0) call ncattr('long_name',vlngnm)
      if (len(trim(vstdnm)).ne.0) call ncattr('standard_name',vstdnm)
      call ncattr('coordinates',
     .  gridid(1:1)//'lon '//gridid(1:1)//'lat')
      call ncattr('cell_measures','area: '//gridid(1:1)//'area')
c
      end subroutine wrtlvl



      subroutine logsrf(pos,gridid,sfac,offs)
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
c$OMP PARALLEL DO 
        do j=1,jj
          do l=1,isu(j)
            do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
              if (physrf(i,j,pos).lt.epsil) then
                physrf(i,j,pos)=0.
              else
                physrf(i,j,pos)=log10(physrf(i,j,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO 
        do j=1,jj
          do l=1,isv(j)
            do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              if (physrf(i,j,pos).lt.epsil) then
                physrf(i,j,pos)=0.
              else
                physrf(i,j,pos)=log10(physrf(i,j,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO 
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (physrf(i,j,pos).lt.epsil) then
                physrf(i,j,pos)=0.
              else
                physrf(i,j,pos)=log10(physrf(i,j,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      else
        write (lp,*) 'cannot identify grid '//gridid//'!'
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
c
      end subroutine logsrf


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
c$OMP PARALLEL DO 
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
c$OMP PARALLEL DO 
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
c$OMP PARALLEL DO 
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
        call xchalt('(diaini)')
               stop '(diaini)'
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
c$OMP PARALLEL DO 
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
c$OMP PARALLEL DO 
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
c$OMP PARALLEL DO 
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
        call xchalt('(diaini)')
               stop '(diaini)'
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
      integer :: i,j,k,l
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
c$OMP PARALLEL DO
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
c$OMP PARALLEL DO
        do j=1,jj
          do i=1,ii
            do k=kmaxu(i,j)+1,ddm
              phylvl(i,j,k,pos)=mskval
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'v') then
c$OMP PARALLEL DO
        do j=1,jj
          do i=1,ii
            do k=kmaxv(i,j)+1,ddm
              phylvl(i,j,k,pos)=mskval
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      elseif (gridid.eq.'p') then
c$OMP PARALLEL DO
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
        call xchalt('(diaini)')
               stop '(diaini)'
      endif
c
      end subroutine msklvl
 
      end module mod_dia
