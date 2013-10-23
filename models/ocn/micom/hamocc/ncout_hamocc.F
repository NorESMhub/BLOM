      subroutine ncwrt_bgc(iogrp)
c
c --- -------------------------------------------
c --- output routine for HAMOCC diagnostic fields
c --- -------------------------------------------
c
      use mod_xc
      use mod_dia, only : diafnm,sigmar1
      use mo_control_bgc, only : dtbgc
      use mo_bgcmean
c      
      implicit none
c      
#include "common_blocks.h"
#include "common_bgc.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_ice.h"
c
      integer iogrp
c
      integer i,j,k,l,nt
      integer nhour,ny,nm,nd,dayfrac,irec(nbgcmax),cmpflg
      character fname(nbgcmax)*80,startdate*20,timeunits*30
      real datenum,rnacc
      logical append2file(nbgcmax)
      data append2file /nbgcmax*.false./
      save fname,irec,append2file
c
c --- set time information
      timeunits=' '
      startdate=' ' 
      write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)')
     . 'days since ',min(1800,nyear0),'-',1,'-',1,' 00:00'
      write(startdate,'(i4.4,a1,i2.2,a1,i2.2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=time-time0-0.5*diagfq_bgc(iogrp)/nstep_in_day
c
c --- get file name
      if (.not.append2file(iogrp)) then
        call diafnm(runid,runid_len,trim(GLB_FNAMETAG(iogrp)),nstep,
     .    filefq_bgc(iogrp)/real(nstep_in_day),filemon_bgc(iogrp),
     .    fileann_bgc(iogrp),fname(iogrp))
        append2file(iogrp)=.true.
        irec(iogrp)=1
      else
        irec(iogrp)=irec(iogrp)+1
      endif
      if (((fileann_bgc(iogrp).and.nday_of_year.eq.1.or.
     .  filemon_bgc(iogrp).and.nday.eq.1).and.mod(nstep,nstep_in_day)
     .  .le.1).or..not.(fileann_bgc(iogrp).or.filemon_bgc(iogrp)).and.
     .  mod(nstep+.5,filefq_bgc(iogrp)).lt.2.) then
        append2file(iogrp)=.false.
      endif
c
c --- prepare output fields 
      if (mnproc.eq.1) then 
        write (lp,'(a,f6.2,a)') ' ncwrt_bgc: fields averaged over ',
     .                          real(nacc_bgc(iogrp)),' steps'
      endif
      rnacc=1./real(nacc_bgc(iogrp))
      cmpflg=GLB_COMPFLAG(iogrp)
c
c --- create output file
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
c --- define spatial and time dimensions
      if (cmpflg.ne.0) then 
        call ncdimc('pcomp',ip,0)
      else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
      endif
      call ncdims('sigma',kdm)
      call ncdims('depth',ddm)
      call ncdims('ks',ks)
      call ncdims('bounds',2)
      call nctime(datenum,calendar,timeunits,startdate)
c
c --- write auxillary dimension information 
      call ncwrt1('sigma','sigma',sigmar1)
      call ncattr('long_name','Potential density')
      call ncattr('standard_name','sea_water_sigma_theta')
      call ncattr('units','kg m-3')
      call ncattr('positive','down')
      call ncwrt1('depth','depth',depthslev)
      call ncattr('long_name','z level')
      call ncattr('units','m')
      call ncattr('positive','down')
      call ncattr('bounds','depth_bnds')
      call ncwrt1('depth_bnds','bounds depth',depthslev_bnds)
c
c --- finalize accumulation 
      call finlyr(jphyto(iogrp),jdp(iogrp))
      call finlyr(jgrazer(iogrp),jdp(iogrp))
      call finlyr(jdoc(iogrp),jdp(iogrp))
      call finlyr(jphosy(iogrp),jdp(iogrp))
      call finlyr(jphosph(iogrp),jdp(iogrp))
      call finlyr(joxygen(iogrp),jdp(iogrp))
      call finlyr(jiron(iogrp),jdp(iogrp))
      call finlyr(jano3(iogrp),jdp(iogrp))
      call finlyr(jalkali(iogrp),jdp(iogrp))
      call finlyr(jsilica(iogrp),jdp(iogrp))
      call finlyr(jdic(iogrp),jdp(iogrp))
      call finlyr(jpoc(iogrp),jdp(iogrp))
      call finlyr(jcalc(iogrp),jdp(iogrp))
      call finlyr(jopal(iogrp),jdp(iogrp))
      call finlyr(jco3(iogrp),jdp(iogrp))
      call finlyr(jph(iogrp),jdp(iogrp))
      call finlyr(jomegac(iogrp),jdp(iogrp))
      call finlyr(jn2o(iogrp),jdp(iogrp))
#ifdef __c_isotopes
      call finlyr(jdic13(iogrp),jdp(iogrp))
      call finlyr(jdic14(iogrp),jdp(iogrp))
#endif  
#ifdef CFC
      call finlyr(jcfc11(iogrp),jdp(iogrp))
      call finlyr(jcfc12(iogrp),jdp(iogrp))
      call finlyr(jsf6(iogrp),jdp(iogrp))
#endif
c
c --- Mask sea floor in level data
      call msklvl(jlvlphyto(iogrp),depths)
      call msklvl(jlvlgrazer(iogrp),depths)
      call msklvl(jlvldoc(iogrp),depths)
      call msklvl(jlvlphosy(iogrp),depths)
      call msklvl(jlvlphosph(iogrp),depths)
      call msklvl(jlvloxygen(iogrp),depths)
      call msklvl(jlvliron(iogrp),depths)
      call msklvl(jlvlano3(iogrp),depths)
      call msklvl(jlvlalkali(iogrp),depths)
      call msklvl(jlvlsilica(iogrp),depths)
      call msklvl(jlvldic(iogrp),depths)
      call msklvl(jlvlpoc(iogrp),depths)
      call msklvl(jlvlcalc(iogrp),depths)
      call msklvl(jlvlopal(iogrp),depths)
      call msklvl(jlvlco3(iogrp),depths)
      call msklvl(jlvlph(iogrp),depths)
      call msklvl(jlvlomegac(iogrp),depths)
      call msklvl(jlvln2o(iogrp),depths)
#ifdef __c_isotopes
      call msklvl(jlvldic13(iogrp),depths)
      call msklvl(jlvldic14(iogrp),depths)
#endif  
#ifdef AGG
      call msklvl(jlvlnos(iogrp),depths)
#endif  
#ifdef CFC
      call msklvl(jlvlcfc11(iogrp),depths)
      call msklvl(jlvlcfc12(iogrp),depths)
      call msklvl(jlvlsf6(iogrp),depths)
#endif
c
c --- Compute log10 of pH
      if (LYR_PH(iogrp).ne.0) call loglyr(jph(iogrp),1.,0.) 
      if (LVL_PH(iogrp).ne.0) call loglvl(jlvlph(iogrp),rnacc,0.) 
c
c --- Store 2d fields
      call wrtsrf(jkwco2(iogrp),SRF_KWCO2(iogrp),rnacc,0.,cmpflg,       
     .  'kwco2',' ',' ',' ')
      call wrtsrf(jpco2(iogrp),SRF_PCO2(iogrp),rnacc,0.,cmpflg,         
     .  'pco2','Surface PCO2',' ','uatm')
      call wrtsrf(jdmsflux(iogrp),SRF_DMSFLUX(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'dmsflux','DMS flux',' ','mol DMS m-2 s-1')
      call wrtsrf(jco2fxd(iogrp),SRF_CO2FXD(iogrp),rnacc*12./dtbgc,0.,
     .  cmpflg,'co2fxd','Downward CO2 flux',' ','kg C m-2 s-1')
      call wrtsrf(jco2fxu(iogrp),SRF_CO2FXU(iogrp),rnacc*12./dtbgc,0.,
     .  cmpflg,'co2fxu','Upward CO2 flux',' ','kg C m-2 s-1')
      call wrtsrf(joxflux(iogrp),SRF_OXFLUX(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'fgo2','Oxygen flux',' ','mol O2 m-2 s-1')
      call wrtsrf(jniflux(iogrp),SRF_NIFLUX(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'fgn2','Nitrogen flux',' ','mol N2 m-2 s-1')
      call wrtsrf(jdms(iogrp),SRF_DMS(iogrp),rnacc,0.,cmpflg,           
     .  'dms','DMS',' ',' ')
      call wrtsrf(jdmsprod(iogrp),SRF_DMSPROD(iogrp),rnacc,0.,cmpflg,   
     .  'dmsprod','DMS production',' ',' ')
      call wrtsrf(jdms_bac(iogrp),SRF_DMS_BAC(iogrp),rnacc,0.,cmpflg,   
     .  'dms_bac',' ',' ',' ')
      call wrtsrf(jdms_uv(iogrp),SRF_DMS_UV(iogrp),rnacc,0.,cmpflg,     
     .  'dms_uv','DMS',' ',' ')
      call wrtsrf(jexport(iogrp),SRF_EXPORT(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'epc100','Export production',' ','mol C m-2 s-1')
      call wrtsrf(jexposi(iogrp),SRF_EXPOSI(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'epsi100','Si export production',' ','mol Si m-2 s-1')
      call wrtsrf(jexpoca(iogrp),SRF_EXPOCA(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'epcalc100','Ca export production',' ','mol Ca m-2 s-1')
#ifdef CFC
      call wrtsrf(jcfc11fx(iogrp),SRF_CFC11(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'cfc11flux','CFC-11 flux',' ','mol CFC12 m-2 s-1')
      call wrtsrf(jcfc12fx(iogrp),SRF_CFC12(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'cfc12flux','CFC-12 flux',' ','mol CFC12 m-2 s-1')
      call wrtsrf(jsf6fx(iogrp),SRF_SF6(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'sf6flux','SF-6 flux',' ','mol SF6 m-2 s-1')
#endif
#if defined(DIFFAT) || defined(CCSMCOUPLED)
      call wrtsrf(jatmco2(iogrp),SRF_ATMCO2(iogrp),rnacc,0.,cmpflg,     
     .  'atmco2','Atmospheric CO2',' ','ppm')
#endif
#ifdef DIFFAT
      call wrtsrf(jatmo2(iogrp),SRF_ATMO2(iogrp),rnacc,0.,cmpflg,       
     .  'atmo2','Atmospheric O2',' ','ppm')
      call wrtsrf(jatmn2(iogrp),SRF_ATMN2(iogrp),rnacc,0.,cmpflg,       
     .  'atmn2','Atmospheric N2',' ','ppm')
#endif 
      call wrtsrf(jn2ofx(iogrp),SRF_N2OFX(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'n2oflux','N2O flux',' ','mol N2O m-2 s-1')
c
c --- Store 3d layer fields
      call wrtlyr(jdp(iogrp),LYR_DP(iogrp),rnacc,0.,cmpflg, 
     .  'pddpo','Layer thickness',' ','m')
      call wrtlyr(jdic(iogrp),LYR_DIC(iogrp),1e3,0.,cmpflg,            
     .  'dissic','Dissolved inorganic carbon',' ','mol C m-3')
      call wrtlyr(jalkali(iogrp),LYR_ALKALI(iogrp),1e3,0.,cmpflg,     
     .  'talk','Alkalinity',' ','eq m-3')
      call wrtlyr(jphosph(iogrp),LYR_PHOSPH(iogrp),1e3,0.,cmpflg,        
     .  'po4','Phosphorus',' ','mol P m-3')
      call wrtlyr(joxygen(iogrp),LYR_OXYGEN(iogrp),1e3,0.,cmpflg,        
     .  'o2','Oxygen',' ','mol O2 m-3')
      call wrtlyr(jano3(iogrp),LYR_ANO3(iogrp),1e3,0.,cmpflg,            
     .  'no3','Nitrate',' ','mol N m-3')
      call wrtlyr(jsilica(iogrp),LYR_SILICA(iogrp),1e3,0.,cmpflg,        
     .  'si','Silicate',' ','mol Si m-3')
      call wrtlyr(jdoc(iogrp),LYR_DOC(iogrp),1e3,0.,cmpflg,              
     .  'dissoc','Dissolved organic carbon',' ','mol P m-3')
      call wrtlyr(jphyto(iogrp),LYR_PHYTO(iogrp),1e3,0.,cmpflg,          
     .  'phyc','Phytoplankton',' ','mol P m-3')
      call wrtlyr(jgrazer(iogrp),LYR_GRAZER(iogrp),1e3,0.,cmpflg,        
     .  'zooc','Zooplankton',' ','mol P m-3')
      call wrtlyr(jpoc(iogrp),LYR_POC(iogrp),1e3,0.,cmpflg,              
     .  'detoc','Detrius',' ','mol P m-3')
      call wrtlyr(jcalc(iogrp),LYR_CALC(iogrp),1e3,0.,cmpflg,            
     .  'calc','CaCO3 shells',' ','mol C m-3')
      call wrtlyr(jopal(iogrp),LYR_OPAL(iogrp),1e3,0.,cmpflg,            
     .  'opal','Opal shells',' ','mol Si m-3')
      call wrtlyr(jiron(iogrp),LYR_IRON(iogrp),1e3,0.,cmpflg,            
     .  'dfe','Dissolved iron',' ','mol Fe m-3')
      call wrtlyr(jphosy(iogrp),LYR_PHOSY(iogrp),1e3/dtbgc,0.,cmpflg,     
     .  'pp','Primary production',' ','mol C m-3 s-1')
      call wrtlyr(jco3(iogrp),LYR_CO3(iogrp),1e3,0.,cmpflg,              
     .  'co3','Carbonate ions',' ','mol C m-3')
      call wrtlyr(jph(iogrp),LYR_PH(iogrp),-1.,0.,cmpflg,             
     .  'ph','pH',' ','-log10([h+])')
      call wrtlyr(jomegac(iogrp),LYR_OMEGAC(iogrp),1.,0.,cmpflg,        
     .  'omegac','OmegaC',' ','1')
      call wrtlyr(jn2o(iogrp),LYR_N2O(iogrp),1e3,0.,cmpflg,            
     .  'n2o','N2O',' ','mol N2O m-3')
#ifdef __c_isotopes
      call wrtlyr(jdic13(iogrp),LYR_DIC13(iogrp),1e3,0.,cmpflg,          
     .  'dissic13','Dissolved C13',' ','mol C m-3')
      call wrtlyr(jdic14(iogrp),LYR_DIC14(iogrp),1e3,0.,cmpflg,          
     .  'dissic14','Dissolved C14',' ','mol C m-3')
#endif
#ifdef AGG
      call wrtlyr(jnos(iogrp),LYR_NOS(iogrp),1.,0.,cmpflg,              
     .  'nos','NOS',' ','???')
#endif
#ifdef CFC
      call wrtlyr(jcfc11(iogrp),LYR_CFC11(iogrp),1e3,0.,cmpflg,            
     .  'cfc11','CFC-11',' ','mol cfc11 m-3')
      call wrtlyr(jcfc12(iogrp),LYR_CFC12(iogrp),1e3,0.,cmpflg,            
     .  'cfc12','CFC-12',' ','mol cfc12 m-3')
      call wrtlyr(jsf6(iogrp),LYR_SF6(iogrp),1e3,0.,cmpflg,            
     .  'sf6','SF-6',' ','mol sf6 m-3')
#endif
c
c --- Store 3d level fields
      call wrtlvl(jlvldic(iogrp),LVL_DIC(iogrp),rnacc*1e3,0.,cmpflg,      
     .  'dissiclvl','Dissolved inorganic carbon',' ','mol C m-3')
      call wrtlvl(jlvlalkali(iogrp),LVL_ALKALI(iogrp),rnacc*1e3,0.,
     .  cmpflg, 'talklvl','Alkalinity',' ','eq m-3')
      call wrtlvl(jlvlphosph(iogrp),LVL_PHOSPH(iogrp),rnacc*1e3,0.,
     .  cmpflg,'po4lvl','Phosphorus',' ','mol P m-3')
      call wrtlvl(jlvloxygen(iogrp),LVL_OXYGEN(iogrp),rnacc*1e3,0.,
     .  cmpflg,'o2lvl','Oxygen',' ','mol O2 m-3')
      call wrtlvl(jlvlano3(iogrp),LVL_ANO3(iogrp),rnacc*1e3,0.,cmpflg,      
     .  'no3lvl','Nitrate',' ','mol N m-3')
      call wrtlvl(jlvlsilica(iogrp),LVL_SILICA(iogrp),rnacc*1e3,0.,
     .  cmpflg, 'silvl','Silicate',' ','mol Si m-3')
      call wrtlvl(jlvldoc(iogrp),LVL_DOC(iogrp),rnacc*1e3,0.,cmpflg,       
     .  'dissoclvl','Dissolved organic carbon',' ','mol P m-3')
      call wrtlvl(jlvlphyto(iogrp),LVL_PHYTO(iogrp),rnacc*1e3,0.,cmpflg,   
     .  'phyclvl','Phytoplankton',' ','mol P m-3')
      call wrtlvl(jlvlgrazer(iogrp),LVL_GRAZER(iogrp),rnacc*1e3,0.,
     .  cmpflg,'zooclvl','Zooplankton',' ','mol P m-3')
      call wrtlvl(jlvlpoc(iogrp),LVL_POC(iogrp),rnacc*1e3,0.,cmpflg,        
     .  'detoclvl','Detrius',' ','mol P m-3')
      call wrtlvl(jlvlcalc(iogrp),LVL_CALC(iogrp),rnacc*1e3,0.,cmpflg,     
     .  'calclvl','CaCO3 shells',' ','mol C m-3')
      call wrtlvl(jlvlopal(iogrp),LVL_OPAL(iogrp),rnacc*1e3,0.,cmpflg,      
     .  'opallvl','Opal shells',' ','mol Si m-3')
      call wrtlvl(jlvliron(iogrp),LVL_IRON(iogrp),rnacc*1e3,0.,cmpflg,      
     .  'dfelvl','Dissolved iron',' ','mol Fe m-3')
      call wrtlvl(jlvlphosy(iogrp),LVL_PHOSY(iogrp),rnacc*1e3/dtbgc,0.,
     .  cmpflg,'pplvl','Primary production',' ','mol C m-3 s-1')
      call wrtlvl(jlvlco3(iogrp),LVL_CO3(iogrp),rnacc*1e3,0.,cmpflg,        
     .  'co3lvl','Carbonate ions',' ','mol C m-3')
      call wrtlvl(jlvlph(iogrp),LVL_PH(iogrp),-1.,0.,cmpflg,          
     .  'phlvl','pH',' ','-log10([h+])')
      call wrtlvl(jlvlomegac(iogrp),LVL_OMEGAC(iogrp),rnacc,0.,cmpflg,  
     .  'omegaclvl','OmegaC',' ','1')
      call wrtlvl(jlvln2o(iogrp),LVL_N2O(iogrp),rnacc*1e3,0.,cmpflg,      
     .  'n2olvl','N2O',' ','mol N2O m-3')
#ifdef __c_isotopes
      call wrtlvl(jlvldic13(iogrp),LVL_DIC13(iogrp),rnacc*1e3,0.,cmpflg,    
     .  'dic13lvl','Dissolved C13',' ','mol C m-3')
      call wrtlvl(jlvldic14(iogrp),LVL_DIC14(iogrp),rnacc*1e3,0.,cmpflg,    
     .  'dic14lvl','Dissolved C14',' ','mol C m-3')
#endif
#ifdef AGG
      call wrtlvl(jlvlnos(iogrp),LVL_NOS(iogrp),rnacc,0.,cmpflg,       
     .  'noslvl','NOS',' ','???')
#endif
#ifdef CFC
      call wrtlvl(jlvlcfc11(iogrp),LVL_CFC11(iogrp),rnacc*1e3,0.,cmpflg,  
     .  'cfc11lvl','CFC-11',' ','mol cfc11 m-3')
      call wrtlvl(jlvlcfc12(iogrp),LVL_CFC12(iogrp),rnacc*1e3,0.,cmpflg,  
     .  'cfc12lvl','CFC-12',' ','mol cfc12 m-3')
      call wrtlvl(jlvlsf6(iogrp),LVL_SF6(iogrp),rnacc*1e3,0.,cmpflg,  
     .  'sf6lvl','SF-6',' ','mol sf6 m-3')
#endif
c
c --- Store sediment fields
      call wrtsdm(jpowaic(iogrp),SDM_POWAIC(iogrp),rnacc*1e3,0.,cmpflg,
     .  'powdic','PoWa DIC',' ','mol C m-3')
      call wrtsdm(jpowaal(iogrp),SDM_POWAAL(iogrp),rnacc*1e3,0.,cmpflg,
     .  'powalk','PoWa alkalinity',' ','eq m-3')
      call wrtsdm(jpowaph(iogrp),SDM_POWAPH(iogrp),rnacc*1e3,0.,cmpflg,
     .  'powpho','PoWa phosphorus',' ','mol P m-3')
      call wrtsdm(jpowaox(iogrp),SDM_POWAOX(iogrp),rnacc*1e3,0.,cmpflg,
     .  'powox','PoWa oxygen',' ','mol O2 m-3')
      call wrtsdm(jpown2(iogrp),SDM_POWN2(iogrp),  rnacc*1e3,0.,cmpflg,
     .  'pown2','PoWa N2',' ','mol N2 m-3')
      call wrtsdm(jpowno3(iogrp),SDM_POWNO3(iogrp),rnacc*1e3,0.,cmpflg,
     .  'powno3','PoWa nitrate',' ','mol N m-3')
      call wrtsdm(jpowasi(iogrp),SDM_POWASI(iogrp),rnacc*1e3,0.,cmpflg,
     .  'powsi','PoWa silicate',' ','mol Si m-3')
      call wrtsdm(jssso12(iogrp),SDM_SSSO12(iogrp),rnacc*1e3,0.,cmpflg,
     .  'ssso12','Sediment detritus',' ','mol P m-3')
      call wrtsdm(jssssil(iogrp),SDM_SSSSIL(iogrp),rnacc*1e3,0.,cmpflg,
     .  'ssssil','Sediment silicate',' ','mol Si m-3')
      call wrtsdm(jsssc12(iogrp),SDM_SSSC12(iogrp),rnacc*1e3,0.,cmpflg,
     .  'sssc12','Sediment CaCO3',' ','mol C m-3')
      call wrtsdm(jssster(iogrp),SDM_SSSTER(iogrp),rnacc*1e3,0.,cmpflg,
     .  'ssster','Sediment clay',' ','mol m-3')
c
c --- Store sediment burial fields
      call wrtbur(jburssso12(iogrp),BUR_SSSO12(iogrp),rnacc*1e3,0.,
     .  cmpflg,'buro12','Burial org carbon',' ','mol P m-2')
      call wrtbur(jbursssc12(iogrp),BUR_SSSC12(iogrp),rnacc*1e3,0.,
     .  cmpflg,'burc12','Burial calcium ',' ','mol C m-2')
      call wrtbur(jburssssil(iogrp),BUR_SSSSIL(iogrp),rnacc*1e3,0.,
     .  cmpflg,'bursil','Burial silicate',' ','mol Si m-2')
      call wrtbur(jburssster(iogrp),BUR_SSSTER(iogrp),rnacc*1e3,0.,
     .  cmpflg,'burter','Burial clay',' ','mol  m-2')
c
c --- close netcdf file
      call ncfcls
c
c --- Initialise fields 
      call inisrf(jkwco2(iogrp),0.)       
      call inisrf(jpco2(iogrp),0.)       
      call inisrf(jdmsflux(iogrp),0.)       
      call inisrf(jco2fxd(iogrp),0.)       
      call inisrf(jco2fxu(iogrp),0.)       
      call inisrf(joxflux(iogrp),0.)       
      call inisrf(jniflux(iogrp),0.)       
      call inisrf(jn2ofx(iogrp),0.)       
      call inisrf(jdms(iogrp),0.)  
      call inisrf(jdmsprod(iogrp),0.)       
      call inisrf(jdms_bac(iogrp),0.)       
      call inisrf(jdms_uv(iogrp),0.)       
      call inisrf(jexport(iogrp),0.)       
      call inisrf(jexposi(iogrp),0.)       
      call inisrf(jexpoca(iogrp),0.)       
#ifdef CFC
      call inisrf(jcfc11fx(iogrp),0.)       
      call inisrf(jcfc12fx(iogrp),0.)       
      call inisrf(jsf6fx(iogrp),0.)       
#endif
#if defined(DIFFAT) || defined(CCSMCOUPLED)
      call inisrf(jatmco2(iogrp),0.)       
#endif   
#ifdef DIFFAT 
      call inisrf(jatmo2(iogrp),0.)       
      call inisrf(jatmn2(iogrp),0.)    
#endif   
c
      call inilyr(jdp(iogrp),0.)       
      call inilyr(jdic(iogrp),0.)       
      call inilyr(jalkali(iogrp),0.)       
      call inilyr(jphosy(iogrp),0.)       
      call inilyr(jphosph(iogrp),0.)       
      call inilyr(joxygen(iogrp),0.)       
      call inilyr(jano3(iogrp),0.)       
      call inilyr(jsilica(iogrp),0.)       
      call inilyr(jdoc(iogrp),0.)       
      call inilyr(jphyto(iogrp),0.)       
      call inilyr(jgrazer(iogrp),0.)       
      call inilyr(jpoc(iogrp),0.)       
      call inilyr(jcalc(iogrp),0.)       
      call inilyr(jopal(iogrp),0.)       
      call inilyr(jiron(iogrp),0.)       
      call inilyr(jco3(iogrp),0.)       
      call inilyr(jph(iogrp),0.)       
      call inilyr(jomegac(iogrp),0.)       
      call inilyr(jn2o(iogrp),0.)       
#ifdef __c_isotopes
      call inilyr(jdic13(iogrp),0.)       
      call inilyr(jdic14(iogrp),0.)       
#endif  
#ifdef AGG
      call inilyr(jnos(iogrp),0.)  
#endif  
#ifdef CFC
      call inilyr(jcfc11(iogrp),0.)       
      call inilyr(jcfc12(iogrp),0.)       
      call inilyr(jsf6(iogrp),0.)       
#endif
c
      call inilvl(jlvldic(iogrp),0.)
      call inilvl(jlvlalkali(iogrp),0.)
      call inilvl(jlvlphosy(iogrp),0.)
      call inilvl(jlvlphosph(iogrp),0.)       
      call inilvl(jlvloxygen(iogrp),0.)
      call inilvl(jlvlano3(iogrp),0.)  
      call inilvl(jlvlsilica(iogrp),0.)
      call inilvl(jlvldoc(iogrp),0.)   
      call inilvl(jlvlphyto(iogrp),0.)
      call inilvl(jlvlgrazer(iogrp),0.)
      call inilvl(jlvlpoc(iogrp),0.)   
      call inilvl(jlvlcalc(iogrp),0.)
      call inilvl(jlvlopal(iogrp),0.)
      call inilvl(jlvliron(iogrp),0.)
      call inilvl(jlvlco3(iogrp),0.)
      call inilvl(jlvlph(iogrp),0.)
      call inilvl(jlvlomegac(iogrp),0.)
      call inilvl(jlvln2o(iogrp),0.)       
#ifdef __c_isotopes
      call inilvl(jlvldic13(iogrp),0.)
      call inilvl(jlvldic14(iogrp),0.)
#endif  
#ifdef AGG
      call inilvl(jlvlnos(iogrp),0.)
#endif  
#ifdef CFC
      call inilvl(jlvlcfc11(iogrp),0.)       
      call inilvl(jlvlcfc12(iogrp),0.)       
      call inilvl(jlvlsf6(iogrp),0.)       
#endif
c
      call inisdm(jpowaic(iogrp),0.)       
      call inisdm(jpowaal(iogrp),0.)       
      call inisdm(jpowaph(iogrp),0.)       
      call inisdm(jpowaox(iogrp),0.)       
      call inisdm(jpown2(iogrp),0.)       
      call inisdm(jpowno3(iogrp),0.)       
      call inisdm(jpowasi(iogrp),0.)       
      call inisdm(jssso12(iogrp),0.)       
      call inisdm(jssssil(iogrp),0.)       
      call inisdm(jsssc12(iogrp),0.)       
      call inisdm(jssster(iogrp),0.)       

      call inibur(jburssso12(iogrp),0.)       
      call inibur(jbursssc12(iogrp),0.)       
      call inibur(jburssssil(iogrp),0.)       
      call inibur(jburssster(iogrp),0.)       
c
      nacc_bgc(iogrp)=0
c
      end
