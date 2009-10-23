      subroutine ncwrt2_bgc(otr_dia2d,nbgc2d,naccbgc2d)
c
c --- ---------------------------------------
c --- output routine for 2d-diagnostic fields
c --- ---------------------------------------
c
      use mod_xc
c      
      implicit none
c      
#include "common_blocks.h"
#include "common_bgc.h"
#include "common_clndr.h"
#include "common_dia.h"
#include "common_forc.h"
#include "common_ice.h"
c
      integer i,j,k,l,nt
      integer ijdm,fid,nhour,ny,nm,nd,dayfrac,nbgc2d,naccbgc2d
      parameter(ijdm=idm*jdm)
      character fname*80,startdate*20,timeunits*30
      real datenum,rnacc_bgc
      real otr_dia2d(idm,jdm,nbgc2d),
     .     tr_pack(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c                     
      INTEGER   jkwco2,jpco2,jdms,jdmsflux,
     .          jdmsprod,jdms_bac,jdms_uv,jco2fxd,
     .          jco2fxu,joxflux,jniflux,jexport,
     .          jexpoca,jexposi
      PARAMETER(jkwco2   =1 ,jpco2    =2 ,jdmsflux =3 ,jco2fxd =4,
     .          jco2fxu  =5 ,joxflux  =6 ,jniflux  =7 ,jdms    =8,
     .          jdmsprod =9 ,jdms_bac =10,jdms_uv  =11,jexport =12,
     .          jexpoca  =13,jexposi  =14)

#ifdef DIFFAT
      INTEGER   jatmco2,jatmo2,jatmn2
      PARAMETER(jatmco2=15,jatmo2=16, jatmn2=17)
#elif CCSMCOUPLED
      INTEGER   jatmco2
      PARAMETER(jatmco2=15)
#endif

#ifdef PCFC
      INTEGER jac14fx,jcfc11fx,jcfc12fx,jpcfc11,jpcfc12
      PARAMETER(jac14fx=22,jcfc11fx=23,jcfc12fx=24,jpcfc11=25,
     .          jpcfc12=26) 
#endif         

#ifdef PCOMPONENT_ANALYSIS
      INTEGER jdpco2_dalk,jdpco2_ddic,jdpco2_dsst,
     .        jdpco2_dsss,jdpco2_dalk_ant,jdpco2_ddic_ant 
      PARAMETER(jdpco2_dalk =27,jdpco2_ddic    =28,jdpco2_dsst    =29,
     .          jdpco2_dsss =30,jdpco2_dalk_ant=31,jdpco2_ddic_ant=32) 
#endif       
c
c --- set time information
      timeunits='hours since 0001-01-01 00:00'
      write (startdate,'(i4,a1,i2,a1,i2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=(time-.5*diag2dfq/nstep_in_day-time0)*24.
c
c --- get file name
      call ncfnam(runid,runid_len,'2b',time,diag2dfq/nstep_in_day,
     .            diag2dmon,fname)
c
c --- prepare output fields 
      if (mnproc.eq.1) then
        write (lp,'(a,f6.2,a)') ' ncwrt2: 2d fields averaged over ',
     .                          real(naccbgc2d),' steps'
      endif
      rnacc_bgc=1./real(naccbgc2d)
c
c     do nt=8,14
cc$OMP PARALLEL DO 
c       do j=1,jj
c         do l=1,isp(j)
c         do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c           otr_dia2d(i,j,nt)=(otr_dia2d(i,j,nt)*onem)
c    .                        /(dp_dia_bgc(i,j,1)+epsil)
c         enddo
c         enddo
c       enddo
cc$OMP END PARALLEL DO
c     enddo
c
c --- create output file
      call ncfopn(path1(1:path1_len)//fname,'w')
      if (mnproc.eq.1) then
        call flush(lp)
      endif
c      
c --- define spatial and time dimensions
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call nctime(datenum,calendar,timeunits,startdate)
c
c --- write data fields
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jkwco2)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('kwco2','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','W m-2') ! ???
      call ncattr('long_name','Heat flux received by ocean') ! ???
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jpco2)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('pco2','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','W m-2') ! ???
      call ncattr('long_name','Surface PCO2')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jco2fxd)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('co2fxd','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol m-2') 
      call ncattr('long_name','Downward CO2 flux')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jco2fxu)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('co2fxu','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol m-2') 
      call ncattr('long_name','Upward CO2 flux')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,joxflux)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('oxflux','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol m-2') 
      call ncattr('long_name','Oxygen flux')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jniflux)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('niflux','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol m-2') 
      call ncattr('long_name','Nitrogen flux')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jexport)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('export','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol C m-2')
      call ncattr('long_name','Export production')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jexposi)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('exposi','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol Si m-2')
      call ncattr('long_name','Si export production')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jexpoca)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('expoca','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol C m-2') 
      call ncattr('long_name','Ca export production')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
#ifdef DIFFAT
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jatmco2)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('atmco2','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','ppm')
      call ncattr('long_name','Atmospheric CO2')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jatmo2)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('atmo2','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','ppm')
      call ncattr('long_name','Atmospheric O2')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jatmn2)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('atmn2','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','ppm')
      call ncattr('long_name','Atmospheric N2')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
#elif CCSMCOUPLED
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tr_pack(i,j)=otr_dia2d(i,j,jatmco2)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('atmco2','x y time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','ppm')
      call ncattr('long_name','Atmospheric CO2')
      call ncattr('coordinates','plon plat')
      call ncattr('cell_measures','area: parea')
c
#endif
c
c     nt=1
cc$OMP PARALLEL DO
c     do j=1,jj
c       do l=1,isp(j)
c       do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c         tr_pack(i,j)=otr_dia2d(i,j,nt)
c       enddo
c       enddo
c     enddo
cc$OMP END PARALLEL DO
c
c     call ncpack('kwrbioz','x y time',tr_pack,ip,1.,0.,1)
c     call ncattr('units','none')
c     call ncattr('long_name','Layer kwrbioz')
c     call ncattr('coordinates','plon plat')
c     call ncattr('cell_measures','area: parea')
c
c     nt=2
cc$OMP PARALLEL DO
c     do j=1,jj
c       do l=1,isp(j)
c       do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c         tr_pack(i,j)=otr_dia2d(i,j,nt)
c       enddo
c       enddo
c     enddo
cc$OMP END PARALLEL DO
c
c     call ncpack('tubioz','x y time',tr_pack,ip,1.,0.,1)
c     call ncattr('units','m')
c     call ncattr('long_name','Depth kwrbioz centre')
c     call ncattr('coordinates','plon plat')
c     call ncattr('cell_measures','area: parea')
c
c     nt=3
cc$OMP PARALLEL DO
c     do j=1,jj
c       do l=1,isp(j)
c       do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c         tr_pack(i,j)=otr_dia2d(i,j,nt)
c       enddo
c       enddo
c     enddo
cc$OMP END PARALLEL DO
c
c     call ncpack('twbioz','x y time',tr_pack,ip,1.,0.,1)
c     call ncattr('units','m')
c     call ncattr('long_name','Depth kwrbioz lower')
c     call ncattr('coordinates','plon plat')
c     call ncattr('cell_measures','area: parea')
c
c --- close netcdf file
      call ncfcls(fid)
c
c --- set accumulation fields to zero
      do nt=1,nbgc2d
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            otr_dia2d(i,j,nt)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      end

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ncwrt3_bgc(ocetr_dia)
c
c --- ---------------------------------------
c --- output routine for 3d-diagnostic fields
c --- ---------------------------------------
c      
      use mod_xc
c
      implicit none
c      
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_dia.h"
c
      integer i,j,k,l,nt
      integer ijdm,ijkdm,fid
      parameter(ijdm=idm*jdm,ijkdm=idm*jdm*kdm)
      character fname*80,startdate*20,timeunits*30
      real datenum
      integer nocetr_dia3d

#ifdef __c_isotopes
      parameter (nocetr_dia3d=19)
#else
      parameter (nocetr_dia3d=17)
#endif
      real ocetr_dia(idm,jdm,kdm,nocetr_dia3d),
     .     tr_pack(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
c                     
      integer isco212,ialkali,iphosph,ioxygen,igasnit,iano3,isilica,
     .        idoc,iphy,izoo,idet,icalc,iopal,isco214,ian2o,idms,
     .        ifdust,iiron
      parameter(isco212=1 ,ialkali= 2,iphosph= 3,ioxygen= 4,igasnit=5 ,
     .          iano3  =6 ,isilica= 7,idoc   = 8,iphy   = 9,izoo   =10,
     .          idet   =11,icalc  =12,iopal  =13,isco214=14,ian2o  =15,
     .          idms   =16,ifdust =17,iiron  =18)

      integer jphyto,jgrazer,jdoc,jphosy,jphosph,joxygen,jiron,jano3,
     .        jalkali,jsilica,jdic,jpoc,jcalc,jopal,jco3,jph,jomegac
      parameter(jphyto = 1,jgrazer= 2,jdoc = 3,jphosy = 4,jphosph= 5,
     .          joxygen= 6,jiron  = 7,jano3= 8,jalkali= 9,jsilica=10,
     .          jdic   =11,jpoc   =12,jcalc=13,jopal  =14,jco3   =15,
     .          jph    =16,jomegac=17)
      
#ifdef __c_isotopes
      integer jdic13,jdic14
      parameter(jdic13=18,jdic14=19)
#endif
#ifdef AGG  
      integer inos,iadust
      parameter (inos=19,iadust=20)
#endif

#ifdef PCFC  
      integer icfc11,icfc12,iantc14
      parameter (icfc11=24,icfc12=25,iantc14=26)
#endif
c
c --- set time information
      timeunits='hours since 0001-01-01 00:00'
      write (startdate,'(i4,a1,i2,a1,i2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=(time-.5*diag3dfq/nstep_in_day-time0)*24.
c
c --- get file name
      call ncfnam(runid,runid_len,'3b',time,diag3dfq/nstep_in_day,
     .            diag3dmon,fname)
c
c --- prepare output fields 
c
      if (nstep.gt.nstep1) then
        do nt=1,nocetr_dia3d
c$OMP PARALLEL DO
          do j=1,jj
            do k=1,kdm
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c               ocetr_dia(i,j,k,nt)=ocetra(i,j,k,nt)
                ocetr_dia(i,j,k,nt)=(ocetr_dia(i,j,k,nt)*(onem*1.027))
     .                              /max(dp_dia_bgc(i,j,k),epsil) 
              enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        enddo
c$OMP PARALLEL DO
          do j=1,jj
            do k=1,kdm
              do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
               if(ocetr_dia(i,j,k,jph).ne.0.0) then
               ocetr_dia(i,j,k,jph)=-log10(ocetr_dia(i,j,k,jph))
               endif
              enddo
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO

      endif
c
      if (mnproc.eq.1) then
        write (lp,'(2a)')
     .    'Writing 3D diagnostics to file: ',path1(1:path1_len)//fname
      endif
c
c --- create output file
      call ncfopn(path1(1:path1_len)//fname,'w')
      if (mnproc.eq.1) then
        call flush(lp)
      endif
c
c --- define spatial and time dimensions
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncdims('z',kdm)
      call nctime(datenum,calendar,timeunits,startdate)    
c
c --- write data fields
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jdic)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('dic','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol C m-3')
      call ncattr('long_name','DIC')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jalkali)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('alkali','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol eq m-3')
      call ncattr('long_name','Alkalinity')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jphosph)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('phosph','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol P m-3')
      call ncattr('long_name','Phosphorus')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,joxygen)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('oxygen','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol O m-3')
      call ncattr('long_name','Oxygen')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jano3)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('ano3','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol N m-3')
      call ncattr('long_name','Nitrate')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jsilica)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('silica','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol Si m-3')
      call ncattr('long_name','Silicate')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jdoc)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('doc','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol C m-3')
      call ncattr('long_name','DOC')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jphyto)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('phy','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol P m-3') ! ???
      call ncattr('long_name','Phytoplankton')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jgrazer)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('zoo','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol P m-3') ! ???
      call ncattr('long_name','Zooplankton')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jpoc)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('det','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol P m-3') ! ???
      call ncattr('long_name','Detritus')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jcalc)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('calc','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol C m-3') ! ???
      call ncattr('long_name','CaCO3 shells')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jopal)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('opal','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol Si m-3') ! ???
      call ncattr('long_name','Opal shells')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
cc$OMP PARALLEL DO
c     do j=1,jj
c       do k=1,kk
c         do l=1,isp(j)
c         do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c           tr_pack(i,j,k)=ocetra(i,j,k,ian2o)
c         enddo
c         enddo
c       enddo
c     enddo
cc$OMP END PARALLEL DO
c
c     call ncpack('an2o','x y z time',tr_pack,ip,1.,0.,1)
c     call ncattr('units','kmol N m-3') ! ???
c     call ncattr('long_name','Nitrous Oxide')
c     call ncattr('coordinates','plon plat sigma')
c     call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jiron)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('iron','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol Fe m-3')
      call ncattr('long_name','Dissolved Iron')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jphosy)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('pripro','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol P m-3') 
      call ncattr('long_name','Primary production')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jco3)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO

      call ncpack('co3','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol C m-3') 
      call ncattr('long_name','Carbonate ions')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c18
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jph)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO

      call ncpack('ph','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','-log10([h+])') 
      call ncattr('long_name','pH')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c19
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jomegac)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO

      call ncpack('omegac','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','no unit') 
      call ncattr('long_name','OmegaC')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
#ifdef __c_isotopes
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jdic13)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('dic13','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol C m-3') ! ???
      call ncattr('long_name','Dissolved C13')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kk
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=ocetr_dia(i,j,k,jdic14)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('dic14','x y z time',tr_pack,ip,1.,0.,1)
      call ncattr('units','kmol C m-3') ! ???
      call ncattr('long_name','Dissolved C14')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
#endif
c
c --- close netcdf file
      call ncfcls
c
c --- set accumulation fields to zero
      do nt=1,nocetr_dia3d
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,kdm
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              ocetr_dia(i,j,k,nt)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,kdm
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp_dia_bgc(i,j,k)=0.0
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      subroutine ncwrt3_sed(sedtr_dia,naccbgc2d)
c
c --- ---------------------------------------
c --- output routine for sediment fields
c --- ---------------------------------------
c
       use mod_xc
c      
      implicit none
c      
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_dia.h"
c
      integer ks
      parameter(ks=12)
      integer i,j,k,l,nt
      integer ijdm,ijksdm,fid
      parameter(ijdm=idm*jdm,ijksdm=idm*jdm*ks)
      character fname*80,startdate*20,timeunits*30
      real datenum,rnacc_bgc
      integer nsedtr_dia3d,naccbgc2d

      parameter (nsedtr_dia3d=11)

      real sedtr_dia(idm,jdm,ks,nsedtr_dia3d),
     .     tr_pack(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ks)
c                     
      INTEGER   jpowaic,jpowaal,jpowaph,jpowaox,jpown2,
     .          jpowno3,jpowasi,jssso12,jssssil,jsssc12,
     .          jssster
      PARAMETER(jpowaic=1,jpowaal=2,jpowaph=3,jpowaox=4,jpown2=5,
     .          jpowno3=6,jpowasi=7,jssso12=8,jssssil=9,jsssc12=10,
     .          jssster=11)
c
c --- set time information
      timeunits='hours since 0001-01-01 00:00'
      write (startdate,'(i4,a1,i2,a1,i2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=(time-.5*diag3dfq/nstep_in_day-time0)*24.
c
c --- get file name
      call ncfnam(runid,runid_len,'3s',time,diag3dfq/nstep_in_day,
     .            diag3dmon,fname)
c
c --- prepare output fields 
      if (mnproc.eq.1) then
        write (lp,'(a,f6.2,a)') ' ncwrt2: 2d fields averaged over ',
     .                          real(naccbgc2d),' steps'
      endif
      rnacc_bgc=1./real(naccbgc2d)
c
      if (mnproc.eq.1) then
        write (lp,'(2a)')
     .    'Writing 3D diagnostics to file: ',path1(1:path1_len)//fname
      endif
c
c --- create output file
      call ncfopn(path1(1:path1_len)//fname,'w')
      if (mnproc.eq.1) then
        call flush(lp)
      endif
c
c --- define spatial and time dimensions
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncdims('z',ks)
      call nctime(datenum,calendar,timeunits,startdate)    
c
c --- write data fields
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpowaic)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('powdic','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol C m-3')
      call ncattr('long_name','PoWa DIC')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpowaal)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('powalk','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol eq m-3')
      call ncattr('long_name','PoWa alkalinity')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpowaph)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('powpho','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol P m-3')
      call ncattr('long_name','PoWa phosphate')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpowaox)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('powox','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol O m-3')
      call ncattr('long_name','PoWa oxygen')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpown2)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('pown2','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol N m-3')
      call ncattr('long_name','PoWa N2')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpowno3)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('powno2','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol N m-3') ! ???
      call ncattr('long_name','PoWa nitrate')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jpowasi)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('powsi','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol Si m-3')
      call ncattr('long_name','PoWa silicate')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jssso12)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('sssol12','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol P m-3')
      call ncattr('long_name','Sediment detritus')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jssssil)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('ssssil','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol Si m-3')
      call ncattr('long_name','Sediment silicate')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jsssc12)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('sssc12','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kmol C m-3')
      call ncattr('long_name','Sediment CaCO3')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c$OMP PARALLEL DO
      do j=1,jj
        do k=1,ks
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tr_pack(i,j,k)=sedtr_dia(i,j,k,jssster)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call ncpack('ssster','x y z time',tr_pack,ip,rnacc_bgc,0.,1)
      call ncattr('units','kg m-3')
      call ncattr('long_name','Sediment clay')
      call ncattr('coordinates','plon plat sigma')
      call ncattr('cell_measures','area: parea')
c
c --- close netcdf file
      call ncfcls
c
c --- set accumulation fields to zero
      do nt=1,nsedtr_dia3d
c$OMP PARALLEL DO
        do j=1,jj
          do k=1,ks
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              sedtr_dia(i,j,k,nt)=0.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c      call diaini3d
c
      end

