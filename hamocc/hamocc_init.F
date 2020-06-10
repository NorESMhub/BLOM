c Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen
c
c This file is part of BLOM/iHAMOCC.
c
c BLOM is free software: you can redistribute it and/or modify it under the
c terms of the GNU Lesser General Public License as published by the Free 
c Software Foundation, either version 3 of the License, or (at your option) 
c any later version. 
c
c BLOM is distributed in the hope that it will be useful, but WITHOUT ANY 
c WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
c FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
c more details. 
c
c You should have received a copy of the GNU Lesser General Public License 
c along with BLOM. If not, see https://www.gnu.org/licenses/.


      subroutine hamocc_init(read_rest,rstfnm_ocn)
c******************************************************************************
c
c  HAMOCC_INIT - initialize HAMOCC and its interface to BLOM.
c
c
c  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-25
c
c
c  Purpose
c  -------
c  - HAMOCC intialization when coupled to BLOM.
c
c
c  Interface to ocean model (parameter list):
c  -----------------------------------------
c  *INTEGER*   *read_rest*  - flag indicating whether to read restart files.
c  *INTEGER*   *rstfnm_ocn* - restart filename.
c
c******************************************************************************
      use mod_xc,        only: ii,jj,kk,idm,jdm,kdm,nbdy,isp,ifp,ilp,
     .                         mnproc,lp,nfu
      USE mo_control_bgc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
      USE mo_carbch,     only: alloc_mem_carbch,ocetra,atm,atm_co2
      USE mo_biomod,     only: alloc_mem_biomod
      USE mo_sedmnt,     only: alloc_mem_sedmnt,sedlay,powtra,burial
      use mo_vgrid,      only: alloc_mem_vgrid,set_vgrid
      use mo_bgcmean,    only: alloc_mem_bgcmean
      use mo_fedep,      only: ini_fedep,fedepfile
      use mo_ndep,       only: ini_ndep,ndepfile
      use mo_Gdata_read, only: inidic,inialk,inipo4,inioxy,inino3,inisil,
     .                         inid13c,inid14c
      use mo_riverinpt,  only: ini_riverinpt,rivinfile
      use mo_intfcblom,  only: alloc_mem_intfcblom,nphys,
     .                         bgc_dx,bgc_dy,bgc_dp,bgc_rho,
     .                         omask,sedlay2,powtra2,burial2
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_geo.h"
c
      integer,          intent(in) :: read_rest
      character(len=*), intent(in) :: rstfnm_ocn

      integer :: i,j,k,l,nt

      namelist /bgcnml/ atm_co2,do_rivinpt,do_ndep,
     .                  ndepfile,fedepfile,rivinfile,
     .                  inidic,inialk,inipo4,inioxy,inino3,inisil,
     .                  inid13c,inid14c

c
c --- Set io units and some control parameters
c
      io_stdo_bgc = lp              !  standard out.
      io_nml = nfu                  !  namelist
      dtbgc = nphys*baclin          !  time step length [sec].
      ndtdaybgc=NINT(86400./dtbgc)  !  time steps per day [No].
      dtb=1./ndtdaybgc              !  time step length [days].
      isac = 1
      ldtbgc = 0
      ldtrunbgc = 0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      WRITE(io_stdo_bgc,*)'********************************************'
      write(io_stdo_bgc,*) 'iHAMOCC: initialisation'
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) 'restart',read_rest
      write(io_stdo_bgc,*) 'dims',idm,jdm,kdm
      write(io_stdo_bgc,*) 'time',nyear,nmonth,nday
      write(io_stdo_bgc,*) 'time step',dtbgc
      endif
c
c --- Read the HAMOCC BGCNML namelist.
c
      open (unit=io_nml,file='ocn_in',status='old',action='read')
      read (unit=io_nml,nml=BGCNML)
      close (unit=io_nml)
      IF (mnproc.eq.1) THEN
        write(io_stdo_bgc,*)
        write(io_stdo_bgc,*) 'iHAMOCC: reading namelist BGCNML'
        write(io_stdo_bgc,nml=BGCNML)
      ENDIF
c                        
c --- Memory allocation 
c
      CALL ALLOC_MEM_INTFCBLOM(idm,jdm,kdm)
      CALL ALLOC_MEM_BGCMEAN(idm,jdm,kdm)
      CALL ALLOC_MEM_VGRID(idm,jdm,kdm)
      CALL ALLOC_MEM_BIOMOD(idm,jdm,kdm)
      CALL ALLOC_MEM_SEDMNT(idm,jdm)
      CALL ALLOC_MEM_CARBCH(idm,jdm,kdm)
c
c --- initialise trc array (two time levels)
c
      do nt=itrbgc,itrbgc+ntrbgc-1
      do k=1,2*kk
      do j=1,jj
      do i=1,ii
        trc(i,j,k,nt)=0.0
      enddo
      enddo
      enddo
      enddo
c
c --- initialise HAMOCC land/ocean mask
c
      do j=1,jj
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        omask(i,j)=1.
      enddo
      enddo
      enddo
c
c --- BLOM to HAMOCC interface
c
      call blom2hamocc(2,1,kk,0)
c
c --- Calculate variables related to the vertical grid
c
      call set_vgrid(idm,jdm,kdm,bgc_dp)
c                        
c --- Initialize sediment layering
c
      CALL BODENSED(idm,jdm,kdm,bgc_dp)
c
c --- Initialize parameters, sediment and ocean tracer.
c 
      CALL BELEG_PARM(idm,jdm)
      CALL BELEG_VARS(read_rest,idm,jdm,kdm,nbdy,bgc_dp,bgc_rho,omask,
     .                plon,plat)
c
c --- Initialise dust input, n-deposition and river input
c
      CALL ini_fedep(idm,jdm,omask)

      CALL ini_ndep(idm,jdm)

      CALL INI_RIVERINPT(idm,jdm,omask)
#ifdef DMSPH
      CALL GET_PI_PH(idm,jdm,kdm,omask,path)
#endif
c
c --- Read restart fields from restart file if requested, otherwise 
c     (at first start-up) copy ocetra and sediment arrays (which are
c     initialised in BELEG_VARS) to both timelevels of their respective
c     two-time-level counterpart
c
      IF(read_rest.eq.1) THEN
         CALL AUFR_BGC(idm,jdm,kdm,ntr,ntrbgc,itrbgc,trc,
     .                 nyear,nmonth,nday,omask,rstfnm_ocn)
      ELSE
         trc(1:idm,1:jdm,1:kdm,      itrbgc:itrbgc+ntrbgc-1) =
     .     ocetra(:,:,:,:)
         trc(1:idm,1:jdm,kdm+1:2*kdm,itrbgc:itrbgc+ntrbgc-1) =
     .     ocetra(:,:,:,:)
#ifndef sedbypass
         sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
         sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:) 
         powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
         powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:) 
         burial2(:,:,1,:)         = burial(:,:,:)
         burial2(:,:,2,:)         = burial(:,:,:) 
#endif
#if defined(BOXATM)
         atm2(:,:,1,:)            = atm(:,:,:)
         atm2(:,:,2,:)            = atm(:,:,:)
#endif
      ENDIF
c
      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      WRITE(io_stdo_bgc,*)'********************************************'
      write(io_stdo_bgc,*) 'iHAMOCC: finished initialisation'
      write(io_stdo_bgc,*)
      endif

      return
c******************************************************************************
      end subroutine hamocc_init
