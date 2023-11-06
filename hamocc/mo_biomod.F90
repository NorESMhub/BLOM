! Copyright (C) 2001  S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, I. Kriest
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.


      MODULE mo_biomod
!******************************************************************************
!
! MODULE mo_biomod - Variables for marine biology.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!     
!     I. Kriest, GEOMAR, 11.08.2016
!     - included T-dependence of cyanobacteria growth
!     - modified stoichiometry for denitrification
! 
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved accumulation of all output fields to seperate SUBROUTINE,
!       new global fields for output defined here
!
!     Purpose
!     -------
!     - declaration and memory allocation.
!
!     Description:
!     ------------
!     Public routines and variable of this MODULE:
!
!     -SUBROUTINE alloc_mem_biomod
!        Allocate memory for biomod variables
!
!
!******************************************************************************
      implicit none

      real, dimension (:,:),   allocatable :: strahl
      real, dimension (:,:),   allocatable :: expoor
      real, dimension (:,:),   allocatable :: expoca
      real, dimension (:,:),   allocatable :: exposi
      real, dimension (:,:),   allocatable :: intphosy
      real, dimension (:,:),   allocatable :: intdnit
      real, dimension (:,:),   allocatable :: intnfix
      real, dimension (:,:),   allocatable :: intdmsprod
      real, dimension (:,:),   allocatable :: intdms_bac
      real, dimension (:,:),   allocatable :: intdms_uv
      real, dimension (:,:),   allocatable :: carflx0100
      real, dimension (:,:),   allocatable :: carflx0500
      real, dimension (:,:),   allocatable :: carflx1000
      real, dimension (:,:),   allocatable :: carflx2000
      real, dimension (:,:),   allocatable :: carflx4000
      real, dimension (:,:),   allocatable :: carflx_bot
      real, dimension (:,:),   allocatable :: bsiflx0100
      real, dimension (:,:),   allocatable :: bsiflx0500
      real, dimension (:,:),   allocatable :: bsiflx1000
      real, dimension (:,:),   allocatable :: bsiflx2000
      real, dimension (:,:),   allocatable :: bsiflx4000
      real, dimension (:,:),   allocatable :: bsiflx_bot
      real, dimension (:,:),   allocatable :: calflx0100
      real, dimension (:,:),   allocatable :: calflx0500
      real, dimension (:,:),   allocatable :: calflx1000
      real, dimension (:,:),   allocatable :: calflx2000
      real, dimension (:,:),   allocatable :: calflx4000
      real, dimension (:,:),   allocatable :: calflx_bot
      real, dimension (:,:,:), allocatable :: phosy3d

      ! Variables for interactive phytoplanktion absorption (use_FB_BGC_OCE=.true.)
      real, dimension (:,:,:), allocatable :: abs_oce
	  
      ! Variables for aggregation scheme (use_AGG=.true.)
      real, dimension (:,:,:), allocatable :: wmass
      real, dimension (:,:,:), allocatable :: wnumb
      real, dimension (:,:,:), allocatable :: eps3d
      real, dimension (:,:,:), allocatable :: asize3d

      ! Variables for bromoform scheme (use_BROMO=.true.)
      real, dimension (:,:),   allocatable :: int_chbr3_prod
      real, dimension (:,:),   allocatable :: int_chbr3_uv

      real :: growth_co2,bifr13_perm
	  

      CONTAINS

      SUBROUTINE ALLOC_MEM_BIOMOD(kpie,kpje,kpke)
!******************************************************************************
! ALLOC_MEM_BIOMOD - Allocate variables in this MODULE
!******************************************************************************
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc
      use mo_control_bgc, only: use_FB_BGC_OCE,use_AGG,use_BROMO

      integer, intent(in) :: kpie,kpje,kpke
      integer             :: errstat
      

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'Memory allocation for marine biology module :'
      write(io_stdo_bgc,*)' '
      endif


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable strahl ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (strahl(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory strahl'
      strahl(:,:) = 0.0

      if (use_FB_BGC_OCE ) then
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable abs_oce'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (abs_oce(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory abs_oce'
         abs_oce(:,:,:) = 0.0
      endif

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable expoor ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (expoor(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory expoor'
      expoor(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable expoca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (expoca(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory expoca'
      expoca(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable exposi ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (exposi(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory exposi'
      exposi(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intphosy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (intphosy(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory intphosy'
      intphosy(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intdnit ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (intdnit(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory intdnit'
      intdnit(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intnfix ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (intnfix(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory intnfix'
      intnfix(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intdmsprod, intdms_bac, intdms_uv ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (intdmsprod(kpie,kpje),stat=errstat)
      ALLOCATE (intdms_bac(kpie,kpje),stat=errstat)
      ALLOCATE (intdms_uv(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory intdmsprod, intdms_bac, intdms_uv'
      intdmsprod(:,:) = 0.0
      intdms_bac(:,:) = 0.0
      intdms_uv(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable carflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (carflx0100(kpie,kpje),stat=errstat)
      ALLOCATE (carflx0500(kpie,kpje),stat=errstat)
      ALLOCATE (carflx1000(kpie,kpje),stat=errstat)
      ALLOCATE (carflx2000(kpie,kpje),stat=errstat)
      ALLOCATE (carflx4000(kpie,kpje),stat=errstat)
      ALLOCATE (carflx_bot(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory carflx*'
      carflx0100(:,:) = 0.0
      carflx0500(:,:) = 0.0
      carflx1000(:,:) = 0.0
      carflx2000(:,:) = 0.0
      carflx4000(:,:) = 0.0
      carflx_bot(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable bsiflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (bsiflx0100(kpie,kpje),stat=errstat)
      ALLOCATE (bsiflx0500(kpie,kpje),stat=errstat)
      ALLOCATE (bsiflx1000(kpie,kpje),stat=errstat)
      ALLOCATE (bsiflx2000(kpie,kpje),stat=errstat)
      ALLOCATE (bsiflx4000(kpie,kpje),stat=errstat)
      ALLOCATE (bsiflx_bot(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory bsiflx*'
      bsiflx0100(:,:) = 0.0
      bsiflx0500(:,:) = 0.0
      bsiflx1000(:,:) = 0.0
      bsiflx2000(:,:) = 0.0
      bsiflx4000(:,:) = 0.0
      bsiflx_bot(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable calflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (calflx0100(kpie,kpje),stat=errstat)
      ALLOCATE (calflx0500(kpie,kpje),stat=errstat)
      ALLOCATE (calflx1000(kpie,kpje),stat=errstat)
      ALLOCATE (calflx2000(kpie,kpje),stat=errstat)
      ALLOCATE (calflx4000(kpie,kpje),stat=errstat)
      ALLOCATE (calflx_bot(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory bsiflx*'
      calflx0100(:,:) = 0.0
      calflx0500(:,:) = 0.0
      calflx1000(:,:) = 0.0
      calflx2000(:,:) = 0.0
      calflx4000(:,:) = 0.0
      calflx_bot(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable phosy3d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      ALLOCATE (phosy3d(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory phosy3d'
      phosy3d(:,:,:) = 0.0
         
      if (use_AGG) then
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable wmass ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (wmass(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory eps3d'
         wmass(:,:,:) = 0.0

         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable wnumb ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (wnumb(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory eps3d'
         wnumb(:,:,:) = 0.0
         
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable eps3d ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (eps3d(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory eps3d'
         eps3d(:,:,:) = 0.0

         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable asize3d ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (asize3d(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory asize3d'
         asize3d(:,:,:) = 0.0
      endif

      if (use_BROMO) then
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable int_chbr3_prod, int_chbr3_uv ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         endif

         ALLOCATE (int_chbr3_prod(kpie,kpje),stat=errstat)
         ALLOCATE (int_chbr3_uv(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory int_chbr3_prod, int_chbr3_uv'
         int_chbr3_prod(:,:) = 0.0
         int_chbr3_uv(:,:) = 0.0
      endif

!******************************************************************************
      END SUBROUTINE ALLOC_MEM_BIOMOD

      END MODULE mo_biomod
