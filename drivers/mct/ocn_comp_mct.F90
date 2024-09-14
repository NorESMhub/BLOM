! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen, Alok Kumar Gupta, Ping-Gin Chiu
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

module ocn_comp_mct

   ! -------------------------------------------------------------------
   ! BLOM interface module for the cesm cpl7 mct system
   ! -------------------------------------------------------------------

   use mct_mod
   use esmf,             only: ESMF_Clock
   use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
   use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata, &
                               seq_infodata_putdata, seq_infodata_start_type_cont, &
                               seq_infodata_start_type_brnch, seq_infodata_start_type_start
   use seq_flds_mod
   use seq_timemgr_mod,  only: seq_timemgr_EClockGetData, &
                               seq_timemgr_RestartAlarmIsOn, &
                               seq_timemgr_EClockDateInSync,seq_timemgr_pauseAlarmIsOn
   use seq_comm_mct,     only: seq_comm_suffix, seq_comm_inst, seq_comm_name
   use shr_file_mod,     only: shr_file_getUnit, shr_file_setIO, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
                               shr_file_setLogUnit, shr_file_setLogLevel, &
                               shr_file_freeUnit
   use shr_cal_mod,      only: shr_cal_date2ymd
   use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
   use shr_const_mod,    only: SHR_CONST_REARTH
   use perf_mod,         only: t_startf, t_stopf

   ! BLOM modules
   use dimensions,       only: idm, jdm, nreg, itdm, jtdm
   use mod_types,        only: r8
   use mod_config,       only: inst_index, inst_name, inst_suffix, resume_flag
   use mod_time,         only: blom_time, nstep, baclin, delt1, dlt
   use mod_cesm,         only: runid_cesm, runtyp_cesm, ocn_cpl_dt_cesm
   use mod_xc,           only: mnproc, mpicom_external, xctilr, lp, nbdy, &
                               ii, jj, kk, i0, j0, nproc, jpr, cplmsk, halo_ps, &
                               ifu, ilu, isv, ifv, ilv, isp, ifp, ilp, isu
   use mod_blom_init,    only: blom_init
   use mod_restart,      only: restart_write, restart_read
   use mod_blom_step,    only: blom_step
   use mod_fill_global,  only: fill_global
   use mod_forcing,      only: sprfac, prfac, flxco2, flxdms, flxbrf, flxn2o, flxnh3
   use mod_constants,    only: L_mks2cgs
   use mod_grid,         only: scp2, plon, plat, scuy, scvx, scuxi, scvyi
   use mod_state,        only: u, v, temp, saln, pbu, pbv, ubflxs, vbflxs, sealv
   use mod_cesm,         only: frzpot
   use blom_cpl_indices

   implicit none
   private

   public  :: ocn_init_mct, ocn_run_mct, ocn_final_mct
   private :: ocn_SetGSMap_mct
   private :: domain_mct
   private :: getprecipfact_mct
   private :: sumsbuff_mct

   integer, dimension(:), allocatable ::  &
      perm              ! Permutation array to reorder points

   real(r8), dimension(:,:,:), allocatable :: &
      sbuff             ! Accumulated sum of send buffer quantities for
                        ! averaging before being sent
   integer :: &
      lsize, &          ! Size of attribute vector
      jjcpl, &          ! y-dimension of local ocean domain to send/receive
                        ! fields to coupler
      nsend, &          ! Number of fields to be sent to coupler
      nrecv             ! Number of fields to be received from coupler

   real(r8) :: &
      tlast_coupled, &  ! Time since last coupling
      precip_fact       ! Correction factor for precipitation and runoff

   logical :: &
      lsend_precip_fact ! Flag for sending precipitation/runoff factor

 contains

   subroutine ocn_init_mct(EClock, cdata_o, x2o_o, o2x_o, NLFilename)

      ! Input/output arguments
      type(ESMF_Clock)            , intent(inout)    :: EClock
      type(seq_cdata)             , intent(inout) :: cdata_o
      type(mct_aVect)             , intent(inout) :: x2o_o, o2x_o
      character (len=*), optional , intent(in)    :: NLFilename ! Namelist filename

      ! Local variables
      type(mct_gsMap), pointer :: gsMap_ocn
      type(mct_gGrid), pointer :: dom_ocn
      type(seq_infodata_type), pointer :: infodata   ! Input init object
      integer :: OCNID, mpicom_ocn, shrlogunit, shrloglev, &
                 start_ymd, start_tod, start_year, start_day, start_month
      logical :: exists
      character(len=32) :: starttype

      ! Set cdata pointers
      call seq_cdata_setptrs(cdata_o, ID = OCNID, mpicom = mpicom_ocn, &
                             gsMap = gsMap_ocn, dom = dom_ocn, &
                             infodata = infodata)

      ! Set communicator to be used by blom
      mpicom_external = mpicom_ocn

      ! Get multiple instance data
      inst_name   = seq_comm_name(OCNID)
      inst_index  = seq_comm_inst(OCNID)
      inst_suffix = seq_comm_suffix(OCNID)

      ! ----------------------------------------------------------------
      ! Initialize the model run
      ! ----------------------------------------------------------------

      call blom_cpl_indices_set()

      call seq_infodata_GetData( infodata, case_name = runid_cesm )

      call seq_infodata_GetData( infodata, start_type = starttype)

      if     (trim(starttype) == trim(seq_infodata_start_type_start)) then
         runtyp_cesm = "initial"
      elseif (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
         runtyp_cesm = "continue"
      elseif (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
         runtyp_cesm = "branch"
      else
         write (lp,*) 'ocn_comp_mct ERROR: unknown starttype'
         call shr_sys_flush(lp)
         call shr_sys_abort()
      endif

      !-----------------------------------------------------------------
      ! Get coupling time interval
      !-----------------------------------------------------------------

      call seq_timemgr_EClockGetData(EClock, dtime = ocn_cpl_dt_cesm)

      ! ----------------------------------------------------------------
      ! Initialize blom
      ! ----------------------------------------------------------------

      call t_startf('blom_init')
      call blom_init
      call t_stopf('blom_init')

      ! ----------------------------------------------------------------
      ! Reset shr logging to my log file
      ! ----------------------------------------------------------------

      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (lp)

      call shr_sys_flush(lp)

      ! ----------------------------------------------------------------
      ! Check for consistency of BLOM calender information and EClock
      ! ----------------------------------------------------------------

      ! This must be completed!

      if (runtyp_cesm == 'initial') then
         call seq_timemgr_EClockGetData(EClock, &
                                        start_ymd = start_ymd, &
                                        start_tod = start_tod)
         call shr_cal_date2ymd(start_ymd, start_year, start_month, start_day)
         if (mnproc == 1) then
            write (lp,'(a,i8,a2,i5,a2,i4.4,a1,i2.2,a1,i2.2)') &
               ' cesm initial date:           ',start_ymd,': ',start_tod,': ', &
               start_year,'.',start_month,'.',start_day
            call shr_sys_flush(lp)
         endif
      endif

      ! ----------------------------------------------------------------
      ! Initialize MCT attribute vectors and indices
      ! ----------------------------------------------------------------

      call t_startf ('blom_mct_init')

      ! Initialize ocn gsMap

      call ocn_SetGSMap_mct(mpicom_ocn, OCNID, gsMap_ocn)

      ! Initialize mct ocn domain (needs ocn initialization info)

      if (mnproc == 1) then
         write (lp, *) 'blom: ocn_init_mct: lsize', lsize
      endif

      call domain_mct(gsMap_ocn, dom_ocn, lsize, perm, jjcpl)

      ! Inialize mct attribute vectors

      call mct_aVect_init(x2o_o, rList = seq_flds_x2o_fields, lsize = lsize)
      call mct_aVect_zero(x2o_o)

      call mct_aVect_init(o2x_o, rList = seq_flds_o2x_fields, lsize = lsize)
      call mct_aVect_zero(o2x_o)

      nsend = mct_avect_nRattr(o2x_o)
      nrecv = mct_avect_nRattr(x2o_o)

      !-----------------------------------------------------------------
      ! Send intial state to driver
      !-----------------------------------------------------------------

      call getprecipfact_mct(lsend_precip_fact, precip_fact)
      if ( lsend_precip_fact )  then
         call seq_infodata_PutData( infodata, precip_fact=precip_fact)
      endif
      allocate(sbuff(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nsend))
      tlast_coupled = 0._r8
      call sumsbuff_mct(nsend, sbuff, tlast_coupled)
      call export_mct(o2x_o, lsize, perm, jjcpl, nsend, sbuff, tlast_coupled)
      if (nreg == 2) then
         call seq_infodata_PutData(infodata, ocn_prognostic = .true., &
                                   ocnrof_prognostic = .true., &
                                   ocn_nx = itdm , ocn_ny = jtdm-1)
      else
         call seq_infodata_PutData(infodata, ocn_prognostic = .true., &
                                   ocnrof_prognostic = .true., &
                                   ocn_nx = itdm , ocn_ny = jtdm)
      endif

      call t_stopf('blom_mct_init')


      if (mnproc == 1) then
        write (lp, *) 'blom: completed initialization!'
      endif

      !-----------------------------------------------------------------
      ! Reset shr logging to original values
      !-----------------------------------------------------------------

      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

      call shr_sys_flush(lp)

   end subroutine ocn_init_mct


   subroutine ocn_run_mct(EClock, cdata_o, x2o_o, o2x_o)

      ! Input/output arguments

      type(ESMF_Clock), intent(inout)    :: EClock
      type(seq_cdata) , intent(inout) :: cdata_o
      type(mct_aVect) , intent(inout) :: x2o_o
      type(mct_aVect) , intent(inout) :: o2x_o

      ! Local variables
      type(seq_infodata_type), pointer :: infodata   ! Input init object
      integer :: shrlogunit, shrloglev, ymd, tod, ymd_sync, tod_sync

      ! ----------------------------------------------------------------
      ! Reset shr logging to my log file
      ! ----------------------------------------------------------------

      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (lp)

      call seq_cdata_setptrs(cdata_o, infodata=infodata)

      if (resume_flag) then
          if (mnproc == 1) then
             call blom_time(ymd, tod)
             write(lp,*)'Resume from restart: ymd=',ymd,' tod= ',tod
          endif
         call restart_read  !! resume_flag is applied
         resume_flag = .false.
      end if
      !-----------------------------------------------------------------
      ! Advance the model in time over a coupling interval
      !-----------------------------------------------------------------

      blom_loop: do

         if (nint(tlast_coupled) == 0) then
            ! Obtain import state from driver
            call import_mct(x2o_o, lsize, perm, jjcpl)
         endif

         ! Advance the model a time step
         call blom_step

         ! Add fields to send buffer sums
         call sumsbuff_mct(nsend, sbuff, tlast_coupled)

         if (nint(ocn_cpl_dt_cesm-tlast_coupled) == 0) then
            ! Return export state to driver and exit integration loop
            call export_mct(o2x_o, lsize, perm, jjcpl, nsend, sbuff, &
                            tlast_coupled)
            exit blom_loop
         endif

         if (mnproc == 1) then
            call shr_sys_flush(lp)
         endif

      enddo blom_loop

      call getprecipfact_mct(lsend_precip_fact, precip_fact)
      if ( lsend_precip_fact ) then
         call seq_infodata_PutData( infodata, precip_fact=precip_fact)
      endif

      !-----------------------------------------------------------------
      ! if requested, write restart file
      !-----------------------------------------------------------------

      if (seq_timemgr_RestartAlarmIsOn(EClock) .or. &
          seq_timemgr_pauseAlarmIsOn(EClock)) then
         call restart_write
      endif
      if (seq_timemgr_pauseAlarmIsOn(EClock)) resume_flag = .true.

      !-----------------------------------------------------------------
      ! check that internal clock is in sync with master clock
      !-----------------------------------------------------------------

      if (mnproc == 1) then
         call blom_time(ymd, tod)
         if (.not. seq_timemgr_EClockDateInSync(EClock, ymd, tod )) then
            call seq_timemgr_EClockGetData(EClock, curr_ymd=ymd_sync, &
               curr_tod=tod_sync )
            write(lp,*)' blom ymd=',ymd     ,'  blom tod= ',tod
            write(lp,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
            call shr_sys_abort( 'ocn_run_mct'// &
               ":: Internal blom clock not in sync with Sync Clock")
         endif
      endif

      !-----------------------------------------------------------------
      ! Reset shr logging to original values
      !-----------------------------------------------------------------

      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

   end subroutine ocn_run_mct


   subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)

      type(ESMF_Clock)            , intent(inout)    :: EClock
      type(seq_cdata)             , intent(inout) :: cdata_o
      type(mct_aVect)             , intent(inout) :: x2o_o
      type(mct_aVect)             , intent(inout) :: o2x_o

      deallocate(perm)
      deallocate(sbuff)

   end subroutine ocn_final_mct


   subroutine ocn_SetGSMap_mct(mpicom_ocn, OCNID, gsMap_ocn)

      ! Input/output arguments

      integer        , intent(in)    :: mpicom_ocn
      integer        , intent(in)    :: OCNID
      type(mct_gsMap), intent(inout) :: gsMap_ocn

      ! Local variables

      integer, allocatable :: gindex(:)
      integer :: i, j, n, gsize

      ! ----------------------------------------------------------------
      ! Build the BLOM grid numbering for MCT
      ! NOTE:  Numbering scheme is: West to East and South to North
      ! starting at south pole.  Should be the same as what's used
      ! in SCRIP
      ! ----------------------------------------------------------------

      if (nreg == 2 .and. nproc == jpr) then
         jjcpl = jj - 1
      else
         jjcpl = jj
      endif

      lsize = ii*jjcpl

      if (nreg == 2) then
         gsize = itdm*(jtdm-1)
      else
         gsize = itdm*jtdm
      endif

      allocate(gindex(lsize))

      n = 0
      do j = 1, jjcpl
         do i = 1, ii
            n = n + 1
            gindex(n) = (j0 + j - 1)*itdm + i0 + i
         enddo
      enddo

      ! ----------------------------------------------------------------
      ! reorder gindex to be in ascending order.
      !  initialize a permutation array and sort gindex in-place
      ! ----------------------------------------------------------------

      allocate(perm(lsize))

      call mct_indexset(perm)
      call mct_indexsort(lsize, perm, gindex)
      call mct_permute(gindex, perm, lsize)
      call mct_gsMap_init(gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize)

      deallocate(gindex)

   end subroutine ocn_SetGSMap_mct

   subroutine domain_mct(gsMap_ocn, dom_ocn, lsize, perm, jjcpl)

     ! Arguments
     type(mct_gsMap)          , intent(in)    :: gsMap_ocn
     type(mct_ggrid)          , intent(inout) :: dom_ocn
     integer                  , intent(in)    :: lsize
     integer, dimension(lsize), intent(in)    :: perm
     integer                  , intent(in)    :: jjcpl

     ! Local variables
     integer, pointer :: idata(:)
     real(r8), pointer :: rdata(:)
     integer i, j, n
     real(r8) :: radius

     ! ----------------------------------------------------------------
     ! Initialize mct domain type
     ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean),
     ! 0 (non-ocean)
     ! ----------------------------------------------------------------

     call mct_gGrid_init(GGrid = dom_ocn, &
                         CoordChars = trim(seq_flds_dom_coord), &
                         OtherChars = trim(seq_flds_dom_other), &
                         lsize = lsize)
     allocate(rdata(lsize))

     ! ----------------------------------------------------------------
     ! Determine global gridpoint number attribute, GlobGridNum, which
     ! is set automatically by MCT
     ! ----------------------------------------------------------------

     call mct_gsMap_orderedPoints(gsMap_ocn, mnproc - 1, idata)
     call mct_gGrid_importIAttr(dom_ocn, 'GlobGridNum', idata, lsize)

     ! ----------------------------------------------------------------
     ! Determine domain (numbering scheme is: West to East and South to
     ! North to South pole)
     ! Initialize attribute vector with special value
     ! ----------------------------------------------------------------

     rdata(:) = -9999.0_r8
     call mct_gGrid_importRAttr(dom_ocn, "lat"  , rdata, lsize)
     call mct_gGrid_importRAttr(dom_ocn, "lon"  , rdata, lsize)
     call mct_gGrid_importRAttr(dom_ocn, "area" , rdata, lsize)
     call mct_gGrid_importRAttr(dom_ocn, "aream", rdata, lsize)
     rdata(:) = 0.0_r8
     call mct_gGrid_importRAttr(dom_ocn, "mask", rdata, lsize)
     call mct_gGrid_importRAttr(dom_ocn, "frac", rdata, lsize)

     ! ----------------------------------------------------------------
     ! Fill in correct values for domain components
     ! ----------------------------------------------------------------

     ! A correction for north pole mapping of velocity fields in the
     ! coupler requires longitudes in the range [0, 360) degrees to
     ! work.
     n = 0
     do j = 1, jjcpl
       do i = 1, ii
         n = n + 1
         rdata(n) = modulo(plon(i,j), 360._r8)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_ocn, "lon", rdata, lsize)

     n = 0
     do j = 1, jjcpl
       do i = 1, ii
         n = n + 1
         rdata(n) = plat(i,j)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_ocn, "lat", rdata, lsize)

     radius = SHR_CONST_REARTH*L_mks2cgs ! Earth's radius in cm

     n = 0
     do j = 1, jjcpl
       do i = 1, ii
         n = n + 1
         rdata(n) = scp2(i,j)/(radius*radius)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_ocn, "area", rdata, lsize)

     n = 0
     do j = 1, jjcpl
       do i = 1, ii
         n = n + 1
         rdata(n) = real(cplmsk(i,j), kind = r8)
       enddo
     enddo
     call mct_gGrid_importRattr(dom_ocn, "mask", rdata, lsize)
     call mct_gGrid_importRattr(dom_ocn, "frac", rdata, lsize)

     !-----------------------------------------------------------------
     ! Permute dom_ocn to have ascending order
     !-----------------------------------------------------------------

     call mct_gGrid_permute(dom_ocn, perm)
     deallocate(rdata)

   end subroutine domain_mct

   subroutine getprecipfact_mct(lsend_precip_fact, precip_fact)
     logical, intent(out)  :: lsend_precip_fact
     real(r8), intent(out) :: precip_fact

     lsend_precip_fact = sprfac
     precip_fact = prfac

   end subroutine getprecipfact_mct

   subroutine sumsbuff_mct(nsend, sbuff, tlast_coupled)

     ! Arguments
     integer, intent(in) :: nsend
     real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nsend), intent(inout) :: sbuff
     real(r8), intent(inout) :: tlast_coupled

     ! Local variables
     integer i, j, l, k, m, n, mm, nn, k1m, k1n

     !-----------------------------------------------------------------
     ! Set send buffer to zero if this is the first call after a
     ! coupling interval
     !-----------------------------------------------------------------

     if (tlast_coupled == 0._r8) then
       do k = 1, nsend
         do j = 1-nbdy, jdm+nbdy
           do i = 1-nbdy, idm+nbdy
             sbuff(i,j,k) = 0._r8
           enddo
         enddo
       enddo
     endif

     !-----------------------------------------------------------------
     ! Accumulate fields in send buffer
     !-----------------------------------------------------------------

     m   = mod(nstep+1,2)+1
     n   = mod(nstep  ,2)+1
     mm  = (m-1)*kk
     nn  = (n-1)*kk
     k1m = 1+mm
     k1n = 1+nn

     call xctilr(sealv, 1,1, 1,1, halo_ps)

     do j = 1, jj
       do l = 1, isu(j)
         do i = max(1,ifu(j,l)), min(ii,ilu(j,l))
           sbuff(i,j,index_o2x_So_u) = sbuff(i,j,index_o2x_So_u) &
                + ( u(i,j,k1n)+ (ubflxs(i,j,m) + ubflxs(i,j,n))*dlt &
                /(pbu(i,j,n)*scuy(i,j)*delt1))*baclin

           sbuff(i,j,index_o2x_So_dhdx) = sbuff(i,j,index_o2x_So_dhdx) &
                + (sealv(i,j)-sealv(i-1,j))*scuxi(i,j)*baclin
         enddo
       enddo
     enddo

     do j = 1, jj
       do l = 1, isv(j)
         do i = max(1,ifv(j,l)), min(ii,ilv(j,l))
           sbuff(i,j,index_o2x_So_v) = sbuff(i,j,index_o2x_So_v) &
                + ( v(i,j,k1n) + (vbflxs(i,j,m) + vbflxs(i,j,n))*dlt &
                / (pbv(i,j,n)*scvx(i,j)*delt1))*baclin
           sbuff(i,j,index_o2x_So_dhdy) = sbuff(i,j,index_o2x_So_dhdy) &
                + (sealv(i,j)-sealv(i,j-1))*scvyi(i,j)*baclin
         enddo
       enddo
     enddo

     do j = 1, jj
       do l = 1, isp(j)
         do i = max(1,ifp(j,l)), min(ii,ilp(j,l))
           sbuff(i,j,index_o2x_So_t) = sbuff(i,j,index_o2x_So_t) + temp(i,j,k1n)*baclin
           sbuff(i,j,index_o2x_So_s) = sbuff(i,j,index_o2x_So_s) + saln(i,j,k1n)*baclin
           sbuff(i,j,index_o2x_Fioo_q) = sbuff(i,j,index_o2x_Fioo_q) + frzpot(i,j)
         enddo
       enddo
     enddo

     if (index_o2x_Faoo_fco2_ocn > 0) then
       do j = 1, jj
         do l = 1, isp(j)
           do i = max(1,ifp(j,l)), min(ii,ilp(j,l))
             sbuff(i,j,index_o2x_Faoo_fco2_ocn) = sbuff(i,j,index_o2x_Faoo_fco2_ocn) &
                  + flxco2(i,j)*baclin
           enddo
         enddo
       enddo
     endif

     if (index_o2x_Faoo_fdms_ocn > 0) then
       do j = 1, jj
         do l = 1, isp(j)
           do i = max(1,ifp(j,l)), min(ii,ilp(j,l))
             sbuff(i,j,index_o2x_Faoo_fdms_ocn) = sbuff(i,j,index_o2x_Faoo_fdms_ocn) &
                  + flxdms(i,j)*baclin
           enddo
         enddo
       enddo
     endif

     if (index_o2x_Faoo_fbrf_ocn > 0) then
       do j = 1, jj
         do l = 1, isp(j)
           do i = max(1,ifp(j,l)), min(ii,ilp(j,l))
             sbuff(i,j,index_o2x_Faoo_fbrf_ocn) = sbuff(i,j,index_o2x_Faoo_fbrf_ocn) &
                  + flxbrf(i,j)*baclin
           enddo
         enddo
       enddo
     endif

      if (index_o2x_Faoo_fn2o_ocn > 0) then
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1,ifp(j,l)), min(ii,ilp(j,l))
               sbuff(i,j,index_o2x_Faoo_fn2o_ocn) = sbuff(i,j,index_o2x_Faoo_fn2o_ocn) &
                    + flxn2o(i,j)*baclin
            enddo
            enddo
         enddo
      endif

      if (index_o2x_Faoo_fnh3_ocn > 0) then
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1,ifp(j,l)), min(ii,ilp(j,l))
               sbuff(i,j,index_o2x_Faoo_fnh3_ocn) = sbuff(i,j,index_o2x_Faoo_fnh3_ocn) &
                    + flxnh3(i,j)*baclin
            enddo
            enddo
         enddo
      endif

     !-----------------------------------------------------------------
     ! Increment time since last coupling
     !-----------------------------------------------------------------

     tlast_coupled = tlast_coupled + baclin

   end subroutine sumsbuff_mct

 end module ocn_comp_mct
