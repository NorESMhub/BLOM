! ------------------------------------------------------------------------------
! Copyright (C) 2022 Mats Bentsen
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

module ocn_comp_nuopc
! ------------------------------------------------------------------------------
! This module contains the NUOPC cap for BLOM.
! ------------------------------------------------------------------------------

   use ESMF ! TODO MOM6 uses "only" statements, while POP and CICE omits this.
   use NUOPC, only: NUOPC_CompDerive, NUOPC_CompSetEntryPoint, &
                    NUOPC_CompSpecialize, NUOPC_CompFilterPhaseMap, &
                    NUOPC_IsUpdated, NUOPC_IsAtTime, NUOPC_CompAttributeGet, &
                    NUOPC_Advertise, NUOPC_SetAttribute, &
                    NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, &
                    NUOPC_IsConnected, NUOPC_Realize
   use NUOPC_Model, only: NUOPC_ModelGet, SetVM, &
                          model_routine_SS           => SetServices, &
                          model_label_Advance        => label_Advance, &
                          model_label_DataInitialize => label_DataInitialize, &
                          model_label_SetRunClock    => label_SetRunClock, &
                          model_label_Finalize       => label_Finalize
   use nuopc_shr_methods, only : ChkErr, set_component_logging, &
                                 get_component_instance, state_setscalar, &
                                 alarmInit
   use shr_file_mod, only: shr_file_getUnit, shr_file_getLogUnit, &
                           shr_file_setLogUnit
   use shr_cal_mod, only : shr_cal_ymd2date
   use mod_nuopc_methods, only: fldlist_type, fldsMax, tlast_coupled, &
                                blom_logwrite, blom_getgindex, blom_checkmesh, &
                                blom_setareacor, blom_getglobdim, &
                                blom_getprecipfact, blom_accflds, &
                                blom_importflds, blom_exportflds, &
                                blom_advertise_imports, blom_advertise_exports, &
                                get_flxdms_from_med
   use mod_xc, only: mpicom_external, lp, nfu
   use mod_cesm, only: runid_cesm, runtyp_cesm, ocn_cpl_dt_cesm
   use mod_config, only: inst_index, inst_name, inst_suffix
   use mod_time, only: blom_time
   use mod_forcing, only : compute_flxdms

   implicit none

   private

   integer, parameter :: cslen = 80  ! Short character string length.
   integer, parameter :: cllen = 265 ! Long character string length.
   character(len=*), parameter :: modname = '(ocn_comp_nuopc)'
   character(len=*), parameter :: u_FILE_u = &
      __FILE__

   integer              :: fldsToOcn_num = 0
   integer              :: fldsFrOcn_num = 0
   type(fldlist_type)   :: fldsToOcn(fldsMax)
   type(fldlist_type)   :: fldsFrOcn(fldsMax)

   character(len=cllen) :: flds_scalar_name = ''
   integer              :: flds_scalar_num = 0
   integer              :: flds_scalar_index_nx = 0
   integer              :: flds_scalar_index_ny = 0
   integer              :: flds_scalar_index_precip_factor = 0

   logical              :: ocn2glc_coupling, flds_dms_med

   integer :: dbug = 0
   logical :: profile_memory = .false.

   public :: SetServices, SetVM

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   subroutine fldlist_realize(state, fldlist_num, fldlist, tag, mesh, rc)
   ! ---------------------------------------------------------------------------
   ! Realize list of import or export fields.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_State)  , intent(inout) :: state
      integer           , intent(in)    :: fldlist_num
      type(fldlist_type), intent(in)    :: fldlist(:)
      character(len=*)  , intent(in)    :: tag
      type(ESMF_Mesh)   , intent(in)    :: mesh
      integer           , intent(inout) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(fldlist_realize)'

      ! Local variables.
      integer                :: n
      type(ESMF_DistGrid)    :: DistGrid
      type(ESMF_Grid)        :: grid
      type(ESMF_Field)       :: field
      character(len=128)     :: stdname
      character(ESMF_MAXSTR) :: msg

      rc = ESMF_SUCCESS

      do n = 1, fldlist_num

         stdname = fldlist(n)%stdname

         if (NUOPC_IsConnected(state, fieldName=stdname)) then

            if (stdname == trim(flds_scalar_name)) then

               ! Create the scalar field.
               call ESMF_LogWrite(subname//trim(tag)//" Field = "// &
                                  trim(stdname)//" is connected on root pe", &
                                  ESMF_LOGMSG_INFO)
               DistGrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
               if (ChkErr(rc, __LINE__, u_FILE_u)) return
               grid = ESMF_GridCreate(DistGrid, rc=rc)
               if (ChkErr(rc, __LINE__, u_FILE_u)) return
               field = ESMF_FieldCreate(name=trim(flds_scalar_name), &
                                        grid=grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        ungriddedLBound=(/1/), &
                                        ungriddedUBound=(/flds_scalar_num/), &
                                        gridToFieldMap=(/2/), rc=rc)
               if (ChkErr(rc, __LINE__, u_FILE_u)) return

            else

               ! Create the field
               if (fldlist(n)%ungridded_lbound > 0 .and. &
                   fldlist(n)%ungridded_ubound > 0) then
                  field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
                             name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                             ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                             ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                             gridToFieldMap=(/2/), rc=rc)
                  if (ChkErr(rc, __LINE__, u_FILE_u)) return
                  write(msg,'(a,i4,2x,i4)') &
                     subname//trim(tag)//" Field = "//trim(stdname)// &
                     " is connected using mesh with lbound, ubound = ", &
                     fldlist(n)%ungridded_lbound, fldlist(n)%ungridded_ubound
                  call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
               else
                  field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, &
                                           name=stdname, &
                                           meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                  if (ChkErr(rc, __LINE__, u_FILE_u)) return
                  write(msg,'(a)') &
                     subname//trim(tag)//" Field = "//trim(stdname)// &
                     " is connected using mesh without ungridded dimension"
                  call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
               endif

            endif

            ! Realize connected field.
            call NUOPC_Realize(state, field=field, rc=rc)
            if (ChkErr(rc, __LINE__, u_FILE_u)) return

         else

            if (stdname /= trim(flds_scalar_name)) then

               call ESMF_LogWrite(subname//trim(tag)//" Field = "// &
                                  trim(stdname)// " is not connected", &
                                  ESMF_LOGMSG_INFO)

               ! Remove a not connected field from state.
               call ESMF_StateRemove(state, (/stdname/), rc=rc)
               if (ChkErr(rc, __LINE__, u_FILE_u)) return

            endif

         endif

      enddo

   end subroutine fldlist_realize

   subroutine ocn_import(importState, rc)
   ! ---------------------------------------------------------------------------
   ! Import data from the mediator to ocean.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_State)     :: importState
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(import)'

      ! Local variables.
      type(ESMF_StateItem_Flag) :: itemType
      type(ESMF_Field) :: field
      integer :: n

      rc = ESMF_SUCCESS

      ! Get data pointers for the fields to be imported.
      do n = 1, fldsToOcn_num
         if (fldsToOcn(n)%stdname == trim(flds_scalar_name)) cycle
         call ESMF_StateGet(importState, trim(fldsToOcn(n)%stdname), &
                            itemType, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         if (itemType == ESMF_STATEITEM_NOTFOUND) then
            fldsToOcn(n)%dataptr => null()
         else
            call ESMF_StateGet(importState, trim(fldsToOcn(n)%stdname), &
                               field=field, rc=rc)
            if (ChkErr(rc, __LINE__, u_FILE_u)) return
            call ESMF_FieldGet(field, farrayPtr=fldsToOcn(n)%dataptr, rc=rc)
            if (ChkErr(rc, __LINE__, u_FILE_u)) return
         endif
      enddo

      ! Import fields from mediator to BLOM arrays.
      call blom_importflds(fldsToOcn_num, fldsToOcn)

   end subroutine ocn_import

   subroutine ocn_export(exportState, rc)
   ! ---------------------------------------------------------------------------
   ! Export data from ocean to the mediator.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_State)     :: exportState
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(export)'

      ! Local variables.
      type(ESMF_StateItem_Flag) :: itemType
      type(ESMF_Field) :: field
      real(ESMF_KIND_R8) :: precip_fact
      integer :: n
      logical :: precip_fact_provided

      rc = ESMF_SUCCESS

      ! Get data pointers for the fields to be exported.
      do n = 1, fldsFrOcn_num
         if (fldsFrOcn(n)%stdname == trim(flds_scalar_name)) cycle
         call ESMF_StateGet(exportState, trim(fldsFrOcn(n)%stdname), &
                            itemType, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         if (itemType == ESMF_STATEITEM_NOTFOUND) then
            fldsFrOcn(n)%dataptr => null()
         else
            call ESMF_StateGet(exportState, trim(fldsFrOcn(n)%stdname), &
                               field=field, rc=rc)
            if (ChkErr(rc, __LINE__, u_FILE_u)) return
            call ESMF_FieldGet(field, farrayPtr=fldsFrOcn(n)%dataptr, rc=rc)
            if (ChkErr(rc, __LINE__, u_FILE_u)) return
         endif
      enddo

      ! Export from BLOM arrays to mediator fields.
      call blom_exportflds(fldsFrOcn_num, fldsFrOcn)

      ! Provide precipitation factor.
      call blom_getprecipfact(precip_fact_provided, precip_fact)
      if (precip_fact_provided) then
         call state_setscalar(precip_fact, &
                              flds_scalar_index_precip_factor, exportState, &
                              flds_scalar_name, flds_scalar_num, rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
      else
         call state_setscalar(1._ESMF_KIND_R8, &
                              flds_scalar_index_precip_factor, exportState, &
                              flds_scalar_name, flds_scalar_num, rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
      endif

   end subroutine ocn_export

   subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
   ! ---------------------------------------------------------------------------
   ! Set which version of the Initialize Phase Definition (IPD) to use.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(InitializeP0)'

      ! Local variables.
      logical :: isPresent, isSet
      character(len=cslen) :: cvalue

      ! Switch to IPDv01 by filtering all other PhaseMap entries
      call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
                                    acceptStringList=(/"IPDv01p"/), rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      profile_memory = .false.
      call NUOPC_CompAttributeGet(gcomp, name="ProfileMemory", value=cvalue, &
                                  isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      if (isPresent .and. isSet) profile_memory = (trim(cvalue) == "true")
      write(cvalue,*) profile_memory
      call ESMF_LogWrite(subname//': ProfileMemory = '//trim(cvalue), &
                         ESMF_LOGMSG_INFO)

   end subroutine InitializeP0

   subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
   ! ---------------------------------------------------------------------------
   ! Called by NUOPC to advertise import and export fields. "Advertise" simply
   ! means that the standard names of all import and export fields are supplied.
   ! The NUOPC layer uses these to match fields between components in the
   ! coupled system.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(InitializeAdvertise)'

      ! local variables
      type(ESMF_VM) :: vm
      type(ESMF_TimeInterval) :: timeStep
      integer :: localPet, nthrds, shrlogunit, n
      character(len=cslen) :: starttype, stdname, cvalue, cname
      character(len=cllen) :: msg
      logical :: isPresent, isSet
      logical :: flds_co2a, flds_co2c

      ! Get debug flag.
      call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=cvalue, &
                                  isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      if (isPresent .and. isSet) read(cvalue,*) dbug
      write(cvalue,*) dbug
      call ESMF_LogWrite(subname//': dbug = '//trim(cvalue), ESMF_LOGMSG_INFO)

      ! Get local MPI communicator and Persistent Execution Thread (PET).
      call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_VMGet(vm, mpiCommunicator=mpicom_external, localPet=localPet, &
                      rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! OpenMP threads
      call ESMF_VMGet(vm, pet=localPet, peCount=nthrds, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      if (nthrds == 1) then
         call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         read(cvalue,*) nthrds
      endif
!$    call omp_set_num_threads(nthrds)

      ! Reset shr logging to components log file.
      call set_component_logging(gcomp, localPet==0, lp, shrlogunit, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! Get generic file unit for master task.
      if (localPet == 0) then
         nfu = shr_file_getUnit()
      else
         nfu = -1
      endif

      ! Get case name.
      call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) runid_cesm

      ! Get start type.
      call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) starttype
      if (trim(starttype) == trim('startup')) then
         runtyp_cesm = "initial"
      else if (trim(starttype) == trim('continue') ) then
         runtyp_cesm = "continue"
      else if (trim(starttype) == trim('branch')) then
         runtyp_cesm = "continue"
      else
         call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg=subname//": unknown starttype - "//trim(starttype), &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
         return
      endif

      ! Get multiple instance data.
      call get_component_instance(gcomp, inst_suffix, inst_index, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      inst_name = "OCN"

      ! Get coupling time interval.
      call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_TimeIntervalGet(timeStep, s=ocn_cpl_dt_cesm, rc=rc )
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! Initialize BLOM.
      ! ------------------------------------------------------------------------

      call blom_init

      ! ------------------------------------------------------------------------
      ! Get ScalarField attributes.
      ! ------------------------------------------------------------------------

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, &
                                  isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      if (isPresent .and. isSet) then
         flds_scalar_name = trim(cvalue)
         call ESMF_LogWrite(subname//': flds_scalar_name = '//trim(cvalue), &
                            ESMF_LOGMSG_INFO)
      else
         call ESMF_LogSetError(ESMF_RC_NOT_SET, &
                               msg=subname//": ScalarFieldName is not set", &
                               line=__LINE__, file=__FILE__, rcToReturn=rc)
         return
      endif

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", &
                                  value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) flds_scalar_num
      write(cvalue,*) flds_scalar_num
      call ESMF_LogWrite(subname//': flds_scalar_num = '//trim(cvalue), &
                         ESMF_LOGMSG_INFO)

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", &
                                  value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) flds_scalar_index_nx
      write(cvalue,*) flds_scalar_index_nx
      call ESMF_LogWrite(subname//': flds_scalar_index_nx = '//trim(cvalue), &
                         ESMF_LOGMSG_INFO)

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", &
                                  value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) flds_scalar_index_ny
      write(cvalue,*) flds_scalar_index_ny
      call ESMF_LogWrite(subname//': flds_scalar_index_ny = '//trim(cvalue), &
                         ESMF_LOGMSG_INFO)

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxPrecipFactor", &
                                  value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) flds_scalar_index_precip_factor
      if ( .not. flds_scalar_index_precip_factor > 0 ) then
         call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, &
            msg=subname//": flds_scalar_index_precip_factor must be > 0", &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
         return
      else
         write(cvalue,*) flds_scalar_index_precip_factor
         call ESMF_LogWrite(subname//': flds_scalar_index_precip_factor = '// &
                            trim(cvalue), ESMF_LOGMSG_INFO)
      endif

      ! Determine if dms flux will be computed in mediator and sent to BLOM
      call NUOPC_CompAttributeGet(gcomp, name='flds_dms_med', value=cvalue, &
           isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! Set the logical flag get_flxdms_from_med in mod_nuopc_methods module
      if (.not. isPresent .and. .not. isSet) then
         get_flxdms_from_med = .false.
      else
         read(cvalue,*) flds_dms_med
         call blom_logwrite(subname//': flds_dms_med = '//trim(cvalue))
         if (flds_dms_med) then
            get_flxdms_from_med = .true.
         else
            get_flxdms_from_med = .false.
         end if
      end if

      ! Set the logical flag compute_flxdms in mod_forcing module
      if (get_flxdms_from_med) then
         compute_flxdms = .false.
      else
         compute_flxdms = .true.
      end if
      write(6,*)'DEBUG: compute_flxdms = ',compute_flxdms

      ! Determine if co2 will be imported from mediator
      call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) flds_co2a
      call blom_logwrite(subname//': flds_co2a = '//trim(cvalue))

      call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) flds_co2c
      call blom_logwrite(subname//': flds_co2c = '//trim(cvalue))

      ! Determine if ocn is sending temperature and salinity data to glc
      ! If data is sent to glc will need to determine number of ocean
      ! levels and ocean level indices
      call NUOPC_CompAttributeGet(gcomp, name="ocn2glc_coupling", value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      read(cvalue,*) ocn2glc_coupling
      write(msg,'(a,l1)') subname//': ocn2glc coupling is ', ocn2glc_coupling
      call blom_logwrite(msg)
      if (ocn2glc_coupling) then
         call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg=subname//": ocn2glc coupling not implemented", &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
         return
      endif

      !NOTE: Nitrogen deposition is always sent from atm now (either CAM or DATM)

      ! ------------------------------------------------------------------------
      ! Advertise import fields.
      ! ------------------------------------------------------------------------

      call blom_advertise_imports(flds_scalar_name, fldsToOcn_num, fldsToOcn, &
           flds_co2a, flds_co2c)

      do n = 1,fldsToOcn_num
         call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, &
                              TransferOfferGeomObject='will provide', rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
      enddo

      ! ------------------------------------------------------------------------
      ! Advertise export fields.
      ! ------------------------------------------------------------------------

      call blom_advertise_exports(flds_scalar_name, fldsFrOcn_num, fldsFrOcn)

      do n = 1,fldsFrOcn_num
         call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, &
                              TransferOfferGeomObject='will provide', rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
      enddo

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine InitializeAdvertise

   subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
   ! ---------------------------------------------------------------------------
   ! Called by NUOPC to realize import and export fields. "Realizing" a field
   ! means that its grid has been defined and an ESMF_Field object has been
   ! created and put into the import or export State.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(InitializeRealize)'

      ! Local variables.
      type(ESMF_DistGrid) :: DistGrid
      type(ESMF_Mesh) :: EMesh
      type(ESMF_Array) :: elemMaskArray
      type(ESMF_Field) :: field
      real(ESMF_KIND_R8), dimension(:), pointer :: &
         ownedElemCoords, lonMesh, latMesh, areaMesh
      integer(ESMF_KIND_I4), dimension(:), pointer :: maskMesh(:)
      integer, allocatable, dimension(:) :: gindex
      integer :: n, spatialDim, numOwnedElements, nx_global, ny_global
      character(len=cslen)  :: cvalue

      if (dbug > 5) call ESMF_LogWrite(subname//': called', ESMF_LOGMSG_INFO)

      ! Get the BLOM global index space for the computational domain.
      call blom_getgindex(gindex)

      ! Create DistGrid from global index array.
      DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! Create the mesh.
      call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=cvalue, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      EMesh = ESMF_MeshCreate(filename=trim(cvalue), &
                              fileformat=ESMF_FILEFORMAT_ESMFMESH, &
                              elementDistgrid=DistGrid, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call blom_logwrite(subname//': mesh file for blom domain is '// &
                         trim(cvalue))

      ! ------------------------------------------------------------------------
      ! Check for consistency of lat, lon and mask between mesh and model grid.
      ! ------------------------------------------------------------------------

      call ESMF_MeshGet(Emesh, spatialDim=spatialDim, &
                        numOwnedElements=numOwnedElements, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      allocate(ownedElemCoords(spatialDim*numOwnedElements), &
               lonMesh(numOwnedElements), latMesh(numOwnedElements), &
               maskMesh(numOwnedElements))

      call ESMF_MeshGet(Emesh, ownedElemCoords=ownedElemCoords, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      do n = 1, numOwnedElements
         lonMesh(n) = ownedElemCoords(2*n-1)
         latMesh(n) = ownedElemCoords(2*n)
      enddo

      elemMaskArray = ESMF_ArrayCreate(Distgrid, maskMesh, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_MeshGet(Emesh, elemMaskArray=elemMaskArray, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call blom_checkmesh(lonMesh, latMesh, maskMesh)

      ! ------------------------------------------------------------------------
      ! Determine flux area correction factors.
      ! ------------------------------------------------------------------------

      field = ESMF_FieldCreate(Emesh, ESMF_TYPEKIND_R8, &
                               meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_FieldRegridGetArea(field, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_FieldGet(field, farrayPtr=areaMesh, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call blom_setareacor(areaMesh, maskMesh)

      ! ------------------------------------------------------------------------
      ! Realize the actively coupled fields.
      ! ------------------------------------------------------------------------

      call fldlist_realize(state=importState, &
                           fldlist_num=fldsToOcn_num, fldlist=fldsToOcn, &
                           tag=subname//':BLOM_Import', mesh=EMesh, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call fldlist_realize(state=exportState, &
                           fldlist_num=fldsFrOcn_num, fldlist=fldsFrOcn, &
                           tag=subname//':BLOM_Export', mesh=EMesh, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! Set scalar data in export state.
      ! ------------------------------------------------------------------------

      call blom_getglobdim(nx_global, ny_global)

      call state_setscalar(real(nx_global, ESMF_KIND_R8), &
                           flds_scalar_index_nx, exportState, &
                           flds_scalar_name, flds_scalar_num, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call state_setscalar(real(ny_global, ESMF_KIND_R8), &
                           flds_scalar_index_ny, exportState, &
                           flds_scalar_name, flds_scalar_num, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine InitializeRealize

   subroutine DataInitialize(gcomp, rc)
   ! ---------------------------------------------------------------------------
   ! Called by NUOPC to do the initial data export from ocean to mediator.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(DataInitialize)'

      ! Local variables.
      type(ESMF_State) :: exportState
      type(ESMF_StateItem_flag) :: itemType

      if (dbug > 5) call ESMF_LogWrite(subname//': called', ESMF_LOGMSG_INFO)

      ! ------------------------------------------------------------------------
      ! Query the Component for its exportState.
      ! ------------------------------------------------------------------------

      call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! TODO
      ! ------------------------------------------------------------------------

      tlast_coupled = 0._ESMF_KIND_R8
      call blom_accflds
      call ocn_export(exportState, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! Check whether all Fields in the exportState are "Updated" TODO
      ! ------------------------------------------------------------------------

      if (NUOPC_IsUpdated(exportState)) then
         call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", &
                                     value="true", rc=rc)
         call ESMF_LogWrite("BLOM - Initialize-Data-Dependency SATISFIED!!!", &
                            ESMF_LOGMSG_INFO)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
      else
         call ESMF_LogWrite("BLOM - Initialize-Data-Dependency NOT SATISFIED!!!", &
                            ESMF_LOGMSG_INFO)
      endif

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine DataInitialize

   subroutine ModelAdvance(gcomp, rc)
   ! ---------------------------------------------------------------------------
   ! Called by NUOPC to advance the model a single timestep.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(ModelAdvance)'

      ! Local variables.
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      type(ESMF_Time) :: currTime
      type(ESMF_Alarm) :: restart_alarm
      integer :: shrlogunit, yr_sync, mon_sync, day_sync, tod_sync, ymd_sync, &
                 ymd, tod
      logical :: first_call = .true., restart_alarm_on
      character(len=cllen) :: msg

      if (dbug > 5) call ESMF_LogWrite(subname//': called', ESMF_LOGMSG_INFO)

      rc = ESMF_SUCCESS

      ! ------------------------------------------------------------------------
      ! Reset shr logging to components log file.
      ! ------------------------------------------------------------------------

      call shr_file_getLogUnit(shrlogunit)
      call shr_file_setLogUnit(lp)

      ! ------------------------------------------------------------------------
      ! Skip first coupling interval for an initial run.
      ! ------------------------------------------------------------------------

      if (first_call) then
         first_call = .false.
         if (runtyp_cesm == 'initial') then
            call blom_logwrite('Returning at first coupling interval')
            return
         endif
      endif

      ! ------------------------------------------------------------------------
      ! Query the Component for its clock, importState and exportState.
      ! ------------------------------------------------------------------------

      call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, &
                          exportState=exportState, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! Check that internal clock is in sync with master clock.
      ! ------------------------------------------------------------------------

      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_TimeGet(currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, &
                        s=tod_sync, rc=rc )
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

      call blom_time(ymd, tod)
      if (ymd /= ymd_sync .or. tod /= tod_sync) then
         write(msg,*) ' blom ymd=',ymd     ,'  blom tod= ',tod
         call blom_logwrite(msg)
         write(msg,*) ' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
         call blom_logwrite(msg)
         call ESMF_LogSetError(ESMF_FAILURE, &
            msg=subname//": Internal blom clock not in sync with Sync Clock", &
            line=__LINE__, file=u_FILE_u, rcToReturn=rc)
         return
      endif

      ! ------------------------------------------------------------------------
      ! Advance the model in time over a coupling interval.
      ! ------------------------------------------------------------------------

      blom_loop: do

         if (nint(tlast_coupled) == 0) then
            ! Obtain import state from driver
            call ocn_import(importState, rc)
            if (ChkErr(rc, __LINE__, u_FILE_u)) return
         endif
      
         ! Advance the model a time step.
         call blom_step

         ! Accumulate BLOM export fields.
         call blom_accflds

         if (nint(ocn_cpl_dt_cesm-tlast_coupled) == 0) then
            ! Return export state to driver and exit integration loop
            call ocn_export(exportState, rc)
            exit blom_loop
         endif

!        if (mnproc == 1) then
!           call shr_sys_flush(lp)
!        endif

      enddo blom_loop

      ! ------------------------------------------------------------------------
      ! If restart alarm is ringing - write restart file. TODO do we need to
      ! consider stop alarm?
      ! ------------------------------------------------------------------------

      call ESMF_ClockGetAlarm(clock, alarmname='restart_alarm', &
                              alarm=restart_alarm, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      restart_alarm_on = ESMF_AlarmIsRinging(restart_alarm, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      if (restart_alarm_on) then

         ! Turn off the alarm
         call ESMF_AlarmRingerOff(restart_alarm, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

         ! Write BLOM restart files.
         call restart_wt

      endif

      ! ------------------------------------------------------------------------
      ! Reset shr logging to original values.
      ! ------------------------------------------------------------------------

      call shr_file_setLogUnit(shrlogunit)

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine ModelAdvance

   subroutine ModelSetRunClock(gcomp, rc)
   ! ---------------------------------------------------------------------------
   ! Synchronize driver and model clock and set restart and stop alarms.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(ModelSetRunClock)'

      ! Local variables.
      type(ESMF_Clock) :: mclock, dclock
      type(ESMF_Time) :: mcurrtime, dcurrtime, mstoptime
      type(ESMF_TimeInterval) :: mtimestep, dtimestep
      type(ESMF_ALARM) :: restart_alarm, stop_alarm
      integer :: restart_n, restart_ymd, stop_n, stop_ymd, alarmcount
      character(len=256) :: cvalue, restart_option, stop_option
      character(len=128) :: name

      if (dbug > 5) call ESMF_LogWrite(subname//': called', ESMF_LOGMSG_INFO)

      ! Query the component for its clocks.

      call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! Force model clock currtime and timestep to match driver and set
      ! stoptime.
      ! ------------------------------------------------------------------------

      mstoptime = mcurrtime + dtimestep
      call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, &
                         stopTime=mstoptime, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! ------------------------------------------------------------------------
      ! Set restart and stop alarms.
      ! ------------------------------------------------------------------------

      call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, &
                                  alarmCount=alarmCount, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      if (alarmCount == 0) then

         call ESMF_GridCompGet(gcomp, name=name, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         call ESMF_LogWrite(subname//': setting alarms for '//trim(name), &
                            ESMF_LOGMSG_INFO)


         ! Restart alarm.

         call NUOPC_CompAttributeGet(gcomp, name="restart_option", &
                                     value=restart_option, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

         call NUOPC_CompAttributeGet(gcomp, name="restart_n", &
                                     value=cvalue, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         read(cvalue,*) restart_n

         call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", &
                                     value=cvalue, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         read(cvalue,*) restart_ymd

         call alarmInit(mclock, restart_alarm, restart_option, &
                        opt_n=restart_n, opt_ymd=restart_ymd, &
                        RefTime=mcurrTime, alarmname='restart_alarm', rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

         call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

         ! Stop alarm.

         call NUOPC_CompAttributeGet(gcomp, name="stop_option", &
                                     value=stop_option, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

         call NUOPC_CompAttributeGet(gcomp, name="stop_n", &
                                     value=cvalue, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         read(cvalue,*) stop_n

         call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", &
                                     value=cvalue, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return
         read(cvalue,*) stop_ymd

         call alarmInit(mclock, stop_alarm, stop_option, &
                        opt_n=stop_n, opt_ymd=stop_ymd, RefTime=mcurrTime, &
                        alarmname='stop_alarm', rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

         call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
         if (ChkErr(rc, __LINE__, u_FILE_u)) return

      endif

      ! ------------------------------------------------------------------------
      ! Advance model clock to trigger alarms then reset model clock back to
      ! currtime.
      ! ------------------------------------------------------------------------

      call ESMF_ClockAdvance(mclock, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, &
                         stopTime=mstoptime, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine ModelSetRunClock

   subroutine ModelFinalize(gcomp, rc)
   ! ---------------------------------------------------------------------------
   ! Called by NUOPC to finalize the model.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(ModelFinalize)'

      if (dbug > 5) call ESMF_LogWrite(subname//': called', ESMF_LOGMSG_INFO)

      rc = ESMF_SUCCESS

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine ModelFinalize

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine SetServices(gcomp, rc)
   ! ---------------------------------------------------------------------------
   ! NUOPC SetService method is the only public entry point. SetServices
   ! registers all of the user-provided subroutines in the module with the NUOPC
   ! layer.
   ! ---------------------------------------------------------------------------

      ! Input/output arguments.
      type(ESMF_GridComp)  :: gcomp ! ESMF_GridComp object.
      integer, intent(out) :: rc    ! Return code.

      ! Local parameters.
      character(len=*), parameter :: &
         subname = modname//':(SetServices)'

      if (dbug > 5) call ESMF_LogWrite(subname//': called', ESMF_LOGMSG_INFO)

      ! The NUOPC gcomp component will register the generic methods.
      call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! Switching to IPD versions.
      call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
                                      userRoutine=InitializeP0, phase=0, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! Set entry point for methods that require specific implementation.
      call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
                                   phaseLabelList=(/"IPDv01p1"/), &
                                   userRoutine=InitializeAdvertise, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
                                   phaseLabelList=(/"IPDv01p3"/), &
                                   userRoutine=InitializeRealize, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      ! Attach specializing method(s).
      call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
                                specRoutine=DataInitialize, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
                                specRoutine=ModelAdvance, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
                                specRoutine=ModelSetRunClock, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

! TODO Method used by POP, not by MOM6 and CICE.
!     call ESMF_MethodRemove(gcomp, label=model_label_CheckImport, rc=rc)
!     if (ChkErr(rc, __LINE__, u_FILE_u)) return
!     call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
!                               specRoutine=ModelCheckImport, rc=rc)
!     if (ChkErr(rc, __LINE__, u_FILE_u)) return

      call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
                                specRoutine=ModelFinalize, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return

      if (dbug > 5) call ESMF_LogWrite(subname//': done', ESMF_LOGMSG_INFO)

   end subroutine SetServices

end module ocn_comp_nuopc
