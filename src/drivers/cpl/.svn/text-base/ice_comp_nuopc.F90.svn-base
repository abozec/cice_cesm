#define FILENAME "ice_comp_nuopc.F90"

module ice_comp_nuopc_mod

#ifdef ESMF_INTERFACE


!-------------------------------------------------------------------------------
!
! Purpose: contains nuopc comp code for ice
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! share code & libs
   !----------------------------------------------------------------------------
   use ESMF
   use NUOPC
   use NUOPC_Model, only: &
      Model_routine_SetServices   => routine_SetServices, &
      Model_label_Advance         => label_Advance

   use ice_comp_esmf
   use seq_flds_mod
   use esmfshr_nuopc_mod
   use seq_infodata_mod
   use shr_kind_mod,      only : SHR_KIND_R8

   implicit none

   private

   public SetServices


!===============================================================================
contains
!===============================================================================


  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc


    ! local variables

    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, Model_routine_SetServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Run routine to execute run functionality
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      userRoutine=routine_Run2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    !-----------------------------------------------------------------

    rc = ESMF_SUCCESS

    !! Setup import-able fields
    call esmfshr_nuopc_advertise_fields( &
      ice_import_fields, importState, tag='CICE import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Setup export-able fields
    call esmfshr_nuopc_advertise_fields( &
      ice_export_fields, exportState, tag='CICE export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeP2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Field)        :: field
    type(ESMF_Mesh)         :: meshIn
    type(ESMF_Mesh)         :: meshOut
    type(ESMF_Array)        :: d2x
    type(ESMF_DistGrid)     :: distgrid
    real (SHR_KIND_R8)      :: nextsw_cday

    rc = ESMF_SUCCESS

    call seq_infodata_Getdata(infodata, nextsw_cday=nextsw_cday)
    call ESMF_AttributeSet(exportState, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ice_init_esmf(gcomp, importState, exportState, ccsm_EClock_i, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    ! create a Mesh object for Fields
    call ESMF_StateGet(exportState, itemName="d2x", array=d2x, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(d2x, distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    meshIn = ESMF_MeshCreate(distgrid, nodalDistGrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    meshOut = meshIn 

    !! Create and Realize Importable fields
    call esmfshr_nuopc_create_fields( &
      ice_import_fields, meshIn, importState, tag='CICE import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Create and Realize Exportable fields
    call esmfshr_nuopc_create_fields( &
      ice_export_fields, meshOut, exportState, tag='CICE export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    ! Copy output from CICE to export State
    call esmfshr_nuopc_copy(ice_export_fields, 'd2x', exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)   :: clock
    type(ESMF_State)   :: importState, exportState
    logical            :: ice_present
    logical            :: icerun_alarm

    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_ClockPrintCurrTime(clock, &
      "------>Advancing ICE from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_ClockPrintStopTime(clock, &
      "--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine


  !-----------------------------------------------------------------------------

  subroutine routine_Run2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)   :: compclock
    logical            :: ice_present
    logical            :: icerun_alarm

    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=compclock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
        name="ice_present", value=ice_present, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
        name="icerun_alarm", value=icerun_alarm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (ice_present  .and.  icerun_alarm) then
      ! Copy import Fields to import data structure before run
      call esmfshr_nuopc_copy(ice_import_fields, importState, 'x2d', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      call ice_run_esmf(gcomp, importState, exportState, ccsm_EClock_i, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      ! Copy export Fields to export data structure after run
      call esmfshr_nuopc_copy(ice_export_fields, 'd2x', exportState, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    end if

  end subroutine

#endif

end module ice_comp_nuopc_mod
