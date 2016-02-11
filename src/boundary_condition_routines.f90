!> \file 
!> \author Ting Yu
!> \brief This module set the boundary conditions for the given equation set
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module handles all boundary conditions routines.
MODULE BOUNDARY_CONDITIONS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE MPI
  USE NODE_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES
  USE LISTS
  USE LINKEDLIST_ROUTINES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_DOFTypes BOUNDARY_CONDITIONS_ROUTINES::DOFTypes
  !> \brief DOF type for boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FREE=0 !<The dof is free. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_MIXED=2 !<The dof is set as a mixed boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}
  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions BOUNDARY_CONDITIONS_ROUTINES::BoundaryConditions
  !> \brief Boundary conditions types. These may be specific to a particular equation type and the solver routines should not need to use these.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE=0 !<The dof is free. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED=1 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INLET=2 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_OUTLET=3 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_WALL=4 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL=5 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE_WALL=6 !<The dof is fixed as a boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT=8 !<The dof is set to a Neumann point boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED=9 !<The dof is set to a Neumann integrated boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DIRICHLET=10 !<The dof is set to a Dirichlet boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CAUCHY=11 !<The dof is set to a Cauchy boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_ROBIN=12 !<The dof is set to a Robin boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INCREMENTED=13 !<The dof is a fixed boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE=14 !<The dof is a surface pressure boundary condition. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE_INCREMENTED=15 !<The dof is a surface pressure boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED=17 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE=18 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_IMPERMEABLE_WALL=19 !<The dof is set such that (via penalty formulation): velocity * normal = 0. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY=20 !<A Neumann integrated boundary condition, and no point values will be integrated over a face or line that includes this dof. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_LINEAR_CONSTRAINT=21 !<The dof is constrained to be a linear combination of other DOFs. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED=22!<A Neumann point boundary condition that is incremented inside a load increment control loop. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_FITTED=23 !<The dof is fixed as a boundary condition to be updated from fitting data \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_NONREFLECTING=24 !<The dof is fixed and set to a non-reflecting type for 1D wave propagation problems. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_CELLML=25 !<The dof is fixed and set to values specified based on the coupled CellML solution at the dof. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_STREE=26 !<The dof is fixed and set to values specified based on the transmission line theory at the dof. \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
  !>@}

  INTEGER(INTG), PARAMETER :: MAX_BOUNDARY_CONDITION_NUMBER=26 !The maximum boundary condition type identifier, used for allocating an array with an entry for each type

  !> \addtogroup BOUNDARY_CONDITIONS_ROUTINES_SparsityTypes BOUNDARY_CONDITIONS_ROUTINES::BoundaryConditions
  !> \brief Storage type for matrices used by boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_SPARSE_MATRICES=1 !<The matrices are stored as sparse matrices.
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FULL_MATRICES=2 !<The matrices are stored as full matrices.
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Adds to the value of the specified local DOF and sets this as a boundary condition on the specified local DOF.
  INTERFACE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF
    MODULE PROCEDURE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1
    MODULE PROCEDURE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS
  END INTERFACE !BOUNDARY_CONDITIONS_ADD_LOCAL_DOF

  !>Sets a boundary condition on the specified local DOF. 
  INTERFACE BOUNDARY_CONDITIONS_SET_LOCAL_DOF
    MODULE PROCEDURE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1
    MODULE PROCEDURE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS
  END INTERFACE !BOUNDARY_CONDITIONS_SET_LOCAL_DOF

  PUBLIC BOUNDARY_CONDITION_DOF_FREE,BOUNDARY_CONDITION_DOF_FIXED,BOUNDARY_CONDITION_DOF_MIXED

  PUBLIC BOUNDARY_CONDITION_FREE,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_FIXED_INLET,&
    & BOUNDARY_CONDITION_FIXED_OUTLET,BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL,BOUNDARY_CONDITION_FREE_WALL,&
    & BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_DIRICHLET,BOUNDARY_CONDITION_NEUMANN_POINT, &
    & BOUNDARY_CONDITION_CAUCHY,BOUNDARY_CONDITION_ROBIN,BOUNDARY_CONDITION_FIXED_INCREMENTED,BOUNDARY_CONDITION_PRESSURE,&
    & BOUNDARY_CONDITION_PRESSURE_INCREMENTED,BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED, &
    & BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,BOUNDARY_CONDITION_IMPERMEABLE_WALL,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
    & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED,BOUNDARY_CONDITION_FIXED_STREE, &
    & BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML

  PUBLIC BOUNDARY_CONDITION_SPARSE_MATRICES,BOUNDARY_CONDITION_FULL_MATRICES

  PUBLIC BOUNDARY_CONDITIONS_CREATE_FINISH,BOUNDARY_CONDITIONS_CREATE_START,BOUNDARY_CONDITIONS_DESTROY
  
  PUBLIC BOUNDARY_CONDITIONS_ADD_CONSTANT,BOUNDARY_CONDITIONS_ADD_LOCAL_DOF,BOUNDARY_CONDITIONS_ADD_ELEMENT, &
    & BOUNDARY_CONDITIONS_ADD_NODE,BOUNDARY_CONDITIONS_VARIABLE_GET

  PUBLIC BOUNDARY_CONDITIONS_SET_CONSTANT,BOUNDARY_CONDITIONS_SET_LOCAL_DOF,BOUNDARY_CONDITIONS_SET_ELEMENT, &
    & BOUNDARY_CONDITIONS_SET_NODE,BoundaryConditions_NeumannIntegrate,BoundaryConditions_NeumannSparsityTypeSet

  PUBLIC BoundaryConditions_ConstrainNodeDofsEqual

  PUBLIC BoundaryConditions_ConstraintsApply

CONTAINS  

  !
  !================================================================================================================================
  !

  !>Finish the creation of boundary conditions.
  SUBROUTINE BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    LOGICAL :: parameterSetRequired(FIELD_NUMBER_OF_VARIABLE_TYPES)
    INTEGER(INTG) :: MPI_IERROR
    INTEGER(INTG) :: variableIdx,parameterSetIdx,conditionTypeIdx
    INTEGER(INTG), POINTER :: conditionTypes(:),dofTypes(:)
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITION_VARIABLE
    TYPE(BoundaryConditionsCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
          IF(COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES>0) THEN
            !Transfer all the boundary conditions to all the computational nodes.
            DO variableIdx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
              BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
                CALL DistributedVector_UpdateStart(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES,err,error,*999)
                CALL DistributedVector_UpdateStart(BOUNDARY_CONDITION_VARIABLE%DOF_TYPES,err,error,*999)
              ELSE
                CALL FlagError("Boundary conditions variable is not associated for variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(variableIdx,"*",ERR,ERROR))//".",ERR,ERROR,*999)
              ENDIF
            ENDDO ! variableIdx
            DO variableIdx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
              BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR
              IF(ASSOCIATED(BOUNDARY_CONDITION_VARIABLE)) THEN
                createValuesCache=>BOUNDARY_CONDITION_VARIABLE%createValuesCache
                IF(ASSOCIATED(createValuesCache)) THEN
                  FIELD_VARIABLE=>BOUNDARY_CONDITION_VARIABLE%VARIABLE
                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                    CALL DistributedVector_UpdateFinish(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES,err,error,*999)
                    CALL DistributedVector_UpdateFinish(BOUNDARY_CONDITION_VARIABLE%DOF_TYPES,err,error,*999)
                    !Count the total number of local dofs
                    NULLIFY(conditionTypes)
                    CALL DistributedVector_DataGet(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES,conditionTypes,err,error,*999)
                    DO conditionTypeIdx=1,MAX_BOUNDARY_CONDITION_NUMBER
                      createValuesCache%dofCounts(conditionTypeIdx)=COUNT(conditionTypes==conditionTypeIdx)
                    END DO !conditionTypeIdx
                    NULLIFY(dofTypes)
                    CALL DistributedVector_DataGet(BOUNDARY_CONDITION_VARIABLE%DOF_TYPES,dofTypes,err,error,*999)

                    ! Update the total number of boundary condition types by summing across all nodes
                    CALL MPI_ALLREDUCE(createValuesCache%dofCounts>0,BOUNDARY_CONDITION_VARIABLE%anyGlobalDofs, &
                      & MAX_BOUNDARY_CONDITION_NUMBER,MPI_LOGICAL,MPI_LOR,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)

                    ! Check that the boundary conditions set are appropriate for equations sets
                    CALL BoundaryConditions_CheckEquations(BOUNDARY_CONDITION_VARIABLE,parameterSetRequired,ERR,ERROR,*999)

                    !Make sure the required parameter sets are created on all computational nodes and begin updating them
                    DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                      IF(parameterSetRequired(parameterSetIdx)) THEN
                        CALL Field_ParameterSetEnsureCreated(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                          & parameterSetIdx,ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_UPDATE_START(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                          & parameterSetIdx,ERR,ERROR,*999)
                      END IF
                    END DO

                    ! Set up pressure incremented condition, if it exists
                    IF(BOUNDARY_CONDITION_VARIABLE%anyGlobalDofs(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)) THEN
                      CALL BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                    END IF

                    ! Set up Neumann condition information if there are any Neumann conditions
                    IF(BOUNDARY_CONDITION_VARIABLE%anyGlobalDofs(BOUNDARY_CONDITION_NEUMANN_POINT).OR. &
                        & BOUNDARY_CONDITION_VARIABLE%anyGlobalDofs(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)) THEN
                      CALL BoundaryConditions_NeumannInitialise(BOUNDARY_CONDITION_VARIABLE,ERR,ERROR,*999)
                    END IF
                    
                    !Finish creating the boundary conditions DOF constraints
                    CALL BoundaryConditions_ConstraintsCreateFinish(boundary_condition_variable,err,error,*999)

                    CALL DistributedVector_DataRestore(BOUNDARY_CONDITION_VARIABLE%CONDITION_TYPES,conditionTypes,err,error,*999)
                    CALL DistributedVector_DataRestore(BOUNDARY_CONDITION_VARIABLE%DOF_TYPES,dofTypes,err,error,*999)

                    ! Finish field update
                    DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                      IF(parameterSetRequired(parameterSetIdx)) THEN
                        CALL FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD_VARIABLE%FIELD,FIELD_VARIABLE%VARIABLE_TYPE, &
                          & parameterSetIdx,ERR,ERROR,*999)
                      END IF
                    END DO
                    !Finalise the create values cache.
                    CALL BoundaryConditions_CreateValuesCacheFinalise(createValuesCache,err,error,*999)
                  ELSE
                    LOCAL_ERROR="Field variable is not associated for variable index "// &
                      & TRIM(NUMBER_TO_VSTRING(variableIdx,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Boundary conditions variable create values cache is not associated for variable index "// &
                      & TRIM(NUMBER_TO_VSTRING(variableIdx,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Boundary conditions variable is not associated for variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(variableIdx,"*",ERR,ERROR))//".",ERR,ERROR,*999)
              ENDIF
            ENDDO ! variableIdx
          ENDIF
          !Set the finished flag
          BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.TRUE.
        ELSE
          CALL FlagError("Boundary conditions variables array is not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Boundary conditions:",ERR,ERROR,*999)
      DO variableIdx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
        BOUNDARY_CONDITION_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",BOUNDARY_CONDITION_VARIABLE%VARIABLE_TYPE, &
            & ERR,ERROR,*999)
      ENDDO !variableIdx
    ENDIF
    
    EXITS("BOUNDARY_CONDITIONS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of boundary conditions for the equation set.
  SUBROUTINE BOUNDARY_CONDITIONS_CREATE_START(SOLVER_EQUATIONS,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to create boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<On exit, a pointer to the created boundary conditions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Boundary conditions are already associated for the solver equations.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
          CALL FlagError("Boundary conditions is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
            !Initialise the boundary conditions
            CALL BOUNDARY_CONDITIONS_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Solver equations solver mapping is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
          !Return the pointer
          BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_CREATE_START")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR)
    RETURN 1

  END SUBROUTINE BOUNDARY_CONDITIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys boundary conditions
  SUBROUTINE BOUNDARY_CONDITIONS_DESTROY(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      CALL BOUNDARY_CONDITIONS_FINALISE(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_DESTROY")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions and deallocate all memory.
  SUBROUTINE BOUNDARY_CONDITIONS_FINALISE(BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx

    ENTERS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
        DO variableIdx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
          IF(ASSOCIATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR, &
                & ERR,ERROR,*999)
          ELSE
            CALL FlagError("Boundary conditions variable number "//TRIM(NUMBER_TO_VSTRING(variableIdx,"*",ERR,ERROR))// &
                  & " is not associated",ERR,ERROR,*999)
          ENDIF
        ENDDO !variableIdx
        NULLIFY(BOUNDARY_CONDITIONS%SOLVER_EQUATIONS%SOLVER%SOLVER_EQUATIONS)
        !BOUNDARY_CONDITIONS%SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED = .FALSE.
        !BOUNDARY_CONDITIONS%SOLVER_EQUATIONS%SOLVER_MAPPING%SOLVER_MAPPING_FINISHED = .FALSE.
        DEALLOCATE(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)
      ENDIF
      DEALLOCATE(BOUNDARY_CONDITIONS)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_FINALISE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the boundary conditions for an equations set.
  SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to initialise the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,equations_set_idx,interface_condition_idx
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: INTERFACE_RHS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Boundary conditions is already associated for these solver equations.",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
          ALLOCATE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate boundary conditions.",ERR,ERROR,*999)
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED=.FALSE.
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES=0
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
          SOLVER_EQUATIONS%BOUNDARY_CONDITIONS%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
          DO equations_set_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_SET
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                  IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
                      EQUATIONS_SET%BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                        & EQUATIONS_MAPPING%DEPENDENT_VARIABLE,ERR,ERROR,*999)
                      RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                      IF(ASSOCIATED(RHS_MAPPING)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                            & RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping has not been finished.",ERR,ERROR,*998)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations equations mapping is not associated.",ERR,ERROR,*998)
                  ENDIF
                ELSE
                  CALL FlagError("Equations has not been finished.",ERR,ERROR,*998)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",ERR,ERROR,*998)
            ENDIF
          ENDDO !equations_set_idx
          DO interface_condition_idx=1,SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
            INTERFACE_CONDITION=>SOLVER_EQUATIONS%SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
              & INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
                  INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                  IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                    IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
                      INTERFACE_CONDITION%BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                      IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                          & INTERFACE_MAPPING%LAGRANGE_VARIABLE,ERR,ERROR,*999)
                        INTERFACE_RHS_MAPPING=>INTERFACE_MAPPING%RHS_MAPPING
                        IF(ASSOCIATED(INTERFACE_RHS_MAPPING)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                              & INTERFACE_RHS_MAPPING%RHS_VARIABLE,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Interface mapping mapping is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Interface mapping has not been finished.",ERR,ERROR,*998)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*998)
                  ENDIF
                ELSE
                  CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*998)
                ENDIF
              ELSE
                CALL FlagError("Interface equations is not associated.",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FlagError("Interface condition not associated.",ERR,ERROR,*998)
            ENDIF
          ENDDO !interface_condition_idx
        ELSE
          CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_FINALISE(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("BOUNDARY_CONDITIONS_INITIALISE",ERR,ERROR)
    RETURN 1

  END SUBROUTINE BOUNDARY_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionAddConstant
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_CONSTANT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("BOUNDARY_CONDITIONS_ADD_CONSTANT",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DEPENDENT_FIELD_VARIABLE)

    !Note: This routine is for constant interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,local_ny,global_ny, &
            & ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,DEPENDENT_FIELD_VARIABLE,ERR,ERROR,*999)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
            & ERR,ERROR,*999)
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
              & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_ADD_CONSTANT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_ADD_CONSTANT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_CONSTANT
  
 !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionsSetConstant
  SUBROUTINE BOUNDARY_CONDITIONS_SET_CONSTANT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("BOUNDARY_CONDITIONS_SET_CONSTANT",ERR,ERROR,*999)

    !Note: This routine is for constant interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,local_ny,global_ny, &
            & ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
              & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_SET_CONSTANT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_SET_CONSTANT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_CONSTANT
  
  !
  !================================================================================================================================
  !
 
  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOF.
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,[DOF_INDEX],[CONDITION],[VALUE], &
        & ERR,ERROR,*999)

    EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOF1
  
  !
  !================================================================================================================================
  !
 
  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOFs.
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:) !<DOF_INDICES(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:) !<CONDITIONS(:). The boundary condition type to set for the i'th dof \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(:). The value of the boundary condition for the i'th dof to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,local_ny
    REAL(DP) :: INITIAL_VALUE
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR,*999)
    NULLIFY(dependent_variable)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          NULLIFY(DEPENDENT_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,DEPENDENT_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
            DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
            IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                & ERR,ERROR,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
                  IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                    DO i=1,SIZE(DOF_INDICES,1)
                      local_ny=DOF_INDICES(i)
                      IF(1<=local_ny.AND.local_ny<=DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                        ! Set boundary condition and dof type, and make sure parameter sets are created
                        CALL BoundaryConditions_SetConditionType(BOUNDARY_CONDITIONS_VARIABLE,local_ny,CONDITIONS(i), &
                          & ERR,ERROR,*999)
                        ! Update field sets by adding boundary condition values
                        SELECT CASE(CONDITIONS(i))
                        CASE(BOUNDARY_CONDITION_FREE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_FIXED)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INLET)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FREE_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUES(i), &
                            & ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
                          ! For increment loops, we need to set the full BC parameter set value by
                          ! getting the current value from the values parameter set
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,INITIAL_VALUE,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,INITIAL_VALUE+VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                          ! For pressure incremented, adding to the values_set parameter value doesn't make sense,
                          ! so just increment the value in the pressure values parameter set
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
                          ! Point value is stored in boundary conditions field set, and is then integrated to
                          ! get the RHS variable value
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
                          ! For integrated Neumann condition, integration is already done, so set the RHS
                          ! dof value directly
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
                          &  BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
                          CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The specified boundary condition type for dof index "// &
                            & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        LOCAL_ERROR="The local dof of  "//&
                          & TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))//" at dof index "// &
                          & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                          & " is invalid. The dof should be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !i
                  ELSE
                    LOCAL_ERROR="The size of the dof indices array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                      & ") does not match the size of the values array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The size of the dof indices array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                    & ") does not match the size of the fixed conditions array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The dependent field variable domain mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated..",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_LOCAL_DOFS
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified DOF.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR,*999)

    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,[DOF_INDEX],[CONDITION],[VALUE], &
      & ERR,ERROR,*999)

    EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOF1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOF1
  
  !
  !================================================================================================================================
  !
 
  !>Sets a boundary condition on the specified DOFs.
  SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:) !<DOF_INDICES(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:) !<CONDITIONS(:). The boundary condition type to set for the i'th dof \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(:). The value of the boundary condition for the i'th dof to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,local_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          NULLIFY(DEPENDENT_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,DEPENDENT_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
            DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
            IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                  & ERR,ERROR,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
                  IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                    DO i=1,SIZE(DOF_INDICES,1)
                      local_ny=DOF_INDICES(i)
                      IF(1<=local_ny.AND.local_ny<=DOMAIN_MAPPING%NUMBER_OF_LOCAL) THEN
                        ! Set boundary condition and dof type
                        CALL BoundaryConditions_SetConditionType(BOUNDARY_CONDITIONS_VARIABLE,local_ny,CONDITIONS(i), &
                          & ERR,ERROR,*999)
                        ! Update field sets with boundary condition value
                        SELECT CASE(CONDITIONS(i))
                        CASE(BOUNDARY_CONDITION_FREE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_FIXED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INLET)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FREE_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny,VALUES(i), &
                            & ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
                          ! No field update
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                        CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
                            & BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & local_ny,VALUES(i),ERR,ERROR,*999)
                          CALL BoundaryConditions_DofConstraintSet(BOUNDARY_CONDITIONS_VARIABLE, &
                            & local_ny,[INTEGER(INTG)::],[REAL(DP)::],err,error,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The specified boundary condition type for dof index "// &
                            & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" is invalid."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        LOCAL_ERROR="The local dof of  "//&
                          & TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))//" at dof index "// &
                          & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                          & " is invalid. The dof should be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !i
                  ELSE
                    LOCAL_ERROR="The size of the dof indices array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                      & ") does not match the size of the values array ("// &
                      & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The size of the dof indices array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                    & ") does not match the size of the fixed conditions array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The dependent field variable domain mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_SET_LOCAL_DOFS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_LOCAL_DOFS

  !
  !================================================================================================================================
  !

  !> Checks the boundary condition type and sets the boundary condition type and dof type for the boundary conditions.
  SUBROUTINE BoundaryConditions_SetConditionType(boundaryConditionsVariable,localDof,condition,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: localDof !<The localDof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dofType

    ENTERS("BoundaryConditions_SetConditionType",err,error,*999)

    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_FIXED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_INLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_MOVED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FREE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_PRESSURE)
      ! Pressure boundary conditions leave the RHS dof as free, as the Neumann terms
      ! are calculated in finite elasticity routines when calculating the element residual
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type for dof number "// &
        & TRIM(NUMBER_TO_VSTRING(localDof,"*",err,error))//" of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT

    !Set the boundary condition type and DOF type
    CALL DistributedVector_ValuesSet(boundaryConditionsVariable%CONDITION_TYPES,localDof,condition,err,error,*999)
    CALL DistributedVector_ValuesSet(boundaryConditionsVariable%DOF_TYPES,localDof,dofType,err,error,*999)

    EXITS("BoundaryConditions_SetConditionType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetConditionType",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_SetConditionType

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsAddElement
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_ELEMENT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_ADD_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
            & local_ny,global_ny,ERR,ERROR,*999)
          NULLIFY(FIELD_VARIABLE)
          NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_ADD_ELEMENT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_ADD_ELEMENT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_ELEMENT
  
  !
  !================================================================================================================================
  !
 
  !> Checks that the specified boundary condition is appropriate for the field variable interpolation type
  SUBROUTINE BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*)

    ! Argument variables
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type being set
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number the boundary condition is set on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: interpolationType
    LOGICAL :: validCondition

    ENTERS("BoundaryConditions_CheckInterpolationType",err,error,*999)

    CALL FIELD_COMPONENT_INTERPOLATION_GET(field,variableType,componentNumber,interpolationType,err,error,*999)

    validCondition=.TRUE.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE, &
        & BOUNDARY_CONDITION_FIXED, &
        & BOUNDARY_CONDITION_FIXED_INCREMENTED)
      ! Valid for all interpolation types
    CASE(BOUNDARY_CONDITION_FIXED_INLET, &
        & BOUNDARY_CONDITION_FIXED_OUTLET)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_FIXED_WALL, &
        & BOUNDARY_CONDITION_MOVED_WALL, &
        & BOUNDARY_CONDITION_FREE_WALL, &
        & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_PRESSURE, &
        & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT, &
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED, &
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
        & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT
    IF(.NOT.validCondition) THEN
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NUMBER_TO_VSTRING(condition,"*",err,error))//" is not valid for the field component "// &
        & "interpolation type of "//TRIM(NUMBER_TO_VSTRING(interpolationType,"*",err,error))//".", &
        & err,error,*999)
    END IF

    EXITS("BoundaryConditions_CheckInterpolationType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckInterpolationType",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_CheckInterpolationType

  !
  !================================================================================================================================
  !

  !> Checks that the applied boundary conditions are supported by the equations sets in the solver equations
  SUBROUTINE BoundaryConditions_CheckEquations(boundaryConditionsVariable,parameterSetRequired,err,error,*)

    ! Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check
    LOGICAL, INTENT(OUT) :: parameterSetRequired(:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    type(varying_string), intent(out) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: boundaryConditionType,equationsSetIdx,specificationSize
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    LOGICAL :: validEquationsSetFound

    ENTERS("BoundaryConditions_CheckEquations",err,error,*999)

    !Get and check pointers we need
    solverEquations=>boundaryConditionsVariable%BOUNDARY_CONDITIONS%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(solverEquations)) THEN
      CALL FlagError("Boundary conditions solver equations are not associated.",err,error,*999)
    END IF
    solverMapping=>solverEquations%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(solverMapping)) THEN
      CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
    END IF

    parameterSetRequired=.FALSE.
    parameterSetRequired(FIELD_VALUES_SET_TYPE)=.TRUE.

    DO boundaryConditionType=1,MAX_BOUNDARY_CONDITION_NUMBER
      !Check if any DOFs have been set to this BC type
      IF(boundaryConditionsVariable%anyGlobalDofs(boundaryConditionType)>0) THEN
        validEquationsSetFound=.FALSE.
        DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
          equationsSet=>solverMapping%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS_SET
          IF(.NOT.ASSOCIATED(equationsSet)) THEN
            CALL FlagError("Solver equations equations set is not associated.",err,error,*999)
          END IF
          IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
            CALL FlagError("Equations set specification is not allocated.",err,error,*999)
          END IF
          specificationSize=SIZE(equationsSet%specification,1)

          SELECT CASE(boundaryConditionType)
          CASE(BOUNDARY_CONDITION_FREE)
            ! Valid for any equations set
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INLET, &
              & BOUNDARY_CONDITION_FIXED_OUTLET)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL, &
              & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED,BOUNDARY_CONDITION_FREE_WALL)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_DARCY_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              ELSE IF(specificationSize==3) THEN
                IF(equationsSet%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSet%specification(2)==EQUATIONS_SET_LAPLACE_EQUATION_TYPE.AND. &
                    & equationsSet%specification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
                  validEquationsSetFound=.TRUE.
                END IF
              END IF
            END IF
            IF(boundaryConditionType==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) &
              & parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
            validEquationsSetFound=.TRUE.
            parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
          CASE(BOUNDARY_CONDITION_PRESSURE,BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
                validEquationsSetFound=.TRUE.
              ELSE IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS .AND. &
                & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE) THEN
                validEquationsSetFound=.TRUE.
              END IF
            ENDIF
            parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
            IF(boundaryConditionType==BOUNDARY_CONDITION_PRESSURE_INCREMENTED) &
              & parameterSetRequired(FIELD_PREVIOUS_PRESSURE_SET_TYPE)=.TRUE.
          CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
            !Not actually used anywhere? So keep it as invalid, although maybe it should be removed?
            validEquationsSetFound=.FALSE.
          CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
            IF(specificationSize>=3) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                  & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
                  & (equationsSet%specification(3)==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE.OR. &
                  & equationsSet%specification(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
                  & equationsSet%specification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            END IF
            parameterSetRequired(FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE)=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
            validEquationsSetFound=.TRUE.
            parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
            parameterSetRequired(FIELD_INTEGRATED_NEUMANN_SET_TYPE)=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_FITTED)
            IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
              & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
            IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                & (equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE DEFAULT
            CALL FlagError("The specified boundary condition type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is invalid.",err,error,*999)
          END SELECT
        END DO
        IF(.NOT.validEquationsSetFound) THEN
            CALL FlagError("The specified boundary condition type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is invalid for the equations sets in the solver equations.",err,error,*999)
        END IF
      END IF
    END DO

    EXITS("BoundaryConditions_CheckEquations")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckEquations",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_CheckEquations

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsSetElement
  SUBROUTINE BOUNDARY_CONDITIONS_SET_ELEMENT(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & CONDITION,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("BOUNDARY_CONDITIONS_SET_ELEMENT",ERR,ERROR,*999)

    !Note: this routine is for element based interpolation
    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
            & local_ny,global_ny,ERR,ERROR,*999)
          NULLIFY(FIELD_VARIABLE)
          NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITION_SET_ELEMENT")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITION_SET_ELEMENT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_ELEMENT
  
  !
  !================================================================================================================================
  !
 
  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user node. \see OPENCMISS_CMISSBoundaryConditionsAddNode
  SUBROUTINE BOUNDARY_CONDITIONS_ADD_NODE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
    & USER_NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR,*999)

    NULLIFY(FIELD_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_NODE(FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
            & USER_NODE_NUMBER,COMPONENT_NUMBER,local_ny,global_ny,ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_ADD_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_ADD_NODE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_ADD_NODE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_ADD_NODE

  !
  !================================================================================================================================
  !

  !>Initialise the Neumann boundary conditions.
  SUBROUTINE BoundaryConditions_NeumannInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,numberOfNonZeros,neumannDofIdx,localDof,adjacentDomainIdx,componentIdx,nodeIdx,localNodeIdx1,localNodeIdx2, &
      & derivativeIdx1,derivativeIdx2,nodeNumber1,nodeNumber2,versionNumber1,versionNumber2,derivativeNumber1,derivativeNumber2, &
      & lineIdx,faceIdx,variableDofIdx,variableDofIdx1,variableDofIdx2,ghostSendIdx,ghostReceiveIdx,numberOfColumns,columnIdx,dummyErr
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:),columnIndices(:),columns(:)
    INTEGER(INTG), POINTER :: conditionTypes(:)
    REAL(DP), POINTER :: pointValuesData(:)
    TYPE(BoundaryConditionsCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannConditions 
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rhsVariableMapping,pointDofMapping
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(LIST_TYPE), POINTER :: ghostSendList,ghostReceiveList 
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_NeumannInitialise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      IF(ASSOCIATED(boundaryConditionsVariable%neumannBoundaryConditions)) THEN
        CALL FlagError("Neumann boundary conditions is already associated for this boundary conditions variable.",err,error,*999)
      END IF
      createValuesCache=>boundaryConditionsVariable%createValuesCache
      IF(.NOT.ASSOCIATED(createValuesCache)) &
        & CALL FlagError("Boundary condition variable create values cache is not associated.",err,error,*999)
      rhsVariable=>boundaryConditionsVariable%variable
      IF(.NOT.ASSOCIATED(rhsVariable)) &
        & CALL FlagError("RHS boundary conditions variable field variable is not associated.",err,error,*999)
      ALLOCATE(boundaryConditionsVariable%neumannBoundaryConditions,stat=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann Boundary Conditions",err,error,*998)
      neumannConditions=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(neumannConditions)) THEN
        NULLIFY(conditionTypes)
        CALL DistributedVector_DataGet(boundaryConditionsVariable%CONDITION_TYPES,conditionTypes,err,error,*999)

        neumannConditions%numberOfNeumann=createValuesCache%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)+ &
          & createValuesCache%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
        NULLIFY(neumannConditions%integrationMatrix)
        NULLIFY(neumannConditions%pointValues)
        NULLIFY(neumannConditions%pointDofMapping)
        ALLOCATE(neumannConditions%setDofs(neumannConditions%numberOfNeumann),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann set DOFs.",err,error,*999)
        neumannDofIdx=0
        DO variableDofIdx=1,rhsVariable%TOTAL_NUMBER_OF_DOFS
          IF(conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
              & conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
            neumannDofIdx=neumannDofIdx+1
            neumannConditions%setDofs(neumannDofIdx)=variableDofIdx
          END IF
        END DO
        !For rows we can re-use the RHS variable row mapping
        rhsVariableMapping=>rhsVariable%DOMAIN_MAPPING
        IF(.NOT.ASSOCIATED(rhsVariableMapping)) &
          & CALL FlagError("RHS field variable mapping is not associated.",err,error,*998)
        !Create a domain mapping for the Neumann point DOFs, required for the distributed matrix columns
        ALLOCATE(neumannConditions%pointDofMapping,stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann DOF domain mapping.",err,error,*999)
        pointDofMapping=>neumannConditions%pointDofMapping
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(pointDofMapping,rhsVariableMapping%NUMBER_OF_DOMAINS,err,error,*999)

        pointDofMapping%NUMBER_OF_LOCAL=COUNT(neumannConditions%setDofs<=rhsVariableMapping%NUMBER_OF_LOCAL)
        pointDofMapping%TOTAL_NUMBER_OF_LOCAL=neumannConditions%numberOfNeumann

        pointDofMapping%NUMBER_OF_ADJACENT_DOMAINS=rhsVariableMapping%NUMBER_OF_ADJACENT_DOMAINS
        ALLOCATE(pointDofMapping%ADJACENT_DOMAINS(pointDofMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
        IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)

        DO adjacentDomainIdx=1,pointDofMapping%NUMBER_OF_ADJACENT_DOMAINS
          CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(pointDofMapping%ADJACENT_DOMAINS(adjacentDomainIdx), &
            & err,error,*999)
          pointDofMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER= &
            & rhsVariableMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
          NULLIFY(ghostSendList)
          CALL LIST_CREATE_START(ghostSendList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(ghostSendList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_INITIAL_SIZE_SET(ghostSendList, &
            & MAX(pointDofMapping%TOTAL_NUMBER_OF_LOCAL-pointDofMapping%NUMBER_OF_LOCAL,1),err,error,*999)
          CALL LIST_CREATE_FINISH(ghostSendList,err,error,*999)
          NULLIFY(ghostReceiveList)
          CALL LIST_CREATE_START(ghostReceiveList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(ghostReceiveList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_INITIAL_SIZE_SET(ghostReceiveList, &
            & MAX(pointDofMapping%TOTAL_NUMBER_OF_LOCAL-pointDofMapping%NUMBER_OF_LOCAL,1),err,error,*999)
          CALL LIST_CREATE_FINISH(ghostReceiveList,err,error,*999)
          DO ghostSendIdx=1,rhsVariableMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS
            localDof=rhsVariableMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES(ghostSendIdx)
            DO neumannDofIdx=1,neumannConditions%numberOfNeumann
              IF(neumannConditions%setDofs(neumannDofIdx)==localDof) THEN
                CALL LIST_ITEM_ADD(ghostSendList,neumannDofIdx,err,error,*999)
                EXIT
              END IF
            END DO !neumannDofIdx
          END DO !ghostSendIdx
          DO ghostReceiveIdx=1,rhsVariableMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
            localDof=rhsVariableMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(ghostReceiveIdx)
            DO neumannDofIdx=1,neumannConditions%numberOfNeumann
              IF(neumannConditions%setDofs(neumannDofIdx)==localDof) THEN
                CALL LIST_ITEM_ADD(ghostReceiveList,neumannDofIdx,err,error,*999)
                EXIT
              END IF
            END DO !neumannDofIdx
          END DO !ghostReceiveIdx
          CALL LIST_REMOVE_DUPLICATES(ghostSendList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(ghostSendList, &
            & pointDofMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS, &
            & pointDofMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES,err,error,*999)
          CALL LIST_REMOVE_DUPLICATES(ghostReceiveList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(ghostReceiveList, &
            & pointDofMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS, &
            & pointDofMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES,err,error,*999)
        ENDDO !adjacentDomainIdx

        !Calculate the local to global map.
        CALL DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE(pointDofMapping,err,error,*999)

        !Set up vector of Neumann point values
        CALL DISTRIBUTED_VECTOR_CREATE_START(pointDofMapping,neumannConditions%pointValues,err,error,*999)
        CALL DISTRIBUTED_VECTOR_CREATE_FINISH(neumannConditions%pointValues,err,error,*999)

        NULLIFY(pointValuesData)
        CALL Field_ParameterSetDataGet(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
          & pointValuesData,err,error,*999)

        !Set point values vector from boundary conditions field parameter set
        CALL DistributedVector_ValuesSet(neumannConditions%pointValues,[(i,i=1,neumannConditions%numberOfNeumann)], &
          & pointValuesData(neumannConditions%setDofs),err,error,*999)
        CALL DISTRIBUTED_VECTOR_UPDATE_START(neumannConditions%pointValues,err,error,*999)

        CALL DISTRIBUTED_MATRIX_CREATE_START(rhsVariableMapping,pointDofMapping,neumannConditions%integrationMatrix,err,error,*999)
        SELECT CASE(boundaryConditionsVariable%BOUNDARY_CONDITIONS%neumannMatrixSparsity)
        CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
          !Work out integration matrix sparsity structure
          !For a single process, compressed column would be more memory efficient, but with
          !multiple processes the number of Neumann point DOFs could be more than the number
          !of local row DOFs, and multiplying a compressed row matrix by a vector is faster,
          !so we will use compressed row storage
          !Allocate lists
          ALLOCATE(columnIndicesLists(rhsVariableMapping%TOTAL_NUMBER_OF_LOCAL),stat=err)
          IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
          !Allocate row indices
          ALLOCATE(rowIndices(rhsVariableMapping%TOTAL_NUMBER_OF_LOCAL),stat=err)
          IF(ERR/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
          rowIndices(1)=1
          !Set up the column indicies lists
          DO variableDofIdx=1,rhsVariableMapping%TOTAL_NUMBER_OF_LOCAL
            NULLIFY(columnIndicesLists(variableDofIdx)%PTR)
            CALL LIST_CREATE_START(columnIndicesLists(variableDofIdx)%PTR,err,error,*999)
            CALL LIST_DATA_TYPE_SET(columnIndicesLists(variableDofIdx)%PTR,LIST_INTG_TYPE,err,error,*999)
            CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(variableDofIdx)%PTR,16,err,error,*999)
            CALL LIST_CREATE_FINISH(columnIndicesLists(variableDofIdx)%PTR,err,error,*999)
          ENDDO !variableDofIdx

          DO componentIdx=1,rhsVariable%NUMBER_OF_COMPONENTS
            SELECT CASE(rhsVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              !Get topology for finding faces/lines
              topology=>rhsVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY
              IF(.NOT.ASSOCIATED(topology)) THEN
                CALL FlagError("Field component topology is not associated.",err,error,*999)
              END IF
              IF(.NOT.ASSOCIATED(topology%NODES)) THEN
                CALL FlagError("Topology nodes are not associated.",err,error,*999)
              END IF
              SELECT CASE(rhsVariable%COMPONENTS(componentIdx)%DOMAIN%NUMBER_OF_DIMENSIONS)
              CASE(1)
                DO nodeIdx=1,topology%nodes%NUMBER_OF_NODES
                  IF(topology%nodes%nodes(nodeIdx)%BOUNDARY_NODE) THEN
                    variableDofIdx=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                    IF(conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      !Find the Neumann condition number
                      DO neumannDofIdx=1,neumannConditions%numberOfNeumann
                        IF(neumannConditions%setDofs(neumannDofIdx)==variableDofIdx) THEN
                          EXIT
                        END IF
                      END DO
                      CALL LIST_ITEM_ADD(columnIndicesLists(variableDofIdx)%PTR,neumannDofIdx,err,error,*999)
                    END IF
                  END IF
                END DO !nodeIdx
              CASE(2)
                IF(.NOT.ASSOCIATED(topology%lines)) THEN
                  CALL FlagError("Topology lines have not been calculated.",err,error,*999)
                END IF
                linesLoop: DO lineIdx=1,topology%lines%NUMBER_OF_LINES
                  IF(topology%lines%lines(lineIdx)%BOUNDARY_LINE) THEN
                    DO localNodeIdx2=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_NODES
                      nodeNumber2=topology%lines%lines(lineIdx)%NODES_IN_LINE(localNodeIdx2)
                      DO derivativeIdx2=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx2)
                        derivativeNumber2=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(1,derivativeIdx2,localNodeIdx2)
                        versionNumber2=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(2,derivativeIdx2,localNodeIdx2)
                        variableDofIdx2=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                          & NODES(nodeNumber2)%DERIVATIVES(derivativeNumber2)%VERSIONS(versionNumber2)
                        IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                          & conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                          !Find the Neumann condition number
                          DO neumannDofIdx=1,neumannConditions%numberOfNeumann
                            IF(neumannConditions%setDofs(neumannDofIdx)==variableDofIdx2) THEN
                              EXIT
                            END IF
                          END DO
                          DO localNodeIdx1=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_NODES
                            nodeNumber1=topology%lines%lines(lineIdx)%NODES_IN_LINE(localNodeIdx1)
                            DO derivativeIdx1=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx1)
                              derivativeNumber1=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(1,derivativeIdx1,localNodeIdx1)
                              versionNumber1=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(2,derivativeIdx1,localNodeIdx1)
                              variableDofIdx1=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                & NODES(nodeNumber1)%DERIVATIVES(derivativeNumber1)%VERSIONS(versionNumber1)
                              CALL LIST_ITEM_ADD(columnIndicesLists(variableDofIdx1)%PTR,neumannDofIdx,err,error,*999)
                            END DO !derivativeIdx1
                          END DO !localNodeIdx1
                        ELSE IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                          CYCLE linesLoop
                        END IF
                      END DO !derivativeIdx2
                    END DO !localNodeIdx2
                  END IF
                END DO linesLoop 
              CASE(3)
                IF(.NOT.ASSOCIATED(topology%faces)) THEN
                  CALL FlagError("Topology faces have not been calculated.",err,error,*999)
                END IF
                facesLoop: DO faceIdx=1,topology%faces%NUMBER_OF_FACES
                  IF(topology%faces%faces(faceIdx)%BOUNDARY_FACE) THEN
                    DO localNodeIdx2=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_NODES
                      nodeNumber2=topology%faces%faces(faceIdx)%NODES_IN_FACE(localNodeIdx2)
                      DO derivativeIdx2=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx2)
                        derivativeNumber2=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(1,derivativeIdx2,localNodeIdx2)
                        versionNumber2=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(2,derivativeIdx2,localNodeIdx2)
                        variableDofIdx2=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                          & NODES(nodeNumber2)%DERIVATIVES(derivativeNumber2)%VERSIONS(versionNumber2)
                        IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                          & conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                          !Find the Neumann condition number
                          DO neumannDofIdx=1,neumannConditions%numberOfNeumann
                            IF(neumannConditions%setDofs(neumannDofIdx)==variableDofIdx2) THEN
                              EXIT
                            END IF
                          END DO
                          DO localNodeIdx1=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_NODES
                            nodeNumber1=topology%faces%faces(faceIdx)%NODES_IN_FACE(localNodeIdx1)
                            DO derivativeIdx1=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx1)
                              derivativeNumber1=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(1,derivativeIdx1,localNodeIdx1)
                              versionNumber1=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(2,derivativeIdx1,localNodeIdx1)
                              variableDofIdx1=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                & NODES(nodeNumber1)%DERIVATIVES(derivativeNumber1)%VERSIONS(versionNumber1)
                              CALL LIST_ITEM_ADD(columnIndicesLists(variableDofIdx1)%PTR,neumannDofIdx,err,error,*999)
                            END DO !derivativeIdx1
                          END DO !localNodeIdx1
                        ELSE IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                          CYCLE facesLoop
                        END IF
                      END DO !derivativeIdx2
                    END DO !localNodeIdx2
                  END IF
                END DO facesLoop
              CASE DEFAULT
                CALL FlagError("The dimension is invalid for point Neumann conditions",err,error,*999)
              END SELECT !number of dimensions
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              CALL FlagError("The interpolation type of "// &
                & TRIM(NUMBER_TO_VSTRING(rhsVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(componentIdx,"*",err,error))//".", &
                & err,error,*999)
            END SELECT
          END DO !componentIdx

          DO variableDofIdx=1,rhsVariableMapping%TOTAL_NUMBER_OF_LOCAL
            CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(variableDofIdx)%PTR,err,error,*999)
            CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(variableDofIdx)%PTR,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(variableDofIdx+1)=numberOfNonZeros+1
          ENDDO !variableDofIdx
          
          !Allocate and setup the column locations
          ALLOCATE(columnIndices(numberOfNonZeros),stat=err)
          IF(ERR/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
          DO variableDofIdx=1,rhsVariableMapping%TOTAL_NUMBER_OF_LOCAL
            CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(variableDofIdx)%PTR,numberOfColumns,columns,err,error,*999)
            DO columnIdx=1,numberOfColumns
              columnIndices(rowIndices(variableDofIdx)+columnIdx-1)=columns(columnIdx)
            ENDDO !columnIdx
            DEALLOCATE(columns)
          ENDDO !variableDofIdx
          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(neumannConditions%integrationMatrix, &
            & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
          CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(neumannConditions%integrationMatrix,numberOfNonZeros,err,error,*999)
          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(neumannConditions%integrationMatrix, &
            & rowIndices,columnIndices,err,error,*999)
          DEALLOCATE(rowIndices)
          DEALLOCATE(columnIndices)
        CASE(BOUNDARY_CONDITION_FULL_MATRICES)
          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(neumannConditions%integrationMatrix, &
            & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
        CASE DEFAULT
          CALL FlagError("The Neumann matrix sparsity type of "// &
              & TRIM(NUMBER_TO_VSTRING(boundaryConditionsVariable%BOUNDARY_CONDITIONS%neumannMatrixSparsity,"*",err,error))// &
              & " is invalid.",err,error,*999)
        END SELECT

        CALL DISTRIBUTED_MATRIX_CREATE_FINISH(neumannConditions%integrationMatrix,err,error,*999)

        CALL Field_ParameterSetDataRestore(rhsVariable%FIELD,rhsVariable%VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
          & pointValuesData,err,error,*999)
        CALL DistributedVector_DataRestore(boundaryConditionsVariable%CONDITION_TYPES,conditionTypes,err,error,*999)
        CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(neumannConditions%pointValues,err,error,*999)
      ELSE
        CALL FlagError("The boundary condition Neumann is not associated",err,error,*998)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    END IF

    EXITS("BoundaryConditions_NeumannInitialise")
    RETURN
999 CALL BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
998 ERRORS("BoundaryConditions_NeumannInitialise",err,error)
    EXITS("BoundaryConditions_NeumannInitialise")
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_NeumannInitialise

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition information for a boundary conditions variable
  SUBROUTINE BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise the Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    ENTERS("BoundaryConditions_NeumannFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        IF(ALLOCATED(boundaryConditionsNeumann%setDofs)) &
          & DEALLOCATE(boundaryConditionsNeumann%setDofs)
        IF(ASSOCIATED(boundaryConditionsNeumann%integrationMatrix)) &
          & CALL DISTRIBUTED_MATRIX_DESTROY(boundaryConditionsNeumann%integrationMatrix,err,error,*999)
        IF(ASSOCIATED(boundaryConditionsNeumann%pointValues)) &
          & CALL DISTRIBUTED_VECTOR_DESTROY(boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(boundaryConditionsNeumann%pointDofMapping,err,error,*999)
        DEALLOCATE(boundaryConditionsNeumann)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannFinalise

  !
  !================================================================================================================================
  !

  !>Calculates integrated Neumann condition values from point values for a boundary conditions variable and
  !>updates the FIELD_INTEGRATED_NEUMANN_SET_TYPE parameter set for the field variable.
  SUBROUTINE BoundaryConditions_NeumannIntegrate(rhsBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(IN) :: rhsBoundaryConditions !<The boundary conditions for the right hand side field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: neumannDofIdx,componentIdx,nodeIdx,localNodeIdx1,localNodeIdx2,derivativeIdx1,derivativeIdx2, &
      & nodeNumber1,nodeNumber2,versionNumber1,versionNumber2,derivativeNumber1,derivativeNumber2, &
      & lineIdx,faceIdx,variableDofIdx,variableDofIdx1,variableDofIdx2,gaussIdx,ms,os
    INTEGER(INTG), POINTER :: conditionTypes(:)
    LOGICAL :: dependentGeometry
    REAL(DP) :: integratedValue,phim,phio
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannConditions
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointMetrics(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: geometricInterpolationParameters(:),rhsInterpolationParameters(:)
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: integratedValues
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme

    ENTERS("BoundaryConditions_NeumannIntegrate",err,error,*999)

    neumannConditions=>rhsBoundaryConditions%neumannBoundaryConditions
    !Check that Neumann conditions are associated, otherwise do nothing
    IF(ASSOCIATED(neumannConditions)) THEN
      rhsVariable=>rhsBoundaryConditions%VARIABLE
      IF(.NOT.ASSOCIATED(rhsVariable)) THEN
        CALL FlagError("Field variable for RHS boundary conditions is not associated.",err,error,*999)
      END IF

      NULLIFY(conditionTypes)
      CALL DistributedVector_DataGet(rhsBoundaryConditions%CONDITION_TYPES,conditionTypes,err,error,*999)

      CALL Field_GeometricGeneralFieldGet(rhsVariable%field,geometricField,dependentGeometry,err,error,*999)

      CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(neumannConditions%integrationMatrix,0.0_DP,err,error,*999)

      !Initialise field interpolation parameters for the geometric field, which are required for the
      !face/line Jacobian and scale factors
      NULLIFY(geometricInterpolationParameters)
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,geometricInterpolationParameters,err,error,*999)
      NULLIFY(rhsInterpolationParameters)
      CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(rhsVariable%field,rhsInterpolationParameters,err,error,*999)
      NULLIFY(interpolatedPoints)
      CALL FIELD_INTERPOLATED_POINTS_INITIALISE(geometricInterpolationParameters,interpolatedPoints,err,error,*999)
      NULLIFY(interpolatedPointMetrics)
      CALL Field_InterpolatedPointsMetricsInitialise(interpolatedPoints,interpolatedPointMetrics,err,error,*999)

      DO componentIdx=1,rhsVariable%NUMBER_OF_COMPONENTS
        SELECT CASE(rhsVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          !Get topology for finding faces/lines
          topology=>rhsVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY
          IF(.NOT.ASSOCIATED(topology)) THEN
            CALL FlagError("Field component topology is not associated.",err,error,*999)
          END IF
          IF(.NOT.ASSOCIATED(topology%NODES)) THEN
            CALL FlagError("Topology nodes are not associated.",err,error,*999)
          END IF
          SELECT CASE(rhsVariable%COMPONENTS(componentIdx)%DOMAIN%NUMBER_OF_DIMENSIONS)
          CASE(1)
            DO nodeIdx=1,topology%nodes%NUMBER_OF_NODES
              IF(topology%nodes%nodes(nodeIdx)%BOUNDARY_NODE) THEN
                variableDofIdx=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                IF(conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                  & conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                  !Find the Neumann condition number
                  DO neumannDofIdx=1,neumannConditions%numberOfNeumann
                    IF(neumannConditions%setDofs(neumannDofIdx)==variableDofIdx) THEN
                      EXIT
                    END IF
                  END DO
                  CALL DistributedMatrix_ValuesSet(neumannConditions%integrationMatrix,variableDofIdx,neumannDofIdx, &
                    & 1.0_DP,err,error,*999)
                END IF
              END IF
            END DO !nodeIdx
          CASE(2)
            IF(.NOT.ASSOCIATED(topology%lines)) THEN
              CALL FlagError("Topology lines have not been calculated.",err,error,*999)
            END IF
            linesLoop: DO lineIdx=1,topology%lines%NUMBER_OF_LINES
              IF(topology%lines%lines(lineIdx)%BOUNDARY_LINE) THEN
                DO localNodeIdx2=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_NODES
                  nodeNumber2=topology%lines%lines(lineIdx)%NODES_IN_LINE(localNodeIdx2)
                  DO derivativeIdx2=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx2)
                    derivativeNumber2=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(1,derivativeIdx2,localNodeIdx2)
                    versionNumber2=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(2,derivativeIdx2,localNodeIdx2)
                    variableDofIdx2=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeNumber2)%DERIVATIVES(derivativeNumber2)%VERSIONS(versionNumber2)
                    IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      !Find the Neumann condition number
                      DO neumannDofIdx=1,neumannConditions%numberOfNeumann
                        IF(neumannConditions%setDofs(neumannDofIdx)==variableDofIdx2) THEN
                          EXIT
                        END IF
                      END DO
                      !Now perform actual integration
                      quadratureScheme=>topology%lines%lines(lineIdx)%basis%QUADRATURE% &
                        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                      IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                        CALL FlagError("Line basis default quadrature scheme is not associated.",err,error,*999)
                      END IF
                      CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,lineIdx, &
                        & geometricInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999, &
                        & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                        CALL Field_InterpolationParametersScaleFactorsLineGet(lineIdx, &
                          & rhsInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                      END IF
                      DO localNodeIdx1=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_NODES
                        nodeNumber1=topology%lines%lines(lineIdx)%NODES_IN_LINE(localNodeIdx1)
                        DO derivativeIdx1=1,topology%lines%lines(lineIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx1)
                          derivativeNumber1=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(1,derivativeIdx1,localNodeIdx1)
                          versionNumber1=topology%lines%lines(lineIdx)%DERIVATIVES_IN_LINE(2,derivativeIdx1,localNodeIdx1)
                          variableDofIdx1=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(nodeNumber1)%DERIVATIVES(derivativeNumber1)%VERSIONS(versionNumber1)
                          ms=topology%lines%lines(lineIdx)%basis%ELEMENT_PARAMETER_INDEX(derivativeIdx1,localNodeIdx1)
                          os=topology%lines%lines(lineIdx)%basis%ELEMENT_PARAMETER_INDEX(derivativeIdx2,localNodeIdx2)
                          integratedValue=0.0_DP
                          !Loop over line gauss points, adding gauss weighted terms to the integral
                          DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                              & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_LINE_TYPE, &
                              & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                            !Get basis function values at gauss points
                            phim=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussIdx)
                            phio=quadratureScheme%GAUSS_BASIS_FNS(os,NO_PART_DERIV,gaussIdx)
                            !Add gauss point value to total line integral
                            integratedValue=integratedValue+phim*phio*quadratureScheme%GAUSS_WEIGHTS(gaussIdx)* &
                              & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian
                          END DO
                          !Multiply by scale factors
                          IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                            integratedValue=integratedValue* &
                              & rhsInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,componentIdx)* &
                              & rhsInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(os,componentIdx)
                          END IF
                          !Add integral term to N matrix
                          CALL DISTRIBUTED_MATRIX_VALUES_ADD(neumannConditions%integrationMatrix,variableDofIdx1,neumannDofIdx, &
                            & integratedValue,err,error,*999)
                        END DO !derivativeIdx1
                      END DO !localNodeIdx1
                    ELSE IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE linesLoop
                    END IF
                  END DO !derivativeIdx2
                END DO !localNodeIdx2
              END IF
            END DO linesLoop
          CASE(3)
            IF(.NOT.ASSOCIATED(topology%faces)) THEN
              CALL FlagError("Topology faces have not been calculated.",err,error,*999)
            END IF
            facesLoop: DO faceIdx=1,topology%faces%NUMBER_OF_FACES
              IF(topology%faces%faces(faceIdx)%BOUNDARY_FACE) THEN
                DO localNodeIdx2=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_NODES
                  nodeNumber2=topology%faces%faces(faceIdx)%NODES_IN_FACE(localNodeIdx2)
                  DO derivativeIdx2=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx2)
                    derivativeNumber2=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(1,derivativeIdx2,localNodeIdx2)
                    versionNumber2=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(2,derivativeIdx2,localNodeIdx2)
                    variableDofIdx2=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeNumber2)%DERIVATIVES(derivativeNumber2)%VERSIONS(versionNumber2)
                    IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      !Find the Neumann condition number
                      DO neumannDofIdx=1,neumannConditions%numberOfNeumann
                        IF(neumannConditions%setDofs(neumannDofIdx)==variableDofIdx2) THEN
                          EXIT
                        END IF
                      END DO
                      !Now perform actual integration
                      quadratureScheme=>topology%faces%faces(faceIdx)%basis%QUADRATURE% &
                        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                      IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                        CALL FlagError("Face basis default quadrature scheme is not associated.",err,error,*999)
                      END IF
                      CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceIdx, &
                        & geometricInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999, &
                        & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                        CALL Field_InterpolationParametersScaleFactorsFaceGet(faceIdx, &
                          & rhsInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                      END IF
                      DO localNodeIdx1=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_NODES
                        nodeNumber1=topology%faces%faces(faceIdx)%NODES_IN_FACE(localNodeIdx1)
                        DO derivativeIdx1=1,topology%faces%faces(faceIdx)%basis%NUMBER_OF_DERIVATIVES(localNodeIdx1)
                          derivativeNumber1=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(1,derivativeIdx1,localNodeIdx1)
                          versionNumber1=topology%faces%faces(faceIdx)%DERIVATIVES_IN_FACE(2,derivativeIdx1,localNodeIdx1)
                          variableDofIdx1=rhsVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(nodeNumber1)%DERIVATIVES(derivativeNumber1)%VERSIONS(versionNumber1)
                          ms=topology%faces%faces(faceIdx)%basis%ELEMENT_PARAMETER_INDEX(derivativeIdx1,localNodeIdx1)
                          os=topology%faces%faces(faceIdx)%basis%ELEMENT_PARAMETER_INDEX(derivativeIdx2,localNodeIdx2)
                          integratedValue=0.0_DP
                          !Loop over face gauss points, adding gauss weighted terms to the integral
                          DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                              & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
                              & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                            !Get basis function values at gauss points
                            phim=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussIdx)
                            phio=quadratureScheme%GAUSS_BASIS_FNS(os,NO_PART_DERIV,gaussIdx)
                            !Add gauss point value to total face integral
                            integratedValue=integratedValue+phim*phio*quadratureScheme%GAUSS_WEIGHTS(gaussIdx)* &
                              & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian
                          END DO
                          !Multiply by scale factors
                          IF(rhsVariable%FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                            integratedValue=integratedValue* &
                              & rhsInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,componentIdx)* &
                              & rhsInterpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr%SCALE_FACTORS(os,componentIdx)
                          END IF
                          !Add integral term to N matrix
                          CALL DISTRIBUTED_MATRIX_VALUES_ADD(neumannConditions%integrationMatrix,variableDofIdx1,neumannDofIdx, &
                            & integratedValue,err,error,*999)
                        END DO !derivativeIdx1
                      END DO !localNodeIdx1
                    ELSE IF(conditionTypes(variableDofIdx2)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE facesLoop
                    END IF
                  END DO !derivativeIdx2
                END DO !localNodeIdx2
              END IF
            END DO facesLoop
          CASE DEFAULT
            CALL FlagError("The dimension is invalid for point Neumann conditions",err,error,*999)
          END SELECT !number of dimensions
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          CALL FlagError("The interpolation type of "// &
            & TRIM(NUMBER_TO_VSTRING(rhsVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
            & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(componentIdx,"*",err,error))//".", &
            & err,error,*999)
        END SELECT
      END DO !componentIdx

      CALL DISTRIBUTED_MATRIX_UPDATE_START(neumannConditions%integrationMatrix,err,error,*999)
      CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(neumannConditions%integrationMatrix,err,error,*999)

      NULLIFY(integratedValues)
      CALL FIELD_PARAMETER_SET_VECTOR_GET(rhsVariable%field,rhsVariable%variable_type,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & integratedValues,err,error,*999)
      CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(integratedValues,0.0_DP,err,error,*999)
      !Perform matrix multiplication, f = N q, to calculate force vector from integration matrix and point values
      CALL DISTRIBUTED_MATRIX_BY_VECTOR_ADD(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
        & neumannConditions%integrationMatrix,neumannConditions%pointValues,integratedValues,err,error,*999)

      CALL DistributedVector_DataRestore(rhsBoundaryConditions%CONDITION_TYPES,conditionTypes,err,error,*999)

      IF(DIAGNOSTICS1) THEN
        IF(dependentGeometry) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Using dependent field geometry",err,error,*999)
        ELSE
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Using undeformed geometry",err,error,*999)
        END IF
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,neumannConditions%numberOfNeumann,6,6,neumannConditions%setDofs, &
          & '("  setDofs:",6(X,I8))', '(10X,6(X,I8))',err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point values",err,error,*999)
        CALL DISTRIBUTED_VECTOR_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%pointValues,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann integration matrix",err,error,*999)
        CALL DISTRIBUTED_MATRIX_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%integrationMatrix,err,error,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Integrated values",err,error,*999)
        CALL DISTRIBUTED_VECTOR_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,integratedValues,err,error,*999)
      END IF
    END IF

    EXITS("BoundaryConditions_NeumannIntegrate")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannIntegrate",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannIntegrate

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the Neumann integration matrices
  SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet(boundaryConditions,sparsityType,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The matrix sparsity type to be set \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions

    ENTERS("BoundaryConditions_NeumannSparsityTypeSet",ERR,ERROR,*999)

    IF(ASSOCIATED(boundaryConditions)) THEN
      SELECT CASE(sparsityType)
      CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
        boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
      CASE(BOUNDARY_CONDITION_FULL_MATRICES)
        boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_FULL_MATRICES
      CASE DEFAULT
        CALL FlagError("The specified Neumann integration matrix sparsity type of "// &
          & TRIM(NUMBER_TO_VSTRING(sparsityType,"*",err,error))//" is invalid.",err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Boundary conditions are not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannSparsityTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannSparsityTypeSet",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user node. \see OPENCMISS_CMISSBoundaryConditionsSetNode
  SUBROUTINE BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
    & USER_NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The boundary condition type to set \see BOUNDARY_CONDITIONS_ROUTINES_BoundaryConditions,BOUNDARY_CONDITIONS_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(FIELD_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_FINISHED) THEN
        CALL FlagError("Boundary conditions have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD)) THEN
          CALL FIELD_COMPONENT_DOF_GET_USER_NODE(FIELD,VARIABLE_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
            & USER_NODE_NUMBER,COMPONENT_NUMBER,local_ny,global_ny,ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
              CALL BoundaryConditions_CheckInterpolationType(CONDITION,FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,ERR,ERROR,*999)
              CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD,VARIABLE_TYPE, &
                & local_ny,CONDITION,VALUE,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The boundary conditions for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("The dependent field variable is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_SET_NODE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_SET_NODE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_SET_NODE

  !
  !================================================================================================================================
  !

  !Applies the constraint to the field variable.
  SUBROUTINE BoundaryConditions_ConstraintsApply(boundaryConditionsVariable,fieldSetType,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(IN) :: boundaryConditionsVariable !<The boundary conditions variable
    INTEGER(INTG), INTENT(IN) :: fieldSetType !The field parameter set type to calculate the constrained values for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfNonZeros
    REAL(DP), POINTER :: constrainedValuesData(:)
    TYPE(BoundaryConditionsConstraintsType), POINTER :: constraints
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: constraintMatrix
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: constraintVector,constrainedValues,variableVector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_ConstraintsApply",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    fieldVariable=>boundaryConditionsVariable%variable
    IF(.NOT.ASSOCIATED(fieldvariable)) &
      & CALL FlagError("Field variable is not associated.",err,error,*999)
    constraints=>boundaryConditionsVariable%constraints
    IF(.NOT.ASSOCIATED(constraints)) &
      & CALL FlagError("Constraints is not associated.",err,error,*999)
    constraintMatrix=>constraints%constraintMatrix
    IF(.NOT.ASSOCIATED(constraintMatrix)) &
      & CALL FlagError("Constraint matrix is not associated.",err,error,*999)
    constraintVector=>constraints%constraintVector
    IF(.NOT.ASSOCIATED(constraintVector)) &
      & CALL FlagError("Constraint vector is not associated.",err,error,*999)
    constrainedValues=>constraints%constrainedValues
    IF(.NOT.ASSOCIATED(constrainedValues)) &
      & CALL FlagError("Constraint values is not associated.",err,error,*999)

    !Check whether there are any linear constraints at all, so we can skip this.
    CALL DistributedMatrix_NumberNonZerosGet(constraintMatrix,numberOfNonZeros,err,error,*999)

    IF(numberOfNonZeros>0) THEN
      !Calculate the constrained values, e.g. u_constr=Pu+c.
      !Copy the non-homogeneous constraints (for example the Dirichlet boundary conditions) to the constrained values.
      CALL DistributedVector_Copy(constraintVector,constrainedValues,1.0_DP,err,error,*999)
      NULLIFY(variableVector)
      CALL Field_ParameterSetVectorGet(fieldVariable%field,fieldVariable%VARIABLE_TYPE,fieldSetType,variableVector,err,error,*999)
      !The constrained dofs can not be coupled to dofs on another domain that are not in this domain's ghost dofs. Therefore, if we
      !include the ghost rows we do not have to update the parameter set later on.
      CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE,1.0_DP,constraintMatrix, &
        & variableVector,constrainedValues,err,error,*999)

      !Update the field parameter set with the constrained values.
      NULLIFY(constrainedValuesData)
      CALL DistributedVector_DataGet(constrainedValues,constrainedValuesData,err,error,*999)
      CALL Field_ParameterSetUpdateLocalDOF(fieldVariable%field,fieldVariable%VARIABLE_TYPE,fieldSetType, &
        & constraints%constrainedDofIndices,constrainedValuesData,err,error,*999)
      CALL DistributedVector_DataRestore(constrainedValues,constrainedValuesData,err,error,*999)
      !No parameter set update needed, as explained above.
    END IF

    EXITS("BoundaryConditions_ConstraintsApply")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstraintsApply",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstraintsApply

  !
  !================================================================================================================================
  !

  !>Constrain multiple equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainDofsEqual(boundaryConditions,fieldVariable,coupledDofs,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: coupledDofs(:) !<The coupled DOFs to be constrained to be equal.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDofs,dofIdx,dofIdx2
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable

    ENTERS("BoundaryConditions_ConstrainDofsEqual",err,error,*999)

    numberOfDofs=SIZE(coupledDofs,1)
    IF(numberOfDofs<2) THEN
      CALL FlagError("Cannot constrain zero or 1 DOF to be equal.",err,error,*999)
    END IF

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDofs
      DO dofIdx2=dofIdx+1,numberOfDofs
        IF(coupledDofs(dofIdx)==coupledDofs(dofIdx2)) THEN
          CALL FlagError("DOF number "//TRIM(NumberToVstring(coupledDofs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOFs constrained to be equal.",err,error,*999)
        END IF
      END DO
    END DO

    NULLIFY(boundaryConditionsVariable)
    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)

    !Add new DOF constraints
    !We set all DOFs except the first to be equal to 1.0 * the first DOF
    !The first DOF is left unconstrained
    DO dofIdx=2,numberOfDofs
      CALL BoundaryConditions_DofConstraintSet(boundaryConditionsVariable,coupledDofs(dofIdx), &
        & [coupledDofs(1)],[1.0_DP],err,error,*999)
      !Set the DOF type and BC type of the constrained DOF
      CALL BoundaryConditions_SetConditionType(boundaryConditionsVariable,coupledDofs(dofIdx), &
        & BOUNDARY_CONDITION_LINEAR_CONSTRAINT,err,error,*999)
    END DO

    EXITS("BoundaryConditions_ConstrainDofsEqual")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstrainDofsEqual",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstrainDofsEqual

  !
  !================================================================================================================================
  !

  !>Constrain multiple nodal equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainNodeDofsEqual(boundaryConditions,field,fieldVariableType,versionNumber, &
      & derivativeNumber,component,nodes,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER, INTENT(IN) :: boundaryConditions !<The solver equations boundary conditions to constrain the DOFs for.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The equations dependent field containing the field DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: fieldVariableType !<The field variable type of the DOFs to be constrained. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version number.
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number.
    INTEGER(INTG), INTENT(IN) :: component !<The field component number of the DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: nodes(:) !<The user numbers of the nodes to be constrained to be equal.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    LOGICAL :: isLocalConstraint,localNodeExists,isGhostNode,isInternalDof
    INTEGER(INTG) :: numberOfNodes,nodeIdx,globalDof,numberOfValidDomains,adjacentDomainIdx,tempDof, &
      & localNodeNumber(SIZE(nodes)),coupledDofs(SIZE(nodes)),MPI_IERROR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_ConstrainNodeDofsEqual",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)

    numberOfNodes=SIZE(nodes,1)
    !Get the field variable for the field
    NULLIFY(fieldVariable)
    CALL FIELD_VARIABLE_GET(field,fieldVariableType,fieldVariable,err,error,*999)
    
    !Make sure that constrained node dofs are all on the same domain. This is to make sure that the off-processor part  
    !of the constraint-matrix is zero (ghost columns are allowed).
    isLocalConstraint=.TRUE.
    numberOfValidDomains=1
    DO nodeIdx=1,numberOfNodes
      CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(fieldVariable%components(component)%domain%topology,nodes(nodeIdx),localNodeExists, &
        & localNodeNumber(nodeIdx),isGhostNode,err,error,*999)
      IF(.NOT.localNodeExists) THEN
        isLocalConstraint=.FALSE.
        numberOfValidDomains=0
        EXIT
      END IF
    END DO
    !Count the number of domains for which all the given nodes are local. If this number is greater than 1 it means that all the
    !constrained nodes are shared with other domains.
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,numberOfValidDomains,1,MPI_INTEGER,MPI_SUM,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,err,error,*999)
    IF(numberOfValidDomains<1) CALL FlagError("All constrained node dofs are not on the same domain.",err,error,*999)
    
    IF(isLocalConstraint) THEN
      !Get field DOFs for nodes
      DO nodeIdx=1,numberOfNodes
        CALL FIELD_COMPONENT_DOF_GET_USER_NODE(field,fieldVariableType,versionNumber,derivativeNumber,nodes(nodeIdx), &
          & component,coupledDofs(nodeIdx),globalDof,err,error,*999)
      END DO !nodeIdx
      IF(numberOfValidDomains==1) THEN
        !Make sure that the constrained node dof (the first one) is not a ghost node dof on another domain that is coupled to
        !another node dof that is not on the same domain as that ghost node dof. If so, take another node dof to be the constrained one.
        domainMapping=>fieldVariable%DOMAIN_MAPPING
        IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Field variable domain mapping is not associated.",err,error,*999)
        DO nodeIdx=1,numberOfNodes
          isInternalDof=.TRUE.
          DO adjacentDomainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
            IF(ANY(domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES==coupledDofs(nodeIdx))) THEN
              isInternalDof=.FALSE.
              EXIT
            END IF
          END DO !adjacentDomainIdx
          IF(isInternalDof) EXIT
        END DO !nodeIdx
      END IF
      !Swap the first dof and the nodeIdx'th dof
      tempDof=coupledDofs(1)
      coupledDofs(1)=coupledDofs(nodeIdx)
      coupledDofs(nodeIdx)=tempDof

      !Now set DOF constraint
      CALL BoundaryConditions_ConstrainDofsEqual(boundaryConditions,fieldVariable,coupledDofs,err,error,*999)
    END IF

    EXITS("BoundaryConditions_ConstrainNodeDofsEqual")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstrainNodeDofsEqual",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstrainNodeDofsEqual

  !
  !================================================================================================================================
  !

  !>Finalise a boundary condition dof constraint.
  SUBROUTINE BoundaryConditions_DofConstraintFinalise(dofConstraint,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables

    ENTERS("BoundaryConditions_DofConstraintFinalise",err,error,*999)

    IF(ASSOCIATED(dofConstraint)) THEN
      IF(ALLOCATED(dofConstraint%coupledDofIndices)) DEALLOCATE(dofConstraint%coupledDofIndices)
      IF(ALLOCATED(dofConstraint%coupledDofCoefficients)) DEALLOCATE(dofConstraint%coupledDofCoefficients)
    END IF

    EXITS("BoundaryConditions_DofConstraintFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DofConstraintFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_DofConstraintFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a boundary condition dof constraint.
  SUBROUTINE BoundaryConditions_DofConstraintInitialise(dofConstraint,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables

    ENTERS("BoundaryConditions_DofConstraintInitialise",err,error,*999)

    IF(ASSOCIATED(dofConstraint)) THEN
      CALL FlagError("Boundary conditions dof constraint is already associated.",err,error,*999)
    END IF

    ALLOCATE(dofConstraint,stat=err)
    IF(err/=0) CALL FlagError("Could not allocate boundary conditions dof constraint.",err,error,*999)
    dofConstraint%constrainedDofIndex=0
    dofConstraint%numberOfCoupledDofs=0

    EXITS("BoundaryConditions_DofConstraintInitialise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DofConstraintInitialise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_DofConstraintInitialise

  !
  !================================================================================================================================
  !

  !>Constrain a DOF to be a linear combination of other DOFs.
  SUBROUTINE BoundaryConditions_DofConstraintSet(boundaryConditionsVariable,constrainedDof,coupledDofIndices, &
      & coupledDofCoefficients,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(IN) :: boundaryConditionsVariable !<The boundary conditions variable for which to constrain the DOF.
    INTEGER(INTG), INTENT(IN) :: constrainedDof !<The constrained DOF to set the constraint on.
    INTEGER(INTG), INTENT(IN) :: coupledDofIndices(:) !<The DOFs that this DOF is constrained to depend on.
    REAL(DP), INTENT(IN) :: coupledDofCoefficients(:) !<The coefficient values in the DOF constraint.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfCoupledDofs,constraintDofIdx,newConstraintDofIdx,dofIdx,dofIdx2
    TYPE(BoundaryConditionsCreateValuesCacheType), POINTER :: createValuesCache 
    TYPE(BoundaryConditionsDofConstraintPtrType), ALLOCATABLE :: newDofConstraints(:)

    ENTERS("BoundaryConditions_DofConstraintSet",err,error,*998)

    !Check pointers for association
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
      CALL FlagError("Boundary conditions variable are not associated.",err,error,*998)
    END IF
    createValuesCache=>boundaryConditionsVariable%createValuesCache
    IF(.NOT.ASSOCIATED(createValuesCache)) THEN
      CALL FlagError("Boundary conditions variable create values cache not associated.",err,error,*998)
    END IF

    numberOfCoupledDofs=SIZE(coupledDofIndices,1)
    IF(numberOfCoupledDofs/=SIZE(coupledDofCoefficients,1)) THEN
      CALL FlagError("Length of coupledDofCoefficients does not match length of dofs array.",err,error,*998)
    END IF

    !Important: when setting constraints between local dofs, this should be done on all domains, even if dofs are ghosts. This way we
    !avoid having to communicate constraints and we can avoid the situation that a dof that is ghost is linearly dependent on a dof
    !that is not local anymore. In other words: a local dof can only be linearly dependent on other local dofs (including ghosts).

    !Check for duplicate dofs
    DO dofIdx=1,numberOfCoupledDofs
      DO dofIdx2=dofIdx+1,numberOfCoupledDofs
        IF(coupledDofIndices(dofIdx)==coupledDofIndices(dofIdx2)) THEN
          CALL FlagError("Dof number "//TRIM(NumberToVstring(coupledDofIndices(dofIdx),"*",err,error))// &
            & " is duplicated in the dof constraint.",err,error,*998)
        END IF
      END DO
    END DO

    !Check that the new constrained dof is not coupled to another constrained dof.
    DO dofIdx=1,numberOfCoupledDofs
      DO constraintDofIdx=1,createValuesCache%numberOfConstraints
        IF(coupledDofIndices(dofIdx)==createValuesCache%dofConstraints(constraintDofIdx)%ptr%constrainedDofIndex) THEN
          CALL FlagError("Dof number "//TRIM(NumberToVstring(coupledDofIndices(dofIdx),"*",err,error))// &
            & " is already constrained.",err,error,*998)
        END IF
      END DO
    END DO

    !Allocate new dof constraints, copy old dof constraints and insert the new one such that the constrained dof indices are in
    !ascending order. 
    ALLOCATE(newDofConstraints(createValuesCache%numberOfConstraints+1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new dof constraints array.",err,error,*998)
    DO constraintDofIdx=1,createValuesCache%numberOfConstraints
      IF(constrainedDof>createValuesCache%dofConstraints(constraintDofIdx)%ptr%constrainedDofIndex) THEN
        newDofConstraints(constraintDofIdx)%ptr=>createValuesCache%dofConstraints(constraintDofIdx)%ptr
      ELSE IF(constrainedDof<createValuesCache%dofConstraints(constraintDofIdx)%ptr%constrainedDofIndex) THEN
        EXIT
      ELSE
        CALL FlagError("Dof number "//TRIM(NumberToVstring(constraintDofIdx,"*",err,error))// &
          & " was already added to the dof constraints.",err,error,*999)
      END IF
    END DO !constraintDofIdx

    !Create the new dof constraint.
    NULLIFY(newDofConstraints(constraintDofIdx)%ptr)
    CALL BoundaryConditions_DofConstraintInitialise(newDofConstraints(constraintDofIdx)%ptr,err,error,*999)
    newDofConstraints(constraintDofIdx)%ptr%constrainedDofIndex=constrainedDof
    newDofConstraints(constraintDofIdx)%ptr%numberOfCoupledDofs=numberOfCoupledDofs
    ALLOCATE(newDofConstraints(constraintDofIdx)%ptr%coupledDofIndices(numberOfCoupledDofs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new dof constraint coupled dof indices array.",err,error,*999)
    newDofConstraints(constraintDofIdx)%ptr%coupledDofIndices=coupledDofIndices
    ALLOCATE(newDofConstraints(constraintDofIdx)%ptr%coupledDofCoefficients(numberOfCoupledDofs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new dof constraint coupled dof coefficients array.",err,error,*999)
    newDofConstraints(constraintDofIdx)%ptr%coupledDofCoefficients=coupledDofCoefficients
    
    !Complete copying the old dof constraints (if necessary).
    newConstraintDofIdx=constraintDofIdx
    DO constraintDofIdx=newConstraintDofIdx,createValuesCache%numberOfConstraints
      newDofConstraints(constraintDofIdx+1)%ptr=>createValuesCache%dofConstraints(constraintDofIdx)%ptr
    END DO !constraintDofIdx

    CALL MOVE_ALLOC(newDofConstraints,createValuesCache%dofConstraints)
    createValuesCache%numberOfConstraints=createValuesCache%numberOfConstraints+1
    
    EXITS("BoundaryConditions_DofConstraintSet")
    RETURN
999 IF(ALLOCATED(newDofConstraints)) THEN
      DO constraintDofIdx=1,SIZE(newDofConstraints)
        CALL BoundaryConditions_DofConstraintFinalise(newDofConstraints(constraintDofIdx)%ptr,err,error,*998)
      END DO !constraintDofIdx
      DEALLOCATE(newDofConstraints)
    END IF
998 ERRORSEXITS("BoundaryConditions_DofConstraintSet",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_DofConstraintSet

  !
  !================================================================================================================================
  !

  !>Finish the creation of the dof constraints for a boundary conditions variable
  SUBROUTINE BoundaryConditions_ConstraintsCreateFinish(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to boundary conditions variable to finish the dof constraints for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintDofIdx,i,adjacentDomainIdx,ghostSendIdx,ghostReceiveIdx,localDof
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:),columnIndices(:)
    TYPE(BoundaryConditionsConstraintsType), POINTER :: constraints
    TYPE(BoundaryConditionsCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: constrainedDofsMapping,variableDomainMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(LIST_TYPE), POINTER :: ghostSendList,ghostReceiveList 

    ENTERS("BoundaryConditions_ConstraintsCreateFinish",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      fieldVariable=>boundaryConditionsVariable%variable
      IF(ASSOCIATED(fieldVariable)) THEN
        constraints=>boundaryConditionsVariable%constraints
        IF(.NOT.ASSOCIATED(constraints)) THEN
          CALL FlagError("Boundary conditions DOF constraints is not associated.",err,error,*998)
        END IF
        createValuesCache=>boundaryConditionsVariable%createValuesCache
        IF(.NOT.ASSOCIATED(createValuesCache)) THEN
          CALL FlagError("Boundary conditions create values cache is not associated.",err,error,*998)
        END IF
        variableDomainMapping=>fieldVariable%domain_mapping
        IF(.NOT.ASSOCIATED(variableDomainMapping)) THEN
          CALL FlagError("Field variable domain mapping is not associated for variable type "// &
            & TRIM(NumberToVstring(fieldVariable%variable_type,"*",err,error))//".",err,error,*998)
        END IF
        
        constraints%numberOfConstraints=createValuesCache%numberOfConstraints

        !Create a domain mapping for the constrained dofs
        ALLOCATE(constraints%constrainedDofsMapping,stat=err)
        IF(err/=0) CALL FlagError("Could not allocate constrained dofs domain mapping.",err,error,*999)
        constrainedDofsMapping=>constraints%constrainedDofsMapping
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(constrainedDofsMapping,variableDomainMapping%NUMBER_OF_DOMAINS, &
          & err,error,*999)

        ALLOCATE(constraints%constrainedDofIndices(constraints%numberOfConstraints),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate constrained dofs domain mapping.",err,error,*999)
        DO constraintDofIdx=1,constraints%numberOfConstraints
          constraints%constrainedDofIndices(constraintDofIdx)=createValuesCache%dofConstraints(constraintDofIdx)%ptr% &
            & constrainedDofIndex
        END DO !constraintDofIdx

        constrainedDofsMapping%NUMBER_OF_LOCAL=COUNT(constraints%constrainedDofIndices<=variableDomainMapping%NUMBER_OF_LOCAL)
        constrainedDofsMapping%TOTAL_NUMBER_OF_LOCAL=constraints%numberOfConstraints

        constrainedDofsMapping%NUMBER_OF_ADJACENT_DOMAINS=variableDomainMapping%NUMBER_OF_ADJACENT_DOMAINS
        ALLOCATE(constrainedDofsMapping%ADJACENT_DOMAINS(constrainedDofsMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
        IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)

        DO adjacentDomainIdx=1,constrainedDofsMapping%NUMBER_OF_ADJACENT_DOMAINS
          CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(constrainedDofsMapping%ADJACENT_DOMAINS(adjacentDomainIdx), &
            & err,error,*999)
          constrainedDofsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER= &
            & variableDomainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
          NULLIFY(ghostSendList)
          CALL LIST_CREATE_START(ghostSendList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(ghostSendList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_INITIAL_SIZE_SET(ghostSendList, &
            & MAX(constrainedDofsMapping%TOTAL_NUMBER_OF_LOCAL-constrainedDofsMapping%NUMBER_OF_LOCAL,1),err,error,*999)
          CALL LIST_CREATE_FINISH(ghostSendList,err,error,*999)
          NULLIFY(ghostReceiveList)
          CALL LIST_CREATE_START(ghostReceiveList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(ghostReceiveList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_INITIAL_SIZE_SET(ghostReceiveList, & 
            & MAX(constrainedDofsMapping%TOTAL_NUMBER_OF_LOCAL-constrainedDofsMapping%NUMBER_OF_LOCAL,1),err,error,*999)
          CALL LIST_CREATE_FINISH(ghostReceiveList,err,error,*999)
          DO ghostSendIdx=1,variableDomainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS
            localDof=variableDomainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES(ghostSendIdx)
            DO constraintDofIdx=1,constraints%numberOfConstraints
              IF(constraints%constrainedDofIndices(constraintDofIdx)==localDof) THEN
                CALL LIST_ITEM_ADD(ghostSendList,constraintDofIdx,err,error,*999)
                EXIT
              END IF
            END DO !constraintDofIdx
          END DO !ghostSendIdx
          DO ghostReceiveIdx=1,variableDomainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
            localDof=variableDomainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(ghostReceiveIdx)
            DO constraintDofIdx=1,constraints%numberOfConstraints
              IF(constraints%constrainedDofIndices(constraintDofIdx)==localDof) THEN
                CALL LIST_ITEM_ADD(ghostReceiveList,constraintDofIdx,err,error,*999)
                EXIT
              END IF
            END DO !constraintDofIdx
          END DO !ghostReceiveIdx
          CALL LIST_REMOVE_DUPLICATES(ghostSendList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(ghostSendList, &
            & constrainedDofsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS, &
            & constrainedDofsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES,err,error,*999)
          CALL LIST_REMOVE_DUPLICATES(ghostReceiveList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(ghostReceiveList, &
            & constrainedDofsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS, &
            & constrainedDofsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES,err,error,*999)
        ENDDO !adjacentDomainIdx

        !Calculate the local to global map.
        CALL DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE(constrainedDofsMapping,err,error,*999)

        !Set up the constraint matrix
        CALL DistributedMatrix_CreateStart(constrainedDofsMapping,variableDomainMapping,constraints%constraintMatrix, &
          & err,error,*999)
        ALLOCATE(rowIndices(constrainedDofsMapping%TOTAL_NUMBER_OF_LOCAL+1),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        ALLOCATE(columnIndices(constrainedDofsMapping%TOTAL_NUMBER_OF_LOCAL),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
        rowIndices(1)=1
        DO constraintDofIdx=1,constraints%numberOfConstraints
          dofConstraint=>createValuesCache%dofConstraints(constraintDofIdx)%ptr
          rowIndices(constraintDofIdx+1)=rowIndices(constraintDofIdx)+dofConstraint%numberOfCoupledDofs
          IF(dofConstraint%numberOfCoupledDofs>0) THEN
            columnIndices(rowIndices(constraintDofIdx):rowIndices(constraintDofIdx+1)-1)=dofConstraint%coupledDofIndices
          END IF
        END DO !constraintDofIdx
        
        CALL DistributedMatrix_StorageTypeSet(constraints%constraintMatrix,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
          & err,error,*999)
        CALL DistributedMatrix_NumberNonZerosSet(constraints%constraintMatrix,rowIndices(constraints%numberOfConstraints+1)-1, &
          & err,error,*999)
        CALL DistributedMatrix_StorageLocationsSet(constraints%constraintMatrix,rowIndices,columnIndices,err,error,*999)
        CALL DistributedMatrix_CreateFinish(constraints%constraintMatrix,err,error,*999)

        !Set the constraint matrix values
        DO constraintDofIdx=1,constraints%numberOfConstraints
          dofConstraint=>createValuesCache%dofConstraints(constraintDofIdx)%ptr
          IF(dofConstraint%numberOfCoupledDofs>0) THEN
            CALL DistributedMatrix_ValuesSet(constraints%constraintMatrix, &
              & [(constraintDofIdx,i=1,dofConstraint%numberOfCoupledDofs)], &
              & dofConstraint%coupledDofIndices,dofConstraint%coupledDofCoefficients,err,error,*999)
          END IF
        END DO !constraintDofIdx

        IF(ALLOCATED(rowIndices)) DEALLOCATE(rowIndices)
        IF(ALLOCATED(columnIndices)) DEALLOCATE(columnIndices)

        !Set up the constraint vector 
        CALL DistributedVector_CreateStart(constrainedDofsMapping,constraints%constraintVector,err,error,*999)
        CALL DistributedVector_CreateFinish(constraints%constraintVector,err,error,*999)
        CALL DistributedVector_AllValuesSet(constraints%constraintVector,0.0_DP,err,error,*999)
        !Set up the constrained values vector 
        CALL DistributedVector_CreateStart(constrainedDofsMapping,constraints%constrainedValues,err,error,*999)
        CALL DistributedVector_CreateFinish(constraints%constrainedValues,err,error,*999)
        CALL DistributedVector_AllValuesSet(constraints%constrainedValues,0.0_DP,err,error,*999)
      ELSE
        CALL FlagError("Field variable is not associated for this boundary conditions variable",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_ConstraintsCreateFinish")
    RETURN
999 CALL BoundaryConditions_ConstraintsFinalise(constraints,err,error,*998)
998 ERRORS("BoundaryConditions_ConstraintsCreateFinish",err,error)
    EXITS("BoundaryConditions_ConstraintsCreateFinish")
    RETURN 1

  END SUBROUTINE BoundaryConditions_ConstraintsCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise the DOF constraints structure
  SUBROUTINE BoundaryConditions_ConstraintsFinalise(constraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsConstraintsType), POINTER :: constraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("BoundaryConditions_ConstraintsFinalise",err,error,*999)

    IF(ASSOCIATED(constraints)) THEN
      IF(ALLOCATED(constraints%constrainedDofIndices)) DEALLOCATE(constraints%constrainedDofIndices)
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(constraints%constrainedDofsMapping,err,error,*999)
      CALL DistributedMatrix_Destroy(constraints%constraintMatrix,err,error,*999)
      CALL DistributedVector_Destroy(constraints%constraintVector,err,error,*999)
      CALL DistributedVector_Destroy(constraints%constrainedValues,err,error,*999)
    ELSE
      CALL FlagError("constraints pointer is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_ConstraintsFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstraintsFinalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_ConstraintsFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the DOF constraints structure
  SUBROUTINE BoundaryConditions_ConstraintsInitialise(constraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsConstraintsType), POINTER :: constraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("BoundaryConditions_ConstraintsInitialise",err,error,*999)

    IF(ASSOCIATED(constraints)) THEN
      constraints%numberOfConstraints=0
      NULLIFY(constraints%constraintMatrix)
      NULLIFY(constraints%constraintVector)
      NULLIFY(constraints%constrainedValues)
    ELSE
      CALL FlagError("constraints pointer is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_ConstraintsInitialise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstraintsInitialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_ConstraintsInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions variable and deallocate all memory.
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE !<A pointer to the boundary conditions variable to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES))  &
        & CALL DistributedVector_Destroy(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES,err,error,*999) 
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES))  &
        & CALL DistributedVector_Destroy(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,err,error,*999) 
      CALL BoundaryConditions_NeumannFinalise(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) &
        & DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)
      IF(ASSOCIATED(boundary_conditions_variable%constraints)) THEN
        CALL BoundaryConditions_ConstraintsFinalise(boundary_conditions_variable%constraints,err,error,*999)
        DEALLOCATE(boundary_conditions_variable%constraints)
      END IF
      IF(ALLOCATED(boundary_conditions_variable%anyGlobalDofs)) DEALLOCATE(boundary_conditions_variable%anyGlobalDofs)
      DEALLOCATE(BOUNDARY_CONDITIONS_VARIABLE)
    ENDIF
       
    EXITS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_VARIABLE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the boundary conditions variable for a variable type if that variable has not already been initialised, otherwise do nothing.
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_INITIALISE(BOUNDARY_CONDITIONS,FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to initialise a variable type for.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to initialise the boundary conditions variable for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variableIdx
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_PTR_TYPE), ALLOCATABLE :: NEW_BOUNDARY_CONDITIONS_VARIABLES(:)
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE

    ENTERS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        VARIABLE_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
        IF(ASSOCIATED(VARIABLE_DOMAIN_MAPPING)) THEN
          NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
          !Check if boundary conditions variable has already been added, if so then we don't do anything as different equations
          !sets can have the same dependent field variables and will both want to add the variable
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
            ALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES(BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new boundary conditions variables array.",ERR,ERROR,*998)
            IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
              DO variableIdx=1,BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES
                NEW_BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR=> &
                    & BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR
              ENDDO
            ENDIF

            ALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES(BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate boundary condition variable.",ERR,ERROR,*998)
            BOUNDARY_CONDITIONS_VARIABLE=>NEW_BOUNDARY_CONDITIONS_VARIABLES( &
                & BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1)%PTR
            BOUNDARY_CONDITIONS_VARIABLE%BOUNDARY_CONDITIONS=>BOUNDARY_CONDITIONS
            BOUNDARY_CONDITIONS_VARIABLE%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            BOUNDARY_CONDITIONS_VARIABLE%VARIABLE=>FIELD_VARIABLE
            CALL DistributedVector_CreateStart(FIELD_VARIABLE%DOMAIN_MAPPING,BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES, &
              & err,error,*999)
            CALL DistributedVector_DataTypeSet(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES,DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE, &
              & err,error,*999)
            CALL DistributedVector_CreateFinish(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES,err,error,*999)
            CALL DistributedVector_AllValuesSet(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES,BOUNDARY_CONDITION_FREE, &
              & err,error,*999)
            CALL DistributedVector_CreateStart(FIELD_VARIABLE%DOMAIN_MAPPING,BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES, &
              & err,error,*999)
            CALL DistributedVector_DataTypeSet(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE, &
              & err,error,*999)
            CALL DistributedVector_CreateFinish(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,err,error,*999)
            CALL DistributedVector_AllValuesSet(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,BOUNDARY_CONDITION_DOF_FREE, &
              & err,error,*999)
            ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%anyGlobalDofs(MAX_BOUNDARY_CONDITION_NUMBER),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate boundary condition DOF counts array.",ERR,ERROR,*999)
            BOUNDARY_CONDITIONS_VARIABLE%anyGlobalDofs=.FALSE.
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions)
            NULLIFY(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)

            CALL MOVE_ALLOC(NEW_BOUNDARY_CONDITIONS_VARIABLES,BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)
            BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES= &
                & BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES+1

            CALL BoundaryConditions_CreateValuesCacheInitialise(boundary_conditions_variable,ERR,ERROR,*999)

            ALLOCATE(boundary_conditions_variable%constraints,stat=err)
            IF(err/=0) CALL FlagError("Could not allocate boundary conditions dof constraints.",err,error,*999)
            CALL BoundaryConditions_ConstraintsInitialise(boundary_conditions_variable%constraints,err,error,*999)
          END IF
        ELSE
          CALL FlagError("Field variable domain mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FlagError("Field variable is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*998)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE")
    RETURN
999 CALL BOUNDARY_CONDITIONS_VARIABLE_FINALISE(BOUNDARY_CONDITIONS_VARIABLE,DUMMY_ERR,DUMMY_ERROR,*998)
    DEALLOCATE(NEW_BOUNDARY_CONDITIONS_VARIABLES)
998 ERRORSEXITS("BOUNDARY_CONDITIONS_VARIABLE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions variable for a given field variable
  SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to initialise a variable type for.
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to initialise the boundary conditions variable for.
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER, INTENT(OUT) :: BOUNDARY_CONDITIONS_VARIABLE !<On return, a pointer to the boundary conditions variable, or NULL if it wasn't found
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    LOGICAL :: VARIABLE_FOUND

    ENTERS("BOUNDARY_CONDITIONS_VARIABLE_GET",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        IF(ALLOCATED(BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES)) THEN
          VARIABLE_FOUND=.FALSE.
          variableIdx=1
          DO WHILE(variableIdx<=BOUNDARY_CONDITIONS%NUMBER_OF_BOUNDARY_CONDITIONS_VARIABLES.AND..NOT.VARIABLE_FOUND)
            VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR%VARIABLE
            IF(ASSOCIATED(VARIABLE)) THEN
              IF(VARIABLE%VARIABLE_TYPE==FIELD_VARIABLE%VARIABLE_TYPE.AND. &
                & VARIABLE%FIELD%USER_NUMBER==FIELD_VARIABLE%FIELD%USER_NUMBER) THEN
                IF(ASSOCIATED(VARIABLE%FIELD%REGION)) THEN
                  IF(VARIABLE%FIELD%REGION%USER_NUMBER==FIELD_VARIABLE%FIELD%REGION%USER_NUMBER) THEN
                    VARIABLE_FOUND=.TRUE.
                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR
                  ENDIF
                ELSEIF(ASSOCIATED(VARIABLE%FIELD%INTERFACE)) THEN
                  IF(VARIABLE%FIELD%INTERFACE%USER_NUMBER==FIELD_VARIABLE%FIELD%INTERFACE%USER_NUMBER) THEN
                    VARIABLE_FOUND=.TRUE.
                    BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLES(variableIdx)%PTR
                  ENDIF
                ENDIF
              ENDIF
              variableIdx=variableIdx+1
            ENDIF
          ENDDO
        ENDIF
      ELSE
        CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_VARIABLE_GET")
    RETURN
999 ERRORSEXITS("BOUNDARY_CONDITIONS_VARIABLE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE BOUNDARY_CONDITIONS_VARIABLE_GET

  !
  !================================================================================================================================
  !

  !>Initialises the pressure incremented boundary condition.
  SUBROUTINE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE(BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*)
    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE !<A pointer to the boundary conditions variable to initialise a boundary conditions pressure incremented type for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: pressureIncrementedDofIdx,variableDofIdx
    INTEGER(INTG), POINTER :: conditionTypes(:)
    TYPE(BoundaryConditionsCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE), POINTER :: pressureIncrementedConditions 
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable

    ENTERS("BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) THEN
        CALL FlagError("Pressure incremented boundary conditions are already associated for this boundary conditions variable." &
           & ,ERR,ERROR,*999)
      ELSE
        createValuesCache=>BOUNDARY_CONDITIONS_VARIABLE%createValuesCache
        IF(.NOT.ASSOCIATED(createValuesCache)) &
          & CALL FlagError("Boundary condition variable create values cache is not associated.",err,error,*999)
        rhsVariable=>BOUNDARY_CONDITIONS_VARIABLE%variable
        IF(.NOT.ASSOCIATED(rhsVariable)) &
          & CALL FlagError("RHS boundary conditions variable field variable is not associated.",err,error,*999)
        ALLOCATE(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate Pressure incremented Boundary Conditions",ERR,ERROR,*999)
        pressureIncrementedConditions=>BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
        pressureIncrementedConditions%numberOfPressureIncrementedDofs= &
          & createValuesCache%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        ALLOCATE(pressureIncrementedConditions%PRESSURE_INCREMENTED_DOF_INDICES( &
          & pressureIncrementedConditions%numberOfPressureIncrementedDofs),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate Pressure incremented DOF indices array",ERR,ERROR,*999)
        
        NULLIFY(conditionTypes)
        CALL DistributedVector_DataGet(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES,conditionTypes,err,error,*999)
        pressureIncrementedDofIdx=0
        DO variableDofIdx=1,rhsVariable%TOTAL_NUMBER_OF_DOFS
          IF(conditionTypes(variableDofIdx)==BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
            pressureIncrementedDofIdx=pressureIncrementedDofIdx+1
            pressureIncrementedConditions%PRESSURE_INCREMENTED_DOF_INDICES(pressureIncrementedDofIdx)= &
              & variableDofIdx
          END IF
        END DO
        CALL DistributedVector_DataRestore(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES,conditionTypes,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE")
    RETURN
!!TODO \todo write BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_FINALISE
999 ERRORSEXITS("BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise a boundary conditions create values cache 
  SUBROUTINE BoundaryConditions_CreateValuesCacheFinalise(createValuesCache,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the boundary conditions create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: constraintDofIdx

    ENTERS("BoundaryConditions_CreateValuesCacheFinalise",err,error,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ALLOCATED(createValuesCache%dofConstraints)) THEN
        DO constraintDofIdx=1,SIZE(createValuesCache%dofConstraints)
          CALL BoundaryConditions_DofConstraintFinalise(createValuesCache%dofConstraints(constraintDofIdx)%ptr,err,error,*999)
        END DO !constraintDofIdx
        DEALLOCATE(createValuesCache%dofConstraints)
      END IF
      IF(ALLOCATED(createValuesCache%dofCounts)) DEALLOCATE(createValuesCache%dofCounts)
      DEALLOCATE(createValuesCache)
    ENDIF
       
    EXITS("BoundaryConditions_CreateValuesCacheFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CreateValuesCacheFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_CreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a boundary conditions create values cache 
  SUBROUTINE BoundaryConditions_CreateValuesCacheInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_CreateValuesCacheInitialise",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable%createValuesCache)) THEN
      CALL FlagError("Boundary conditions create values cache is already associated.",err,error,*998)
    ELSE
      ALLOCATE(boundaryConditionsVariable%createValuesCache,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate boundary conditions create values cache.",err,error,*999)
      boundaryConditionsVariable%createValuesCache%numberOfConstraints=0
      ALLOCATE(boundaryConditionsVariable%createValuesCache%dofCounts(MAX_BOUNDARY_CONDITION_NUMBER),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate boundary conditions create values cache dof counts.",err,error,*999)
      boundaryConditionsVariable%createValuesCache%dofCounts=0
    ENDIF
       
    EXITS("BoundaryConditions_CreateValuesCacheInitialise")
    RETURN
999 CALL BoundaryConditions_CreateValuesCacheFinalise(boundaryConditionsVariable%createValuesCache,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditions_CreateValuesCacheInitialise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !

END MODULE BOUNDARY_CONDITIONS_ROUTINES
