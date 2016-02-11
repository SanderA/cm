!> \file
!> \author Chris Bradley
!> \brief This module handles all solver mapping routines.
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

!> This module handles all solver mapping routines.
MODULE SOLVER_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE COMP_ENVIRONMENT
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SOLVER_MAPPING_EquationsTypes SOLVER_MAPPING::EquationsTypes
  !> \brief Equations Matrix types
  !> \see SOLVER_MAPPING
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET=1 !<The equations in the solver mapping is from an equations set \see SOLVER_MAPPING_EquationsTypes,SOLVER_MAPPING
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION=2 !<The equations in the solver mapping is from an interface condition \see SOLVER_MAPPING_EquationsTypes,SOLVER_MAPPING
  !>@}
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION

  PUBLIC SOLVER_MAPPING_CREATE_FINISH,SOLVER_MAPPING_CREATE_START

  PUBLIC SOLVER_MAPPING_DESTROY
  
  PUBLIC SOLVER_MAPPING_EQUATIONS_SET_ADD

  PUBLIC SOLVER_MAPPING_INTERFACE_CONDITION_ADD

  PUBLIC SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
CONTAINS

  !
  !=================================================================================================================================
  !
  
  !>Calculates the solver mappings
  SUBROUTINE SOLVER_MAPPING_CALCULATE(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,equationsMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx,localDof,variableIdx, &
      & numberOfLocalSolverDofs,totalNumberOfLocalSolverDofs,constraintIdx,componentIdx,adjacentDomainIdx, &
      & solverMatrixIdx,solverVariableIdx,solverVariableIdx2,ghostSendIdx,ghostReceiveIdx
    INTEGER(INTG), POINTER :: dofTypes(:)
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING,ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(LIST_TYPE), POINTER :: ghostSendList,ghostReceiveList
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("SOLVER_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(SOLVER_MAPPING)) &
      & CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
    SOLVER_EQUATIONS=>SOLVER_MAPPING%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(SOLVER_EQUATIONS)) & 
      & CALL FlagError("The solver mapping solver equations are not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) &
      & CALL FlagError("The solver mapping create values cache is not associated.",ERR,ERROR,*999)
    BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
    IF(.NOT.ASSOCIATED(BOUNDARY_CONDITIONS)) &
      & CALL FlagError("The solver equations boundary conditions are not associated.",ERR,ERROR,*999)

    SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS=SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS
    SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS=SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_CONDITIONS
    SOLVER_MAPPING%NUMBER_OF_VARIABLES=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS+SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS

    ALLOCATE(SOLVER_MAPPING%VARIABLES(SOLVER_MAPPING%NUMBER_OF_VARIABLES),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate solver mapping variables.",ERR,ERROR,*999)

    !Allocate and initialise all solver variables
    variableIdx=0
    DO equationsSetIdx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
      variableIdx=variableIdx+1
      CALL SOLVER_MAPPING_VARIABLE_INITIALISE(SOLVER_MAPPING%VARIABLES(variableIdx),ERR,ERROR,*999)
      EQUATIONS_SET=>SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_SETS(equationsSetIdx)%PTR
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
            VARIABLE=>EQUATIONS_MAPPING%DEPENDENT_VARIABLE
            IF(ASSOCIATED(VARIABLE)) THEN
              NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE=>VARIABLE
                SOLVER_MAPPING%VARIABLES(variableIdx)%boundaryConditionsVariable=>BOUNDARY_CONDITIONS_VARIABLE
                SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE_TYPE=VARIABLE%VARIABLE_TYPE
                SOLVER_MAPPING%VARIABLES(variableIdx)%EQUATION_TYPE=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                SOLVER_MAPPING%VARIABLES(variableIdx)%EQUATION_INDEX=equationsSetIdx
                ALLOCATE(SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(VARIABLE% &
                  & TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap=0
                IF(ERR/=0) CALL FlagError("Could not allocate variable dof to solver dof map.",ERR,ERROR,*999)
              ELSE
                CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Dependent field variable is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set equations equations mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDDO !equationsSetIdx
    DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
      variableIdx=variableIdx+1
      CALL SOLVER_MAPPING_VARIABLE_INITIALISE(SOLVER_MAPPING%VARIABLES(variableIdx),ERR,ERROR,*999)
      INTERFACE_CONDITION=>SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
      IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
        SELECT CASE(INTERFACE_CONDITION%METHOD)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
            IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
              VARIABLE=>INTERFACE_MAPPING%LAGRANGE_VARIABLE
              IF(ASSOCIATED(VARIABLE)) THEN
                NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                  SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE=>VARIABLE
                  SOLVER_MAPPING%VARIABLES(variableIdx)%boundaryConditionsVariable=>BOUNDARY_CONDITIONS_VARIABLE
                  SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE_TYPE=VARIABLE%VARIABLE_TYPE
                  SOLVER_MAPPING%VARIABLES(variableIdx)%EQUATION_TYPE=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                  SOLVER_MAPPING%VARIABLES(variableIdx)%EQUATION_INDEX=interfaceConditionIdx
                  ALLOCATE(SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap( &
                    & VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                  SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap=0
                  IF(ERR/=0) CALL FlagError("Could not allocate variable dof to solver dof map.",ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Lagrange field variable is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interace equations equations mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface condition interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The interface condition method of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDDO !interfaceConditionIdx

    !Loop over all the local (excluding ghost) dofs of all variables
    numberOfLocalSolverDofs=0
    totalNumberOfLocalSolverDofs=0
    DO variableIdx=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
      VARIABLE=>SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE
      BOUNDARY_CONDITIONS_VARIABLE=>SOLVER_MAPPING%VARIABLES(variableIdx)%boundaryConditionsVariable
      NULLIFY(dofTypes)
      CALL DistributedVector_DataGet(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,dofTypes,err,error,*999)
      !Loop over the local (excluding ghost) dofs for this variable.
      constraintIdx=0
      DO localDof=1,VARIABLE%NUMBER_OF_DOFS
        IF(dofTypes(localDof)==BOUNDARY_CONDITION_DOF_FREE) THEN
          numberOfLocalSolverDofs=numberOfLocalSolverDofs+1
          totalNumberOfLocalSolverDofs=totalNumberOfLocalSolverDofs+1
          !Set the mappings
          SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof)=totalNumberOfLocalSolverDofs
        ELSE
          !The dof is constrained, so we map it to the negative constraintIdx. This way we know that the variable dof is
          !constrained and we immediately know which constraint it is.
          constraintIdx=constraintIdx+1
          SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof)=-constraintIdx
        END IF
      ENDDO !localDof
      CALL DistributedVector_DataRestore(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,dofTypes,err,error,*999)
    ENDDO !variableIdx
    !Loop over all the ghost dofs of all variables
    DO variableIdx=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
      VARIABLE=>SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE
      BOUNDARY_CONDITIONS_VARIABLE=>SOLVER_MAPPING%VARIABLES(variableIdx)%boundaryConditionsVariable
      NULLIFY(dofTypes)
      CALL DistributedVector_DataGet(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,dofTypes,err,error,*999)
      constraintIdx=BOUNDARY_CONDITIONS_VARIABLE%constraints%constrainedDofsMapping%NUMBER_OF_LOCAL
      !Loop over the ghost dofs for this variable.
      DO localDof=VARIABLE%NUMBER_OF_DOFS+1,VARIABLE%TOTAL_NUMBER_OF_DOFS
        IF(dofTypes(localDof)==BOUNDARY_CONDITION_DOF_FREE) THEN
          totalNumberOfLocalSolverDofs=totalNumberOfLocalSolverDofs+1
          !Set the mappings
          SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof)=totalNumberOfLocalSolverDofs
        ELSE
          !The dof is constrained, so we map it to the negative constraintIdx. This way we know that the variable dof is
          !constrained and we immediately know which constraint it is.
          constraintIdx=constraintIdx+1
          SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof)=-constraintIdx
        END IF
      ENDDO !localDof
      CALL DistributedVector_DataRestore(BOUNDARY_CONDITIONS_VARIABLE%DOF_TYPES,dofTypes,err,error,*999)
    ENDDO !variableIdx

    !Setup the adjacent domains information

    !How are interface variables and variables that belong to other regions handled?

    !Get the adjacent domain numbers. Do this by searching for the first valid domain type that we can find from the solver
    !variables.
    solverVariablesLoop: DO variableIdx=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
      VARIABLE=>SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE
      DO componentIdx=1,VARIABLE%NUMBER_OF_COMPONENTS
        DOMAIN=>VARIABLE%COMPONENTS(componentIdx)%DOMAIN
        !Are the domains of other regions also associated on this computational node?
        IF(ASSOCIATED(DOMAIN)) EXIT solverVariablesLoop
      END DO !componentIdx
    END DO solverVariablesLoop
    !Get the information about adjacent domains from the elements mapping
    ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
    !There has to be at least one valid elements mapping.
    IF(.NOT.ASSOCIATED(ELEMENTS_MAPPING)) THEN
      CALL FlagError("Elements mapping is not associated.",ERR,ERROR,*999)
    END IF

    ALLOCATE(SOLVER_MAPPING%DOFS_MAPPING,STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate solver mapping dofs mapping.",ERR,ERROR,*999)
    !!TODO: what is the real number of domains for a solver???
    CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLVER_MAPPING%DOFS_MAPPING, &
      & COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)

    !Allocate the solver dofs domain mapping
    DOFS_MAPPING=>SOLVER_MAPPING%DOFS_MAPPING
    DOFS_MAPPING%NUMBER_OF_LOCAL=numberOfLocalSolverDofs
    DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL=totalNumberOfLocalSolverDofs

    DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
    ALLOCATE(DOFS_MAPPING%ADJACENT_DOMAINS(DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",ERR,ERROR,*999)

    DO adjacentDomainIdx=1,DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
      CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(DOFS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx),ERR,ERROR,*999)
      DOFS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER=ELEMENTS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
        & DOMAIN_NUMBER
      NULLIFY(ghostSendList)
      CALL List_CreateStart(ghostSendList,ERR,ERROR,*999)
      CALL List_DataTypeSet(ghostSendList,LIST_INTG_TYPE,ERR,ERROR,*999)
      CALL List_InitialSizeSet(ghostSendList, &
        & MAX(totalNumberOfLocalSolverDofs-numberOfLocalSolverDofs,1),ERR,ERROR,*999)
      CALL List_CreateFinish(ghostSendList,ERR,ERROR,*999)
      NULLIFY(ghostReceiveList)
      CALL List_CreateStart(ghostReceiveList,ERR,ERROR,*999)
      CALL List_DataTypeSet(ghostReceiveList,LIST_INTG_TYPE,ERR,ERROR,*999)
      CALL List_InitialSizeSet(ghostReceiveList, &
        & MAX(totalNumberOfLocalSolverDofs-numberOfLocalSolverDofs,1),ERR,ERROR,*999)
      CALL List_CreateFinish(ghostReceiveList,ERR,ERROR,*999)
      DO variableIdx=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
        VARIABLE=>SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE
        DO ghostSendIdx=1,VARIABLE%DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS
          localDof=VARIABLE%DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES(ghostSendIdx)
          IF(SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof)>0) THEN
            CALL LIST_ITEM_ADD(ghostSendList,SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof), &
              & ERR,ERROR,*999)
          END IF
        END DO !ghostSendIdx
        DO ghostReceiveIdx=1,VARIABLE%DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
          localDof=VARIABLE%DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(ghostReceiveIdx)
          IF(SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof)>0) THEN
            CALL LIST_ITEM_ADD(ghostReceiveList,SOLVER_MAPPING%VARIABLES(variableIdx)%variableDofToSolverDofMap(localDof), &
              & ERR,ERROR,*999)
          END IF
        END DO !ghostReceiveIdx
      END DO !variableIdx
      CALL List_RemoveDuplicates(ghostSendList,ERR,ERROR,*999)
      CALL List_DetachAndDestroy(ghostSendList, &
        & DOFS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS, &
        & DOFS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES,ERR,ERROR,*999)
      CALL List_RemoveDuplicates(ghostReceiveList,ERR,ERROR,*999)
      CALL List_DetachAndDestroy(ghostReceiveList, &
        & DOFS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS, &
        & DOFS_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES,ERR,ERROR,*999)
    ENDDO !adjacentDomainIdx

    !Calculate the local to global map.
    CALL DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE(DOFS_MAPPING,err,error,*999)

    SOLVER_MAPPING%NUMBER_OF_DOFS=DOFS_MAPPING%NUMBER_OF_LOCAL
    SOLVER_MAPPING%TOTAL_NUMBER_OF_DOFS=DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
    SOLVER_MAPPING%NUMBER_OF_GLOBAL_DOFS=DOFS_MAPPING%NUMBER_OF_GLOBAL

    solverVariableIdx=0 
    ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate solver mapping equations set to solver map.",ERR,ERROR,*999)      
    DO equationsSetIdx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
      solverVariableIdx=solverVariableIdx+1
      EQUATIONS_SET=>SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_SETS(equationsSetIdx)%PTR
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
            CALL SolverMapping_EquationsSetToSolverMapInitialise(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
              & equationsSetIdx),ERR,ERROR,*999)
            SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS_SET_INDEX=equationsSetIdx
            SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%SOLVER_MAPPING=>SOLVER_MAPPING
            SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS_SET=>EQUATIONS_SET
            SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS=>EQUATIONS
            SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%solverVariableMap=> &
              & SOLVER_MAPPING%VARIABLES(solverVariableIdx)
            ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
             & SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations to solver matrix maps.",ERR,ERROR,*999)
            IF(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES==1) THEN
              solverMatrixIdx=1
              SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
               & solverMatrixIdx)%SOLVER_MATRIX_NUMBER=solverMatrixIdx
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                  & solverMatrixIdx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=DYNAMIC_MAPPING% &
                  & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                  & solverMatrixIdx)%dynamicMatrixToSolverVariableMap(DYNAMIC_MAPPING% &
                  & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate dynamic equations dof to solver dof maps .",ERR,ERROR,*999)
                DO equationsMatrixIdx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  DO solverVariableIdx2=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
                    IF(ASSOCIATED(SOLVER_MAPPING%VARIABLES(solverVariableIdx2)%VARIABLE, &
                      & DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equationsMatrixIdx)%VARIABLE)) THEN
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                        & solverMatrixIdx)%dynamicMatrixToSolverVariableMap(equationsMatrixIdx)%ptr=> &
                        & SOLVER_MAPPING%VARIABLES(solverVariableIdx2)
                      EXIT
                    END IF
                  END DO !solverVariableIdx2
                END DO !equationsMatrixIdx
              END IF
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                  & solverMatrixIdx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=LINEAR_MAPPING% &
                  & NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                  & solverMatrixIdx)%linearMatrixToSolverVariableMap(LINEAR_MAPPING% &
                  & NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate linear equations dof to solver dof maps .",ERR,ERROR,*999)
                DO equationsMatrixIdx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                  DO solverVariableIdx2=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
                    IF(ASSOCIATED(SOLVER_MAPPING%VARIABLES(solverVariableIdx2)%VARIABLE, &
                      & LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equationsMatrixIdx)%VARIABLE)) THEN
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                        & solverMatrixIdx)%linearMatrixToSolverVariableMap(equationsMatrixIdx)%ptr=> &
                        & SOLVER_MAPPING%VARIABLES(solverVariableIdx2)
                      EXIT
                    END IF
                  END DO !solverVariableIdx2
                END DO !equationsMatrixIdx
              END IF
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                  & solverMatrixIdx)%NUMBER_OF_EQUATIONS_JACOBIANS=NONLINEAR_MAPPING%NUMBER_OF_JACOBIANS
                ALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                  & solverMatrixIdx)%jacobianToSolverVariableMap(NONLINEAR_MAPPING%NUMBER_OF_JACOBIANS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate nonlinear equations dof to solver dof maps .",ERR,ERROR,*999)
                DO equationsMatrixIdx=1,NONLINEAR_MAPPING%NUMBER_OF_JACOBIANS
                  DO solverVariableIdx2=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
                    IF(ASSOCIATED(SOLVER_MAPPING%VARIABLES(solverVariableIdx2)%VARIABLE, &
                      & NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(equationsMatrixIdx)%VARIABLE)) THEN
                      SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps( &
                        & solverMatrixIdx)%jacobianToSolverVariableMap(equationsMatrixIdx)%ptr=> &
                        & SOLVER_MAPPING%VARIABLES(solverVariableIdx2)
                      EXIT
                    END IF
                  END DO !solverVariableIdx2
                END DO !equationsMatrixIdx
              END IF
            !More than 1 solver matrix is only appropriate for eigen analysis. This would than typically be a second block
            !diagonal solver matrix consisting of mass matrices corresponding to the solver variables. So perhaps put another
            !component "massMatrix" in the equations matrices type?
            ELSE IF(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES>1) THEN
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          ELSE
            CALL FlagError("Equations set equations mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDDO !equationsSetIdx
    
    ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate solver mapping interface condition to solver map.",ERR,ERROR,*999)
    DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
      solverVariableIdx=solverVariableIdx+1
      INTERFACE_CONDITION=>SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
      IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
        INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          CALL SolverMapping_InterfaceToSolverMapInitialise(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
            & interfaceConditionIdx),ERR,ERROR,*999)
          SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%INTERFACE_CONDITION_INDEX= &
            & interfaceConditionIdx
          SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%SOLVER_MAPPING=>SOLVER_MAPPING
          SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%INTERFACE_CONDITION=>INTERFACE_CONDITION
          SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
          SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%solverVariableMap=> &
            & SOLVER_MAPPING%VARIABLES(solverVariableIdx)
          ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
            & interfaceToSolverVariablesMaps(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate interface to solver matrix maps.",ERR,ERROR,*999)      
          IF(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES==1) THEN
            solverMatrixIdx=1
            SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%interfaceToSolverVariablesMaps( &
              & solverMatrixIdx)%SOLVER_MATRIX_NUMBER=solverMatrixIdx
            INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
            IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
              SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%interfaceToSolverVariablesMaps( &
                & solverMatrixIdx)%NUMBER_OF_INTERFACE_MATRICES=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
              ALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(equationsSetIdx)%interfaceToSolverVariablesMaps( &
                & solverMatrixIdx)%interfaceMatrixToSolverVariableMap(INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES), &
                & STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate interface dof to solver dof maps .",ERR,ERROR,*999)
              DO interfaceMatrixIdx=1,INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                DO solverVariableIdx2=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
                  IF(ASSOCIATED(SOLVER_MAPPING%VARIABLES(solverVariableIdx2)%VARIABLE, &
                    & INTERFACE_MAPPING%INTERFACE_MATRIX_TO_VAR_MAPS(interfaceMatrixIdx)%VARIABLE)) THEN
                    SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
                      & interfaceToSolverVariablesMaps(solverMatrixIdx)% &
                      & interfaceMatrixToSolverVariableMap(interfaceMatrixIdx)%ptr=> &
                      & SOLVER_MAPPING%VARIABLES(solverVariableIdx2)
                    EXIT
                  END IF
                END DO !solverVariableIdx2
              END DO !interfaceMatrixIdx
            END IF
          !More than 1 solver matrix is only appropriate for eigen analysis. This would than typically be a second block
          !diagonal solver matrix consisting of mass matrices corresponding to the solver variables. So perhaps put another
          !component "massMatrix" in the interface matrices type?
          ELSE IF(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES>1) THEN
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          END IF
        ELSE
          CALL FlagError("Interface condition interface equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDDO !interfaceConditionIdx

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Solver mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of solver matrices = ",SOLVER_MAPPING% &
        & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)               
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Solver variables :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",SOLVER_MAPPING%NUMBER_OF_VARIABLES, &
        & ERR,ERROR,*999)        
      DO variableIdx=1,SOLVER_MAPPING%NUMBER_OF_VARIABLES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Variable index: ",variableIdx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",SOLVER_MAPPING%VARIABLES(variableIdx)% &
          & VARIABLE_TYPE,ERR,ERROR,*999)        
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations type = ",SOLVER_MAPPING%VARIABLES(variableIdx)% &
          & EQUATION_TYPE,ERR,ERROR,*999)        
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Equations index = ",SOLVER_MAPPING%VARIABLES(variableIdx)% &
          & EQUATION_INDEX,ERR,ERROR,*999) 
        DO dofIdx=1,SOLVER_MAPPING%VARIABLES(variableIdx)%VARIABLE%TOTAL_NUMBER_OF_DOFS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable local dof : ",dofIdx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver local dof = ",SOLVER_MAPPING%VARIABLES(variableIdx)% &
            & variableDofToSolverDofMap(dofIdx),ERR,ERROR,*999)        
        ENDDO !dofIdx
      ENDDO !variableIdx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets to solver mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of equations sets = ",SOLVER_MAPPING% &
        & NUMBER_OF_EQUATIONS_SETS,ERR,ERROR,*999)
      DO equationsSetIdx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS_SET
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number        = ",EQUATIONS_SET%REGION%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",EQUATIONS_SET%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set to solver variable maps:",ERR,ERROR,*999)
        DO solverMatrixIdx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix : ",solverMatrixIdx,ERR,ERROR,*999)
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Dynamic mappings:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps(solverMatrixIdx)% &
              & NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES,ERR,ERROR,*999)
          END IF
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Linear mappings:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of linear equations matrices = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps(solverMatrixIdx)% &
              & NUMBER_OF_LINEAR_EQUATIONS_MATRICES,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Nonlinear mappings:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of equations jacobians = ",SOLVER_MAPPING% &
              & EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equationsToSolverVariablesMaps(solverMatrixIdx)% &
              & NUMBER_OF_EQUATIONS_JACOBIANS,ERR,ERROR,*999)
          ENDIF
        ENDDO !solverMatrixIdx
      END DO !equationsSetIdx
      IF(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS>0) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions to solver mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of interface conditions = ",SOLVER_MAPPING% &
          & NUMBER_OF_INTERFACE_CONDITIONS,ERR,ERROR,*999)
        DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
          INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
          & INTERFACE_CONDITION
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition index : ",interfaceConditionIdx, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Parent region user number       = ", &
            & INTERFACE_CONDITION%INTERFACE%PARENT_REGION%USER_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition user number = ", &
            & INTERFACE_CONDITION%USER_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition to solver variable maps:",ERR,ERROR,*999)
          DO solverMatrixIdx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,ERR,ERROR,*999)        
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of interface matrices = ",SOLVER_MAPPING% &
              & INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%interfaceToSolverVariablesMaps(solverMatrixIdx)% &
              & NUMBER_OF_INTERFACE_MATRICES,ERR,ERROR,*999)
          ENDDO !solverMatrixIdx        
        ENDDO !interfaceConditionIdx
      ENDIF
    ENDIF
    
    EXITS("SOLVER_MAPPING_CALCULATE")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver mapping
  SUBROUTINE SOLVER_MAPPING_CREATE_FINISH(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("SOLVER_MAPPING_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FlagError("Solver mapping has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
          CALL SOLVER_MAPPING_CALCULATE(SOLVER_MAPPING,ERR,ERROR,*999)
          CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLVER_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
          SOLVER_MAPPING%SOLVER_MAPPING_FINISHED=.TRUE.
        ELSE
          CALL FlagError("Solver mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    EXITS("SOLVER_MAPPING_CREATE_FINISH")
    RETURN
999 CALL SOLVER_MAPPING_FINALISE(SOLVER_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SOLVER_MAPPING_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%EQUATIONS_SETS)) DEALLOCATE(CREATE_VALUES_CACHE%EQUATIONS_SETS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%INTERFACE_CONDITIONS)) DEALLOCATE(CREATE_VALUES_CACHE%INTERFACE_CONDITIONS)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    EXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solver mapping create values cache 
  SUBROUTINE SolverMapping_CreateValuesCacheInitialise(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("SolverMapping_CreateValuesCacheInitialise",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FlagError("Solver mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate solver mapping create values cache.",ERR,ERROR,*999)
        SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS=0
        SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_CONDITIONS=0
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("SolverMapping_CreateValuesCacheInitialise")
    RETURN
999 CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLVER_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SolverMapping_CreateValuesCacheInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SolverMapping_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !
  !>Finishes the process of creating a solver mapping for a problem solver
  SUBROUTINE SOLVER_MAPPING_CREATE_START(SOLVER_EQUATIONS,SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to create the solver mapping on.
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<On return, a pointer to the solver mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("SOLVER_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FlagError("Solver equations has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          CALL FlagError("Solver mapping is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(SOLVER_MAPPING)
          CALL SOLVER_MAPPING_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MAPPING%CREATE_VALUES_CACHE)
          CALL SolverMapping_CreateValuesCacheInitialise(SOLVER_EQUATIONS%SOLVER_MAPPING,ERR,ERROR,*999)
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("SOLVER_MAPPING_CREATE_START")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_CREATE_START",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy a solver mapping.
  SUBROUTINE SOLVER_MAPPING_DESTROY(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SOLVER_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      CALL SOLVER_MAPPING_FINALISE(SOLVER_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Solver mapping is not associated",ERR,ERROR,*999)
    ENDIF
        
    EXITS("SOLVER_MAPPING_DESTROY")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_DESTROY",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a solver mapping
  SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_ADD(SOLVER_MAPPING,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On exit, the index of the equations set in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    LOGICAL :: addEquationsSet 
    TYPE(EQUATIONS_SET_PTR_TYPE), ALLOCATABLE :: NEW_EQUATIONS_SETS(:)
    
    ENTERS("SOLVER_MAPPING_EQUATIONS_SET_ADD",ERR,ERROR,*999)

    EQUATIONS_SET_INDEX=0
    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FlagError("Solver mapping has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_SET)) THEN
          IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
            IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
              IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_SETS)) THEN
                DO equationsSetIdx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS
                  IF(ASSOCIATED(EQUATIONS_SET,SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_SETS(equationsSetIdx)%PTR)) THEN
                    addEquationsSet=.FALSE.
                    EXIT
                  END IF
                END DO !equationsSetIdx
              END IF
              IF(addEquationsSet) THEN
                ALLOCATE(NEW_EQUATIONS_SETS(SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new equations sets.",ERR,ERROR,*999)
                DO equationsSetIdx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS
                  NEW_EQUATIONS_SETS(equationsSetIdx)%PTR=>SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                    & EQUATIONS_SETS(equationsSetIdx)%PTR
                END DO !equationsSetIdx
                CALL MOVE_ALLOC(NEW_EQUATIONS_SETS,SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_SETS)
                SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS= &
                  & SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS+1
                EQUATIONS_SET_INDEX=SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS
              END IF
            ELSE
              CALL FlagError("Solvers mapping create values cache is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    EXITS("SOLVER_MAPPING_EQUATIONS_SET_ADD")
    RETURN
999 IF(ALLOCATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
    ERRORSEXITS("SOLVER_MAPPING_EQUATIONS_SET_ADD",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_EQUATIONS_SET_ADD

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver map and deallocates all memory.
  SUBROUTINE SolverMapping_EquationsSetToSolverMapFinalise(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solverMatrixIdx
    
    ENTERS("SolverMapping_EquationsSetToSolverMapFinalise",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_SET_TO_SOLVER_MAP%equationsToSolverVariablesMaps)) THEN
      DO solverMatrixIdx=1,SIZE(EQUATIONS_SET_TO_SOLVER_MAP%equationsToSolverVariablesMaps,1)
        CALL SolverMapping_EquatsToSolVarMapsFinalise(EQUATIONS_SET_TO_SOLVER_MAP% &
          & equationsToSolverVariablesMaps(solverMatrixIdx),ERR,ERROR,*999)
      ENDDO !solverMatrixIdx
      DEALLOCATE(EQUATIONS_SET_TO_SOLVER_MAP%equationsToSolverVariablesMaps)
    ENDIF
        
    EXITS("SolverMapping_EquationsSetToSolverMapFinalise")
    RETURN
999 ERRORS("SolverMapping_EquationsSetToSolverMapFinalise",ERR,ERROR)    
    EXITS("SolverMapping_EquationsSetToSolverMapFinalise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetToSolverMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver map.
  SUBROUTINE SolverMapping_EquationsSetToSolverMapInitialise(EQUATIONS_SET_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TO_SOLVER_MAP_TYPE) :: EQUATIONS_SET_TO_SOLVER_MAP !<The equations set to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SolverMapping_EquationsSetToSolverMapInitialise",ERR,ERROR,*999)

    EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_SET_INDEX=0
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%EQUATIONS_SET)
    NULLIFY(EQUATIONS_SET_TO_SOLVER_MAP%solverVariableMap)

    EXITS("SolverMapping_EquationsSetToSolverMapInitialise")
    RETURN
999 ERRORS("SolverMapping_EquationsSetToSolverMapInitialise",ERR,ERROR)    
    EXITS("SolverMapping_EquationsSetToSolverMapInitialise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetToSolverMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrix map sm and deallocates all memory.
  SUBROUTINE SolverMapping_EquatsToSolVarMapsFinalise(EQUATIONS_TO_SOLVER_VARIABLE_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_VARIABLES_MAPS_TYPE) :: EQUATIONS_TO_SOLVER_VARIABLE_MAPS !<The equations set to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("SolverMapping_EquatsToSolVarMapsFinalise",ERR,ERROR,*999)

    IF(ALLOCATED(EQUATIONS_TO_SOLVER_VARIABLE_MAPS%dynamicMatrixToSolverVariableMap)) THEN
      DEALLOCATE(EQUATIONS_TO_SOLVER_VARIABLE_MAPS%dynamicMatrixToSolverVariableMap)
    ENDIF
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_VARIABLE_MAPS%linearMatrixToSolverVariableMap)) THEN
      DEALLOCATE(EQUATIONS_TO_SOLVER_VARIABLE_MAPS%linearMatrixToSolverVariableMap)
    ENDIF
    IF(ALLOCATED(EQUATIONS_TO_SOLVER_VARIABLE_MAPS%jacobianToSolverVariableMap)) THEN
      DEALLOCATE(EQUATIONS_TO_SOLVER_VARIABLE_MAPS%jacobianToSolverVariableMap)
    ENDIF

    EXITS("SolverMapping_EquatsToSolVarMapsFinalise")
    RETURN
999 ERRORSEXITS("SolverMapping_EquatsToSolVarMapsFinalise",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquatsToSolVarMapsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an equations to solver matrix maps sm.
  SUBROUTINE SolverMapping_EquatsToSolVarMapsInitialise(EQUATIONS_TO_SOLVER_VARIABLE_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TO_SOLVER_VARIABLES_MAPS_TYPE) :: EQUATIONS_TO_SOLVER_VARIABLE_MAPS !<The equations to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SolverMapping_EquatsToSolVarMapsInitialise",ERR,ERROR,*999)

    EQUATIONS_TO_SOLVER_VARIABLE_MAPS%SOLVER_MATRIX_NUMBER=0
    EQUATIONS_TO_SOLVER_VARIABLE_MAPS%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES=0
    EQUATIONS_TO_SOLVER_VARIABLE_MAPS%NUMBER_OF_LINEAR_EQUATIONS_MATRICES=0
    EQUATIONS_TO_SOLVER_VARIABLE_MAPS%NUMBER_OF_EQUATIONS_JACOBIANS=0

    EXITS("SolverMapping_EquatsToSolVarMapsInitialise")
    RETURN
999 ERRORS("SolverMapping_EquatsToSolVarMapsInitialise",ERR,ERROR)    
    EXITS("SolverMapping_EquatsToSolVarMapsInitialise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquatsToSolVarMapsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_FINALISE(SOLVER_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,variableIdx

    ENTERS("SOLVER_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(ALLOCATED(SOLVER_MAPPING%VARIABLES)) THEN
        DO variableIdx=1,SIZE(SOLVER_MAPPING%VARIABLES,1)
          CALL SOLVER_MAPPING_VARIABLE_FINALISE(SOLVER_MAPPING%VARIABLES(variableIdx),ERR,ERROR,*999)
        ENDDO !variableIdx
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)) THEN
        DO equationsSetIdx=1,SIZE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP,1)
          CALL SolverMapping_EquationsSetToSolverMapFinalise(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
            & equationsSetIdx),ERR,ERROR,*999)
        ENDDO !equationsSetIdx
        DEALLOCATE(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP)
      ENDIF
      IF(ALLOCATED(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP)) THEN
        DO interfaceConditionIdx=1,SIZE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP,1)
          CALL SolverMapping_InterfConditionToSolverMapFinalise(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
            & interfaceConditionIdx),ERR,ERROR,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP)
      ENDIF
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(SOLVER_MAPPING%DOFS_MAPPING,ERR,ERROR,*999)
      CALL SOLVER_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLVER_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(SOLVER_MAPPING)
    ENDIF
       
    EXITS("SOLVER_MAPPING_FINALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to initialise the solver mapping on.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("SOLVER_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) THEN
        CALL FlagError("Solver equations solver mapping is already associated",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER_EQUATIONS%SOLVER_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate solver equations solver mapping",ERR,ERROR,*999)
        SOLVER_EQUATIONS%SOLVER_MAPPING%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
        SOLVER_EQUATIONS%SOLVER_MAPPING%SOLVER_MAPPING_FINISHED=.FALSE.
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES=1
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_DOFS=0
        SOLVER_EQUATIONS%SOLVER_MAPPING%TOTAL_NUMBER_OF_DOFS=0
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_GLOBAL_DOFS=0
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS=0
        SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS=0
        NULLIFY(SOLVER_EQUATIONS%SOLVER_MAPPING%DOFS_MAPPING)
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF
    
    EXITS("SOLVER_MAPPING_INITIALISE")
    RETURN
999 CALL SOLVER_MAPPING_FINALISE(SOLVER_EQUATIONS%SOLVER_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SOLVER_MAPPING_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds an interface condition to a solver mapping
  SUBROUTINE SOLVER_MAPPING_INTERFACE_CONDITION_ADD(SOLVER_MAPPING,INTERFACE_CONDITION,INTERFACE_CONDITION_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer the solver mapping to add the interface condition to
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to add
    INTEGER(INTG), INTENT(OUT) :: INTERFACE_CONDITION_INDEX !<On exit, the index of the interface condition in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx
    LOGICAL :: EQUATIONS_SET_FOUND,addInterfaceCondition
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_CONDITION_PTR_TYPE), ALLOCATABLE :: NEW_INTERFACE_CONDITIONS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD",ERR,ERROR,*999)

    INTERFACE_CONDITION_INDEX=0
    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FlagError("Solver mapping has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
            IF(ASSOCIATED(SOLVER_MAPPING%CREATE_VALUES_CACHE)) THEN
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                  !Check that the interface variables are already part of an added equations set.
                  SELECT CASE(INTERFACE_CONDITION%METHOD)
                  CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                    interfaceMatrixIdx=1
                    EQUATIONS_SET=>INTERFACE_MAPPING%INTERFACE_MATRIX_TO_VAR_MAPS(interfaceMatrixIdx)%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      EQUATIONS_SET_FOUND=.FALSE.
                      DO equationsSetIdx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS
                        IF(ASSOCIATED(EQUATIONS_SET,SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                          & EQUATIONS_SETS(equationsSetIdx)%PTR)) THEN
                          EQUATIONS_SET_FOUND=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO !equationsSetIdx
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for interface matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(interfaceMatrixIdx,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The interface condition method of "// &
                      & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  IF(EQUATIONS_SET_FOUND) THEN
                    addInterfaceCondition=.TRUE.
                    IF(ALLOCATED(SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_CONDITIONS)) THEN
                      DO interfaceConditionIdx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_CONDITIONS
                        IF(ASSOCIATED(INTERFACE_CONDITION,SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                          & INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR)) THEN
                          addInterfaceCondition=.FALSE.
                          EXIT
                        END IF
                      END DO !interfaceConditionIdx
                    END IF
                    IF(addInterfaceCondition) THEN
                      ALLOCATE(NEW_INTERFACE_CONDITIONS(SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_CONDITIONS+1), &
                        & STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate new interface conditions.",ERR,ERROR,*999)
                      DO interfaceConditionIdx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_CONDITIONS
                        NEW_INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR=> &
                          & SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
                      END DO !interfaceConditionIdx
                      CALL MOVE_ALLOC(NEW_INTERFACE_CONDITIONS,SOLVER_MAPPING%CREATE_VALUES_CACHE%INTERFACE_CONDITIONS)
                      SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS=SOLVER_MAPPING%CREATE_VALUES_CACHE% &
                        & NUMBER_OF_INTERFACE_CONDITIONS+1
                      INTERFACE_CONDITION_INDEX=SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_CONDITIONS
                    END IF
                  ELSE
                    LOCAL_ERROR="The equations set associated with interface "// &
                      & "matrix number "//TRIM(NUMBER_TO_VSTRING(interfaceMatrixIdx,"*",ERR,ERROR))// &
                      & " has not been added to the solver equations."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface equations mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition interface equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solvers mapping create values cache is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface condition has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      END IF
    ELSE
      CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    EXITS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD")
    RETURN
999 IF(ALLOCATED(NEW_INTERFACE_CONDITIONS)) DEALLOCATE(NEW_INTERFACE_CONDITIONS)
    ERRORSEXITS("SOLVER_MAPPING_INTERFACE_CONDITION_ADD",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MAPPING_INTERFACE_CONDITION_ADD

  !
  !================================================================================================================================
  !

  !>Finalises an interface condition to solver map and deallocates all memory.
  SUBROUTINE SolverMapping_InterfConditionToSolverMapFinalise(INTERFACE_CONDITION_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE) :: INTERFACE_CONDITION_TO_SOLVER_MAP !<The interface condition to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solverMatrixIdx
    
    ENTERS("SolverMapping_InterfConditionToSolverMapFinalise",ERR,ERROR,*999)

    INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%SOLVER_MAPPING)
    IF(ALLOCATED(INTERFACE_CONDITION_TO_SOLVER_MAP%interfaceToSolverVariablesMaps)) THEN
      DO solverMatrixIdx=1,SIZE(INTERFACE_CONDITION_TO_SOLVER_MAP%interfaceToSolverVariablesMaps,1)
        CALL SolverMapping_InterfToSolVarMapsFinalise(INTERFACE_CONDITION_TO_SOLVER_MAP% &
          & interfaceToSolverVariablesMaps(solverMatrixIdx),ERR,ERROR,*999)
      ENDDO !solverMatrixIdx
      DEALLOCATE(INTERFACE_CONDITION_TO_SOLVER_MAP%interfaceToSolverVariablesMaps)
    ENDIF
        
    EXITS("SolverMapping_InterfConditionToSolverMapFinalise")
    RETURN
999 ERRORS("SolverMapping_InterfConditionToSolverMapFinalise",ERR,ERROR)    
    EXITS("SolverMapping_InterfConditionToSolverMapFinalise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfConditionToSolverMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition to solver map.
  SUBROUTINE SolverMapping_InterfaceToSolverMapInitialise(INTERFACE_CONDITION_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TO_SOLVER_MAP_TYPE) :: INTERFACE_CONDITION_TO_SOLVER_MAP !<The interface condition to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SolverMapping_InterfaceToSolverMapInitialise",ERR,ERROR,*999)

    INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_CONDITION_INDEX=0
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%SOLVER_MAPPING)
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_CONDITION)
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%INTERFACE_EQUATIONS)
    NULLIFY(INTERFACE_CONDITION_TO_SOLVER_MAP%solverVariableMap)
    
    EXITS("SolverMapping_InterfaceToSolverMapInitialise")
    RETURN
999 ERRORS("SolverMapping_InterfaceToSolverMapInitialise",ERR,ERROR)    
    EXITS("SolverMapping_InterfaceToSolverMapInitialise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceToSolverMapInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface to solver matrix map sm and deallocates all memory.
  SUBROUTINE SolverMapping_InterfToSolVarMapsFinalise(interfaceToSolverVariablesMaps,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_VARIABLES_MAPS_TYPE) :: interfaceToSolverVariablesMaps !<The interface to solver matrix maps sm to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("SolverMapping_InterfToSolVarMapsFinalise",ERR,ERROR,*999)

    interfaceToSolverVariablesMaps%SOLVER_MATRIX_NUMBER=0
    interfaceToSolverVariablesMaps%NUMBER_OF_INTERFACE_MATRICES=0
    IF(ALLOCATED(interfaceToSolverVariablesMaps%interfaceMatrixToSolverVariableMap)) THEN
      DEALLOCATE(interfaceToSolverVariablesMaps%interfaceMatrixToSolverVariableMap)
    ENDIF
   
    EXITS("SolverMapping_InterfToSolVarMapsFinalise")
    RETURN
999 ERRORSEXITS("SolverMapping_InterfToSolVarMapsFinalise",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfToSolVarMapsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface to solver matrix maps sm.
  SUBROUTINE SolverMapping_InterfToSolVarMapsInitialise(interfaceToSolverVariablesMaps,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TO_SOLVER_VARIABLES_MAPS_TYPE) :: interfaceToSolverVariablesMaps !<The interface to solver matrix maps sm to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SolverMapping_InterfToSolVarMapsInitialise",ERR,ERROR,*999)

    interfaceToSolverVariablesMaps%SOLVER_MATRIX_NUMBER=0
    interfaceToSolverVariablesMaps%NUMBER_OF_INTERFACE_MATRICES=0
        
    EXITS("SolverMapping_InterfToSolVarMapsInitialise")
    RETURN
999 ERRORS("SolverMapping_InterfToSolVarMapsInitialise",ERR,ERROR)    
    EXITS("SolverMapping_InterfToSolVarMapsInitialise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfToSolVarMapsInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solver matrices for the solver mapping
  SUBROUTINE SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SOLVER_MATRICES !<The number of solver matrices for the solver.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
      IF(SOLVER_MAPPING%SOLVER_MAPPING_FINISHED) THEN
        CALL FlagError("Solver mappings has been finished",ERR,ERROR,*999)
      ELSE
        MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES=1
        DO equationsSetIdx=1,SOLVER_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_EQUATIONS_SETS
          EQUATIONS_SET=>SOLVER_MAPPING%CREATE_VALUES_CACHE%EQUATIONS_SETS(equationsSetIdx)%PTR
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES=MAX(LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES, &
            & MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES)
          ENDIF
        ENDDO !equationsSetIdx
        !Check number of matrices to set is valid
        IF(NUMBER_OF_SOLVER_MATRICES>0.AND.NUMBER_OF_SOLVER_MATRICES<=MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES) THEN
          SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES=NUMBER_OF_SOLVER_MATRICES
        ELSE
          LOCAL_ERROR="The specified number of solver matrices of "// &
            & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
            & " is invalid. The number must be >= 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(MAXIMUM_NUMBER_OF_EQUATIONS_MATRICES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variable and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLE_FINALISE(SOLVER_MAPPING_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE) :: SOLVER_MAPPING_VARIABLE !<The solver mapping variable to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SOLVER_MAPPING_VARIABLE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_MAPPING_VARIABLE%variableDofToSolverDofMap)) &
      & DEALLOCATE(SOLVER_MAPPING_VARIABLE%variableDofToSolverDofMap)
       
    EXITS("SOLVER_MAPPING_VARIABLE_FINALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_VARIABLE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_VARIABLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SOLVER_MAPPING_VARIABLE_INITIALISE(SOLVER_MAPPING_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE) :: SOLVER_MAPPING_VARIABLE !<The solver mapping variable to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("SOLVER_MAPPING_VARIABLE_INITIALISE",ERR,ERROR,*999)

    NULLIFY(SOLVER_MAPPING_VARIABLE%VARIABLE)
    NULLIFY(SOLVER_MAPPING_VARIABLE%boundaryConditionsVariable)
    SOLVER_MAPPING_VARIABLE%VARIABLE_TYPE=0
    SOLVER_MAPPING_VARIABLE%EQUATION_TYPE=0
    SOLVER_MAPPING_VARIABLE%EQUATION_INDEX=0
    SOLVER_MAPPING_VARIABLE%EQUATION_TYPE=0
    
    EXITS("SOLVER_MAPPING_VARIABLE_INITIALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MAPPING_VARIABLE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MAPPING_VARIABLE_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE SOLVER_MAPPING_ROUTINES
