!> \file
!> \author Chris Bradley
!> \brief This module contains all constraint matrices routines.
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
!> Contributor(s): Sander Arens
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

!>This module contains all constraint matrices routines.
MODULE CONSTRAINT_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE CONSTRAINT_CONDITIONS_CONSTANTS
  USE CONSTRAINT_MATRICES_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes CONSTRAINT_MATRICES_ROUTINES::ConstraintMatrixStructureTypes
  !> \brief Constraint matrices structure (sparsity) types
  !> \see CONSTRAINT_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes,CONSTRAINT_MATRICES_ROUTINES 
  !>@}

  !> \addtogroup CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes CONSTRAINT_MATRICES_ROUTINES::ConstraintMatricesSparsityTypes
  !> \brief Constraint matrices sparsity types
  !> \see CONSTRAINT_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_SPARSE_MATRICES=1 !<Use sparse constraint matrices \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_FULL_MATRICES=2 !<Use fully populated constraint matrices \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes,CONSTRAINT_MATRICES_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Constraints

  PUBLIC CONSTRAINT_MATRIX_NO_STRUCTURE,CONSTRAINT_MATRIX_FEM_STRUCTURE

  PUBLIC CONSTRAINT_MATRICES_SPARSE_MATRICES,CONSTRAINT_MATRICES_FULL_MATRICES

  PUBLIC CONSTRAINT_MATRICES_CREATE_FINISH,CONSTRAINT_MATRICES_CREATE_START

  PUBLIC CONSTRAINT_MATRICES_DESTROY

  PUBLIC CONSTRAINT_MATRICES_ELEMENT_ADD

  PUBLIC CONSTRAINT_MATRICES_ELEMENT_CALCULATE

  PUBLIC CONSTRAINT_MATRICES_ELEMENT_FINALISE,ConstraintMatrices_ElementInitialise

  PUBLIC CONSTRAINT_MATRICES_OUTPUT

  PUBLIC CONSTRAINT_MATRICES_STORAGE_TYPE_SET

  PUBLIC CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET

  PUBLIC CONSTRAINT_MATRICES_VALUES_INITIALISE
  
  PUBLIC ConstraintMatrix_TimeDependenceTypeSet,ConstraintMatrix_TimeDependenceTypeGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalise a constraint matrix and deallocate all memory
  SUBROUTINE CONSTRAINT_MATRIX_FINALISE(CONSTRAINT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX !<A pointer to the constraint matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
      IF(ASSOCIATED(CONSTRAINT_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE)) CALL DISTRIBUTED_MATRIX_DESTROY(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE, &
        & ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      DEALLOCATE(CONSTRAINT_MATRIX)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an constraint matrix.
  SUBROUTINE CONSTRAINT_MATRIX_INITIALISE(CONSTRAINT_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the constraint matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The matrix number in the constraint matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRIX_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES) THEN
        CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          IF(ASSOCIATED(CONSTRAINT_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
            LOCAL_ERROR="Constraint matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
              & " is already associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          ELSE
            ALLOCATE(CONSTRAINT_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix.",ERR,ERROR,*999)
            CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
            CONSTRAINT_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
            CONSTRAINT_MATRIX%CONSTRAINT_MATRICES=>CONSTRAINT_MATRICES
            CONSTRAINT_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
            CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
            CONSTRAINT_MATRIX%UPDATE_MATRIX=.TRUE.
            CONSTRAINT_MATRIX%FIRST_ASSEMBLY=.TRUE.
            CONSTRAINT_MATRIX%HAS_TRANSPOSE=CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%HAS_TRANSPOSE
            CONSTRAINT_MATRIX%NUMBER_OF_ROWS=CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_ROWS
            CONSTRAINT_MATRIX%TOTAL_NUMBER_OF_ROWS=CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)% &
              & TOTAL_NUMBER_OF_ROWS
            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%CONSTRAINT_MATRIX=>CONSTRAINT_MATRIX
            NULLIFY(CONSTRAINT_MATRIX%MATRIX)
            NULLIFY(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE)
            NULLIFY(CONSTRAINT_MATRIX%TEMP_VECTOR)
            NULLIFY(CONSTRAINT_MATRIX%TEMP_TRANSPOSE_VECTOR)
            CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified constraint matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRIX_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRIX_FINALISE(CONSTRAINT_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRIX_INITIALISE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds the element matrices into the constraint matrices.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_ELEMENT_ADD()")
#endif

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      !Add the element matrices
      DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
        CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
          IF(CONSTRAINT_MATRIX%UPDATE_MATRIX) THEN
            !Add the element matrix into the distributed constraint equations matrix
            CALL DISTRIBUTED_MATRIX_VALUES_ADD(CONSTRAINT_MATRIX%MATRIX,CONSTRAINT_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
              & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
              & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
              & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
              & ERR,ERROR,*999)
            !If the constraint matrix has a transpose add it
            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,CONSTRAINT_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),TRANSPOSE(CONSTRAINT_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS)), &
                & ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="Constraint matrix for constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        IF(RHS_VECTOR%UPDATE_VECTOR) THEN
          !Add the rhs element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(RHS_VECTOR%RHS_VECTOR,RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS(1: &
            & RHS_VECTOR%ELEMENT_VECTOR%NUMBER_OF_ROWS),RHS_VECTOR%ELEMENT_VECTOR%VECTOR(1:RHS_VECTOR% &
            & ELEMENT_VECTOR%NUMBER_OF_ROWS),ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_ELEMENT_ADD()")
#endif
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_ADD")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_ADD",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_ADD")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Calculate the positions of the element matrices in the constraint matrices. 
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_CALCULATE(constraintMatrices,constraintElementNumber,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: constraintMatrices !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: constraintElementNumber !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,rowsElementNumber,rowsMeshIdx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable,rowsFieldVariable
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: constraintCondition
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: constraintEquations
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: constraintMapping
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: constraintMatrix
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: rhsVector
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: meshConnectivity
    TYPE(ConstraintPointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(VARYING_STRING) :: localError

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE",err,error,*999)

    IF(ASSOCIATED(constraintMatrices)) THEN
      constraintMapping=>constraintMatrices%CONSTRAINT_MAPPING
      IF(ASSOCIATED(constraintMapping)) THEN
        constraintEquations=>constraintMapping%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(constraintEquations)) THEN
          constraintCondition=>constraintEquations%CONSTRAINT_CONDITION
          IF(ASSOCIATED(constraintCondition)) THEN
            constraint=>constraintCondition%CONSTRAINT
            IF(ASSOCIATED(constraint)) THEN
              SELECT CASE(constraintCondition%integrationType)
              CASE(CONSTRAINT_CONDITION_GAUSS_INTEGRATION)
                meshConnectivity=>constraint%MESH_CONNECTIVITY
                IF(ASSOCIATED(meshConnectivity)) THEN
                  IF(ALLOCATED(meshConnectivity%ELEMENT_CONNECTIVITY)) THEN
                    !Calculate the row and columns for the constraint equations matrices
                    DO matrixIdx=1,constraintMatrices%NUMBER_OF_CONSTRAINT_MATRICES
                      constraintMatrix=>constraintMatrices%MATRICES(matrixIdx)%PTR
                      IF(ASSOCIATED(constraintMatrix)) THEN
                        rowsFieldVariable=>constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                        colsFieldVariable=>constraintMapping%LAGRANGE_VARIABLE !\todo: TEMPORARY: Needs generalising
                        rowsMeshIdx=constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%MESH_INDEX
                        IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
                          ! If the rows and column variables are both the Lagrange variable (this is the diagonal matrix)
                          rowsElementNumber=ConstraintElementNumber
                        ELSE
                          rowsElementNumber=meshConnectivity%ELEMENT_CONNECTIVITY(ConstraintElementNumber,rowsMeshIdx)% &
                            & COUPLED_MESH_ELEMENT_NUMBER
                        ENDIF
                        CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(constraintMatrix%ELEMENT_MATRIX, &
                          & constraintMatrix%UPDATE_MATRIX,[rowsElementNumber],[constraintElementNumber],rowsFieldVariable, &
                          & colsFieldVariable,err,error,*999)
                      ELSE
                        localError="Constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                          & " is not associated."
                        CALL FLAG_ERROR(localError,err,error,*999)
                      ENDIF
                    ENDDO !matrixIdx
                  ELSE
                    CALL FLAG_ERROR("Constraint element connectivity is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Constraint mesh connectivity is not associated.",err,error,*999)              
                ENDIF
              CASE(CONSTRAINT_CONDITION_DATA_POINTS_INTEGRATION)
                pointsConnectivity=>constraint%pointsConnectivity
                IF(ASSOCIATED(pointsConnectivity)) THEN
                  IF(ALLOCATED(pointsConnectivity%coupledElements)) THEN
                    DO matrixIdx=1,constraintMatrices%NUMBER_OF_CONSTRAINT_MATRICES
                      constraintMatrix=>constraintMatrices%MATRICES(matrixIdx)%PTR
                      IF(ASSOCIATED(constraintMatrix)) THEN 
                        IF(constraintCondition%METHOD==CONSTRAINT_CONDITION_PENALTY_METHOD .AND. &
                            matrixIdx==constraintMatrices%NUMBER_OF_CONSTRAINT_MATRICES) THEN
                          rowsFieldVariable=>constraintMapping%LAGRANGE_VARIABLE
                          colsFieldVariable=>constraintMapping%LAGRANGE_VARIABLE
                          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(constraintMatrix%ELEMENT_MATRIX, &
                            & constraintMatrix%UPDATE_MATRIX,[ConstraintElementNumber],[ConstraintElementNumber], &
                            & rowsFieldVariable,colsFieldVariable,err,error,*999)
                        ELSE
                          rowsFieldVariable=>constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                          colsFieldVariable=>constraintMapping%LAGRANGE_VARIABLE !\todo: TEMPORARY: Needs generalising
                          rowsMeshIdx=constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%MESH_INDEX
                          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(constraintMatrix%ELEMENT_MATRIX, &
                            & constraintMatrix%UPDATE_MATRIX,pointsConnectivity%coupledElements(ConstraintElementNumber, &
                            & rowsMeshIdx)%elementNumbers,[ConstraintElementNumber],rowsFieldVariable,colsFieldVariable, &
                            & err,error,*999)
                        ENDIF
                      ELSE
                        localError="Constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                          & " is not associated."
                        CALL FLAG_ERROR(localError,err,error,*999)
                      ENDIF
                    ENDDO !matrixIdx
                  ELSE
                    CALL FLAG_ERROR("Constraint points connectivity coupled elements is not allocated.",err,error,*999)             
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)              
                ENDIF
              CASE DEFAULT
                localError="The constraint condition integration type of "// &
                  & TRIM(NUMBER_TO_VSTRING(constraintCondition%integrationType,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(localError,ERR,ERROR,*999)
              END SELECT
              !RHS element matrix dofs are the same for both mesh and points connectivity, right now
              rhsVector=>constraintMatrices%RHS_VECTOR
              IF(ASSOCIATED(rhsVector)) THEN
                rhsMapping=>constraintMapping%RHS_MAPPING
                IF(ASSOCIATED(rhsMapping)) THEN
                  !Calculate the rows  for the equations RHS
                  rowsFieldVariable=>rhsMapping%RHS_VARIABLE
                  CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(rhsVector%ELEMENT_VECTOR,rhsVector%UPDATE_VECTOR, &
                    & constraintElementNumber,rowsFieldVariable,err,error,*999)
                ELSE
                  CALL FLAG_ERROR("Constraint mapping rhs mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint condition constraint is not associated.",err,error,*999)            
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint mapping is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not allocated",err,error,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_ELEMENT_CALCULATE()")
#endif
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE",err,error)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information for constraint matrices and deallocate all memory
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<The constraint matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
        CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="Constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        !Finalise the constraint element vector
        RHS_VECTOR%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)
      ENDIF

      ENDDO !matrix_idx
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the constraint matrices
  SUBROUTINE ConstraintMatrices_ElementInitialise(constraintMatrices,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: constraintMatrices !The constraint matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,rowsMeshIdx
    INTEGER(INTG) :: rowsNumberOfElements,colsNumberOfElements !Number of elements in the row and col variables whose dofs are present in constraint element matrix
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: constraintMapping
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: constraintCondition
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: constraintEquations
    TYPE(ConstraintPointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: constraintMatrix
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: rhsVector
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable,rowsFieldVariable
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("ConstraintMatrices_ElementInitialise",err,error,*999)

    IF(ASSOCIATED(constraintMatrices)) THEN
      constraintMapping=>constraintMatrices%CONSTRAINT_MAPPING
      IF(ASSOCIATED(constraintMapping)) THEN
        constraintEquations=>constraintMapping%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(constraintEquations)) THEN
          constraintCondition=>constraintEquations%CONSTRAINT_CONDITION
          IF(ASSOCIATED(constraintCondition)) THEN
            SELECT CASE(constraintCondition%integrationType)
              CASE(CONSTRAINT_CONDITION_GAUSS_INTEGRATION)
              DO matrixIdx=1,constraintMatrices%NUMBER_OF_CONSTRAINT_MATRICES
                constraintMatrix=>constraintMatrices%MATRICES(matrixIdx)%PTR
                IF(ASSOCIATED(constraintMatrix)) THEN
                  rowsNumberOfElements=1
                  colsNumberOfElements=1
                  rowsFieldVariable=>constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                  colsFieldVariable=>constraintMapping%LAGRANGE_VARIABLE !TEMPORARY: Needs generalising
                  CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(constraintMatrix%ELEMENT_MATRIX,rowsFieldVariable, &
                    & colsFieldVariable,rowsNumberOfElements,colsNumberOfElements,err,error,*999)
                ELSE
                  localError="Constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                    & " is not associated."
                  CALL FLAG_ERROR(localError,err,error,*999)
                ENDIF
              ENDDO !matrixIdx
            CASE(CONSTRAINT_CONDITION_DATA_POINTS_INTEGRATION) 
              constraint=>constraintCondition%CONSTRAINT
              IF(ASSOCIATED(constraint))THEN
                pointsConnectivity=>constraint%pointsConnectivity
                IF(ASSOCIATED(pointsConnectivity)) THEN
                  IF(ALLOCATED(pointsConnectivity%coupledElements)) THEN
                    DO matrixIdx=1,constraintMatrices%NUMBER_OF_CONSTRAINT_MATRICES !\todo: Need to separate the case for penalty matrix
                      constraintMatrix=>constraintMatrices%MATRICES(matrixIdx)%PTR
                      IF(ASSOCIATED(constraintMatrix)) THEN
                        colsNumberOfElements=1
                        rowsFieldVariable=>constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                        colsFieldVariable=>constraintMapping%LAGRANGE_VARIABLE !TEMPORARY: Needs generalising
                        rowsMeshIdx=constraintMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%MESH_INDEX
                        CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(constraintMatrix%ELEMENT_MATRIX,rowsFieldVariable, &
                          & colsFieldVariable,pointsConnectivity%maxNumberOfCoupledElements(rowsMeshIdx), &
                          & colsNumberOfElements,err,error,*999)
                      ELSE
                        localError="Constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                          & " is not associated."
                        CALL FLAG_ERROR(localError,err,error,*999)
                      ENDIF
                    ENDDO !matrixIdx
                  ELSE
                    CALL FLAG_ERROR("Constraint points connectivity coupled elements is not allocated.",err,error,*999)             
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)              
                ENDIF
              ELSE
                CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The constraint condition integration type of "// &
                & TRIM(NUMBER_TO_VSTRING(constraintCondition%integrationType,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(localError,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Constraint condition is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations is not associated.",err,error,*999)
        ENDIF
        rhsVector=>constraintMatrices%RHS_VECTOR
        IF(ASSOCIATED(rhsVector)) THEN
          !Initialise the RHS element vector
          rhsMapping=>constraintMapping%RHS_MAPPING
          IF(ASSOCIATED(rhsMapping)) THEN
            rowsFieldVariable=>rhsMapping%RHS_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(rhsVector%ELEMENT_VECTOR,rowsFieldVariable,err,error,*999)
          ELSE
            CALL FLAG_ERROR("RHS mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices mapping is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintMatrices_ElementInitialise")
    RETURN
999 CALL ERRORS("ConstraintMatrices_ElementInitialise",err,error)
    CALL EXITS("ConstraintMatrices_ElementInitialise")
    RETURN 1
  END SUBROUTINE ConstraintMatrices_ElementInitialise

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for an constraint matrix.
  SUBROUTINE CONSTRAINT_MATRIX_STRUCTURE_CALCULATE(CONSTRAINT_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
    & TRANSPOSE_ROW_INDICES,TRANSPOSE_COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX !<A pointer to the constraint matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return, the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return, a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return, a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: TRANSPOSE_ROW_INDICES(:) !<On return, if the constraint matrix has a transpose a pointer to transpose row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: TRANSPOSE_COLUMN_INDICES(:) !<On return, if the constraint matrix has a transpose a pointer to the transpose column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   INTEGER(INTG) :: column_version,column_derivative,column_idx,column_component_idx,column_local_derivative_idx, &
     & column_local_node_idx, column_node,DUMMY_ERR,domain_element,global_column,global_row,constraint_element_idx, &
     & CONSTRAINT_MESH_INDEX,local_column,local_row,MATRIX_NUMBER,NUMBER_OF_COLUMNS,NUMBER_OF_ROWS,row_component_idx, &
     & row_version,row_derivative,row_local_derivative_idx,row_idx,row_local_node_idx,row_node,TRANSPOSE_NUMBER_OF_NON_ZEROS
    INTEGER(INTG), ALLOCATABLE :: COLUMNS(:),TRANSPOSE_COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: COLUMN_BASIS,ROW_BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOFS_DOMAIN_MAPPING,ROW_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: COLUMN_DOMAIN_ELEMENTS,ROW_DOMAIN_ELEMENTS
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: MESH_CONNECTIVITY
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: COLUMN_DOFS_PARAM_MAPPING,ROW_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLUMN_VARIABLE,ROW_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: TRANSPOSE_COLUMN_INDICES_LISTS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*999)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
      IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
        IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
          IF(.NOT.ASSOCIATED(TRANSPOSE_ROW_INDICES)) THEN
            IF(.NOT.ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) THEN
              MATRIX_NUMBER=CONSTRAINT_MATRIX%MATRIX_NUMBER
              SELECT CASE(CONSTRAINT_MATRIX%STRUCTURE_TYPE)
              CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
                CALL FLAG_ERROR("There is no structure to calculate for a matrix with no structure.",ERR,ERROR,*998)
              CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
                SELECT CASE(CONSTRAINT_MATRIX%STORAGE_TYPE)
                CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                  CONSTRAINT_MATRICES=>CONSTRAINT_MATRIX%CONSTRAINT_MATRICES
                  IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
                    CONSTRAINT_EQUATIONS=>CONSTRAINT_MATRICES%CONSTRAINT_EQUATIONS
                    IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
                      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
                        CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
                        IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
                          CONSTRAINT=>CONSTRAINT_CONDITION%CONSTRAINT
                          IF(ASSOCIATED(CONSTRAINT)) THEN
                            MESH_CONNECTIVITY=>CONSTRAINT%MESH_CONNECTIVITY
                            IF(ASSOCIATED(MESH_CONNECTIVITY)) THEN
                              ROW_VARIABLE=>CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%VARIABLE
                              CONSTRAINT_MESH_INDEX=CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%MESH_INDEX
                              IF(ASSOCIATED(ROW_VARIABLE)) THEN
                                COLUMN_VARIABLE=>CONSTRAINT_MAPPING%LAGRANGE_VARIABLE
                                IF(ASSOCIATED(COLUMN_VARIABLE)) THEN
                                  ROW_DOFS_DOMAIN_MAPPING=>ROW_VARIABLE%DOMAIN_MAPPING
                                  IF(ASSOCIATED(ROW_DOFS_DOMAIN_MAPPING)) THEN
                                    COLUMN_DOFS_DOMAIN_MAPPING=>COLUMN_VARIABLE%DOMAIN_MAPPING
                                    IF(ASSOCIATED(COLUMN_DOFS_DOMAIN_MAPPING)) THEN
                                      ROW_DOFS_PARAM_MAPPING=>ROW_VARIABLE%DOF_TO_PARAM_MAP
                                      IF(ASSOCIATED(ROW_DOFS_PARAM_MAPPING)) THEN
                                        COLUMN_DOFS_PARAM_MAPPING=>COLUMN_VARIABLE%DOF_TO_PARAM_MAP
                                        IF(ASSOCIATED(COLUMN_DOFS_PARAM_MAPPING)) THEN
                                          !Allocate lists
                                          ALLOCATE(COLUMN_INDICES_LISTS(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",ERR,ERROR,*999)
                                          DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            !Set up list
                                            NULLIFY(COLUMN_INDICES_LISTS(local_row)%PTR)
                                            CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                            CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_row)%PTR,LIST_INTG_TYPE, &
                                              & ERR,ERROR,*999)
                                            CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_row)%PTR,50,ERR,ERROR,*999)
                                            CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                          ENDDO !local_row
                                          !Allocate row indices
                                          ALLOCATE(ROW_INDICES(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",ERR,ERROR,*999)
                                          IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                            !Allocate transpose lists
                                            ALLOCATE(TRANSPOSE_COLUMN_INDICES_LISTS(COLUMN_DOFS_DOMAIN_MAPPING% &
                                              & TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate transpose column indices lists.", &
                                              & ERR,ERROR,*999)
                                            DO local_column=1,COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                              !Set up list
                                              NULLIFY(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR)
                                              CALL LIST_CREATE_START(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & ERR,ERROR,*999)
                                              CALL LIST_DATA_TYPE_SET(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & LIST_INTG_TYPE,ERR,ERROR,*999)
                                              CALL LIST_INITIAL_SIZE_SET(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR,50, &
                                                & ERR,ERROR,*999)
                                              CALL LIST_CREATE_FINISH(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & ERR,ERROR,*999)
                                            ENDDO !local_column
                                            !Allocate transpose row indices
                                            ALLOCATE(TRANSPOSE_ROW_INDICES(COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1), &
                                              & STAT=ERR)
                                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate transpose row indices.",ERR,ERROR,*999)
                                          ENDIF
                                          !Loop over the number of components in the Lagrange multipler variable
                                          DO column_component_idx=1,COLUMN_VARIABLE%NUMBER_OF_COMPONENTS
                                            IF(COLUMN_VARIABLE%COMPONENTS(column_component_idx)%INTERPOLATION_TYPE== &
                                              & FIELD_NODE_BASED_INTERPOLATION) THEN
                                              !Loop over the elements in the constraint mesh
                                              COLUMN_DOMAIN_ELEMENTS=>COLUMN_VARIABLE%COMPONENTS(column_component_idx)%DOMAIN% &
                                                & TOPOLOGY%ELEMENTS
                                              DO constraint_element_idx=1,COLUMN_DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                                                COLUMN_BASIS=>COLUMN_DOMAIN_ELEMENTS%ELEMENTS(constraint_element_idx)%BASIS
                                                !Loop over the column DOFs in the element
                                                DO column_local_node_idx=1,COLUMN_BASIS%NUMBER_OF_NODES
                                                  column_node=COLUMN_DOMAIN_ELEMENTS%ELEMENTS(constraint_element_idx)% &
                                                    & ELEMENT_NODES(column_local_node_idx)
                                                  DO column_local_derivative_idx=1,COLUMN_BASIS% &
                                                    & NUMBER_OF_DERIVATIVES(column_local_node_idx)
                                                    column_derivative=COLUMN_DOMAIN_ELEMENTS%ELEMENTS(constraint_element_idx)% &
                                                      & ELEMENT_DERIVATIVES(column_local_derivative_idx,column_local_node_idx)
                                                    column_version=COLUMN_DOMAIN_ELEMENTS%ELEMENTS(constraint_element_idx)% &
                                                      & elementVersions(column_local_derivative_idx,column_local_node_idx)
                                                    local_column=COLUMN_VARIABLE%COMPONENTS(column_component_idx)% &
                                                      & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(column_node)% &
                                                      & DERIVATIVES(column_derivative)%VERSIONS(column_version)
                                                    global_column=COLUMN_DOFS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                                    !Loop over the components in the dependent variable
                                                    DO row_component_idx=1,ROW_VARIABLE%NUMBER_OF_COMPONENTS
                                                      SELECT CASE(ROW_VARIABLE%COMPONENTS(row_component_idx)%INTERPOLATION_TYPE)
                                                      CASE(FIELD_CONSTANT_INTERPOLATION)
                                                        local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                          & CONSTANT_PARAM2DOF_MAP
                                                        CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                          & ERR,ERROR,*999)
                                                        IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                                          global_row=ROW_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                          CALL LIST_ITEM_ADD(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                            & global_row,ERR,ERROR,*999)
                                                        ENDIF
                                                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                                        domain_element=MESH_CONNECTIVITY% &
                                                          & ELEMENT_CONNECTIVITY(constraint_element_idx,CONSTRAINT_MESH_INDEX)% &
                                                          & COUPLED_MESH_ELEMENT_NUMBER
                                                        local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                          & ELEMENT_PARAM2DOF_MAP%ELEMENTS(domain_element)
                                                        CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                          & ERR,ERROR,*999)
                                                        IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                                          global_row=ROW_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                          CALL LIST_ITEM_ADD(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                            & global_row,ERR,ERROR,*999)
                                                        ENDIF
                                                      CASE(FIELD_NODE_BASED_INTERPOLATION)
                                                        ROW_DOMAIN_ELEMENTS=>ROW_VARIABLE%COMPONENTS(row_component_idx)%DOMAIN% &
                                                          & TOPOLOGY%ELEMENTS
                                                        domain_element=MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY( &
                                                          & constraint_element_idx,CONSTRAINT_MESH_INDEX)%COUPLED_MESH_ELEMENT_NUMBER
                                                        ROW_BASIS=>ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)%BASIS
                                                        !Loop over the row DOFs in the domain mesh element
                                                        DO row_local_node_idx=1,ROW_BASIS%NUMBER_OF_NODES
                                                          row_node=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                            & ELEMENT_NODES(row_local_node_idx)
                                                          DO row_local_derivative_idx=1,ROW_BASIS% &
                                                            & NUMBER_OF_DERIVATIVES(row_local_node_idx)
                                                            row_derivative=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                              & ELEMENT_DERIVATIVES(row_local_derivative_idx,row_local_node_idx)
                                                            row_version=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                              & elementVersions(row_local_derivative_idx,row_local_node_idx)
                                                            local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)% &
                                                              & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(row_node)% &
                                                              & DERIVATIVES(row_derivative)%VERSIONS(row_version)
                                                            CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                              & ERR,ERROR,*999)
                                                            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                                              global_row=ROW_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                              CALL LIST_ITEM_ADD(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                                & global_row,ERR,ERROR,*999)
                                                            ENDIF
                                                          ENDDO !row_local_derivative_idx
                                                        ENDDO !row_local_node_idx
                                                      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                                      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                                      CASE DEFAULT
                                                        LOCAL_ERROR="The row variable interpolation type of "// &
                                                          & TRIM(NUMBER_TO_VSTRING(ROW_VARIABLE%COMPONENTS(row_component_idx)% &
                                                          INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid."
                                                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                                      END SELECT
                                                    ENDDO !row_component_idx
                                                  ENDDO !column_local_derivative_idx
                                                ENDDO !column_local_node_idx
                                              ENDDO !constraint_element_idx
                                            ELSE
                                              CALL FLAG_ERROR("Only node based fields implemented.",ERR,ERROR,*999)
                                            ENDIF
                                          ENDDO !column_component_idx
                                          ROW_INDICES(1)=1
                                          DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                            CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                              & ERR,ERROR,*999)
                                            NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                            ROW_INDICES(local_row+1)=NUMBER_OF_NON_ZEROS+1
                                          ENDDO !local_row
                                          IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                            TRANSPOSE_NUMBER_OF_NON_ZEROS=0
                                            TRANSPOSE_ROW_INDICES(1)=1
                                            DO local_column=1,COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                              CALL LIST_REMOVE_DUPLICATES(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & ERR,ERROR,*999)
                                              CALL LIST_NUMBER_OF_ITEMS_GET(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                                              TRANSPOSE_NUMBER_OF_NON_ZEROS=TRANSPOSE_NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                              TRANSPOSE_ROW_INDICES(local_column+1)=TRANSPOSE_NUMBER_OF_NON_ZEROS+1
                                            ENDDO !local_column
                                            !Sanity check - the number of non-zeros should be the same
                                            IF(TRANSPOSE_NUMBER_OF_NON_ZEROS/=NUMBER_OF_NON_ZEROS) THEN
                                              LOCAL_ERROR="Invalid number of non-zeros. The number of non-zeros in the "// &
                                                & "transposed matrix ("//TRIM(NUMBER_TO_VSTRING(TRANSPOSE_NUMBER_OF_NON_ZEROS, &
                                                & "*",ERR,ERROR))//") does not match the number of non-zeros in the constraint "// &
                                                & "matrix ("//TRIM(NUMBER_TO_VSTRING(TRANSPOSE_NUMBER_OF_NON_ZEROS,"*",ERR, &
                                                & ERROR))//")."
                                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                            ENDIF
                                          ENDIF
                                          !Allocate and setup the column locations
                                          ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",ERR,ERROR,*999)
                                          DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                              & COLUMNS,ERR,ERROR,*999)
                                            DO column_idx=1,NUMBER_OF_COLUMNS
                                              COLUMN_INDICES(ROW_INDICES(local_row)+column_idx-1)=COLUMNS(column_idx)
                                            ENDDO !column_idx
                                            DEALLOCATE(COLUMNS)
                                          ENDDO !local_row
                                          IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                            !Allocate and setup the column locations
                                            ALLOCATE(TRANSPOSE_COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                            IF(ERR/=0) &
                                              & CALL FLAG_ERROR("Could not allocate transpose column indices.",ERR,ERROR,*999)
                                            DO local_column=1,COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                              CALL LIST_DETACH_AND_DESTROY(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & NUMBER_OF_ROWS,TRANSPOSE_COLUMNS,ERR,ERROR,*999)
                                              DO row_idx=1,NUMBER_OF_ROWS
                                                TRANSPOSE_COLUMN_INDICES(TRANSPOSE_ROW_INDICES(local_column)+row_idx-1)= &
                                                  & TRANSPOSE_COLUMNS(row_idx)
                                              ENDDO !row_idx
                                              DEALLOCATE(TRANSPOSE_COLUMNS)
                                            ENDDO !local_column
                                          ENDIF
                                          IF(DIAGNOSTICS1) THEN
                                            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Constraint matrix structure:",ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Constraint matrix number : ", &
                                              & MATRIX_NUMBER,ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                              & ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                              & COLUMN_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                              & NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                                            IF(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                                              & COLUMN_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL/=0) THEN
                                              SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(ROW_DOFS_DOMAIN_MAPPING% &
                                                & TOTAL_NUMBER_OF_LOCAL*COLUMN_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,DP)*100.0_DP
                                              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY, &
                                                & "F6.2",ERR,ERROR,*999)
                                            ENDIF
                                            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ROW_DOFS_DOMAIN_MAPPING% &
                                              & TOTAL_NUMBER_OF_LOCAL+1,5,5,ROW_INDICES, &
                                              & '("  Row indices              :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
                                            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                              & COLUMN_INDICES,'("  Column indices           :",5(X,I13))','(28X,5(X,I13))', &
                                              & ERR,ERROR,*999)
                                            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN 
                                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,COLUMN_DOFS_DOMAIN_MAPPING% &
                                                & TOTAL_NUMBER_OF_LOCAL+1,5,5,TRANSPOSE_ROW_INDICES, &
                                                & '("  Transpose row indices    :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
                                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                                & TRANSPOSE_COLUMN_INDICES,'("  Transpose column indices :",5(X,I13))', &
                                                & '(28X,5(X,I13))',ERR,ERROR,*999)
                                            ENDIF
                                          ENDIF
                                        ELSE
                                          CALL FLAG_ERROR("Column dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Row dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Column dofs domain mapping is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Row dofs domain mapping is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Column field variable is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Row field variable is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Constraint mesh connectivity is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Constraint condition constraint is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE DEFAULT
                  LOCAL_ERROR="The matrix storage type of "// &
                    & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE DEFAULT
                LOCAL_ERROR="The matrix structure type of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MATRIX%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Transpose column indices is already associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Transpose row indieces is already associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Column indices is already associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indieces is already associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*998)
    ENDIF
     
    CALL EXITS("CONSTRAINT_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(TRANSPOSE_ROW_INDICES)) DEALLOCATE(TRANSPOSE_ROW_INDICES)
    IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) DEALLOCATE(TRANSPOSE_COLUMN_INDICES)
    IF(ALLOCATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(TRANSPOSE_COLUMNS)) DEALLOCATE(TRANSPOSE_COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO local_row=1,SIZE(COLUMN_INDICES_LISTS,1)
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(local_row)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(local_row)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_row
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
    IF(ALLOCATED(TRANSPOSE_COLUMN_INDICES_LISTS)) THEN
      DO local_column=1,SIZE(TRANSPOSE_COLUMN_INDICES_LISTS,1)
        IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR)) &
          & CALL LIST_DESTROY(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_row
      DEALLOCATE(TRANSPOSE_COLUMN_INDICES_LISTS)
    ENDIF
998 CALL ERRORS("CONSTRAINT_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRIX_STRUCTURE_CALCULATE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRIX_STRUCTURE_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the constraint matrices for the constraint equations
  SUBROUTINE CONSTRAINT_MATRICES_CREATE_FINISH(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<The pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:),TRANSPOSE_ROW_INDICES(:),TRANSPOSE_COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)
    NULLIFY(TRANSPOSE_ROW_INDICES)
    NULLIFY(TRANSPOSE_COLUMN_INDICES)

    NULLIFY(ROW_DOMAIN_MAP)
    NULLIFY(COLUMN_DOMAIN_MAP)

    CALL ENTERS("CONSTRAINT_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have already been finished.",ERR,ERROR,*998)
      ELSE
        CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          COLUMN_DOMAIN_MAP=>CONSTRAINT_MAPPING%COLUMN_DOFS_MAPPING
          IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
            !Now create the individual constraint matrices
            DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
              CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                ROW_DOMAIN_MAP=>CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING
                IF(ASSOCIATED(ROW_DOMAIN_MAP)) THEN
                  !Create the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,CONSTRAINT_MATRICES% &
                    & MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(CONSTRAINT_MATRIX%MATRIX,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(CONSTRAINT_MATRIX%MATRIX,CONSTRAINT_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                  IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                    CALL DISTRIBUTED_MATRIX_CREATE_START(COLUMN_DOMAIN_MAP,ROW_DOMAIN_MAP,CONSTRAINT_MATRICES% &
                      & MATRICES(matrix_idx)%PTR%MATRIX_TRANSPOSE,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,CONSTRAINT_MATRIX%STORAGE_TYPE, &
                      & ERR,ERROR,*999)
                  ENDIF
                  !Calculate and set the matrix structure/sparsity pattern
                  IF(CONSTRAINT_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                    & CONSTRAINT_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                    CALL CONSTRAINT_MATRIX_STRUCTURE_CALCULATE(CONSTRAINT_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                      & TRANSPOSE_ROW_INDICES,TRANSPOSE_COLUMN_INDICES,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(CONSTRAINT_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(CONSTRAINT_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                      & ERR,ERROR,*999)
                    IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                      CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,NUMBER_OF_NON_ZEROS, &
                        & ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,TRANSPOSE_ROW_INDICES, &
                        & TRANSPOSE_COLUMN_INDICES,ERR,ERROR,*999)
                     ENDIF
                    IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                    IF(ASSOCIATED(TRANSPOSE_ROW_INDICES)) DEALLOCATE(TRANSPOSE_ROW_INDICES)
                    IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) DEALLOCATE(TRANSPOSE_COLUMN_INDICES)
                  ENDIF
                  CALL DISTRIBUTED_MATRIX_CREATE_FINISH(CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
                  IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                    CALL DISTRIBUTED_MATRIX_CREATE_FINISH(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Row domain map for constraint matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Constraint matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
            RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              !Set up the constraint RHS vector
              CALL DISTRIBUTED_VECTOR_CREATE_START(COLUMN_DOMAIN_MAP,CONSTRAINT_MATRICES%RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(RHS_VECTOR%RHS_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
            ENDIF
            !Finish up
            CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Column domain map is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(TRANSPOSE_ROW_INDICES)) DEALLOCATE(TRANSPOSE_ROW_INDICES)
    IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) DEALLOCATE(TRANSPOSE_COLUMN_INDICES)
    CALL CONSTRAINT_MATRICES_FINALISE(CONSTRAINT_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of the constraint matrices and rhs for the constraint equations
  SUBROUTINE CONSTRAINT_MATRICES_CREATE_START(CONSTRAINT_EQUATIONS,CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<The pointer to the constraint equations to create the constraint equations matrices for
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<On return, a pointer to the constraint matrices being created. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables

    CALL ENTERS("CONSTRAINT_MATRICES_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN      
      IF(CONSTRAINT_EQUATIONS%CONSTRAINT_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
          CALL FLAG_ERROR("Constraint matrices is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(CONSTRAINT_MATRICES)
          !Initialise the constraint matrices
          CALL CONSTRAINT_MATRICES_INITIALISE(CONSTRAINT_EQUATIONS,ERR,ERROR,*999)
          !Return the pointer
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_CREATE_START")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_CREATE_START",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_CREATE_START")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRICES_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_DESTROY(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer the constraint matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CALL CONSTRAINT_MATRICES_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("CONSTRAINT_MATRICES_DESTROY")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_DESTROY",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MATRICES_DESTROY")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MATRICES_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the constraint matrices and deallocate all memory.
  SUBROUTINE CONSTRAINT_MATRICES_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
   
    CALL ENTERS("CONSTRAINT_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(ALLOCATED(CONSTRAINT_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(CONSTRAINT_MATRICES%MATRICES,1)
          CALL CONSTRAINT_MATRIX_FINALISE(CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(CONSTRAINT_MATRICES%MATRICES)
      ENDIF
      CALL CONSTRAINT_MATRICES_RHS_FINALISE(CONSTRAINT_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(CONSTRAINT_MATRICES)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the constraint matrices for the constraint equations.
  SUBROUTINE CONSTRAINT_MATRICES_INITIALISE(CONSTRAINT_EQUATIONS,ERR,ERROR,*)
    
     !Argument variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<A pointer to the constraint equations to initialise the constraint matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
      IF(ASSOCIATED(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES)) THEN
        CALL FLAG_ERROR("Constraint matrices is already associated for this constraint equations.",ERR,ERROR,*998)
      ELSE
        CONSTRAINT_MAPPING=>CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
            ALLOCATE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint equations constraint matrices.",ERR,ERROR,*999)
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%CONSTRAINT_EQUATIONS=>CONSTRAINT_EQUATIONS
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED=.FALSE.
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%CONSTRAINT_MAPPING=>CONSTRAINT_MAPPING
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%SOLVER_MAPPING)
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NUMBER_OF_COLUMNS=CONSTRAINT_MAPPING%NUMBER_OF_COLUMNS
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%TOTAL_NUMBER_OF_COLUMNS=CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NUMBER_OF_GLOBAL_COLUMNS=CONSTRAINT_MAPPING%NUMBER_OF_GLOBAL_COLUMNS
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%RHS_VECTOR)
            !Allocate and initialise the matrices
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES=CONSTRAINT_MAPPING%NUMBER_OF_CONSTRAINT_MATRICES
            ALLOCATE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%MATRICES(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES% &
              & NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
              NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL CONSTRAINT_MATRIX_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
            CALL CONSTRAINT_MATRICES_RHS_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Constraint mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations constraint mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MATRICES_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRICES_FINALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Outputs the constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_OUTPUT(ID,CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: RHS_VECTOR
    
    CALL ENTERS("CONSTRAINT_MATRICES_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"Constraint matrices:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Number of constraint matrices = ",CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES, &
          & ERR,ERROR,*999)
        DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
          CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
            CALL WRITE_STRING_VALUE(ID,"Constraint matrix : ",matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING(ID,"Standard matrix:",ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_OUTPUT(ID,CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
              CALL WRITE_STRING(ID,"Transposed matrix:",ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,ERR,ERROR,*999)
            ENDIF
         ELSE
            CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
        RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          CALL WRITE_STRING(ID,"Constraint RHS vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_OUTPUT")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_OUTPUT",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_OUTPUT")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRICES_OUTPUT
  
  !
  !================================================================================================================================
  !

  !>Finalises the constraint matrices RHS vector and deallocates all memory
  SUBROUTINE CONSTRAINT_MATRICES_RHS_FINALISE(RHS_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: RHS_VECTOR !<A pointer to the equation matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("CONSTRAINT_MATRICES_RHS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_VECTOR)) THEN
      IF(ASSOCIATED(RHS_VECTOR%RHS_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(RHS_VECTOR)
    ENDIF      
     
    CALL EXITS("CONSTRAINT_MATRICES_RHS_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_RHS_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_RHS_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_RHS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the constraint matrices RHS vector
  SUBROUTINE CONSTRAINT_MATRICES_RHS_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the equation matrices to initialise the rhs vector for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_RHS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
        RHS_MAPPING=>CONSTRAINT_MAPPING%RHS_MAPPING
        IF(ASSOCIATED(RHS_MAPPING)) THEN
          IF(ASSOCIATED(CONSTRAINT_MATRICES%RHS_VECTOR)) THEN
            CALL FLAG_ERROR("Constraint matrices RHS vector is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(CONSTRAINT_MATRICES%RHS_VECTOR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices RHS vector.",ERR,ERROR,*999)
            CONSTRAINT_MATRICES%RHS_VECTOR%UPDATE_VECTOR=.TRUE.
            CONSTRAINT_MATRICES%RHS_VECTOR%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(CONSTRAINT_MATRICES%RHS_VECTOR%RHS_VECTOR)
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(CONSTRAINT_MATRICES%RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices equation mapping is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_RHS_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRICES_RHS_FINALISE(CONSTRAINT_MATRICES%RHS_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_RHS_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_RHS_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_RHS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_STORAGE_TYPE_SET(CONSTRAINT_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th inteface matrices. \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes,CONSTRAINT_MATRICES_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STORAGE_TYPE,1)==CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES) THEN
          DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
            CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
              SELECT CASE(STORAGE_TYPE(matrix_idx))
              CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
              CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                CONSTRAINT_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
              CASE DEFAULT
                LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                  & " for constraint matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of constraint matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_STORAGE_TYPE_SET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRICES_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the constraint matrices.
  SUBROUTINE CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET(CONSTRAINT_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The structure type for the  matrix_idx'th constraint matrix \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes,CONSTRAINT_MATRICES_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STRUCTURE_TYPE,1)==CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES) THEN
          DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
            CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
              SELECT CASE(STRUCTURE_TYPE(matrix_idx))
              CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
                CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
              CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
                CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_FEM_STRUCTURE
              CASE DEFAULT
                LOCAL_ERROR="The specified strucutre type of "// &
                  & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for constraint matrix number "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of constraint matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Initialise the values of the constraint matrices to the given value e.g., 0.0_DP
  SUBROUTINE CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the values for
    REAL(DP), INTENT(IN) :: VALUE !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(CONSTRAINT_RHS_TYPE), POINTER :: RHS_VECTOR
    
    CALL ENTERS("CONSTRAINT_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      DO matrix_idx=1,CONSTRAINT_MATRICES%NUMBER_OF_CONSTRAINT_MATRICES
        CONSTRAINT_MATRIX=>CONSTRAINT_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
          IF(CONSTRAINT_MATRIX%UPDATE_MATRIX) THEN
            CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(CONSTRAINT_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,VALUE,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        IF(RHS_VECTOR%UPDATE_VECTOR) THEN
          CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(RHS_VECTOR%RHS_VECTOR,VALUE,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_MATRICES_VALUES_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_VALUES_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_VALUES_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_VALUES_INITIALISE

  !
  !================================================================================================================================
  !
  
  SUBROUTINE ConstraintMatrix_TimeDependenceTypeSet(ConstraintCondition, &
    & constraintMatrixIndex,IsTranspose,TimeDependenceType,Err,Error,*)
    
    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: ConstraintCondition
    INTEGER(INTG), INTENT(IN) :: ConstraintMatrixIndex
    LOGICAL, INTENT(IN) :: IsTranspose
    INTEGER(INTG), INTENT(IN) :: TimeDependenceType
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: ConstraintEquations
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: ConstraintMatrices
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: ConstraintMatrix
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ConstraintMatrix_TimeDependenceTypeSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(ConstraintCondition)) THEN
      ConstraintEquations=>ConstraintCondition%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(ConstraintEquations)) THEN
        ConstraintMatrices=>ConstraintEquations%CONSTRAINT_MATRICES
        IF(ASSOCIATED(ConstraintMatrices)) THEN
          ConstraintMatrix=>ConstraintMatrices%MATRICES(ConstraintMatrixIndex)%PTR
          IF(ASSOCIATED(ConstraintMatrix)) THEN
            IF(.NOT.IsTranspose) THEN
              ConstraintMatrix%CONSTRAINT_MATRIX_TIME_DEPENDENCE_TYPE=TimeDependenceType
            ELSE
              IF(ConstraintMatrix%HAS_TRANSPOSE) THEN
                ConstraintMatrix%CONSTRAINT_MATRIX_TRANSPOSE_TIME_DEPENDENCE_TYPE=TimeDependenceType
              ELSE
                LOCAL_ERROR="Constraint matrices has_transpose flag is .false. but constraint matrix type is transpose."
                CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint matrix is not associated",Err,Error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint equations not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",Err,Error,*999)
    ENDIF
    
    CALL EXITS("ConstraintMatrix_TimeDependenceTypeSet")
    RETURN
999 CALL ERRORS("ConstraintMatrix_TimeDependenceTypeSet",Err,Error)
    CALL EXITS("ConstraintMatrix_TimeDependenceTypeSet")
    RETURN 1
  END SUBROUTINE ConstraintMatrix_TimeDependenceTypeSet

  !
  !================================================================================================================================
  !
  
  SUBROUTINE ConstraintMatrix_TimeDependenceTypeGet(ConstraintCondition, &
    & constraintMatrixIndex,IsTranspose,TimeDependenceType,Err,Error,*)
    
    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: ConstraintCondition
    INTEGER(INTG), INTENT(IN) :: ConstraintMatrixIndex
    LOGICAL, INTENT(IN) :: IsTranspose
    INTEGER(INTG), INTENT(OUT) :: TimeDependenceType
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: ConstraintEquations
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: ConstraintMatrices
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: ConstraintMatrix
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ConstraintMatrix_TimeDependenceTypeGet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(ConstraintCondition)) THEN
      ConstraintEquations=>ConstraintCondition%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(ConstraintEquations)) THEN
        ConstraintMatrices=>ConstraintEquations%CONSTRAINT_MATRICES
        IF(ASSOCIATED(ConstraintMatrices)) THEN
          ConstraintMatrix=>ConstraintMatrices%MATRICES(ConstraintMatrixIndex)%PTR
          IF(ASSOCIATED(ConstraintMatrix)) THEN
            IF(.NOT.IsTranspose) THEN
              TimeDependenceType=ConstraintMatrix%CONSTRAINT_MATRIX_TIME_DEPENDENCE_TYPE
            ELSE
              IF(ConstraintMatrix%HAS_TRANSPOSE) THEN
                TimeDependenceType=ConstraintMatrix%CONSTRAINT_MATRIX_TRANSPOSE_TIME_DEPENDENCE_TYPE
              ELSE
                LOCAL_ERROR="Constraint matrices has_transpose flag is .false. but constraint matrix type is transpose."
                CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
              ENDIF
            ENDIF
            !Sanity check
            IF(.NOT.(TimeDependenceType>0.AND.TimeDependenceType<=NUMBER_OF_CONSTRAINT_MATRIX_TYPES)) THEN
              LOCAL_ERROR="Invalid time dependence type of "//TRIM(NUMBER_TO_VSTRING(TimeDependenceType,"*",ERR,ERROR))// &
                & ". Must be > 0 and <= "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_CONSTRAINT_MATRIX_TYPES,"*",ERR,ERROR))// &
                & "."
              CALL FLAG_ERROR(LOCAL_ERROR,Err,Error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint matrix is not associated",Err,Error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint equations not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",Err,Error,*999)
    ENDIF
    
    CALL EXITS("ConstraintMatrix_TimeDependenceTypeGet")
    RETURN
999 CALL ERRORS("ConstraintMatrix_TimeDependenceTypeGet",Err,Error)
    CALL EXITS("ConstraintMatrix_TimeDependenceTypeGet")
    RETURN 1
  END SUBROUTINE ConstraintMatrix_TimeDependenceTypeGet

  !
  !================================================================================================================================
  !


END MODULE CONSTRAINT_MATRICES_ROUTINES
