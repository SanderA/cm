!> \file
!> \author Chris Bradley
!> \brief This module handles all constraint matrix and rhs routines.
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

!> This module handles all constraint matrix and rhs routines.
MODULE CONSTRAINT_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE CONSTRAINT_CONDITION_CONSTANTS
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES
  USE LINKEDLIST_ROUTINES
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes CONSTRAINT_MATRICES_ROUTINES::ConstraintMatrixStructureTypes
  !> \brief Constraint matrices structure (sparsity) types
  !> \see CONSTRAINT_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE=3 !<Diagonal matrix structure. \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatrixStructureTypes,CONSTRAINT_MATRICES_ROUTINES
  !>@}

  !> \addtogroup CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes CONSTRAINT_MATRICES_ROUTINES::ConstraintMatricesSparsityTypes
  !> \brief Constraint matrices sparsity types
  !> \see CONSTRAINT_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_SPARSE_MATRICES=1 !<Use sparse constraint matrices \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_FULL_MATRICES=2 !<Use fully populated constraint matrices \see CONSTRAINT_MATRICES_ROUTINES_ConstraintMatricesSparsityTypes,CONSTRAINT_MATRICES_ROUTINES
  !>@}

  !> \addtogroup CONSTRAINT_MATRICES_ROUTINES_SelectMatricesTypes CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes
  !> \brief The types of selection available for the constraint matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_ALL=1 !<Select all the constraint matrices and vectors \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic constraint matrices and vectors \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_LINEAR_ONLY=3 !<Select only the linear constraint matrices and vectors \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear constraint matrices and vectors \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian constraint matrix \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual constraint vector \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_RHS_ONLY=7 !<Select only the RHS constraint vector \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_RHS_RESIDUAL_ONLY=9 !<Select only the RHS and residual constraint vectors \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_MATRICES_VECTORS_ONLY=12 !<Assemble only the constraint vectors \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
  !>@}

  !> \addtogroup CONSTRAINT_MATRICES_ROUTINES_JacobianCalculationTypes CONSTRAINT_MATRICES_ROUTINES:JacobianCalculationTypes
  !> \brief Jacobian calculation types
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_JACOBIAN_FINITE_DIFFERENCE_CALCULATED=1 !<Use finite differencing to calculate the Jacobian
  INTEGER(INTG), PARAMETER :: CONSTRAINT_JACOBIAN_ANALYTIC_CALCULATED=2 !<Use an analytic Jacobian evaluation
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  INTERFACE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET
    MODULE PROCEDURE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0
    MODULE PROCEDURE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1
  END INTERFACE !CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET

  INTERFACE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET
    MODULE PROCEDURE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0
    MODULE PROCEDURE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1
  END INTERFACE !CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET

  PUBLIC CONSTRAINT_MATRIX_NO_STRUCTURE,CONSTRAINT_MATRIX_FEM_STRUCTURE,CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE

  PUBLIC CONSTRAINT_MATRICES_SPARSE_MATRICES,CONSTRAINT_MATRICES_FULL_MATRICES

  PUBLIC CONSTRAINT_JACOBIAN_FINITE_DIFFERENCE_CALCULATED,CONSTRAINT_JACOBIAN_ANALYTIC_CALCULATED

  PUBLIC CONSTRAINT_MATRICES_CREATE_FINISH,CONSTRAINT_MATRICES_CREATE_START,CONSTRAINT_MATRICES_DESTROY

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC CONSTRAINT_MATRICES_ELEMENT_ADD,CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD,CONSTRAINT_MATRICES_ELEMENT_CALCULATE, &
    & CONSTRAINT_MATRICES_ELEMENT_INITIALISE,CONSTRAINT_MATRICES_ELEMENT_FINALISE,CONSTRAINT_MATRICES_VALUES_INITIALISE

  PUBLIC CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE,CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE, &
    & CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE,CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP

  PUBLIC ConstraintMatrices_JacobianTypesSet

  PUBLIC CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE,CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE, &
    & CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE,CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP

  PUBLIC CONSTRAINT_MATRICES_ALL,CONSTRAINT_MATRICES_LINEAR_ONLY,CONSTRAINT_MATRICES_NONLINEAR_ONLY,CONSTRAINT_MATRICES_JACOBIAN_ONLY, &
    & CONSTRAINT_MATRICES_RESIDUAL_ONLY,CONSTRAINT_MATRICES_RHS_ONLY,CONSTRAINT_MATRICES_RHS_RESIDUAL_ONLY, &
    & CONSTRAINT_MATRICES_VECTORS_ONLY
  
  PUBLIC CONSTRAINT_MATRICES_OUTPUT,CONSTRAINT_MATRICES_JACOBIAN_OUTPUT

  PUBLIC CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET,CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET

  PUBLIC CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET,CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET

  PUBLIC CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET,CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalise the constraint Jacobian and deallocate all memory
  SUBROUTINE CONSTRAINT_JACOBIAN_FINALISE(CONSTRAINT_JACOBIAN,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: CONSTRAINT_JACOBIAN !<A pointer to the constraint Jacobian to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_JACOBIAN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_JACOBIAN)) THEN
      IF(ASSOCIATED(CONSTRAINT_JACOBIAN%JACOBIAN)) CALL DISTRIBUTED_MATRIX_DESTROY(CONSTRAINT_JACOBIAN%JACOBIAN,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE)) CALL DISTRIBUTED_MATRIX_DESTROY(CONSTRAINT_JACOBIAN%JACOBIAN_TRANSPOSE,ERR,ERROR,*999)
      CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE(CONSTRAINT_JACOBIAN%ELEMENT_JACOBIAN,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_JACOBIAN_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_JACOBIAN_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_JACOBIAN_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_JACOBIAN_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the constraint Jacobian.
  SUBROUTINE CONSTRAINT_JACOBIAN_INITIALISE(NONLINEAR_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer to the constraint matrices nonlinear matrices to initialise the Jacobian for
    INTEGER(INTG), INTENT(IN) :: MATRIX_NUMBER !<The index of the Jacobian matrix to initialise for the nonlinear matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("CONSTRAINT_JACOBIAN_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
      CONSTRAINT_MATRICES=>NONLINEAR_MATRICES%CONSTRAINT_MATRICES
      IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
        CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            IF(ALLOCATED(NONLINEAR_MATRICES%JACOBIANS)) THEN
              IF(ASSOCIATED(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR)) THEN
                CALL FLAG_ERROR("Nonlinear matrices Jacobian is already associated.",ERR,ERROR,*998)
              ELSE
                ALLOCATE(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint Jacobian.",ERR,ERROR,*999)
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN_NUMBER=MATRIX_NUMBER
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%NONLINEAR_MATRICES=>NONLINEAR_MATRICES
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%NUMBER_OF_COLUMNS= &
                    & NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%UPDATE_JACOBIAN=.TRUE.
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%FIRST_ASSEMBLY=.TRUE.
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%HAS_TRANSPOSE=NONLINEAR_MAPPING% &
                  & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(MATRIX_NUMBER)%HAS_TRANSPOSE
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%NUMBER_OF_ROWS=NONLINEAR_MAPPING% &
                  & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(MATRIX_NUMBER)%NUMBER_OF_ROWS
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%TOTAL_NUMBER_OF_ROWS=NONLINEAR_MAPPING% &
                  & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(MATRIX_NUMBER)%TOTAL_NUMBER_OF_ROWS
                NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(MATRIX_NUMBER)%CONSTRAINT_JACOBIAN=> &
                  & NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)
                NULLIFY(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN)
                NULLIFY(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN_TRANSPOSE)
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR% &
                    & ELEMENT_JACOBIAN,ERR,ERROR,*999)
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN_CALCULATION_TYPE= &
                  & CONSTRAINT_JACOBIAN_FINITE_DIFFERENCE_CALCULATED
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint matrices nonlinear matrieces Jacobian is not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Nonlinear matrices constraint matrices is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Nonlinear matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_JACOBIAN_INITIALISE")
    RETURN
999 CALL CONSTRAINT_JACOBIAN_FINALISE(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_JACOBIAN_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_JACOBIAN_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_JACOBIAN_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the constraint matrices and RHS for the the equations
  SUBROUTINE CONSTRAINT_MATRICES_CREATE_FINISH(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<The pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    type(LinkedList),pointer :: list(:) 
    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)

    CALL ENTERS("CONSTRAINT_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have already been finished.",ERR,ERROR,*998)
      ELSE
        CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          ROW_DOMAIN_MAP=>CONSTRAINT_MAPPING%ROW_DOFS_MAPPING
          IF(ASSOCIATED(ROW_DOMAIN_MAP)) THEN
            DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
            IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
              !Dynamic matrices
              DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                !Now create the individual dynamic constraint matrices
                DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
                  CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
                  IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                    COLUMN_DOMAIN_MAP=>DYNAMIC_MAPPING%CONSTRAINT_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
                      !Create the distributed constraint matrix
                      CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,CONSTRAINT_MATRICES% &
                           & DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(CONSTRAINT_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(CONSTRAINT_MATRIX%MATRIX,CONSTRAINT_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                      !Calculate and set the matrix structure/sparsity pattern
                      IF(CONSTRAINT_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                        & CONSTRAINT_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                        CALL ConstraintMatrix_StructureCalculate(CONSTRAINT_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                          & list,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_LINKLIST_SET(CONSTRAINT_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(CONSTRAINT_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(CONSTRAINT_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                        IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                      ENDIF
                      CALL DISTRIBUTED_MATRIX_CREATE_FINISH(CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Column domain map for dynamic matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Constraint matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                      & " is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx                
              ELSE
                CALL FLAG_ERROR("Constraint mapping dynamic mapping is not associated.",ERR,ERROR,*999)                
              ENDIF
            ENDIF
            LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
            IF(ASSOCIATED(LINEAR_MATRICES)) THEN
              !Linear matrices
              LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
              IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                !Now create the individual linear constraint matrices
                DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                  CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
                  IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                    COLUMN_DOMAIN_MAP=>LINEAR_MAPPING%CONSTRAINT_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
                      !Create the distributed constraint matrix
                      CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,CONSTRAINT_MATRICES% &
                           & LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(CONSTRAINT_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(CONSTRAINT_MATRIX%MATRIX,CONSTRAINT_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                      !Calculate and set the matrix structure/sparsity pattern
                      IF(CONSTRAINT_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                        & CONSTRAINT_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                        CALL ConstraintMatrix_StructureCalculate(CONSTRAINT_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                          & list,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_LINKLIST_SET(CONSTRAINT_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(CONSTRAINT_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(CONSTRAINT_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                        IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                      ENDIF
                      CALL DISTRIBUTED_MATRIX_CREATE_FINISH(CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Column domain map for linear matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Constraint matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                      & " is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
              ELSE
                CALL FLAG_ERROR("Constraint mapping linear mapping is not associated.",ERR,ERROR,*999)                
              ENDIF
            ENDIF
            NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
            IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
              !Nonlinear matrices
              NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
              IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                !Set up the Jacobian matrices
                DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
                  JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
                  IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                    COLUMN_DOMAIN_MAP=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(matrix_idx)%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
  !!TODO: Set the distributed matrix not to allocate the data if the Jacobian is not calculated.
                      !Create the distributed Jacobian matrix
                      CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(JACOBIAN_MATRIX%JACOBIAN,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(JACOBIAN_MATRIX%JACOBIAN,JACOBIAN_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                      !Calculate and set the matrix structure/sparsity pattern
                      IF(JACOBIAN_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                        & JACOBIAN_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                        CALL JacobianMatrix_StructureCalculate(JACOBIAN_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(JACOBIAN_MATRIX%JACOBIAN,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(JACOBIAN_MATRIX%JACOBIAN,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                        IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                      ENDIF
                      CALL DISTRIBUTED_MATRIX_CREATE_FINISH(JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Column domain map is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Jacobian matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO
                !Set up the residual vector                
                CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,CONSTRAINT_MATRICES%NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NONLINEAR_MATRICES%RESIDUAL,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_CREATE_FINISH(NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
                !Initialise the residual vector to zero for time dependent problems so that the previous residual is set to zero
                CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(NONLINEAR_MATRICES%RESIDUAL,0.0_DP,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Constraint mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDIF
            RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              !Set up the equations RHS vector          
              CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,CONSTRAINT_MATRICES%RHS_VECTOR%VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(RHS_VECTOR%VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(RHS_VECTOR%VECTOR,ERR,ERROR,*999)
            ENDIF
            !Finish up
            CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Row domain map is not associated.",ERR,ERROR,*999)
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
    CALL CONSTRAINT_MATRICES_FINALISE(CONSTRAINT_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of the constraint matrices and rhs for the the equations
  SUBROUTINE CONSTRAINT_MATRICES_CREATE_START(EQUATIONS,CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<The pointer to the equations to create the constraint matrices for
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<On return, a pointer to the constraint matrices being created.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR    

    CALL ENTERS("CONSTRAINT_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN      
      IF(CONSTRAINT_EQUATIONS%EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
          CALL FLAG_ERROR("Constraint matrices is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(CONSTRAINT_MATRICES)
          !Initialise the constraint matrices
          CALL CONSTRAINT_MATRICES_INITIALISE(EQUATIONS,ERR,ERROR,*999)
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_CREATE_START")
    RETURN
999 CALL CONSTRAINT_MATRICES_FINALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_CREATE_START",ERR,ERROR)
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
      CALL FLAG_ERROR("Constraint matrices is not associated",ERR,ERROR,*999)
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

  !>Calculate the positions in the constraint matrices of the element matrix. Old CMISS name MELGE.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE(ELEMENT_MATRIX,UPDATE_MATRIX,ROW_ELEMENT_NUMBERS,COLUMN_ELEMENT_NUMBERS, &
    & ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix to calculate
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if the element matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: ROW_ELEMENT_NUMBERS(:) !<The row element number to calculate
    INTEGER(INTG), INTENT(IN) :: COLUMN_ELEMENT_NUMBERS(:) !<The column element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLS_FIELD_VARIABLE !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,global_ny,local_ny,node,node_idx,version,dataPointIdx, &
      & localDataPointNumber,elementIdx,rowElementNumber,colElementNumber
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      IF(ASSOCIATED(COLS_FIELD_VARIABLE)) THEN
        ELEMENT_MATRIX%NUMBER_OF_ROWS=0
        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
        IF(UPDATE_MATRIX) THEN
          IF(ASSOCIATED(ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE)) THEN
            !Row and columns variable is the same.
            DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              DO elementIdx=1,SIZE(ROW_ELEMENT_NUMBERS)
                rowElementNumber=ROW_ELEMENT_NUMBERS(elementIdx)
                IF(rowElementNumber>=1.AND.rowElementNumber<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                      & ELEMENTS(rowElementNumber)
                    global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%BASIS
                    DO node_idx=1,BASIS%NUMBER_OF_NODES
                      node=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_NODES(node_idx)
                      DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                        derivative=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                        version=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%elementVersions(derivative_idx,node_idx)
                        local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                          & DERIVATIVES(derivative)%VERSIONS(version)
                        global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                        ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                        ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                      ENDDO !derivative_idx
                    ENDDO !node_idx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                    decompositionData=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
                    DO dataPointIdx=1,decompositionData%elementDataPoint(rowElementNumber)%numberOfProjectedData
                      localDataPointNumber=decompositionData%elementDataPoint(rowElementNumber)% &
                        & dataIndices(dataPointIdx)%localNumber
                      local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                        & DATA_POINTS(localDataPointNumber)
                      global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "// &
                      & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                      & " of rows field variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                  END SELECT
                ELSE
                  LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(rowElementNumber,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of rows field variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elementIdx
            ENDDO !component_idx
          ELSE
            !Row and column variables are different
            !Row mapping
            DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              DO elementIdx=1,SIZE(ROW_ELEMENT_NUMBERS)
                rowElementNumber=ROW_ELEMENT_NUMBERS(elementIdx)
                IF(rowElementNumber>=1.AND.rowElementNumber<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                      & ELEMENTS(rowElementNumber)
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%BASIS
                    DO node_idx=1,BASIS%NUMBER_OF_NODES
                      node=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_NODES(node_idx)
                      DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                        derivative=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                        version=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%elementVersions(derivative_idx,node_idx)
                        local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                          & DERIVATIVES(derivative)%VERSIONS(version)
                        ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                        ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                      ENDDO !derivative_idx
                    ENDDO !node_idx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                    decompositionData=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
                    DO dataPointIdx=1,decompositionData%elementDataPoint(colElementNumber)%numberOfProjectedData
                      localDataPointNumber=decompositionData%elementDataPoint(colElementNumber)% &
                        & dataIndices(dataPointIdx)%localNumber
                      local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                        & DATA_POINTS(localDataPointNumber)
                      global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "// &
                      & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                      & " of rows field variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                  END SELECT
                ELSE
                  LOCAL_ERROR="Row element number "//TRIM(NUMBER_TO_VSTRING(rowElementNumber,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of rows field variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elementIdx
            ENDDO !component_idx
            !Column mapping
            DO component_idx=1,COLS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              DO elementIdx=1,SIZE(COLUMN_ELEMENT_NUMBERS)
                colElementNumber=COLUMN_ELEMENT_NUMBERS(elementIdx)
                IF(colElementNumber>=1.AND.colElementNumber<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  SELECT CASE(COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                      & ELEMENTS(colElementNumber)
                    global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%BASIS
                    DO node_idx=1,BASIS%NUMBER_OF_NODES
                      node=ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%ELEMENT_NODES(node_idx)
                      DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                        derivative=ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                        version=ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%elementVersions(derivative_idx,node_idx)
                        local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                          & DERIVATIVES(derivative)%VERSIONS(version)
                        global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                        ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                      ENDDO !derivative_idx
                    ENDDO !node_idx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                    decompositionData=>COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
                    DO dataPointIdx=1,decompositionData%elementDataPoint(colElementNumber)%numberOfProjectedData
                      localDataPointNumber=decompositionData%elementDataPoint(colElementNumber)% &
                        & dataIndices(dataPointIdx)%localNumber
                      local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                        & DATA_POINTS(localDataPointNumber)
                      global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "// &
                      & TRIM(NUMBER_TO_VSTRING(COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                      & " of column field variable type "// &
                      & TRIM(NUMBER_TO_VSTRING(COLS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                  END SELECT
                ELSE
                  LOCAL_ERROR="Column element number "//TRIM(NUMBER_TO_VSTRING(colElementNumber,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of column field variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(COLS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elementIdx
            ENDDO !component_idx
          ENDIF
          ELEMENT_MATRIX%MATRIX=0.0_DP
        ENDIF
      ELSE
        CALL FLAG_ERROR("Columns field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise an element matrix and deallocate all memory
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE):: ELEMENT_MATRIX !<The element matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR,*999)

    ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
    IF(ALLOCATED(ELEMENT_MATRIX%ROW_DOFS)) DEALLOCATE(ELEMENT_MATRIX%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%COLUMN_DOFS)) DEALLOCATE(ELEMENT_MATRIX%COLUMN_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%MATRIX)) DEALLOCATE(ELEMENT_MATRIX%MATRIX)
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element matrix.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !The element matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE",ERR,ERROR,*999)

    ELEMENT_MATRIX%CONSTRAINT_MATRIX_NUMBER=0
    ELEMENT_MATRIX%NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
       
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets up the element matrix for the row and column field variables.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP(elementMatrix,rowsFieldVariable,columnsFieldVariable, &
    & rowsNumberOfElements,colsNumberOfElements,err,error,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: elementMatrix !<The element matrix to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: columnsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(IN)  :: rowsNumberOfElements !<Number of elements in the row variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(IN)  :: colsNumberOfElements !<Number of elements in the col variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr, componentIdx
    TYPE(VARYING_STRING) :: dummyError

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP",err,error,*998)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      IF(ASSOCIATED(columnsFieldVariable)) THEN
        elementMatrix%MAX_NUMBER_OF_ROWS = 0
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          elementMatrix%MAX_NUMBER_OF_ROWS=elementMatrix%MAX_NUMBER_OF_ROWS+ &
            & rowsFieldVariable%COMPONENTS(componentIdx)%maxNumberElementInterpolationParameters
        ENDDO
        elementMatrix%MAX_NUMBER_OF_ROWS=elementMatrix%MAX_NUMBER_OF_ROWS*rowsNumberOfElements
        elementMatrix%MAX_NUMBER_OF_COLUMNS = 0
        DO componentIdx=1,columnsFieldVariable%NUMBER_OF_COMPONENTS
          elementMatrix%MAX_NUMBER_OF_COLUMNS=elementMatrix%MAX_NUMBER_OF_COLUMNS+ &
            & columnsFieldVariable%COMPONENTS(componentIdx)%maxNumberElementInterpolationParameters
        ENDDO
        elementMatrix%MAX_NUMBER_OF_COLUMNS=elementMatrix%MAX_NUMBER_OF_COLUMNS*colsNumberOfElements
        IF(ALLOCATED(elementMatrix%ROW_DOFS)) THEN
          CALL FLAG_ERROR("Element matrix row dofs already allocated.",err,error,*999)
        ELSE
          ALLOCATE(elementMatrix%ROW_DOFS(elementMatrix%MAX_NUMBER_OF_ROWS),STAT=err)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix row dofs.",err,error,*999)
        ENDIF
        IF(ALLOCATED(elementMatrix%COLUMN_DOFS)) THEN
          CALL FLAG_ERROR("Element matrix column dofs already allocated.",err,error,*999)
        ELSE
          ALLOCATE(elementMatrix%COLUMN_DOFS(elementMatrix%MAX_NUMBER_OF_COLUMNS),STAT=err)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix column dofs.",err,error,*999)
        ENDIF
        IF(ALLOCATED(elementMatrix%MATRIX)) THEN
          CALL FLAG_ERROR("Element matrix already allocated.",err,error,*999)
        ELSE
          ALLOCATE(elementMatrix%MATRIX(elementMatrix%MAX_NUMBER_OF_ROWS,elementMatrix%MAX_NUMBER_OF_COLUMNS),STAT=err)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Columns field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP")
    RETURN
999 CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE(elementMatrix,dummyErr,dummyError,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP",err,error)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations rhs of the element rhs vector. Old CMISS name MELGE.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE(ELEMENT_VECTOR,UPDATE_VECTOR,ELEMENT_NUMBER,ROWS_FIELD_VARIABLE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element vector to calculate.
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the element vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,local_ny,node,node_idx,version,dataPointIdx,localDataPointNumber
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      !Calculate the rows for the element vector
      ELEMENT_VECTOR%NUMBER_OF_ROWS=0
      IF(UPDATE_VECTOR) THEN
        DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
          IF(ELEMENT_NUMBER>=1.AND.ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
            SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
              ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                & ELEMENTS(ELEMENT_NUMBER)
              ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
              ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO node_idx=1,BASIS%NUMBER_OF_NODES
                node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                  derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                  version=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%elementVersions(derivative_idx,node_idx)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & DERIVATIVES(derivative)%VERSIONS(version)
                  ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                  ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
                ENDDO !derivative_idx
              ENDDO !node_idx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              decompositionData=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
              DO dataPointIdx=1,decompositionData%elementDataPoint(ELEMENT_NUMBER)%numberOfProjectedData
                localDataPointNumber=decompositionData%elementDataPoint(ELEMENT_NUMBER)% &
                  & dataIndices(dataPointIdx)%localNumber
                local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                  & DATA_POINTS(localDataPointNumber)
                ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
              ENDDO
            CASE DEFAULT
              LOCAL_ERROR="The interpolation type of "// &
                & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                & " is invalid for component number "// &
                & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " of rows field variable type "// &
                & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
            END SELECT
          ELSE
            LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
              & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " of rows field variable type "//TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
              & ". The element number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !component_idx
        ELEMENT_VECTOR%VECTOR=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise an element vector and deallocate all memory
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE):: ELEMENT_VECTOR !<The element vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(ELEMENT_VECTOR%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_VECTOR%VECTOR)) DEALLOCATE(ELEMENT_VECTOR%VECTOR)
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE")

    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element vector
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !The element vector to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE",ERR,ERROR,*999)

    ELEMENT_VECTOR%NUMBER_OF_ROWS=0
    ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
       
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets up the element vector for the row field variables.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP(elementVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: elementVector !<The element vector to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,componentIdx
    TYPE(VARYING_STRING) :: dummyError

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP",err,error,*998)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      elementVector%MAX_NUMBER_OF_ROWS = 0
      DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
        elementVector%MAX_NUMBER_OF_ROWS=elementVector%MAX_NUMBER_OF_ROWS+ &
          & rowsFieldVariable%COMPONENTS(componentIdx)%maxNumberElementInterpolationParameters
      ENDDO
      IF(ALLOCATED(elementVector%ROW_DOFS)) THEN
        CALL FLAG_ERROR("Element vector row dofs is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(elementVector%ROW_DOFS(elementVector%MAX_NUMBER_OF_ROWS),STAT=err)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector row dofs.",err,error,*999)
      ENDIF
      IF(ALLOCATED(elementVector%VECTOR)) THEN
        CALL FLAG_ERROR("Element vector vector already allocated.",err,error,*999)
      ELSE
        ALLOCATE(elementVector%VECTOR(elementVector%MAX_NUMBER_OF_ROWS),STAT=err)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector vector.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP")
    RETURN
999 CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE(elementVector,DUMMY_ERR,dummyError,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP",err,error)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP

  !
  !================================================================================================================================
  !

  !>Adds the element matrices and rhs vector into the constraint matrices and rhs vector.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,matrix_idx,row_idx
    REAL(DP) :: SUM
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_ELEMENT_ADD()")
#endif

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
      IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
        !Add the element matrices
        DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
          CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
            IF(CONSTRAINT_MATRIX%UPDATE_MATRIX) THEN
              !Add the element matrice into the distributed constraint matrix
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(CONSTRAINT_MATRIX%MATRIX,CONSTRAINT_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
                & ERR,ERROR,*999)
            ENDIF
            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,CONSTRAINT_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),TRANSPOSE(CONSTRAINT_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS)), &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Constraint matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
      IF(ASSOCIATED(LINEAR_MATRICES)) THEN
        !Add the element matrices
        DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
          CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
            IF(CONSTRAINT_MATRIX%UPDATE_MATRIX) THEN
              !Add the element matrice into the distributed constraint matrix
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(CONSTRAINT_MATRIX%MATRIX,CONSTRAINT_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
                & ERR,ERROR,*999)
            ENDIF
            !If the constraint matrix has a transpose add it
            IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE,CONSTRAINT_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),CONSTRAINT_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),TRANSPOSE(CONSTRAINT_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                & CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:CONSTRAINT_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS)), &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Constraint matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
          !Add the residual element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(NONLINEAR_MATRICES%RESIDUAL,NONLINEAR_MATRICES%ELEMENT_RESIDUAL%ROW_DOFS(1: &
            & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%NUMBER_OF_ROWS),NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(1:NONLINEAR_MATRICES% &
            & ELEMENT_RESIDUAL%NUMBER_OF_ROWS),ERR,ERROR,*999)
        ENDIF
      ENDIF
      RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        IF(RHS_VECTOR%UPDATE_VECTOR) THEN
          !Add the rhs element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(RHS_VECTOR%VECTOR,RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS(1: &
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

  !>Calculate the positions in the constraint matrices and rhs of the element matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
        DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          !Calculate the row and columns for the dynamic constraint matrices
          DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                ROW_FIELD_VARIABLE=>DYNAMIC_MAPPING%CONSTRAINT_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE
                COL_FIELD_VARIABLE=>DYNAMIC_MAPPING%LAGRANGE_VARIABLE
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,CONSTRAINT_MATRIX%UPDATE_MATRIX, &
                  & [ELEMENT_NUMBER],[ELEMENT_NUMBER],ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Constraint matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Constraint mapping dynamic mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          !Calculate the row and columns for the linear constraint matrices
          LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                ROW_FIELD_VARIABLE=>LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE
                COL_FIELD_VARIABLE=>LINEAR_MAPPING%LAGRANGE_VARIABLE
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,CONSTRAINT_MATRIX%UPDATE_MATRIX, &
                  & [ELEMENT_NUMBER],[ELEMENT_NUMBER],ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Constraint matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Constraint mapping linear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          !Calculate the rows and columns of the Jacobian
          NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                ROW_FIELD_VARIABLE=>NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE
                COL_FIELD_VARIABLE=>NONLINEAR_MAPPING%LAGRANGE_VARIABLE
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_CALCULATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN,JACOBIAN_MATRIX%UPDATE_JACOBIAN, &
                  & [ELEMENT_NUMBER],[ELEMENT_NUMBER],ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
            !Calculate the rows of the equations residual
            ROW_FIELD_VARIABLE=>NONLINEAR_MAPPING%LAGRANGE_VARIABLE
            CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,NONLINEAR_MATRICES% &
              & UPDATE_RESIDUAL,ELEMENT_NUMBER,ROW_FIELD_VARIABLE,ERR,ERROR,*999)
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL_CALCULATED=0
          ELSE
            CALL FLAG_ERROR("Constraint mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          RHS_MAPPING=>CONSTRAINT_MAPPING%RHS_MAPPING
          IF(ASSOCIATED(RHS_MAPPING)) THEN
            !Calculate the rows  for the equations RHS
            ROW_FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_CALCULATE(RHS_VECTOR%ELEMENT_VECTOR,RHS_VECTOR%UPDATE_VECTOR,ELEMENT_NUMBER, &
              & ROW_FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Constraint mapping rhs mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_ELEMENT_CALCULATE()")
#endif
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information and deallocate all memory
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<The constraint matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
      IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
        !Finalise the dynamic element matrices
        DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
          CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
            CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Constraint matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
      IF(ASSOCIATED(LINEAR_MATRICES)) THEN
        !Finalise the linear element matrices
        DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
          CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
            CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Constraint matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MAX_NUMBER_OF_ROWS=0
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MAX_NUMBER_OF_COLUMNS=0
            IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS)
            IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS)
            IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX)
          ELSE
            CALL FLAG_ERROR("Nonlinear matrices Jacobian number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                & " is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO
        NONLINEAR_MATRICES%ELEMENT_RESIDUAL%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%ROW_DOFS)) DEALLOCATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%ROW_DOFS)
        IF(ALLOCATED(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR)) DEALLOCATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR)
      ENDIF
      RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        !Finalise the element vector
        RHS_VECTOR%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)
      ENDIF
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
  SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !The constraint matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    INTEGER(INTG) :: rowsNumberOfElements,colsNumberOfElements !Number of elements in the row and col variables whose dofs are present in the element matrix
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_ELEMENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      rowsNumberOfElements=1
      colsNumberOfElements=1
      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
        DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          !Initialise the dynamic element matrices
          DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                ROW_FIELD_VARIABLE=>DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE
                COL_FIELD_VARIABLE=>DYNAMIC_MAPPING%LAGRANGE_VARIABLE
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP(CONSTRAINT_MATRIX%ELEMENT_MATRIX, &
                  & ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE,rowsNumberOfElements,colsNumberOfElements,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Constraint dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Constraint mapping dynamic mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          !Initialise the linear element matrices
          LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                ROW_FIELD_VARIABLE=>LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE
                COL_FIELD_VARIABLE=>LINEAR_MAPPING%LAGRANGE_VARIABLE
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP(CONSTRAINT_MATRIX%ELEMENT_MATRIX, &
                  & ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE,rowsNumberOfElements,colsNumberOfElements,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Constraint linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Constraint mapping linear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          !Initialise the Jacobian element matrices
          NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                ROW_FIELD_VARIABLE=>NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE
                COL_FIELD_VARIABLE=>NONLINEAR_MAPPING%LAGRANGE_VARIABLE
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_SETUP(JACOBIAN_MATRIX%ELEMENT_JACOBIAN, &
                  & ROW_FIELD_VARIABLE,COL_FIELD_VARIABLE,rowsNumberOfElements,colsNumberOfElements,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
            ROW_FIELD_VARIABLE=>NONLINEAR_MAPPING%LAGRANGE_VARIABLE
            CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ROW_FIELD_VARIABLE,ERR,ERROR,*999)
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL_CALCULATED=0
          ELSE
            CALL FLAG_ERROR("Constraint mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          !Initialise the RHS element vector
          RHS_MAPPING=>CONSTRAINT_MAPPING%RHS_MAPPING
          IF(ASSOCIATED(RHS_MAPPING)) THEN
            ROW_FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_SETUP(RHS_VECTOR%ELEMENT_VECTOR,ROW_FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("RHS mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_ELEMENT_INITIALISE

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
      CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_FINALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT_MATRIX%TEMP_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(CONSTRAINT_MATRIX%TEMP_VECTOR,ERR,ERROR,*999)
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

  !>Initialise the dynamic constraint matrix.
  SUBROUTINE CONSTRAINT_MATRIX_DYNAMIC_INITIALISE(DYNAMIC_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES !<A pointer to the dynamic matrices to initialise the dynamic constraint matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The dynamic matrix number in the dynamic constraint matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRIX_DYNAMIC_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
        CONSTRAINT_MATRICES=>DYNAMIC_MATRICES%CONSTRAINT_MATRICES
        IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
          CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
          IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
            DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
            IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
              IF(ASSOCIATED(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
                LOCAL_ERROR="Constraint matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is already associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                ALLOCATE(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix.",ERR,ERROR,*999)
                CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
                CONSTRAINT_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
                CONSTRAINT_MATRIX%DYNAMIC_MATRICES=>DYNAMIC_MATRICES
                NULLIFY(CONSTRAINT_MATRIX%LINEAR_MATRICES)
                CONSTRAINT_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
                CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
                CONSTRAINT_MATRIX%UPDATE_MATRIX=.TRUE.
                CONSTRAINT_MATRIX%FIRST_ASSEMBLY=.TRUE.
                CONSTRAINT_MATRIX%HAS_TRANSPOSE=DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%HAS_TRANSPOSE
                CONSTRAINT_MATRIX%NUMBER_OF_ROWS=DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_ROWS
                CONSTRAINT_MATRIX%TOTAL_NUMBER_OF_ROWS=DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)% &
                & TOTAL_NUMBER_OF_ROWS
                DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%CONSTRAINT_MATRIX=>CONSTRAINT_MATRIX
                NULLIFY(CONSTRAINT_MATRIX%MATRIX)
                NULLIFY(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE)
                NULLIFY(CONSTRAINT_MATRIX%TEMP_VECTOR)
                NULLIFY(CONSTRAINT_MATRIX%TEMP_TRANSPOSE_VECTOR)
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint mapping dynamic mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dynamic matrices constraint matrices is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified dynamic matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Dynamic matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRIX_DYNAMIC_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRIX_FINALISE(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRIX_DYNAMIC_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRIX_DYNAMIC_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRIX_DYNAMIC_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialise the linear constraint matrix.
  SUBROUTINE CONSTRAINT_MATRIX_LINEAR_INITIALISE(LINEAR_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES !<A pointer to the linear matrices to initialise the linear constraint matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The linear matrix number in the linear constraint matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRIX_LINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(LINEAR_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
        CONSTRAINT_MATRICES=>LINEAR_MATRICES%CONSTRAINT_MATRICES
        IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
          CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
          IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
            LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
              IF(ASSOCIATED(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
                LOCAL_ERROR="Constraint matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is already associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                ALLOCATE(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix.",ERR,ERROR,*999)
                CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
                CONSTRAINT_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
                NULLIFY(CONSTRAINT_MATRIX%DYNAMIC_MATRICES)
                CONSTRAINT_MATRIX%LINEAR_MATRICES=>LINEAR_MATRICES
                CONSTRAINT_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
                CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
                CONSTRAINT_MATRIX%UPDATE_MATRIX=.TRUE.
                CONSTRAINT_MATRIX%FIRST_ASSEMBLY=.TRUE.
                CONSTRAINT_MATRIX%HAS_TRANSPOSE=LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%HAS_TRANSPOSE
                CONSTRAINT_MATRIX%NUMBER_OF_ROWS=LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_ROWS
                CONSTRAINT_MATRIX%TOTAL_NUMBER_OF_ROWS=LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%
                & TOTAL_NUMBER_OF_ROWS
                LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%CONSTRAINT_MATRIX=>CONSTRAINT_MATRIX
                NULLIFY(CONSTRAINT_MATRIX%MATRIX)
                NULLIFY(CONSTRAINT_MATRIX%MATRIX_TRANSPOSE)
                NULLIFY(CONSTRAINT_MATRIX%TEMP_VECTOR)
                NULLIFY(CONSTRAINT_MATRIX%TEMP_TRANSPOSE_VECTOR)
                CALL CONSTRAINT_MATRICES_ELEMENT_MATRIX_INITIALISE(CONSTRAINT_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint mapping linear mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Linear matrices constraint matrices is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified linear matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRIX_LINEAR_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRIX_FINALISE(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRIX_LINEAR_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRIX_LINEAR_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRIX_LINEAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the constraint matrices dynamic matrices and deallocates all memory
  SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_FINALISE(DYNAMIC_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES !<A pointer to the constraint matrices dynamic matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
     
    CALL ENTERS("CONSTRAINT_MATRICES_DYNAMIC_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
      IF(ALLOCATED(DYNAMIC_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(DYNAMIC_MATRICES%MATRICES,1)
          CALL CONSTRAINT_MATRIX_FINALISE(DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(DYNAMIC_MATRICES%MATRICES)
      ENDIF
      IF(ASSOCIATED(DYNAMIC_MATRICES%TEMP_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(DYNAMIC_MATRICES%TEMP_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(DYNAMIC_MATRICES)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_DYNAMIC_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the constraint matrices dynamic matrices
  SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
     
    CALL ENTERS("CONSTRAINT_MATRICES_DYNAMIC_INITIALISE",ERR,ERROR,*998)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
        DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
          IF(ASSOCIATED(CONSTRAINT_MATRICES%DYNAMIC_MATRICES)) THEN
            CALL FLAG_ERROR("Constraint matrices dynamic matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(CONSTRAINT_MATRICES%DYNAMIC_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices dynamic matrices.",ERR,ERROR,*999)
            CONSTRAINT_MATRICES%DYNAMIC_MATRICES%CONSTRAINT_MATRICES=>CONSTRAINT_MATRICES
            CONSTRAINT_MATRICES%DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
            ALLOCATE(CONSTRAINT_MATRICES%DYNAMIC_MATRICES%MATRICES(DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices dynamic matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              NULLIFY(CONSTRAINT_MATRICES%DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL CONSTRAINT_MATRIX_DYNAMIC_INITIALISE(CONSTRAINT_MATRICES%DYNAMIC_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
            NULLIFY(CONSTRAINT_MATRICES%DYNAMIC_MATRICES%TEMP_VECTOR)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices equations mapping is not associated.",ERR,ERROR,*998)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRICES_DYNAMIC_FINALISE(CONSTRAINT_MATRICES%DYNAMIC_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_DYNAMIC_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Adds the Jacobain matrices into the constraint Jacobian.
  SUBROUTINE CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobian_matrix_idx
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD()")
#endif

    CALL ENTERS("CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        DO jacobian_matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(jacobian_matrix_idx)%PTR
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
              !Add in Jacobian element matrices
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(JACOBIAN_MATRIX%JACOBIAN,JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS,1:JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS), &
                & ERR,ERROR,*999)
            ENDIF
            IF(JACOBIAN_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(JACOBIAN_MATRIX%JACOBIAN_TRANSPOSE,JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS),TRANSPOSE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS,1:JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS)), &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Jacobian matrix for Jacobian matrix index "// &
              & TRIM(NUMBER_TO_VSTRING(jacobian_matrix_idx,"*",ERR,ERROR))//" is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !jacobian_matrix_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD()")
#endif
    
    CALL EXITS("CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Outputs the constraint Jacobian matrices
  SUBROUTINE CONSTRAINT_MATRICES_JACOBIAN_OUTPUT(ID,CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint Jacobian matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobian_matrix_idx
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_JACOBIAN_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"",ERR,ERROR,*999)
          CALL WRITE_STRING(ID,"Jacobian matrices:",ERR,ERROR,*999)
          DO jacobian_matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(jacobian_matrix_idx)%PTR
            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN        
              CALL WRITE_STRING(ID,"Jacobian matrix:",ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="Jacobian matrix for Jacobian matrix index "// &
                & TRIM(NUMBER_TO_VSTRING(jacobian_matrix_idx,"*",ERR,ERROR))//" is not associated."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !jacobian_matrix_idx
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_JACOBIAN_OUTPUT")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_JACOBIAN_OUTPUT",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_JACOBIAN_OUTPUT")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_JACOBIAN_OUTPUT
  
  !
  !================================================================================================================================
  !

  !>Sets the Jacobian calculation types of the residual variables
  SUBROUTINE ConstraintMatrices_JacobianTypesSet(constraintMatrices,jacobianTypes,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: constraintMatrices !<A pointer to the constraint matrices to set the Jacobian type for.
    INTEGER(INTG), INTENT(IN) :: jacobianTypes(:) !<jacobianTypes(matrix_idx). The Jacobian calculation type for the matrix_idx'th Jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    INTEGER(INTG) :: matrixIdx,numberOfjacobians,jacobianType
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("ConstraintMatrices_JacobianTypesSet",err,error,*999)

    IF(ASSOCIATED(constraintMatrices)) THEN
      IF(constraintMatrices%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",err,error,*999)
      ELSE
        nonlinearMatrices=>constraintMatrices%NONLINEAR_MATRICES
        IF(ASSOCIATED(nonlinearMatrices)) THEN
          numberOfJacobians=SIZE(jacobianTypes,1)
          IF(numberOfJacobians==nonlinearMatrices%NUMBER_OF_JACOBIANS) THEN
            DO matrixIdx=1,numberOfJacobians
              jacobianType=jacobianTypes(matrixIdx)
              SELECT CASE(jacobianType)
              CASE(CONSTRAINT_JACOBIAN_FINITE_DIFFERENCE_CALCULATED, &
                  & CONSTRAINT_JACOBIAN_ANALYTIC_CALCULATED)
                nonlinearMatrices%JACOBIANS(matrixIdx)%PTR%JACOBIAN_CALCULATION_TYPE=jacobianType
              CASE DEFAULT
                localError="Invalid Jacobian calculation type of " &
                  & //TRIM(NUMBER_TO_VSTRING(jacobianType,"*",err,error))//"."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            END DO
          ELSE
            localError="Invalid number of Jacobian calculation types. The number of types " &
              & //TRIM(NUMBER_TO_VSTRING(numberOfJacobians,"*",err,error)) &
              & //" should be "//TRIM(NUMBER_TO_VSTRING(nonlinearMatrices%NUMBER_OF_JACOBIANS,"*",err,error))
            CALL FLAG_ERROR(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices nonlinear matrices are not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices are not associated",err,error,*999)
    ENDIF

    CALL EXITS("ConstraintMatrices_JacobianTypesSet")
    RETURN
999 CALL ERRORS("ConstraintMatrices_JacobianTypesSet",err,error)
    CALL EXITS("ConstraintMatrices_JacobianTypesSet")
    RETURN 1
  END SUBROUTINE ConstraintMatrices_JacobianTypesSet

  !
  !================================================================================================================================
  !

  !>Finalises the constraint matrices linear matrices and deallocates all memory
  SUBROUTINE CONSTRAINT_MATRICES_LINEAR_FINALISE(LINEAR_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES !<A pointer to the constraint matrices linear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
     
    CALL ENTERS("CONSTRAINT_MATRICES_LINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_MATRICES)) THEN
      IF(ALLOCATED(LINEAR_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(LINEAR_MATRICES%MATRICES,1)
          CALL CONSTRAINT_MATRIX_FINALISE(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(LINEAR_MATRICES%MATRICES)
      ENDIF
      DEALLOCATE(LINEAR_MATRICES)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_LINEAR_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_LINEAR_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the constraint matrices linear matrices
  SUBROUTINE CONSTRAINT_MATRICES_LINEAR_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the linear matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
     
    CALL ENTERS("CONSTRAINT_MATRICES_LINEAR_INITIALISE",ERR,ERROR,*998)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
        LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
        IF(ASSOCIATED(LINEAR_MAPPING)) THEN
          IF(ASSOCIATED(CONSTRAINT_MATRICES%LINEAR_MATRICES)) THEN
            CALL FLAG_ERROR("Constraint matrices linear matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(CONSTRAINT_MATRICES%LINEAR_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices linear matrices.",ERR,ERROR,*999)
            CONSTRAINT_MATRICES%LINEAR_MATRICES%CONSTRAINT_MATRICES=>CONSTRAINT_MATRICES
            CONSTRAINT_MATRICES%LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES=LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
            ALLOCATE(CONSTRAINT_MATRICES%LINEAR_MATRICES%MATRICES(LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices linear matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
              NULLIFY(CONSTRAINT_MATRICES%LINEAR_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL CONSTRAINT_MATRIX_LINEAR_INITIALISE(CONSTRAINT_MATRICES%LINEAR_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices equations mapping is not associated.",ERR,ERROR,*998)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRICES_LINEAR_FINALISE(CONSTRAINT_MATRICES%LINEAR_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_LINEAR_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_LINEAR_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the constraint matrices nonlinear matrices and deallocates all memory
  SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_FINALISE(NONLINEAR_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer to the constraint matrices nonlinear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
     
    CALL ENTERS("CONSTRAINT_MATRICES_NONLINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
      IF(ALLOCATED(NONLINEAR_MATRICES%JACOBIANS)) THEN
        DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          CALL CONSTRAINT_JACOBIAN_FINALISE(NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO
        DEALLOCATE(NONLINEAR_MATRICES%JACOBIANS)
      ENDIF
      IF(ASSOCIATED(NONLINEAR_MATRICES%RESIDUAL)) CALL DISTRIBUTED_VECTOR_DESTROY(NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
      CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ERR,ERROR,*999)
      DEALLOCATE(NONLINEAR_MATRICES)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_NONLINEAR_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the constraint matrices nonlinear matrices
  SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the nonlinear matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_NONLINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CONSTRAINT_MAPPING=>CONSTRAINT_MATRICES%CONSTRAINT_MAPPING
      IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
        NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
        IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
          IF(ASSOCIATED(CONSTRAINT_MATRICES%NONLINEAR_MATRICES)) THEN
            CALL FLAG_ERROR("Constraint matrices nonlinear matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(CONSTRAINT_MATRICES%NONLINEAR_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices nonlinear matrices.",ERR,ERROR,*999)
            CONSTRAINT_MATRICES%NONLINEAR_MATRICES%CONSTRAINT_MATRICES=>CONSTRAINT_MATRICES
            CONSTRAINT_MATRICES%NONLINEAR_MATRICES%UPDATE_RESIDUAL=.TRUE.
            CONSTRAINT_MATRICES%NONLINEAR_MATRICES%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(CONSTRAINT_MATRICES%NONLINEAR_MATRICES%RESIDUAL)
            CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE(CONSTRAINT_MATRICES%NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ERR,ERROR,*999)
            CONSTRAINT_MATRICES%NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS=NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
            ALLOCATE(CONSTRAINT_MATRICES%NONLINEAR_MATRICES%JACOBIANS(NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrices Jacobian matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
              NULLIFY(CONSTRAINT_MATRICES%NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR)
              CALL CONSTRAINT_JACOBIAN_INITIALISE(CONSTRAINT_MATRICES%NONLINEAR_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint matrices equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MATRICES_NONLINEAR_FINALISE(CONSTRAINT_MATRICES%NONLINEAR_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MATRICES_NONLINEAR_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_INITIALISE
  
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
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    
    CALL ENTERS("CONSTRAINT_MATRICES_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"",ERR,ERROR,*999)
        CALL WRITE_STRING(ID,"Constraint matrices:",ERR,ERROR,*999)
        DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Dynamic matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of dynamic matrices = ",DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,ERR,ERROR,*999)
          DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
            CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
              CALL WRITE_STRING_VALUE(ID,"Constraint matrix : ",matrix_idx,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
        LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Linear matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of linear matrices = ",LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,ERR,ERROR,*999)
          DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
            CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
              CALL WRITE_STRING_VALUE(ID,"Constraint matrix : ",matrix_idx,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,CONSTRAINT_MATRIX%MATRIX,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Nonlinear vectors:",ERR,ERROR,*999)
          IF(ASSOCIATED(NONLINEAR_MATRICES%RESIDUAL)) THEN
            CALL WRITE_STRING(ID,"Residual vector:",ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Nonlinear matrices residual is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          CALL WRITE_STRING(ID,"RHS vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,RHS_VECTOR%VECTOR,ERR,ERROR,*999)
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
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR !<A pointer to the constraint matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("CONSTRAINT_MATRICES_RHS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_VECTOR)) THEN
      IF(ASSOCIATED(RHS_VECTOR%VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(RHS_VECTOR%VECTOR,ERR,ERROR,*999)
      CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_FINALISE(RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
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
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the rhs vector for
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
            NULLIFY(CONSTRAINT_MATRICES%RHS_VECTOR%VECTOR)
            CALL CONSTRAINT_MATRICES_ELEMENT_VECTOR_INITIALISE(CONSTRAINT_MATRICES%RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
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

  !>Sets the storage type (sparsity) of the dynamic constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET(CONSTRAINT_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th dynamic constraint matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have already been finished.",ERR,ERROR,*999)
      ELSE
        DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          IF(SIZE(STORAGE_TYPE,1)==DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
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
                    & " for the dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the linear constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET(CONSTRAINT_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th linear constraint matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          IF(SIZE(STORAGE_TYPE,1)==LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
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
                    & " for the linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of linear matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices linear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_LINEAR_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0(CONSTRAINT_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th Jacobian constraint matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          IF(SIZE(STORAGE_TYPE,1)==NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                SELECT CASE(STORAGE_TYPE(matrix_idx))
                CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
                CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                  JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
                CASE DEFAULT
                  LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the Jacobian matrix is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of Jacobian matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of all nonlinear (Jacobian) constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1(CONSTRAINT_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE !<STORAGE_TYPE. The storage type for all Jacobian constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: STORAGE_TYPES(:)
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES

    CALL ENTERS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          ALLOCATE(STORAGE_TYPES(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate storage types.",ERR,ERROR,*999)
          STORAGE_TYPES=STORAGE_TYPE
          CALL CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_0(CONSTRAINT_MATRICES,STORAGE_TYPES,ERR,ERROR,*999)
          DEALLOCATE(STORAGE_TYPES)
        ELSE
          CALL FLAG_ERROR("Constraint matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STORAGE_TYPE_SET_1

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the dynamic constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(CONSTRAINT_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The storage type for the matrix_idx'th dynamic constraint matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
           DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                SELECT CASE(STRUCTURE_TYPE(matrix_idx))
                CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
                  CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
                CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
                  CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_FEM_STRUCTURE
                CASE(CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE)
                  CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for dynamic matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the linear constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET(CONSTRAINT_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The storage type for the matrix_idx'th linear constraint matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
                SELECT CASE(STRUCTURE_TYPE(matrix_idx))
                CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
                  CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
                CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
                  CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_FEM_STRUCTURE
                CASE(CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE)
                  CONSTRAINT_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for linear matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of linear matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices linear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_LINEAR_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0(CONSTRAINT_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The structure type for the matrix_idx'th Jacobian constraint matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                SELECT CASE(STRUCTURE_TYPE(matrix_idx))
                CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_NO_STRUCTURE
                CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_FEM_STRUCTURE
                CASE(CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for the Jacobian matrix is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of Jacobian matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of all nonlinear (Jacobian) constraint matrices
  SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1(CONSTRAINT_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE !<The structure type for all Jacobian constraint matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: STRUCTURE_TYPES(:)
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES

    CALL ENTERS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Constraint matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          ALLOCATE(STRUCTURE_TYPES(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate storage types.",ERR,ERROR,*999)
          STRUCTURE_TYPES=STRUCTURE_TYPE
          CALL CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_0(CONSTRAINT_MATRICES,STRUCTURE_TYPES,ERR,ERROR,*999)
          DEALLOCATE(STRUCTURE_TYPES)
        ELSE
          CALL FLAG_ERROR("Constraint matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1")
    RETURN
999 CALL ERRORS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET_1

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
   
    CALL ENTERS("CONSTRAINT_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      CALL CONSTRAINT_MATRICES_DYNAMIC_FINALISE(CONSTRAINT_MATRICES%DYNAMIC_MATRICES,ERR,ERROR,*999)
      CALL CONSTRAINT_MATRICES_LINEAR_FINALISE(CONSTRAINT_MATRICES%LINEAR_MATRICES,ERR,ERROR,*999)
      CALL CONSTRAINT_MATRICES_NONLINEAR_FINALISE(CONSTRAINT_MATRICES%NONLINEAR_MATRICES,ERR,ERROR,*999)
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

  !>Initialise the constraint matrices for the equations.
  SUBROUTINE CONSTRAINT_MATRICES_INITIALISE(EQUATIONS,ERR,ERROR,*)
    
     !Argument variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<A pointer to the equations to initialise the constraint matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("CONSTRAINT_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
      IF(ASSOCIATED(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES)) THEN
        CALL FLAG_ERROR("Constraint matrices is already associated for this equations.",ERR,ERROR,*998)
      ELSE
        CONSTRAINT_MAPPING=>CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
            ALLOCATE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations constraint matrices.",ERR,ERROR,*999)
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%CONSTRAINT_EQUATIONSEQUATIONS
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%CONSTRAINT_MATRICES_FINISHED=.FALSE.
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%CONSTRAINT_MAPPING=>CONSTRAINT_MAPPING
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%SOLVER_MAPPING)
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NUMBER_OF_ROWS=CONSTRAINT_MAPPING%NUMBER_OF_ROWS
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%TOTAL_NUMBER_OF_ROWS=CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_ROWS
            CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NUMBER_OF_GLOBAL_ROWS=CONSTRAINT_MAPPING%NUMBER_OF_GLOBAL_ROWS
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%DYNAMIC_MATRICES)
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%LINEAR_MATRICES)
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%NONLINEAR_MATRICES)
            NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES%RHS_VECTOR)
            CALL CONSTRAINT_MATRICES_DYNAMIC_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,ERR,ERROR,*999)            
            CALL CONSTRAINT_MATRICES_LINEAR_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,ERR,ERROR,*999)            
            CALL CONSTRAINT_MATRICES_NONLINEAR_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,ERR,ERROR,*999)            
            CALL CONSTRAINT_MATRICES_RHS_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES,ERR,ERROR,*999)            
          ELSE
            CALL FLAG_ERROR("Constraint mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
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

  !>Initialise the values of the constraint matrices and vectors to the given value e.g., 0.0_DP
  SUBROUTINE CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,SELECTION_TYPE,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES !<A pointer to the constraint matrices to initialise the values for
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The selection of constraint matrices to be initialised \see CONSTRAINT_MATRICES_ROUTINES::SelectMatricesTypes,CONSTRAINT_MATRICES_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: CONSTRAINT_MATRIX
    
    CALL ENTERS("CONSTRAINT_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
      IF(SELECTION_TYPE==CONSTRAINT_MATRICES_ALL.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_NONLINEAR_ONLY) THEN
        DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
            CONSTRAINT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
              IF(CONSTRAINT_MATRIX%UPDATE_MATRIX) THEN
                CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(CONSTRAINT_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==CONSTRAINT_MATRICES_ALL.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_NONLINEAR_ONLY) THEN
        LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
            CONSTRAINT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(CONSTRAINT_MATRIX)) THEN
              IF(CONSTRAINT_MATRIX%UPDATE_MATRIX) THEN
                CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(CONSTRAINT_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==CONSTRAINT_MATRICES_ALL.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_JACOBIAN_ONLY) THEN
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
              IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
                CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(JACOBIAN_MATRIX%JACOBIAN,VALUE,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==CONSTRAINT_MATRICES_ALL.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_RESIDUAL_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_RHS_RESIDUAL_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_VECTORS_ONLY) THEN
        NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
            CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(NONLINEAR_MATRICES%RESIDUAL,VALUE,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==CONSTRAINT_MATRICES_ALL.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_RHS_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_RHS_RESIDUAL_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_VECTORS_ONLY) THEN    
        RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          IF(RHS_VECTOR%UPDATE_VECTOR) THEN
            CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(RHS_VECTOR%VECTOR,VALUE,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==CONSTRAINT_MATRICES_ALL.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==CONSTRAINT_MATRICES_VECTORS_ONLY) THEN    
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

  !>Caclulates the matrix structure (sparsity) for a constraint matrix.
  SUBROUTINE ConstraintMatrix_StructureCalculate(constraintMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRIX_TYPE), POINTER :: constraintMatrix !<A pointer to the constraint matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    type(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) ::  columnIdx,dummyErr,elementIdx,globalColumn,localColumn,local_ny,matrixNumber,mk,mp,ne,nh,nh2,nn,nnk,np
    INTEGER(INTG) ::  numberOfColumns,nyy,nyyg,npg,nhg,local_cols,local_dof,mv
    INTEGER(INTG) :: dofIdx,nodeIdx,componentIdx,localDofIdx
    INTEGER(INTG) :: versionIdx,derivativeIdx,numberOfDerivatives,numberOfVersions
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: dependentDofsDomainMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: equations
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: constraintMapping
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: constraintMatrices
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: constraintCondition
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: dependentDofsParamMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError
    CALL ENTERS("ConstraintMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    IF(ASSOCIATED(constraintMatrix)) THEN
      IF(.NOT.ASSOCIATED(rowIndices)) THEN
        IF(.NOT.ASSOCIATED(columnIndices)) THEN
          matrixNumber=constraintMatrix%MATRIX_NUMBER
          SELECT CASE(constraintMatrix%STRUCTURE_TYPE)
          CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
            CALL FLAG_ERROR("There is no structure to calculate for a matrix with no structure.",err,error,*998)
          CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
            SELECT CASE(constraintMatrix%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              linearMatrices=>constraintMatrix%LINEAR_MATRICES
              dynamicMatrices=>constraintMatrix%DYNAMIC_MATRICES
              IF(ASSOCIATED(dynamicMatrices).OR.ASSOCIATED(linearMatrices)) THEN
                IF(ASSOCIATED(dynamicMatrices)) THEN
                  constraintMatrices=>dynamicMatrices%CONSTRAINT_MATRICES
                ELSE
                  constraintMatrices=>linearMatrices%CONSTRAINT_MATRICES
                ENDIF
                IF(ASSOCIATED(constraintMatrices)) THEN
                  constraintEquations=>constraintMatrices%EQUATIONS
                  IF(ASSOCIATED(constraintEquations)) THEN
                    constraintMapping=>constraintMatrices%CONSTRAINT_MAPPING
                    IF(ASSOCIATED(constraintMapping)) THEN
                      dynamicMapping=>constraintMapping%DYNAMIC_MAPPING
                      linearMapping=>constraintMapping%LINEAR_MAPPING
                      IF(ASSOCIATED(dynamicMapping).OR.ASSOCIATED(linearMapping)) THEN
                        constraintCondition=>constraintEquations%CONSTRAINT_CONDITION
                        IF(ASSOCIATED(constraintCondition)) THEN
                          IF(ASSOCIATED(dynamicMatrices)) THEN
                            rowVariable=>dynamicMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrixNumber)%VARIABLE
                            columnvariable=>dynamicMapping%LAGRANGE_VARIABLE
                          ELSE
                            rowVariable=>linearMapping%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrixNumber)%VARIABLE
                            columnvariable=>linearMapping%LAGRANGE_VARIABLE
                          ENDIF
                          IF(ASSOCIATED(rowVariable)) THEN
                            rowDofsDomainMapping=>rowVariable%DOMAIN_MAPPING
                            IF(ASSOCIATED(columnVariable)) THEN
                              columnDofsDomainMapping=>columnVariable%DOMAIN_MAPPING
                              IF(ASSOCIATED(rowDofsDomainMapping)) THEN
                                rowDofsParamMapping=>rowVariable%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(columnDofsDomainMapping)) THEN
                                  columnDofsParamMapping=>columnVariable%DOF_TO_PARAM_MAP
                                  IF(ASSOCIATED(rowDofsParamMapping)) THEN
                                    !Allocate lists
                                    ALLOCATE(columnIndicesLists(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",err,error,*999)
                                    DO local_row=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                      !Set up list
                                      NULLIFY(columnIndicesLists(local_row)%PTR)
                                      CALL LIST_CREATE_START(columnIndicesLists(local_row)%PTR,ERR,ERROR,*999)
                                      CALL LIST_DATA_TYPE_SET(columnIndicesLists(local_row)%PTR,LIST_INTG_TYPE, &
                                        & ERR,ERROR,*999)
                                      CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(local_row)%PTR,50,ERR,ERROR,*999)
                                      CALL LIST_CREATE_FINISH(columnIndicesLists(local_row)%PTR,ERR,ERROR,*999)
                                    ENDDO !local_row
                                    !Allocate row indices
                                    ALLOCATE(rowIndices(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",ERR,ERROR,*999)
                                    IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                      !Allocate transpose lists
                                      ALLOCATE(transposeColumnIndicesLists(columnDofsDomainMapping% &
                                        & TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate transpose column indices lists.", &
                                        & ERR,ERROR,*999)
                                      DO local_column=1,columnDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        !Set up list
                                        NULLIFY(transposeColumnIndicesLists(local_column)%PTR)
                                        CALL LIST_CREATE_START(transposeColumnIndicesLists(local_column)%PTR, &
                                          & ERR,ERROR,*999)
                                        CALL LIST_DATA_TYPE_SET(transposeColumnIndicesLists(local_column)%PTR, &
                                          & LIST_INTG_TYPE,ERR,ERROR,*999)
                                        CALL LIST_INITIAL_SIZE_SET(transposeColumnIndicesLists(local_column)%PTR,50, &
                                          & ERR,ERROR,*999)
                                        CALL LIST_CREATE_FINISH(transposeColumnIndicesLists(local_column)%PTR, &
                                          & ERR,ERROR,*999)
                                      ENDDO !local_column
                                      !Allocate transpose row indices
                                      ALLOCATE(transposeRowIndices(columnDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1), &
                                        & STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate transpose row indices.",ERR,ERROR,*999)
                                    ENDIF
                                    !First, loop over the Lagrange multiplier column dofs and calculate the number of non-zeros
                                    numberOfNonZeros=0
                                    DO local_column=1,columnDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                      IF(columnDofsParamMapping%DOF_TYPE(1,local_column)==FIELD_NODE_DOF_TYPE) THEN
                                        column_field_dof_idx=columnDofsParamMapping%DOF_TYPE(2,local_column)!value for a particular field dof (local_ny)
                                        column_node_idx=columnDofsParamMapping%NODE_DOF2PARAM_MAP(3,column_field_dof_idx)
                                        column_component_idx=columnDofsParamMapping%NODE_DOF2PARAM_MAP(4,column_field_dof_idx)
                                        domainNodes=>columnVariable%COMPONENTS(column_component_idx)%DOMAIN%TOPOLOGY%NODES
                                        
                                        !Loop over all elements containing the dof
                                        DO element_idx=1,domainNodes%NODES(column_node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                                          row_local_element_idx=domainNodes%NODES(column_node_idx)%SURROUNDING_ELEMENTS(elementIdx)
                                          DO n=1,fieldVariable%NUMBER_OF_COMPONENTS
                                            domainElements=>fieldVariable%COMPONENTS(nh2)%DOMAIN%TOPOLOGY%ELEMENTS
                                            basis=>domainElements%ELEMENTS(ne)%BASIS
                                            DO nn=1,basis%NUMBER_OF_NODES
                                              mp=domainElements%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                              DO nnk=1,basis%NUMBER_OF_DERIVATIVES(nn)
                                                mk=domainElements%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                                mv=domainElements%ELEMENTS(ne)%elementVersions(nnk,nn)
                                                !Find the local and global column and add the global column to the indices list
                                                localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                                  & NODES(mp)%DERIVATIVES(mk)%VERSIONS(mv)
                                                globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                            
                                                CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn,err,error,*999)
                                                  
                                              ENDDO !mk
                                            ENDDO !nn
                                          ENDDO !nh2
                                        ENDDO !elementIdx
                                        CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                        CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_ny)%PTR,numberOfColumns, &
                                          & err,error,*999)
                                        numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                        rowIndices(local_ny+1)=numberOfNonZeros+1
                                      ELSE
                                        localError="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",err,error))// &
                                          & " is not a node based dof."
                                        CALL FLAG_ERROR(localError,err,error,*999)
                                      ENDIF
                                    ENDDO !local_ny
                                    !Loop over the number of components in the Lagrange multipler variable
                                    DO column_component_idx=1,columnVariable%NUMBER_OF_COMPONENTS
                                      IF(columnVariable%COMPONENTS(column_component_idx)%INTERPOLATION_TYPE== &
                                        & FIELD_NODE_BASED_INTERPOLATION) THEN
                                        !Loop over the elements in the constraint mesh
                                        columnDomainElements=>columnVariable%COMPONENTS(column_component_idx)%DOMAIN% &
                                          & TOPOLOGY%ELEMENTS
                                        DO element_idx=1,columnDomainElements%TOTAL_NUMBER_OF_ELEMENTS
                                          columnBasis=>columnDomainElements%ELEMENTS(element_idx)%BASIS
                                          !Loop over the column DOFs in the element
                                          DO column_local_node_idx=1,columnBasis%NUMBER_OF_NODES
                                            column_node=columnDomainElements%ELEMENTS(element_idx)% &
                                              & ELEMENT_NODES(column_local_node_idx)
                                            DO column_local_derivative_idx=1,columnBasis% &
                                                & NUMBER_OF_DERIVATIVES(column_local_node_idx)
                                              column_derivative=columnDomainElements%ELEMENTS(element_idx)% &
                                                & ELEMENT_DERIVATIVES(column_local_derivative_idx,column_local_node_idx)
                                              column_version=columnDomainElements%ELEMENTS(element_idx)% &
                                                & elementVersions(column_local_derivative_idx,column_local_node_idx)
                                              local_column=columnVariable%COMPONENTS(column_component_idx)% &
                                                & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(column_node)% &
                                                & DERIVATIVES(column_derivative)%VERSIONS(column_version)
                                              global_column=columnDofsDomainMapping%LOCAL_TO_GLOBAL_MAP(local_column)
                                              !Loop over the components in the dependent variable
                                              DO row_component_idx=1,rowVariable%NUMBER_OF_COMPONENTS
                                                SELECT CASE(rowVariable%COMPONENTS(row_component_idx)%INTERPOLATION_TYPE)
                                                CASE(FIELD_CONSTANT_INTERPOLATION)
                                                  local_row=rowVariable%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                    & CONSTANT_PARAM2DOF_MAP
                                                  CALL LIST_ITEM_ADD(columnIndicesLists(local_row)%PTR,global_column, &
                                                    & ERR,ERROR,*999)
                                                  IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                                    global_row=rowVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                    CALL LIST_ITEM_ADD(transposeColumnIndicesLists(local_column)%PTR, &
                                                      & global_row,ERR,ERROR,*999)
                                                  ENDIF
                                                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                                  domain_element=MESH_CONNECTIVITY% &
                                                    & ELEMENT_CONNECTIVITY(element_idx,CONSTRAINT_MESH_INDEX)% &
                                                    & COUPLED_MESH_ELEMENT_NUMBER
                                                  domain_element=rowVariable%COMPONENTS(row_component_idx)%DOMAIN%TOPOLOGY%ELEMENTS

                                                  local_row=rowVariable%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                    & ELEMENT_PARAM2DOF_MAP%ELEMENTS(domain_element)
                                                  CALL LIST_ITEM_ADD(columnIndicesLists(local_row)%PTR,global_column, &
                                                    & ERR,ERROR,*999)
                                                  IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                                    global_row=rowVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                    CALL LIST_ITEM_ADD(transposeColumnIndicesLists(local_column)%PTR, &
                                                      & global_row,ERR,ERROR,*999)
                                                  ENDIF
                                                CASE(FIELD_NODE_BASED_INTERPOLATION)
                                                  rowDomainElements=>rowVariable%COMPONENTS(row_component_idx)%DOMAIN% &
                                                    & TOPOLOGY%ELEMENTS
                                                  domain_element=MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY( &
                                                    & element_idx,CONSTRAINT_MESH_INDEX)%COUPLED_MESH_ELEMENT_NUMBER
                                                  rowBasis=>rowDomainElements%ELEMENTS(domain_element)%BASIS
                                                  !Loop over the row DOFs in the domain mesh element
                                                  DO row_local_node_idx=1,rowBasis%NUMBER_OF_NODES
                                                    row_node=rowDomainElements%ELEMENTS(domain_element)% &
                                                      & ELEMENT_NODES(row_local_node_idx)
                                                    DO row_local_derivative_idx=1,rowBasis% &
                                                        & NUMBER_OF_DERIVATIVES(row_local_node_idx)
                                                      row_derivative=rowDomainElements%ELEMENTS(domain_element)% &
                                                        & ELEMENT_DERIVATIVES(row_local_derivative_idx,row_local_node_idx)
                                                      row_version=rowDomainElements%ELEMENTS(domain_element)% &
                                                        & elementVersions(row_local_derivative_idx,row_local_node_idx)
                                                      local_row=rowVariable%COMPONENTS(row_component_idx)% &
                                                        & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(row_node)% &
                                                        & DERIVATIVES(row_derivative)%VERSIONS(row_version)
                                                      CALL LIST_ITEM_ADD(columnIndicesLists(local_row)%PTR,global_column, &
                                                        & ERR,ERROR,*999)
                                                      IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                                        global_row=rowVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                        CALL LIST_ITEM_ADD(transposeColumnIndicesLists(local_column)%PTR, &
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
                                                    & TRIM(NUMBER_TO_VSTRING(rowVariable%COMPONENTS(row_component_idx)% &
                                                    INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid."
                                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                                END SELECT
                                              ENDDO !row_component_idx
                                            ENDDO !column_local_derivative_idx
                                          ENDDO !column_local_node_idx
                                        ENDDO !element_idx
                                      ELSE
                                        CALL FLAG_ERROR("Only node based fields implemented.",ERR,ERROR,*999)
                                      ENDIF
                                    ENDDO !column_component_idx
                                    rowIndices(1)=1
                                    DO local_row=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                      CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_row)%PTR,ERR,ERROR,*999)
                                      CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                        & ERR,ERROR,*999)
                                      NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                      rowIndices(local_row+1)=NUMBER_OF_NON_ZEROS+1
                                    ENDDO !local_row
                                    IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                      transposeNUMBER_OF_NON_ZEROS=0
                                      transposeRowIndices(1)=1
                                      DO local_column=1,columnDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        CALL LIST_REMOVE_DUPLICATES(transposeColumnIndicesLists(local_column)%PTR, &
                                          & ERR,ERROR,*999)
                                        CALL LIST_NUMBER_OF_ITEMS_GET(transposeColumnIndicesLists(local_column)%PTR, &
                                          & NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                                        transposeNUMBER_OF_NON_ZEROS=transposeNUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                        transposeRowIndices(local_column+1)=transposeNUMBER_OF_NON_ZEROS+1
                                      ENDDO !local_column
                                      !Sanity check - the number of non-zeros should be the same
                                      IF(transposeNUMBER_OF_NON_ZEROS/=NUMBER_OF_NON_ZEROS) THEN
                                        LOCAL_ERROR="Invalid number of non-zeros. The number of non-zeros in the "// &
                                          & "transposed matrix ("//TRIM(NUMBER_TO_VSTRING(transposeNUMBER_OF_NON_ZEROS, &
                                          & "*",ERR,ERROR))//") does not match the number of non-zeros in the constraint "// &
                                          & "matrix ("//TRIM(NUMBER_TO_VSTRING(transposeNUMBER_OF_NON_ZEROS,"*",ERR, &
                                          & ERROR))//")."
                                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                      ENDIF
                                    ENDIF
                                    !Allocate and setup the column locations
                                    ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",ERR,ERROR,*999)
                                    DO local_row=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                      CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                        & COLUMNS,ERR,ERROR,*999)
                                      DO column_idx=1,NUMBER_OF_COLUMNS
                                        COLUMN_INDICES(rowIndices(local_row)+column_idx-1)=COLUMNS(column_idx)
                                      ENDDO !column_idx
                                      DEALLOCATE(COLUMNS)
                                    ENDDO !local_row
                                    IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN
                                      !Allocate and setup the column locations
                                      ALLOCATE(transposeColumnIndices(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                      IF(ERR/=0) &
                                        & CALL FLAG_ERROR("Could not allocate transpose column indices.",ERR,ERROR,*999)
                                      DO local_column=1,columnDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        CALL LIST_DETACH_AND_DESTROY(transposeColumnIndicesLists(local_column)%PTR, &
                                          & NUMBER_OF_ROWS,transposeColumns,ERR,ERROR,*999)
                                          DO row_idx=1,NUMBER_OF_ROWS
                                          transposeColumnIndices(transposeRowIndices(local_column)+row_idx-1)= &
                                            & transposeColumns(row_idx)
                                          ENDDO !row_idx
                                        DEALLOCATE(transposeColumns)
                                      ENDDO !local_column
                                    ENDIF
                                    IF(DIAGNOSTICS1) THEN
                                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Constraint matrix structure:",ERR,ERROR,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Constraint matrix number : ", &
                                        & MATRIX_NUMBER,ERR,ERROR,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                        & rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                        & columnDofsDomainMapping%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                        & NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                                      IF(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                        & columnDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                      SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(rowDofsDomainMapping% &
                                        & TOTAL_NUMBER_OF_LOCAL*columnDofsDomainMapping%NUMBER_OF_GLOBAL,DP)*100.0_DP
                                      CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY, &
                                        & "F6.2",ERR,ERROR,*999)
                                    ENDIF
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,rowDofsDomainMapping% &
                                      & TOTAL_NUMBER_OF_LOCAL+1,5,5,rowIndices, &
                                      & '("  Row indices              :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                      & COLUMN_INDICES,'("  Column indices           :",5(X,I13))','(28X,5(X,I13))', &
                                      & ERR,ERROR,*999)
                                    IF(CONSTRAINT_MATRIX%HAS_TRANSPOSE) THEN 
                                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,columnDofsDomainMapping% &
                                        & TOTAL_NUMBER_OF_LOCAL+1,5,5,transposeRowIndices, &
                                        & '("  Transpose row indices    :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
                                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                        & transposeColumnIndices,'("  Transpose column indices :",5(X,I13))', &
                                        & '(28X,5(X,I13))',ERR,ERROR,*999)
                                    ENDIF
                                  ENDIF


!BLABLABLABLA

                                    
                                    
                                    !Allocate and setup the column locations
                                    ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)

                                    ALLOCATE(list(dependentDofsDomainMapping%NUMBER_OF_GLOBAL))

                                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",err,error,*999)
                                    DO local_ny=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                     
                                      CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(local_ny)%PTR,numberOfColumns,columns, &
                                        & err,error,*999)        
                                      DO columnIdx=1,numberOfColumns
                                        !COLUMNS store the list of nonzero column indices for each local row (local_ny)
                                        columnIndices(rowIndices(local_ny)+columnIdx-1)=columns(columnIdx) 

                                        ! global to local columns
                                         IF(ASSOCIATED(linearMapping).OR.ASSOCIATED(dynamicMapping)) THEN
                                           IF(ASSOCIATED(dynamicMatrices)) THEN
                                             local_cols=constraintMatrices%equations_mapping%dynamic_mapping &
                                               & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                               & (columns(columnIdx))%LOCAL_NUMBER(1)
                                             local_dof = local_cols
                                             ! Column to dof mapping?
                                             !local_dof=constraintMatrices%equations_mapping%dynamic_mapping% &
                                              ! & equations_matrix_to_var_maps(1)%column_to_dof_map(local_cols)
                                           ELSE
                                             local_cols=constraintMatrices%equations_mapping%linear_mapping &
                                               & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                               & (columns(columnIdx))%LOCAL_NUMBER(1)
                                             local_dof = local_cols
                                           ENDIF
                                         ENDIF
                                         nyyg=dependentDofsParamMapping%DOF_TYPE(2,local_dof)
                                         npg=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,nyyg)
                                         nhg=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,nyyg)
                                         domainNodes=>fieldVariable%COMPONENTS(nhg)%DOMAIN%TOPOLOGY%NODES
                              
                                        ! Check whether boundary node    
                                        IF(domainNodes%NODES(npg)%BOUNDARY_NODE)THEN
                                          CALL LinkedList_Add(list(columns(columnIdx)),local_ny)
                                        ENDIF
                                      
                                      ENDDO !columnIdx
                                      DEALLOCATE(columns)                                    
                                    ENDDO !local_ny

                                   
                                    IF(DIAGNOSTICS1) THEN
                                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Constraint matrix structure:",err,error,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Constraint matrix number : ",matrixNumber, &
                                        & err,error,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                        & dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                        & dependentDofsDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
                                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                        & numberOfNonZeros,err,error,*999)
                                      IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                        & dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                        sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
                                          & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
                                        CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ", &
                                          & sparsity,"F6.2",err,error,*999)
                                      ENDIF
                                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,dependentDofsDomainMapping% &
                                        & TOTAL_NUMBER_OF_LOCAL+1,8,8,rowIndices,'("  Row indices    :",8(X,I13))', &
                                        & '(18X,8(X,I13))',err,error,*999)
                                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
                                        & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Dependent dofs parameter mapping is not associated.",err,error,*999)
                                  ENDIF
                              ELSE
                                CALL FLAG_ERROR("Dependent dofs domain mapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Dependent field variable is not associated.",err,error,*999)
                            ENDIF
                        ELSE
                          CALL FLAG_ERROR("Constraint condition is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Either equations mapping dynamic mapping or linear mapping is not associated.", &
                          & err,error,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Constraint mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dynamic or linear matrices constraint matrices is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Either constraint matrix dynamic or linear matrices is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(constraintMatrix%STORAGE_TYPE,"*",err,error))//" is invalid."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE(CONSTRAINT_MATRIX_DIAGONAL_STRUCTURE)
            CALL FLAG_ERROR("There is not structure to calculate for a diagonal matrix.",err,error,*998)
          CASE DEFAULT
            localError="The matrix structure type of "// &
              & TRIM(NUMBER_TO_VSTRING(constraintMatrix%STRUCTURE_TYPE,"*",err,error))//" is invalid."
            CALL FLAG_ERROR(localError,err,error,*998)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Column indices is already associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indices is already associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrix is not associated.",err,error,*999)
    ENDIF
      
    CALL EXITS("ConstraintMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%PTR)) &
          & CALL LIST_DESTROY(columnIndicesLists(localDofIdx)%PTR,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 CALL ERRORS("ConstraintMatrix_StructureCalculate",err,error)
    CALL EXITS("ConstraintMatrix_StructureCalculate")
    RETURN 1
  END SUBROUTINE ConstraintMatrix_StructureCalculate

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a Jacobian matrix.
  SUBROUTINE JacobianMatrix_StructureCalculate(jacobianMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_JACOBIAN_TYPE), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) ::  columnIdx,dummyErr,elementIdx,globalColumn,localColumn,local_ny,mk,mp,ne,nh,nh2,nn,nnk,np,mv, &
      & numberOfColumns,nyy,matrixNumber
    INTEGER(INTG) :: dofIdx,nodeIdx,componentIdx,versionIdx,derivativeIdx,numberOfVersions,numberOfDerivatives
    INTEGER(INTG) :: localDofIdx
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: dependentDofsDomainMapping,rowDofsDomainMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: equations
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: constraintMapping
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: constraintMatrices
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: nonLinearMatrices
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: constraintCondition
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: dependentDofsParamMapping,rowDofsParamMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,rowVariable
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    CALL ENTERS("JacobianMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    IF(ASSOCIATED(jacobianMatrix)) THEN
      matrixNumber=jacobianMatrix%JACOBIAN_NUMBER
      IF(.NOT.ASSOCIATED(rowIndices)) THEN
        IF(.NOT.ASSOCIATED(columnIndices)) THEN
          SELECT CASE(jacobianMatrix%STRUCTURE_TYPE)
          CASE(CONSTRAINT_MATRIX_NO_STRUCTURE)
            CALL FLAG_ERROR("Not implemented.",err,error,*998)
          CASE(CONSTRAINT_MATRIX_FEM_STRUCTURE)
            SELECT CASE(jacobianMatrix%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              nonlinearMatrices=>jacobianMatrix%NONLINEAR_MATRICES
              IF(ASSOCIATED(nonlinearMatrices)) THEN
                constraintMatrices=>nonlinearMatrices%CONSTRAINT_MATRICES
                IF(ASSOCIATED(constraintMatrices)) THEN
                  constraintEquations=>constraintMatrices%EQUATIONS
                  IF(ASSOCIATED(constraintEquations)) THEN
                    constraintMapping=>constraintMatrices%CONSTRAINT_MAPPING
                    IF(ASSOCIATED(constraintMapping)) THEN
                      nonlinearMapping=>constraintMapping%NONLINEAR_MAPPING
                      IF(ASSOCIATED(nonlinearMapping)) THEN
                        constraintCondition=>constraintEquations%CONSTRAINT_CONDITION
                        IF(ASSOCIATED(constraintCondition)) THEN
                          dependentField=>constraintCondition%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixNumber)%VARIABLE
                            IF(ASSOCIATED(fieldVariable)) THEN
                              dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
                              IF(ASSOCIATED(dependentDofsDomainMapping)) THEN
                                dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(dependentDofsParamMapping)) THEN
                                  !If RHS variable exists, use this for row DOFs, else use the first nonlinear variable
                                  IF(ASSOCIATED(constraintMapping%RHS_MAPPING)) THEN
                                    rowVariable=>constraintMapping%RHS_MAPPING%RHS_VARIABLE
                                  ELSE
                                    rowVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
                                  ENDIF
                                  IF(ASSOCIATED(rowVariable)) THEN
                                    rowDofsDomainMapping=>rowVariable%DOMAIN_MAPPING
                                    rowDofsParamMapping=>rowVariable%DOF_TO_PARAM_MAP
                                  ELSE
                                    CALL FLAG_ERROR("RHS or first nonlinear variable is not associated",err,error,*999)
                                  ENDIF
                                  IF(ASSOCIATED(rowDofsDomainMapping)) THEN
                                    IF(ASSOCIATED(rowDofsParamMapping)) THEN
                                      !Allocate lists
                                      ALLOCATE(columnIndicesLists(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",err,error,*999)
                                      !Allocate row indices
                                      ALLOCATE(rowIndices(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",err,error,*999)
                                      rowIndices(1)=1
                                      !First, loop over the rows and calculate the number of non-zeros
                                      numberOfNonZeros=0
                                      DO local_ny=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        SELECT CASE(rowDofsParamMapping%DOF_TYPE(1,local_ny))
                                        CASE(FIELD_CONSTANT_INTERPOLATION)
                                          CALL FLAG_ERROR("Constant interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_NODE_DOF_TYPE)
                                          nyy=rowDofsParamMapping%DOF_TYPE(2,local_ny)
                                          np=rowDofsParamMapping%NODE_DOF2PARAM_MAP(3,nyy) !node number
                                          nh=rowDofsParamMapping%NODE_DOF2PARAM_MAP(4,nyy) !component number
                                          domainNodes=>rowVariable%COMPONENTS(nh)%DOMAIN%TOPOLOGY%NODES
                                          !Set up list
                                          NULLIFY(columnIndicesLists(local_ny)%PTR)
                                          CALL LIST_CREATE_START(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_DATA_TYPE_SET(columnIndicesLists(local_ny)%PTR,LIST_INTG_TYPE,err,error,*999)
                                          CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(local_ny)%PTR,domainNodes%NODES(np)% &
                                            & NUMBER_OF_SURROUNDING_ELEMENTS*rowVariable%COMPONENTS(nh)% &
                                            & maxNumberElementInterpolationParameters,err,error,*999)
                                          CALL LIST_CREATE_FINISH(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          !Loop over all elements containing the dof
                                          DO elementIdx=1,domainNodes%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                                            ne=domainNodes%NODES(np)%SURROUNDING_ELEMENTS(elementIdx)
                                            DO nh2=1,fieldVariable%NUMBER_OF_COMPONENTS
                                              SELECT CASE(fieldVariable%COMPONENTS(nh2)%INTERPOLATION_TYPE)
                                              CASE(FIELD_CONSTANT_INTERPOLATION)
                                                ! do nothing? this will probably never be encountered...?
                                              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                                localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                  & ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                                                globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                                CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn,err,error,*999)
                                              CASE(FIELD_NODE_BASED_INTERPOLATION)
                                                domainElements=>fieldVariable%COMPONENTS(nh2)%DOMAIN%TOPOLOGY%ELEMENTS
                                                basis=>domainElements%ELEMENTS(ne)%BASIS
                                                DO nn=1,basis%NUMBER_OF_NODES
                                                  mp=domainElements%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                                  DO nnk=1,basis%NUMBER_OF_DERIVATIVES(nn)
                                                    mk=domainElements%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                                    mv=domainElements%ELEMENTS(ne)%elementVersions(nnk,nn)
                                                    !Find the local and global column and add the global column to the indices list
                                                    localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                      & NODE_PARAM2DOF_MAP%NODES(mp)%DERIVATIVES(mk)%VERSIONS(mv)
                                                    globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                                    CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn, &
                                                      & err,error,*999)
                                                  ENDDO !mk
                                                ENDDO !nn
                                              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                                CALL FLAG_ERROR("Grid point based interpolation is not implemented yet.",& 
                                                  & err,error,*999)
                                              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                                CALL FLAG_ERROR("Gauss point based interpolation is not implemented yet.",&
                                                  & err,error,*999)
                                              CASE DEFAULT
                                                localError="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",err,error))// &
                                                  & " has invalid interpolation type."
                                                CALL FLAG_ERROR(localError,err,error,*999)
                                              END SELECT
                                            ENDDO !nh2
                                          ENDDO !elementIdx
                                          CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_ny)%PTR,numberOfColumns, &
                                            & err,error,*999)
                                          numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                          rowIndices(local_ny+1)=numberOfNonZeros+1
                                        CASE(FIELD_ELEMENT_DOF_TYPE)
                                          ! row corresponds to a variable that's element-wisely interpolated
                                          nyy=rowDofsParamMapping%DOF_TYPE(2,local_ny)          ! nyy = index in ELEMENT_DOF2PARAM_MAP
                                          ne=rowDofsParamMapping%ELEMENT_DOF2PARAM_MAP(1,nyy)   ! current element (i.e. corresponds to current dof)
                                          nh=rowDofsParamMapping%ELEMENT_DOF2PARAM_MAP(2,nyy)   ! current variable component
                                          domainElements=>rowVariable%COMPONENTS(nh)%DOMAIN%TOPOLOGY%ELEMENTS
                                          basis=>domainElements%ELEMENTS(ne)%BASIS
                                          !Set up list
                                          NULLIFY(columnIndicesLists(local_ny)%PTR)
                                          CALL LIST_CREATE_START(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_DATA_TYPE_SET(columnIndicesLists(local_ny)%PTR,LIST_INTG_TYPE,err,error,*999)
                                          CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(local_ny)%PTR, &
                                            & rowVariable%COMPONENTS(nh)%maxNumberElementInterpolationParameters+1, &
                                            & err,error,*999) ! size = all nodal dofs + itself
                                          CALL LIST_CREATE_FINISH(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          DO nh2=1,fieldVariable%NUMBER_OF_COMPONENTS
                                            SELECT CASE(fieldVariable%COMPONENTS(nh2)%INTERPOLATION_TYPE)
                                            CASE(FIELD_CONSTANT_INTERPOLATION)
                                              CALL FLAG_ERROR("Constant interpolation is not implemented yet.",err,error,*999)
                                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                              ! it's assumed that element-based variables arne't directly coupled
                                              ! put a diagonal entry
                                              localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                & ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                                              globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                              CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn,err,error,*999)
                                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                                              ! loop over all nodes in the element (and dofs belonging to them)
                                              DO nn=1,basis%NUMBER_OF_NODES
                                                mp=domainElements%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                                DO nnk=1,basis%NUMBER_OF_DERIVATIVES(nn)
                                                  mk=domainElements%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                                  mv=domainElements%ELEMENTS(ne)%elementVersions(nnk,nn)
                                                  !Find the local and global column and add the global column to the indices list
                                                  localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                                    & NODES(mp)%DERIVATIVES(mk)%VERSIONS(mv)
                                                  globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                                  CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn, &
                                                    & err,error,*999)
                                                ENDDO !mk
                                              ENDDO !nn
                                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                              CALL FLAG_ERROR("Grid point based interpolation is not implemented yet.", &
                                                & err,error,*999)
                                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                              CALL FLAG_ERROR("Gauss point based interpolation is not implemented yet.", &
                                                & err,error,*999)
                                            CASE DEFAULT
                                              localError="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",err,error))// &
                                                & " has invalid interpolation type."
                                              CALL FLAG_ERROR(localError,err,error,*999)
                                            END SELECT
                                          ENDDO !nh2
                                          ! clean up the list
                                          CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_ny)%PTR,numberOfColumns, &
                                            & err,error,*999)
                                          numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                          rowIndices(local_ny+1)=numberOfNonZeros+1
                                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                          CALL FLAG_ERROR("Grid point based interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                          CALL FLAG_ERROR("Gauss point based interpolation is not implemented yet.",err,error,*999)
                                        CASE DEFAULT
                                          localError="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",err,error))// &
                                            & " has an invalid type."
                                          CALL FLAG_ERROR(localError,err,error,*999)
                                        END SELECT
                                      ENDDO !local_ny
                                      !Allocate and setup the column locations
                                      ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",err,error,*999)
                                      DO local_ny=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(local_ny)%PTR,numberOfColumns,columns, &
                                          & err,error,*999)
                                        DO columnIdx=1,numberOfColumns
                                          columnIndices(rowIndices(local_ny)+columnIdx-1)=columns(columnIdx)
                                        ENDDO !columnIdx
                                        DEALLOCATE(columns)
                                      ENDDO !local_ny
                                      IF(DIAGNOSTICS1) THEN
                                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Jacobian matrix structure:",err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                          & dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                          & dependentDofsDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                          & numberOfNonZeros,err,error,*999)
                                        IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                          & dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                          sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
                                            & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
                                          CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ", &
                                            & sparsity,"F6.2",err,error,*999)
                                        ENDIF
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,rowDofsDomainMapping% &
                                          & TOTAL_NUMBER_OF_LOCAL+1,8,8,rowIndices,'("  Row indices    :",8(X,I13))', &
                                          & '(18X,8(X,I13))',err,error,*999)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices,&
                                          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Row dofs parameter mapping is not associated.",err,error,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Row dofs domain mapping is not associated.",err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Dependent dofs parameter mapping is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Dependent dofs domain mapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Dependent field variable is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Constraint condition dependent field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Constraint condition is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Constraint mapping nonlinear mapping is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Constraint mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Nonlinear matrices constraint matrices is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Constraint matrix nonlinear matrices is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(jacobianMatrix%STORAGE_TYPE,"*",err,error))//" is invalid."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The matrix structure type of "// &
              & TRIM(NUMBER_TO_VSTRING(jacobianMatrix%STRUCTURE_TYPE,"*",err,error))//" is invalid."
            CALL FLAG_ERROR(localError,err,error,*998)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Column indices is already associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indices is already associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Jacobian matrix is not associated.",err,error,*999)
    ENDIF
      
    CALL EXITS("JacobianMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%PTR)) &
          & CALL LIST_DESTROY(columnIndicesLists(localDofIdx)%PTR,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 CALL ERRORS("JacobianMatrix_StructureCalculate",err,error)
    CALL EXITS("JacobianMatrix_StructureCalculate")
    RETURN 1
  END SUBROUTINE JacobianMatrix_StructureCalculate

  !
  !================================================================================================================================
  !
 
END MODULE CONSTRAINT_MATRICES_ROUTINES